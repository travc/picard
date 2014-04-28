/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam.markduplicates;

import picard.PicardException;
import htsjdk.samtools.util.Histogram;
import picard.sam.DuplicationMetrics;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeEOFException;
import htsjdk.samtools.util.SortingCollection;

import java.util.*;

/**
 * This will iterate through a coordinate sorted SAM file (iterator) and either mark or
 * remove duplicates as appropriate.  This class relies on the coordinate sort order as
 * well as the mate cigar (MC) optional SAM tag.
 */
public class MarkDuplicatesWithMateCigarIterator implements SAMRecordIterator {

    /**
     * Ideas to make things faster:
     * - run the iterator in a different thread
     * - we could try to maintain a minimum buffer of reads that have been duplicate marked
     * - object pool for ReadEndsMC
     */

    private SAMFileHeader header = null;
    private PeekableIterator<SAMRecord> backingIterator = null;

    private boolean removeDuplicates = false;
    private boolean skipPairsWithNoMateCigar = true;

    private int toMarkQueueMinimumDistance = -1;

    private final Map<String,DuplicationMetrics> metricsByLibrary = new HashMap<String,DuplicationMetrics>();
    private final Map<String,Short> libraryIds = new HashMap<String,Short>();
    private short nextLibraryId = 1;

    private int numRecordsWithNoMateCigar = 0;
    private int numDuplicates = 0;

    // Variables used for optical duplicate detection and tracking
    private final Histogram<Short> opticalDupesByLibraryId = new Histogram<Short>();

    private boolean foundUnmappedEOFReads = false;
    private int referenceIndex = 0;

    private final PriorityQueue<ReadEndsMC> toMarkQueue = new PriorityQueue<ReadEndsMC>(100, new ReadEndsMCComparator()); // sorted by 5' start unclipped position
    private final ArrayDeque<SAMRecord> alignmentStartSortedBuffer;
    private final Map<String, Boolean> isDuplicateMarkedMap = new HashMap<String, Boolean>();

    private SAMRecord nextRecord = null;
    private OpticalDuplicateFinder opticalDuplicateFinder = null;

    private final SAMRecordCoordinateComparator sortComparator = new SAMRecordCoordinateComparator();

    private ScoringStrategy scoringStrategy = ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH_THEN_MAPQ_THEN_READ_NAME;

    enum ScoringStrategy {
        SUM_OF_BASE_QUALITIES,
        TOTAL_MAPPED_REFERENCE_LENGTH_THEN_MAPQ_THEN_READ_NAME
    }

    // Convenience variables to be set and used by processChunkSubset and markTheDuplicates.  They are used
    // for when marking fragment reads to check if paired reads were previously found for the same orientation.
    private boolean hasPairsR, hasPairsF; // for convenience, no static variables in Java (go C)

    /**
     * Initializes the mark duplicates iterator
     * @param header the SAM header
     * @param iterator an iterator over the SAM records to consider
     * @param opticalDuplicateFinder the algorithm for optical duplicate detection
     * @param toMarkQueueMinimumDistance minimum distance for which to buf
     * @param removeDuplicates true to remove duplicates, false to mark duplicates
     * @param skipPairsWithNoMateCigar true to not return mapped pairs with no mate cigar, false otherwise
     * @throws PicardException if the inputs are not in coordinate sort order
     */
    public MarkDuplicatesWithMateCigarIterator(final SAMFileHeader header,
                                               final CloseableIterator<SAMRecord> iterator,
                                               final OpticalDuplicateFinder opticalDuplicateFinder,
                                               final int toMarkQueueMinimumDistance,
                                               final boolean removeDuplicates,
                                               final ScoringStrategy scoringStrategy,
                                               final boolean skipPairsWithNoMateCigar,
                                               final int maxRecordsInRam) throws PicardException {
        if (header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new PicardException(this.getClass().getName() + " expects the input to be in coordinate sort order.");
        }

        this.header = header;
        this.backingIterator = new PeekableIterator<SAMRecord>(iterator);
        this.alignmentStartSortedBuffer = new ArrayDeque<SAMRecord>();
//                SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(header), new SAMRecordQueryNameComparator(), maxRecordsInRam);

        this.removeDuplicates = removeDuplicates;
        this.skipPairsWithNoMateCigar = skipPairsWithNoMateCigar;
        this.opticalDuplicateFinder = opticalDuplicateFinder;
        this.scoringStrategy = scoringStrategy;

        // set up metrics
        for(final SAMReadGroupRecord readGroup : this.header.getReadGroups()) {
            final String library = readGroup.getLibrary();
            DuplicationMetrics metrics = this.metricsByLibrary.get(library);
            if (metrics == null) {
                metrics = new DuplicationMetrics();
                metrics.LIBRARY = library;
                this.metricsByLibrary.put(library, metrics);
            }
        }

        this.toMarkQueueMinimumDistance = toMarkQueueMinimumDistance;

        // get the first record
        this.nextRecord = this.markDuplicatesAndGetTheNextAvailable(); // get one directly, or null
    }

    /**
     * Set the scoring strategy to select which records should be called duplicate within a group of comparable records.
     */
    public void setScoringStrategy(final ScoringStrategy scoringStrategy) {
        this.scoringStrategy = scoringStrategy;
    }

    /**
     * Establishes that records returned by this iterator are expected to
     * be in the specified sort order.  If this method has been called,
     * then implementers must throw an IllegalStateException from next()
     * when a record is read that violates the sort order.  This method
     * may be called multiple times over the course of an iteration,
     * changing the expected sort, if desired -- from the time it is called,
     * it validates whatever sort is set, or stops validating if it
     * is set to null or SAMFileHeader.SortOrder.unsorted.  If this method
     * is not called, then no validation of the iterated records is done.
     *
     * @param sortOrder The order in which records are expected to be returned
     * @return  This SAMRecordIterator
     */
    @Override
    public SAMRecordIterator assertSorted(final SAMFileHeader.SortOrder sortOrder) {
        if (sortOrder != SAMFileHeader.SortOrder.coordinate) {
            throw new IllegalStateException("Cannot assort " + sortOrder + " when expecting coordinate sorted input");
        }
        return this;
    }

    @Override
    public boolean hasNext() {

        // fast fail
        if (!this.backingIterator.hasNext() && this.alignmentStartSortedBuffer.isEmpty() && null == this.nextRecord) return false;

        // return if we have next record
        return (null != this.nextRecord);
    }

    @Override
    public SAMRecord next() throws PicardException {
        final SAMRecord toReturn = this.nextRecord; // save for return
        if (hasNext()) { // call hasNext, since we may need it to update this.nextRecord
            this.nextRecord = this.markDuplicatesAndGetTheNextAvailable(); // get one more, if possible
        }
        else {
            this.nextRecord = null;
        }

        // should always return an element
        if (null == toReturn) {
            throw new NoSuchElementException();
        }

        // check for sorted order
        if (null != this.nextRecord && 0 < this.sortComparator.fileOrderCompare(toReturn, nextRecord)) {
            System.err.print("Previous: " + toReturn.getSAMString());
            System.err.print("Current:" + nextRecord.getSAMString());
            throw new PicardException("Records were not found coordinate sort order");
        }

        return toReturn;
    }

    /**
     * This tries to get a record that has been evaluated for duplicate marking.  It does this by first seeing if there
     * are any records that have been through duplicate marking.  If none are available, it will try to get more records
     * from the input iterator until there are reads available that have been duplicate marked.  If there are no more
     * records available from the input iterator, it will duplicate mark the final chunk of records.  Finally, if there
     * are no more records, it will return null;
     */
    private SAMRecord markDuplicatesAndGetTheNextAvailable() {

        // Check if there are any we can flush...
        {
            final SAMRecord record = this.flush();
            if (null != record) return record;
        }

        // Check if there are any more records to read in
        if (!this.backingIterator.hasNext()) { // no more records to read in

            // Check if there are any more to mark
            if (this.toMarkQueue.isEmpty()) {
                if (this.alignmentStartSortedBuffer.isEmpty()) return null; // no need to flush; no records in either queue or buffer
            }
            else {
                // force marking duplicates on the remaining records
                while (!toMarkQueue.isEmpty()) {
                    chunkAndMarkTheDuplicates();
                }
            }

            // update our coordinate to past the end of the reference
            this.referenceIndex = this.header.getSequenceDictionary().getSequences().size();
            return this.markDuplicatesAndGetTheNextAvailable(); // try again, we should either hit flush, or enter in this "if" again and return null
        }

        // We need to retrieve more records from the input iterator and duplicate mark, until we can return one that
        // has been through duplicate marking.
        while (this.backingIterator.hasNext()) {
            // NB: we could get rid of this if we made this.nextRecord into a list...
            SAMRecord record = this.backingIterator.peek(); // peek: used for unmapped reads

            ReadEndsMC readEnds = null;
            boolean performedChunkAndMarkTheDuplicates = false;

            // remove duplicate information
            record.setDuplicateReadFlag(false);

            // ignore/except-on paired records with mapped mate and no mate cigar
            if (record.getReadPairedFlag() &&
                    !record.getMateUnmappedFlag() && null == SAMUtils.getMateCigar(record)) { // paired with one end unmapped and no mate cigar

                // NB: we are not truly examining these records. Do we want to count them?
                if (!record.isSecondaryOrSupplementary()) {
                    // update metrics
                    final DuplicationMetrics metrics = getMetrics(record);
                    if (record.getReadUnmappedFlag()) {
                        ++metrics.UNMAPPED_READS;
                    }
                    else if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()) {
                        ++metrics.UNPAIRED_READS_EXAMINED;
                    }
                    else {
                        ++metrics.READ_PAIRS_EXAMINED;
                    }
                }

                if (this.skipPairsWithNoMateCigar) { // pseudo-silently ignores them
                    // NB: need to add/set-flag as chunking/flushing of the toMarkQueue may need to occur
                    this.add(record); // now wrapped record will be stored in alignmentStartSortedBuffer for return
                    this.setHasBeenMarkedFlag(getReadKey(record), true); // indicate the present wrapped record is available for return
                    this.numRecordsWithNoMateCigar++;
                    this.backingIterator.next(); // remove it, since we called this.backingIterator.peek()
                    continue;
                }
                else {
                    throw new PicardException("Read " + record.getReadName() + " was mapped and had a mapped mate, but no mate cigar (\"MC\") tag.");
                }
            }

            // check for an unmapped read
            if (record.getReadUnmappedFlag()) {

                // unmapped reads at the end of the file!
                if (-1 == record.getReferenceIndex()) {
                    // when we find unmapped reads with -1 as their reference index, there should be no mapped reads later in the file.
                    if (this.foundUnmappedEOFReads) { // previously found unmapped reads at the end of the file
                        final SAMRecord unmappedRecord = this.backingIterator.next(); // since we called this.backingIterator.peek()

                        if (!record.isSecondaryOrSupplementary()) {
                            // update metrics
                            final DuplicationMetrics metrics = getMetrics(record);
                            ++metrics.UNMAPPED_READS;
                        }

                        // We should have no more in the queue
                        if (!this.alignmentStartSortedBuffer.isEmpty()) {
                            throw new PicardException("Encountered unmapped reads at the end of the file, but the alignment start buffer was not empty.");
                        }
                        return unmappedRecord;
                    }
                    else {
                        this.foundUnmappedEOFReads = true;
                        // move past all mapped reads
                        this.referenceIndex = this.header.getSequenceDictionary().getSequences().size();

                        // do the final round of duplicate marking
                        while (!this.toMarkQueue.isEmpty()) {
                            chunkAndMarkTheDuplicates();
                        }

                        return this.markDuplicatesAndGetTheNextAvailable(); // this should flush the buffer
                    }
                }
                if (!record.isSecondaryOrSupplementary()) {
                    // update metrics
                    final DuplicationMetrics metrics = getMetrics(record);
                    ++metrics.UNMAPPED_READS;
                }
                // we will check for unmapped reads later so as not to add them to mark queue
            }
            else {
                if (-1 == this.toMarkQueueMinimumDistance) {
                    // use the first read's length + 1
                    this.toMarkQueueMinimumDistance = Math.max(record.getReadBases().length + 1, 100);

                    // use twice the first read's length
//                    this.toMarkQueueMinimumDistance = Math.max(2*record.getReadBases().length, 100);
                }

                // build a read end for use in the to-mark queue
                readEnds = new ReadEndsMC(header, record);

                // Check that we are not incorrectly performing any duplicate marking, by having too few of the records.  This
                // can happen if the alignment start is increasing but 5' soft-clipping is increasing such that we miss reads with
                // the same 5' unclipped alignment start.
                if (!toMarkQueue.isEmpty()) {
                    final ReadEndsMC end = toMarkQueue.peek();
                    if (end.read1Sequence == readEnds.read1Sequence && toMarkQueueMinimumDistance <= end.read1Coordinate - readEnds.read1Coordinate) {
                        if (checkCigarForSkips(end.getRecord().getCigar())) {
                            throw new PicardException("Found a record with sufficiently large code length that we may have\n"
                                    + " missed including it in an early duplicate marking iteration.  Alignment contains skipped"
                                    + " reference bases (N's). If this is an\n RNAseq aligned bam, please use MarkDuplicates instead,"
                                    + " as this tool does not work well with spliced reads.\n Minimum distance set to "
                                    + this.toMarkQueueMinimumDistance + " but " + (end.read1Coordinate - readEnds.read1Coordinate - 1)
                                    + " would be required.\n" + "Record was: " + end.getRecord().getSAMString());
                        }
                        else {
                            throw new PicardException("Found a record with sufficiently large clipping that we may have\n"
                                    + " missed including it in an early duplicate marking iteration.  Please increase the"
                                    + " minimum distance to at least " + (end.read1Coordinate - readEnds.read1Coordinate - 1)
                                    + "bp\nto ensure it is considered (was " + this.toMarkQueueMinimumDistance + ").\n"
                                    + "Record was: " + end.getRecord().getSAMString());
                        }
                    }
                }

                // do duplicate marking on the available records
                while (!toMarkQueue.isEmpty() &&
                        (this.referenceIndex != readEnds.read1Sequence ||
                                toMarkQueueMinimumDistance < readEnds.read1Coordinate - toMarkQueue.peek().read1Coordinate)) {
                    chunkAndMarkTheDuplicates();
                    performedChunkAndMarkTheDuplicates = true; // indicates we can perhaps find a record if we flush
                }
            }

            this.backingIterator.next(); // remove the record, since we called this.backingIterator.peek()

            // now wrapped record will be tracked by alignmentStartSortedBuffer until it has been duplicate marked
            this.add(record);

            if (record.isSecondaryOrSupplementary() || record.getReadUnmappedFlag()) { // do not consider these
                  this.setHasBeenMarkedFlag(getReadKey(record), true); // indicate the present wrapped record is available for return
            }
            else {
                // update metrics
                final DuplicationMetrics metrics = getMetrics(record);

                // Bring the simple metrics up to date
                if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()) {
                    ++metrics.UNPAIRED_READS_EXAMINED;
                }
                else {
                    ++metrics.READ_PAIRS_EXAMINED; // will need to be divided by 2 at the end
                }

                // add for duplicate marking
                toMarkQueue.add(readEnds);
            }

            // Check if there are any we can flush, which happens if we just performed duplicate marking
            if (performedChunkAndMarkTheDuplicates) {
                record = this.flush();
                if (null != record) return record;
            }
        }

        // try again, as we may have records we can flush, or we want to see if we are at the EOF
        return this.markDuplicatesAndGetTheNextAvailable();
    }

    private String getReadKey(final SAMRecord record) {return record.getReadName() +  (record.getReadPairedFlag() ? record.getFirstOfPairFlag(): 0);}

    private void setHasBeenMarkedFlag(final String readName, final boolean wasDuplicateMarked) {
        final Boolean flag = this.isDuplicateMarkedMap.get(readName);
        if (flag == null) {
            throw new RuntimeException("Error: Attempting to set duplicate marking flag for un-tracked readName " + readName);
        }
        this.isDuplicateMarkedMap.put(readName, wasDuplicateMarked);
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void close() {
        this.backingIterator.close();
    }

    /**
     * Checks a Cigar for the presence of N operators. Reads with skipped bases may be spliced RNAseq reads
     * @param cigar
     */
    private boolean checkCigarForSkips(final Cigar cigar) {
        final List<CigarElement> elements = cigar.getCigarElements();
        for (final CigarElement el : elements)
            if (el.getOperator() == CigarOperator.N) return true;
        return false;
    }

    /** Useful for statistics after the iterator has completed */
    public int getNumRecordsWithNoMateCigar() { return this.numRecordsWithNoMateCigar; }
    public int getNumDuplicates() { return this.numDuplicates; }
    public Map<String,Short> getLibraryIds() { return this.libraryIds; }
    public Histogram<Short> getOpticalDupesByLibraryId() { return this.opticalDupesByLibraryId; }
    public Map<String,DuplicationMetrics> getMetricsByLibrary() { return this.metricsByLibrary; }

    /**
     * Gets a SAMRecord if one is available after marking.  This enforces that we return records in the original
     * coordinate sort order in a stable fashion.
     */
    private SAMRecord flush() {
        while (!this.alignmentStartSortedBuffer.isEmpty()) {
            final SAMRecord record = this.alignmentStartSortedBuffer.getFirst(); // get the first in the buffer
            // check if we can proceed
            if (this.isDuplicateMarkedMap.get(getReadKey(record))) {
                this.alignmentStartSortedBuffer.removeFirst(); // remove the record, since we called "get"
                this.isDuplicateMarkedMap.remove(getReadKey(record)); // remove the key from the marking tracker
                // remove or just mark?
                if (!this.removeDuplicates || !record.getDuplicateReadFlag()) {
                    // we do not want to remove or it is not marked as a duplicate
                    return record; // return the underlying SAMRecord for output from the iterator
                }
            }
            else { // nope, we have not duplicate marked this read
                break;
            }
        }
        return null;
    }

    /**
     * Adds a record to the alignment start buffer
     * @param record - a SAMRecord wrapped alongside a boolean tracking whether it has been through duplicate marking
     * @throws PicardException
     */
    private void add(final SAMRecord record) throws PicardException {
        final int recordReferenceIndex = record.getReferenceIndex();
        if (recordReferenceIndex < this.referenceIndex) {
            throw new PicardException("Records out of order: " + recordReferenceIndex + " < " + this.referenceIndex);
        }
        else if (this.referenceIndex < recordReferenceIndex) {
            // new reference, so we need to mark duplicates on the current ones
            while (!this.toMarkQueue.isEmpty()) {
                chunkAndMarkTheDuplicates();
            }
            // update genomic coordinate
            this.referenceIndex = recordReferenceIndex;
        }

        // add the wrapped record to the list
        this.alignmentStartSortedBuffer.add(record);
        this.isDuplicateMarkedMap.put(getReadKey(record), false);
    }

    /**
     * This performs duplicate marking on a subset of chunk: [startChunkIdx, endChunkIdx] (inclusive).
     */
    private void markTheDuplicates(final List<ReadEndsMC> chunk,
                                   final boolean isPaired) {
        boolean prevWasBest = false; // for pairs, this is the default

        if (chunk.isEmpty()) return; // should never happen

        final ReadEndsMC firstReadEndsMC = chunk.get(0);

        // we mark all as duplicates if a valid pair read was found previously
        if (!isPaired) {
            prevWasBest = (ReadEndsMC.F == firstReadEndsMC.orientation && hasPairsF)
                    || (ReadEndsMC.R == firstReadEndsMC.orientation && hasPairsR);
        }

        // must be more than one to mark duplicates
        if (prevWasBest || 1 < chunk.size()) {
            ReadEndsMC best = null;
            int totalNumberToConsider = 0;

            if (!prevWasBest) { // should we find a 'best'?
                // check that we do not have all one paired end with both ends mapped to the same 5' location

                // get the best scoring end
                for (final ReadEndsMC end : chunk) {
                    // ignore unmapped reads
                    if (end.getRecord().getReadUnmappedFlag()) continue;

                    // get the best score end
                    if (null == best || compareRecordsByScoringStrategy(best.getRecord(), end.getRecord()) < 0) {
                        best = end;
                    }
                    totalNumberToConsider++;
                }
            }

            // do not mark duplicates if we have only one or all the read names are the same
            if (prevWasBest || (1 < totalNumberToConsider)) {
                // annotate the record(s) and decrement the coordinate counts
                for (final ReadEndsMC end : chunk) {
                    // to duplicate mark, or not to duplicate mark, that is the question
                    // NB: do not mark records that have the same read name as the "best"
                    if (null != best && (best == end || best.getRecordReadName().equals(end.getRecordReadName()))) {
                        // not a duplicate
                        end.getRecord().setDuplicateReadFlag(false);
                    }
                    else if (!end.getRecord().getReadUnmappedFlag()) { // do not mark unmapped reads as duplicates
                        // duplicate
                        end.getRecord().setDuplicateReadFlag(true);
                    }
                }
            }
        }

        // optical duplicates and set flags for if future fragments should all be duplicates
        if (isPaired) { // only for pairs
            // optical duplicates
            if (1 < chunk.size()) {
                final List<ReadEndsMC> firstOfPairs = new LinkedList<ReadEndsMC>();
                for (final ReadEndsMC end : chunk) {
                    if (end.isPaired() && end.getRecord().getFirstOfPairFlag()) firstOfPairs.add(end);
                }
                AbstractMarkDuplicateFindingAlgorithm.trackOpticalDuplicates(firstOfPairs,
                        this.opticalDuplicateFinder, this.opticalDupesByLibraryId);
            }
            // for future fragment reads
            if (firstReadEndsMC.getRecord().getReadNegativeStrandFlag()) this.hasPairsR = true;
            else this.hasPairsF = true;
        }
    }

    /**
     * We only mark duplicates on groups of reads with the same library, genomic coordinate (5' unclipped), and orientation
     * (mapped strand).  In the case of pairs, we in addition consider the mate's genomic coordinate (5' unclipped).  We
     * also provide the ability to ignore orientation as well as the mate's genomic coordinate.
     */
    private boolean areComparableForDuplicates(final ReadEndsMC lhs, final ReadEndsMC rhs, final boolean withOrientation, final boolean withSecondEnd) {
        boolean retval = lhs.libraryId == rhs.libraryId &&
                lhs.read1Sequence == rhs.read1Sequence &&
                lhs.read1Coordinate == rhs.read1Coordinate &&
                (!withOrientation || lhs.orientation == rhs.orientation);
        if (withSecondEnd && lhs.record.getReadPairedFlag() && rhs.record.getReadPairedFlag()) {
            retval = retval &&
                    lhs.read2Sequence   == rhs.read2Sequence &&
                    lhs.read2Coordinate == rhs.read2Coordinate;
        }
        return retval;
    }

    /**
     * This will split the records based on orientation, and then run mark duplicates.  We assumed this is called twice
     * in succession, once for a chunk of paired reads, and once for a chunk of fragment/unpaired reads.  The chunk
     * provided are not divided by orientation (with areComparableForDuplicates).      */
    private void processChunkSubset(final List<ReadEndsMC> chunk,
                                    final boolean isPaired) {
        ReadEndsMC start = null;
        int currentChunkIdx = 0;
        int startChunkIdx = 0;

        // go through the chunk.  We wish to split further split them by orientation.
        for (final ReadEndsMC next : chunk) {
            // update
            if (null == start) {
                start = next;
                startChunkIdx = currentChunkIdx;
            }
            else if (!areComparableForDuplicates(start, next, true, isPaired)) { // now include orientation
                // mark duplicates
                markTheDuplicates(chunk.subList(startChunkIdx, currentChunkIdx), isPaired);

                // move to the next chunk
                startChunkIdx = currentChunkIdx;
                start = next;
            }
            // next
            currentChunkIdx++;
        }
        // mark the last chunk, always
        markTheDuplicates(chunk.subList(startChunkIdx, currentChunkIdx), isPaired);
    }

    /**
     * Chunk the records and mark the duplicates.
     */
    private void chunkAndMarkTheDuplicates()
    {
        final List<ReadEndsMC> chunk = new ArrayList<ReadEndsMC>();
        int lastPairBothMappedIndex = -1; // zero-based

        if (!toMarkQueue.isEmpty() && this.alignmentStartSortedBuffer.isEmpty()) {
            throw new PicardException("0 < toMarkQueue && !alignmentStartSortedBuffer.iterator().hasNext()");
        }

        // Get all records at this position.
        // We do not care about the orientation or second end yet.
        while (!toMarkQueue.isEmpty()) {
            final ReadEndsMC next;

            if (chunk.isEmpty() || areComparableForDuplicates(chunk.get(0), toMarkQueue.peek(), false, false)) {
                next = toMarkQueue.poll();
            }
            else {
                break;
            }
            chunk.add(next);
            this.setHasBeenMarkedFlag(getReadKey(next.getRecord()), true); // track that this record has been through duplicate marking

            // optimizing for later loops
            if (next.isPaired() && !next.getRecord().getReadUnmappedFlag() && !next.getRecord().getMateUnmappedFlag()) {
                lastPairBothMappedIndex++;
            }
        }

        // no duplicate marking necessary on one record
        if (1 < chunk.size()) {

            // reset
            this.hasPairsF = this.hasPairsR = false;

            // process the pairs (both ends mapped) first
            if (0 <= lastPairBothMappedIndex) {
                processChunkSubset(chunk.subList(0, lastPairBothMappedIndex+1), true);
            }

            // next come the unpaired and fragments
            if (lastPairBothMappedIndex+1 < chunk.size()) {
                processChunkSubset(chunk.subList(lastPairBothMappedIndex+1, chunk.size()), false);
            }
        }

        // count the duplicate metrics
        for (final ReadEndsMC end : chunk) {
            final DuplicationMetrics metrics = getMetrics(end.getRecord());

            if (end.getRecord().getDuplicateReadFlag()) {
                // duplicate
                this.numDuplicates++;

                // Update the duplication metrics
                if (!end.getRecord().getReadPairedFlag() || end.getRecord().getMateUnmappedFlag()) {
                    ++metrics.UNPAIRED_READ_DUPLICATES;
                }
                else {
                    ++metrics.READ_PAIR_DUPLICATES;// will need to be divided by 2 at the end
                }
            }
        }
    }

    /** Get the duplication metrics for the library associated with end. */
    private DuplicationMetrics getMetrics(final SAMRecord record) {
        final String library = AbstractMarkDuplicateFindingAlgorithm.getLibraryName(this.header, record);
        DuplicationMetrics metrics = this.metricsByLibrary.get(library);
        if (metrics == null) {
            metrics = new DuplicationMetrics();
            metrics.LIBRARY = library;
            this.metricsByLibrary.put(library, metrics);
        }
        return metrics;
    }

    /** Get the library ID for the given SAM record. */
    private short getLibraryId(final SAMFileHeader header, final SAMRecord rec) {
        final String library = AbstractMarkDuplicateFindingAlgorithm.getLibraryName(header, rec);
        Short libraryId = this.libraryIds.get(library);

        if (libraryId == null) {
            libraryId = this.nextLibraryId++;
            this.libraryIds.put(library, libraryId);
        }

        return libraryId;
    }

    /** We allow different scoring strategies. Larger is better. */
    private int compareRecordsByScoringStrategy(final SAMRecord rec1, final SAMRecord rec2) {
        int cmp = 0;
        int referenceLength1, referenceLength2;

        // always prefer paired over non-paired
        if (rec1.getReadPairedFlag() != rec2.getReadPairedFlag()) return rec1.getReadPairedFlag() ? 1 : -1;

        switch (this.scoringStrategy) {
            case SUM_OF_BASE_QUALITIES:
                // NB: we could cache this somewhere if necessary
                cmp = getSumOfBaseQualities(rec1) - getSumOfBaseQualities(rec2);
                break;
            case TOTAL_MAPPED_REFERENCE_LENGTH_THEN_MAPQ_THEN_READ_NAME:
                // get the reference length, both ends if paired
                referenceLength1 = rec1.getCigar().getReferenceLength();
                referenceLength2 = rec2.getCigar().getReferenceLength();
                if (rec1.getReadPairedFlag() && !rec1.getMateUnmappedFlag()) referenceLength1 += SAMUtils.getMateCigar(rec1).getReferenceLength();
                if (rec2.getReadPairedFlag() && !rec2.getMateUnmappedFlag()) referenceLength2 += SAMUtils.getMateCigar(rec2).getReferenceLength();
                // compare
                cmp = referenceLength1 - referenceLength2;
                if (0 == cmp) cmp = rec1.getMappingQuality() - rec2.getMappingQuality();
                if (0 == cmp) cmp = rec1.getReadName().compareTo(rec2.getReadName());
                break;
        }
        return cmp;
    }

    /** Calculates a score for the read which is the sum of scores over Q15. */
    private short getSumOfBaseQualities(final SAMRecord rec) {
        short score = 0;
        for (final byte b : rec.getBaseQualities()) {
            if (b >= 15) score += b;
        }

        return score;
    }


    // NB: could refactor to use ReadEndsMC.java
    private class ReadEndsMC extends ReadEnds {
        // to see if we either end is unmapped
        byte hasUnmapped = 0;

        // we need this reference so we can access the mate cigar among other things
        private SAMRecord record = null;


        /** Builds a read ends object that represents a single read. */
        public ReadEndsMC(final SAMFileHeader header, final SAMRecord rec) {
            this.readGroup = -1;
            this.tile = -1;
            this.x = this.y = -1;
            this.read2Sequence = this.read2Coordinate = -1;
            this.hasUnmapped = 0;

            this.record = rec;

            this.read1Sequence = this.record.getReferenceIndex();
            this.read1Coordinate = this.record.getReadNegativeStrandFlag() ? this.record.getUnclippedEnd() : this.record.getUnclippedStart();
            if (this.record.getReadUnmappedFlag()) {
                throw new PicardException("Found an unexpected unmapped read");
            }

            if (this.record.getReadPairedFlag() && !this.record.getReadUnmappedFlag() && !this.record.getMateUnmappedFlag()) {
                this.read2Sequence = this.record.getMateReferenceIndex();
                this.read2Coordinate = this.record.getMateNegativeStrandFlag() ? SAMUtils.getMateUnclippedEnd(this.record) : SAMUtils.getMateUnclippedStart(this.record);

                // set orientation
                this.orientation = ReadEnds.getOrientationByte(this.record.getReadNegativeStrandFlag(), this.record.getMateNegativeStrandFlag());
            }
            else {
                this.orientation = this.record.getReadNegativeStrandFlag() ? ReadEndsMC.R : ReadEndsMC.F;
            }

            // Fill in the library ID
            this.libraryId = getLibraryId(header, this.record);

            // Is this unmapped or its mate?
            if (this.record.getReadUnmappedFlag() || (this.record.getReadPairedFlag() && this.record.getMateNegativeStrandFlag())) {
                this.hasUnmapped = 1;
            }

            // Fill in the location information for optical duplicates
            if (opticalDuplicateFinder.addLocationInformation(this.record.getReadName(), this)) {
                // calculate the RG number (nth in list)
                // NB: could this be faster if we used a hash?
                this.readGroup = 0;
                final String rg = (String) this.record.getAttribute("RG");
                final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
                if (rg != null && readGroups != null) {
                    for (final SAMReadGroupRecord readGroup : readGroups) {
                        if (readGroup.getReadGroupId().equals(rg)) break;
                        else this.readGroup++;
                    }
                }
            }
        }

        // track that underlying SAMRecord has been through duplicate marking
//        public void setHasBeenMarkedFlag() { super.setHasBeenMarkedFlag(record.getReadName(), true);}
        public SAMRecord getRecord() { return this.record; }
        public SAMRecord setRecord(final SAMRecord record) { return this.record = record; }
        public String getRecordReadName() { return this.record.getReadName(); }
    }

    /** Comparator for ReadEndsMC that orders by read1 position then pair orientation then read2 position. */
    // Could be a Singleton, but no static variables in inner classes.
    class ReadEndsMCComparator implements Comparator<ReadEndsMC> {
        public int compare(final ReadEndsMC lhs, final ReadEndsMC rhs) {
            int retval = lhs.libraryId - rhs.libraryId;
            if (retval == 0) retval = lhs.read1Sequence - rhs.read1Sequence;
            if (retval == 0) retval = lhs.read1Coordinate - rhs.read1Coordinate;
            if (retval == 0) retval = rhs.orientation - lhs.orientation; // IMPORTANT: reverse the order to get pairs first
            if (retval == 0 && lhs.isPaired() != rhs.isPaired()) return lhs.isPaired() ? -1 : 1; // unpaired goes first...
            if (retval == 0) retval = lhs.hasUnmapped - rhs.hasUnmapped;
            if (retval == 0) retval = lhs.read2Sequence   - rhs.read2Sequence;
            if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
            if (retval == 0) retval = lhs.getRecordReadName().compareTo(rhs.getRecordReadName());

            return retval;
        }
    }
}
