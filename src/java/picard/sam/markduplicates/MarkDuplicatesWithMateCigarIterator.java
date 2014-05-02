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
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.DiskBackedQueue;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.io.File;
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
    private final List<File> tmpDirs;

    private boolean removeDuplicates = false;
    private boolean skipPairsWithNoMateCigar = true;

    private final Map<String,DuplicationMetrics> metricsByLibrary = new HashMap<String,DuplicationMetrics>();
    private final Map<String,Short> libraryIds = new HashMap<String,Short>();
    private short nextLibraryId = 1;

    private int numRecordsWithNoMateCigar = 0;

    // Variables used for optical duplicate detection and tracking
    private final Histogram<Short> opticalDupesByLibraryId = new Histogram<Short>();

    private boolean foundUnmappedEOFReads = false;
    private int referenceIndex = 0;

    private final DiskBackedQueue<SAMRecord> alignmentStartSortedBuffer;
    private final Set<String> isDuplicateMarkedSet = new HashSet<String>();
    private final MarkQueue toMarkQueue = new MarkQueue();

    private SAMRecord nextRecord = null;
    private OpticalDuplicateFinder opticalDuplicateFinder = null;

    private final SAMRecordCoordinateComparator sortComparator = new SAMRecordCoordinateComparator();

    private ScoringStrategy scoringStrategy = ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH_THEN_MAPQ_THEN_READ_NAME;

    enum ScoringStrategy {
        SUM_OF_BASE_QUALITIES,
        TOTAL_MAPPED_REFERENCE_LENGTH_THEN_MAPQ_THEN_READ_NAME
    }

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
                                               final int maxRecordsInRam,
                                               final List<File> tmpDirs) throws PicardException {
        if (header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new PicardException(this.getClass().getName() + " expects the input to be in coordinate sort order.");
        }

        this.header = header;
        this.backingIterator = new PeekableIterator<SAMRecord>(iterator);
        this.tmpDirs = tmpDirs;
        this.alignmentStartSortedBuffer = DiskBackedQueue.newInstance(new BAMRecordCodec(header), maxRecordsInRam, tmpDirs);

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

        this.toMarkQueue.setToMarkQueueMinimumDistance(toMarkQueueMinimumDistance);

        // get the first record
        this.nextRecord = this.markDuplicatesAndGetTheNextAvailable(); // get one directly, or null
    }

    public String getRecordKey(final SAMRecord record) {return (record.getReadPairedFlag()) ? (record.getReadName() + record.getFirstOfPairFlag() + record.getNotPrimaryAlignmentFlag()) : record.getReadName();}
//    private String getRecordKey(final SAMRecord record) {return record.getReadName() +  (record.getReadPairedFlag() ? record.getFirstOfPairFlag(): 0);}

    public void logMemoryStats(final Log log) {
        System.gc();
        final Runtime runtime = Runtime.getRuntime();
        log.info("freeMemory: " + runtime.freeMemory() +
                "; totalMemory: " + runtime.totalMemory() +
                "; maxMemory: " + runtime.maxMemory() +
                "; buffer size: " + alignmentStartSortedBuffer.size() +
                "; queue size: " + toMarkQueue.size()
        );
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
     * then implementers must throw an IllegalStateException from tmpReadEnds()
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
        if (!this.backingIterator.hasNext() && !this.alignmentStartSortedBuffer.hasNext() && null == this.nextRecord) return false;

        // return if we have tmpReadEnds record
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
                if (!this.alignmentStartSortedBuffer.hasNext()) return null; // no need to flush; no records in either queue or buffer
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
                    this.isDuplicateMarkedSet.add(getRecordKey(record)); // indicate the present wrapped record is available for return
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
                        if (this.alignmentStartSortedBuffer.hasNext()) {
                            System.out.println("ASSB size: " + this.alignmentStartSortedBuffer.size());
                            System.out.println("First record: " + this.alignmentStartSortedBuffer.peek().getSAMString());
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
                if (-1 == this.toMarkQueue.getToMarkQueueMinimumDistance()) {
                    // use twice the first read's length
                    this.toMarkQueue.setToMarkQueueMinimumDistance(Math.max(2*record.getReadBases().length, 100));
                }

                // build a read end for use in the to-mark queue
                readEnds = new ReadEndsMC(header, record);

                // Check that we are not incorrectly performing any duplicate marking, by having too few of the records.  This
                // can happen if the alignment start is increasing but 5' soft-clipping is increasing such that we miss reads with
                // the same 5' unclipped alignment start.
                if (!toMarkQueue.isEmpty()) {
                    final ReadEndsMC end = toMarkQueue.peek();
                    if (end.read1Sequence == readEnds.read1Sequence && this.toMarkQueue.getToMarkQueueMinimumDistance() <= end.read1Coordinate - readEnds.read1Coordinate) {
                        if (checkCigarForSkips(end.getRecord().getCigar())) {
                            throw new PicardException("Found a record with sufficiently large code length that we may have\n"
                                    + " missed including it in an early duplicate marking iteration.  Alignment contains skipped"
                                    + " reference bases (N's). If this is an\n RNAseq aligned bam, please use MarkDuplicates instead,"
                                    + " as this tool does not work well with spliced reads.\n Minimum distance set to "
                                    + this.toMarkQueue.getToMarkQueueMinimumDistance() + " but " + (end.read1Coordinate - readEnds.read1Coordinate - 1)
                                    + " would be required.\n" + "Record was: " + end.getRecord().getSAMString());
                        }
                        else {
                            System.err.println("end: " + end.getRecord().getSAMString());
                            System.err.println("readEnds: " + readEnds.getRecord().getSAMString());
                            throw new PicardException("Found a record with sufficiently large clipping that we may have\n"
                                    + " missed including it in an early duplicate marking iteration.  Please increase the"
                                    + " minimum distance to at least " + (end.read1Coordinate - readEnds.read1Coordinate - 1)
                                    + "bp\nto ensure it is considered (was " + this.toMarkQueue.getToMarkQueueMinimumDistance() + ").\n"
                                    + "Record was: " + end.getRecord().getSAMString());
                        }
                    }
                }

                // do duplicate marking on the available records
                while (!toMarkQueue.isEmpty() &&
                        (this.referenceIndex != readEnds.read1Sequence ||
                                this.toMarkQueue.getToMarkQueueMinimumDistance() < readEnds.read1Coordinate - toMarkQueue.peek().read1Coordinate)) {
                    chunkAndMarkTheDuplicates();
                    performedChunkAndMarkTheDuplicates = true; // indicates we can perhaps find a record if we flush
                    // greedily exit this look if we could flush!
                    if (this.isDuplicateMarkedSet.contains(getRecordKey(this.alignmentStartSortedBuffer.peek()))) break;
                }
            }

            this.backingIterator.next(); // remove the record, since we called this.backingIterator.peek()

            // now wrapped record will be tracked by alignmentStartSortedBuffer until it has been duplicate marked
            this.add(record);

            if (record.isSecondaryOrSupplementary() || record.getReadUnmappedFlag()) { // do not consider these
                  this.isDuplicateMarkedSet.add(getRecordKey(record)); // indicate the present wrapped record is available for return
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
    public int getNumDuplicates() { return this.toMarkQueue.getNumDuplicates(); }
    public Map<String,Short> getLibraryIds() { return this.libraryIds; }
    public Histogram<Short> getOpticalDupesByLibraryId() { return this.opticalDupesByLibraryId; }
    public Map<String,DuplicationMetrics> getMetricsByLibrary() { return this.metricsByLibrary; }

    /**
     * Gets a SAMRecord if one is available after marking.  This enforces that we return records in the original
     * coordinate sort order in a stable fashion.
     */
    private SAMRecord flush() {
        while (this.alignmentStartSortedBuffer.hasNext()) {
            final SAMRecord record = this.alignmentStartSortedBuffer.peek(); // get the first in the buffer
            // check if we can proceed
            if (this.isDuplicateMarkedSet.contains(getRecordKey(record))) {
                this.alignmentStartSortedBuffer.remove(); // remove the record, since we called "peek"
                this.isDuplicateMarkedSet.remove(getRecordKey(record)); // remove the key from the marking tracker
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
}

    /**
     * Chunk the records and mark the duplicates.
     */
    private void chunkAndMarkTheDuplicates()
    {
        if (this.toMarkQueue.isEmpty()) return;

        if (!toMarkQueue.isEmpty() && !this.alignmentStartSortedBuffer.hasNext()) {
            throw new PicardException("0 < toMarkQueue && !alignmentStartSortedBuffer.iterator().hasNext()");
        }

        final ReadEndsMC next = this.toMarkQueue.poll(); // get the first one!
        // it is not marked as a duplicate :)
        this.isDuplicateMarkedSet.add(getRecordKey(next.getRecord())); // track that this record has been through duplicate marking

        // track optical duplicates
        if (next.isPaired() && 0 == next.hasUnmapped) {
            final Set<PhysicalLocationMC> locations = this.toMarkQueue.getLocations(next);
            if (!locations.isEmpty()) {
                    AbstractMarkDuplicateFindingAlgorithm.trackOpticalDuplicates(new ArrayList<PhysicalLocationMC>(locations),
                            this.opticalDuplicateFinder, this.opticalDupesByLibraryId);
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
    private short getLibraryIdFromRecord(final SAMFileHeader header, final SAMRecord rec) {
        final String library = AbstractMarkDuplicateFindingAlgorithm.getLibraryName(header, rec);
        Short libraryId = this.libraryIds.get(library);

        if (libraryId == null) {
            libraryId = this.nextLibraryId++;
            this.libraryIds.put(library, libraryId);
        }

        return libraryId;
    }

    /** We allow different scoring strategies. We return <0 if rec1 has a better strategy than rec2. */
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
        return -cmp;
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
            this.libraryId = getLibraryIdFromRecord(header, this.record);

            // Is this unmapped or its mate?
            if (this.record.getReadUnmappedFlag() || (this.record.getReadPairedFlag() && this.record.getMateUnmappedFlag())) {
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

        // Track that underlying SAMRecord has been through duplicate marking
//        public void setHasBeenMarkedFlag() { super.setHasBeenMarkedFlag(record.getReadName(), true);}
        public SAMRecord getRecord() { return this.record; }
        public SAMRecord setRecord(final SAMRecord record) { return this.record = record; }
        public String getRecordReadName() { return this.record.getReadName(); }

        @Override
        public boolean isPaired() { return this.getRecord().getReadPairedFlag(); }
    }

    /** Stores the minimal information needed for optical duplicate detection. */
    private class PhysicalLocationMC implements OpticalDuplicateFinder.PhysicalLocation {

        // Information used to detect optical dupes
        short readGroup = -1;
        short tile = -1;
        short x = -1, y = -1;
        short libraryId;

        public PhysicalLocationMC(final OpticalDuplicateFinder.PhysicalLocation rec) {
            this.setReadGroup(rec.getReadGroup());
            this.setTile(rec.getTile());
            this.setX(rec.getX());
            this.setY(rec.getY());
            this.setLibraryId(rec.getLibraryId());
        }

        public short getReadGroup() { return this.readGroup; }
        public void  setReadGroup(final short rg) { this.readGroup = rg; }
        public short getTile() { return this.tile; }
        public void  setTile(final short tile) { this.tile = tile; }
        public short getX() { return this.x; }
        public void  setX(final short x) { this.x = x; }
        public short getY() { return this.y; }
        public void  setY(final short y) { this.y = y;}
        public short getLibraryId() { return this.libraryId; }
        public void  setLibraryId(final short libraryId) { this.libraryId = libraryId; }
    }

    /**
     * This is the mark queue.
     *
     * This stores a current set of read ends that need to be duplicate marked.  It only stores internally the "best" read end for a given
     * possible duplicate location, preferring to perform duplicate marking as read ends come in, rather than wait for all "comparable"
     * read ends to arrive.
     */
    private class MarkQueue {

        /** Comparator to order the mark queue set.  The set of all the read ends that have are compared to be the same should
         * be used for duplicate marking. */
        private class MarkQueueComparator implements Comparator<ReadEndsMC> {
            public int compare(final ReadEndsMC lhs, final ReadEndsMC rhs) {
                int retval = lhs.libraryId - rhs.libraryId;
                if (retval == 0) retval = lhs.read1Sequence - rhs.read1Sequence;
                if (retval == 0) retval = lhs.read1Coordinate - rhs.read1Coordinate;
                if (retval == 0) retval = rhs.orientation - lhs.orientation; // to get pairs first
                if (retval == 0) retval = lhs.read2Sequence - rhs.read2Sequence;
                if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
                return retval;
            }
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
                if (retval == 0) retval = lhs.read2Sequence - rhs.read2Sequence;
                if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
                if (retval == 0) retval = compareRecordsByScoringStrategy(lhs.getRecord(), rhs.getRecord());
                if (retval == 0) retval = lhs.getRecordReadName().compareTo(rhs.getRecordReadName());

                return retval;
            }
        }

        private int toMarkQueueMinimumDistance = -1;

        private int numDuplicates = 0;

        /** The set of all read ends sorted by 5' start unclipped position.  Some read ends in this set may eventually be duplicates. */
        private final TreeSet<ReadEndsMC> set = new TreeSet<ReadEndsMC>(new MarkQueueComparator());

        /** Physical locations used for optical duplicate tracking.  This is only stored for paired end reads where both ends are mapped,
         * and when we see the first mate.
         */
        private final Map<ReadEndsMC, Set<PhysicalLocationMC>> locations = new HashMap<ReadEndsMC, Set<PhysicalLocationMC>>();

        /** If we have two items that are the same with respect to being in the "set", then we must choose one.  The "one" will
         * eventually be the end that is not marked as a duplicate in most cases (see poll() for the exceptions).
         */
        private final Comparator<ReadEndsMC> comparator = new ReadEndsMCComparator();

        /** temporary so we do not need to create many objects */
        private ReadEndsMC tmpReadEnds = null;

        public MarkQueue() {
        }

        public int getNumDuplicates() { return this.numDuplicates; }

        public int size() {
            return this.set.size();
        }

        public boolean isEmpty() {
            return this.set.isEmpty();
        }

        public void setToMarkQueueMinimumDistance(final int toMarkQueueMinimumDistance) {
            this.toMarkQueueMinimumDistance = toMarkQueueMinimumDistance;
        }

        public int getToMarkQueueMinimumDistance() { return this.toMarkQueueMinimumDistance; }

        /** For tracking optical duplicates */
        public Set<PhysicalLocationMC> getLocations(final ReadEndsMC current) {
            // NB: only needed for pairs!!!
            if (!current.isPaired() || 0 != current.hasUnmapped) throw new NotImplementedException();
            return this.locations.remove(current);
        }

        public ReadEndsMC peek() {
            return this.set.first();
        }

        /**
         * The poll method will return the read end that is *not* the duplicate of all comparable read ends that
         * have been seen.  All comparable read ends and the returned read end will have their seen duplicate flag set.
         *
         * NB: we must remove all fragments or unpaireds if this is a mapped pair
         * NB: we must remove all fragments if this is an unpaired
         */
        public ReadEndsMC poll() {
            final ReadEndsMC current = this.set.pollFirst();

            // If we are a paired read end, we need to make sure we remove unpaired (if we are not also unpaired), as
            // well as fragments from the set, as they should all be duplicates.
            if (current.isPaired()) {
                // NB: only care about read1Sequence, read1Coordinate, and orientation in the set
                if (null == this.tmpReadEnds) { // initialize
                    this.tmpReadEnds = new ReadEndsMC(header, current.getRecord());
                    this.tmpReadEnds.read2Sequence = this.tmpReadEnds.read2Coordinate = -1;
                    this.tmpReadEnds.record = null;
                }
                else {
                    this.tmpReadEnds.read1Sequence = current.read1Sequence;
                    this.tmpReadEnds.read1Coordinate = current.read1Coordinate;
                }

                // remove from the set fragments and unpaired, which only have two possible orientations
                for (final byte orientation : new byte[]{ReadEnds.F, ReadEnds.R}) { // go through non-paired (both mapped) orientations
                    this.tmpReadEnds.orientation = orientation;
                    if (this.set.contains(this.tmpReadEnds)) { // found in the set
                        // get the duplicate read end
                        final SortedSet<ReadEndsMC> sortedSet = this.set.subSet(this.tmpReadEnds, true, this.tmpReadEnds, true);
                        if (1 != sortedSet.size()) throw new PicardException("SortedSet should have size one (has size " + sortedSet.size() + " )");
                        final ReadEndsMC duplicate = sortedSet.first();

                        // mark as duplicate and set that it has been through duplicate marking
                        duplicate.getRecord().setDuplicateReadFlag(true);
                        isDuplicateMarkedSet.add(getRecordKey(duplicate.getRecord())); // track that this record has been through duplicate marking

                        // remove from the set
                        this.set.remove(this.tmpReadEnds);

                        // count the duplicate metrics
                        final DuplicationMetrics metrics = getMetrics(duplicate.getRecord());
                        // Update the duplication metrics
                        if (!duplicate.getRecord().getReadPairedFlag() || duplicate.getRecord().getMateUnmappedFlag()) {
                            ++metrics.UNPAIRED_READ_DUPLICATES;
                        }
                        else {
                            ++metrics.READ_PAIR_DUPLICATES;// will need to be divided by 2 at the end
                        }
                        this.numDuplicates++;
                    }
                }
            }

            // this read end is now ok to be outputted
            isDuplicateMarkedSet.add(getRecordKey(current.getRecord())); // track that this record has been through duplicate marking

            return current;
        }

        public void add(final ReadEndsMC other) {
            Set<PhysicalLocationMC> locationSet = null;
            ReadEndsMC duplicate = null;

            // 1. check the queue if there exists one at the location, if so compare and keep the best.
            // 2. add physical location info if paired
            if (this.set.contains(other)) { // already in the set
                final SortedSet<ReadEndsMC> sortedSet = this.set.subSet(other, true, other, true);
                if (1 != sortedSet.size()) throw new PicardException("SortedSet should have size one (has size " + sortedSet.size() + " )");
                final ReadEndsMC current = sortedSet.first();

                // compare "current: with "other"
                final int comparison = this.comparator.compare(current, other); // if we are to re-add, then other should make this > 0

                if (0 < comparison) { // re-add
                    if (current.isPaired() && 0 == current.hasUnmapped) {
                        locationSet = this.locations.remove(current);
                    }
                    else {
                        locationSet = new HashSet<PhysicalLocationMC>();
                    }
                    this.locations.put(other, locationSet);
                    this.set.remove(current);
                    this.set.add(other);

                    // current is a now duplicate :/
                    current.getRecord().setDuplicateReadFlag(true);
                    duplicate = current;
                    isDuplicateMarkedSet.add(getRecordKey(current.getRecord())); // track that this record has been through duplicate marking
                }
                else { // keep current
                    if (current.isPaired() && 0 == current.hasUnmapped) {
                        locationSet = this.locations.get(current);
                    }
                    // NB: else is technically not needed, since if this was not paired and the other one was, we would enter here and add it later

                    if (0 != comparison) {
                        // other is a duplicate :/
                        other.getRecord().setDuplicateReadFlag(true);
                        duplicate = other;
                    }
                    isDuplicateMarkedSet.add(getRecordKey(other.getRecord())); // track that this record has been through duplicate marking
                }
            }
            else { // not in the set
                if (other.isPaired() && 0 == other.hasUnmapped) {
                    locationSet = new HashSet<PhysicalLocationMC>();
                    this.locations.put(other, locationSet);
                }
                this.set.add(other);
            }

            // add to the physical locations
            final SAMRecord record = other.getRecord();
            if (record.getReadPairedFlag()
                    && !record.getReadUnmappedFlag()
                    && !record.getMateUnmappedFlag()
                    && record.getFirstOfPairFlag()) { // only first of pairs!
                if (null != locationSet) locationSet.add(new PhysicalLocationMC(other)); // the null check is not needed.
            }

            if (null != duplicate) {
                // count the duplicate metrics
                final DuplicationMetrics metrics = getMetrics(duplicate.getRecord());

                if (duplicate.getRecord().getDuplicateReadFlag()) {
                    // Update the duplication metrics
                    if (!duplicate.getRecord().getReadPairedFlag() || duplicate.getRecord().getMateUnmappedFlag()) {
                        ++metrics.UNPAIRED_READ_DUPLICATES;
                    }
                    else {
                        ++metrics.READ_PAIR_DUPLICATES;// will need to be divided by 2 at the end
                    }
                }
                this.numDuplicates++;
            }
        }
    }

    private class DuplicateMarkingBuffer implements Iterator<SAMRecord> {
        private int maxRecordsInMemory;
        private int availableRecordsInMemory;
        private int blockSize;
        private int queueHeadRecordIndex;
        private int queueTailRecordIndex;
        private Deque<BufferBlock> blocks;

        public DuplicateMarkingBuffer(final int maxRecordsInMemory, final int blockSize) {
            this.maxRecordsInMemory = maxRecordsInMemory;
            this.availableRecordsInMemory = maxRecordsInMemory;
            this.queueHeadRecordIndex = 0;
            this.queueTailRecordIndex = 0;
            this.blockSize = blockSize;
            this.blocks = new ArrayDeque<BufferBlock>();
        }

        public void add(final SAMRecord record, final int recordIndex) {
            // If necessary, create a new block, using as much ram as available up to its total size
            if (this.blocks.size() == 0 || !this.blocks.getLast().isAvailableForWriting()) {
                // once ram is given to a block, we can't give it to another block (until some is recovered from the head of the queue)
                final int blockRam = Math.min(blockSize, availableRecordsInMemory);
                this.availableRecordsInMemory = this.availableRecordsInMemory - blockRam;
                final BufferBlock block = new BufferBlock(this.blockSize, blockRam, tmpDirs);
                this.blocks.add(block);
            }
            this.blocks.getLast().add(record, recordIndex);
        }

        @Override
        public boolean hasNext() { return (blocks.size() != 0); }     // NB: (&& blocks.peekFirst().hasNext()) ??

        public boolean canEmit() { return (this.blocks.size() !=0 && this.blocks.getFirst().canEmit()); }

        @Override
        public SAMRecord next() {
            //TODO reclaim ram if possible

            if (this.blocks.size() == 0)
                throw new NoSuchElementException("Attempting to remove an element from an empty DuplicateMarkingBuffer");
            final BufferBlock head = this.blocks.peekFirst();
            if (!head.hasNext())
                throw new PicardException("DuplicateMarkingBuffer attempted to remove a record from an empty block.");
            if (!head.canEmit())
                throw new PicardException("Attempting to get a record from the DuplicateMarkingBuffer that has not been through " +
                                          "duplicate marking. canEmit() must return true in order to call next()");
            final SAMRecord record = head.next();
            if (!head.hasNext())
                blocks.poll(); // remove the block as it is now empty
            return record;
        }

        @Override
        public void remove() { this.next(); }

        /**
         * Return the total number of elements in the queue, both in memory and on disk
         */
        public int size() { return this.queueTailRecordIndex - this.queueHeadRecordIndex; }

        public void setDuplicateMarkingFlags(final int recordIndex) {
            // TODO
        }

        private class BufferBlock implements Iterator<SAMRecord> {
            private DiskBackedQueue diskBackedQueue;
            private boolean isAvailableForWriting;
            private int maxBlockSize;
            private int startIndex;
            private int endIndex;
            private byte[] wasCheckedIndexes;
            private byte[] isDuplicateIndexes;
            private int maxBlockRecordsInMemory;

            public BufferBlock(final int maxBlockSize, final int maxBlockRecordsInMemory, final List<File> tmpDirs) {
                this.diskBackedQueue = DiskBackedQueue.newInstance(new BAMRecordCodec(header), maxBlockRecordsInMemory, tmpDirs);
                this.maxBlockSize = maxBlockSize;
                this.maxBlockRecordsInMemory = maxBlockRecordsInMemory;
                this.isAvailableForWriting = false;
                this.startIndex = -1;
                this.endIndex = -1;
            }

            public boolean isAvailableForWriting() { return isAvailableForWriting; }

            public void add(final SAMRecord record, final int recordIndex) {
                // TODO
            }

            @Override
            public boolean hasNext() {
                // TODO
                return false;
            }

            @Override
            public SAMRecord next() {
                // TODO
                return null;
            }

            @Override
            public void remove() { this.next(); }

            /**
             * Return the total number of elements in the queue, both in memory and on disk
             */
            public int size() { return this.endIndex - this.startIndex + 1; }

            public void setDuplicateMarkingFlags(final int recordIndex) {
                // TODO
            }

            public boolean canEmit() {
                // TODO
                return false; }
        }
    }
}
