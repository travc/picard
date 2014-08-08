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
import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.util.CloseableIterator;
import picard.sam.markduplicates.util.DuplicateMarkingBuffer;
import picard.sam.markduplicates.util.*;

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
    private int backingIteratorRecordIndex = 0;

    private boolean removeDuplicates = false;
    private boolean skipPairsWithNoMateCigar = true;

    private final LibraryIdGenerator libraryIdGenerator;

    private int numRecordsWithNoMateCigar = 0;

    private boolean foundUnmappedEOFReads = false;
    private int referenceIndex = 0;

    private DuplicateMarkingBuffer alignmentStartSortedBuffer = null;
    private final Set<String> isDuplicateMarkedSet = new HashSet<String>();
    private final MarkQueue toMarkQueue;

    private SAMRecord nextRecord = null;
    private OpticalDuplicateFinder opticalDuplicateFinder = null;

    private final SAMRecordCoordinateComparator sortComparator = new SAMRecordCoordinateComparator();

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
                                               final ScoringStrategy duplicateScoringStrategy,
                                               final int toMarkQueueMinimumDistance,
                                               final boolean removeDuplicates,
                                               final boolean skipPairsWithNoMateCigar,
                                               final int maxRecordsInRam,
                                               final int blockSize,
                                               final List<File> tmpDirs) throws PicardException {
        if (header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new PicardException(this.getClass().getName() + " expects the input to be in coordinate sort order.");
        }

        this.header = header;
        this.backingIterator = new PeekableIterator<SAMRecord>(iterator);
        this.alignmentStartSortedBuffer = new DuplicateMarkingBuffer(maxRecordsInRam, blockSize, tmpDirs, header);

        this.removeDuplicates = removeDuplicates;
        this.skipPairsWithNoMateCigar = skipPairsWithNoMateCigar;
        this.opticalDuplicateFinder = opticalDuplicateFinder;

        // Check for supported scoring strategies
        switch (duplicateScoringStrategy) {
            case SUM_OF_BASE_QUALITIES:
                throw new PicardException("SUM_OF_BASE_QUALITIES not supported as this may cause inconsistencies across ends in a pair.  Please use a different scoring strategy.");
            default:
                break;
        }

        toMarkQueue = new MarkQueue(duplicateScoringStrategy);
        libraryIdGenerator = new LibraryIdGenerator(header);

        // set up metrics
        for(final SAMReadGroupRecord readGroup : this.header.getReadGroups()) {
            final String library = readGroup.getLibrary();
            DuplicationMetrics metrics = this.libraryIdGenerator.getMetricsByLibrary(library);
            if (metrics == null) {
                metrics = new DuplicationMetrics();
                metrics.LIBRARY = library;
                this.libraryIdGenerator.addMetricsByLibrary(library, metrics);
            }
        }

        this.toMarkQueue.setToMarkQueueMinimumDistance(toMarkQueueMinimumDistance);

        // get the first samRecordIndex
        this.nextRecord = this.markDuplicatesAndGetTheNextAvailable(); // get one directly, or null
    }

//    public String getRecordKey(final SAMRecord samRecordIndex) {return (samRecordIndex.getReadPairedFlag()) ? (samRecordIndex.getReadName() + samRecordIndex.getFirstOfPairFlag() + samRecordIndex.getNotPrimaryAlignmentFlag()) : samRecordIndex.getReadName();}
//    private String getRecordKey(final SAMRecord samRecordIndex) {return samRecordIndex.getReadName() +  (samRecordIndex.getReadPairedFlag() ? samRecordIndex.getFirstOfPairFlag(): 0);}

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
     * Establishes that records returned by this iterator are expected to
     * be in the specified sort order.  If this method has been called,
     * then implementers must throw an IllegalStateException from tmpReadEnds()
     * when a samRecordIndex is read that violates the sort order.  This method
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

        // fast fail      TODO- why are we checking this.nextRecord twice?
        if (!this.backingIterator.hasNext() && this.alignmentStartSortedBuffer.isEmpty() && null == this.nextRecord) return false;

        // return if we have tmpReadEnds samRecordIndex
        return (null != this.nextRecord);
    }

    @Override
    public SAMRecord next() throws PicardException {
        final SAMRecord toReturn = this.nextRecord; // save for return
        if (hasNext()) { // call hasNext, since we may need it to update this.nextRecord
            this.nextRecord = this.markDuplicatesAndGetTheNextAvailable(); // get one more, if possible
        } else {
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

        // Check if there are any we can flush from the mark duplicates queue
        {
            final SAMRecord record = this.flush();
            if (null != record) return record;
        }

        // Check if there are any more records to read in
        if (!this.backingIterator.hasNext()) { // no more records to read in

            // Check if there are any more to mark
            if (this.toMarkQueue.isEmpty()) {
                if (this.alignmentStartSortedBuffer.isEmpty()) {
                    return null;
                } // no need to flush; no records in either queue or buffer
            }
            else {
                // force marking duplicates on the remaining records  TODO- are these all non-duplicates???
                while (!toMarkQueue.isEmpty()) {
                    chunkAndMarkTheDuplicates();
                }
            }
            /*
            System.err.println("toMarkQueue.isEmpty()=" + toMarkQueue.isEmpty());
            System.err.println("this.alignmentStartSortedBuffer.isEmpty()=" + this.alignmentStartSortedBuffer.isEmpty());
            System.err.println("this.alignmentStartSortedBuffer.canEmit()=" + this.alignmentStartSortedBuffer.canEmit());
            */

            // update our coordinate to past the end of the reference
            this.referenceIndex = this.header.getSequenceDictionary().getSequences().size();
            return this.markDuplicatesAndGetTheNextAvailable(); // try again, we should either hit flush, or enter in this "if" again and return null
        }

        // We need to retrieve more records from the input iterator and duplicate mark, until we can return one that
        // has been through duplicate marking.
        while (this.backingIterator.hasNext()) {
            // NB: we could get rid of this if we made this.nextRecord into a list...
            SAMRecord record = this.backingIterator.peek(); // peek: used for unmapped reads
            final SamRecordIndex samRecordIndex = new SamRecordIndex(record, this.backingIteratorRecordIndex);

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
                    this.add(samRecordIndex); // now samRecordIndex will be stored in alignmentStartSortedBuffer for return
                    this.backingIteratorRecordIndex++;
                    this.alignmentStartSortedBuffer.setDuplicateMarkingFlags(samRecordIndex, false); // indicate the present wrapped samRecordIndex is available for return
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
                        return unmappedRecord; // unmapped end of file records can simply be emitted - no need to duplicate mark them
                    }
                    else {
                        this.foundUnmappedEOFReads = true;
                        // move past all mapped reads
                        this.referenceIndex = this.header.getSequenceDictionary().getSequences().size();

                        // do the final round of duplicate marking
                        while (!this.toMarkQueue.isEmpty()) {
//                            System.out.println("chunking and marking at EOF reads");
                            chunkAndMarkTheDuplicates();
                        }
                        // NB: we do not call next here since we will recurse and perhaps hit the flush, or re-enter the if with unmapped EOF reads
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
                readEnds = new ReadEndsMC(header, samRecordIndex, opticalDuplicateFinder, libraryIdGenerator.getLibraryIdFromRecord(samRecordIndex.getRecord()));

                // Check that we are not incorrectly performing any duplicate marking, by having too few of the records.  This
                // can happen if the alignment start is increasing but 5' soft-clipping is increasing such that we miss reads with
                // the same 5' unclipped alignment start.
                if (!toMarkQueue.isEmpty()) {
                    final ReadEndsMC end = toMarkQueue.peek();
                    if (end.read1Sequence == readEnds.read1Sequence && this.toMarkQueue.getToMarkQueueMinimumDistance() <= end.read1Coordinate - readEnds.read1Coordinate) {
                        if (checkCigarForSkips(end.getRecord().getCigar())) {
                            throw new PicardException("Found a samRecordIndex with sufficiently large code length that we may have\n"
                                    + " missed including it in an early duplicate marking iteration.  Alignment contains skipped"
                                    + " reference bases (N's). If this is an\n RNAseq aligned bam, please use MarkDuplicates instead,"
                                    + " as this tool does not work well with spliced reads.\n Minimum distance set to "
                                    + this.toMarkQueue.getToMarkQueueMinimumDistance() + " but " + (end.read1Coordinate - readEnds.read1Coordinate - 1)
                                    + " would be required.\n" + "Record was: " + end.getRecord().getSAMString());
                        }
                        else {
                            System.err.println("end: " + end.getRecord().getSAMString());
                            System.err.println("readEnds: " + readEnds.getRecord().getSAMString());
                            throw new PicardException("Found a samRecordIndex with sufficiently large clipping that we may have\n"
                                    + " missed including it in an early duplicate marking iteration.  Please increase the"
                                    + " minimum distance to at least " + (end.read1Coordinate - readEnds.read1Coordinate - 1)
                                    + "bp\nto ensure it is considered (was " + this.toMarkQueue.getToMarkQueueMinimumDistance() + ").\n"
                                    + "Record was: " + end.getRecord().getSAMString());
                        }
                    }
                }

                // do duplicate marking on the available records            TODO- NOTE- no next in here either
                while (!toMarkQueue.isEmpty() &&
                        (this.referenceIndex != readEnds.read1Sequence ||
                                this.toMarkQueue.getToMarkQueueMinimumDistance() < readEnds.read1Coordinate - toMarkQueue.peek().read1Coordinate)) {
//                    System.out.println("chunking and marking bc new record at " + readEnds.read1Coordinate + " is > " + this.toMarkQueue.toMarkQueueMinimumDistance + " away from head record in TMQ at " + toMarkQueue.peek().read1Coordinate);
                    chunkAndMarkTheDuplicates();
                    performedChunkAndMarkTheDuplicates = true; // indicates we can perhaps find a samRecordIndex if we flush
                    // greedily exit this look if we could flush!
//                    if  (this.alignmentStartSortedBuffer.canEmit()) break;
//                    if (this.isDuplicateMarkedSet.contains(getRecordKey(this.alignmentStartSortedBuffer.peek()))) break;
                }
            }

            this.backingIterator.next(); // remove the samRecordIndex, since we called this.backingIterator.peek()

            // now wrapped samRecordIndex will be tracked by alignmentStartSortedBuffer until it has been duplicate marked
            this.add(samRecordIndex);
            this.backingIteratorRecordIndex++;

            // We do not consider these. Indicate the present samRecordIndex is available for return
            if (record.isSecondaryOrSupplementary() || record.getReadUnmappedFlag()) {
                this.alignmentStartSortedBuffer.setDuplicateMarkingFlags(samRecordIndex, false);
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
                toMarkQueue.add(readEnds, alignmentStartSortedBuffer, getMetrics(readEnds.getRecord()));
//                System.out.println("adding record to TMQ at " + readEnds.read1Coordinate);
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
        this.alignmentStartSortedBuffer.close();
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
    // TODO: enforce that these cannot be accessed until the iterator has been closed
    public int getNumRecordsWithNoMateCigar() { return this.numRecordsWithNoMateCigar; }
    public int getNumDuplicates() { return this.toMarkQueue.getNumDuplicates(); }
    public LibraryIdGenerator getLibraryIdGenerator() { return this.libraryIdGenerator; }
    public Map<String,Short> getLibraryIds() { return this.libraryIdGenerator.getLibraryIdsMap(); }
    public Histogram<Short> getOpticalDupesByLibraryId() { return this.libraryIdGenerator.getOpticalDupesByLibraryIdMap(); }
    public Map<String,DuplicationMetrics> getMetricsByLibrary() { return this.libraryIdGenerator.getMetricsByLibraryMap(); }

    /**
     * Gets a SAMRecord if one is available after marking.  This enforces that we return records in the original
     * coordinate sort order in a stable fashion.
     *
     * @return record representing the head of the alignment-start sorted buffer, or null if the head record has not yet been duplicate marked
     */
    private SAMRecord flush() {
        // Check that there is at least one samRecordIndex in the coordinate-sorted buffer, and that the head record has been through duplicate-marking
        while (!this.alignmentStartSortedBuffer.isEmpty() && this.alignmentStartSortedBuffer.canEmit()) {
            // the buffer contains wrapped SAMRecords, which we want to unwrap
            final SAMRecord record = this.alignmentStartSortedBuffer.next().getRecord();

            // If this read is a duplicate, do we want to remove it (continue the loop) or return it for emission?
            if (!this.removeDuplicates || !record.getDuplicateReadFlag()) {
                return record;
            }
        }
        return null;
    }

    /**
     * Adds a samRecordIndex to the alignment start buffer
     * @param samRecordIndex - a SAMRecord wrapped alongside a boolean tracking whether it has been through duplicate marking
     * @throws PicardException
     */
    private void add(final SamRecordIndex samRecordIndex) throws PicardException {
        final int recordReferenceIndex = samRecordIndex.getRecord().getReferenceIndex();
        if (recordReferenceIndex < this.referenceIndex) {
            throw new PicardException("Records out of order: " + recordReferenceIndex + " < " + this.referenceIndex);
        }
        else if (this.referenceIndex < recordReferenceIndex) {
            // new reference, so we need to mark duplicates on the current ones
            while (!this.toMarkQueue.isEmpty()) {
//                System.out.println("chunking and marking from record add, because we're on a new chromosome");
                chunkAndMarkTheDuplicates();           // TODO- seems like this is just going to process remaining records out as non-duplicates. Are we missing inter-chromosomal duplicates here????
            }
            // update genomic coordinate
            this.referenceIndex = recordReferenceIndex;
        }

        // add the wrapped samRecordIndex to the list
        this.alignmentStartSortedBuffer.add(samRecordIndex);
}

    /**
     * Chunk the records and mark the duplicates.      //TODO - rename this function. really we are just processing out reads we're calling non-duplicates
     */
    private void chunkAndMarkTheDuplicates()
    {   //TODO: are we always calling this in a loop????? If so, move the looping logic in here.
        if (this.toMarkQueue.isEmpty()) return;

        if (!toMarkQueue.isEmpty() && this.alignmentStartSortedBuffer.isEmpty()) {
            throw new PicardException("0 < toMarkQueue && alignmentStartSortedBuffer.isEmpty()");
        }
        // Poll will track that this samRecordIndex has been through duplicate marking. It is not marked as a duplicate :)
        final ReadEndsMC next = this.toMarkQueue.poll(alignmentStartSortedBuffer, header, opticalDuplicateFinder, libraryIdGenerator); // get the first one!

        // track optical duplicates using only those reads that are the first end...
        if (this.toMarkQueue.shouldBeInLocations(next) && next.getRecord().getFirstOfPairFlag()) {
            final Set<ReadEnds> locations = this.toMarkQueue.getLocations(next);

            if (!locations.isEmpty()) {
                AbstractMarkDuplicateFindingAlgorithm.trackOpticalDuplicates(new ArrayList<ReadEnds>(locations),
                        this.opticalDuplicateFinder, this.libraryIdGenerator);
            }
        }
    }

    /** Get the duplication metrics for the library associated with end. */
    private DuplicationMetrics getMetrics(final SAMRecord record) {
        final String library = this.libraryIdGenerator.getLibraryName(this.header, record);
        DuplicationMetrics metrics = this.libraryIdGenerator.getMetricsByLibrary(library);
        if (metrics == null) {
            metrics = new DuplicationMetrics();
            metrics.LIBRARY = library;
            this.libraryIdGenerator.addMetricsByLibrary(library, metrics);
        }
        return metrics;
    }

}
