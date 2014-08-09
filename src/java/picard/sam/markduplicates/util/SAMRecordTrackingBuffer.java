package picard.sam.markduplicates.util;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.DiskBackedQueue;
import picard.PicardException;

import java.io.File;
import java.util.ArrayDeque;
import java.util.BitSet;
import java.util.Deque;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * This class stores SAMRecords for return.  The purpose of this class is to buffer records that need to be modified or processed in some
 * fashion, and only return (or emit) them when they have been recorded as being fully examined.  If we have too many records in RAM,
 * we can spill over to disk.  The order in which they are given (via SamRecordIndex) determines their order of being returned.  It is the
 * responsibility of the user of this class to make sure all records have unique index and are added in order.
 *
 * When a record is examined, we also store a state for that examination.  This is currently a boolean to reduce on memory and disk footprint.
 *
 * We store groups of records in blocks and the size of these blocks can be controlled.  If we have too many records in RAM, we start
 * spilling blocks to disk.
 *
 * Users should check isEmpty() to see if any records are still being tracked.  If so, they should check canEmit() to see if the
 * next record can be returned.  If so, they can call next() to get that record.
 *
 * When users are done with this structure, call close().
 *
 * @author bradtaylor
 */
public class SAMRecordTrackingBuffer {
    private int availableRecordsInMemory; // how many more records can we store in memory
    private final int blockSize; // the size of each block
    private final List<File> tmpDirs; // the list of temporary directories to use
    private int queueHeadRecordIndex; // the index of the head of the buffer
    private int queueTailRecordIndex; // the index of the tail of the buffer
    private final Deque<BufferBlock> blocks; // the queue of blocks, in which records are contained
    private final SAMFileHeader header;

    private final Class<? extends SamRecordIndex> clazz; // the class to create

    /**
     * @param maxRecordsInRam how many records to buffer before spilling to disk
     * @param blockSize the number of records in a given block
     * @param tmpDirs the temporary directories to use when spilling to disk
     * @param header the header
     * @param clazz the class that extends SamRecordIndex
     */
    public SAMRecordTrackingBuffer(final int maxRecordsInRam, final int blockSize, final List<File> tmpDirs, final SAMFileHeader header, final Class<? extends SamRecordIndex> clazz) {
        this.availableRecordsInMemory = maxRecordsInRam;
        this.blockSize = blockSize;
        this.tmpDirs = tmpDirs;
        this.queueHeadRecordIndex = -1;
        this.queueTailRecordIndex = -1;
        this.blocks = new ArrayDeque<BufferBlock>();
        this.header = header;
        this.clazz = clazz;
    }

    /** Returns true if we are tracking no records, false otherwise */
    public boolean isEmpty() { return (blocks.size() == 0  || this.blocks.getFirst().isEmpty()); }

    /** Returns true if we can return the next record (it has been examined). */
    public boolean canEmit() { return (this.blocks.size() !=0 && this.blocks.getFirst().canEmit()); }

    /**
     * Add the given SAMRecordIndex to the buffer.  The records must be added in order.
     * @param samRecordIndex The samRecordIndex to be added
     */
    public void add(final SamRecordIndex samRecordIndex) {
        if (this.isEmpty()) {
            this.queueHeadRecordIndex = samRecordIndex.getRecordIndex();
            this.queueTailRecordIndex = samRecordIndex.getRecordIndex() - 1;
        }
        this.queueTailRecordIndex++;
        if (samRecordIndex.getRecordIndex() != this.queueTailRecordIndex) {
            throw new PicardException("The records were added out of order");
        }
        // If necessary, create a new block, using as much ram as available up to its total size
        if (this.blocks.size() == 0 || !this.blocks.getLast().canAdd()) {
            // once ram is given to a block, we can't give it to another block (until some is recovered from the head of the queue)
            final int blockRam = Math.min(this.blockSize, this.availableRecordsInMemory);
            this.availableRecordsInMemory = this.availableRecordsInMemory - blockRam;
            final BufferBlock block = new BufferBlock(this.blockSize, blockRam, this.tmpDirs, this.header, samRecordIndex.getRecordIndex());
            this.blocks.addLast(block);
        }
        this.blocks.getLast().add(samRecordIndex);
    }

    /**
     * Returns the next element in the iteration.
     *
     * @return The next element in the iteration.
     * @throws java.util.NoSuchElementException if the buffer is empty.
     * @throws picard.PicardException if the buffer is not competent to emit (canEmit returns false)
     */
    public SamRecordIndex next() {
        if (this.isEmpty())
            throw new NoSuchElementException("Attempting to remove an element from an empty SAMRecordTrackingBuffer");
        final BufferBlock headBlock = this.blocks.getFirst();
        if (!headBlock.canEmit())
            throw new PicardException("Attempting to get a samRecordIndex from the SAMRecordTrackingBuffer that has not been through " +
                    "marked as examined. canEmit() must return true in order to call next()");

        // If the samRecordIndex was stored in memory, reclaim its ram for use in additional blocks at tail of queue
        // NB: this must be checked before calling next(), as that method updates the block-head
        if (!headBlock.headRecordIsFromDisk()) {
            this.availableRecordsInMemory++;
        }
        final SamRecordIndex samRecordIndex = headBlock.next();
        if (headBlock.hasBeenDrained()) {
            blocks.poll(); // remove the block as it is now empty
            headBlock.clear(); // free any disk io resources associated with empty block
        }
        this.queueHeadRecordIndex++;
        return samRecordIndex;
    }

    /** Removes the next record from this buffer */
    public void remove() { this.next(); }

    /**
     * Return the total number of elements in the queue, both in memory and on disk
     */
    public int size() { return this.queueTailRecordIndex - this.queueHeadRecordIndex + 1; }

    /** Returns the block that holds the sam record at the given index, null if no such block exists */
    private BufferBlock getBlock(final SamRecordIndex samRecordIndex) {
        for (final BufferBlock block : this.blocks) {
            if (block.getStartIndex() <= samRecordIndex.getRecordIndex() && block.getEndIndex() >= samRecordIndex.getRecordIndex()) {
                return block;
            }
        }
        return null;
    }

    /** Returns true if this buffer contains the record at the given index, false otherwise */
    public boolean contains(final SamRecordIndex samRecordIndex) {
        return (null != getBlock(samRecordIndex));
    }

    /**
     * Mark the current samRecordIndex as having been examined.
     *
     * @param samRecordIndex The samRecordIndex to be marked
     * @param examinedState Boolean flag indicating whether the state of the examination.
     * @throws PicardException if the provided recordIndex is not found within the SAMRecordTrackingBuffer
     */
    public void setExamined(final SamRecordIndex samRecordIndex, final boolean examinedState) {
        final BufferBlock block = getBlock(samRecordIndex);
        if (null == block) {
            throw new PicardException("Attempted to set examined information on a samRecordIndex whose index is not found " +
                    "in the SAMRecordTrackingBuffer. recordIndex: " + samRecordIndex.getRecordIndex());
        }
        block.setExamined(samRecordIndex, examinedState);
    }

    /**
     * Close IO resources associated with each underlying BufferBlock
     */
    public void close() {
        while (!blocks.isEmpty()) {
            final BufferBlock block = blocks.pollFirst();
            block.clear();
        }
    }

    /**
     * This stores blocks of records, either in memory or on disk, or both!
     */
    private class BufferBlock {
        private final DiskBackedQueue<SAMRecord> recordsQueue;
        private final int maxBlockSize;
        private int currentStartIndex;
        private final int originalStartIndex;
        private int endIndex;

        private final BitSet wasExaminedIndexes;
        private final BitSet examinedStateIndexes;

        /** Creates an empty block buffer, with an allowable # of records in RAM */
        public BufferBlock(final int maxBlockSize, final int maxBlockRecordsInMemory, final List<File> tmpDirs,
                           final SAMFileHeader header,
                           final int originalStartIndex) {
            this.recordsQueue = DiskBackedQueue.newInstance(new BAMRecordCodec(header), maxBlockRecordsInMemory, tmpDirs);
            this.maxBlockSize = maxBlockSize;
            this.currentStartIndex = 0;
            this.endIndex = -1;
            this.wasExaminedIndexes = new BitSet(maxBlockSize);
            this.examinedStateIndexes = new BitSet(maxBlockSize);
            this.originalStartIndex = originalStartIndex;
        }

        /**
         * Check that the tail of the block has not grown past the maximum block size (even if records were popped) and that the underlying queue can be added to.
         * TODO - reimplement with a circular byte array buffer PROVIDED RECORDS ARE IN MEMORY
         * @return
         */
        public boolean canAdd() { return (this.endIndex - this.originalStartIndex + 1) < this.maxBlockSize && this.recordsQueue.canAdd(); }

        /** Returns true if the record at the front of the buffer is on disk */
        public boolean headRecordIsFromDisk() { return this.recordsQueue.headRecordIsFromDisk(); }

        /**
         * Check whether we have read all possible records from this block (and it is available to be destroyed)
         * @return true if we have read the last /possible/ record (ie the block size, or if !canAdd the end index)
         */
        public boolean hasBeenDrained() {
            final int maximalIndex = (this.canAdd()) ? (this.originalStartIndex + this.maxBlockSize) : this.endIndex;
            return this.currentStartIndex > maximalIndex;       //NB: watch out for an off by one here
        }

        /** Gets the index of the first record in this block */
        public int getStartIndex() { return this.currentStartIndex; }

        /** Gets the index of the last record in this block */
        public int getEndIndex() { return this.endIndex; }

        /** Add a record to this block */
        public void add(final SamRecordIndex samRecordIndex) {
            if (this.recordsQueue.canAdd()) {
                if (this.recordsQueue.isEmpty()) {
                    this.currentStartIndex = samRecordIndex.getRecordIndex();
                    this.endIndex = samRecordIndex.getRecordIndex() - 1;
                }
                this.recordsQueue.add(samRecordIndex.getRecord());
                this.endIndex++;
            } else {
                throw new IllegalStateException("Cannot add to DiskBackedQueue whose canAdd() method returns false");
            }
        }


        /**
         * Mark the current samRecordIndex as having been examined.
         *
         * @param samRecordIndex The samRecordIndex to be marked
         * @param examinedState Boolean flag indicating whether the state of the examination.
         *
         * This assumes that this record index does not fall out of range.
         */
        public void setExamined(final SamRecordIndex samRecordIndex, final boolean examinedState) {
            // find the correct byte array index and update both metadata byte arrays
            this.wasExaminedIndexes.set(samRecordIndex.getRecordIndex() - this.originalStartIndex, true);
            this.examinedStateIndexes.set(samRecordIndex.getRecordIndex() - this.originalStartIndex, examinedState);
        }

        public boolean isEmpty() {
            return (this.recordsQueue.isEmpty());
        }

        public boolean canEmit() {
            // TODO: what if isEmpty() == true?
            return this.wasExaminedIndexes.get(this.currentStartIndex - this.originalStartIndex);
        }

        public SamRecordIndex next() throws IllegalStateException {
            if (this.canEmit()) {
                try {
                    // create a wrapped record for the head of the queue, and set the underlying record's examined information appropriately
                    final SamRecordIndex samRecordIndex = clazz.newInstance();
                    samRecordIndex.setRecord(this.recordsQueue.poll());
                    samRecordIndex.setRecordIndex(this.currentStartIndex);
                    samRecordIndex.setExaminedState(this.examinedStateIndexes.get(this.currentStartIndex - this.originalStartIndex));
                    this.currentStartIndex++;
                    return samRecordIndex;
                } catch (final Exception e) {
                    throw new RuntimeException(e);
                }
            } else {
                throw new IllegalStateException("Cannot call next() on a buffer block where canEmit() is false!");
            }
        }

        /**
         * Remove, but do not return, the next samRecordIndex in the iterator
         */
        public void remove() { this.next(); }

        /**
         * Return the total number of elements in the block, both in memory and on disk
         */
        public int size() { return this.endIndex - this.currentStartIndex + 1; }

        /**
         * Close disk IO resources associated with the underlying records queue.
         * This must be called when a block is no longer needed in order to prevent memory leaks.
         */
        public void clear() { this.recordsQueue.clear(); }
    }
}
