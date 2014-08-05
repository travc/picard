package picard.sam.markduplicates.util;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.DiskBackedQueue;
import picard.PicardException;

import java.io.File;
import java.util.ArrayDeque;
import java.util.Deque;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * @author nhomer
 */
public class DuplicateMarkingBuffer {
    private int availableRecordsInMemory;
    private final int blockSize;
    private final List<File> tmpDirs;
    private int queueHeadRecordIndex;
    private int queueTailRecordIndex;
    private final Deque<BufferBlock> blocks;
    private final SAMFileHeader header;

    public DuplicateMarkingBuffer(final int maxRecordsInMemory, final int blockSize, final List<File> tmpDirs, final SAMFileHeader header) {
        this.availableRecordsInMemory = maxRecordsInMemory;
        this.blockSize = blockSize;
        this.tmpDirs = tmpDirs;
        this.queueHeadRecordIndex = -1;
        this.queueTailRecordIndex = -1;
        this.blocks = new ArrayDeque<BufferBlock>();
//            System.out.println("creating the buffer");
//            System.out.println("Max block size: " + this.blockSize);
//            System.out.println("max records in ram: " + this.availableRecordsInMemory);
        this.header = header;
    }

    public boolean isEmpty() { return (blocks.size() == 0  || this.blocks.getFirst().isEmpty()); }

    /**
     * Returns true if the head samRecordIndex in the DuplicateMarkingBuffer is annotated as having been through duplicate marking
     *
     * @return true if the head samRecordIndex in the buffer has been through duplicate marking
     */
    public boolean canEmit() { return (this.blocks.size() !=0 && this.blocks.getFirst().canEmit()); }

    /**
     * Add the provided samRecordIndex to the tail of this DuplicateMarkingBuffer
     * @param samRecordIndex The samRecordIndex to be added
     */
    public void add(final SamRecordIndex samRecordIndex) {
        if (this.isEmpty()) {
            this.queueHeadRecordIndex = samRecordIndex.getCoordinateSortedIndex();
            this.queueTailRecordIndex = samRecordIndex.getCoordinateSortedIndex() - 1;
        }
        this.queueTailRecordIndex++;
        //TODO - OOOOOOOOOK. block OriginalRecordIndex getting set whenever queue is empty, which means that the wrong byte arrays are bing accessed and that records are being flushed prematurely. Only want to set OriginalRecordIndex WHEN BLOCK IS CREATED
        // If necessary, create a new block, using as much ram as available up to its total size
        if (this.blocks.size() == 0 || !this.blocks.getLast().canAdd()) {
            // once ram is given to a block, we can't give it to another block (until some is recovered from the head of the queue)
            final int blockRam = Math.min(this.blockSize, this.availableRecordsInMemory);
            this.availableRecordsInMemory = this.availableRecordsInMemory - blockRam;
            final BufferBlock block = new BufferBlock(this.blockSize, blockRam, this.tmpDirs, this.header);
            block.setOriginalRecordIndex(samRecordIndex.getCoordinateSortedIndex());
            this.blocks.addLast(block);
//                System.out.println("\nNumber of blocks is " + this.blocks.size());
//                System.out.println("size of queue on addition: " + this.size());
        }
        this.blocks.getLast().add(samRecordIndex);  // TODO- need to catch BufferBlock and DiskBackedQueue exceptions here?
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
            throw new NoSuchElementException("Attempting to remove an element from an empty DuplicateMarkingBuffer");
        final BufferBlock headBlock = this.blocks.getFirst();
        if (!headBlock.canEmit())
            throw new PicardException("Attempting to get a samRecordIndex from the DuplicateMarkingBuffer that has not been through " +
                    "duplicate marking. canEmit() must return true in order to call next()");

        // If the samRecordIndex was stored in memory, reclaim its ram for use in additional blocks at tail of queue
        // NB: this must be checked before calling next(), as that method updates the block-head
        if (!headBlock.headRecordIsFromDisk())
            this.availableRecordsInMemory++;
        final SamRecordIndex samRecordIndex = headBlock.next();
        if (headBlock.hasBeenDrained()) { // TODO- need to catch BufferBlock and DiskBackedQueue exceptions here?
            blocks.poll(); // remove the block as it is now empty
            headBlock.clear(); // free any disk io resources associated with empty block
//                System.out.println("removing block! New num blocks: " + this.blocks.size());
//                System.out.println("Queue size at removal time: " + this.size());
        }
        this.queueHeadRecordIndex++;
        return samRecordIndex;
    }

    public void remove() { this.next(); }

    /**
     * Return the total number of elements in the queue, both in memory and on disk
     */
    public int size() { return this.queueTailRecordIndex - this.queueHeadRecordIndex + 1; }

    private BufferBlock getBlock(final SamRecordIndex samRecordIndex) {
        for (final BufferBlock block : this.blocks) {
            if (block.getStartIndex() <= samRecordIndex.getCoordinateSortedIndex() && block.getEndIndex() >= samRecordIndex.getCoordinateSortedIndex()) {
                return block;
            }
        }
        return null;
    }

    /**
     * TODO: document
     */
    public boolean contains(final SamRecordIndex samRecordIndex) {
        return (null != getBlock(samRecordIndex));
    }

    /**
     * Mark the current samRecordIndex as having been through duplicate-marking, and whether it is a duplicate
     *
     * @param samRecordIndex The samRecordIndex to be marked
     * @param isDuplicate Boolean flag indicating whether this is a duplicate samRecordIndex
     * @throws PicardException if the provided recordIndex is not found within the DuplicateMarkingBuffer
     */
    public void setDuplicateMarkingFlags(final SamRecordIndex samRecordIndex, final boolean isDuplicate) {
        BufferBlock block = getBlock(samRecordIndex);
        if (null == block) {
            throw new PicardException("Attempted to set duplicate-marking information on a samRecordIndex whose index is not found " +
                    "in the DuplicateMarkingBuffer. recordIndex: " + samRecordIndex.getCoordinateSortedIndex());
        }
        block.setDuplicateMarkingIndexes(samRecordIndex, isDuplicate);
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
     * TODO - document this class
     */
    private class BufferBlock {
        private final DiskBackedQueue<SAMRecord> recordsQueue;
        private final int maxBlockSize;
        private int currentStartIndex;
        private int originalStartIndex;
        private int endIndex;
        private byte[] wasCheckedIndexes = null;
        private byte[] isDuplicateIndexes = null;

        // TODO- Method javadocs

        public BufferBlock(final int maxBlockSize, final int maxBlockRecordsInMemory, final List<File> tmpDirs, final SAMFileHeader header) {
            this.recordsQueue = DiskBackedQueue.newInstance(new BAMRecordCodec(header), maxBlockRecordsInMemory, tmpDirs);
            this.maxBlockSize = maxBlockSize;
            this.currentStartIndex = 0;
            this.originalStartIndex = 0;
            this.endIndex = -1;
            this.wasCheckedIndexes = new byte[maxBlockSize];
            this.isDuplicateIndexes = new byte[maxBlockSize];

//                System.out.println("creating a new buffer block!");
        }

        /**
         * Check that the tail of the block has not grown past the maximum block size (even if records were popped) and that the underlying queue can be added to.
         * TODO - reimplement with a circular byte array buffer PROVIDED RECORDS ARE IN MEMORY
         * @return
         */
        public boolean canAdd() { return (this.endIndex - this.originalStartIndex + 1) < this.maxBlockSize && this.recordsQueue.canAdd(); }

        public boolean headRecordIsFromDisk() { return this.recordsQueue.headRecordIsFromDisk(); }

        /**
         * Check whether we have read all possible records from this block (and it is available to be destroyed)
         * @return true if we have read the last /possible/ record (ie the block size, or if !canAdd the end index)
         */
        public boolean hasBeenDrained() {
            final int maximalIndex = (this.canAdd()) ? (this.originalStartIndex + this.maxBlockSize) : this.endIndex;
            return this.currentStartIndex > maximalIndex;       //TODO- watch out for an off by one here
        }

        public int getStartIndex() { return this.currentStartIndex; }

        public int getEndIndex() { return this.endIndex; }

        public void add(final SamRecordIndex samRecordIndex) {
            if (this.recordsQueue.canAdd()) {
                if (this.recordsQueue.isEmpty()) {
                    this.currentStartIndex = samRecordIndex.getCoordinateSortedIndex();
//                        this.originalStartIndex = samRecordIndex.getCoordinateSortedIndex();
                    this.endIndex = samRecordIndex.getCoordinateSortedIndex() - 1;
//                        System.out.println("This block's queue was empty. Adding first record, index: " + samRecordIndex.getCoordinateSortedIndex());
                }
                this.recordsQueue.add(samRecordIndex.getRecord());
                this.endIndex++;
            } else {
                throw new IllegalStateException("Cannot add to DiskBackedQueue whose canAdd() method returns false");
            }
        }

        public void setDuplicateMarkingIndexes(final SamRecordIndex samRecordIndex, final boolean isDuplicate) {
            // find the correct byte array index and update both metadata byte arrays
            this.wasCheckedIndexes[samRecordIndex.getCoordinateSortedIndex() - this.originalStartIndex] = 1;
            this.isDuplicateIndexes[samRecordIndex.getCoordinateSortedIndex() - this.originalStartIndex] = (isDuplicate) ? (byte)1 : 0; //NB: why the need to cast here?
        }


        public boolean isEmpty() {
            return (this.recordsQueue.isEmpty());
        }

        public boolean canEmit() {
            return (this.wasCheckedIndexes[this.currentStartIndex - this.originalStartIndex] == 1);
        }

        public SamRecordIndex next() throws IllegalStateException {
            if (this.canEmit()) {
                // create a wrapped record for the head of the queue, and set the underlying record's dup information appropriately
                final SamRecordIndex samRecordIndex = new SamRecordIndex(this.recordsQueue.poll(), this.currentStartIndex);
                samRecordIndex.getRecord().setDuplicateReadFlag(this.isDuplicateIndexes[this.currentStartIndex - this.originalStartIndex] == 1);
                this.currentStartIndex++;
                return samRecordIndex;
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

        public void setOriginalRecordIndex(final int recordIndex) {
            this.originalStartIndex = recordIndex;
        }
    }
}
