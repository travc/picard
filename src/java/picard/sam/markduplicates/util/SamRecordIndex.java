package picard.sam.markduplicates.util;

import htsjdk.samtools.SAMRecord;

/**
 * @author nhomer
 */
public class SamRecordIndex {
    private final SAMRecord record;
    private final int coordinateSortedIndex;

    public SamRecordIndex(final SAMRecord record, final int recordIndex) {
        this.record = record;
        this.coordinateSortedIndex = recordIndex;
    }

    public SAMRecord getRecord() { return this.record; }
    public int getCoordinateSortedIndex() { return this.coordinateSortedIndex; }
}
