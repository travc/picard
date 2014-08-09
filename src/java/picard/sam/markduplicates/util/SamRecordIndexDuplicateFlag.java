package picard.sam.markduplicates.util;

import htsjdk.samtools.SAMRecord;

/**
 * @author nhomer
 */
public class SamRecordIndexDuplicateFlag extends SamRecordIndex {

    public SamRecordIndexDuplicateFlag() {
        super();
    }

    public SamRecordIndexDuplicateFlag(final SAMRecord record, final int recordIndex) {
        super(record, recordIndex);
    }

    @Override
    public void setExaminedState(final boolean examinedState) {
        this.getRecord().setDuplicateReadFlag(examinedState);
    }
}
