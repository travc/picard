package picard.sam.markduplicates.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import picard.PicardException;

import java.util.List;

/**
 * @author nhomer
 */
// NB: could refactor to use ReadEndsMC.java
public class ReadEndsMC extends ReadEnds {
    // to see if we either end is unmapped
    byte hasUnmapped = 0;

    // we need this reference so we can access the mate cigar among other things
    public SamRecordIndex samRecordIndex = null;

    /** Builds a read ends object that represents a single read. */
    public ReadEndsMC(final SAMFileHeader header, final SamRecordIndex samRecordIndex, final OpticalDuplicateFinder opticalDuplicateFinder, final short libraryId) {
        this.readGroup = -1;
        this.tile = -1;
        this.x = this.y = -1;
        this.read2Sequence = this.read2Coordinate = -1;
        this.hasUnmapped = 0;

        this.samRecordIndex = samRecordIndex;

        this.read1Sequence = this.samRecordIndex.getRecord().getReferenceIndex();
        this.read1Coordinate = this.samRecordIndex.getRecord().getReadNegativeStrandFlag() ? this.samRecordIndex.getRecord().getUnclippedEnd() : this.samRecordIndex.getRecord().getUnclippedStart();
        if (this.samRecordIndex.getRecord().getReadUnmappedFlag()) {
            throw new PicardException("Found an unexpected unmapped read");
        }

        if (this.samRecordIndex.getRecord().getReadPairedFlag() && !this.samRecordIndex.getRecord().getReadUnmappedFlag() && !this.samRecordIndex.getRecord().getMateUnmappedFlag()) {
            this.read2Sequence = this.samRecordIndex.getRecord().getMateReferenceIndex();
            this.read2Coordinate = this.samRecordIndex.getRecord().getMateNegativeStrandFlag() ? SAMUtils.getMateUnclippedEnd(this.samRecordIndex.getRecord()) : SAMUtils.getMateUnclippedStart(this.samRecordIndex.getRecord());

            // set orientation
            this.orientation = ReadEnds.getOrientationByte(this.samRecordIndex.getRecord().getReadNegativeStrandFlag(), this.samRecordIndex.getRecord().getMateNegativeStrandFlag());
        }
        else {
            this.orientation = this.samRecordIndex.getRecord().getReadNegativeStrandFlag() ? ReadEndsMC.R : ReadEndsMC.F;
        }

        // Fill in the library ID
        this.libraryId = libraryId;

        // Is this unmapped or its mate?
        if (this.samRecordIndex.getRecord().getReadUnmappedFlag() || (this.samRecordIndex.getRecord().getReadPairedFlag() && this.samRecordIndex.getRecord().getMateUnmappedFlag())) {
            this.hasUnmapped = 1;
        }

        // Fill in the location information for optical duplicates
        if (opticalDuplicateFinder.addLocationInformation(this.samRecordIndex.getRecord().getReadName(), this)) {
            // calculate the RG number (nth in list)
            // NB: could this be faster if we used a hash?
            this.readGroup = 0;
            final String rg = (String) this.samRecordIndex.getRecord().getAttribute("RG");
            final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
            if (rg != null && readGroups != null) {
                for (final SAMReadGroupRecord readGroup : readGroups) {
                    if (readGroup.getReadGroupId().equals(rg)) break;
                    else this.readGroup++;
                }
            }
        }
    }

    public ReadEndsMC(final ReadEndsMC other, final SamRecordIndex samRecordIndex) {
        this.readGroup = other.readGroup;
        this.tile = other.tile;
        this.x = other.x;
        this.y = other.y;
        this.read1Sequence = other.read1Sequence;
        this.read1Coordinate = other.read1Coordinate;
        this.read2Sequence = other.read2Sequence;
        this.read2Coordinate = other.read2Coordinate;
        this.hasUnmapped = other.hasUnmapped;
        this.samRecordIndex = samRecordIndex;
        this.orientation = other.orientation;
        this.libraryId = other.libraryId;
    }

    // Track that underlying SAMRecord has been through duplicate marking
//        public void setHasBeenMarkedFlag() { super.setHasBeenMarkedFlag(samRecordIndex.getReadName(), true);}
    public SamRecordIndex getsamRecordIndex() { return this.samRecordIndex; }
    public SAMRecord getRecord() { return this.samRecordIndex.getRecord(); }
    public String getRecordReadName() { return this.samRecordIndex.getRecord().getReadName(); }
    public int getRecordIndex() { return this.samRecordIndex.getCoordinateSortedIndex(); }

    @Override
    public boolean isPaired() { return this.getRecord().getReadPairedFlag(); }
}
