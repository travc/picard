package picard.sam.markduplicates;

/**
 *  Little struct-like class to hold read pair (and fragment) end data for MarkDuplicates
 * @author Nils Homer
 */
class ReadEndsMarkDuplicates extends ReadEnds {
    public static final int SIZE_OF = (1*1) + (2*1) + (4*4) + (8*2) + 2 + 1 + 2 + 2
            + 8 + // last 8 == reference overhead
            13; // This is determined experimentally with JProfiler

    short score = 0;
    long read1IndexInFile = -1;
    long read2IndexInFile = -1;
}
