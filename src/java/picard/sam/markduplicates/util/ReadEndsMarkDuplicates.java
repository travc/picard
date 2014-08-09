package picard.sam.markduplicates.util;

/**
 *  Little struct-like class to hold read pair (and fragment) end data for MarkDuplicatesWithMateCigar
 * @author Nils Homer
 */
public class ReadEndsMarkDuplicates extends ReadEnds {
    public static final int SIZE_OF = (1*1) + (2*1) + (4*4) + (8*2) + 2 + 1 + 2 + 2
            + 8 + // last 8 == reference overhead
            13; // This is determined experimentally with JProfiler

    public short score = 0;
    public long read1IndexInFile = -1;
    public long read2IndexInFile = -1;
    public String name = "";
}
