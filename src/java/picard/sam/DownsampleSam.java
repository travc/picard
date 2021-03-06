package picard.sam;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.Usage;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 * Class to randomly downsample a BAM file while respecting that we should either get rid
 * of both ends of a pair or neither end of the pair!
 */
public class DownsampleSam extends CommandLineProgram {
    @Usage public final String USAGE = getStandardUsagePreamble() + " Randomly down-sample a SAM or BAM file to retain " +
            "a random subset of the reads. Mate-pairs are either both kept or both discarded. Reads marked as not primary " +
            "alignments are all discarded. Each read is given a probability P of being retained - results with the exact " +
            "same input in the same order and with the same value for RANDOM_SEED will produce the same results.";

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input SAM or BAM file to downsample.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output, downsampled, SAM or BAM file to write.")
    public File OUTPUT;

    @Option(shortName="R", doc="Random seed to use if reproducibilty is desired.  " +
            "Setting to null will cause multiple invocations to produce different results.")
    public Long RANDOM_SEED = 1L;

    @Option(shortName="P", doc="The probability of keeping any individual read, between 0 and 1.")
    public double PROBABILITY = 1;

    private final Log log = Log.getInstance(DownsampleSam.class);

    public static void main(final String[] args) {
        new DownsampleSam().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final Random r = RANDOM_SEED == null ? new Random() : new Random(RANDOM_SEED);
        final SAMFileReader in = new SAMFileReader(INPUT);
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(), true, OUTPUT);
        final Map<String,Boolean> decisions = new HashMap<String,Boolean>();

        long total = 0;
        long kept  = 0;
        
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Read");

        for (final SAMRecord rec : in) {
            if (rec.isSecondaryOrSupplementary()) continue;
            ++total;

            final String key = rec.getReadName();
            final Boolean previous = decisions.remove(key);
            final boolean keeper;

            if (previous == null) {
                keeper = r.nextDouble() <= PROBABILITY;
                if (rec.getReadPairedFlag()) decisions.put(key, keeper);
            }
            else {
                keeper = previous;
            }

            if (keeper) {
                out.addAlignment(rec);
                ++kept;
            }

            progress.record(rec);
        }

        out.close();
        log.info("Finished! Kept " + kept + " out of " + total + " reads.");

        return 0;
    }
}
