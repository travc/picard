package picard.sam.markduplicates;


import picard.cmdline.CommandLineProgram;
import picard.sam.testers.SamFileTester;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.TestUtil;
import org.testng.Assert;

import java.io.File;

/**
 * This class is an extension of SamFileTester used to test MarkDuplicates with SAM files generated on the fly.
 */
abstract public class AbstractMarkDuplicateFindingAlgorithmTester extends SamFileTester {

    public AbstractMarkDuplicateFindingAlgorithmTester() {
        super(50, true);

        final File metrics = new File(getOutputDir(), "metrics.txt");
        addArg("METRICS_FILE=" + metrics);
    }

    @Override
    public void test() {
        try {
            // Read the output and check the duplicate flag
            final SAMFileReader reader = new SAMFileReader(getOutput());
            for (final SAMRecord record : reader) {
                final String key = samRecordToDuplicatesFlagsKey(record);
                if (!this.duplicateFlags.containsKey(key)) {
                    System.err.println("DOES NOT CONTAIN KEY: " + key);
                }
                Assert.assertTrue(this.duplicateFlags.containsKey(key));
                final boolean value = this.duplicateFlags.get(key);
                this.duplicateFlags.remove(key);
                if (value != record.getDuplicateReadFlag()) {
                    System.err.println("Mismatching read:");
                    System.err.print(record.getSAMString());
                }
                Assert.assertEquals(record.getDuplicateReadFlag(), value);
            }
            reader.close();
        } finally {
            TestUtil.recursiveDelete(getOutputDir());
        }
    }

    @Override
    abstract protected CommandLineProgram getProgram();
}

