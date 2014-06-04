package picard.sam.markduplicates;

import java.io.File;

import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Created by bradt on 6/4/14.
 */
public class DuplicateMarkingBufferTest {
    // Create a separate directory for files so it is possible to confirm that the directory is emptied
    protected final File tmpDir = new File(System.getProperty("java.io.tmpdir") + "/" + System.getProperty("user.name"),
            "DuplicateMarkingBufferTest");
    @BeforeTest
    void setup() {
        // Clear out any existing files if the directory exists
        if (tmpDir.exists()) {
            for (final File f : tmpDir.listFiles()) {
                f.delete();
            }
        }
        tmpDir.mkdirs();
    }

    @AfterTest
    void tearDown() {
        System.err.println("In SortingCollectionTest.tearDown.  tmpDir: " + tmpDir);
        for (final File f : tmpDir.listFiles()) {
            f.delete();
        }
        tmpDir.delete();
    }

    private void fillRecordSet(final SAMRecordSetBuilder samRecordSetBuilder) {
        for (int iii = 0; iii < 100; iii++) {
            // Add a fairly generic fragment (positive strand, mapped, no cigar, no qual string, default base q, primary alignment)
            samRecordSetBuilder.addFrag("record"+iii, 1, 0, false, false, null, null,
                                        AbstractMarkDuplicateFindingAlgorithmTest.DEFAULT_BASE_QUALITY, false);
        }
    }

    private MarkDuplicatesWithMateCigarIterator.DuplicateMarkingBuffer fillDuplicateMarkingBuffer(final Iterator<SAMRecord> recordIterator,
                                                              final int maxRecordsInMemory,
                                                              final int blockSize,
                                                              final List<File> tmpDirs) {
        final MarkDuplicatesWithMateCigarIterator blorpth = new MarkDuplicatesWithMateCigarIterator()
    }

    @DataProvider (name = "duplicateMarkingBufferProvider")
    public Object[][] createTestData() {
        final SAMRecordSetBuilder samRecords = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.coordinate);
        fillRecordSet(samRecords);
        return new Object[][] {
                {"all in memory", samRecords.iterator(), 25, 1000},
                {"all on disk", samRecords.iterator(), 25, 0}
        };
    }

    @Test(dataProvider = "duplicateMarkingBufferProvider")
    public void testPositive(final String testName,
                             final Iterator<SAMRecord> recordIterator,
                             final int blockSize,
                             final int maxRecordsInMemory) {
        while (recordIterator.hasNext()) {
            Assert.assertTrue(recordIterator.next().getAlignmentStart() == 0);
        }
    }
}
