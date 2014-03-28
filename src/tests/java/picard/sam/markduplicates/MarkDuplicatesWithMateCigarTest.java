package picard.sam.markduplicates;

import picard.PicardException;
import org.testng.annotations.Test;
/**
 * @author nhomer@broadinstitute.org
 */
public class MarkDuplicatesWithMateCigarTest extends AbstractMarkDuplicateFindingAlgorithmTest {
    protected AbstractMarkDuplicateFindingAlgorithmTester getTester() {
        return new MarkDuplicatesWithMateCigarTester();
    }

    // TODO: test program record chaining, including failures. Use MarkDuplicate's facility.
    // TODO: check if one mate is dup, the other is as well, only if both are mapped

    // NB: this test should return different results than MarkDuplicatesWithMateCigar, as we have the mate cigar
    @Test
    public void testTwoMappedPairsWithSoftClippingFirstOfPairOnly() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        // NB: no duplicates
        // 5'1: 2, 5'2:46+73M=118
        // 5'1: 2, 5'2:51+68M=118
        tester.addMappedPair(0, 12, 46, false, false, "6S42M28S", "3S73M", true, 50); // only add the first one
        // NB: this next record should not be a duplicate in MarkDuplicates, but is here, because have the mate cigar
        tester.addMappedPair(0, 12, 51, true, true, "6S42M28S", "8S68M", true, 50); // only add the first one
        tester.runTest();
    }

    @Test
    public void testTwoFragmentsLargeSoftClipWithMinimumDistanceOK() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addArg("MINIMUM_DISTANCE=990");
        tester.addMappedFragment(0, 1000, false, "100M", DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(0, 2000, false, "10S100M", DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(0, 3000, true, "2000S100M", DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test(expectedExceptions = PicardException.class)
    public void testTwoFragmentsLargeSoftClipWithMinimumDistanceFailure() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addArg("MINIMUM_DISTANCE=989");
        tester.addMappedFragment(0, 1000, false, "100M", DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(0, 2000, false, "10S100M", DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(0, 3000, true, "2000S100M", DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test(expectedExceptions = PicardException.class)
    public void testTwoFragmentsLargeSoftClip() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedFragment(0, 1000, false, "100M", DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(0, 2000, false, "10S100M", DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(0, 3000, true, "2000S100M", DEFAULT_BASE_QUALITY);
        tester.runTest();
    }
}
