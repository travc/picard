/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.sam.markduplicates;

import picard.PicardException;
import org.testng.annotations.Test;

public abstract class AbstractMarkDuplicateFindingAlgorithmTest {

    protected abstract AbstractMarkDuplicateFindingAlgorithmTester getTester();

    protected final static int DEFAULT_BASE_QUALITY = 10;

    //TODO - Add a test case that includes reads from multiple chromosomes, with some expected duplicates from each chromosome

    @Test
    public void testSingleUnmappedFragment() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addUnmappedFragment(-1, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testSingleUnmappedPair() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addUnmappedPair(-1, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }


    @Test
    public void testSingleMappedFragment() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedFragment(1, 1, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoMappedFragments() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedFragment(0, 1, false, DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(0, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testSingleMappedPair() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testSingleMappedFragmentAndSingleMappedPair() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedFragment(1, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairs() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testThreeMappedPairs() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testSingleMappedFragmentAndTwoMappedPairs() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedFragment(1, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsMatesSoftClipped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 10022, 10051, false, false, "76M", "8S68M", false, true, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 10022, 10063, false, false, "76M", "5S71M", false, true, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithSoftClipping() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        // NB: no duplicates
        // 5'1: 2, 5'2:46+73M=118
        // 5'1: 2, 5'2:51+68M=118
        tester.addMappedPair(1, 2, 46, false, false, "6S42M28S", "3S73M", false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 2, 51, true, true, "6S42M28S", "8S68M", false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithSoftClippingFirstOfPairOnlyNoMateCigar() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.setNoMateCigars(true);
        // NB: no duplicates
        // 5'1: 2, 5'2:46+73M=118
        // 5'1: 2, 5'2:51+68M=118
        tester.addMappedPair(1, 12, 46, false, false, "6S42M28S", null, true, DEFAULT_BASE_QUALITY); // only add the first one
        tester.addMappedPair(1, 12, 51, false, false, "6S42M28S", null, true, DEFAULT_BASE_QUALITY); // only add the first one
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithSoftClippingBoth() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 10046, 10002, false, false, "3S73M", "6S42M28S", true, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 10051, 10002, true, true, "8S68M", "6S48M22S", true, false, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testMatePairFirstUnmapped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMatePair(1, 10049, 10049, false, true, false, false, "11M2I63M", null, false, false, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testMatePairSecondUnmapped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMatePair(1, 10056, 10056, true, false, false, false, null, "54M22S", false, false, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }


    @Test
    public void testMappedFragmentAndMatePairOneUnmapped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMatePair(1, 10049, 10049, false, true, false, false, "11M2I63M", null, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(1, 10049, true, DEFAULT_BASE_QUALITY); // duplicate
        tester.runTest();
    }

    @Test
    public void testMappedPairAndMatePairOneUnmapped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMatePair(1, 10040, 10040, false, true, true, false, "76M", null, false, false, false, DEFAULT_BASE_QUALITY); // first a duplicate,
        // second end unmapped
        tester.addMappedPair(1, 10189, 10040, false, false, "41S35M", "65M11S", true, false, false, DEFAULT_BASE_QUALITY); // mapped OK
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithOppositeOrientations() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 10182, 10038, false, false, "32S44M", "66M10S", true, false, false, DEFAULT_BASE_QUALITY); // -/+
        tester.addMappedPair(1, 10038, 10182, true, true, "70M6S", "32S44M", false, true, false, DEFAULT_BASE_QUALITY); // +/-, both are duplicates
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithOppositeOrientationsNumberTwo() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 10038, 10182, false, false, "70M6S", "32S44M", false, true, false, DEFAULT_BASE_QUALITY); // +/-, both are duplicates
        tester.addMappedPair(1, 10182, 10038, true, true, "32S44M", "66M10S", true, false, false, DEFAULT_BASE_QUALITY); // -/+
        tester.runTest();
    }

    @Test
    public void testThreeMappedPairsWithMatchingSecondMate() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        // Read0 and Read2 are duplicates
        // 10181+35=10216, 10058
        tester.addMappedPair(1, 10181, 10058, false, false, "41S35M", "47M29S", true, false, false, DEFAULT_BASE_QUALITY); // -/+
        // 10181+37=10218, 10058
        tester.addMappedPair(1, 10181, 10058, false, false, "37S39M", "44M32S", true, false, false, DEFAULT_BASE_QUALITY); // -/+
        // 10180+36=10216, 10058
        tester.addMappedPair(1, 10180, 10058, true, true, "36S40M", "50M26S", true, false, false, DEFAULT_BASE_QUALITY); // -/+, both are duplicates
        tester.runTest();
    }

    @Test
    public void testMappedPairWithSamePosition() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 4914, 4914, false, false, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairWithSamePosition() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(0, 5604914, 5604914, false, false, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.addMappedPair(0, 5604914, 5604914, true, true, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.runTest();
    }

    @Test
    public void testMappedPairWithFirstEndSamePositionAndOther() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(0, 5604914, 5605914, false, false, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.addMappedPair(0, 5604914, 5604914, false, false, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.runTest();
    }

    @Test
    public void testTwoGroupsOnDifferentChromosomesOfTwoFragments() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedFragment(0, 1, false, DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(0, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedFragment(1, 1, false, DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(1, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }


    @Test
    public void testTwoGroupsOnDifferentChromosomesOfTwoMappedPairs() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(0, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testTwoGroupsOnDifferentChromosomesOfThreeMappedPairs() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(0, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testThreeGroupsOnDifferentChromosomesOfThreeMappedPairs() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(0, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(2, 1, 100, false, false , DEFAULT_BASE_QUALITY);
        tester.addMappedPair(2, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(2, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testBulkFragmentsNoDuplicates() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        for(int position = 1; position <= 10000; position += 1) {
            tester.addMappedFragment(0, position, false, "100M", DEFAULT_BASE_QUALITY);
        }
        tester.runTest();
    }

    @Test
    public void testBulkFragmentsWithDuplicates() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        for(int position = 1; position <= 10000; position += 1) {
            tester.addMappedFragment(0, position, false, "100M", DEFAULT_BASE_QUALITY);
            tester.addMappedFragment(0, position, true, "100M", DEFAULT_BASE_QUALITY);
            tester.addMappedFragment(0, position, true, "100M", DEFAULT_BASE_QUALITY);
            tester.addMappedFragment(0, position, true, "100M", DEFAULT_BASE_QUALITY);
            tester.addMappedFragment(0, position, true, "100M", DEFAULT_BASE_QUALITY);
        }
        tester.runTest();
    }
}
