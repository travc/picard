package picard.sam.markduplicates;

import picard.cmdline.CommandLineProgram;
import picard.sam.markduplicates.MarkDuplicates;

public class MarkDuplicatesTester extends AbstractMarkDuplicateFindingAlgorithmTester {

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicates(); }
}
