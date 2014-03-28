package picard.sam.markduplicates;

import picard.cmdline.CommandLineProgram;
import picard.sam.markduplicates.MarkDuplicatesWithMateCigar;

public class MarkDuplicatesWithMateCigarTester extends AbstractMarkDuplicateFindingAlgorithmTester {

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicatesWithMateCigar(); }
}
