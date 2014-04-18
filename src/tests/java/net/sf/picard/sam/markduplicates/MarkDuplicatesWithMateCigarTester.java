package net.sf.picard.sam.markduplicates;

import net.sf.picard.cmdline.CommandLineProgram;

public class MarkDuplicatesWithMateCigarTester extends AbstractMarkDuplicateFindingAlgorithmTester {

    public MarkDuplicatesWithMateCigarTester(final String scoringStrategy) {
        super();

        addArg("SCORING_STRATEGY=" + scoringStrategy);
    }

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicatesWithMateCigar(); }
}
