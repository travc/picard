package picard.sam.markduplicates;

import picard.cmdline.CommandLineProgram;
import picard.sam.markduplicates.MarkDuplicatesWithMateCigar;

public class MarkDuplicatesWithMateCigarTester extends AbstractMarkDuplicateFindingAlgorithmTester {

    public MarkDuplicatesWithMateCigarTester(final String scoringStrategy) {
        super();

        addArg("SCORING_STRATEGY=" + scoringStrategy);
    }

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicatesWithMateCigar(); }
}
