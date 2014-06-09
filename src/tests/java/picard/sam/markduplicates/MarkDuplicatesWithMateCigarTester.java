package picard.sam.markduplicates;

import picard.cmdline.CommandLineProgram;
import picard.sam.markduplicates.MarkDuplicatesWithMateCigar;

public class MarkDuplicatesWithMateCigarTester extends AbstractMarkDuplicateFindingAlgorithmTester {

    public MarkDuplicatesWithMateCigarTester(final String scoringStrategy) {
        super();

        addArg("SCORING_STRATEGY=" + scoringStrategy);
        addArg("MAX_RECORDS_IN_RAM=1000");
        addArg("BLOCK_SIZE=250");
    }

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicatesWithMateCigar(); }
}
