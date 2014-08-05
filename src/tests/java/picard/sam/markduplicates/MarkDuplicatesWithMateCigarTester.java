package picard.sam.markduplicates;

import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import picard.cmdline.CommandLineProgram;

public class MarkDuplicatesWithMateCigarTester extends AbstractMarkDuplicateFindingAlgorithmTester {

    public MarkDuplicatesWithMateCigarTester() {
        // NB: to be equivalent to MarkDuplicates we need to use SUM_OF_BASE_QUALITIES
        super(ScoringStrategy.SUM_OF_BASE_QUALITIES);

        addArg("MAX_RECORDS_IN_RAM=1000");
        addArg("BLOCK_SIZE=250");
    }

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicatesWithMateCigar(); }
}
