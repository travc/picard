package picard.sam.markduplicates;

import htsjdk.samtools.DuplicateScoringStrategy;
import picard.cmdline.CommandLineProgram;
import picard.sam.markduplicates.MarkDuplicates;

public class MarkDuplicatesTester extends AbstractMarkDuplicateFindingAlgorithmTester {

    public MarkDuplicatesTester() {
        super(DuplicateScoringStrategy.ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH);
    }

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicates(); }
}
