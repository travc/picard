package net.sf.picard.sam.markduplicates;

import net.sf.picard.cmdline.CommandLineProgram;

public class MarkDuplicatesTester extends AbstractMarkDuplicateFindingAlgorithmTester {

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicates(); }
}
