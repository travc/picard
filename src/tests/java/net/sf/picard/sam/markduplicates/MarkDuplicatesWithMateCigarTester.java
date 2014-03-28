package net.sf.picard.sam.markduplicates;

import net.sf.picard.cmdline.CommandLineProgram;

public class MarkDuplicatesWithMateCigarTester extends AbstractMarkDuplicateFindingAlgorithmTester {

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicatesWithMateCigar(); }
}
