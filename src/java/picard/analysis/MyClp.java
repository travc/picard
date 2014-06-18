package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.PicardException;
import picard.analysis.directed.InsertSizeMetricsCollector;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.Unrollable;
import picard.cmdline.Usage;
import picard.util.RExecutor;

import java.io.File;
import java.util.ArrayList;
import java.util.Set;

/**
 * @author mccowan
 */
public class MyClp extends CommandLineProgram {
    private static class StringBin extends ArrayList<String> {
        
    }
    
    final Log logger = Log.getInstance(MyClp.class);
            
    
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Testing!";

    @Option
    public StringBin STRINGS = new StringBin();
    
    @Option
    public Unrollable<Integer> I = Unrollable.of(Integer.class);

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new MyClp().instanceMainWithExit(argv);
    }

    @Override
    protected int doWork() {
        for (final Integer integer : I.unrolled()) {
            logger.info(integer);   
        }
        
        return 0;
    }
}