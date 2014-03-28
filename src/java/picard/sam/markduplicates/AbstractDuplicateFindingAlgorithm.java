package picard.sam.markduplicates;

import picard.cmdline.*;
import htsjdk.samtools.util.Log;

/**
 * Abstract class that holds parameters and methods common to classes that perform duplicate
 * detection and/or marking within SAM/BAM files.
 *
 * @author Tim Fennell
 */
public abstract class AbstractDuplicateFindingAlgorithm extends CommandLineProgram {
    protected static Log LOG = Log.getInstance(AbstractDuplicateFindingAlgorithm.class);

    @Option(doc="Regular expression that can be used to parse read names in the incoming SAM file. Read names are " +
            "parsed to extract three variables: tile/region, x coordinate and y coordinate. These values are used " +
            "to estimate the rate of optical duplication in order to give a more accurate estimated library size. " +
            "The regular expression should contain three capture groups for the three variables, in order. " +
            "It must match the entire read name. " +
            "Note that if the default regex is specified, a regex match is not actually done, but instead the read name " +
            " is split on colon character and the 2nd, 3rd and 4th elements are assumed to be tile, x and y values.")
    public String READ_NAME_REGEX = OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX;

    @Option(doc="The maximum offset between two duplicte clusters in order to consider them optical duplicates. This " +
            "should usually be set to some fairly small number (e.g. 5-10 pixels) unless using later versions of the " +
            "Illumina pipeline that multiply pixel values by 10, in which case 50-100 is more normal.")
    public int OPTICAL_DUPLICATE_PIXEL_DISTANCE = OpticalDuplicateFinder.DEFAULT_OPTICAL_DUPLICATE_DISTANCE;

    // The tool with which to find optical duplicates
    protected OpticalDuplicateFinder opticalDuplicateFinder = null;

    @Override
    protected String[] customCommandLineValidation() {
        this.opticalDuplicateFinder = new OpticalDuplicateFinder(READ_NAME_REGEX, OPTICAL_DUPLICATE_PIXEL_DISTANCE, LOG);
        return super.customCommandLineValidation();
    }
}
