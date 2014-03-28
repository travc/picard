package picard.sam.markduplicates;

import picard.cmdline.Option;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Contains methods for finding optical duplicates.
 * @author Tim Fennell
 * @author Nils Homer
 */
public class OpticalDuplicateFinder {

    public static final String DEFAULT_READ_NAME_REGEX = "[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*".intern();

    public static final int DEFAULT_OPTICAL_DUPLICATE_DISTANCE = 100;

    public String readNameRegex;

    public int opticalDuplicatePixelDistance;

    private Pattern readNamePattern;

    private boolean warnedAboutRegexNotMatching = false;

    private Log log;

    public OpticalDuplicateFinder() {
        this(DEFAULT_READ_NAME_REGEX, DEFAULT_OPTICAL_DUPLICATE_DISTANCE);
    }

    public OpticalDuplicateFinder(int opticalDuplicatePixelDistance) {
        this(DEFAULT_READ_NAME_REGEX, opticalDuplicatePixelDistance);
    }

    public OpticalDuplicateFinder(String readNameRegex) {
        this(readNameRegex, DEFAULT_OPTICAL_DUPLICATE_DISTANCE);
    }

    public OpticalDuplicateFinder(String readNameRegex, int opticalDuplicatePixelDistance) {
        this(readNameRegex, opticalDuplicatePixelDistance, null);
    }

    public OpticalDuplicateFinder(String readNameRegex, int opticalDuplicatePixelDistance, Log log) {
        this.readNameRegex = readNameRegex;
        this.opticalDuplicatePixelDistance = opticalDuplicatePixelDistance;
        this.log = log;
    }

    /**
     * Small interface that provides access to the physical location information about a cluster.
     * All values should be defaulted to -1 if unavailable.  ReadGroup and Tile should only allow
     * non-zero positive integers, x and y coordinates may be negative.
     */
    public static interface PhysicalLocation {
        short getReadGroup();
        void  setReadGroup(short rg);
        short  getTile();
        void  setTile(short tile);
        short getX();
        void  setX(short x);
        short getY();
        void  setY(short y);
    }

    /**
     * Method used to extract tile/x/y from the read name and add it to the PhysicalLocation so that it
     * can be used later to determine optical duplication
     *
     * @param readName the name of the read/cluster
     * @param loc the object to add tile/x/y to
     * @return true if the read name contained the information in parsable form, false otherwise
     */
    private final int[] tmpLocationFields = new int[10];
    boolean addLocationInformation(final String readName, final PhysicalLocation loc) {
        // Optimized version if using the default read name regex (== used on purpose):
        if (this.readNameRegex == this.DEFAULT_READ_NAME_REGEX) {
            final int fields = getRapidDefaultReadNameRegexSplit(readName, ':', tmpLocationFields);
            if (fields < 5) {
                if (null != log && !this.warnedAboutRegexNotMatching) {
                    this.log.warn(String.format("Default READ_NAME_REGEX '%s' did not match read name '%s'.  " +
                            "You may need to specify a READ_NAME_REGEX in order to correctly identify optical duplicates.  " +
                            "Note that this message will not be emitted again even if other read names do not match the regex.",
                            this.readNameRegex, readName));
                    this.warnedAboutRegexNotMatching = true;
                }
                return false;
            }
            loc.setTile((short) tmpLocationFields[2]);
            loc.setX((short) tmpLocationFields[3]);
            loc.setY((short) tmpLocationFields[4]);
            return true;
        }
        else if (this.readNameRegex == null) {
            return false;
        }
        else {
            // Standard version that will use the regex
            if (this.readNamePattern == null) this.readNamePattern = Pattern.compile(this.readNameRegex);

            final Matcher m = this.readNamePattern.matcher(readName);
            if (m.matches()) {
                loc.setTile((short) Integer.parseInt(m.group(1)));
                loc.setX((short) Integer.parseInt(m.group(2)));
                loc.setY((short) Integer.parseInt(m.group(3)));
                return true;
            }
            else {
                if (null != log && !this.warnedAboutRegexNotMatching) {
                    this.log.warn(String.format("READ_NAME_REGEX '%s' did not match read name '%s'.  Your regex may not be correct.  " +
                            "Note that this message will not be emitted again even if other read names do not match the regex.",
                            this.readNameRegex, readName));
                    warnedAboutRegexNotMatching = true;
                }
                return false;
            }
        }
    }


    /**
     * Single pass method to parse the read name for the default regex.  This will only insert the 2nd to the 4th
     * tokens (inclusive).  It will also stop after the fifth token has been successfully parsed.
     */
    protected int getRapidDefaultReadNameRegexSplit(String readName, char delim, int[] tokens) {
        int tokensIdx = 0;
        int prevIdx = 0;
        for (int i = 0; i < readName.length(); i++) {
            if (readName.charAt(i) == delim) {
                if (1 < tokensIdx && tokensIdx < 5) tokens[tokensIdx] = rapidParseInt(readName.substring(prevIdx, i)); // only fill in 2-4 inclusive
                tokensIdx++;
                if (4 < tokensIdx) return tokensIdx; // early return, only consider the first five tokens
                prevIdx = i+1;
            }
        }
        if (prevIdx < readName.length()) {
            if (1 < tokensIdx && tokensIdx < 5) tokens[tokensIdx] = rapidParseInt(readName.substring(prevIdx, readName.length())); // only fill in 2-4 inclusive
            tokensIdx++;
        }
        return tokensIdx;
    }

    /**
     * Very specialized method to rapidly parse a sequence of digits from a String up until the first
     * non-digit character.
     */
    protected final int rapidParseInt(final String input) {
        final int len = input.length();
        int val = 0;
        int i = 0;
        boolean isNegative = false;

        if (0 < len  && '-' == input.charAt(0)) {
            i = 1;
            isNegative = true;
        }

        for (;i < len; ++i) {
            final char ch = input.charAt(i);
            if (Character.isDigit(ch)) {
                val = (val*10) + (ch-48);
            }
            else {
                break;
            }
        }

        if (isNegative) val = -val;

        return val;
    }

    /**
     * Finds which reads within the list of duplicates are likely to be optical duplicates of
     * one another.
     *
     * Note: this method will perform a sort() of the list; if it is imperative that the list be
     * unmodified a copy of the list should be passed to this method.
     *
     * @param list a list of reads that are determined to be duplicates of one another
     * @return a boolean[] of the same length as the incoming list marking which reads are optical duplicates
     */
    boolean[] findOpticalDuplicates(final List<? extends PhysicalLocation> list) {
        final int length = list.size();
        final boolean[] opticalDuplicateFlags = new boolean[length];

        Collections.sort(list, new Comparator<PhysicalLocation>() {
            public int compare(final PhysicalLocation lhs, final PhysicalLocation rhs) {
                int retval = lhs.getReadGroup() - rhs.getReadGroup();
                if (retval == 0) retval = lhs.getTile() - rhs.getTile();
                if (retval == 0) retval = lhs.getX() - rhs.getX();
                if (retval == 0) retval = lhs.getY() - rhs.getY();
                return retval;
            }
        });

        outer: for (int i=0; i<length; ++i) {
            PhysicalLocation lhs = list.get(i);
            if (lhs.getTile() < 0) continue;

            for (int j=i+1; j<length; ++j) {
                PhysicalLocation rhs = list.get(j);

                if (lhs.getReadGroup() != rhs.getReadGroup()) continue outer;
                if (lhs.getTile() != rhs.getTile()) continue outer;
                if (rhs.getX() > lhs.getX() + this.opticalDuplicatePixelDistance) continue outer;

                if (Math.abs(lhs.getY()  - rhs.getY()) <= this.opticalDuplicatePixelDistance) {
                    opticalDuplicateFlags[j] = true;
                }
            }
        }
        return opticalDuplicateFlags;
    }
}
