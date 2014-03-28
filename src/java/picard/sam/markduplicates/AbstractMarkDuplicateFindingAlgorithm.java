package picard.sam.markduplicates;

import picard.PicardException;
import picard.cmdline.CommandLineParser;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.Usage;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.MergingSamRecordIterator;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.util.*;

/**
 * Abstract class that holds parameters and methods common to classes that perform duplicate
 * detection and/or marking within SAM/BAM files.
 *
 * @author Nils Homer
 */
public abstract class AbstractMarkDuplicateFindingAlgorithm extends AbstractDuplicateFindingAlgorithm {

    @Usage
    public final String USAGE =
            CommandLineParser.getStandardUsagePreamble(getClass()) +
                    "Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules. " +
                    "All records are then written to the output file with the duplicate records flagged.";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc="One or more input SAM or BAM files to analyze. Must be coordinate sorted.")
    public List<File> INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc="The output file to write marked records to")
    public File OUTPUT;

    @Option(shortName="M",
            doc="File to write duplication metrics to")
    public File METRICS_FILE;

    @Option(shortName=StandardOptionDefinitions.PROGRAM_RECORD_ID_SHORT_NAME,
            doc="The program record ID for the @PG record(s) created by this program. Set to null to disable " +
                    "PG record creation.  This string may have a suffix appended to avoid collision with other " +
                    "program record IDs.",
            optional=true)
    public String PROGRAM_RECORD_ID = "MarkDuplicates";

    @Option(shortName="PG_VERSION",
            doc="Value of VN tag of PG record to be created. If not specified, the version will be detected automatically.",
            optional=true)
    public String PROGRAM_GROUP_VERSION;

    @Option(shortName="PG_COMMAND",
            doc="Value of CL tag of PG record to be created. If not supplied the command line will be detected automatically.",
            optional=true)
    public String PROGRAM_GROUP_COMMAND_LINE;

    @Option(shortName="PG_NAME",
            doc="Value of PN tag of PG record to be created.")
    public String PROGRAM_GROUP_NAME = getClass().getSimpleName();

    @Option(shortName="CO",
            doc="Comment(s) to include in the output file's header.",
            optional=true)
    public List<String> COMMENT = new ArrayList<String>();

    @Option(doc="If true do not write duplicates to the output file instead of writing them with appropriate flags set.")
    public boolean REMOVE_DUPLICATES = false;

    @Option(shortName=StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME,
            doc="If true, assume that the input file is coordinate sorted even if the header says otherwise.")
    public boolean ASSUME_SORTED = false;

    // Gather all PG IDs seen in merged input files in first pass.  These are gathered for two reasons:
    // - to know how many different PG records to create to represent this program invocation.
    // - to know what PG IDs are already used to avoid collisions when creating new ones.
    // Note that if there are one or more records that do not have a PG tag, then a null value
    // will be stored in this set.
    protected final Set<String> pgIdsSeen = new HashSet<String>();

    protected Map<String, String> getChainedPgIds(final SAMFileHeader outputHeader) {
        final Map<String, String> chainedPgIds;
        // Generate new PG record(s)
        if (PROGRAM_RECORD_ID != null) {
            final PgIdGenerator pgIdGenerator = new PgIdGenerator(outputHeader);
            if (PROGRAM_GROUP_VERSION == null) {
                PROGRAM_GROUP_VERSION = this.getVersion();
            }
            if (PROGRAM_GROUP_COMMAND_LINE == null) {
                PROGRAM_GROUP_COMMAND_LINE = this.getCommandLine();
            }
            chainedPgIds = new HashMap<String, String>();
            for (final String existingId : this.pgIdsSeen) {
                final String newPgId = pgIdGenerator.getNonCollidingId(PROGRAM_RECORD_ID);
                chainedPgIds.put(existingId, newPgId);
                final SAMProgramRecord programRecord = new SAMProgramRecord(newPgId);
                programRecord.setProgramVersion(PROGRAM_GROUP_VERSION);
                programRecord.setCommandLine(PROGRAM_GROUP_COMMAND_LINE);
                programRecord.setProgramName(PROGRAM_GROUP_NAME);
                programRecord.setPreviousProgramGroupId(existingId);
                outputHeader.addProgramRecord(programRecord);
            }
        } else {
            chainedPgIds = null;
        }
        return chainedPgIds;
    }

    protected void writeMetrics(final Map<String,DuplicationMetrics> metricsByLibrary,
                                final Histogram<Short> opticalDupesByLibraryId,
                                final Map<String,Short> libraryIds) {
        // Write out the metrics
        final MetricsFile<DuplicationMetrics,Double> file = getMetricsFile();
        for (final Map.Entry<String,DuplicationMetrics> entry : metricsByLibrary.entrySet()) {
            final String libraryName = entry.getKey();
            final DuplicationMetrics metrics = entry.getValue();

            metrics.READ_PAIRS_EXAMINED = metrics.READ_PAIRS_EXAMINED / 2;
            metrics.READ_PAIR_DUPLICATES = metrics.READ_PAIR_DUPLICATES / 2;

            // Add the optical dupes to the metrics
            final Short libraryId = libraryIds.get(libraryName);
            if (libraryId != null) {
                final Histogram<Short>.Bin bin = opticalDupesByLibraryId.get(libraryId);
                if (bin != null) {
                    metrics.READ_PAIR_OPTICAL_DUPLICATES = (long) bin.getValue();
                }
            }
            metrics.calculateDerivedMetrics();
            file.addMetric(metrics);
        }

        if (metricsByLibrary.size() == 1) {
            file.setHistogram(metricsByLibrary.values().iterator().next().calculateRoiHistogram());
        }

        file.write(METRICS_FILE);
    }

    static class PgIdGenerator {
        private int recordCounter;

        private final Set<String> idsThatAreAlreadyTaken = new HashSet<String>();

        PgIdGenerator(final SAMFileHeader header) {
            for (final SAMProgramRecord pgRecord : header.getProgramRecords()) {
                idsThatAreAlreadyTaken.add(pgRecord.getProgramGroupId());
            }
            recordCounter = idsThatAreAlreadyTaken.size();
        }

        String getNonCollidingId(final String recordId) {
            if(!idsThatAreAlreadyTaken.contains(recordId)) {
                // don't remap 1st record. If there are more records
                // with this id, they will be remapped in the 'else'.
                idsThatAreAlreadyTaken.add(recordId);
                ++recordCounter;
                return recordId;
            } else {
                String newId;
                // Below we tack on one of roughly 1.7 million possible 4 digit base36 at random. We do this because
                // our old process of just counting from 0 upward and adding that to the previous id led to 1000s of
                // calls idsThatAreAlreadyTaken.contains() just to resolve 1 collision when merging 1000s of similarly
                // processed bams.
                while(idsThatAreAlreadyTaken.contains(newId = recordId + "." + SamFileHeaderMerger.positiveFourDigitBase36Str(recordCounter++)));

                idsThatAreAlreadyTaken.add( newId );
                return newId;
            }

        }
    }

    /** Little class used to package up a header and an iterable/iterator. */
    protected static final class SamHeaderAndIterator {
        final SAMFileHeader header;
        final CloseableIterator<SAMRecord> iterator;

        protected SamHeaderAndIterator(final SAMFileHeader header, final CloseableIterator<SAMRecord> iterator) {
            this.header = header;
            this.iterator = iterator;
        }
    }

    /**
     * Since this may read it's inputs more than once this method does all the opening
     * and checking of the inputs.
     */
    protected SamHeaderAndIterator openInputs() {
        final List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>(INPUT.size());
        final List<SAMFileReader> readers = new ArrayList<SAMFileReader>(INPUT.size());

        for (final File f : INPUT) {
            final SAMFileReader reader = new SAMFileReader(f);
            final SAMFileHeader header = reader.getFileHeader();

            if (!ASSUME_SORTED && header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new PicardException("Input file " + f.getAbsolutePath() + " is not coordinate sorted.");
            }

            headers.add(header);
            readers.add(reader);
        }

        if (headers.size() == 1) {
            return new SamHeaderAndIterator(headers.get(0), readers.get(0).iterator());
        }
        else {
            final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, headers, false);
            final MergingSamRecordIterator iterator = new MergingSamRecordIterator(headerMerger, readers, ASSUME_SORTED);
            return new SamHeaderAndIterator(headerMerger.getMergedHeader(), iterator);
        }
    }

    /**
     * Looks through the set of reads and identifies how many of the duplicates are
     * in fact optical duplicates, and stores the data in the instance level histogram.
     */
    protected static void trackOpticalDuplicates(final List<? extends ReadEnds> list,
                                                 final OpticalDuplicateFinder opticalDuplicateFinder,
                                                 final Histogram<Short> opticalDupesByLibraryId) {

        final boolean[] opticalDuplicateFlags = opticalDuplicateFinder.findOpticalDuplicates(list);

        int opticalDuplicates = 0;
        for (final boolean b: opticalDuplicateFlags) if (b) ++opticalDuplicates;
        if (opticalDuplicates > 0) {
            opticalDupesByLibraryId.increment(list.get(0).libraryId, opticalDuplicates);
        }
    }

    /**
     * Gets the library name from the header for the record. If the RG tag is not present on
     * the record, or the library isn't denoted on the read group, a constant string is
     * returned.
     */
    public static String getLibraryName(final SAMFileHeader header, final SAMRecord rec) {
        final String readGroupId = (String) rec.getAttribute("RG");

        if (readGroupId != null) {
            final SAMReadGroupRecord rg = header.getReadGroup(readGroupId);
            if (rg != null) {
                final String libraryName = rg.getLibrary();
                if (null != libraryName) return libraryName;
            }
        }

        return "Unknown Library";
    }
}
