/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam.markduplicates;

import picard.cmdline.Option;
import htsjdk.samtools.util.Histogram;
import picard.sam.DuplicationMetrics;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.*;

import java.io.*;
import java.util.*;

/**
 * An even better duplication marking algorithm that handles all cases including clipped
 * and gapped alignments.
 *
 * This tool differs with MarkDuplicates as it may break ties differently.  Furthermore,
 * as it is a one-pass algorithm, it cannot know the program records contained in the file
 * that should be chained in advance. Therefore, it will only create chain from program records
 * that are have no previous program record in any other program record (it is the last in the
 * given chain).  If a read is encountered without a program record, or not one as previously
 * defined, it will not be updated.
 *
 * @author Nils Homer
 */
public class MarkDuplicatesWithMateCigar extends AbstractMarkDuplicateFindingAlgorithm {
    private final Log log = Log.getInstance(MarkDuplicatesWithMateCigar.class);

    @Option(doc="The minimum distance to buffer records to account for clipping on the 5' end of the records." +
            "Set this number to -1 to use twice the first read's read length (or 100, whichever is smaller).", optional=true)
    public int MINIMUM_DISTANCE = -1;

    @Option(doc="The scoring strategy to select which record should be not called a duplicate among comparable"
            + " (potential duplicate) records.", optional = true)
    public MarkDuplicatesWithMateCigarIterator.ScoringStrategy SCORING_STRATEGY = MarkDuplicatesWithMateCigarIterator.ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH_THEN_MAPQ_THEN_READ_NAME;

    @Option(doc="Skip record pairs with no mate cigar and include them in the output.")
    boolean SKIP_PAIRS_WITH_NO_MATE_CIGAR = true;

    private boolean warnedNullProgramRecords = false;
    private boolean warnedMissingProgramRecords = false;

    /** Stock main method. */
    public static void main(final String[] args) {
        System.exit(new MarkDuplicatesWithMateCigar().instanceMain(args));
    }

    /**
     * Main work method.
     */
    protected int doWork() {
        for (final File f : INPUT) IOUtil.assertFileIsReadable(f);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(METRICS_FILE);

        // Open the inputs
        final SamHeaderAndIterator headerAndIterator = openInputs();
        final SAMFileHeader header = headerAndIterator.header;

        // Create the output header
        final SAMFileHeader outputHeader = header.clone();
        outputHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (final String comment : COMMENT) outputHeader.addComment(comment);

        // Since this is one-pass, unlike MarkDuplicates, we cannot only chain together program
        // group records we have seen, we have to assume all of them may be seen.  We can perhaps
        // filter out any program groups which have been referenced previously.
        setPGIdsSeen(outputHeader);
        // Key: previous PG ID on a SAM Record (or null).  Value: New PG ID to replace it.
        final Map<String, String> chainedPgIds = getChainedPgIds(outputHeader);

        // Open the output
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(outputHeader,
                true,
                OUTPUT);

        // Create the mark duplicate iterator
        final MarkDuplicatesWithMateCigarIterator iterator = new MarkDuplicatesWithMateCigarIterator(headerAndIterator.header,
                headerAndIterator.iterator,
                this.opticalDuplicateFinder,
                this.MINIMUM_DISTANCE,
                this.REMOVE_DUPLICATES,
                this.SCORING_STRATEGY,
                this.SKIP_PAIRS_WITH_NO_MATE_CIGAR);

        // progress logger!
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");

        // Go through the records
        for (final SAMRecord record : new IterableAdapter<SAMRecord>(iterator)) {
            progress.record(record);

            // Update the program record if necessary
            if (PROGRAM_RECORD_ID != null) {
                final String pgId = record.getStringAttribute(SAMTag.PG.name());
                if (null == pgId) {
                    if (!warnedNullProgramRecords) {
                        warnedNullProgramRecords = true;
                        log.warn("Encountered a record with no program record, program group chaining will not occur for this read: " + record);
                    } // else already warned!
                }
                else if (!chainedPgIds.containsKey(pgId)) {
                    if (!warnedMissingProgramRecords) {
                        warnedMissingProgramRecords = true;
                        log.warn("Encountered a record with a intermediate program record, program group chaining will not occur for this read: " + record);
                    } // else already warned!
                }
                else {
                    record.setAttribute(SAMTag.PG.name(), chainedPgIds.get(pgId));
                }
            }

            out.addAlignment(record);
        }

        iterator.close();
        out.close();

        // For convenience to reference
        final Map<String,Short> libraryIds = iterator.getLibraryIds();
        final Histogram<Short> opticalDupesByLibraryId = iterator.getOpticalDupesByLibraryId();
        final Map<String,DuplicationMetrics> metricsByLibrary = iterator.getMetricsByLibrary();

        // Log info
        log.info("Found " + iterator.getNumRecordsWithNoMateCigar() + " records with no mate cigar optional tag.");
        log.info("Marking " + iterator.getNumDuplicates() + " records as duplicates.");
        log.info("Found " + ((long) opticalDupesByLibraryId.getSumOfValues()) + " optical duplicate clusters.");

        // Write out the metrics
        writeMetrics(metricsByLibrary, opticalDupesByLibraryId, libraryIds);

        return 0;
    }

    private void setPGIdsSeen(final SAMFileHeader header) {
        final Set<String> pgIdsSeenAsPrevious = new HashSet<String>();

        // get all program record ids that are mentioned as previously seen
        for (final SAMProgramRecord samProgramRecord : header.getProgramRecords()) {
            final String previousProgramGroupID = samProgramRecord.getPreviousProgramGroupId();
            if (null != previousProgramGroupID) pgIdsSeenAsPrevious.add(previousProgramGroupID);
        }

        // ignore those that were previously seen
        for (final SAMProgramRecord samProgramRecord : header.getProgramRecords()) {
            final String pgId = samProgramRecord.getId();
            if (!pgIdsSeenAsPrevious.contains(pgId)) this.pgIdsSeen.add(pgId);
        }
    }
}
