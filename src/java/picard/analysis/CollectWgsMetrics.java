package picard.analysis;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SamLocusIterator;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.Usage;
import picard.util.MathUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing experiments.
 *
 * @author tfennell
 */
public class CollectWgsMetrics extends CommandLineProgram {

    @Usage
    public final String usage = "Computes a number of metrics that are useful for evaluating coverage and performance of " +
            "whole genome sequencing experiments.";

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input SAM or BAM file.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output metrics file.")
    public File OUTPUT;

    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="The reference sequence fasta aligned to.")
    public File REFERENCE_SEQUENCE;

    @Option(shortName="MQ", doc="Minimum mapping quality for a read to contribute coverage.")
    public int MINIMUM_MAPPING_QUALITY = 20;

    @Option(shortName="Q", doc="Minimum base quality for a base to contribute coverage.")
    public int MINIMUM_BASE_QUALITY = 20;

    @Option(shortName="CAP", doc="Treat bases with coverage exceeding this value as if they had coverage at this value.")
    public int COVERAGE_CAP = 250;

    @Option(doc="For debugging purposes, stop after processing this many genomic bases.")
    public long STOP_AFTER = -1;


    @Option(doc="The number of threads to use")
    public int NUM_THREADS = 1;

    @Option(doc="The size of the genome within a chromosome to use when chunking")
    public int GENOMIC_CHUNK_SIZE = 25000000;

    @Option(doc="The time to wait in ms to poll threads when they are all busy")
    public int SLEEP_TIME = 100;

    private final Log log = Log.getInstance(CollectWgsMetrics.class);

    /** Metrics for evaluating the performance of whole genome sequencing experiments. */
    public static class WgsMetrics extends MetricBase {
        /** The number of non-N bases in the genome reference over which coverage will be evaluated. */
        public long GENOME_TERRITORY;
        /** The mean coverage in bases of the genome territory, after all filters are applied. */
        public double MEAN_COVERAGE;
        /** The standard deviation of coverage of the genome after all filters are applied. */
        public double SD_COVERAGE;
        /** The median coverage in bases of the genome territory, after all filters are applied. */
        public double MEDIAN_COVERAGE;
        /** The median absolute deviation of coverage of the genome after all filters are applied. */
        public double MAD_COVERAGE;

        /** The fraction of aligned bases that were filtered out because they were in reads with low mapping quality (default is < 20). */
        public double PCT_EXC_MAPQ;
        /** The fraction of aligned bases that were filtered out because they were in reads marked as duplicates. */
        public double PCT_EXC_DUPE;
        /** The fraction of aligned bases that were filtered out because they were in reads without a mapped mate pair. */
        public double PCT_EXC_UNPAIRED;
        /** The fraction of aligned bases that were filtered out because they were of low base quality (default is < 20). */
        public double PCT_EXC_BASEQ;
        /** The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads. */
        public double PCT_EXC_OVERLAP;
        /** The fraction of aligned bases that were filtered out because they would have raised coverage above the capped value (default cap = 250x). */
        public double PCT_EXC_CAPPED;
        /** The total fraction of aligned bases excluded due to all filters. */
        public double PCT_EXC_TOTAL;

        /** The fraction of bases that attained at least 5X sequence coverage in post-filtering bases. */
        public double PCT_5X;
        /** The fraction of bases that attained at least 10X sequence coverage in post-filtering bases. */
        public double PCT_10X;
        /** The fraction of bases that attained at least 15X sequence coverage in post-filtering bases. */
        public double PCT_15X;
        /** The fraction of bases that attained at least 20X sequence coverage in post-filtering bases. */
        public double PCT_20X;
        /** The fraction of bases that attained at least 25X sequence coverage in post-filtering bases. */
        public double PCT_25X;
        /** The fraction of bases that attained at least 30X sequence coverage in post-filtering bases. */
        public double PCT_30X;
        /** The fraction of bases that attained at least 40X sequence coverage in post-filtering bases. */
        public double PCT_40X;
        /** The fraction of bases that attained at least 50X sequence coverage in post-filtering bases. */
        public double PCT_50X;
        /** The fraction of bases that attained at least 60X sequence coverage in post-filtering bases. */
        public double PCT_60X;
        /** The fraction of bases that attained at least 70X sequence coverage in post-filtering bases. */
        public double PCT_70X;
        /** The fraction of bases that attained at least 80X sequence coverage in post-filtering bases. */
        public double PCT_80X;
        /** The fraction of bases that attained at least 90X sequence coverage in post-filtering bases. */
        public double PCT_90X;
        /** The fraction of bases that attained at least 100X sequence coverage in post-filtering bases. */
        public double PCT_100X;
    }

    public static void main(final String[] args) {
        new CollectWgsMetrics().instanceMainWithExit(args);
    }

    /**
     * This class are for metrics collectors to collect metrics over an interval using SamLocusIterator, such as per-locus
     * information.  This thread should be able to return a class via getCollection that returns an intermediate
     * collection of metrics used to compute the final set of metrics.
     *
     * At some point we should have this class implement the locus iteration directly.
     */
    abstract class CollectMetricsWithSamLocusIteratorThread<T extends MetricsCollection> extends Thread {
        protected final IntervalList intervalList;

        public CollectMetricsWithSamLocusIteratorThread() {
            throw new PicardException("Not implemented");
        }

        public CollectMetricsWithSamLocusIteratorThread(final IntervalList intervalList) {
            this.intervalList = intervalList;
        }

        /** This should only be called after the method has completed */
        public abstract T getCollection();

    }

    /**
     * The specific thread to compute the intermediate WGS metrics over a given interval. WE ass
     */
    class CollectWgsMetricsThread extends CollectMetricsWithSamLocusIteratorThread<WgsMetricsCollection> {

        private final WgsMetricsCollection collection;
        private final ProgressLogger progress;

        public CollectWgsMetricsThread(final IntervalList intervalList, final ProgressLogger progress) {
            super(intervalList);
            this.collection = new WgsMetricsCollection();
            this.progress = progress;
        }

        /**
         * This method assumes we are using only one chromosome.
         * At some point, this method should be put into a better form within CollectMetricsWithSamLocusIteratorThread */
        @Override
        public void run() {
            // Setup all the inputs
            final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
            final SAMFileReader in        = new SAMFileReader(INPUT);
            final SAMFileHeader header = in.getFileHeader();

            // Check that the input file and the reference are the same
            header.getSequenceDictionary().assertSameDictionary(refWalker.getSequenceDictionary());

            final SamLocusIterator iterator = new SamLocusIterator(in, this.intervalList);
            final List<SamRecordFilter> filters   = new ArrayList<SamRecordFilter>();
            filters.add(collection.mapqFilter);
            filters.add(collection.dupeFilter);
            filters.add(collection.pairFilter);
            filters.add(new SecondaryAlignmentFilter()); // Not a counting filter because we never want to count reads twice
            iterator.setSamFilters(filters);
            iterator.setEmitUncoveredLoci(true);
            iterator.setMappingQualityScoreCutoff(0); // Handled separately because we want to count bases
            iterator.setQualityScoreCutoff(0);        // Handled separately because we want to count bases
            iterator.setIncludeNonPfReads(false);

            final int max = COVERAGE_CAP;

            // Loop through all the loci
            int counter = 0; // count the # of loci seen
            int lastPosition = 0;
            int lastSequenceIndex = 0;
            while (iterator.hasNext()) {
                final SamLocusIterator.LocusInfo info = iterator.next();

                // Check that the reference is not N
                final ReferenceSequence ref = refWalker.get(info.getSequenceIndex());
                final byte base = ref.getBases()[info.getPosition()-1];
                if (base == 'N') continue;

                lastSequenceIndex = info.getSequenceIndex();

                // Figure out the coverage while not counting overlapping reads twice, and excluding various things
                final HashSet<String> readNames = new HashSet<String>(info.getRecordAndPositions().size());
                for (final SamLocusIterator.RecordAndOffset recs : info.getRecordAndPositions()) {
                    if (recs.getBaseQuality() < MINIMUM_BASE_QUALITY)                   { ++collection.basesExcludedByBaseq;   continue; }
                    if (!readNames.add(recs.getRecord().getReadName()))                 { ++collection.basesExcludedByOverlap; continue; }
                }

                final int depth = Math.min(readNames.size(), max);
                if (depth < readNames.size()) collection.basesExcludedByCapping += readNames.size() - max;
                collection.histogramArray[depth]++;

                lastPosition = info.getPosition();

                counter++;
            }
            // Record progress
            if (0 < counter) {
                progress.record(header.getSequence(lastSequenceIndex).getSequenceName(), lastPosition, counter);
            }

            iterator.close();
            in.close();

            // Construct and write the outputs
            for (int i=0; i<collection.histogramArray.length; ++i) {
                collection.histogram.increment(i, collection.histogramArray[i]);
            }
        }

        public WgsMetricsCollection getCollection() { return this.collection; }
    }

    /** This is the base class for all intermediate metrics collections.  This is used to get all necessary intermediate values
     * per thread, and then can be combined into the final set of metrics by the parent thread. */
    abstract class MetricsCollection<T extends MetricsCollection> {
        public abstract void add(final T collection);
    }

    /** This stores intermediate values of WgsMetrics, for us to scatter and gather */
    class WgsMetricsCollection extends MetricsCollection<WgsMetricsCollection> {
        final public Histogram<Integer> histogram;
        final public long[] histogramArray;
        final public CountingFilter dupeFilter;
        final public CountingFilter mapqFilter;
        final public CountingPairedFilter pairFilter;
        public long basesExcludedByBaseq   = 0;
        public long basesExcludedByOverlap = 0;
        public long basesExcludedByCapping = 0;

        public WgsMetricsCollection() {
            histogram = new Histogram<Integer>("coverage", "count");
            histogramArray   = new long[COVERAGE_CAP+1];
            dupeFilter       = new CountingDuplicateFilter();
            mapqFilter       = new CountingMapQFilter(MINIMUM_MAPPING_QUALITY);
            pairFilter       = new CountingPairedFilter();
        }

        @Override
        public void add(final WgsMetricsCollection collection) {
            this.histogram.addHistogram(collection.histogram);
            for (int i = 0; i < this.histogramArray.length; i++) {
                this.histogramArray[i] += collection.histogramArray[i];
            }
            this.dupeFilter.combineFilter(collection.dupeFilter);
            this.mapqFilter.combineFilter(collection.mapqFilter);
            this.pairFilter.combineFilter(collection.pairFilter);
            this.basesExcludedByBaseq += collection.basesExcludedByBaseq;
            this.basesExcludedByOverlap += collection.basesExcludedByOverlap;
            this.basesExcludedByCapping += collection.basesExcludedByCapping;
        }
    }

    public WgsMetricsCollection computeWgsMetricsCollection() {
        final SAMFileReader in = new SAMFileReader(INPUT);
        final SAMFileHeader header = in.getFileHeader();
        final SAMSequenceDictionary dict = header.getSequenceDictionary();
        in.close();

        // Check that the input file and the reference are the same
        final ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
        dict.assertSameDictionary(referenceSequenceFile.getSequenceDictionary());
        CloserUtil.close(referenceSequenceFile);

        final ProgressLogger progress = new ProgressLogger(log, 10000000, "Processed", "loci");

        // Create the list of intervals on which we process separately
        final List<Interval> intervals = new ArrayList<Interval>();
        int counter = 0;
        for (int i = 0; i < dict.getSequences().size(); i++) {
            final SAMSequenceRecord sequenceRecord = dict.getSequence(i);
            int startCoordinate = 0;
            final int sequenceLength = sequenceRecord.getSequenceLength();
            final String sequenceName = sequenceRecord.getSequenceName();
            while (startCoordinate < sequenceLength) {
                int endCoordinate = startCoordinate + GENOMIC_CHUNK_SIZE - 1;
                if (sequenceLength <= endCoordinate) endCoordinate = sequenceLength - 1;
                intervals.add(new Interval(sequenceName, startCoordinate+1, endCoordinate+1));
                startCoordinate += GENOMIC_CHUNK_SIZE;
                counter += intervals.size();
            }
            if (0 < STOP_AFTER && STOP_AFTER <= counter) break;
        }

        // Create the threads and run them
        final CollectWgsMetricsThread[] threads = new CollectWgsMetricsThread[NUM_THREADS];
        for (int i = 0; i < threads.length; i++) {
            threads[i] = null;
        }

        final WgsMetricsCollection collection = new WgsMetricsCollection();
        for (final Interval interval : intervals) {
            final IntervalList intervalList = new IntervalList(header);
            intervalList.add(interval);

            // find a thread that is done, then run another thread
            while (true) {
                boolean threadLaunched = false;
                for (int i = 0; i < threads.length; i++) {
                    if (null == threads[i] || !threads[i].isAlive()) {
                        if (null != threads[i]) { // the thread is dead, so combine the metrics
                            collection.add(threads[i].getCollection());
                        }
                        // create a new thread and start it
                        threads[i] = new CollectWgsMetricsThread(intervalList, progress);
                        threads[i].start();
                        threadLaunched = true;
                        break;
                    }
                }
                if (threadLaunched) break;
                // Sleep for a little while if we did not launch any threads
                try {
                    Thread.sleep(SLEEP_TIME);
                } catch (final InterruptedException e) {
                    throw new RuntimeException(e);
                }
            }
        }

        return collection;
    }

    private void finalizeMetrics(final WgsMetricsCollection wgs) {
        final WgsMetrics metrics = new WgsMetrics();
        metrics.GENOME_TERRITORY = (long) wgs.histogram.getSumOfValues();
        metrics.MEAN_COVERAGE    = wgs.histogram.getMean();
        metrics.SD_COVERAGE      = wgs.histogram.getStandardDeviation();
        metrics.MEDIAN_COVERAGE  = wgs.histogram.getMedian();
        metrics.MAD_COVERAGE     = wgs.histogram.getMedianAbsoluteDeviation();

        final long basesExcludedByDupes   = wgs.dupeFilter.getFilteredBases();
        final long basesExcludedByMapq    = wgs.mapqFilter.getFilteredBases();
        final long basesExcludedByPairing = wgs.pairFilter.getFilteredBases();
        final long basesExcludedByBaseq   = wgs.basesExcludedByBaseq;
        final long basesExcludedByOverlap = wgs.basesExcludedByOverlap;
        final long basesExcludedByCapping = wgs.basesExcludedByCapping;
        final double total             = wgs.histogram.getSum();
        final double totalWithExcludes = total + basesExcludedByDupes + basesExcludedByMapq + basesExcludedByPairing + basesExcludedByBaseq + basesExcludedByOverlap + basesExcludedByCapping;
        metrics.PCT_EXC_DUPE     = basesExcludedByDupes   / totalWithExcludes;
        metrics.PCT_EXC_MAPQ     = basesExcludedByMapq    / totalWithExcludes;
        metrics.PCT_EXC_UNPAIRED = basesExcludedByPairing / totalWithExcludes;
        metrics.PCT_EXC_BASEQ    = basesExcludedByBaseq   / totalWithExcludes;
        metrics.PCT_EXC_OVERLAP  = basesExcludedByOverlap / totalWithExcludes;
        metrics.PCT_EXC_CAPPED   = basesExcludedByCapping / totalWithExcludes;
        metrics.PCT_EXC_TOTAL    = (totalWithExcludes - total) / totalWithExcludes;

        metrics.PCT_5X     = MathUtil.sum(wgs.histogramArray, 5, wgs.histogramArray.length)   / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_10X    = MathUtil.sum(wgs.histogramArray, 10, wgs.histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_15X    = MathUtil.sum(wgs.histogramArray, 15, wgs.histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_20X    = MathUtil.sum(wgs.histogramArray, 20, wgs.histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_25X    = MathUtil.sum(wgs.histogramArray, 25, wgs.histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_30X    = MathUtil.sum(wgs.histogramArray, 30, wgs.histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_40X    = MathUtil.sum(wgs.histogramArray, 40, wgs.histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_50X    = MathUtil.sum(wgs.histogramArray, 50, wgs.histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_60X    = MathUtil.sum(wgs.histogramArray, 60, wgs.histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_70X    = MathUtil.sum(wgs.histogramArray, 70, wgs.histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_80X    = MathUtil.sum(wgs.histogramArray, 80, wgs.histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_90X    = MathUtil.sum(wgs.histogramArray, 90, wgs.histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_100X   = MathUtil.sum(wgs.histogramArray, 100, wgs.histogramArray.length) / (double) metrics.GENOME_TERRITORY;

        final MetricsFile<WgsMetrics, Integer> out = getMetricsFile();
        out.addMetric(metrics);
        out.addHistogram(wgs.histogram);
        out.write(OUTPUT);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        // Compute the metrics and write them out
        finalizeMetrics(computeWgsMetricsCollection());

        return 0;
    }
}

/**
 * A SamRecordFilter that counts the number of aligned bases in the reads which it filters out. Abstract and designed
 * to be subclassed to implement the desired filter.
 */
abstract class CountingFilter implements SamRecordFilter {
    private long filteredRecords = 0;
    private long filteredBases = 0;

    /** Gets the number of records that have been filtered out thus far. */
    public long getFilteredRecords() { return this.filteredRecords; }

    /** Gets the number of bases that have been filtered out thus far. */
    public long getFilteredBases() { return this.filteredBases; }

    @Override public final boolean filterOut(final SAMRecord record) {
        final boolean filteredOut = reallyFilterOut(record);
        if (filteredOut) {
            ++filteredRecords;
            for (final AlignmentBlock block : record.getAlignmentBlocks()) {
                this.filteredBases += block.getLength();
            }
        }
        return filteredOut;
    }

    abstract public boolean reallyFilterOut(final SAMRecord record);

    @Override public boolean filterOut(final SAMRecord first, final SAMRecord second) {
        throw new UnsupportedOperationException();
    }

    public void combineFilter(final CountingFilter filter) {
        this.filteredRecords += filter.filteredRecords;
        this.filteredBases += filter.filteredBases;
    }
}

/** Counting filter that discards reads that have been marked as duplicates. */
class CountingDuplicateFilter extends CountingFilter {
    @Override public boolean reallyFilterOut(final SAMRecord record) { return record.getDuplicateReadFlag(); }
}

/** Counting filter that discards reads below a configurable mapping quality threshold. */
class CountingMapQFilter extends CountingFilter {
    private final int minMapq;
    CountingMapQFilter(final int minMapq) { this.minMapq = minMapq; }
    @Override public boolean reallyFilterOut(final SAMRecord record) { return record.getMappingQuality() < minMapq; }
}

/** Counting filter that discards reads that are unpaired in sequencing and paired reads who's mates are not mapped. */
class CountingPairedFilter extends CountingFilter {
    @Override public boolean reallyFilterOut(final SAMRecord record) { return !record.getReadPairedFlag() || record.getMateUnmappedFlag(); }
}

