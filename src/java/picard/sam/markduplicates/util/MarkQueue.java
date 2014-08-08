package picard.sam.markduplicates.util;

import htsjdk.samtools.DuplicateScoringStrategy;
import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import picard.PicardException;
import picard.sam.DuplicationMetrics;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;
import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;

import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * This is the mark queue.
 *
 * This stores a current set of read ends that need to be duplicate marked.  It only stores internally the "best" read end for a given
 * possible duplicate location, preferring to perform duplicate marking as read ends come in, rather than wait for all "comparable"
 * read ends to arrive.  This reduces the memory footprint of this data structure.
 */
public class MarkQueue {

    /** Comparator to order the mark queue set.  The set of all the read ends that have are compared to be the same should
     * be used for duplicate marking. */
    private class MarkQueueComparator implements Comparator<ReadEndsMC> {
        public int compare(final ReadEndsMC lhs, final ReadEndsMC rhs) {
            int retval = lhs.libraryId - rhs.libraryId;
            if (retval == 0) retval = lhs.read1Sequence - rhs.read1Sequence;
            if (retval == 0) retval = lhs.read1Coordinate - rhs.read1Coordinate;
            if (retval == 0) retval = rhs.orientation - lhs.orientation; // to get pairs first
            if (retval == 0) retval = lhs.read2Sequence - rhs.read2Sequence;
            if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
            return retval;
        }
    }

    /** Comparator for ReadEndsMC that orders by read1 position then pair orientation then read2 position. */
    // Could be a Singleton, but no static variables in inner classes.
    class ReadEndsMCComparator implements Comparator<ReadEndsMC> {
        private final ScoringStrategy duplicateScoringStrategy;

        public ReadEndsMCComparator(final ScoringStrategy duplicateScoringStrategy) {
            this.duplicateScoringStrategy = duplicateScoringStrategy;
        }

        public int compare(final ReadEndsMC lhs, final ReadEndsMC rhs) {
            int retval = lhs.libraryId - rhs.libraryId;
            if (retval == 0) retval = lhs.read1Sequence - rhs.read1Sequence;
            if (retval == 0) retval = lhs.read1Coordinate - rhs.read1Coordinate;
            if (retval == 0) retval = rhs.orientation - lhs.orientation; // IMPORTANT: reverse the order to get pairs first
            if (retval == 0 && lhs.isPaired() != rhs.isPaired()) return lhs.isPaired() ? -1 : 1; // unpaired goes first...
            if (retval == 0) retval = lhs.hasUnmapped - rhs.hasUnmapped;
            if (retval == 0) retval = lhs.read2Sequence - rhs.read2Sequence;
            if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
            // TODO: cache the scores?
            if (retval == 0) retval = DuplicateScoringStrategy.compare(lhs.getRecord(), rhs.getRecord(), this.duplicateScoringStrategy, true);
            if (retval == 0) retval = lhs.getRecordReadName().compareTo(rhs.getRecordReadName());

            return retval;
        }
    }

    private int toMarkQueueMinimumDistance = -1;

    private int numDuplicates = 0;

    /** The set of all read ends sorted by 5' start unclipped position.  Some read ends in this set may eventually be duplicates. */
    private final TreeSet<ReadEndsMC> set = new TreeSet<ReadEndsMC>(new MarkQueueComparator());

    /** Reads in the main set may occasionally have mates with the same chromosome, coordinate, and orientation, causing collisions
     * We store the 'best' end of the mate pair in the main set, and the other end in this set.  We only remove from this.pairSet when
     * we remove something from this.set.
     */
    private final TreeSet<ReadEndsMC> pairSet = new TreeSet<ReadEndsMC>(new MarkQueueComparator());

    /** Physical locations used for optical duplicate tracking.  This is only stored for paired end reads where both ends are mapped,
     * and when we see the first mate.
     */
    private final Map<ReadEndsMC, LocationSet> locations = new HashMap<ReadEndsMC, LocationSet>();

    /** If we have two items that are the same with respect to being in the "set", then we must choose one.  The "one" will
     * eventually be the end that is not marked as a duplicate in most cases (see poll() for the exceptions).
     */
    private final Comparator<ReadEndsMC> comparator;

    /** temporary so we do not need to create many objects */
    private ReadEndsMC tmpReadEnds = null;

    public MarkQueue(final ScoringStrategy duplicateScoringStrategy) {
        comparator = new ReadEndsMCComparator(duplicateScoringStrategy);
    }

    public int getNumDuplicates() { return this.numDuplicates; }

    public int size() {
        return this.set.size();
    }

    public boolean isEmpty() {
        return this.set.isEmpty();
    }

    public void setToMarkQueueMinimumDistance(final int toMarkQueueMinimumDistance) {
        this.toMarkQueueMinimumDistance = toMarkQueueMinimumDistance;
    }

    public int getToMarkQueueMinimumDistance() { return this.toMarkQueueMinimumDistance; }

    public boolean shouldBeInLocations(final ReadEndsMC current) {
        return (current.isPaired() && 0 == current.hasUnmapped);
    }

    /** For tracking optical duplicates. */
    public Set<ReadEnds> getLocations(final ReadEndsMC current) {
        // NB: only needed for pairs!!!
        if (!shouldBeInLocations(current)) throw new NotImplementedException();
        final Set<ReadEnds> locationSet = this.locations.remove(current).getReadEnds();
        if (null == locationSet) throw new PicardException("Locations was empty: unexpected error");
        return locationSet;
    }

    public ReadEndsMC peek() {
        return this.set.first();
    }

    private void updateMetrics(final ReadEndsMC duplicate, final DuplicationMetrics metrics) {
        // Update the duplication metrics
        if (!duplicate.getRecord().getReadPairedFlag() || duplicate.getRecord().getMateUnmappedFlag()) {
            ++metrics.UNPAIRED_READ_DUPLICATES;
        }
        else {
            ++metrics.READ_PAIR_DUPLICATES;// will need to be divided by 2 at the end
        }
        this.numDuplicates++;
    }

    /**
     * The poll method will return the read end that is *not* the duplicate of all comparable read ends that
     * have been seen.  All comparable read ends and the returned read end will have their seen duplicate flag set.
     *
     * NB: we must remove all fragments or unpaireds if this is a mapped pair
     * NB: we must remove all fragments if this is an unpaired
     */
    public ReadEndsMC poll(final DuplicateMarkingBuffer alignmentStartSortedBuffer,
                           final SAMFileHeader header,
                           final OpticalDuplicateFinder opticalDuplicateFinder,
                           final LibraryIdGenerator libraryIdGenerator) {
        final ReadEndsMC current = this.set.pollFirst();

        // Remove this record's comparable pair, if present.
        if (this.pairSet.contains(current)) { // the pair of this end is not a duplicate, if found
            final ReadEndsMC pair = this.pairSet.subSet(current, true, current, true).first();
            alignmentStartSortedBuffer.setDuplicateMarkingFlags(pair.getsamRecordIndex(), false);
            this.pairSet.remove(current);
            // NB: do not need to update metrics since this record is not a duplicate
        }

        // If we are a paired read end, we need to make sure we remove unpaired (if we are not also unpaired), as
        // well as fragments from the set, as they should all be duplicates.
        if (current.isPaired()) {

            // NB: only care about read1Sequence, read1Coordinate, and orientation in the set
            if (null == this.tmpReadEnds) { // initialize
                this.tmpReadEnds = new ReadEndsMC(header, current.getsamRecordIndex(), opticalDuplicateFinder, current.libraryId);
                this.tmpReadEnds.read2Sequence = this.tmpReadEnds.read2Coordinate = -1;
                this.tmpReadEnds.samRecordIndex = null;
            }
            else {
                this.tmpReadEnds.read1Sequence = current.read1Sequence;
                this.tmpReadEnds.read1Coordinate = current.read1Coordinate;
            }

            // We should search for one of F/R
            if (current.orientation == ReadEnds.FF || current.orientation == ReadEnds.FR || current.orientation == ReadEnds.F) {
                this.tmpReadEnds.orientation = ReadEnds.F;
            }
            else {
                this.tmpReadEnds.orientation = ReadEnds.R;
            }

            // remove from the set fragments and unpaired, which only have two possible orientations
            //this.tmpReadEnds.orientation = orientation;
            if (this.set.contains(this.tmpReadEnds)) { // found in the set
                // get the duplicate read end
                final SortedSet<ReadEndsMC> sortedSet = this.set.subSet(this.tmpReadEnds, true, this.tmpReadEnds, true);
                if (1 != sortedSet.size()) throw new PicardException("SortedSet should have size one (has size " + sortedSet.size() + " )");
                final ReadEndsMC duplicate = sortedSet.first();

                /** mark as duplicate and set that it has been through duplicate marking
                 * duplicate.getRecord().setDuplicateReadFlag(true); HANDLED BY THE METHOD CALL BELOW*/
                alignmentStartSortedBuffer.setDuplicateMarkingFlags(duplicate.getsamRecordIndex(), true);

                // remove from the set
                this.set.remove(this.tmpReadEnds);

                // update the metrics
                updateMetrics(duplicate, libraryIdGenerator.getMetricsByLibrary(libraryIdGenerator.getLibraryName(header, duplicate.getRecord())));
            }
        }

        // this read end is now ok to be emitted. track that it has been through duplicate marking
        alignmentStartSortedBuffer.setDuplicateMarkingFlags(current.getsamRecordIndex(), false);

        return current;
    }

    /**
     * Add a record to the mark queue.
     */
    private boolean debug = false; // TODO: remove me
    public void add(final ReadEndsMC other,
                    final DuplicateMarkingBuffer alignmentStartSortedBuffer,
                    final DuplicationMetrics metrics) {
        LocationSet locationSet = null;
        boolean addToLocationSet = true; // only false if we have mates mapped to the same position
        ReadEndsMC duplicate = null;

        if (debug) System.err.print("TMQ add: " + other.getRecord().getSAMString());

        // 1. check the queue to see if there exists a comparable record at the location, if so compare and keep the best.
        // 2. add physical location info if paired
        if (this.set.contains(other)) { // the "other" record already in the set
            if (debug) System.err.print("TMQ add: other contains " + other.getRecord().getSAMString());

            // Get the subset of records that are comparable, which should be of size one
            final SortedSet<ReadEndsMC> sortedSet = this.set.subSet(other, true, other, true);
            if (1 != sortedSet.size()) throw new PicardException("SortedSet should have size one (has size " + sortedSet.size() + " )");
            final ReadEndsMC current = sortedSet.first();

            // when checking for a read's pair, must check that read group (if present) AND read name match
            final String otherName = SAMUtils.getCanonicalRecordName(other.getRecord());
            final String currentName = SAMUtils.getCanonicalRecordName(current.getRecord());

            // check for the read's pair
            if (currentName.equals(otherName)) { // "other" is paired-end mate of "current". We need to choose the best end to store in the main set.
                if (debug) System.err.print("TMQ add: other contains : name equals " + other.getRecord().getSAMString());
                final int comparison = this.comparator.compare(current, other);
                // TODO: should this ever happen?  We should probably throw an exception in this case.
                if (0 < comparison) { // other is the best end. Swap for current.
                    if (debug) System.err.print("TMQ add: other contains : name equals : swap : " + other.getRecord().getSAMString());
                    // Swap the set
                    this.set.remove(current);
                    this.set.add(other);
                    this.pairSet.add(current);
                    // Swap "current" and "other" in the locations
                    if (shouldBeInLocations(other)) {
                        locationSet = this.locations.remove(current);
                        locationSet.replace(current, other);
                        this.locations.put(other, locationSet);
                        addToLocationSet = false;
                    }
                } else { // other is less desirable. Store it in the pair set.
                    if (debug) System.err.print("TMQ add: other contains : name equals : no swap : " + other.getRecord().getSAMString());
                    this.pairSet.add(other);
                    if (shouldBeInLocations(current)) {
                        locationSet = this.locations.get(current);
                        addToLocationSet = false;
                    }
                }
                if (debug) System.err.print("TMQ add: other contains : name equals : DONE : " + other.getRecord().getSAMString());
            } else { // "other" is a unique record at the same location and must be compared against "current"
                if (debug) System.err.print("TMQ add: other contains : name not equals " + other.getRecord().getSAMString());
                final int comparison = this.comparator.compare(current, other); // if we are to re-add, then other should make this > 0

                if (0 < comparison) { // remove the current, and add the other in its place
                    if (debug) System.err.print("TMQ add: other contains : name not equals : 0<comparison " + other.getRecord().getSAMString());
                    if (shouldBeInLocations(current)) { // was this in the location set?
                        // NB we could also just check if locationSet == null after remove?
                        locationSet = this.locations.remove(current);
                    }
                    else {
                        locationSet = new LocationSet();
                    }
                    this.locations.put(other, locationSet);
                    // remove current and add the other
                    this.set.remove(current);
                    this.set.add(other);

                    // update the pair set in case current's pair is in that set
                    if (this.pairSet.contains(current)) {
                        final ReadEndsMC pair = this.pairSet.subSet(current, true, current, true).first();
                        this.pairSet.remove(current);
                        alignmentStartSortedBuffer.setDuplicateMarkingFlags(pair.getsamRecordIndex(), true); // track that this samRecordIndex has been through duplicate marking
                        updateMetrics(pair, metrics);
                    }

                    duplicate = current;
                }
                else { // keep the current record, and the "other" is now a duplicate
                    if (shouldBeInLocations(current)) { // Get the location set
                        locationSet = this.locations.get(current);
                    }
                    // NB: else is technically not needed, since if this was not paired and the other one was, we would enter here and add it later

                    // other is a duplicate :/
                    duplicate = other;
                }
            }
        } else { // 'other' ReadEndMC is not in the main set, thus the first record at this location. Store it for now.
            if (debug) System.err.print("TMQ add: other not contains " + other.getRecord().getSAMString());

            if (shouldBeInLocations(other)) {
                locationSet = new LocationSet();
                this.locations.put(other, locationSet);
            }
            this.set.add(other);
        }
        //System.err.println("\tSET SIZE=" + this.set.size());

        // add to the physical locations
        final SAMRecord record = other.getRecord();
        if (debug && record.getReadPairedFlag()) System.err.print("TESTING FOR LOCATION SET (" + record.getFirstOfPairFlag() + ") : " + record.getSAMString());
        if (record.getReadPairedFlag()
                && !record.getReadUnmappedFlag()
                && !record.getMateUnmappedFlag()
                && addToLocationSet) {
            if (null == locationSet) throw new PicardException("location set was null: " + record.getSAMString());
            locationSet.add(other);
            if (debug) System.err.print("TESTING FOR LOCATION SET OK (" + locationSet.size() + " " +
            other.getTile() + ":" + other.getX() + ":" + other.getY() + ") : " + record.getSAMString());
            if (debug) {
                for (ReadEnds loc : locationSet.getReadEnds()) {
                    System.err.println("LOCATION SET loc: + " + loc.getX() + ":" + loc.getY());
                }
            }
        }

        // NB: locationSet can be empty, presumably when the second end of a pair is added first
        //if (null != locationSet && locationSet.isEmpty()) throw new PicardException("location set was unexpectedly empty");

        if (null != duplicate) {
            alignmentStartSortedBuffer.setDuplicateMarkingFlags(duplicate.getsamRecordIndex(), true);
            // count the duplicate metrics
            updateMetrics(duplicate, metrics);
        }
    }

    protected class LocationSet {
        private final Set<ReadEnds> readEnds = new HashSet<ReadEnds>();
        private final Set<PhysicalLocationMC> physicalLocations = new HashSet<PhysicalLocationMC>();

        // PhysicalLocationMC
        public LocationSet() {}

        public void add(final ReadEndsMC end) {
            final PhysicalLocationMC location = new PhysicalLocationMC(end);
            if (!physicalLocations.contains(location)) {
                readEnds.add(end);
                physicalLocations.add(new PhysicalLocationMC(location));
            }
        }

        public int size() { return physicalLocations.size(); }

        public void remove(final ReadEndsMC end) {
            final PhysicalLocationMC location = new PhysicalLocationMC(end);
            if (physicalLocations.contains(location)) {
                readEnds.remove(end);
                physicalLocations.remove(location);
            }
        }

        public Set<ReadEnds> getReadEnds() { return this.readEnds; }

        public void replace(final ReadEndsMC current, final ReadEndsMC other) {
            final PhysicalLocationMC location = new PhysicalLocationMC(current);
            if (!physicalLocations.contains(location)) {
                throw new PicardException("Trying to replace something not in the set");
            }
            this.remove(current);
            this.add(other);
        }
    }
}
