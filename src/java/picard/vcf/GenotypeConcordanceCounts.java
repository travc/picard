package picard.vcf;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import picard.PicardException;
import picard.vcf.GenotypeConcordanceStates.*;

/**
 * A class to store the counts for various truth and call state classifications relative to a reference.  With these counts and a provided
 * scheme, summary metrics can be returned.
 * @author nhomer
 */
public class GenotypeConcordanceCounts {

    /**
     * Pre-defined sets based on if the caller wishes to return the sensitivity given the common homozygous reference, heterozygous, and homozygous variant cases.
     */
    public static final Set<TruthState> HOM_REF_TRUTH_STATE_SET = new HashSet<TruthState>() {{ add(TruthState.HOM_REF); }};
    public static final Set<TruthState> HET_TRUTH_STATE_SET     = new HashSet<TruthState>() {{ add(TruthState.HET_REF_VAR1); add(TruthState.HET_VAR1_VAR2);}};
    public static final Set<TruthState> HOM_VAR_TRUTH_STATE_SET = new HashSet<TruthState>() {{ add(TruthState.HOM_VAR1); }};

    /**
     * Pre-defined sets based on if the caller wishes to return the PPV given the common homozygous reference, heterozygous, and homozygous variant cases.
     */
    public static final Set<CallState> HOM_REF_CALL_STATE_SET = new HashSet<CallState>() {{ add(CallState.HOM_REF); }};
    public static final Set<CallState> HET_CALL_STATE_SET     = new HashSet<CallState>() {{ addAll(Arrays.asList(CallState.HET_REF_VAR1, CallState.HET_REF_VAR2, CallState.HET_REF_VAR3, CallState.HET_VAR1_VAR2, CallState.HET_VAR1_VAR3, CallState.HET_VAR3_VAR4)); }};
    public static final Set<CallState> HOM_VAR_CALL_STATE_SET = new HashSet<CallState>() {{ add(CallState.HOM_VAR1); add(CallState.HOM_VAR2); add(CallState.HOM_VAR3); }};


    /** The underling counts table */
    private final Map<TruthAndCallStates, Integer> counts = new HashMap<TruthAndCallStates, Integer>();

    /**
     * Creates a new counts table initializing the counts to zero.
     */
    public GenotypeConcordanceCounts() {
        // initialize the counts to zero
        for (final TruthState truthState : TruthState.values()) {
            for (final CallState callState : CallState.values()) {
                this.counts.put(new TruthAndCallStates(truthState, callState), 0);
            }
        }
    }

    /**
     * Increments a count for the truth/call state tuple.
     * @param truthState
     * @param callState
     */
    public void increment(final TruthState truthState, final CallState callState) {
        // TODO: so many object creations!
        final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
        this.counts.put(truthAndCallStates, this.counts.get(truthAndCallStates) + 1);
    }

    /**
     * Validates that there are no counts for NA states in the underlying scheme
     */
    public void validateCountsAgainstScheme(final GenotypeConcordanceScheme scheme) {
        for (final TruthState truthState : TruthState.values()) {
            for (final CallState callState : CallState.values()) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                if (0 < this.counts.get(truthAndCallStates) && scheme.getConcordanceStateSet(truthAndCallStates).containsAll(GenotypeConcordanceScheme.NA)) {
                    throw new PicardException(String.format("Found counts for an illegal set of states: [%s, %s]", truthState.name(), callState.name()));
                }
            }
        }
    }

    /**
     * Returns the sensitivity defined by the scheme across the subset of truth states.
     */
    public double getSensitivity(final GenotypeConcordanceScheme scheme, final Set<TruthState> truthStateSet) {
        /**
         * Sensitivity is the TP / P = TP / (TP + FN)
         */
        double numerator = 0.0;
        double denominator = 0.0;

        scheme.validateScheme();

        for (final TruthState truthState : truthStateSet) {
            for (final CallState callState : CallState.values()) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                int count = this.counts.get(truthAndCallStates);
                for (final ContingencyState contingencyState : scheme.getConcordanceStateSet(truthAndCallStates)) {
                    if (ContingencyState.TP == contingencyState) {
                        numerator += count;
                        denominator += count;
                    }
                    else if (ContingencyState.FN == contingencyState) {
                        denominator += count;
                    }
                }
            }
        }

        return (numerator / denominator);
    }

    /**
     * Returns the sensitivity defined by the scheme across the all truth states.
     */
    public double getSensitivity(final GenotypeConcordanceScheme scheme) {
        return getSensitivity(scheme, new HashSet<TruthState>(Arrays.asList(TruthState.values())));

    }

    /**
     * Returns the PPV defined by the scheme across the subset of truth states.
     */
    public double getPPV(final GenotypeConcordanceScheme scheme, final Set<CallState> callStateSet) {
        /**
         * PPV is the TP / (TP + FP)
         */
        double numerator = 0.0;
        double denominator = 0.0;

        scheme.validateScheme();

        for (final CallState callState : callStateSet) {
            for (final TruthState truthState : TruthState.values()) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                int count = this.counts.get(truthAndCallStates);
                for (final ContingencyState contingencyState : scheme.getConcordanceStateSet(truthAndCallStates)) {
                    if (ContingencyState.TP == contingencyState) {
                        numerator += count;
                        denominator += count;
                    }
                    else if (ContingencyState.FP == contingencyState) {
                        denominator += count;
                    }
                }
            }
        }

        return (numerator / denominator);
    }

    /**
     * Returns the PPV defined by the scheme across the all truth states.
     */
    public double getPPV(final GenotypeConcordanceScheme scheme) {
        return getPPV(scheme, new HashSet<CallState>(Arrays.asList(CallState.values())));
    }

    /**
     * Returns the count defined by the truth state set and call state set.
     */
    public int getCount(final TruthState truthState, final CallState callState) {
        return this.counts.get(new TruthAndCallStates(truthState, callState));
    }

    /**
     * Returns the sum of the all pairs of tuples defined by the truth state set and call state set.
     */
    public int getSum(final Set<TruthState> truthStateSet, final Set<CallState> callStateSet) {
        int count = 0;
        for (final TruthState truthState : truthStateSet) {
            for (final CallState callState : callStateSet) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                count += this.counts.get(truthAndCallStates);
            }
        }
        return count;
    }

    /**
     * Returns the sum of the all pairs of tuples defined by the truth state set and call state set.
     */
    public int getSum() {
        return getSum(new HashSet<TruthState>(Arrays.asList(TruthState.values())), new HashSet<CallState>(Arrays.asList(CallState.values())));
    }
}
