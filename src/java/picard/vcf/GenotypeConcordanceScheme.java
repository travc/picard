package picard.vcf;

import picard.PicardException;
import picard.vcf.GenotypeConcordanceStates.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * This defines for each valid TruthState and CallState tuple, the set of contingency table entries that to which the tuple should contribute.
 * @author nhomer
 */
public class GenotypeConcordanceScheme {

    /** The underling scheme */
    protected final Map<TruthAndCallStates, Set<ContingencyState>> scheme = new HashMap<TruthAndCallStates, Set<ContingencyState>>();

    /** These are convenience variables for defining a scheme.  NA means that such a tuple should never be observed. */
    public static final Set<ContingencyState>    NA       = new HashSet<ContingencyState>() {{ add(ContingencyState.NA); }};
    protected static final Set<ContingencyState> EMPTY    = new HashSet<ContingencyState>();
    protected static final Set<ContingencyState> TP_ONLY  = new HashSet<ContingencyState>() {{ add(ContingencyState.TP); }};
    protected static final Set<ContingencyState> FP_ONLY  = new HashSet<ContingencyState>() {{ add(ContingencyState.FP); }};
    protected static final Set<ContingencyState> TN_ONLY  = new HashSet<ContingencyState>() {{ add(ContingencyState.TN); }};
    protected static final Set<ContingencyState> FN_ONLY  = new HashSet<ContingencyState>() {{ add(ContingencyState.FN); }};
    protected static final Set<ContingencyState> TP_FN    = new HashSet<ContingencyState>() {{ add(ContingencyState.TP); add(ContingencyState.FN); }};
    protected static final Set<ContingencyState> TP_FP    = new HashSet<ContingencyState>() {{ add(ContingencyState.TP); add(ContingencyState.FP); }};
    protected static final Set<ContingencyState> FP_FN    = new HashSet<ContingencyState>() {{ add(ContingencyState.FP); add(ContingencyState.FN); }};
    protected static final Set<ContingencyState> TP_FP_FN = new HashSet<ContingencyState>() {{ add(ContingencyState.TP); add(ContingencyState.FP); add(ContingencyState.FN);}};

    /** Has this boolean been previously validated */
    private boolean isValidated = false;

    /**
     * Adds a row to the scheme
     * @param callState the call state (row)
     * @param concordanceStateSets the concordance state sets for each truth value, in order
     */
    protected void addRow(final CallState callState, final Set<ContingencyState>... concordanceStateSets) {
        if (concordanceStateSets.length != TruthState.values().length) {
            throw new PicardException("Length mismatch between concordanceStateSets and TruthState.values()");
        }
        for (int i = 0; i < concordanceStateSets.length; i++) {
            scheme.put(new TruthAndCallStates(TruthState.values()[i], callState), concordanceStateSets[i]);
            i++;
        }
    }

    /**
     * The scheme is defined in the constructor.
     */
    public GenotypeConcordanceScheme() {

        /**          ROW STATE            HOM_REF       HET_REF_VAR1       HET_VAR1_VAR2        HOM_VAR1        NO_CALL        LOW_GQ        LOW_DP        IS_MIXED    **/
        addRow(CallState.HOM_REF,         TN_ONLY,      FN_ONLY,           FN_ONLY,             FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HET_REF_VAR1,    FP_ONLY,      TP_ONLY,           TP_FN,               TP_FN,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HET_REF_VAR2,    NA,           TP_FP_FN,          TP_FN,               FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HET_REF_VAR3,    NA,           NA,                FP_FN,               NA,             EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HET_VAR1_VAR2,   FP_ONLY,      TP_FP,             TP_ONLY,             TP_FP_FN,       EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HET_VAR1_VAR3,   NA,           NA,                TP_FP_FN,            NA,             EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HET_VAR3_VAR4,   FP_ONLY,      FP_FN,             FP_FN,               FP_FN,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HOM_VAR1,        FP_ONLY,      TP_FP,             TP_FN,               TP_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HOM_VAR2,        NA,           FP_FN,             TP_FN,               FP_FN,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HOM_VAR3,        NA,           NA,                FP_FN,               NA,             EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.NO_CALL,         EMPTY,        EMPTY,             EMPTY,               EMPTY,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.FILTERED,        EMPTY,        EMPTY,             EMPTY,               EMPTY,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.LOW_GQ,          EMPTY,        EMPTY,             EMPTY,               EMPTY,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.LOW_DP,          EMPTY,        EMPTY,             EMPTY,               EMPTY,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.IS_MIXED,        EMPTY,        EMPTY,             EMPTY,               EMPTY,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY);

        validateScheme();
    }{}

    /**
     * Get the concordance state set associate with the given truth state and call state tuple.
     */
    public Set<ContingencyState> getConcordanceStateSet(final TruthState truthState, final CallState callState) {
        return this.getConcordanceStateSet(new TruthAndCallStates(truthState, callState));
    }

    /**
     * Get the concordance state set associate with the given truth state and call state tuple.
     */
    public Set<ContingencyState> getConcordanceStateSet(final TruthAndCallStates truthAndCallStates) {
        return this.scheme.get(truthAndCallStates);
    }

    /**
     * Check that all cells in the scheme exist.
     * @throws PicardException if a missing tuple was found.
     */
    public void validateScheme() throws PicardException {
        if (!isValidated) {
            for (final TruthState truthState : TruthState.values()) {
                for (final CallState callState : CallState.values()) {
                    if (!scheme.containsKey(new TruthAndCallStates(truthState, callState))) {
                        throw new PicardException(String.format("Missing scheme tuple: [%s, %s]", truthState.name(), callState.name()));
                    }
                }
            }
        }

        isValidated = true;
    }
}
