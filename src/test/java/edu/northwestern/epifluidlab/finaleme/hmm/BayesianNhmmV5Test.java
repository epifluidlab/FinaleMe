package edu.northwestern.epifluidlab.finaleme.hmm;

import static org.junit.jupiter.api.Assertions.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.math3.util.Pair;
import edu.northwestern.epifluidlab.finaleme.utils.ObservationVector;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

/**
 * Unit tests for BayesianNhmmV5 HMM model and related components.
 */
public class BayesianNhmmV5Test {

    private static final int NB_STATES = 2;
    private static final int NB_CPG_DIST_STATES = 10;
    private static final double BAYESIAN_FACTOR = 0.5;
    private static final double EPSILON = 1e-10;

    private BayesianNhmmV5<ObservationVector> hmm;

    @BeforeEach
    void setUp() {
        ArrayList<Integer> mixNumberInFeature = new ArrayList<>();
        mixNumberInFeature.add(1);
        mixNumberInFeature.add(1);
        mixNumberInFeature.add(1);

        OpdfMultiMixtureGaussianFactory factory =
            new OpdfMultiMixtureGaussianFactory(3, mixNumberInFeature,
                new double[]{0.0, 0.0, 0.0}, new double[]{1.0, 1.0, 1.0});

        hmm = new BayesianNhmmV5<>(NB_STATES, NB_CPG_DIST_STATES, factory, BAYESIAN_FACTOR);
    }

    @Test
    void testBasicProperties() {
        assertEquals(NB_STATES, hmm.nbStates());
        assertEquals(NB_CPG_DIST_STATES, hmm.nbCpgDistState());
        assertEquals(BAYESIAN_FACTOR, hmm.getBayesianFactor(), EPSILON);
    }

    @Test
    void testInitialPiValues() {
        // Initial pi should be uniform: 1/nbStates for all states and dist indices
        for (int r = 0; r <= NB_CPG_DIST_STATES; r++) {
            double sum = 0;
            for (int i = 0; i < NB_STATES; i++) {
                assertEquals(0.5, hmm.getPri(r, i), EPSILON,
                    "Pi should be 1/nbStates initially for r=" + r + " state=" + i);
                sum += hmm.getPri(r, i);
            }
            assertEquals(1.0, sum, EPSILON, "Pi values should sum to 1 for r=" + r);
        }
    }

    @Test
    void testInitialTransitionValues() {
        // Initial transitions should be uniform
        for (int r = 0; r <= NB_CPG_DIST_STATES; r++) {
            for (int i = 0; i < NB_STATES; i++) {
                double sum = 0;
                for (int j = 0; j < NB_STATES; j++) {
                    assertEquals(0.5, hmm.getArij(r, i, j), EPSILON);
                    sum += hmm.getArij(r, i, j);
                }
                assertEquals(1.0, sum, EPSILON, "Transition row should sum to 1");
            }
        }
    }

    @Test
    void testSetGetPi() {
        hmm.setPri(5, 0, 0.7);
        hmm.setPri(5, 1, 0.3);
        assertEquals(0.7, hmm.getPri(5, 0), EPSILON);
        assertEquals(0.3, hmm.getPri(5, 1), EPSILON);
    }

    @Test
    void testSetGetTransition() {
        hmm.setArij(3, 0, 0, 0.9);
        hmm.setArij(3, 0, 1, 0.1);
        assertEquals(0.9, hmm.getArij(3, 0, 0), EPSILON);
        assertEquals(0.1, hmm.getArij(3, 0, 1), EPSILON);
    }

    @Test
    void testClone() throws CloneNotSupportedException {
        hmm.setPri(2, 0, 0.8);
        hmm.setArij(1, 0, 1, 0.6);

        BayesianNhmmV5<ObservationVector> clone = hmm.clone();

        // Values should be equal
        assertEquals(0.8, clone.getPri(2, 0), EPSILON);
        assertEquals(0.6, clone.getArij(1, 0, 1), EPSILON);

        // Modifying clone should NOT affect original
        clone.setPri(2, 0, 0.1);
        assertEquals(0.8, hmm.getPri(2, 0), EPSILON);
        assertEquals(0.1, clone.getPri(2, 0), EPSILON);

        clone.setArij(1, 0, 1, 0.99);
        assertEquals(0.6, hmm.getArij(1, 0, 1), EPSILON);
    }

    @Test
    void testToString() {
        String s = hmm.toString();
        assertNotNull(s);
        assertTrue(s.contains("HMM with " + NB_STATES + " state(s)"));
        assertTrue(s.contains("State 0"));
        assertTrue(s.contains("State 1"));
        assertTrue(s.contains("Pi:"));
        assertTrue(s.contains("Aij:"));
    }

    @Test
    void testOpdfNotNull() {
        for (int i = 0; i < NB_STATES; i++) {
            assertNotNull(hmm.getOpdf(i), "Opdf for state " + i + " should not be null");
        }
    }

    @Test
    void testForwardBackwardProbability() {
        // Create a simple observation sequence with CpG distance state info
        Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>> obsSeqPair = createTestSequence(5);

        // Probability should be a finite positive number
        double lnProb = hmm.lnProbability(obsSeqPair);
        assertTrue(Double.isFinite(lnProb), "Log probability should be finite, got: " + lnProb);
    }

    @Test
    void testViterbiDecoding() {
        Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>> obsSeqPair = createTestSequence(5);

        int[] stateSeq = hmm.mostLikelyStateSequence(obsSeqPair, 1, 0.5);

        assertNotNull(stateSeq);
        assertEquals(5, stateSeq.length);
        for (int s : stateSeq) {
            assertTrue(s >= 0 && s < NB_STATES, "State should be valid: " + s);
        }
    }

    @Test
    void testViterbiConsistency() {
        // Create HMM with clear state separation
        ArrayList<Integer> mixNumberInFeature = new ArrayList<>();
        mixNumberInFeature.add(1);
        mixNumberInFeature.add(1);
        mixNumberInFeature.add(1);

        // State 0: low values, State 1: high values
        BayesianNhmmV5<ObservationVector> clearHmm = new BayesianNhmmV5<>(
            2, NB_CPG_DIST_STATES,
            new OpdfMultiMixtureGaussianFactory(3, mixNumberInFeature,
                new double[]{-2.0, -2.0, -2.0}, new double[]{0.1, 0.1, 0.1}),
            BAYESIAN_FACTOR);

        // Set strong self-transitions
        for (int r = 0; r <= NB_CPG_DIST_STATES; r++) {
            clearHmm.setArij(r, 0, 0, 0.9);
            clearHmm.setArij(r, 0, 1, 0.1);
            clearHmm.setArij(r, 1, 0, 0.1);
            clearHmm.setArij(r, 1, 1, 0.9);
        }

        Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>> obsSeqPair = createTestSequence(10);

        int[] stateSeq = clearHmm.mostLikelyStateSequence(obsSeqPair, 1, 0.5);
        assertNotNull(stateSeq);
        assertEquals(10, stateSeq.length);
    }

    @Test
    void testBaumWelchSingleIteration() {
        BaumWelchBayesianNhmmV5ScaledLearner learner = new BaumWelchBayesianNhmmV5ScaledLearner();

        List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> sequences = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            sequences.add(createTestSequence(5 + i));
        }

        BayesianNhmmV5<ObservationVector> updated = learner.iterate(hmm, sequences);

        assertNotNull(updated);
        assertEquals(NB_STATES, updated.nbStates());

        // Pi and transitions should still be valid probabilities
        for (int r = 0; r <= NB_CPG_DIST_STATES; r++) {
            for (int i = 0; i < NB_STATES; i++) {
                double pi = updated.getPri(r, i);
                assertTrue(pi >= 0 && pi <= 1.0 + EPSILON, "Pi should be in [0,1]: " + pi);

                double rowSum = 0;
                for (int j = 0; j < NB_STATES; j++) {
                    double aij = updated.getArij(r, i, j);
                    assertTrue(Double.isFinite(aij), "Transition should be finite");
                    rowSum += aij;
                }
                // Transition rows should approximately sum to 1
                if (rowSum > 0) {
                    assertEquals(1.0, rowSum, 0.01, "Transition row should sum to ~1 for r=" + r + " i=" + i);
                }
            }
        }
    }

    @Test
    void testKLDivergenceZeroForIdenticalModels() throws CloneNotSupportedException {
        BayesianNhmmV5<ObservationVector> hmm2 = hmm.clone();

        List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> matrix = new ArrayList<>();
        for (int i = 0; i < 20; i++) {
            matrix.add(createTestSequence(15));
        }

        KullbackLeiblerDistanceBayesianNhmmV5Calculator<ObservationVector> klCalc =
            new KullbackLeiblerDistanceBayesianNhmmV5Calculator<>(matrix);

        double dist = klCalc.distance(hmm, hmm2, true);
        assertEquals(0.0, dist, 0.01, "KL distance between identical models should be ~0");
    }

    /**
     * Creates a test observation sequence with CpG distance state information.
     */
    private Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>> createTestSequence(int length) {
        HashMap<Integer, Pair<Integer, Double>> cpgDistState = new HashMap<>();
        List<ObservationVector> observations = new ArrayList<>();

        for (int t = 0; t < length; t++) {
            int dist = Math.min(t + 1, NB_CPG_DIST_STATES);
            double methyPrior = 0.5;
            cpgDistState.put(t, new Pair<>(dist, methyPrior));
            observations.add(new ObservationVector(new double[]{
                Math.sin(t * 0.5),
                Math.cos(t * 0.3),
                0.5 + 0.1 * t
            }));
        }

        return new Pair<>(cpgDistState, observations);
    }
}
