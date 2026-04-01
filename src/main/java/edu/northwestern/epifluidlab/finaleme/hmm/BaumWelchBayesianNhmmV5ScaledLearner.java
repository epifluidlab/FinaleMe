package edu.northwestern.epifluidlab.finaleme.hmm;

import java.util.*;
import java.util.concurrent.*;

import org.apache.commons.math3.util.Pair;

import edu.northwestern.epifluidlab.finaleme.utils.CcInferenceUtils;
import be.ac.ulg.montefiore.run.jahmm.*;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;


/**
 * An implementation of the Baum-Welch learning algorithm.  It uses a
 * scaling mechanism so as to avoid underflows.
 * <p>
 * For more information on the scaling procedure, read <i>Rabiner</i> and 
 * <i>Juang</i>'s <i>Fundamentals of speech recognition</i> (Prentice Hall,
 * 1993).
 */
public class BaumWelchBayesianNhmmV5ScaledLearner
extends BaumWelchScaledLearner
{	
	/**
	 * Number of iterations performed by the {@link #learn} method.
	 */
	private int nbIterations = 9;
	private int threadCount = -1;
	
	/**
	 * Initializes a Baum-Welch algorithm implementation.
	 */
	public BaumWelchBayesianNhmmV5ScaledLearner()
	{
	}

	public BaumWelchBayesianNhmmV5ScaledLearner(int threadCount)
	{
		this.threadCount = threadCount;
	}

	public void setThreadCount(int threadCount)
	{
		this.threadCount = threadCount;
	}

	private int resolveThreadCount()
	{
		int resolved = threadCount > 0 ? threadCount : Runtime.getRuntime().availableProcessors();
		return Math.max(1, resolved);
	}
	
	/**
	 * Performs one iteration of the Baum-Welch algorithm for non homogenous HMM.
	 * In one iteration, a new non homogenous HMM is computed using a previously estimated
	 * non homogenous HMM.
	 *
	 * @param hmm A previously estimated HMM.
	 * @param sequences The observation sequences on which the learning is
	 *         based.  Each sequence must have a length higher or equal to
         *         2.
	 * @return A new, updated HMM.
	 */
	public <O extends Observation> BayesianNhmmV5<O>
	iterate(BayesianNhmmV5<O> hmm, List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>>> sequences)
	{		
		BayesianNhmmV5<O> nhmm;
		try {
			nhmm = hmm.clone();
		} catch(CloneNotSupportedException e) {
			throw new InternalError();
		}
			
		/* a[i][j] = aijNum[i][j] / aijDen[i]
		 * aijDen[i] = expected number of transitions from state i
		 * aijNum[i][j] = expected number of transitions from state i to j
		 */
		double aijNum[][] = new double[hmm.nbStates()][hmm.nbStates()];
		double aijDen[] = new double[hmm.nbStates()];
		HashMap<Integer,Pair<double[][], double[]>> arij = new HashMap<Integer,Pair<double[][], double[]>>();

			Arrays.fill(aijDen, 0.);

			for (int j = 0; j < hmm.nbStates(); j++)
				Arrays.fill(aijNum[j], 0.);

		// Phase 1: Parallel forward-backward and xi computation
		int nThreads = resolveThreadCount();
		ExecutorService executor = Executors.newFixedThreadPool(nThreads);

		// Each task computes xi (gamma is derived from xi on-the-fly during accumulation)
		List<Future<Object[]>> futures = new ArrayList<Future<Object[]>>(sequences.size());
		for (int idx = 0; idx < sequences.size(); idx++) {
			final int seqIdx = idx;
			final Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> obsSeqPair = sequences.get(seqIdx);
			futures.add(executor.submit(new Callable<Object[]>() {
				@Override
				public Object[] call() {
					ForwardBackwardBayesianNhmmV5ScaledCalculator fbc =
						generateForwardBackwardCalculator(obsSeqPair, hmm);
					double[][][] xi = estimateXi(obsSeqPair, fbc, hmm);
					return new Object[]{xi, seqIdx};
				}
			}));
		}

		// Collect parallel results
		double[][][][] allXi = new double[sequences.size()][][][];
		try {
			for (Future<Object[]> future : futures) {
				Object[] result = future.get();
				double[][][] xi = (double[][][]) result[0];
				int seqIdx = (Integer) result[1];
				allXi[seqIdx] = xi;
			}
		} catch (Exception e) {
			throw new RuntimeException("Error in parallel forward-backward computation", e);
		}
		executor.shutdown();

		// Pre-compute total observations for PDF weight arrays
		int totalObs = 0;
		for (Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> seq : sequences)
			totalObs += seq.getSecond().size();

		// Memory-efficient accumulators: firstGamma for pi, pdfWeights for opdf fitting
		double[][] firstGamma = new double[sequences.size()][hmm.nbStates()];
		double[][] pdfWeights = new double[hmm.nbStates()][totalObs];
		double[] pdfWeightSums = new double[hmm.nbStates()];

		// Sequential accumulation into aijNum, aijDen, arij, firstGamma, pdfWeights
		int obsOffset = 0;
		for (int g = 0; g < sequences.size(); g++) {
			double[][][] xi = allXi[g];
			Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> obsSeqPair = sequences.get(g);
			List<? extends O> obsSeq = obsSeqPair.getSecond();
			HashMap<Integer, Pair<Integer, Double>> cpgDistState = obsSeqPair.getFirst();

			// Recompute gamma from xi on-the-fly (avoids storing allGamma)
			double[][] gamma = estimateGamma(xi, (ForwardBackwardBayesianNhmmV5Calculator) null);

			// Store first-timestep gamma for pi computation
			System.arraycopy(gamma[0], 0, firstGamma[g], 0, hmm.nbStates());

			// Accumulate PDF weights
			for (int i = 0; i < hmm.nbStates(); i++) {
				for (int t = 0; t < obsSeq.size(); t++) {
					pdfWeights[i][obsOffset + t] = gamma[t][i];
					pdfWeightSums[i] += gamma[t][i];
					if(Double.isNaN(gamma[t][i]) || Double.isInfinite(gamma[t][i])){
						System.err.println(gamma[t][i] + "\t" + pdfWeightSums[i] + "\t" + g + "\t" + t + "\t" + i);
						System.exit(1);
					}
				}
			}
			obsOffset += obsSeq.size();

			for (int i = 0; i < hmm.nbStates(); i++)
				for (int t = 0; t < obsSeq.size() - 1; t++) {
						aijDen[i] += gamma[t][i];
						if(arij.containsKey(cpgDistState.get(t+1).getFirst())){
							Pair<double[][], double[]> tmp = arij.get(cpgDistState.get(t+1).getFirst());
							double[] denTmp = tmp.getSecond();
							double[][] numTmp = tmp.getFirst();
							denTmp[i] += gamma[t][i];
							for (int j = 0; j < hmm.nbStates(); j++)
								numTmp[i][j] += xi[t][i][j];
							arij.put(cpgDistState.get(t+1).getFirst(), new Pair<double[][], double[]>(numTmp, denTmp));

						}else{
							double numTmp[][] = new double[hmm.nbStates()][hmm.nbStates()];
							double denTmp[] = new double[hmm.nbStates()];
							Arrays.fill(denTmp, 0.);
							for (int j = 0; j < hmm.nbStates(); j++)
								Arrays.fill(numTmp[j], 0.);
							denTmp[i] += gamma[t][i];
							for (int j = 0; j < hmm.nbStates(); j++)
								numTmp[i][j] += xi[t][i][j];
							arij.put(cpgDistState.get(t+1).getFirst(), new Pair<double[][], double[]>(numTmp, denTmp));
						}
						for (int j = 0; j < hmm.nbStates(); j++)
							aijNum[i][j] += xi[t][i][j];

				}
		}
		allXi = null; // allow GC

		for (int r = 0; r <= hmm.nbCpgDistState(); r++) {
			if(arij.containsKey(r)){
				Pair<double[][], double[]> tmp = arij.get(r);
				double[] denTmp = tmp.getSecond();
				double[][] numTmp = tmp.getFirst();

				for (int i = 0; i < hmm.nbStates(); i++) {
					if (denTmp[i] == 0.){
						for (int j = 0; j < hmm.nbStates(); j++)
							nhmm.setArij(r, i, j, hmm.getArij(r, i, j));
					}else{
						for (int j = 0; j < hmm.nbStates(); j++)
							nhmm.setArij(r, i, j, numTmp[i][j] / denTmp[i]);

					}

				}
			}else{
				for (int i = 0; i < hmm.nbStates(); i++) {
					if (aijDen[i] == 0.) // State i is not reachable
						for (int j = 0; j < hmm.nbStates(); j++)
							nhmm.setArij(r, i, j, hmm.getArij(r, i, j));
					else
						for (int j = 0; j < hmm.nbStates(); j++)
							nhmm.setArij(r, i, j, aijNum[i][j] / aijDen[i]);
				}
			}
		}


		/* pi computation */
		for (int r = 1; r < hmm.nbCpgDistState(); r++) {
			for (int i = 0; i < hmm.nbStates(); i++)
				nhmm.setPri(r, i, 0.);
		}

		for (int o = 0; o < sequences.size(); o++){
			Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> obsSeqPair = sequences.get(o);
			Integer r = obsSeqPair.getFirst().get(0).getFirst();
				for (int i = 0; i < hmm.nbStates(); i++){
					nhmm.setPri(r, i,
								nhmm.getPri(r, i) + firstGamma[o][i] / sequences.size());
				}
		}
		firstGamma = null; // allow GC

		/* rescale pi */
		HashMap<Integer, Double> sumPi = new HashMap<Integer, Double>();

		for (int r = 0; r <= hmm.nbCpgDistState(); r++) {
			for (int i = 0; i < hmm.nbStates(); i++){
				if(sumPi.containsKey(r)){
					sumPi.put(r, (sumPi.get(r) + nhmm.getPri(r, i)));
				}else{
					sumPi.put(r, nhmm.getPri(r, i));
				}

			}
		}

		for (int r = 0; r <= hmm.nbCpgDistState(); r++) {
			for (int i = 0; i < hmm.nbStates(); i++){
				if(sumPi.get(r) == 0){
					nhmm.setPri(r, i, 0.5);
				}else{
					nhmm.setPri(r, i, nhmm.getPri(r, i)/sumPi.get(r));
				}

			}
		}
		/* pdfs computation -- weights already accumulated during sequential phase */
		List<O> observations = CcInferenceUtils.flatPair(sequences);
		for (int i = 0; i < hmm.nbStates(); i++) {
			double sum = pdfWeightSums[i];
			double[] weights = pdfWeights[i];

			for (int j = weights.length - 1; j >= 0; j--){
				weights[j] /= sum;
				if(Double.isNaN(weights[j]) || Double.isInfinite(weights[j])){
					System.err.println(weights[j] + "\t" + sum + "\t" +  j);
					System.exit(1);
				}
			}

			Opdf<O> opdf = nhmm.getOpdf(i);
			opdf.fit(observations, weights);
		}

		return nhmm;
	}

	/**
	 * Does a fixed number of iterations (see {@link #getNbIterations}) of the
	 * Baum-Welch algorithm.
	 * 
	 * @param initialNhmm An initial estimation of the expected HMM.  This
	 *         estimate is critical as the Baum-Welch algorithm only find
	 *         local minima of its likelihood function.
	 * @param sequences The observation sequences on which the learning is
	 *         based.  Each sequence must have a length higher or equal to 2.
	 * @return The HMM that best matches the set of observation sequences given
	 *         (according to the Baum-Welch algorithm).
	 */
	public <O extends Observation> BayesianNhmmV5<O>
	learn(BayesianNhmmV5<O> initialNhmm, List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>>> sequences)
	{
		BayesianNhmmV5<O> hmm = initialNhmm;
		
		for (int i = 0; i < nbIterations; i++)
			hmm = iterate(hmm, sequences);
		
		return hmm;
	}
	
	protected <O extends Observation> ForwardBackwardBayesianNhmmV5ScaledCalculator
	generateForwardBackwardCalculator(Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> sequence,
			BayesianNhmmV5<O> hmm)
	{
		return new ForwardBackwardBayesianNhmmV5ScaledCalculator(sequence, hmm, 
				EnumSet.allOf(ForwardBackwardBayesianNhmmV5ScaledCalculator.Computation.class));
	}
	
	
	/* Here, the xi (and, thus, gamma) values are not divided by the
	 probability of the sequence because this probability might be
	 too small and induce an underflow. xi[t][i][j] still can be
	 interpreted as P[q_t = i and q_(t+1) = j | obsSeq, hmm] because
	 we assume that the scaling factors are such that their product
	 is equal to the inverse of the probability of the sequence. */
	protected <O extends Observation> double[][][]
	estimateXi(Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> obsSeqPair, ForwardBackwardBayesianNhmmV5Calculator fbc,
			BayesianNhmmV5<O> hmm)
	{	
		List<? extends O> sequence = obsSeqPair.getSecond();
		HashMap<Integer, Pair<Integer, Double>> cpgDistState = obsSeqPair.getFirst();
		
		if (sequence.size() <= 1)
			throw new IllegalArgumentException("Observation sequence too " + 
			"short");
		
		double xi[][][] = 
			new double[sequence.size() - 1][hmm.nbStates()][hmm.nbStates()];
		
		Iterator<? extends O> seqIterator = sequence.iterator();
		seqIterator.next();
		
		for (int t = 0; t < sequence.size() - 1; t++) {
			O observation = seqIterator.next();
			
			
				
				if(cpgDistState.containsKey(t)){
					
					for (int i = 0; i < hmm.nbStates(); i++){
						for (int j = 0; j < hmm.nbStates(); j++){
							
							xi[t][i][j] = fbc.alphaElement(t, i) *
									hmm.getArij(cpgDistState.get(t+1).getFirst(), i, j) *
									hmm.getOpdfProb(j,observation) *
									fbc.betaElement(t + 1, j);
							if(Double.isNaN(xi[t][i][j])){
								System.err.println(t + "\t" + i + "\t" + j + "\t" + cpgDistState.get(t) + "\t" + observation + "\t" +  hmm.getArij(cpgDistState.get(t).getFirst(), i, j));
								System.err.println(fbc.alphaElement(t, i) + "\t" + hmm.getArij(cpgDistState.get(t).getFirst(), i, j) + "\t" + hmm.getOpdf(j).probability(observation) + "\t" + fbc.betaElement(t + 1, j));
							}
						}
							
					}
						
				}else{
					for (int i = 0; i < hmm.nbStates(); i++)
						for (int j = 0; j < hmm.nbStates(); j++)
							xi[t][i][j] = 0.5;
				}
				
		}
		
		return xi;
	}
	
	/* gamma[][] could be computed directly using the alpha and betacd 
	 * arrays, but this (slower) method is prefered because it doesn't
	 * change if the xi array has been scaled (and should be changed with
	 * the scaled alpha and beta arrays).
	 */
	protected double[][]
	estimateGamma(double[][][] xi, ForwardBackwardBayesianNhmmV5Calculator fbc)
	{
		double[][] gamma = new double[xi.length + 1][xi[0].length];
		
		for (int t = 0; t < xi.length + 1; t++)
				Arrays.fill(gamma[t], 0.);
		
		for (int t = 0; t < xi.length; t++)
				for (int i = 0; i < xi[0][0].length; i++)
					for (int j = 0; j < xi[0][0].length; j++)
						gamma[t][i] += xi[t][i][j];
		
		
			for (int j = 0; j < xi[0][0].length; j++)
				for (int i = 0; i < xi[0][0].length; i++)
					gamma[xi.length][j] += xi[xi.length - 1][i][j];
		
		return gamma;
	}
	
}
