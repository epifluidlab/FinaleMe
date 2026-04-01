/**
 * KullbackLeiblerDistanceBayesianNhmmV5Calculator.java
 * Apr 28, 2016
 * 11:16:57 AM
 * yaping    lyping1986@gmail.com
 */
package edu.northwestern.epifluidlab.finaleme.hmm;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.*;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.Pair;

import be.ac.ulg.montefiore.run.jahmm.Observation;


/**
 *
 */
public class KullbackLeiblerDistanceBayesianNhmmV5Calculator<O extends Observation>{

	private int sequencesLength = 10;
	private int nbSequences = 10000;
	private int threadCount = -1;
	//private ArrayList<Pair<Integer, Double>> cpgDistFreq; //maybe just fit cpg distance distribution with gaussian/poisson mixture model for the simplicity first?
	private List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>>> matrix;
	private MersenneTwister randomEngine;
	public KullbackLeiblerDistanceBayesianNhmmV5Calculator(List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>>> matrix, MersenneTwister randomEngine){ //cpg distance is a poisson distribution
		this.matrix = matrix;
		this.randomEngine = randomEngine;
		this.nbSequences = Math.max(this.nbSequences, matrix.size()/100);
	}
	
	public KullbackLeiblerDistanceBayesianNhmmV5Calculator(List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>>> matrix){ //cpg distance is a poisson distribution
		this.matrix = matrix;
	}

	public KullbackLeiblerDistanceBayesianNhmmV5Calculator(List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>>> matrix, int threadCount){ //cpg distance is a poisson distribution
		this.matrix = matrix;
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
	 * Computes the Kullback-Leibler distance between two HMMs.
	 *
	 * @param hmm1 The first HMM against which the distance is computed.
	 *             The distance is mesured with regard to this HMM (this must
	 *             be defined since the Kullback-Leibler distance is not
	 *             symetric).
	 * @param hmm2 The second HMM against which the distance is computed.
	 * @return The distance between <code>hmm1</code> and <code>hmm2</code> with
	 *      regard to <code>hmm1</code>
	 */
	public double 
	distance(BayesianNhmmV5<O> hmm1, BayesianNhmmV5<O> hmm2)
	{			
		double distance = 0.;
		//int tried = 0;
		for (int i = 0; i < nbSequences; i++) {
			//Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair = new BayesianNhmmV5MarkovGenerator<O>(hmm1, cpgDistFreq).
			//observationSequence(sequencesLength);
			Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair = matrix.get(randomEngine.nextInt(matrix.size()));
			while(oseqPair.getSecond().size()<sequencesLength){
				oseqPair = matrix.get(randomEngine.nextInt(matrix.size()));
			}
			//for(int s = 0; s < oseqPair.getSecond().size(); s++){
			//	System.err.print(oseqPair.getSecond().get(s));
			//}
			//System.err.println();
			//for(int s = 0; s < oseqPair.getSecond().size(); s++){
			//	System.err.print(oseqPair.getFirst().get(s) + "\t");
			//}
			//System.err.println();
			double prob1 = new ForwardBackwardBayesianNhmmV5ScaledCalculator(oseqPair, hmm1).
					lnProbability();
			double prob2 = new ForwardBackwardBayesianNhmmV5ScaledCalculator(oseqPair, hmm2).
					lnProbability();
			//if(Double.isNaN(prob1) || Double.isNaN(prob2) || Double.isInfinite(prob1) || Double.isInfinite(prob2)){
			//	i--;
			//	tried++;
			//	if(tried > 20){
			//		distance +=  (prob1 - prob2) / sequencesLength;
			//		i++;
			//	}
			//	continue;
			//}else{
				distance +=  (prob1 - prob2) / oseqPair.getSecond().size();
			//	tried = 0;
			//}
			
			
			//System.err.println(distance + "\t" + new ForwardBackwardBayesianNhmmV5ScaledCalculator(oseqPair, hmm1).
			//		lnProbability() + "\t" + new ForwardBackwardBayesianNhmmV5ScaledCalculator(oseqPair, hmm2).
			//		lnProbability() + "\t" + sequencesLength);
			//if(Double.isNaN(distance)){
			//	break;
			//}
			
		}
		
		return distance / nbSequences;
	}

	public double
	distance(BayesianNhmmV5<O> hmm1, BayesianNhmmV5<O> hmm2, boolean allSites)
	{
		// Filter sequences that meet minimum length requirement
		List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>>> validSeqs =
			new ArrayList<Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>>>();
		for (Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair : matrix) {
			if (oseqPair.getSecond().size() >= sequencesLength) {
				validSeqs.add(oseqPair);
			}
		}

		if (validSeqs.isEmpty()) {
			return 0.;
		}

		// Parallel forward-backward computation across all valid sequences
		int nThreads = resolveThreadCount();
		ExecutorService executor = Executors.newFixedThreadPool(nThreads);
		List<Future<Double>> futures = new ArrayList<Future<Double>>(validSeqs.size());

		for (final Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair : validSeqs) {
			futures.add(executor.submit(new Callable<Double>() {
				@Override
				public Double call() {
					double prob1 = new ForwardBackwardBayesianNhmmV5ScaledCalculator(oseqPair, hmm1).
							lnProbability();
					double prob2 = new ForwardBackwardBayesianNhmmV5ScaledCalculator(oseqPair, hmm2).
							lnProbability();
					return (prob1 - prob2) / oseqPair.getSecond().size();
				}
			}));
		}

		double distance = 0.;
		try {
			for (Future<Double> future : futures) {
				distance += future.get();
			}
		} catch (Exception e) {
			throw new RuntimeException("Error in parallel KL-divergence computation", e);
		}
		executor.shutdown();

		return distance / validSeqs.size();
	}
	
	

	/**
	 * Returns the number of sequences generated to estimate a distance.
	 * 
	 * @return The number of generated sequences.
	 */
	public int getNbSequences()
	{
		return nbSequences;
	}


	/**
	 * Sets the number of sequences generated to estimate a distance.
	 * 
	 * @param nb The number of generated sequences.
	 */
	public void setNbSequences(int nb)
	{
		this.nbSequences = nb;
	}
	
	
	/**
	 * Returns the length of sequences generated to estimate a distance.
	 * 
	 * @return The sequences length.
	 */
	public int getSequencesLength()
	{
		return sequencesLength;
	}


	/**
	 * Sets the length of sequences generated to estimate a distance.
	 * 
	 * @param length The sequences length.
	 */
	public void setSequencesLength(int length)
	{
		this.sequencesLength = length;
	}

}
