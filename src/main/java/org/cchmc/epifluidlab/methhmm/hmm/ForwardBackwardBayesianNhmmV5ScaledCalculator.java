/**
 * ForwardBackwardNhmmScaledCalculator.java
 * Apr 27, 2016
 * 1:31:23 PM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.methhmm.hmm;

import java.util.Arrays;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.util.Pair;

import java.math.BigDecimal;

import be.ac.ulg.montefiore.run.jahmm.Observation;

/**
 *
 */
public class ForwardBackwardBayesianNhmmV5ScaledCalculator extends ForwardBackwardBayesianNhmmV5Calculator {

	/*
	 * Warning, the semantic of the alpha and beta elements are changed;
	 * in this class, they have their value scaled.
	 */
	// Scaling factors
	private double[] ctFactors;
	private double[] ctFactorsTmp;
	private double lnProbability;
	
	
	/**
	 * Computes the probability of occurence of an observation sequence
	 * given a Hidden Markov Model.  The algorithms implemented use scaling
	 * to avoid underflows.
	 *
	 * @param hmm A Hidden Markov Model;
	 * @param oseq An observations sequence.
	 * @param flags How the computation should be done. See the
	 *              {@link ForwardBackwardCalculator.Computation}.
	 *              The alpha array is always computed.
	 */
	public <O extends Observation> 
	ForwardBackwardBayesianNhmmV5ScaledCalculator(Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair,
			BayesianNhmmV5<O> hmm, EnumSet<Computation> flags)
	{
		List<? extends O> oseq = oseqPair.getSecond();
		
		if (oseq.isEmpty())
			throw new IllegalArgumentException();
		
		ctFactors = new double[oseq.size()];
		Arrays.fill(ctFactors, 0.);
		ctFactorsTmp = new double[oseq.size()];
		Arrays.fill(ctFactorsTmp, 0.);
		
		computeAlpha(hmm, oseqPair);
		
		if (flags.contains(Computation.BETA))
			computeBeta(hmm, oseqPair);
		
		computeProbability(oseq, hmm, flags);
	}
	
	
	/**
	 * Computes the probability of occurence of an observation sequence
	 * given a Hidden Markov Model.  This computation computes the scaled
	 * <code>alpha</code> array as a side effect.
	 * @see #ForwardBackwardScaledCalculator(List, Hmm, EnumSet)
	 */
	public <O extends Observation>
	ForwardBackwardBayesianNhmmV5ScaledCalculator(Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair,BayesianNhmmV5<O> hmm)
	{
		this(oseqPair, hmm, EnumSet.of(Computation.ALPHA));
	}
	
	
	/* Computes the content of the scaled alpha array */
	protected <O extends Observation> void
	computeAlpha(BayesianNhmmV5<? super O> hmm, Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair)
	{
		List<? extends O> oseq = oseqPair.getSecond();
		HashMap<Integer, Pair<Integer, Double>> cpgDistState = oseqPair.getFirst();
		
		alpha = new double[oseq.size()][hmm.nbStates()];
		alphaTmp = new double[oseq.size()][hmm.nbStates()];
		
		for (int i = 0; i < hmm.nbStates(); i++)
			computeAlphaInit(hmm, oseq.get(0), cpgDistState.get(0), i, oseq.size());
		scale(ctFactors, alpha, 0, hmm);
		//scale(ctFactorsTmp, alphaTmp, 0);
		
		Iterator<? extends O> seqIterator = oseq.iterator();
		if (seqIterator.hasNext())
			seqIterator.next();
		
		for (int t = 1; t < oseq.size(); t++) {
			O observation = seqIterator.next();
			
			for (int i = 0; i < hmm.nbStates(); i++)
				computeAlphaStep(hmm, observation, cpgDistState.get(t), t, i, oseq.size());
			scale(ctFactors, alpha, t, hmm);
			//scale(ctFactorsTmp, alphaTmp, t);
		}
	}
	
	
	/* Computes the content of the scaled beta array.  The scaling factors are
	 those computed for alpha. */
	protected <O extends Observation> void 
	computeBeta(BayesianNhmmV5<? super O> hmm, Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair)
	{	
		List<? extends O> oseq = oseqPair.getSecond();
		HashMap<Integer, Pair<Integer, Double>> cpgDistState = oseqPair.getFirst();
		
		beta = new double[oseq.size()][hmm.nbStates()];
		for (int i = 0; i < hmm.nbStates(); i++)
			beta[oseq.size()-1][i] = 1. / ctFactors[oseq.size()-1];
		
		for (int t = oseq.size() - 2; t >= 0; t--)
			for (int i = 0; i < hmm.nbStates(); i++) {
				computeBetaStep(hmm, oseq.get(t+1), cpgDistState.get(t+1), t, i, oseq.size());
				beta[t][i] /= ctFactors[t];
				if(Double.isNaN(beta[t][i]) || Double.isInfinite(beta[t][i])){
					//beta[t][0]=0.5;
					//beta[t][1]=0.5;
					System.err.println("beta\t" + ctFactors[t] + "\t" + i + "\t" + t);
				}
			}
	/*
		for (int i = 0; i < hmm.nbStates(); i++){
			if(ctFactors[oseq.size()-1]==0){
				beta[oseq.size()-1][i] = 0;
			}else{
				beta[oseq.size()-1][i] = 1. / ctFactors[oseq.size()-1];
			}
		}
			
		
		for (int t = oseq.size() - 2; t >= 0; t--)
			for (int i = 0; i < hmm.nbStates(); i++) {
				computeBetaStep(hmm, oseq.get(t+1), cpgDistState.get(t+1), t, i);
				if(ctFactors[t]==0){
					beta[t][i]=0;
				}else{
					beta[t][i] /= ctFactors[t];
				}
				
			}
			*/
	}
	
	
	/* Normalize alpha[t] and put the normalization factor in ctFactors[t] */
	private <O extends Observation> void scale(double[] ctFactors, double[][] array, int t, BayesianNhmmV5<? super O> hmm)
	{
		double[] table = array[t];
		double sum = 0.;
		
		for (int i = 0; i < table.length; i++)
			sum += table[i];
		
		
		
		if(Double.isNaN(table[0]/sum) || Double.isNaN(table[1]/sum) || Double.isInfinite(table[0]/sum) || Double.isInfinite(table[1]/sum)){
			//System.err.println("scale\t" + table[0] + "\t" + table[1] + "\t" + sum + "\t" + table[0]/sum + "\t" + table[1]/sum);
			//System.err.println(Math.pow(10, Math.log10(table[0])-Math.log10(sum)) + "\t" + Math.pow(10, Math.log10(table[1])-Math.log10(sum)));
			//System.err.println(new BigDecimal(table[0]).divide(new BigDecimal(sum)) + "\t" + new BigDecimal(table[1]).divide(new BigDecimal(sum)));
			
		}
		
		ctFactors[t] = sum;
		for (int i = 0; i < table.length; i++){
			//if(Double.compare(sum, 0.0) == 0 || Double.isInfinite(sum) || Double.isNaN(sum) || Double.isInfinite(table[i]/sum)) {
			if(Double.compare(sum, 0.0) == 0){
		//	if(Double.compare(sum, 0.0) == 0 || (t >=10 && t % 10 == 0)){
				table[i] = 0.5;
				ctFactors[t] = 1;
			}else{
				table[i] /= sum;
			}
		}
		
		//rescale by the number of cpg in the fragment, to avoid too much weight from previous cpg
		/*
		if(t>10){
			sum = 0.;
			for (int i = 0; i < table.length; i++){
				table[i] = Math.exp(Math.log(table[i])/(t-10) + Math.log(0.5));
				sum += table[i];
			}
			ctFactors[t] = sum;
			for (int i = 0; i < table.length; i++){
				table[i] /= sum;
			}
		}
		*/
			
		if(Double.isNaN(table[0]) || Double.isNaN(table[1]) || Double.isInfinite(table[0]) || Double.isInfinite(table[1])){
			System.err.println("scale\t" + table[0] + "\t" + table[1] + "\t" + sum + "\t" + table[0]/sum + "\t" + table[1]/sum);
			System.err.println(Math.pow(10, Math.log10(table[0])-Math.log10(sum)) + "\t" + Math.pow(10, Math.log10(table[1])-Math.log10(sum)));
			
				System.err.println("alpha\t" + t);
				System.err.println( hmm.getOpdf(0));
				System.err.println(hmm.getOpdf(1));
				System.err.println(alpha[0][0] + "\t" + alpha[0][1]);
				System.err.println(alpha[t][0] + "\t" + alpha[t][1]);
				System.exit(1);
			
			System.err.println(new BigDecimal(table[0]).divide(new BigDecimal(sum)) + "\t" + new BigDecimal(table[1]).divide(new BigDecimal(sum)));
			System.exit(1);
		}
		
		
	}
	
	
	private <O extends Observation> void
	computeProbability(List<O> oseq, BayesianNhmmV5<? super O> hmm, 
			EnumSet<Computation> flags)
	{	
		lnProbability = 0.;
		
		for (int t = 0; t < oseq.size(); t++){
			lnProbability += Math.log(ctFactors[t]);
			//if(Double.compare(ctFactorsTmp[t], 0.0) == 0){
			//	continue;
			//}
			//lnProbability += Math.log(ctFactorsTmp[t]);
		}
		
		probability = Math.exp(lnProbability);
	}
	
	
	/**
	 * Return the neperian logarithm of the probability of the sequence that
	 * generated this object.
	 *
	 * @return The probability of the sequence of interest's neperian logarithm.
	 */
	public double lnProbability()
	{
		return lnProbability;
	}
	
	

}
