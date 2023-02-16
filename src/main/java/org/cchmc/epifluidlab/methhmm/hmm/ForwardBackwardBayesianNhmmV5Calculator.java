/**
 * ForwardBackwardNhmmCalculator.java
 * Apr 27, 2016
 * 2:19:18 PM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.methhmm.hmm;

import java.util.EnumSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.util.Pair;
















import be.ac.ulg.montefiore.run.jahmm.Observation;


/**
 *
 */
public class ForwardBackwardBayesianNhmmV5Calculator {

	/**
	 * Flags used to explain how the observation sequence probability
	 * should be computed (either forward, using the alpha array, or backward,
	 * using the beta array).
	 */
	public static enum Computation { ALPHA, BETA };
	
	
	/* alpha[t][i] = P(O(1), O(2),..., O(t+1), i(t+1) = i+1 | hmm), that is the
	 probability of the beginning of the state sequence (up to time t+1)
	 with the (t+1)th state being i+1. */
	protected double[][] alpha = null;
	protected double[][] beta = null;
	protected double probability;
	protected double[][] alphaTmp = null;
	
	
	protected ForwardBackwardBayesianNhmmV5Calculator()
	{
	};
	
	
	/**
	 * Computes the probability of occurence of an observation sequence
	 * given a Hidden Markov Model.
	 *
	 * @param hmm A Hidden Markov Model;
	 * @param oseq An observation sequence.
	 * @param flags How the computation should be done. See the
	 *              {@link Computation Computation} enum.
	 */
	public <O extends Observation>
	ForwardBackwardBayesianNhmmV5Calculator(Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair,
			BayesianNhmmV5<O> hmm, EnumSet<Computation> flags)
	{
		List<? extends O> oseq = oseqPair.getSecond();
		//HashMap<Integer, Integer> cpgDistState = oseqPair.getFirst();
		if (oseq.isEmpty())
			throw new IllegalArgumentException("Invalid empty sequence");
		
		if (flags.contains(Computation.ALPHA))
			computeAlpha(hmm, oseqPair);
		
		if (flags.contains(Computation.BETA))
			computeBeta(hmm, oseqPair);
		
		computeProbability(oseqPair, hmm, flags);
	}
	
	
	/**
	 * Computes the probability of occurence of an observation sequence
	 * given a Hidden Markov Model.  This computation computes the
	 * <code>alpha</code> array as a side effect.
	 * @see #ForwardBackwardCalculator(List, Hmm, EnumSet)
	 */
	public <O extends Observation>
	ForwardBackwardBayesianNhmmV5Calculator(Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair,BayesianNhmmV5<O> hmm)
	{
		this(oseqPair, hmm, EnumSet.of(Computation.ALPHA));
	}
	
	
	/* Computes the content of the alpha array */
	protected <O extends Observation> void
	computeAlpha(BayesianNhmmV5<? super O> hmm, Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair)
	{
		List<? extends O> oseq = oseqPair.getSecond();
		HashMap<Integer, Pair<Integer, Double>> cpgDistState = oseqPair.getFirst();
		
		alpha = new double[oseq.size()][hmm.nbStates()];
		alphaTmp = new double[oseq.size()][hmm.nbStates()];
		for (int i = 0; i < hmm.nbStates(); i++)
			computeAlphaInit(hmm, oseq.get(0), cpgDistState.get(0), i, oseq.size());
		
		if((Double.compare(alpha[0][0], 0.0) == 0 && Double.compare(alpha[0][1], 0.0) == 0)){
			//alpha[0][0]=0.5;
			//alpha[0][1]=0.5;
			System.err.println("alpha\t" + cpgDistState.get(0) + "\t" + cpgDistState.get(0).getSecond());
			System.err.println(Double.compare(alpha[0][0], alpha[0][1]));
			System.err.println(hmm.getPri(cpgDistState.get(0).getFirst(), 0) + "\t" + hmm.getOpdfProb(0,oseq.get(0)) );
			System.err.println(hmm.getPri(cpgDistState.get(0).getFirst(), 1) + "\t" + hmm.getOpdfProb(1,oseq.get(0)) );

			System.exit(1);
		}
		//double unmethy = alpha[0][0]/(alpha[0][0]+alpha[0][1]);
		//double methy = alpha[0][1]/(alpha[0][0]+alpha[0][1]);
		//unmethy = unmethy * (1-cpgDistState.get(0).getSecond());
		//methy = methy * (cpgDistState.get(0).getSecond());
		//alpha[0][0] = unmethy/(unmethy + methy);
		//alpha[0][1] = methy/(unmethy + methy);
		
		Iterator<? extends O> seqIterator = oseq.iterator();
		if (seqIterator.hasNext())
			seqIterator.next();
		
		for (int t = 1; t < oseq.size(); t++) {
			O observation = seqIterator.next();
			
			for (int i = 0; i < hmm.nbStates(); i++)
				computeAlphaStep(hmm, observation, cpgDistState.get(t), t, i, oseq.size());
			
			//double unmethyS = alpha[t][0]/(alpha[t][0]+alpha[t][1]);
			//double methyS = alpha[t][1]/(alpha[t][0]+alpha[t][1]);
			//unmethyS = unmethyS * (1-cpgDistState.get(t).getSecond());
			//methyS = methyS * (cpgDistState.get(t).getSecond());
			//alpha[t][0] = unmethyS/(unmethyS + methyS);
			//alpha[t][1] = methyS/(unmethyS + methyS);
		}
	}
	
	
	/* Computes alpha[0][i] */
	protected <O extends Observation> void
	computeAlphaInit(BayesianNhmmV5<? super O> hmm, O o, Pair<Integer, Double> r, int i, int numCpg)
	{
		//System.err.println(r + "\t" + i);
		//System.err.println(hmm.getPri(r.getFirst(), i));
		//alpha[0][i] = hmm.getPri(r.getFirst(), i) * hmm.getOpdf(i).probability(o);
		//alpha[0][i] = hmm.getPri(r.getFirst(), i) * hmm.getOpdfBayesianProb(i, r,o, numCpg);
		alpha[0][i] = hmm.getPri(r.getFirst(), i) * hmm.getOpdfProb(i,o);
		//alpha[0][i] = hmm.getBayesianPri(r, i) * hmm.getOpdfProb(i,o);
		
		
		//if(Double.compare(alpha[0][i], 0.0) == 0){
		//	alpha[0][i] += 0.001;
		//}
		if(Double.isNaN(alpha[0][i]) || (Double.compare(alpha[0][0], 0.0) == 0 && Double.compare(alpha[0][1], 0.0) == 0)){
			//alpha[0][0]=0.5;
			//alpha[0][1]=0.5;
		//	System.err.println("alpha\t" + r + "\t" + i + "\t" + r.getSecond());
		//	System.err.println(Double.compare(alpha[0][0], alpha[0][1]));
		//	System.err.println(hmm.getPri(r.getFirst(), 0) + "\t" + hmm.getOpdfProb(0,o) );
		//	System.err.println(hmm.getPri(r.getFirst(), 1) + "\t" + hmm.getOpdfProb(1,o) );
		//	System.err.println(hmm.getOpdf(i));
		//	System.exit(1);
		}
	}
	
	
	/* Computes alpha[t][j] (t > 0) */
	protected <O extends Observation> void 
	computeAlphaStep(BayesianNhmmV5<? super O> hmm, O o, Pair<Integer, Double> r, int t, int j, int numCpg)
	{
		double sum = 0.;
		
		for (int i = 0; i < hmm.nbStates(); i++)
			sum += alpha[t-1][i] * hmm.getArij(r.getFirst(), i, j);		

		//alpha[t][j] = sum * hmm.getOpdf(j).probability(o);
		//alpha[t][j] = sum * hmm.getOpdfBayesianProb(j, r,o, numCpg);
		alpha[t][j] = sum * hmm.getOpdfProb(j,o);
		
		if(Double.isNaN(alpha[t][j]) || (Double.compare(alpha[t][0], 0.0) == 0 && Double.compare(alpha[t][1], 0.0) == 0)){
			//alpha[t][0]=0.5;
			//alpha[t][1]=0.5;
			//System.err.println("xxalpha\t" + r + "\t" + j + "\t" + t + "\t" + sum);
			//System.err.println(Double.compare(alpha[t][0], alpha[t][1]));
			//System.err.println(alpha[t-1][0] + "\t" + alpha[t-1][1] + "\t" + hmm.getArij(r.getFirst(), 0, j) + "\t" + hmm.getArij(r.getFirst(), 1, j) + "\t" + hmm.getOpdf(j).probability(o)  + "\t" + o + "\t" + r.getFirst());
			//System.err.println(hmm.getOpdf(j));
			//for (int s = 0; s < t; s++) {
			//	System.err.println(alpha[s][0] + "\t" + alpha[s][1]);
			//}
			//System.exit(1);
		}
	}
	
	
	/* Computes the content of the beta array.  Needs a O(1) access time
	 to the elements of oseq to get a theoretically optimal algorithm. */
	protected <O extends Observation> void 
	computeBeta(BayesianNhmmV5<? super O> hmm, Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair)
	{
		List<? extends O> oseq = oseqPair.getSecond();
		HashMap<Integer, Pair<Integer, Double>> cpgDistState = oseqPair.getFirst();
		
		beta = new double[oseq.size()][hmm.nbStates()];
		
		for (int i = 0; i < hmm.nbStates(); i++)
			beta[oseq.size()-1][i] = 1.;
		//System.err.println(beta[1][0] + "\t" + beta[1][1] + "\t" + (oseq.size()-1));
		for (int t = oseq.size()-2; t >= 0; t--)
			for (int i = 0; i < hmm.nbStates(); i++)
				computeBetaStep(hmm, oseq.get(t+1), cpgDistState.get(t+1), t, i, oseq.size());
	}
	
	
	/* Computes beta[t][i] (t < obs. seq.le length - 1) */
	protected <O extends Observation> void 
	computeBetaStep(BayesianNhmmV5<? super O> hmm, O o, Pair<Integer, Double> r, int t, int i, int numCpg)
	{
		double sum = 0.;
		
		for (int j = 0; j < hmm.nbStates(); j++){
			//if(Double.isNaN(beta[t+1][j]) || Double.isInfinite(beta[t+1][j])){
			//	beta[t+1][0] = 0.5;
			//	beta[t+1][1] = 0.5;
			//}
					
			sum += beta[t+1][j] * hmm.getArij(r.getFirst(), i, j) * 
					//hmm.getOpdfBayesianProb(j, r,o, numCpg);
					hmm.getOpdfProb(j, o);
		}
			
			//hmm.getOpdf(j).probability(o);
		
		//double unmethyS = beta[t+1][0] * hmm.getArij(r.getFirst(), i, 0) * hmm.getOpdf(0).probability(o);
		//double methyS = beta[t+1][1] * hmm.getArij(r.getFirst(), i, 1) * hmm.getOpdf(1).probability(o);
		//double unmethy = unmethyS/(unmethyS + methyS);
		//double methy = methyS/(unmethyS + methyS);
		//unmethy = unmethy * (1-r.getSecond());
		//methy = methy * (r.getSecond());
		//beta[t][i] = unmethy + methy;
		
		beta[t][i] = sum;
		if(Double.isNaN(beta[t][i]) || Double.isInfinite(beta[t][i])){
		//	beta[t][0]=0.5;
		//	beta[t][1]=0.5;
		//	System.err.println("beta\t" + r + "\t" + i + "\t" + t);
		//	System.err.println(beta[t+1][0] + "\t" + beta[t+1][1] + "\t" +  hmm.getArij(r.getFirst(), i, 0) + "\t" + hmm.getArij(r.getFirst(), i, 1) + "\t" + hmm.getOpdf(0).probability(o)  + "\t" + hmm.getOpdf(1).probability(o));
		//	System.err.println(beta[t][0] + "\t" + beta[t][1]);
		
		//	System.err.println(hmm.getOpdf(i));
		//if(Double.isNaN(beta[t][i])){
			//System.exit(1);
		}
	}
	
	
	/**
	 * Returns an element of the <i>alpha</i> array.
	 * 
	 * @param t The temporal argument of the array (positive but strictly
	 *          smaller than the length of the sequence that helped generating
	 *          the array).
	 * @param i A state index of the HMM that helped generating the array.
	 * @throws {@link UnsupportedOperationException 
	 *          UnsupportedOperationException} if alpha array has not been
	 *          computed.
	 * @return The <i>alpha</i> array (t, i) element.
	 */ 
	public double alphaElement(int t, int i)
	{
		if (alpha == null)
			throw new UnsupportedOperationException("Alpha array has not " +
					"been computed");
		//if(Double.compare(alpha[t][i], 0.0) == 0)
		//	System.err.println(alpha[t][i] + "\t" + alpha[t][1-i] + "\t" + t + "\t" + i +  "\talpha");
		return alpha[t][i];
	}
	
	
	/**
	 * Returns an element of the <i>beta</i> array.
	 * 
	 * @param t The temporal argument of the array (positive but smaller than
	 *          the length of the sequence that helped generating the array).
	 * @param i A state index of the HMM that helped generating the array.
	 * @throws {@link UnsupportedOperationException 
	 *          UnsupportedOperationException} if beta array has not been
	 *          computed.
	 * @return The <i>beta</i> beta (t, i) element.
	 */ 
	public double betaElement(int t, int i)
	{
		if (beta == null)
			throw new UnsupportedOperationException("Beta array has not " +
					"been computed");
		
		return beta[t][i];
	}
	
	
	private <O extends Observation> void 
	computeProbability(Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair, BayesianNhmmV5<? super O> hmm, 
			EnumSet<Computation> flags)
	{
		probability = 0.;
		List<? extends O> oseq = oseqPair.getSecond();
		HashMap<Integer, Pair<Integer, Double>> cpgDistState = oseqPair.getFirst();
		
		if (flags.contains(Computation.ALPHA))
			for (int i = 0; i < hmm.nbStates(); i++) 
				probability += alpha[oseq.size()-1][i];
		else
			for (int i = 0; i < hmm.nbStates(); i++){
				//double unmethy = Math.log(hmm.getOpdf(i).probability(oseq.get(0)))+ Math.log(1-cpgDistState.get(i).getSecond() )*hmm.getBayesianFactor();
				//double methy = Math.log(hmm.getOpdf(i).probability(oseq.get(0))) + Math.log(cpgDistState.get(i).getSecond())*hmm.getBayesianFactor();
				//unmethy = Math.exp(unmethy);
				//methy = Math.exp(methy);
				//probability += hmm.getPri(cpgDistState.get(0).getFirst(), i) * (i % 2 == 0 ? unmethy : methy)* beta[0][i];
				probability +=
						//hmm.getBayesianPri(cpgDistState.get(0),i) *
								hmm.getPri(cpgDistState.get(0).getFirst(),i) *
					//hmm.getOpdfBayesianProb(i, cpgDistState.get(0),oseq.get(0), oseq.size()) * beta[0][i];
					hmm.getOpdfProb(i, oseq.get(0)) * beta[0][i];
				//	hmm.getOpdf(i).probability(oseq.get(0)) * beta[0][i];
			}
				
	}
	
	
	/**
	 * Return the probability of the sequence that generated this object.
	 * For long sequences, this probability might be very small, or even
	 * meaningless because of underflows.
	 *
	 * @return The probability of the sequence of interest.
	 */
	public double probability()
	{
		return probability;
	}

}
