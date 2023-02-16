/**
 * ViterbiNhmmCalculator.java
 * Apr 27, 2016
 * 2:20:06 PM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.methhmm.hmm;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.util.Pair;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ViterbiCalculator;

/**
 *
 */
public class ViterbiBayesianNhmmV5Calculator {

	/*
	 * The psy and delta values, as described in Rabiner and Juand classical
	 * papers.
	 */
	private double[][] delta; 
	private double[][] deltaWithMethyPrior; 
	private int[][] psy;
	private int[] stateSequence;
	private double lnProbability;
	private double lnProbabilityWithMethyPrior;
	private double pCriteria = 0;
	private int methylatedState = 1;

	
	/**
	 * Computes the most likely state sequence matching an observation
	 * sequence given an HMM.
	 *
	 * @param hmm A Hidden Markov Model;
	 * @param oseq An observations sequence.
	 */
	public <O extends Observation> 
	ViterbiBayesianNhmmV5Calculator(Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseqPair, BayesianNhmmV5<O> hmm, int methylatedState, double pCriteria)
	{
		List<? extends O> oseq = oseqPair.getSecond();
		HashMap<Integer, Pair<Integer, Double>> cpgDistState = oseqPair.getFirst();
		if (oseq.isEmpty())
			throw new IllegalArgumentException("Invalid empty sequence");
		this.pCriteria = pCriteria;
		delta = new double[oseq.size()][hmm.nbStates()];
		deltaWithMethyPrior = new double[oseq.size()][hmm.nbStates()];
		
		psy = new int[oseq.size()][hmm.nbStates()];
		stateSequence = new int[oseq.size()];
		this.methylatedState = methylatedState;
		for (int i = 0; i < hmm.nbStates(); i++) {
			delta[0][i] = -Math.log(hmm.getPri(cpgDistState.get(0).getFirst(), i)) -
					Math.log(hmm.getOpdf(i).probability(oseq.get(0)));
			//delta[0][i] = -Math.log(hmm.getBayesianPri(cpgDistState.get(0), i)) -
			//				Math.log(hmm.getOpdf(i).probability(oseq.get(0)));
					//Math.log(hmm.getOpdfBayesianProb(i, cpgDistState.get(0), oseq.get(0)));
					//Math.log(hmm.getOpdf(i).probability(oseq.get(0)))  - hmm.getBayesianFactor()*Math.log((i % 2 == 0 ? 1-cpgDistState.get(0).getSecond() : cpgDistState.get(0).getSecond()));
			psy[0][i] = 0;
			deltaWithMethyPrior[0][i] = delta[0][i];
		}
		
		Iterator<? extends O> oseqIterator = oseq.iterator();
		O prevO = null;
		if (oseqIterator.hasNext()){
			prevO = oseqIterator.next();
		}
			
		
		int t = 1;
		while (oseqIterator.hasNext()) {
			O observation = oseqIterator.next();
			
			for (int i = 0; i < hmm.nbStates(); i++)
				computeStep(hmm, observation, cpgDistState.get(t), t, i,cpgDistState.get(t-1), prevO, oseq.size());
			prevO = observation;
			t++;
		}
		
		lnProbability = Double.MAX_VALUE;
		lnProbabilityWithMethyPrior = Double.MAX_VALUE;
		/*
		for (int i = 0; i < hmm.nbStates(); i++) {
			double thisProbability = delta[oseq.size()-1][i];
			
			if (lnProbability > thisProbability) {
				lnProbability = thisProbability;
				stateSequence[oseq.size() - 1] = i;
			}
		}
		lnProbability = -lnProbability;
		*/
		
		
		
		double unmethyDelta = delta[oseq.size()-1][1-methylatedState];
		double methyDelta = delta[oseq.size()-1][methylatedState];
		//if(methylatedState == 0){
		//	unmethyDelta = delta[oseq.size()-1][1];
		//	methyDelta = delta[oseq.size()-1][0];
		//}
		
		unmethyDelta = Math.exp(0-unmethyDelta);
		methyDelta = Math.exp(0-methyDelta);
		/*
		double factor = hmm.getBayesianFactor();
		double priorUpperBound = 0.5 + factor/2;
		double priorLowerBound = 0.5 - factor/2;
		double methyPrior = priorLowerBound + cpgDistState.get(t-1).getSecond() * factor;
		double unmethyPrior = priorUpperBound - cpgDistState.get(t-1).getSecond() * factor;
		double methyPosterior = methyDelta*methyPrior/(methyDelta*methyPrior+unmethyDelta*unmethyPrior);
		double unmethyPosterior = unmethyDelta*unmethyPrior/(methyDelta*methyPrior+unmethyDelta*unmethyPrior);
		*/
		
		if((methyDelta/(unmethyDelta+methyDelta)) > (unmethyDelta/(unmethyDelta+methyDelta)+pCriteria)){
			lnProbability = -delta[oseq.size()-1][methylatedState];
			stateSequence[oseq.size() - 1] = methylatedState;
			lnProbabilityWithMethyPrior = -deltaWithMethyPrior[oseq.size()-1][methylatedState];
		}else{
			lnProbability = -delta[oseq.size()-1][1-methylatedState];
			stateSequence[oseq.size() - 1] = 1-methylatedState;
			lnProbabilityWithMethyPrior = -deltaWithMethyPrior[oseq.size()-1][1-methylatedState];
		}
		
		
		for (int t2 = oseq.size() - 2; t2 >= 0; t2--)
			stateSequence[t2] = psy[t2+1][stateSequence[t2+1]];
	}
	
	
	/*
	 * Computes delta and psy[t][j] (t > 0) 
	 */
	private <O extends Observation> void
	computeStep(BayesianNhmmV5<O> hmm, O o, Pair<Integer, Double> r, int t, int j, Pair<Integer, Double> s, O prevO, int numCpg) 
	{
		double minDelta = Double.MAX_VALUE;
		int min_psy = 0;

		/*
		for (int i = 0; i < hmm.nbStates(); i++) {
			
					
			//if(i==1){
			//	thisDelta -= Math.log(pCriteria);
			//}
			if (minDelta > thisDelta) {
				minDelta = thisDelta;
				min_psy = i;
			}
		}
		*/
		
		//double unmethyDelta = delta[t-1][0] + Math.log(hmm.getOpdfProb(0, prevO)) - Math.log(hmm.getOpdfBayesianProb(0, s, prevO)) - Math.log(hmm.getArij(r.getFirst(), 0, j));
		//double methyDelta = delta[t-1][1]  + Math.log(hmm.getOpdfProb(1, prevO)) - Math.log(hmm.getOpdfBayesianProb(1, s, prevO)) - Math.log(hmm.getArij(r.getFirst(), 1, j));
		/*
		double unmethyDelta = delta[t-1][1-methylatedState] + Math.log(hmm.getOpdfProb(1-methylatedState, prevO)) - Math.log(hmm.getOpdfBayesianProb(1-methylatedState, s, prevO)) - Math.log(hmm.getArij(r.getFirst(), 1-methylatedState, j));
		double methyDelta = delta[t-1][methylatedState]  + Math.log(hmm.getOpdfProb(methylatedState, prevO)) - Math.log(hmm.getOpdfBayesianProb(methylatedState, s, prevO)) - Math.log(hmm.getArij(r.getFirst(), methylatedState, j));
		if(noMethyPrior){
			unmethyDelta = delta[t-1][1-methylatedState]  - Math.log(hmm.getArij(r.getFirst(), 1-methylatedState, j));
			methyDelta = delta[t-1][methylatedState] - Math.log(hmm.getArij(r.getFirst(), methylatedState, j));
			
		}
		
		unmethyDelta = Math.exp(0-unmethyDelta);
		methyDelta = Math.exp(0-methyDelta);
				
		if((methyDelta/(unmethyDelta+methyDelta)) > (unmethyDelta/(unmethyDelta+methyDelta)+pCriteria)){
			minDelta = delta[t-1][methylatedState] - Math.log(hmm.getArij(r.getFirst(), methylatedState, j));
			min_psy = methylatedState;
			
		}else{
			minDelta = delta[t-1][1-methylatedState] - Math.log(hmm.getArij(r.getFirst(), 1-methylatedState, j));
			min_psy = 1-methylatedState;
		}
		delta[t][j] = minDelta - Math.log(hmm.getOpdfProb(j,o));
		*/

			/*
			double unmethyLikelihood = delta[t-1][1-methylatedState]  - Math.log(hmm.getArij(r.getFirst(), 1-methylatedState, j));
			double methyLikelihood = delta[t-1][methylatedState]  - Math.log(hmm.getArij(r.getFirst(), methylatedState, j));


			//double unmethyPosterior = Math.exp(0-unmethyLikelihood);
			//double methyPosterior = Math.exp(0-methyLikelihood);
			unmethyLikelihood = Math.exp(0-unmethyLikelihood);
			methyLikelihood = Math.exp(0-methyLikelihood);
			double factor = hmm.getBayesianFactor();
			double priorUpperBound = 0.5 + factor/2;
			double priorLowerBound = 0.5 - factor/2;
			double methyPrior = priorLowerBound + s.getSecond() * factor;
			double unmethyPrior = priorUpperBound - s.getSecond() * factor;
			double methyPosterior = methyLikelihood*methyPrior/(methyLikelihood*methyPrior+unmethyLikelihood*unmethyPrior);
			double unmethyPosterior = unmethyLikelihood*unmethyPrior/(methyLikelihood*methyPrior+unmethyLikelihood*unmethyPrior);

			//double methyPrior = s.getSecond();
			//double methyPosterior = methyLikelihood*methyPrior/(methyLikelihood*methyPrior+unmethyLikelihood*(1-methyPrior));
			//double unmethyPosterior = unmethyLikelihood*(1-methyPrior)/(methyLikelihood*methyPrior+unmethyLikelihood*(1-methyPrior));
			if((methyPosterior/(methyPosterior+unmethyPosterior)) > (unmethyPosterior/(methyPosterior+unmethyPosterior)+pCriteria)){
				minDelta = delta[t-1][methylatedState] - Math.log(hmm.getArij(r.getFirst(), methylatedState, j));
				min_psy = methylatedState;

			}else{
				minDelta = delta[t-1][1-methylatedState] - Math.log(hmm.getArij(r.getFirst(), 1-methylatedState, j));
				min_psy = 1-methylatedState;
			}
			delta[t][j] = minDelta - Math.log(hmm.getOpdfProb(j,o));
			*/


			double unmethyDelta = delta[t-1][1-methylatedState] + Math.log(hmm.getOpdfProb(1-methylatedState, prevO)) - Math.log(hmm.getOpdfBayesianProb(1-methylatedState, s, prevO, numCpg)) - Math.log(hmm.getArij(r.getFirst(), 1-methylatedState, j));
			double methyDelta = delta[t-1][methylatedState]  + Math.log(hmm.getOpdfProb(methylatedState, prevO)) - Math.log(hmm.getOpdfBayesianProb(methylatedState, s, prevO, numCpg)) - Math.log(hmm.getArij(r.getFirst(), methylatedState, j));
			
			
			unmethyDelta = Math.exp(0-unmethyDelta);
			methyDelta = Math.exp(0-methyDelta);
					
			if((methyDelta/(unmethyDelta+methyDelta)) > (unmethyDelta/(unmethyDelta+methyDelta)+pCriteria)){
				minDelta = delta[t-1][methylatedState] - Math.log(hmm.getArij(r.getFirst(), methylatedState, j));
				min_psy = methylatedState;
				
			}else{
				minDelta = delta[t-1][1-methylatedState] - Math.log(hmm.getArij(r.getFirst(), 1-methylatedState, j));
				min_psy = 1-methylatedState;
			}
			delta[t][j] = minDelta - Math.log(hmm.getOpdfProb(j,o));


		psy[t][j] = min_psy;
	}
	
	
	/**
	 * Returns the neperian logarithm of the probability of the given
	 * observation sequence on the most likely state sequence of the given
	 * HMM.
	 *
	 * @return <code>ln(P[O,S|H])</code> where <code>O</code> is the given
	 *         observation sequence, <code>H</code> the given HMM and 
	 *         <code>S</code> the most likely state sequence of this observation
	 *         sequence given this HMM.
	 */
	public double lnProbability(){
		return lnProbability;
	}
	
	
	public double lnProbability(boolean withMethyPrior){
		return lnProbabilityWithMethyPrior;
	}
	
	/**
	 * Returns a (clone of) the array containing the computed most likely
	 * state sequence.
	 *
	 * @return The state sequence; the i-th value of the array is the index
	 *         of the i-th state of the state sequence.
	 */
	public int[] stateSequence(){
		return stateSequence.clone();
	}

}
