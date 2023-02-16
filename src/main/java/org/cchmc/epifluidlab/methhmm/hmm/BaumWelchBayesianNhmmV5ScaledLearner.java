package org.cchmc.epifluidlab.methhmm.hmm;

//import java.math.BigDecimal;
//import java.math.RoundingMode;
import java.util.*;

import org.apache.commons.math3.util.Pair;

import org.cchmc.epifluidlab.methhmm.utils.CcInferenceUtils;
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
//for those of cpg distance without enough instance to estimate, i will use average Aij
public class BaumWelchBayesianNhmmV5ScaledLearner
extends BaumWelchScaledLearner
{	
	/**
	 * Number of iterations performed by the {@link #learn} method.
	 */
	private int nbIterations = 9;
	
	/**
	 * Initializes a Baum-Welch algorithm implementation.
	 */
	public BaumWelchBayesianNhmmV5ScaledLearner()
	{
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
			
		/* gamma and xi arrays are those defined by Rabiner and Juang */
		/* allGamma[n] = gamma array associated to observation sequence n */
		double allGamma[][][] = new double[sequences.size()][][];
	//	double allGammaWithMethy[][][] = new double[sequences.size()][][];
		
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
		
			
			
		//System.err.println("estimate xi and gamma");
		int g = 0;
		for (Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> obsSeqPair : sequences) {	    
			ForwardBackwardBayesianNhmmV5ScaledCalculator fbc = 
				generateForwardBackwardCalculator(obsSeqPair, hmm);
			
			double xi[][][] = estimateXi(obsSeqPair, fbc, hmm);
		//	double xiWithMethy[][][] = estimateXiWithMethy(obsSeqPair, fbc, hmm);
			//System.err.println(hmm.nbCpgDistStates() + "\t" + xi.length + "\t" + xi[0].length + "\t" + xi[0][0].length + "\t" + xi[0][0][0].length);
			
			//double gamma[][] = allGamma[g++] = estimateGamma(xi, fbc, obsSeqPair, hmm);
			double gamma[][] = allGamma[g] = estimateGamma(xi, fbc);
			//allGammaWithMethy[g] = estimateGamma(xiWithMethy, fbc);
			g++;
			//System.err.println(hmm.nbCpgDistStates() + "\t" + allGamma.length + "\t" + allGamma[0].length + "\t" + allGamma[0][0].length + "\t" + allGamma[0][0][0].length);
			
			//double gamma[][][] = allGamma[g];
			//System.err.println(gamma.length + "\t" + gamma[0].length + "\t" + gamma[0][0].length);

			
			List<? extends O> obsSeq = obsSeqPair.getSecond();
			HashMap<Integer, Pair<Integer, Double>> cpgDistState = obsSeqPair.getFirst();
			/*
			for (int i = 0; i < hmm.nbStates(); i++)
				for (int t = 0; t < obsSeq.size() - 1; t++) {
					double unmethy = xi[t][i][0]/(xi[t][i][0]+xi[t][i][1]);
					double methy = xi[t][i][1]/(xi[t][i][0]+xi[t][i][1]);
						unmethy = unmethy * (1-cpgDistState.get(t).getSecond());
						methy = methy * (cpgDistState.get(t).getSecond());
						if(arij.containsKey(cpgDistState.get(t+1).getFirst())){
							Pair<double[][], double[]> tmp = arij.get(cpgDistState.get(t+1).getFirst());
							double[] denTmp = tmp.getSecond();
							double[][] numTmp = tmp.getFirst();
							
							for (int j = 0; j < hmm.nbStates(); j++){
								denTmp[i] += (j % 2 == 0 ? unmethy : methy);
								numTmp[i][j] += (j % 2 == 0 ? unmethy : methy);
							}
								
							arij.put(cpgDistState.get(t+1).getFirst(), new Pair<double[][], double[]>(numTmp, denTmp));
							
						}else{
							double numTmp[][] = new double[hmm.nbStates()][hmm.nbStates()];
							double denTmp[] = new double[hmm.nbStates()];
							Arrays.fill(denTmp, 0.);
							for (int j = 0; j < hmm.nbStates(); j++)
								Arrays.fill(numTmp[j], 0.);
							for (int j = 0; j < hmm.nbStates(); j++){
								denTmp[i] += (j % 2 == 0 ? unmethy : methy);
								numTmp[i][j] += (j % 2 == 0 ? unmethy : methy);
							}
							arij.put(cpgDistState.get(t+1).getFirst(), new Pair<double[][], double[]>(numTmp, denTmp));
						}
						for (int j = 0; j < hmm.nbStates(); j++){
							aijDen[i] += (j % 2 == 0 ? unmethy : methy);
							aijNum[i][j] += (j % 2 == 0 ? unmethy : methy);
						}
							
					
				}
			*/
			
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
		
		
		
		//System.err.println("estimate Arij");
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
		
		
		//System.err.println("estimate Pri");
		/* pi computation */
		for (int r = 1; r < hmm.nbCpgDistState(); r++) {
			for (int i = 0; i < hmm.nbStates(); i++)
				nhmm.setPri(r, i, 0.);
		}
			
		
		
		for (int o = 0; o < sequences.size(); o++){
			Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> obsSeqPair = sequences.get(o);
			Integer r = obsSeqPair.getFirst().get(0).getFirst();
				for (int i = 0; i < hmm.nbStates(); i++){
					//System.err.println(o + "\t" + r + "\t" + i + "\t" + allGamma.length + "\t" + allGamma[0].length + "\t" + allGamma[0][0].length + "\t" + allGamma[o][0][i]);
					
				//	if(Double.isNaN(allGamma[o][0][i])){
					//	nhmm.setPri(r, i,
						//		nhmm.getPri(r, i) );
				//	}else{
						nhmm.setPri(r, i,
								nhmm.getPri(r, i) + allGamma[o][0][i] / sequences.size());
				//	}
						//System.err.println("PRi: " + nhmm.getPri(r,i) + "\t" + r + "\t" + i);
				}

			
				
		}
		
		
		
		
		//System.err.println(nhmm.getPri(0,0) + "\t" + nhmm.getPri(0,1));
		//System.err.println(nhmm.getPri(50,0) + "\t" + nhmm.getPri(50,1));
		//rescale pi
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
		//System.err.println(nhmm.getPri(1,0) + "\t" + nhmm.getPri(1,1));
		//System.err.println(nhmm.getPri(50,0) + "\t" + nhmm.getPri(50,1));
		//System.err.println(sumPi);
			
		//System.err.println("estimate pdfs");
		/* pdfs computation */
		
		for (int i = 0; i < hmm.nbStates(); i++) {
			List<O> observations = CcInferenceUtils.flatPair(sequences);
			double[] weights = new double[observations.size()];
			double sum = 0.;
			
			int j = 0;
			
			int o = 0;
			for (Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> obsSeqPair : sequences) {
				List<? extends O> obsSeq = obsSeqPair.getSecond();
				for (int t = 0; t < obsSeq.size(); t++, j++){

						sum += allGamma[o][t][i];
						weights[j] += allGamma[o][t][i];
						if(Double.isNaN(weights[j]) || Double.isInfinite(weights[j])){
							
							System.err.println(weights[j] + "\t" + sum + "\t" + allGamma[o][t][i] + "\t" + o + "\t" + t + "\t" + i);
							System.exit(1);
						}
				}
					
				o++;
			}
			
			for (j--; j >= 0; j--){
				weights[j] /= sum;
				if(Double.isNaN(weights[j]) || Double.isInfinite(weights[j])){
					
					System.err.println(weights[j] + "\t" + sum + "\t" +  j);
					System.exit(1);
				}
			}
				
			
			Opdf<O> opdf = nhmm.getOpdf(i);
			opdf.fit(observations, weights);
			
		}
		
		/*
		for (int i = 0; i < hmm.nbStates(); i++) {
			List<O> observations = CcInferenceUtils.flatPair(sequences);
			double[] weights = new double[observations.size()];
			double sum = 0.;
			BigDecimal[] weightsBg = new BigDecimal[observations.size()];
			BigDecimal sumBg = BigDecimal.ZERO;
			for (int t = 0; t < observations.size(); t++){
				weightsBg[t] = BigDecimal.ZERO;
			}
			int j = 0;
			
			int o = 0;
			for (Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> obsSeqPair : sequences) {
				List<? extends O> obsSeq = obsSeqPair.getSecond();
				for (int t = 0; t < obsSeq.size(); t++, j++){

						sumBg = sumBg.add(new BigDecimal(allGamma[o][t][i]));
						weightsBg[j] = weightsBg[j].add(new BigDecimal(allGamma[o][t][i]));
						
				}
					
				o++;
			}
			
			for (j--; j >= 0; j--){
				weights[j]=weightsBg[j].divide(sumBg, 20, RoundingMode.HALF_UP).doubleValue();
				
			}
			
			
			
			
			Opdf<O> opdf = nhmm.getOpdf(i);
			opdf.fit(observations, weights);
			
		}
		*/
//System.err.println(hmm);
//System.err.println(nhmm.getPri(1,0) + "\t" + nhmm.getPri(1,1));
//System.err.println(nhmm.getPri(50,0) + "\t" + nhmm.getPri(50,1));
//System.err.println(nhmm.getArij(1,0,0) + "\t" + nhmm.getArij(1,0,1) + "\t" + nhmm.getArij(1,1,0) + "\t" + nhmm.getArij(1,1,1));
//System.err.println(nhmm.getArij(50,0,0) + "\t" + nhmm.getArij(50,0,1) + "\t" + nhmm.getArij(50,1,0) + "\t" + nhmm.getArij(50,1,1));
//System.err.println(nhmm);
		
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
									//hmm.getOpdf(j).probability(observation) *
									//hmm.getOpdfBayesianProb(j, cpgDistState.get(t+1),observation, sequence.size()) *
									hmm.getOpdfProb(j,observation) *
									fbc.betaElement(t + 1, j);
							//if(Double.isNaN(xi[t][i][j]) || Double.compare(xi[t][i][j], 0.0) == 0){
							if(Double.isNaN(xi[t][i][j])){
								System.err.println(t + "\t" + i + "\t" + j + "\t" + cpgDistState.get(t) + "\t" + observation + "\t" +  hmm.getArij(cpgDistState.get(t).getFirst(), i, j));
								System.err.println(fbc.alphaElement(t, i) + "\t" + hmm.getArij(cpgDistState.get(t).getFirst(), i, j) + "\t" + hmm.getOpdf(j).probability(observation) + "\t" + fbc.betaElement(t + 1, j));
								//xi[t][i][j] = 0.;
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
	
	protected <O extends Observation> double[][][]
	estimateXiWithMethy(Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> obsSeqPair, ForwardBackwardBayesianNhmmV5Calculator fbc,
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
									//hmm.getOpdf(j).probability(observation) *
									//hmm.getOpdfBayesianProb(j, cpgDistState.get(t+1),observation, sequence.size()) *
									hmm.getOpdfProb(j,observation) *
									fbc.betaElement(t + 1, j);
							//if(Double.isNaN(xi[t][i][j]) || Double.compare(xi[t][i][j], 0.0) == 0){
							if(Double.isNaN(xi[t][i][j])){
								System.err.println(t + "\t" + i + "\t" + j + "\t" + cpgDistState.get(t) + "\t" + observation + "\t" +  hmm.getArij(cpgDistState.get(t).getFirst(), i, j));
								System.err.println(fbc.alphaElement(t, i) + "\t" + hmm.getArij(cpgDistState.get(t).getFirst(), i, j) + "\t" + hmm.getOpdf(j).probability(observation) + "\t" + fbc.betaElement(t + 1, j));
								//xi[t][i][j] = 0.;
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
	
	protected <O extends Observation> double[][]
	estimateGamma(double[][][] xi, ForwardBackwardBayesianNhmmV5Calculator fbc, Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> obsSeqPair, BayesianNhmmV5<O> hmm)
	{
		HashMap<Integer, Pair<Integer, Double>> cpgDistState = obsSeqPair.getFirst();
		
		double[][] gamma = new double[xi.length + 1][xi[0].length];
		
		for (int t = 0; t < xi.length + 1; t++)
				Arrays.fill(gamma[t], 0.);
		
		for (int t = 0; t < xi.length; t++){
			for (int i = 0; i < xi[0][0].length; i++)
				for (int j = 0; j < xi[0][0].length; j++){
					double unmethy = xi[t][i][0]/(xi[t][i][0]+xi[t][i][1]);
					double methy = xi[t][i][1]/(xi[t][i][0]+xi[t][i][1]);
						unmethy = unmethy * (1-cpgDistState.get(t).getSecond());
						methy = methy * (cpgDistState.get(t).getSecond());
						
					gamma[t][i] +=  (j % 2 == 0 ? unmethy : methy);
					if(Double.isNaN(gamma[t][i])){
						System.err.println(t + "\t" + i + "\t" + j + "\t" + cpgDistState.get(t) + "\t" + xi[t][i][0] + "\t" +  xi[t][i][1]);
						System.err.println(unmethy + "\t" + methy);
						//xi[t][i][j] = 0.;
					}
					//gamma[t][i] += xi[t][i][j];
				}
			
		}
				
						
		
		
			for (int j = 0; j < xi[0][0].length; j++)
				for (int i = 0; i < xi[0][0].length; i++){
					double unmethy = xi[xi.length-1][i][0]/(xi[xi.length-1][i][0]+xi[xi.length-1][i][1]);
					double methy = xi[xi.length-1][i][1]/(xi[xi.length-1][i][0]+xi[xi.length-1][i][1]);
						unmethy = unmethy * (1-cpgDistState.get(xi.length-1).getSecond());
						methy = methy * (cpgDistState.get(xi.length-1).getSecond());
					gamma[xi.length][j] += (j % 2 == 0 ? unmethy : methy);
					if(Double.isNaN(gamma[xi.length][j])){
						System.err.println( i + "\t" + j + "\t" + cpgDistState.get(xi.length-1) + "\t" + xi[xi.length-1][i][0] + "\t" +  xi[xi.length-1][i][1]);
						System.err.println(unmethy + "\t" + methy);
						//xi[t][i][j] = 0.;
					}
					//gamma[xi.length][j] += xi[xi.length - 1][i][j];
				}
					
		
		return gamma;
	}
}

