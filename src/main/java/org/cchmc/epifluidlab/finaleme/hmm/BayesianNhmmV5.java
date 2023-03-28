/**
 * Nhmm.java
 * Apr 27, 2016
 * 12:43:46 PM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.finaleme.hmm;

import java.io.Serializable;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import org.apache.commons.math3.util.Pair;

import be.ac.ulg.montefiore.run.jahmm.Observation;
import org.cchmc.epifluidlab.finaleme.utils.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfFactory;


/**
 *
 */
public class BayesianNhmmV5<O extends Observation> 
implements Serializable, Cloneable {

	
	private TreeMap<Integer, Double[]> pi;
	private TreeMap<Integer, Double[][]>  a;
	private ArrayList<Opdf<O>> opdfs;
	private int nbStates;
	private int nbCpgDistStates;
	private double bayesianFactor;
	private int methyState = 1;
	private double maxCpgInFrag;
	private double minCpgInFrag;
	//private GammaDistribution gamma;
	//private ExponentialDistribution expDist;
	//private WeibullDistribution dist;
	/**
	 * 
	 */
	private static final long serialVersionUID = -1203330572216045699L;


	public BayesianNhmmV5(int nbStates, int nbCpgDistStates, double bayesianFactor) {
		
		if (nbStates <= 0)
			throw new IllegalArgumentException("Number of states must be " +
			"positive");
		this.nbStates = nbStates;
		this.nbCpgDistStates = nbCpgDistStates;
		pi = new TreeMap<Integer, Double[]>();
		this.bayesianFactor = bayesianFactor;
		a = new TreeMap<Integer, Double[][]>();
		for(int i = 0; i<= nbCpgDistStates; i++){
			a.put(i, new Double[nbStates][nbStates]);
			pi.put(i, new Double[nbStates]);
		}
		opdfs = new ArrayList<Opdf<O>>(nbStates);
		for (int i = 0; i < nbStates; i++)
			opdfs.add(null);
		//gamma = new GammaDistribution(0.5,1.5);
		//expDist = new ExponentialDistribution(1);
		//dist = new WeibullDistribution(bayesianFactor,1);
	}
	
	public BayesianNhmmV5(int nbStates, int nbCpgDistStates, OpdfFactory<? extends Opdf<O>> opdfFactory, double bayesianFactor)
	{
		if (nbStates <= 0)
			throw new IllegalArgumentException("Number of states must be " +
			"strictly positive");
		this.nbStates = nbStates;
		this.nbCpgDistStates = nbCpgDistStates;
		pi = new TreeMap<Integer, Double[]>();
		a = new TreeMap<Integer, Double[][]>();
		opdfs = new ArrayList<Opdf<O>>(nbStates);
		this.bayesianFactor = bayesianFactor;
			
			for(int n = 0; n<= nbCpgDistStates; n++){
				Double[][] t = new Double[nbStates][nbStates];
				Double[] pit = new Double[nbStates];
				for (int i = 0; i < nbStates; i++) {
					pit[i] = 1. / ((double) nbStates);
					for (int j = 0; j < nbStates; j++)
						t[i][j] = 1. / ((double) nbStates);
				}
				a.put(n, t);
				pi.put(n, pit);
			}
		
		for (int i = 0; i < nbStates; i++) {
			opdfs.add(opdfFactory.factor());
		}
		//gamma = new GammaDistribution(0.5,1.5);
		//expDist = new ExponentialDistribution(1);
		//dist = new WeibullDistribution(bayesianFactor,1);
	}

	
	
	/**
	 * Returns the <i>pi</i> value associated with a given state.
	 *
	 * @param stateNb A state number such that
	 *                <code>0 &le; stateNb &lt; nbStates()</code>
	 * @return The <i>pi</i> value associated to <code>stateNb</code>.
	 */
	//public double getPi(int stateNb)
	//{
	//	return pi[stateNb];
	//}
	
	
	
	/**
	 * Returns the <i>pi</i> value associated with a given state.
	 *
	 * @param stateNb A state number such that
	 *                <code>0 &le; stateNb &lt; nbStates()</code>
	 * @return The <i>pi</i> value associated to <code>stateNb</code>.
	 */
	public double getBayesianPri(Pair<Integer, Double> nbCpgDistStates, int stateNb)//even states represent unmethylated, odd states represent methylated, states number could not be odd asnymore
	{
		double methyPrior = nbCpgDistStates.getSecond();
		int dist = nbCpgDistStates.getFirst();
		//double state1 = stateNb==0 ? pi.get(dist)[stateNb]*(100-methyPrior) : pi.get(dist)[stateNb]*methyPrior; //0 is unmethylated, 1 is methlated
		//double state2 = (1-stateNb)==0 ? pi.get(dist)[1-stateNb]*(100-methyPrior) : pi.get(dist)[1-stateNb]*methyPrior; //0 is unmethylated, 1 is methlated
		/*
		double state1 = pi.get(dist)[stateNb]*(100-methyPrior)/100;
		double state2 = pi.get(dist)[1-stateNb]*methyPrior/100;
		if(stateNb!=0){
			state1 = pi.get(dist)[stateNb]*methyPrior/100;
			state2 = pi.get(dist)[1-stateNb]*(100-methyPrior)/100;
		}
		*/
		/*
		if(pi.get(dist) == null){
			System.err.println(dist + "\t" + stateNb);
			System.err.println(pi.get(dist).length);
		}
		double methyPriorScale = methyPrior/(100*nbStates()/2);
		double unmethyPriorScale = (100-methyPrior)/(100*nbStates()/2);
		double stateSum = 0.;
		for(int i = 0; i < nbStates(); i +=2 ){
			stateSum = stateSum + pi.get(dist)[i] + unmethyPriorScale + pi.get(dist)[i+1] + methyPriorScale;
			
		}
		double state;
		if(stateNb % 2 ==0){
			state = pi.get(dist)[stateNb] + unmethyPriorScale;
		}else{
			state = pi.get(dist)[stateNb] + methyPriorScale;
		}
		//System.err.println(methyPrior + "\t" + methyPriorScale + "\t" + unmethyPriorScale + "\t" + stateSum + "\t" + pi.get(dist)[stateNb] + "\t" + stateNb);
		*/
		double methyPriorScale = pi.get(dist)[1]*methyPrior;
		double unmethyPriorScale = pi.get(dist)[0]*(1-methyPrior);
		//double methyPriorScale = methyPrior+pi.get(dist)[1];
		//double unmethyPriorScale = 1-methyPrior+pi.get(dist)[0];
		//double methyPriorScale = 1-methyPrior+pi.get(dist)[1];
		//double unmethyPriorScale = methyPrior+pi.get(dist)[0];
		//double methyPriorScale = 0+pi.get(dist)[1];
		//double unmethyPriorScale = 1+pi.get(dist)[0];
		//System.err.println(methyPrior + "\t" + methyPriorScale + "\t" + unmethyPriorScale  + "\t" + pi.get(dist)[1]  + "\t" + pi.get(dist)[0] + "\t" + stateNb);
		if(stateNb % 2 ==0){
			return unmethyPriorScale/(methyPriorScale + unmethyPriorScale);
			//return unmethyPriorScale;
		}else{
			return methyPriorScale/(methyPriorScale + unmethyPriorScale);
			//return methyPriorScale;
		}
		
		//return state/stateSum; 
	}
	
	public double getPri(int nbCpgDistStates, int stateNb)
	{
		//System.err.println(nbCpgDistStates + "\t" + stateNb + "\t" + pi.length + "\t" + pi[0].length + "\t" + pi[nbCpgDistStates][stateNb]);
		return pi.get(nbCpgDistStates)[stateNb]; //TODO: should be 
	}
	
	public void setPri(int nbCpgDistStates, int stateNb, double value)
	{
		//System.err.println(nbCpgDistStates + "\t" + stateNb + "\t" + pi.length + "\t" + pi[0].length + "\t" + pi[nbCpgDistStates][stateNb]);
		Double[] t = pi.get(nbCpgDistStates);
		t[stateNb]=value;
		pi.put(nbCpgDistStates, t); 
	}
	
	public double getBayesianFactor()
	{
		return bayesianFactor;
	}
	
	public void setBayesianFactor(double value)
	{
		this.bayesianFactor = value;
	}
	
	public void setMaxCpgNum(double v){
		this.maxCpgInFrag = v;
	}

	public void setMinCpgNum(double v){
		this.minCpgInFrag = v;
	}

	/**
	 * Sets the <i>pi</i> value associated with a given state.
	 *
	 * @param stateNb A state number such that
	 *                <code>0 &le; stateNb &lt; nbStates()</code>.
	 * @param value The <i>pi</i> value to associate to state number
	 *              <code>stateNb</code>
	 */
	//public void setPi(int stateNb, double value)
	//{
	//	System.err.print(nbCpgDistStates + "\t" + stateNb + "\t" + value + "\t" + pi[nbCpgDistStates][stateNb] + "\t");
	//	pi[stateNb] = value;
	//	System.err.print(pi[nbCpgDistStates][stateNb] + "\n");
	//}
	
	/**
	 * Returns the probability associated with the transition going from
	 * state <i>i</i> to state <i>j</i> (<i>a<sub>i,j</sub></i>).
	 *
	 *@param r The input variable that could affect transition probability
	 *
	 * @param i The first state number such that
	 *        <code>0 &le; i &lt; nbStates()</code>.
	 * @param j The second state number such that
	 *        <code>0 &le; j &lt; nbStates()</code>.
	 * @return The probability associated to the transition going from
	 *         <code>i</code> to state <code>j</code>.
	 */
	public double getBayesianArij(Pair<Integer, Double> nbCpgDistStates, int i, int j)
	{
		double methyPrior = nbCpgDistStates.getSecond();
		int dist = nbCpgDistStates.getFirst();
		
		/*
		double methyPriorScale = methyPrior/(100*nbStates()/2);
		double unmethyPriorScale = (100-methyPrior)/(100*nbStates()/2);
		double stateSum = 0.;
		for(int z = 0; z < nbStates(); z +=2 ){
			stateSum = stateSum + a.get(dist)[i][z] + unmethyPriorScale + a.get(dist)[i][z+1] + methyPriorScale;
			
		}
		double state;
		if(j % 2 ==0){
			state = a.get(dist)[i][j] + unmethyPriorScale;
		}else{
			state = a.get(dist)[i][j] + methyPriorScale;
		}
		
		return state/stateSum;
		*/
		double methyPriorScale = 1-methyPrior+a.get(dist)[i][1];
		double unmethyPriorScale = methyPrior+a.get(dist)[i][0];
		//double methyPriorScale = 1.5+a.get(dist)[i][1];
		//double unmethyPriorScale = 0+a.get(dist)[i][0];
		//System.err.println(methyPrior + "\t" + methyPriorScale + "\t" + unmethyPriorScale  + "\t" + pi.get(dist)[1]  + "\t" + pi.get(dist)[0] + "\t" + stateNb);
		if(j % 2 ==0){
			return unmethyPriorScale/(methyPriorScale + unmethyPriorScale);
		}else{
			return methyPriorScale/(methyPriorScale + unmethyPriorScale);
		}
		
	}
	
	public double getArij(int nbCpgDistStates, int i, int j)
	{
		return a.get(nbCpgDistStates)[i][j];
	}
	
	public void setArij(int nbCpgDistStates, int i, int j, double value)
	{
		Double t[][] = a.get(nbCpgDistStates);
		t[i][j]=value;
		a.put(nbCpgDistStates, t);
	}
	
	
	/**
	 * Sets the probability associated to the transition going from
	 * state <i>i</i> to state <i>j</i> (<i>A<sub>i,j</sub></i>).
	 *
	 *@param r The input variable that could affect transition probability
	 *
	 * @param i The first state number such that
	 *        <code>0 &le; i &lt; nbStates()</code>.
	 * @param j The second state number such that
	 *        <code>0 &le; j &lt; nbStates()</code>.
	 * @param value The value of <i>A<sub>i,j</sub></i>.
	 */
//	public void setAij(int i, int j, double value)
//	{
//		a.put(key, value)
//		a[i][j] = value;
//	}
	
//	public double getAij(int i, int j)
//	{
//		return a[i][j];
//	}
	
	public int[] mostLikelyStateSequence(Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseq, int methylatedState, double pCriteria)
	{
		return (new ViterbiBayesianNhmmV5Calculator(oseq, this, methylatedState, pCriteria)).stateSequence();
	}
	
	

	public double probability(Pair<HashMap<Integer, Pair<Integer, Double>>, List<O>> oseq)
	{
		return (new ForwardBackwardBayesianNhmmV5Calculator(oseq, this)).probability();
	}
	
	

	public double lnProbability(Pair<HashMap<Integer, Pair<Integer, Double>>, List< O>> oseq)
	{
		return (new ForwardBackwardBayesianNhmmV5ScaledCalculator(oseq, this)).lnProbability();
	}
	
	
	/**
	 * Returns the probability of an observation sequence along a state
	 * sequence given this HMM.
	 *
	 * @param oseq A non-empty observation sequence.
	 * @param sseq An array containing a sequence of state numbers.  The length
	 *             of this array must be equal to the length of
	 *             <code>oseq</code>
	 * @return The probability P[oseq,sseq|H], where H is this HMM.
	 */
	public double probability(Pair<HashMap<Integer, Pair<Integer, Double>>, List<? extends O>> oseqPair, int[] sseq)
	{
		
		List<? extends O> oseq = oseqPair.getSecond();
		HashMap<Integer, Pair<Integer, Double>> cpgDistState = oseqPair.getFirst();
		
		if (oseq.size() != sseq.length || oseq.isEmpty())
			throw new IllegalArgumentException();
		
		double probability = getPri(cpgDistState.get(0).getFirst(),sseq[0]);
		
		Iterator<? extends O> oseqIterator = oseq.iterator();
		
		for (int i = 1; i < sseq.length-1; i++)
			probability *= 
				getOpdf(sseq[i]).probability(oseqIterator.next()) *
				getArij(cpgDistState.get(i).getFirst(), sseq[i], sseq[i+1]); //TODO: need double check here
		
		return probability * getOpdf(sseq[sseq.length-1]).
		probability(oseq.get(sseq.length-1));
	}
	
	public double getProbability(Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>> oseqPair, int[] sseq)
	{
		
		List<ObservationVector> oseq = oseqPair.getSecond();
		HashMap<Integer, Pair<Integer, Double>> cpgDistState = oseqPair.getFirst();
		
		if (oseq.size() != sseq.length || oseq.isEmpty())
			throw new IllegalArgumentException();
		
		double probability = getPri(cpgDistState.get(0).getFirst(),sseq[0]);
		
		Iterator<ObservationVector> oseqIterator = oseq.iterator();
		
		for (int i = 1; i < sseq.length-1; i++)
			probability *= 
				getOpdf(sseq[i]).probability((O)oseqIterator.next()) *
				getArij(cpgDistState.get(i).getFirst(), sseq[i], sseq[i+1]); //TODO: need double check here
		
		return probability * getOpdf(sseq[sseq.length-1]).
		probability((O)oseq.get(sseq.length-1));
	}

	/**
	 * Returns the number of states of this HMM.
	 *
	 * @return The number of states of this HMM.
	 */
	public int nbStates()
	{
		return this.nbStates;
	}
	
	public int nbCpgDistState()
	{
		return this.nbCpgDistStates;
	}
	
	/**
	 * Returns the opdf associated with a given state.
	 *
	 * @param stateNb A state number such that
	 *                <code>0 &le; stateNb &lt; nbStates()</code>.
	 * @return The opdf associated to state <code>stateNb</code>.
	 */
	public Opdf<O> getOpdf(int stateNb)
	{
		//System.err.println(opdfs.size() + "\t" + stateNb);
		return opdfs.get(stateNb);
	}
	
	/*
	public double getOpdfBayesianProb(int stateNb, Pair<Integer, Double> nbCpgDistStates, O o)
	{
		double methyPrior = nbCpgDistStates.getSecond()*getBayesianFactor();
		double unmethyLikelihood = opdfs.get(0).probability(o);
		double methyLikelihood = opdfs.get(1).probability(o);
		//System.err.println(stateNb + "\t" + methyPrior + "\t" + unmethyLikelihood + "\t" + methyLikelihood + "\t" + unmethyLikelihood*(1-methyPrior)/(unmethyLikelihood*(1-methyPrior) + methyLikelihood*methyPrior) + "\t" 
		//		+ methyLikelihood*methyPrior/(unmethyLikelihood*(1-methyPrior) + methyLikelihood*methyPrior));
		return stateNb % 2 ==0 ? unmethyLikelihood*(1-methyPrior)/(unmethyLikelihood*(1-methyPrior) + methyLikelihood*methyPrior) : methyLikelihood*methyPrior/(unmethyLikelihood*(1-methyPrior) + methyLikelihood*methyPrior);
		//return opdfs.get(stateNb);
	}
	*/


	public double getOpdfBayesianProb(int stateNb, Pair<Integer, Double> nbCpgDistStates, O o, double cpgNum)
	{
		
		//double cpgScaled = (double)(cpgNum-minCpgInFrag)/(double)(maxCpgInFrag-minCpgInFrag);
		//double cpgScaled = (double)(cpgNum-1)/(double)(50-1);
		//cpgScaled = cpgScaled > 1.0 ? 1.0 : cpgScaled;
		//double factor = gamma.density(cpgScaled);
		//double factor = expDist.density(cpgScaled);
		//double factor = dist.density(cpgScaled);
		//double factor = 1.0;
		//factor = factor > getBayesianFactor() ? getBayesianFactor() : factor;
		//factor = factor > 1.0 ? 1.0 : factor;
		//double factor = cpgNum >maxCpgInFrag ? getBayesianFactor() : 1;
		double factor = getBayesianFactor();
		double priorUpperBound = 0.5 + factor/2;
		double priorLowerBound = 0.5 - factor/2;
		double methyPrior = priorLowerBound + nbCpgDistStates.getSecond() * factor;
		double unmethyPrior = priorUpperBound - nbCpgDistStates.getSecond() * factor;
		double methyLikelihood = opdfs.get(methyState).probability(o);
		double unmethyLikelihood = opdfs.get(1-methyState).probability(o);
		
		double unmethyLikelihoodScale = unmethyLikelihood/(unmethyLikelihood + methyLikelihood) ;
		double methyLikelihoodScale = methyLikelihood/(unmethyLikelihood + methyLikelihood);
				
		//System.err.println(stateNb  + "\t" + nbCpgDistStates.getSecond()+  "\t" + maxCpgInFrag + "\t" + minCpgInFrag + "\t" +  cpgScaled + "\t" +factor + "\t" + methyPrior  + "\t" + unmethyPrior + "\t" + unmethyLikelihood + "\t" + methyLikelihood  + "\t"  + methyLikelihoodScale + "\t" + methyLikelihoodScale*methyPrior/(methyLikelihoodScale*methyPrior+unmethyLikelihoodScale*unmethyPrior)
		//		+ "\t" + unmethyLikelihoodScale*unmethyPrior/(methyLikelihoodScale*methyPrior+unmethyLikelihoodScale*unmethyPrior));
		
		return stateNb == methyState ? methyLikelihoodScale*methyPrior/(methyLikelihoodScale*methyPrior+unmethyLikelihoodScale*unmethyPrior) : 
			unmethyLikelihoodScale*unmethyPrior/(methyLikelihoodScale*methyPrior+unmethyLikelihoodScale*unmethyPrior);
	}
	
	/*
	public double getOpdfBayesianProb(int stateNb, Pair<Integer, Double> nbCpgDistStates, O o, double factor)
	{
		double priorUpperBound = 0.5 + factor/2;
		double priorLowerBound = 0.5 - factor/2;
		double methyPrior = priorLowerBound + nbCpgDistStates.getSecond() * factor;
		double unmethyPrior = priorUpperBound - nbCpgDistStates.getSecond() * factor;
		double methyLikelihood = opdfs.get(methyState).probability(o);
		double unmethyLikelihood = opdfs.get(1-methyState).probability(o);
		
		double unmethyLikelihoodScale = unmethyLikelihood/(unmethyLikelihood + methyLikelihood) ;
		double methyLikelihoodScale = methyLikelihood/(unmethyLikelihood + methyLikelihood);
				
		//System.err.println(stateNb  + "\t" + nbCpgDistStates.getSecond()+ "\t" + methyPrior  + "\t" + unmethyPrior + "\t" + unmethyLikelihood + "\t" + methyLikelihood  + "\t"  + methyLikelihoodScale + "\t" + methyLikelihoodScale*methyPrior/(methyLikelihoodScale*methyPrior+unmethyLikelihoodScale*unmethyPrior)
		//		+ "\t" + unmethyLikelihoodScale*unmethyPrior/(methyLikelihoodScale*methyPrior+unmethyLikelihoodScale*unmethyPrior));
		
		return stateNb == methyState ? methyLikelihoodScale*methyPrior/(methyLikelihoodScale*methyPrior+unmethyLikelihoodScale*unmethyPrior) : 
			unmethyLikelihoodScale*unmethyPrior/(methyLikelihoodScale*methyPrior+unmethyLikelihoodScale*unmethyPrior);
		
		/*
		double methyPrior = Math.log(nbCpgDistStates.getSecond()) * getBayesianFactor();
		double unmethyPrior = Math.log(1-nbCpgDistStates.getSecond()) * getBayesianFactor();
		double methyLikelihood = opdfs.get(methyState).probability(o);
		double unmethyLikelihood = opdfs.get(1-methyState).probability(o);
		
		double unmethyLikelihoodScale = unmethyLikelihood/(unmethyLikelihood + methyLikelihood) ;
		double methyLikelihoodScale = methyLikelihood/(unmethyLikelihood + methyLikelihood);
		
		double methyPosterior = Math.exp(Math.log(methyLikelihoodScale) + methyPrior);
		double unmethyPosterior = Math.exp(Math.log(unmethyLikelihoodScale) + unmethyPrior);
		
		//System.err.println(stateNb  + "\t" + nbCpgDistStates.getSecond()+ "\t" + methyPrior  + "\t" + unmethyLikelihood + "\t" + methyLikelihood  + "\t"  + methyLikelihoodScale + "\t" + methyPosterior);
		
		return stateNb == methyState ? methyPosterior/(methyPosterior+unmethyPosterior) : unmethyPosterior/(unmethyPosterior+methyPosterior);
		*/
		/*
		double methyPrior = nbCpgDistStates.getSecond() * getBayesianFactor();
		double unmethyPrior = (1-nbCpgDistStates.getSecond()) * 1/getBayesianFactor();
		double unmethyLikelihood = opdfs.get(1-methyState).probability(o);
		double methyLikelihood = opdfs.get(methyState).probability(o);
		double unmethyLikelihoodScale = unmethyLikelihood/(unmethyLikelihood + methyLikelihood) ;
		double methyLikelihoodScale = methyLikelihood/(unmethyLikelihood + methyLikelihood);
		//double unmethyLikelihoodScale = unmethyLikelihood/(unmethyLikelihood + methyLikelihood) + 1- methyPrior;
		//double methyLikelihoodScale = methyLikelihood/(unmethyLikelihood + methyLikelihood) + methyPrior;
		//double sum = unmethyLikelihood + methyLikelihood;
		GammaDistribution gamma = new GammaDistribution(1.0,1.0);
		double gammashape = gamma.density(Math.abs(methyLikelihoodScale-unmethyLikelihoodScale));
		if(Double.compare((unmethyLikelihood + methyLikelihood),0.0)==0){
			//	unmethyLikelihoodScale = gammashape*unmethyPrior/(methyPrior+unmethyPrior);
			//	methyLikelihoodScale = gammashape*methyPrior/(methyPrior+unmethyPrior);
				unmethyLikelihoodScale = 0.5;
				methyLikelihoodScale = 0.5;

			}
	//	if(Math.max(methyLikelihoodScale/unmethyLikelihoodScale, unmethyLikelihoodScale/methyLikelihoodScale) >= 2){
	//		return stateNb == methyState ? methyLikelihoodScale/(unmethyLikelihoodScale + methyLikelihoodScale) : unmethyLikelihoodScale/(unmethyLikelihoodScale + methyLikelihoodScale) ;
	//	}else{
			

	//	}
		
		if(methyLikelihoodScale - unmethyLikelihoodScale > 0.3){
			System.err.println(stateNb  + "\t" + gammashape+ "\t" + nbCpgDistStates.getSecond()+ "\t" + methyPrior + "\t" + unmethyPrior + "\t" + unmethyLikelihood + "\t" + methyLikelihood + "\t" + unmethyLikelihoodScale + "\t"  + methyLikelihoodScale + "\t" + unmethyLikelihoodScale*unmethyPrior/(unmethyLikelihoodScale*unmethyPrior + methyLikelihoodScale*methyPrior) + "\t" 
				+ methyLikelihoodScale*methyPrior/(unmethyLikelihoodScale*unmethyPrior + methyLikelihoodScale*methyPrior));
		}
		//System.err.println(stateNb  + "\t" + nbCpgDistStates.getSecond()+ "\t" + methyPrior  + "\t" + unmethyLikelihood + "\t" + methyLikelihood + "\t" + unmethyLikelihoodScale + "\t"  + methyLikelihoodScale + "\t" + unmethyLikelihoodScale/(unmethyLikelihoodScale + methyLikelihoodScale) + "\t" 
		//		+ methyLikelihoodScale/(unmethyLikelihoodScale + methyLikelihoodScale));

		return stateNb == methyState ? gammashape*methyLikelihoodScale*methyPrior/(unmethyLikelihoodScale*unmethyPrior + methyLikelihoodScale*methyPrior) : gammashape*unmethyLikelihoodScale*unmethyPrior/(unmethyLikelihoodScale*unmethyPrior + methyLikelihoodScale*methyPrior) ;
		//return stateNb == methyState ? methyLikelihoodScale/(unmethyLikelihoodScale + methyLikelihoodScale) : unmethyLikelihoodScale/(unmethyLikelihoodScale + methyLikelihoodScale) ;
		
		//return stateNb == methyState ? methyLikelihoodScale*methyPrior/(unmethyLikelihoodScale*unmethyPrior + methyLikelihoodScale*methyPrior) : unmethyLikelihoodScale*unmethyPrior/(unmethyLikelihoodScale*unmethyPrior + methyLikelihoodScale*methyPrior) ;
		//return opdfs.get(stateNb).probability(o);
	}
	*/

	public double getOpdfProb(int stateNb, O o)
	{
		
		return opdfs.get(stateNb).probability(o);
	}

	
	/**
	 * Sets the opdf associated with a given state.
	 *
	 * @param stateNb A state number such that
	 *                <code>0 &le; stateNb &lt; nbStates()</code>.
	 * @param opdf An observation probability function.
	 */
	public void setOpdf(int stateNb, Opdf<O> opdf)
	{
		opdfs.set(stateNb, opdf);
		//System.err.println(opdfs.size() + "\t" + stateNb);
	}
	
	
	public void setMethyState(int stateNb)
	{
		this.methyState = stateNb;
		//System.err.println(opdfs.size() + "\t" + stateNb);
	}
	
	public int getMethyState(){
		int methyState = 0;
		double coverage = Double.NEGATIVE_INFINITY;
		for(int z = 0; z < nbStates(); z++){
			if(((OpdfMultiMixtureGaussian)getOpdf(z)).mean()[2] > coverage){
				coverage = ((OpdfMultiMixtureGaussian)getOpdf(z)).mean()[2];
				methyState = z;
			}
		}
		
		return methyState;
	}
	
	public int getMethyState(boolean lowcoverage){
		int methyState = 0;
		double coverage = Double.NEGATIVE_INFINITY;
		for(int z = 0; z < nbStates(); z++){
			if(((OpdfMultiMixtureGaussian)getOpdf(z)).mean()[lowcoverage ? 1 : 2] > coverage){
				coverage = ((OpdfMultiMixtureGaussian)getOpdf(z)).mean()[lowcoverage ? 1 : 2];
				methyState = z;
			}
		}
		
		return methyState;
	}
	
	/**
	 * Gives a description of this HMM.
	 * 
	 * @param nf A number formatter used to print numbers (e.g. Aij values).
	 * @return A textual description of this HMM.
	 */
	public String toString(NumberFormat nf)
	{
		//String s = "HMM with " + nbStates() + " state(s)\n" + pi.keySet().toArray(new Integer[pi.keySet().size()])[pi.keySet().size()/2] + "\n";
		String s = "HMM with " + nbStates() + " state(s)\n";
		//s += "  Poisson lambda: " + poisson.getMean() + "\n";
		
		for (int i = 0; i < nbStates(); i++) {
			s += "\nState " + i + "\n";
			//s += "  Pi: " + getPri(1,i) + "\n";
			s += "  Pi: " + nf.format(getPri(pi.keySet().toArray(new Integer[pi.keySet().size()])[pi.keySet().size()/2],i)) + " (" + nf.format(getPri(pi.firstKey(),i))+ ", " + nf.format(getPri(pi.lastKey(),i)) + ") " + "\n";
			
			
			s += "  Aij:";
			
			for (int j = 0; j < nbStates(); j++){
			//	s += " " + nf.format(getArij(1,i,j));
				s += " " + nf.format(getArij(pi.keySet().toArray(new Integer[pi.keySet().size()])[pi.keySet().size()/2],i,j)) + " (" + nf.format(getArij(pi.firstKey(),i,j)) + ", " + nf.format(getArij(pi.lastKey(),i,j)) + ") ";
			}
				
			s += "\n";
			
			s += "  Opdf: " + ((Opdf<O>) getOpdf(i)).toString(nf) + "\n";
		}
			
		return s;
	}
	
	public String toString()
	{
		return toString(NumberFormat.getInstance());
	}
	
	public BayesianNhmmV5<O> clone()throws CloneNotSupportedException{

				BayesianNhmmV5<O> hmm = new BayesianNhmmV5<O>(nbStates(), nbCpgDistStates, bayesianFactor);
				hmm.pi = (TreeMap<Integer, Double[]>) pi.clone();
				for (Integer key : a.keySet()){
						hmm.pi.put(key, pi.get(key).clone());
					
				}
				
				hmm.a = (TreeMap<Integer, Double[][]>) a.clone();
				for (Integer key : a.keySet()){
					Double[][] t = hmm.a.get(key);
					for (int i = 0; i < a.get(key).length; i++){
						
						t[i] = a.get(key)[i].clone();
					}
					hmm.a.put(key, t);
				}
					
				
				for (int i = 0; i < hmm.opdfs.size(); i++)
					hmm.opdfs.set(i, opdfs.get(i).clone());
				hmm.setMaxCpgNum(this.maxCpgInFrag);
				hmm.setMinCpgNum(this.minCpgInFrag);
				return hmm;
	}

}
