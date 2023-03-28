package org.cchmc.epifluidlab.finaleme.hmm;

import java.util.*;

import org.apache.commons.math3.ml.clustering.CentroidCluster;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;
import org.apache.commons.math3.ml.clustering.MultiKMeansPlusPlusClusterer;
import org.apache.commons.math3.ml.distance.EuclideanDistance;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.Pair;

import org.cchmc.epifluidlab.finaleme.utils.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.*;


/**
 * An implementation of the K-Means learning algorithm.
 */
public class KMeansPlusLearner
{	
	private ClustersPlus clusters;
	private int nbStates;
	private List<List<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>> obsSeqs;
	private ArrayList<HashMap<Integer, Pair<Integer, Double>>> nbCpgDistStates;
	private List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>>> matrixObj;
	private OpdfMultiMixtureGaussianFactory opdfFactory;
	private int nbCpgDist;
	private int features;
	private double bayesianFactor;
	private ArrayList<Integer> mixNumberInFeature;
	private MersenneTwister randomEngine;
	private int methyState = -1;
	private boolean terminated = false;
	private double tol= 0.005;
	private int iteration = 20;
	private double decayRate = 0.01;
	private double minCpgNum;
	private double maxCpgNum;
	private boolean lowCov = false;

	
	/**
	 * Initializes a K-Means algorithm implementation.  This algorithm
	 * finds a HMM that models a set of observation sequences.
	 *
	 * @param nbStates  The number of states the resulting HMM will be made of.
	 * @param opdfFactory A class that builds the observation probability
	 *                    distributions associated to the states of the HMM.
	 * @param sequences A vector of observation sequences.  Each observation
	 *                sequences is a vector of
	 *                {@link be.ac.ulg.montefiore.run.jahmm.Observation
	 *                observations} compatible with the 
	 *                {@link be.ac.ulg.montefiore.run.jahmm.CentroidFactory
	 *                k-means algorithm}.
	 */
	public KMeansPlusLearner(int nbStates, OpdfMultiMixtureGaussianFactory opdfFactory, List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>>> matrixObj, int nbCpgDist, int features, ArrayList<Integer> mixNumberInFeature, double bayesianFactor,
                             MersenneTwister randomEngine, double tolKmeans, double decayKmeans, double maxCpgNum, double minCpgNum, boolean lowCov)
	{	


		this.opdfFactory = opdfFactory;
		this.nbStates = nbStates;
		this.nbCpgDist = nbCpgDist;
		this.features =features;
		this.mixNumberInFeature = mixNumberInFeature;
		this.bayesianFactor = bayesianFactor;
		this.obsSeqs = new ArrayList<List<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>>();
		this.matrixObj = matrixObj;
		this.randomEngine = randomEngine;
		this.tol = tolKmeans;
		this.decayRate = decayKmeans;
		this.maxCpgNum = maxCpgNum;
		this.minCpgNum = minCpgNum;
		this.lowCov = lowCov;
		this.nbCpgDistStates = new ArrayList<HashMap<Integer, Pair<Integer, Double>>>();
		for( Pair<HashMap<Integer, Pair<Integer, Double>>, List<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>> row : matrixObj){
			List<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> tmp = new ArrayList<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>();
			for(int j = 0; j < row.getSecond().size();j++){
				tmp.add(row.getSecond().get(j));
			}
			obsSeqs.add(tmp);
			nbCpgDistStates.add(row.getFirst());
			
		}
		List<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> observations = flat(obsSeqs);
		clusters = new ClustersPlus(nbStates, observations, randomEngine);
	}
	
	

	
	public BayesianNhmmV5<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> iterate()
	{	
		BayesianNhmmV5<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> hmm = new BayesianNhmmV5(nbStates,nbCpgDist, (OpdfFactory<? extends Opdf<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>>) new OpdfMultiMixtureGaussianFactory(features, mixNumberInFeature), bayesianFactor);
		
		hmm.setMaxCpgNum(maxCpgNum);
		hmm.setMinCpgNum(minCpgNum);
		learnPi(hmm);
		learnAij(hmm);
		learnOpdf(hmm);
		if(this.methyState < 0){
			//this.methyState = getMethyState(hmm, false);
			hmm.setBayesianFactor(bayesianFactor);
			this.methyState = getMethyState(hmm);
			//this.methyState = getMethyState();
			System.out.println(getMethyState(hmm, false) + "\t" + getMethyState(hmm, true) + "\t" + getMethyState(hmm) + "\t" + getMethyState());
					System.out.println("Methy states here is:" + this.methyState);
					hmm.setMethyState(this.methyState);
					
		}
		 
		
		terminated = optimizeCluster(hmm);
		//terminated = true;
		return hmm;
	}
	
	
	/**
	 * Returns <code>true</code> if the algorithm has reached a fix point,
	 * else returns <code>false</code>.
	 */
	public boolean isTerminated()
	{
		return terminated;
	}
	
	
	/**
	 * Does iterations of the K-Means algorithm until a fix point is reached.
	 * 
	 * @return The HMM that best matches the set of observation sequences given
	 *         (according to the K-Means algorithm).
	 */
	public BayesianNhmmV5<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> learn()
	{	
		BayesianNhmmV5<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> hmm = null;
		BayesianNhmmV5<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> prevHmm = null;
		KullbackLeiblerDistanceBayesianNhmmV5Calculator klc = 
				new KullbackLeiblerDistanceBayesianNhmmV5Calculator(matrixObj);
		int i = 0;
		double distance = Double.MAX_VALUE;
		double distancePre = 0.01;
		do{
			if(hmm != null){
				try {
					prevHmm = hmm.clone();
				} catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			System.out.println("k-means iterations");
			hmm = iterate();
			System.out.println(hmm);
			if(prevHmm != null){
				distancePre = distance;
				distance = klc.distance(prevHmm, hmm, true);
				System.out.println("Distance in Kmeans at iteration " + i + ": " +
						distance + "\t" + Math.abs((Math.abs(distance)-Math.abs(distancePre))/distancePre));
				if(Math.abs((Math.abs(distance)-Math.abs(distancePre))/distancePre) < decayRate || Math.abs(distance) < tol || i >= iteration){
					terminated = true;
				}
			}
			i++;
		}while(!isTerminated());
		//this.methyState = getMethyState(hmm, false);
		this.methyState = getMethyState(hmm);
		//this.methyState = getMethyState(hmm, true);
		//this.methyState = getMethyState();
		System.out.println(getMethyState(hmm, false) + "\t" + getMethyState(hmm, true) + "\t" + getMethyState(hmm) + "\t" + getMethyState());
		System.out.println("Methystates final is:" + this.methyState);
		return hmm;
	}
	
	//return the states number that represent methylated status, need to be called after learn()
	public int getMethyState(){
		HashMap<Integer, Pair<Double, Double>> stateMethy = new HashMap<Integer, Pair<Double, Double>>();
		int z = 0;
		for (List<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> obsSeq : obsSeqs) {
			
				
			HashMap<Integer, Pair<Integer, Double>> cpgDistPairs = nbCpgDistStates.get(z);
			
			for (int i = 0; i < obsSeq.size(); i++) {
				int state = clusters.clusterNb(obsSeq.get(i));
				double cpgMethy = cpgDistPairs.get(i).getSecond();
				if(stateMethy.containsKey(state)){
					Pair<Double, Double> tmp = stateMethy.get(state);
					stateMethy.put(state, new Pair<Double, Double>(tmp.getFirst() + cpgMethy, tmp.getSecond() + 1.0));
				}else{
					stateMethy.put(state, new Pair<Double, Double>(cpgMethy, 1.0));
				}
			}
			z++;
		}
		double maxMethy = 0;
		int methyState = 0;
		for(int i : stateMethy.keySet()){
			double methy = stateMethy.get(i).getFirst()/stateMethy.get(i).getSecond();
			if(methy > maxMethy){
				maxMethy = methy;
				methyState = i;
			}
		}
		return methyState;
	}
	
	public int getMethyState(BayesianNhmmV5<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> hmm){
		int methyState = 0;
		double coverage = Double.NEGATIVE_INFINITY;
		for(int z = 0; z < hmm.nbStates(); z++){
			if(((OpdfMultiMixtureGaussian)hmm.getOpdf(z)).mean()[(lowCov ? 1 : 2)] > coverage){
				coverage = ((OpdfMultiMixtureGaussian)hmm.getOpdf(z)).mean()[(lowCov ? 1 : 2)];
				methyState = z;
			}
		}
		
		return methyState;
	}
	
	//based on which state is assigned to methylated state give the least number of state changes
	public int getMethyState(BayesianNhmmV5<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> hmm, boolean notUseMethyPrior){
		long numberLeastChanged = Long.MAX_VALUE;
		int methyState = 0;
	//	double bicSmallest = Double.MAX_VALUE;
		for(int z = 0; z < hmm.nbStates(); z++){
			long numberChanged = 0L;
		//	double likelihood = 0;
		//	double likelihoodWithMethy = 0;
		//	long number = 0L;
			for (Pair<HashMap<Integer, Pair<Integer, Double>>, List<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>> obsSeq : matrixObj) {

				ViterbiBayesianNhmmV5Calculator vc = new ViterbiBayesianNhmmV5Calculator(obsSeq, hmm, z, 0.0);
				//ViterbiBayesianNhmmV5Calculator vb = new ViterbiBayesianNhmmV5Calculator(obsSeq, hmm, z, 0.0,false);
				int states[] = vc.stateSequence();
			//	double prob = vb.lnProbability();
			//	if(Double.isNaN(prob) || Double.isInfinite(prob)){
					
			//	}else{
			///		likelihood += prob/states.length;
			//	}
			//	double probWithMethy = vb.lnProbability(true);
			//	if(Double.isNaN(probWithMethy) || Double.isInfinite(probWithMethy)){
					
			//	}else{
			//		likelihoodWithMethy += probWithMethy/states.length;
			//	}
				//System.err.println("optimize cluster: " + j);
				for (int i = 0; i < states.length; i++) {
					org.cchmc.epifluidlab.finaleme.utils.ObservationVector o = obsSeq.getSecond().get(i);
					//number++;
					if (clusters.clusterNb(o) != states[i]) {
						numberChanged++;
					}
				}
			}
		//	System.out.println("hmm state " + z + "\tLikelihood is:" + likelihood);	
		//	System.out.println("LikelihoodWithMethy is:" + likelihoodWithMethy);	
		//	double bic = Math.log(number)*(hmm.nbCpgDistState()*3+3+6)-likelihood;
		//	double bicWithMethy = Math.log(number)*(hmm.nbCpgDistState()*3+3+6)-likelihoodWithMethy;
		//	System.out.println("BIC is:" + bic);
		//	System.out.println("BICWithMethy is:" + bicWithMethy);
			if(numberChanged < numberLeastChanged){
				numberLeastChanged = numberChanged;
				methyState = z;
			}
			//if(bic < bicSmallest){
			//	methyState = z;
			//}
		
		}
		return methyState;
	}

	public int returnMethyState(){
		return this.methyState;
	}
	
	private void learnPi(BayesianNhmmV5<?> hmm)
	{	
		TreeMap<Integer, Double[]> pi = new TreeMap<Integer, Double[]>();
		for(int i = 0; i<= nbCpgDist; i++){
			Double[] piTmp = new Double[nbStates];
			for(int j = 0; j < nbStates; j++){
				piTmp[j] = 0.;
			}
			pi.put(i, piTmp);
		}
		
		int z = 0;
		for (List<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> sequence : obsSeqs){
			int state = clusters.clusterNb(sequence.get(0));
			int cpgDist = nbCpgDistStates.get(z).get(0).getFirst();
			Double[] piTmp = pi.get(cpgDist);
			piTmp[state]+=1;
			z++;
		}
		for(int i = 0; i<= nbCpgDist; i++){
				Double[] piTmp = pi.get(i);
				double sum = 0.;
				for(int j = 0; j < nbStates; j++){	
					sum += piTmp[j];
				}
				if(Double.compare(sum,0.0)==0){
					piTmp[0]=0.5;
					piTmp[1]=0.5;
				}else{
					piTmp[0]/=sum;
					piTmp[1]/=sum;
				}
				hmm.setPri (i, 0 ,piTmp[0]);
				hmm.setPri (i, 1 ,piTmp[1]);
		}
		
	}
	
	
	private void learnAij(BayesianNhmmV5<?> hmm)
	{	
		
		TreeMap<Integer, Double[][]>  arij = new TreeMap<Integer, Double[][]>();
		for(int r = 0; r <= nbCpgDist; r++){
			for (int i = 0; i < hmm.nbStates(); i++){
				for(int j = 0; j < hmm.nbStates(); j++){
					hmm.setArij(r, i, j, 0.);
				}
			}
		}
		
		int z = 0;
		for (List<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> obsSeq : obsSeqs) {
			if (obsSeq.size() < 2){
				z++;
				continue;
			}
				
			HashMap<Integer, Pair<Integer, Double>> cpgDistPairs = nbCpgDistStates.get(z);
			int first_state;
			int second_state = clusters.clusterNb(obsSeq.get(0));
			for (int i = 1; i < obsSeq.size(); i++) {
				first_state = second_state;
				second_state =
					clusters.clusterNb(obsSeq.get(i));
				int cpgDist = cpgDistPairs.get(i).getFirst();
				hmm.setArij(cpgDist, first_state, second_state,
						hmm.getArij(cpgDist, first_state, second_state)+1.);
			}
			z++;
		}
		
		/* Normalize Aij array */
		for(int r = 0; r <= nbCpgDist; r++){
			for (int i = 0; i < hmm.nbStates(); i++){
				double sum = 0;
				for(int j = 0; j < hmm.nbStates(); j++){
					sum += hmm.getArij(r, i, j);
					//System.err.println(r + "\t" + i + "\t" + j + "\t" + hmm.getArij(r, i, j));
				}
				if(Double.compare(sum,0.0)==0){
					for (int j = 0; j < hmm.nbStates(); j++) 
						hmm.setArij(r, i, j, 1. / hmm.nbStates());
				}else{
					for (int j = 0; j < hmm.nbStates(); j++) 
						hmm.setArij(r, i, j, hmm.getArij(r, i, j) / sum);
				}
				
			}
		}
		
	}
	
	
	private void learnOpdf(BayesianNhmmV5<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> hmm)
	{
		for (int i = 0; i < hmm.nbStates(); i++) {
			Collection<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> clusterObservations = clusters.cluster(i);
			
			if (clusterObservations.isEmpty())
				hmm.setOpdf(i, (Opdf<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>) opdfFactory.factor());
			else
				hmm.getOpdf(i).fit(clusterObservations);
		}
	}
	
	/* Return true if no modification */
	private boolean optimizeCluster(BayesianNhmmV5<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> hmm)
	{	
		boolean modif = false;
		int j = 0;
		for (Pair<HashMap<Integer, Pair<Integer, Double>>, List<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>> obsSeq : matrixObj) {

			ViterbiBayesianNhmmV5Calculator vc = new ViterbiBayesianNhmmV5Calculator(obsSeq, hmm, this.methyState, 0.0);
			int states[] = vc.stateSequence();
			//System.err.println("optimize cluster: " + j);
			for (int i = 0; i < states.length; i++) {
				ObservationVector o = obsSeq.getSecond().get(i);
				
				if (clusters.clusterNb(o) != states[i]) {
					modif = true;
					clusters.remove(o, clusters.clusterNb(o));
					clusters.put(o, states[i]);
				}
			}
			j++;
		}
		
		return !modif;
	}

	
	static <T> List<T> flat(List<? extends List<? extends T>> lists)
	{	
		List<T> v = new ArrayList<T>();
		
		for (List<? extends T> list : lists)
			v.addAll(list);
		
		return v;
	}
	
	
}


/*
 * This class holds the matching between observations and clusters.
 */
class ClustersPlus
{	
	class ValuePlus
	{
		private int clusterNb;
		
		ValuePlus(int clusterNb)
		{
			this.clusterNb = clusterNb;
		}
		
		void setClusterNb(int clusterNb)
		{
			this.clusterNb = clusterNb;
		}
		
		int getClusterNb()
		{
			return clusterNb;
		}
	}
	
	
	private Hashtable<org.cchmc.epifluidlab.finaleme.utils.ObservationVector,ValuePlus> clustersHash;
	//private ArrayList<Collection<ObservationVector>> clusters;
	private ArrayList<Hashtable<org.cchmc.epifluidlab.finaleme.utils.ObservationVector,Integer>> clusters;
	
	public ClustersPlus(int k, List<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> observations, MersenneTwister randomEngine)
	{
		
		clustersHash = new Hashtable<org.cchmc.epifluidlab.finaleme.utils.ObservationVector,ValuePlus>();
		clusters = new ArrayList<Hashtable<org.cchmc.epifluidlab.finaleme.utils.ObservationVector,Integer>>();

		KMeansPlusPlusClusterer<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> kmc = new KMeansPlusPlusClusterer<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>(k, 10000, new EuclideanDistance(), randomEngine);
		
		MultiKMeansPlusPlusClusterer<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> mkmc = new MultiKMeansPlusPlusClusterer<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>(kmc, 20);
		List<CentroidCluster<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>> clusterResults = mkmc.cluster( (Collection<org.cchmc.epifluidlab.finaleme.utils.ObservationVector>) observations);
		
		for (int i = 0; i < k; i++) {
			Collection<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> cluster = clusterResults.get(i).getPoints();
			
			Hashtable<org.cchmc.epifluidlab.finaleme.utils.ObservationVector,Integer> tmp = new Hashtable<org.cchmc.epifluidlab.finaleme.utils.ObservationVector,Integer>();
			for (org.cchmc.epifluidlab.finaleme.utils.ObservationVector element : cluster){
				clustersHash.put(element, new ValuePlus(i));
				tmp.put(element, i);
			}
			clusters.add(tmp);	
		}
	}
	
	
	public boolean isInCluster(Observation o, int clusterNb)
	{
		return clusterNb(o) == clusterNb;
	}
	
	
	public int clusterNb(Observation o)
	{
		return clustersHash.get(o).getClusterNb();
	}
	
	
	public Collection<org.cchmc.epifluidlab.finaleme.utils.ObservationVector> cluster(int clusterNb)
	{
		return clusters.get(clusterNb).keySet();
	}
	
	
	public void remove(Observation o, int clusterNb)
	{
		clustersHash.get(o).setClusterNb(-1);
		clusters.get(clusterNb).remove(o);
	}
	
	
	public void put(org.cchmc.epifluidlab.finaleme.utils.ObservationVector o, int clusterNb)
	{
		clustersHash.get(o).setClusterNb(clusterNb);
		clusters.get(clusterNb).put(o, clusterNb);
	}
}