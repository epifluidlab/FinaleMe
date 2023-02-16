package org.cchmc.epifluidlab.methhmm.hmm;

import java.util.*;

import org.apache.commons.math3.distribution.MixtureMultivariateNormalDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.fitting.MultivariateNormalMixtureExpectationMaximization;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.Pair;

import org.cchmc.epifluidlab.methhmm.utils.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.*;


/**
 * An implementation of the K-Means learning algorithm.
 */
public class GMMLearner
{	
	private ClustersGMM clusters;
	private int nbStates;
	private List<List<org.cchmc.epifluidlab.methhmm.utils.ObservationVector>> obsSeqs;
	private ArrayList<HashMap<Integer, Pair<Integer, Double>>> nbCpgDistStates;
	private List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<org.cchmc.epifluidlab.methhmm.utils.ObservationVector>>> matrixObj;
	List<org.cchmc.epifluidlab.methhmm.utils.ObservationVector> observations;
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
	public GMMLearner(int nbStates, OpdfMultiMixtureGaussianFactory opdfFactory, List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<org.cchmc.epifluidlab.methhmm.utils.ObservationVector>>> matrixObj, int nbCpgDist, int features, ArrayList<Integer> mixNumberInFeature, double bayesianFactor,
                      MersenneTwister randomEngine, double tolKmeans, double decayKmeans, double maxCpgNum, double minCpgNum, boolean lowCov)
	{	


		this.opdfFactory = opdfFactory;
		this.nbStates = nbStates;
		this.nbCpgDist = nbCpgDist;
		this.features =features;
		this.mixNumberInFeature = mixNumberInFeature;
		this.bayesianFactor = bayesianFactor;
		this.obsSeqs = new ArrayList<List<org.cchmc.epifluidlab.methhmm.utils.ObservationVector>>();
		this.matrixObj = matrixObj;
		this.randomEngine = randomEngine;
		this.tol = tolKmeans;
		this.decayRate = decayKmeans;
		this.maxCpgNum = maxCpgNum;
		this.minCpgNum = minCpgNum;
		this.lowCov = lowCov;
		this.nbCpgDistStates = new ArrayList<HashMap<Integer, Pair<Integer, Double>>>();
		
		for( Pair<HashMap<Integer, Pair<Integer, Double>>, List<org.cchmc.epifluidlab.methhmm.utils.ObservationVector>> row : matrixObj){
			List<org.cchmc.epifluidlab.methhmm.utils.ObservationVector> tmp = new ArrayList<org.cchmc.epifluidlab.methhmm.utils.ObservationVector>();
			for(int j = 0; j < row.getSecond().size();j++){
				tmp.add(row.getSecond().get(j));
			}
			obsSeqs.add(tmp);
			nbCpgDistStates.add(row.getFirst());
			
		}
		this.observations = flat(obsSeqs);
		clusters = new ClustersGMM(nbStates, observations, features);
	}
	
	

	
	public BayesianNhmmV5<org.cchmc.epifluidlab.methhmm.utils.ObservationVector> iterate()
	{	
		BayesianNhmmV5<org.cchmc.epifluidlab.methhmm.utils.ObservationVector> hmm = new BayesianNhmmV5(nbStates,nbCpgDist, (OpdfFactory<? extends Opdf<org.cchmc.epifluidlab.methhmm.utils.ObservationVector>>) new OpdfMultiMixtureGaussianFactory(features, mixNumberInFeature), bayesianFactor);
		
		hmm.setMaxCpgNum(maxCpgNum);
		hmm.setMinCpgNum(minCpgNum);
		learnPi(hmm);
		learnAij(hmm);
		learnOpdf(hmm);
		if(this.methyState < 0){
			hmm.setBayesianFactor(bayesianFactor);
			this.methyState = getMethyState(hmm);
					System.out.println("Methy states here is:" + this.methyState);
					hmm.setMethyState(this.methyState);
					
		}

		terminated = true;
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
	public BayesianNhmmV5<org.cchmc.epifluidlab.methhmm.utils.ObservationVector> learn()
	{	
		BayesianNhmmV5<org.cchmc.epifluidlab.methhmm.utils.ObservationVector> hmm = null;
		BayesianNhmmV5<org.cchmc.epifluidlab.methhmm.utils.ObservationVector> prevHmm = null;
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
			//System.out.println("k-means iterations");
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
		this.methyState = getMethyState(hmm);
		System.out.println("Methystates final is:" + this.methyState);
		return hmm;
	}
	

	
	public int getMethyState(BayesianNhmmV5<org.cchmc.epifluidlab.methhmm.utils.ObservationVector> hmm){
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
		for (List<org.cchmc.epifluidlab.methhmm.utils.ObservationVector> sequence : obsSeqs){
			ArrayList<Double> stateWeight = clusters.clusterHash.get(sequence.get(0));
			int cpgDist = nbCpgDistStates.get(z).get(0).getFirst();
			Double[] piTmp = pi.get(cpgDist);
			for(int j = 0; j < hmm.nbStates(); j++){
				piTmp[j]+=stateWeight.get(j);
			}
			pi.put(cpgDist, piTmp);
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
		for (List<org.cchmc.epifluidlab.methhmm.utils.ObservationVector> obsSeq : obsSeqs) {
			if (obsSeq.size() < 2){
				z++;
				continue;
			}
				
			HashMap<Integer, Pair<Integer, Double>> cpgDistPairs = nbCpgDistStates.get(z);
			ArrayList<Double> first_state;
			ArrayList<Double> second_state = clusters.clusterHash.get(obsSeq.get(0));
			for (int i = 1; i < obsSeq.size(); i++) {
				first_state = second_state;
				second_state = clusters.clusterHash.get(obsSeq.get(i));
				int cpgDist = cpgDistPairs.get(i).getFirst();
				for(int x = 0; x < nbStates; x++){
					for(int y = 0; y < nbStates; y++){
						hmm.setArij(cpgDist, x, y,
								hmm.getArij(cpgDist, x, y)+first_state.get(x)*second_state.get(y));
					}
				}
				
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
	
	
	private void learnOpdf(BayesianNhmmV5<org.cchmc.epifluidlab.methhmm.utils.ObservationVector> hmm)
	{
		for (int i = 0; i < hmm.nbStates(); i++) {
			double sum = 0;
			for(ObservationVector o : observations){
				sum += clusters.clusters.get(i).get(o);
			}
			double[] weight = new double[observations.size()];
			for(int j = 0; j < observations.size(); j++){
				weight[j] = clusters.clusters.get(i).get(observations.get(j))/sum;
			}
			
			hmm.getOpdf(i).fit(observations, weight);
		}
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
class ClustersGMM
{	
	
	
	//private ArrayList<Collection<ObservationVector>> clusters;
	public ArrayList<Hashtable<org.cchmc.epifluidlab.methhmm.utils.ObservationVector,Double>> clusters;
	public Hashtable<org.cchmc.epifluidlab.methhmm.utils.ObservationVector,ArrayList<Double>> clusterHash;
	
	public ClustersGMM(int k, List<org.cchmc.epifluidlab.methhmm.utils.ObservationVector> observations, int feature)
	{
		clusters = new ArrayList<Hashtable<org.cchmc.epifluidlab.methhmm.utils.ObservationVector,Double>>();
		clusterHash = new Hashtable<org.cchmc.epifluidlab.methhmm.utils.ObservationVector,ArrayList<Double>>();
		
		double[][] data = new double[observations.size()][feature];
		for(int i = 0; i < observations.size(); i++){
			for(int j = 0; j < feature; j++){
				data[i][j] = observations.get(i).value(j);
			}
			
		}
		MixtureMultivariateNormalDistribution distribution = MultivariateNormalMixtureExpectationMaximization.estimate(data, k);
		
		for (int i = 0; i < k; i++) {
			Hashtable<org.cchmc.epifluidlab.methhmm.utils.ObservationVector,Double> tmp = new Hashtable<org.cchmc.epifluidlab.methhmm.utils.ObservationVector,Double>();
			clusters.add(tmp);
		}
		
		for (int i = 0; i < observations.size(); i++){
			org.cchmc.epifluidlab.methhmm.utils.ObservationVector element = observations.get(i);
			double[] row = data[i];
			double sum = 0;
			for (int j = 0; j < k; j++) {
				MultivariateNormalDistribution d = distribution.getComponents().get(j).getSecond();
				sum += d.density(row);
			}
			ArrayList<Double> list = new ArrayList<Double>();
			for (int j = 0; j < k; j++) {
				MultivariateNormalDistribution d = distribution.getComponents().get(j).getSecond();
				double v = d.density(row)/sum;
				Hashtable<org.cchmc.epifluidlab.methhmm.utils.ObservationVector,Double> tmp = clusters.get(j);
				tmp.put(element, v);
				clusters.set(j,tmp);
				list.add(v);
			}
			clusterHash.put(element, list);
		}
		
		
	}
	

	
	
	public double clusterWeight(Observation o, int clusterNb)
	{
		return clusterHash.get(o).get(clusterNb);
	}
	
	
	public Collection<org.cchmc.epifluidlab.methhmm.utils.ObservationVector> cluster(int clusterNb)
	{
		return clusters.get(clusterNb).keySet();
	}
	
}