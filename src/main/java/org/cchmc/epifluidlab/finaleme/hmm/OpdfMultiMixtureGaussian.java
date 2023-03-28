/* jahmm package - v0.6.1 */

/*
  *  Copyright (c) 2004-2006, Jean-Marc Francois.
 *
 *  This file is part of Jahmm.
 *  Jahmm is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Jahmm is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Jahmm; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 */

package org.cchmc.epifluidlab.finaleme.hmm;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import org.apache.commons.lang.ArrayUtils;

import be.ac.ulg.montefiore.run.distributions.GaussianDistribution;
import org.cchmc.epifluidlab.finaleme.utils.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.Opdf;


/**
 * This class represents a multivariate gaussian distribution function.
 */
public class OpdfMultiMixtureGaussian
implements Opdf<ObservationVector>
{	
	private MultiMixtureGaussianDistribution distribution;
	
	
	/**
	 * Builds a new gaussian probability distribution with zero mean and
	 * identity covariance matrix.
	 *
	 * @param dimension The dimension of the vectors.
	 */
	
	public OpdfMultiMixtureGaussian(int dimension, ArrayList<Integer> mixNumberInFeature)
	{
		distribution = new MultiMixtureGaussianDistribution(dimension,mixNumberInFeature);
	}
	
	public OpdfMultiMixtureGaussian(int dimension, ArrayList<Integer> mixNumberInFeature, double[] mean, double[] variance)
	{
		distribution = new MultiMixtureGaussianDistribution(dimension,mixNumberInFeature, mean, variance);
	}

	
	/**
	 * Builds a new gaussian probability distribution with a given mean and
	 * covariance matrix.
	 *
	 * @param mean The distribution's mean.
	 * @param covariance The distribution's covariance matrix.
	 */
	public OpdfMultiMixtureGaussian(double[] mean, double[][] covariance, ArrayList<Integer> mixNumberInFeature, ArrayList<ArrayList<Double>> meanInEachGaussian, 
			ArrayList<ArrayList<Double>> varianceInEachGaussian, ArrayList<ArrayList<Double>> propInEachGaussian)
	{		
		if (covariance.length == 0 || mean.length != covariance.length ||
				covariance.length != covariance[0].length)
			throw new IllegalArgumentException();
		
		distribution = new MultiMixtureGaussianDistribution(mean, covariance, mixNumberInFeature, meanInEachGaussian, varianceInEachGaussian, propInEachGaussian);
	}
	
	
	/**
	 * Returns (a copy of) this distribution's mean vector.
	 *
	 * @return The mean vector.
	 */
	public double[] mean()
	{
		return distribution.mean();
	}
	
	
	/**
	 * Returns (a copy of) this distribution's covariance matrix.
	 *
	 * @return The covariance matrix.
	 */
	public double[][] covariance()
	{
		return distribution.covariance();
	}
	
	
	/**
	 * Returns the dimension of the vectors handled by this distribution.
	 *
	 * @return The dimension of the vectors handled by this distribution.
	 */
	public int dimension()
	{
		return distribution.dimension();
	}
	
	public ArrayList<Integer> mixtureNumber()
	{
		return distribution.mixtureNumber();
	}
	
	
	public double probability(ObservationVector o)
	{
		if (o.dimension() != distribution.dimension())
			throw new IllegalArgumentException("Vector has a wrong " +
			"dimension");
		
		return distribution.probability(o.values());
	}
	
	/*
	public double probability(ObservationVector o, double bayesianFactor)
	{
		if (o.dimension() != distribution.dimension())
			throw new IllegalArgumentException("Vector has a wrong " +
			"dimension");
		
		return distribution.probability(o.values(), bayesianFactor);
	}
	*/
	public ObservationVector generate()
	{
		return new ObservationVector(distribution.generate());
	}
	
	
	public void fit(ObservationVector... oa)
	{
		fit(Arrays.asList(oa));
	}
	
	
	public void fit(Collection<? extends ObservationVector> co)
	{
		if (co.isEmpty())
			throw new IllegalArgumentException("Empty observation set");
		
		double[] weights = new double[co.size()];
		Arrays.fill(weights, 1. / co.size());
		
		fit(co, weights);
	}
	
	
	public void fit(ObservationVector[] o, double[] weights)
	{
		fit(Arrays.asList(o), weights);
	}
	
	public void fit(Collection<? extends ObservationVector> co, 
			double[] weights)
	{
		if (co.isEmpty() || co.size() != weights.length)
			throw new IllegalArgumentException();
		
		double[] mean = new double[dimension()];
		double[] variance = new double[dimension()];
		
		ObservationVector[] allO = co.toArray(new ObservationVector[co.size()]);
		
		ArrayList<ArrayList<Double>> meanInEachGaussian = new ArrayList<ArrayList<Double>>();
		ArrayList<ArrayList<Double>> varianceInEachGaussian = new ArrayList<ArrayList<Double>>();
		ArrayList<ArrayList<Double>> propInEachGaussian = new ArrayList<ArrayList<Double>>();
		
		ArrayList<Integer> mixtureNumber = mixtureNumber();
		for(int z = 0; z < dimension(); z++){
			ArrayList<Double> meanInEachGaussianTmp = new ArrayList<Double>();
			ArrayList<Double> varianceInEachGaussianTmp = new ArrayList<Double>();
			ArrayList<Double> propInEachGaussianTmp = new ArrayList<Double>();
			
			
			if(mixtureNumber.get(z) == 1){
				int i = 0;
				for (ObservationVector o : co){
					mean[z] += o.values()[z] * weights[i++];
				}
				
				
				i = 0;
				for (ObservationVector o : co){
					variance[z] += (o.values()[z]-mean[z])*(o.values()[z]-mean[z]) * weights[i++];
				}
				meanInEachGaussianTmp.add(mean[z]);
				varianceInEachGaussianTmp.add(variance[z]);
				propInEachGaussianTmp.add(1.0);
				
			}else{
				
				double[] o = new double[co.size()];
				for (int i = 0; i < co.size(); i++){
					o[i] = allO[i].value(z);
				}
				double[][] delta = getDelta(o, z);
				//System.err.println(o[0] + "\t" + delta[0][0]);
				double[] newMixingProportions = computeNewMixingProportions(delta, o, weights, z);
				double[] newMeans = computeNewMeans(delta, o, weights,z);
				double[] newVariances = computeNewVariances(delta, o, weights, z);
				for(int j = 0; j < mixtureNumber.get(z); j++){
					mean[z] += newMixingProportions[j]*newMeans[j];
					variance[z] += newMixingProportions[j]*newVariances[j];
				}
				
				//System.err.println(newMeans[0] + "\t" + newMeans[1]);
				meanInEachGaussianTmp.addAll(Arrays.asList(ArrayUtils.toObject(newMeans)));
				varianceInEachGaussianTmp.addAll(Arrays.asList(ArrayUtils.toObject(newVariances)));
				propInEachGaussianTmp.addAll(Arrays.asList(ArrayUtils.toObject(newMixingProportions)));
			}
			meanInEachGaussian.add(meanInEachGaussianTmp);
			varianceInEachGaussian.add(varianceInEachGaussianTmp);
			propInEachGaussian.add(propInEachGaussianTmp);
		}
		
		double[][] covariance = new double[dimension()][dimension()];
		int i = 0;
		for (ObservationVector o : co) {
			double[] obs = o.values();
			double[] omm = new double[obs.length];
			
			for (int j = 0; j < obs.length; j++)
				omm[j] = obs[j] - mean[j];
			
			for (int r = 0; r < dimension(); r++){
				for (int c = 0; c < dimension(); c++){
					covariance[r][c] += omm[r] * omm[c] * weights[i];
					
				}
					
			}
				
			
			i++;
		}
		distribution =  new MultiMixtureGaussianDistribution(mean, covariance, mixtureNumber, meanInEachGaussian, varianceInEachGaussian, propInEachGaussian);
		
	}
	
	private double[][] getDelta(double[] o, int z)
	{
		double[][] delta = new double[distribution.nbGaussians(z)][o.length];
		
		for (int i = 0; i < distribution.nbGaussians(z); i++) {
			ArrayList<Double> proportions = distribution.proportions(z);
			ArrayList<GaussianDistribution> distributions =
				distribution.distributions(z);
			
			for (int t = 0; t < o.length; t++){
				delta[i][t] = proportions.get(i) *
						distributions.get(i).probability(o[t]) / probability(distributions, proportions, o[t]);
				if(Double.isNaN(delta[i][t])){
					delta[i][t] = 0.0;
					//System.err.println(o[t] + "\t" + proportions.get(i) + "\t" + distributions.get(i).probability(o[t]) + "\t" + probability(distributions, proportions, o[t]) + "\t" + distributions.get(i).mean() + "\t" + distributions.get(i).variance());
				}
			}
				
		}
			
		return delta;
	}
	
	public double probability(ArrayList<GaussianDistribution> distributions, ArrayList<Double> proportions, double n)
	{
		double sum = 0.;
		
		for (int i = 0; i < distributions.size(); i++)
			sum += distributions.get(i).probability(n) * proportions.get(i);
		
		return sum;
	}

	private double[] computeNewMixingProportions(double[][] delta, 
			double[] o, double[] weights, int z)
	{
		double[] num = new double[distribution.nbGaussians(z)];
		double sum = 0.0;
		
		Arrays.fill(num, 0.0);
		
		for (int i = 0; i < distribution.nbGaussians(z); i++){
			for (int t = 0; t < weights.length; t++) {
				num[i] += weights[t] * delta[i][t];
				sum += weights[t] * delta[i][t];
			}
			if(Double.compare(num[i], 0.0) == 0){
				num[i] += 0.0001;
				sum += 0.0001;
			}
			
		}
		
		double[] newMixingProportions = new double[distribution.nbGaussians(z)];
		for (int i = 0; i < distribution.nbGaussians(z); i++) {
			newMixingProportions[i] = num[i]/sum;
			
		}
			
		
		return newMixingProportions;
	}
	
	private double[] computeNewMeans(double[][] delta, double[] o,
			double[] weights, int z)
	{
		double[] num = new double[distribution.nbGaussians(z)];
		double[] sum = new double[distribution.nbGaussians(z)];
		
		Arrays.fill(num, 0.0);
		Arrays.fill(sum, 0.0);
		
		for (int i = 0; i < distribution.nbGaussians(z); i++)
			for (int t = 0; t < o.length; t++) {
				num[i] += weights[t] * delta[i][t] * o[t];
				sum[i] += weights[t] * delta[i][t];
				//System.err.println(weights[t] * delta[i][t] * o[t] + "\t" + weights[t] + "\t" + o[t] + "\t" + delta[i][t]);
				//System.err.println(num[1] + "\t" + sum[1]);
			}
		//System.err.println(num[0] + "\t" + sum[0]);
		//System.err.println(num[1] + "\t" + sum[1]);
		double[] newMeans = new double[distribution.nbGaussians(z)];
		for (int i = 0; i < distribution.nbGaussians(z); i++)
			newMeans[i] = num[i] / sum[i];
		
		return newMeans;
	}
	
	private double[] computeNewVariances(double[][] delta, double[] o,
			double[] weights, int z)
	{
		double[] num = new double[distribution.nbGaussians(z)];
		double[] sum = new double[distribution.nbGaussians(z)];
		
		Arrays.fill(num, 0.);
		Arrays.fill(sum, 0.);
		
		for (int i = 0; i < distribution.nbGaussians(z); i++) {
			ArrayList<GaussianDistribution> distributions = distribution.distributions(z);
			
			for (int t = 0; t < o.length; t++) {
				num[i] += weights[t] * delta[i][t] *
				(o[t] - distributions.get(i).mean()) *
				(o[t] - distributions.get(i).mean());
				sum[i] += weights[t] * delta[i][t];
			}
		}
		
		double[] newVariances = new double[distribution.nbGaussians(z)];
		for (int i = 0; i < distribution.nbGaussians(z); i++) 
			newVariances[i] = num[i] / sum[i];
		
		return newVariances;
	}

	

	
	public OpdfMultiMixtureGaussian clone()
	{
		try {
			return (OpdfMultiMixtureGaussian) super.clone();
		} catch(CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
	}
	
	
	public String toString()
	{
		return toString(NumberFormat.getInstance());
	}
	
	
	public String toString(NumberFormat numberFormat)
	{
		String s = "Multi-variate Mixture Gaussian distribution --- Mean: [ ";
		double[] mean = distribution.mean();
		
		for (int i = 0; i < mean.length; i++){
			int mixed = mixtureNumber().get(i);
			if(mixed == 1){
				s += numberFormat.format(mean[i]) + " ";
			}else{
				s += numberFormat.format(mean[i]) + " (";
				for(int z = 0; z < mixed; z++){
					s +=  numberFormat.format(distribution.proportions(i).get(z)) + ":" + numberFormat.format(distribution.distributions(i).get(z).mean()) + ", ";
					//s +=  distribution.proportions(i).get(z) + ":" + distribution.distributions(i).get(z).mean() + ", ";
				}
				s += ") ";
			}
			
		}
			
		
		return s + "]";
	}


	private static final long serialVersionUID = 1L;
}
