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

package org.cchmc.epifluidlab.methhmm.hmm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import be.ac.ulg.montefiore.run.distributions.GaussianDistribution;
import be.ac.ulg.montefiore.run.distributions.MultiRandomDistribution;



/**
 * This class implements a multi-variate Gaussian distribution.
 */
public class MultiMixtureGaussianDistribution
implements MultiRandomDistribution
{
	final private int dimension;
	private ArrayList<Integer> mixNumberInFeature;
	private ArrayList<ArrayList<Double>> meanInEachGaussian;
	private ArrayList<ArrayList<Double>> varianceInEachGaussian;
	private ArrayList<ArrayList<Double>> propInEachGaussian;
	private ArrayList<ArrayList<GaussianDistribution>> GaussianDistributions;
	final private double[] mean;
	final private double[][] covariance;
	private double[][] covarianceL = null; // covariance' Cholesky decomposition
	private double[][] covarianceInv = null;
	private double covarianceDet;
	private final static Random randomGenerator = new Random();
	//private double bayesianFactor = 1;
	
	
	/**
	 * Creates a new pseudo-random, multivariate gaussian distribution.
	 *
	 * @param mean The mean vector of the generated numbers.  This array is
	 *             copied.
	 * @param covariance The covariance of the generated numbers.  This array
	 *                   is copied.  <code>covariance[r][c]</code> is the
	 *                   element at row <code>r</code> and column
	 *                   <code>c</code>.
	 */
	public MultiMixtureGaussianDistribution(double[] mean, double[][] covariance, ArrayList<Integer> mixNumberInFeature, ArrayList<ArrayList<Double>> meanInEachGaussian, 
			ArrayList<ArrayList<Double>> varianceInEachGaussian, ArrayList<ArrayList<Double>> propInEachGaussian)
	{	
		if (!SimpleMatrix.isSquare(covariance))
			throw new IllegalArgumentException("Covariance must be a square " +
			"matrix");
		
		dimension = SimpleMatrix.nbRows(covariance);
		if (mean.length != dimension)
			throw new IllegalArgumentException("mean and covariance " +
			"dimensions don't match");
		
		this.mean = SimpleMatrix.vector(mean);
		this.covariance = SimpleMatrix.matrix(covariance);
		this.mixNumberInFeature = mixNumberInFeature;
		this.meanInEachGaussian = meanInEachGaussian;
		this.varianceInEachGaussian = varianceInEachGaussian;
		this.propInEachGaussian = propInEachGaussian;
		
		this.GaussianDistributions = new ArrayList<ArrayList<GaussianDistribution>>();
		for(int i = 0; i < dimension; i++){
			ArrayList<GaussianDistribution> GaussianDistributionsTmp = new ArrayList<GaussianDistribution>();
			for(int j = 0; j < mixNumberInFeature.get(i); j++){
				GaussianDistributionsTmp.add(new GaussianDistribution(meanInEachGaussian.get(i).get(j), varianceInEachGaussian.get(i).get(j)));
			}
			this.GaussianDistributions.add(GaussianDistributionsTmp);
		}
		
	}
	
	
	/**
	 * Creates a new pseudo-random, multivariate gaussian distribution with
	 * zero mean and identity covariance.
	 *
	 * @param dimension This distribution dimension.
	 */
	
	public MultiMixtureGaussianDistribution(int dimension, ArrayList<Integer> mixNumberInFeature)
	{
		if (dimension <= 0)
			throw new IllegalArgumentException();
		
		this.dimension = dimension;
		
		mean = SimpleMatrix.vector(dimension);
		covariance = SimpleMatrix.matrixIdentity(dimension);
		this.mixNumberInFeature = mixNumberInFeature;
		if(mixNumberInFeature.size() != dimension)
			throw new IllegalArgumentException();
		
		this.meanInEachGaussian = new ArrayList<ArrayList<Double>>();
		this.varianceInEachGaussian = new ArrayList<ArrayList<Double>>();
		this.propInEachGaussian = new ArrayList<ArrayList<Double>>();
		this.GaussianDistributions = new ArrayList<ArrayList<GaussianDistribution>>();
		for(int i = 0; i < dimension; i++){
			ArrayList<Double> meanInEachGaussianTmp = new ArrayList<Double>();
			ArrayList<Double> varianceInEachGaussianTmp = new ArrayList<Double>();
			ArrayList<Double> propInEachGaussianTmp = new ArrayList<Double>();
			if(mixNumberInFeature.get(i) < 1){
				throw new IllegalArgumentException();
			}else if(mixNumberInFeature.get(i) == 1){
				meanInEachGaussianTmp.add(mean[i]);
				varianceInEachGaussianTmp.add(covariance[i][i]);
				propInEachGaussianTmp.add(1.0);
			}else{
				for(int j = 0; j < mixNumberInFeature.get(i); j++){
					meanInEachGaussianTmp.add((1. + 2. * (double) j) / (2. * (double) mixNumberInFeature.get(i)));
					varianceInEachGaussianTmp.add(1.0);
					propInEachGaussianTmp.add(1. / ((double) mixNumberInFeature.get(i)));
				}
			}
			ArrayList<GaussianDistribution> GaussianDistributionsTmp = new ArrayList<GaussianDistribution>();
			for(int j = 0; j < mixNumberInFeature.get(i); j++){
				GaussianDistributionsTmp.add(new GaussianDistribution(meanInEachGaussianTmp.get(j), varianceInEachGaussianTmp.get(j)));
			}
			this.meanInEachGaussian.add(meanInEachGaussianTmp);
			this.varianceInEachGaussian.add(varianceInEachGaussianTmp);
			this.propInEachGaussian.add(propInEachGaussianTmp);
			this.GaussianDistributions.add(GaussianDistributionsTmp);
		}
		
	
	}
	
	
	public MultiMixtureGaussianDistribution(int dimension, ArrayList<Integer> mixNumberInFeature, double[] mean, double[] variance)
	{
		if (dimension <= 0)
			throw new IllegalArgumentException();
		
		this.dimension = dimension;
		
		this.mean = mean.clone();
		covariance = SimpleMatrix.matrix(dimension, dimension);
		for (int i = 0; i < dimension; i++)
			covariance[i][i] = variance[i];
		
		this.mixNumberInFeature = mixNumberInFeature;
		if(mixNumberInFeature.size() != dimension)
			throw new IllegalArgumentException();
		
		this.meanInEachGaussian = new ArrayList<ArrayList<Double>>();
		this.varianceInEachGaussian = new ArrayList<ArrayList<Double>>();
		this.propInEachGaussian = new ArrayList<ArrayList<Double>>();
		this.GaussianDistributions = new ArrayList<ArrayList<GaussianDistribution>>();
		for(int i = 0; i < dimension; i++){
			ArrayList<Double> meanInEachGaussianTmp = new ArrayList<Double>();
			ArrayList<Double> varianceInEachGaussianTmp = new ArrayList<Double>();
			ArrayList<Double> propInEachGaussianTmp = new ArrayList<Double>();
			if(mixNumberInFeature.get(i) < 1){
				throw new IllegalArgumentException();
			}else if(mixNumberInFeature.get(i) == 1){
				meanInEachGaussianTmp.add(mean[i]);
				varianceInEachGaussianTmp.add(covariance[i][i]);
				propInEachGaussianTmp.add(1.0);
			}else{
				for(int j = 0; j < mixNumberInFeature.get(i); j++){
					meanInEachGaussianTmp.add(mean[i] + 2. * (double) j);
					varianceInEachGaussianTmp.add(covariance[i][i]);
					//varianceInEachGaussianTmp.add(1.0);
					propInEachGaussianTmp.add(1. / ((double) mixNumberInFeature.get(i)));
				}
			}
			ArrayList<GaussianDistribution> GaussianDistributionsTmp = new ArrayList<GaussianDistribution>();
			for(int j = 0; j < mixNumberInFeature.get(i); j++){
				GaussianDistributionsTmp.add(new GaussianDistribution(meanInEachGaussianTmp.get(j), varianceInEachGaussianTmp.get(j)));
			}
			this.meanInEachGaussian.add(meanInEachGaussianTmp);
			this.varianceInEachGaussian.add(varianceInEachGaussianTmp);
			this.propInEachGaussian.add(propInEachGaussianTmp);
			this.GaussianDistributions.add(GaussianDistributionsTmp);
		}
		
	
	}
	
	public int dimension()
	{
		return dimension;
	}

	public ArrayList<Integer> mixtureNumber()
	{
		return this.mixNumberInFeature;
	}
	
	public ArrayList<Double> proportions(int i) 
	{
		return (ArrayList<Double>)propInEachGaussian.get(i).clone();
	}

	public ArrayList<GaussianDistribution> distributions(int i)
	{
		return (ArrayList<GaussianDistribution>)GaussianDistributions.get(i).clone();
	}

	public int nbGaussians(int i)
	{
		return this.mixNumberInFeature.get(i);
	}


	
	/**
	 * Returns (a copy of) this distribution's mean vector.
	 *
	 * @return This distribution's mean vector.
	 */
	public double[] mean()
	{
		return (double[]) mean.clone();
	}
	
	
	/**
	 * Returns (a copy of) this distribution's covariance matrix.
	 *
	 * @return This distribution's covariance matrix.
	 */
	public double[][] covariance()
	{
		return SimpleMatrix.matrix(covariance);
	}
	
	
	private double[][] covarianceL()
	{
		if (covarianceL == null) {
			covarianceL = SimpleMatrix.decomposeCholesky(covariance);
			covarianceDet = SimpleMatrix.determinantCholesky(covarianceL);
		}
		
		return covarianceL;
	}
	
	
	private double[][] covarianceInv()
	{
		if (covarianceInv == null)
			covarianceInv = SimpleMatrix.inverseCholesky(covarianceL());
		
		return covarianceInv;
	}
	
	
	/**
	 * Returns the covariance matrix determinant.
	 *
	 * @return The covariance matrix determinant.
	 */
	public double covarianceDet()
	{
		covarianceL();
		
		return covarianceDet;
	}
	
	
	/**
	 * Generates a pseudo-random vector according to this distribution.
	 * The vectors are generated using the Cholesky decomposition of the
	 * covariance matrix.
	 *
	 * @return A pseudo-random vector.
	 */
	public double[] generate()
	{
		double[] d = SimpleMatrix.vector(dimension);
		
		for (int i = 0; i < dimension; i++){
				if(mixNumberInFeature.get(i) == 1){
					d[i] = randomGenerator.nextGaussian();
				}else{
					double r = randomGenerator.nextDouble();
					double sum = 0.;	
					for (int z = 0; z < propInEachGaussian.get(i).size(); z++) {
						sum += propInEachGaussian.get(i).get(z);
						
						if (r <= sum){
							d[i] =  (GaussianDistributions.get(i).get(z).generate()-GaussianDistributions.get(i).get(z).mean())/Math.sqrt(GaussianDistributions.get(i).get(z).variance());
							break;
						}
							
						
					}
				}
		}
		
		
		
		return SimpleMatrix.plus(SimpleMatrix.times(covarianceL(), d), mean);
	}
	
	/*
	public double probability(double[] v)
	{
		if (v.length != dimension)
			throw new IllegalArgumentException("Argument array size is not " +
					"compatible with this distribution");
		
		double[][] vmm = SimpleMatrix.matrix(SimpleMatrix.minus(v, mean));
		
		double expArg =
			(SimpleMatrix.times(SimpleMatrix.transpose(vmm),
					SimpleMatrix.times(covarianceInv(), vmm))[0][0]) * -.5;
		
		return Math.exp(expArg) / 
		(Math.pow(2. * Math.PI, ((double) dimension) / 2.) * 
				Math.pow(covarianceDet(), .5)); 
	}
	*/
	
	public double probability(double[] v)
	{
		if (v.length != dimension)
			throw new IllegalArgumentException("Argument array size is not " +
					"compatible with this distribution");
		
		double sumAll = 1.;
		
		for (int i = 0; i < GaussianDistributions.size(); i++){
			ArrayList<GaussianDistribution> GaussianDistributionsTmp = GaussianDistributions.get(i);
			double value = v[i];
			double sum = 0.;
			ArrayList<Double> propInEachGaussianTmp = propInEachGaussian.get(i);
			for (int j = 0; j < GaussianDistributionsTmp.size(); j++){
				sum += GaussianDistributionsTmp.get(j).probability(value) * propInEachGaussianTmp.get(j);
				if(Double.isNaN(sum)){
					System.err.println(GaussianDistributionsTmp.get(j).mean() + "\t" + GaussianDistributionsTmp.get(j).variance() + "\t" + GaussianDistributionsTmp.get(j).probability(value) + "\t" + propInEachGaussianTmp.get(j));
					System.err.println("Need to reduce the number of mixture in feature " + (i+1));
					System.exit(1);
				}
			}
			sumAll *= sum;
			//if(i==0){
			//	sumAll += Math.log(sum)*bayesianFactor;
			//}else{
			//	sumAll += Math.log(sum);
			//}
			
		}
			
		
		return Math.exp(sumAll);
	}
	/*
	public double probability(double[] v, double bayesianFactor)
	{
		if (v.length != dimension)
			throw new IllegalArgumentException("Argument array size is not " +
					"compatible with this distribution");
		
		double sumAll = 0.;
		
		for (int i = 0; i < GaussianDistributions.size(); i++){
			ArrayList<GaussianDistribution> GaussianDistributionsTmp = GaussianDistributions.get(i);
			double value = v[i];
			double sum = 0.;
			ArrayList<Double> propInEachGaussianTmp = propInEachGaussian.get(i);
			for (int j = 0; j < GaussianDistributionsTmp.size(); j++){
				sum += GaussianDistributionsTmp.get(j).probability(value) * propInEachGaussianTmp.get(j);
				if(Double.isNaN(sum)){
					System.err.println(GaussianDistributionsTmp.get(j).mean() + "\t" + GaussianDistributionsTmp.get(j).variance() + "\t" + GaussianDistributionsTmp.get(j).probability(value) + "\t" + propInEachGaussianTmp.get(j));
					System.err.println("Need to reduce the number of mixture in feature " + (i+1));
					System.exit(1);
				}
			}
			if(i==0){
				sumAll += Math.log(sum)*bayesianFactor;
			}else{
				sumAll += Math.log(sum);
			}
			
		}
			
		
		return Math.exp(sumAll);
	}
	*/
	
	
	private static final long serialVersionUID = -2438571303843585271L;



}
