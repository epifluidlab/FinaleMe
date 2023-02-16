/**
 * CcNhmm.java
 * Apr 27, 2016
 * 10:25:30 AM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.methhmm.hmm;


import htsjdk.samtools.util.IntervalTree;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.Logger;
import org.cchmc.epifluidlab.methhmm.utils.ObservationVector;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;


import be.ac.ulg.montefiore.run.jahmm.io.FileFormatException;


/**
 * use cpg distance dependant transition probability matrix for the non homogenous HMM
 */
public class MethHMM {

	@Option(name="-states",usage="number of states in HMM model. ##current only allow even states number##. default: 2")
	public int states = 2;

	@Option(name="-features",usage="number of features for each observation, default: 3")
	public int features = 3;
	
	@Option(name="-tol",usage="tolerence level for the converge, default: 1e-4")
	public double tol = 1e-4;
	
	@Option(name="-decayRate",usage="distance changes less than decayRate, default: 0.01")
	public double decayRate = 0.01;

	@Option(name="-tolKmeans",usage="tolerence level for the K-means part, default: 0.005")
	public double tolKmeans = 0.005;

	@Option(name="-decayKmeans",usage="distance changes less than decayRate, default: 0.01")
	public double decayKmeans = 0.01;
	
	@Option(name="-bin",usage="base pair to bin to calculate non-homogenous Pri and Arij matrix, when number of data point small, need to increase the siz of bin. default: 1")
	public int bin = 1;

	@Option(name="-covOutlier",usage="the X sd from the mean of the coverage will be filtered out, <0 means no filter. default: -1")
	public double covOutlier = -1;

	
	@Option(name="-mixNumberInFeature",usage="the number of gaussian is mixed in each feature. could provide multiple option, default is no mixture for each of features, default: null")
	public ArrayList<Integer> mixNumberInFeature = null;
	
	@Option(name="-iteration",usage="maximum number of iteration for HMM model converge, default: 50")
	public int iteration = 50;

	@Option(name="-decodeP",usage="criteria used in decoding step, -1.0 to 1.0 to identify methylated points. default: 0")
	public double decodeP = 0.;
	
	@Option(name="-miniDataPoints",usage="minimum CpG data point per fragments, default: 1")
	public int miniDataPoints = 1;

	@Option(name="-maxFragLen",usage="maximum fragment length allowed for the model, also it represents the number of possible cpg distance. default: 500")
	public int maxFragLen = 500;
	
	@Option(name="-bayesianFactor",usage="bayesian factor used to weight prior. default: 0.9")
	public double bayesianFactor = 0.9;

	@Option(name="-minFragLen",usage="minimum fragment length allowed for the model. default: 30")
	public int minFragLen = 30;

	@Option(name="-maxCpgDist",usage="maximum cpg distance allowed for the model. default: 500")
	public int maxCpgDist = 500;
	
	@Option(name="-maxCpgs",usage="maximum number of Cpg. default: 1000")
	public int maxCpgs = 1000;

	@Option(name="-methylatedState",usage="which states in HMM model indicate the methylated states, only allow 0 or 1. default: 1")
	public int methylatedState = 1;

//	@Option(name="-notAutomateIdentifyMethyState",usage="by default, the program will mask -methylatedState and identify the methylation state by using the most abundant state. default: false")
//	public boolean notAutomateIdentifyMethyState = false;


	@Option(name="-lowCoverage",usage="low coverage mode, differnt filter for extream coverage region. default: false")
	public boolean lowCoverage = false;

	@Option(name="-randomPerm",usage="use the prior methylation to randomly assign m or um point, instead of trained HMM. default: false")
	public boolean randomPerm = false;

	@Option(name="-aucMode",usage="calcualte auc in each threshold. default: false")
	public boolean aucMode = false;
	
	@Option(name="-decodeModeOnly",usage="only decode, no training. default: false")
	public boolean decodeModeOnly = false;
	
	@Option(name="-wgbs",usage="if the input is from wgbs, then use pre-estimated Opdf and Adj will reduce the training time and converge much faster. default: false")
	public boolean wgbs = false;

	@Option(name="-gmm",usage="use GMM to initialize the hmm model (not good speed for large number of data points...). default: false")
	public boolean gmm = false;

	@Option(name="-cpgNumClip",usage="maximum number of Cpg used to scale the number of cpg in the fragment in HMM. -1 means use maximum cpg in each dataset itself. default: 500")
	public int cpgNumClip = 500;

	

	@Option(name="-seed",usage="seed for randomness.when less than 0, it will be random and not repeatable. default: 2017")
	public int seed = 2017;
	
	@Option(name="-region",usage="only look at data points within these regions. need to be bed format. default: null")
	public String region = null;

	@Option(name="-exclude",usage="exclude data points within these regions. need to be bed format. default: null")
	public String exclude = null;

	
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "CcNhmm [opts] model input_matrix.txt[.gz] [prediction.txt.gz]";
	
	private static Logger log = Logger.getLogger(MethHMM.class);

	private static long startTime = -1;
	private static long points = 0;
	private MersenneTwister randomEngine;
	private double maxCpgNum = Double.NEGATIVE_INFINITY;
	//private double minCpgNum = Double.POSITIVE_INFINITY;
	//private int methylatedState = 1;



	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		MethHMM cnh = new MethHMM();
		//BasicConfigurator.configure();
		cnh.doMain(args);
	}
	

	public void doMain(String[] args)
			throws Exception {

					CmdLineParser parser = new CmdLineParser(this);
					//parser.setUsageWidth(80);
					try
					{
						if(help || args.length < 2) throw new CmdLineException(parser, USAGE, new Throwable());
						parser.parseArgument(args);
						
					
					}
					catch (CmdLineException e)
					{
						System.err.println(e.getMessage());
						// print the list of available options
						parser.printUsage(System.err);
						System.err.println();
						return;
					}

					//read input bed file, for each row,
					//String intervalFile = arguments.get(0);
					String modelFile = arguments.get(0);
					String inputFile = arguments.get(1);
					String outputFile = arguments.get(2);
					//if(states % 2 != 0){
					//	throw new IllegalArgumentException("states number could not be odd asnymore");
					//}
					initiate();
					//System.out.println("new3:\n" );
					//input matrix processing
					
					
					//test cpg distance distributipn
					//FileOutputStream output = new FileOutputStream(outputFile);
					//OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
					//for(Integer cpgDist : matrixObj.cpgDistFreq){
					//	writer.write(cpgDist + "\n");
					//}
					//writer.close();
					//output.close();
					//if(noMethyPrior){
					//	bayesianFactor=0;
					//}
					//training
					MatrixObj matrixObj = processMatrixFile(inputFile);
					//System.err.println(miniDataPoints + "\t" + matrixObj.matrix.size());
					if(aucMode){
						
						aucMode(matrixObj, modelFile, outputFile);
					}else{
						
						if(!decodeModeOnly){
							int miniDataPointsPre = miniDataPoints;
							
							if(miniDataPoints < 2){
								miniDataPoints = 2;
								MatrixObj matrixObj2 = processMatrixFile(inputFile);
								//System.err.println(miniDataPoints + "\t" + matrixObj2.matrix.size() + "\t" + matrixObj.matrix.size());
								trainHmm(matrixObj2, modelFile);
							}else{
								trainHmm(matrixObj, modelFile);
							}
							miniDataPoints = miniDataPointsPre;
						}
						//System.err.println(miniDataPoints + "\t" + matrixObj.matrix.size());
						decodeHmm(matrixObj, modelFile, outputFile, inputFile, false);
					}
					
					
					finish();
					

	}
	
	private MatrixObj processMatrixFile(String matrixFile) throws FileNotFoundException, IOException, FileFormatException{
		ArrayList<Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>>> matrix = new ArrayList<Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>>>();
		HashMap<String, TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>>> matrixProcess = new HashMap<String, TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>>>();
		ArrayList<Pair<Integer, Double>> cpgDistFreq = new ArrayList<Pair<Integer, Double>>();
		HashMap<String,IntervalTree<Integer>> overlapLoc = null;
		//HashMap<String,Integer> fragLen = new HashMap<String,Integer>();
		if(region!=null ){
			log.info("Loading overalpped regions... ");
			overlapLoc = new HashMap<String,IntervalTree<Integer>>();
			GZIPInputStream gzipInputStream = null;
			BufferedReader br;
			if(region.endsWith(".gz")){
				gzipInputStream = new GZIPInputStream(new FileInputStream(region));
				br = new BufferedReader(new InputStreamReader(gzipInputStream));
				
			}else{
				br = new BufferedReader(new FileReader(region));
			}
			String line;
			while( (line = br.readLine()) != null){
				if(line.startsWith("#"))
					continue;
				String[] splitLines = line.split("\t");
				if(splitLines.length<3){
					continue;
				}
				String chr = splitLines[0];
				int start = Integer.parseInt(splitLines[1]);
				int end = Integer.parseInt(splitLines[2]);
				IntervalTree<Integer> tree;
				if(overlapLoc.containsKey(chr)){
					tree = overlapLoc.get(chr);
				}else{
					tree = new IntervalTree<Integer>();
				}
				tree.put(start, end, 1);
				overlapLoc.put(chr, tree);
					
			}
			if(region.endsWith(".gz")){
				gzipInputStream.close();
			}
			br.close();

		}
		
		HashMap<String,IntervalTree<Integer>> excludeLoc = null;
		if(exclude != null){
			log.info("Loading excluded regions... ");
			excludeLoc = new HashMap<String,IntervalTree<Integer>>();
			GZIPInputStream gzipInputStream = null;
			BufferedReader br;
			if(exclude.endsWith(".gz")){
				gzipInputStream = new GZIPInputStream(new FileInputStream(exclude));
				br = new BufferedReader(new InputStreamReader(gzipInputStream));
				
			}else{
				br = new BufferedReader(new FileReader(exclude));
			}
			String line;
			while( (line = br.readLine()) != null){
				if(line.startsWith("#"))
					continue;
				String[] splitLines = line.split("\t");
				if(splitLines.length<3){
					continue;
				}
				String chr = splitLines[0];
				int start = Integer.parseInt(splitLines[1]);
				int end = Integer.parseInt(splitLines[2]);
				IntervalTree<Integer> tree;
				if(excludeLoc.containsKey(chr)){
					tree = excludeLoc.get(chr);
				}else{
					tree = new IntervalTree<Integer>();
				}
				tree.put(start, end, 1);
				excludeLoc.put(chr, tree);
					
			}
			if(exclude.endsWith(".gz")){
				gzipInputStream.close();
			}
			br.close();
		}
		
		
		GZIPInputStream gzipInputStream = null;
		BufferedReader br;
		if(matrixFile.endsWith(".gz")){
			gzipInputStream = new GZIPInputStream(new FileInputStream(matrixFile));
			br = new BufferedReader(new InputStreamReader(gzipInputStream));
			
		}else{
			br = new BufferedReader(new FileReader(matrixFile));
		}
			
			String line;
			SummaryStatistics[] stats = new SummaryStatistics[3];
			//HashMap<String, SummaryStatistics[]> stats = new HashMap<String, SummaryStatistics[]>();
			for(int i = 0; i < 3; i++){
				stats[i] = new SummaryStatistics();
			}

			while( (line = br.readLine()) != null){
				if(line.startsWith("#"))
					continue;
				String[] splitLines = line.split("\t");
				if(splitLines.length< (features + 4) || splitLines[1].equalsIgnoreCase("start") || Integer.parseInt(splitLines[4]) >= maxFragLen || Integer.parseInt(splitLines[4])  <= minFragLen || Double.parseDouble(splitLines[8]) <= 5){
					continue;
				}
				String chr = splitLines[0];
				int start = Integer.parseInt(splitLines[1]);
				int end = Integer.parseInt(splitLines[2]);
				
				if(region!=null ){
					if(overlapLoc.containsKey(chr)){
						if(overlapLoc.get(chr).minOverlapper(start, end)==null){
							continue;
						}
						
					}else{
						continue;
					}
				}
				
				if(exclude !=null ){
					if(excludeLoc.containsKey(chr)){
						if(excludeLoc.get(chr).minOverlapper(start, end)!=null){
							//System.err.println(chr + "\t" + start + "\t" + end + "\t" + excludeLoc.size() + "\t" + excludeLoc.get(chr).minOverlapper(start, end));
							continue;
						}
					}
				}
				
				Integer offset = Integer.parseInt(splitLines[9]);
				if(offset < 0){
					continue;
				}
				Double methyPrior = Double.parseDouble(splitLines[11]);
				if(Double.isNaN(methyPrior)){
					continue;
				}
				double DistToCenter = Double.parseDouble(splitLines[4])/2-Double.parseDouble(splitLines[10])+0.5;
				stats[0].addValue(Double.parseDouble(splitLines[4]));
				stats[1].addValue(Double.parseDouble(splitLines[7]));
				stats[2].addValue(DistToCenter);
				//if(lowCoverage){
				//	stats[1].addValue(Double.parseDouble(splitLines[10]));
					
				//}else{
				//	stats[1].addValue(Double.parseDouble(splitLines[7]));
				//	stats[2].addValue(Double.parseDouble(splitLines[10]));
				//}
				
				//stats[3].addValue(Double.parseDouble(splitLines[12])*Double.parseDouble(splitLines[4]));
				//if(stats.containsKey(chr)){
				//	SummaryStatistics[] statsTmp = stats.get(chr);
					
					//statsTmp[0].addValue(Double.parseDouble(splitLines[4]));
					//statsTmp[1].addValue(Double.parseDouble(splitLines[7]));
					//statsTmp[2].addValue(Double.parseDouble(splitLines[10]));
				//	stats.put(chr,statsTmp);
			//	}else{
				//	SummaryStatistics[] statsTmp = new SummaryStatistics[features];
			//		for(int i = 0; i < features; i++){
			//			statsTmp[i] = new SummaryStatistics();
						
			//		}
			//		statsTmp[0].addValue(Double.parseDouble(splitLines[4]));
			//		statsTmp[1].addValue(Double.parseDouble(splitLines[7]));
			//		statsTmp[2].addValue(Double.parseDouble(splitLines[10]));
			//		stats.put(chr,statsTmp);
			//	}
				
				
			}
			if(matrixFile.endsWith(".gz")){
				gzipInputStream.close();
			}
			br.close();
		
		
		if(matrixFile.endsWith(".gz")){
			gzipInputStream = new GZIPInputStream(new FileInputStream(matrixFile));
			br = new BufferedReader(new InputStreamReader(gzipInputStream));
			
		}else{
			br = new BufferedReader(new FileReader(matrixFile));
		}
		for(int i = 0; i < 3; i++){
			log.info("Feature " + i + ": " + stats[i]);
		}
			

			while( (line = br.readLine()) != null){
				if(line.startsWith("#"))
					continue;
				String[] splitLines = line.split("\t");
				if(splitLines.length< (features + 4) || splitLines[1].equalsIgnoreCase("start") || Integer.parseInt(splitLines[4]) >= maxFragLen || Integer.parseInt(splitLines[4])  <= minFragLen || Double.parseDouble(splitLines[8]) <= 5){
					continue;
				}
				String chr = splitLines[0];
				int start = Integer.parseInt(splitLines[1]);
				int end = Integer.parseInt(splitLines[2]);
				String loc = chr + ":" + start + ":" + end;
				String readName = splitLines[3];
				
				if(region!=null ){
					if(overlapLoc.containsKey(chr)){
						if(overlapLoc.get(chr).minOverlapper(start, end)==null){
							continue;
						}
						
					}else{
						continue;
					}
				}
				if(exclude !=null ){
					if(excludeLoc.containsKey(chr)){
						if(excludeLoc.get(chr).minOverlapper(start, end)!=null){
							continue;
						}
					}
				}
				
				Integer offset = Integer.parseInt(splitLines[9]);
				if(offset < 0){
					continue;
				}
				Double methyPrior = Double.parseDouble(splitLines[11]);
				//if(noMethyPrior){
				//	methyPrior = 50.0;
				//}else{
					if(Double.compare(methyPrior, 100.0)==0){
						methyPrior -= 0.01;
					}else if(Double.compare(methyPrior, 0.0)==0){
						methyPrior += 0.01;
					}
				//}

				methyPrior /= 100;
				if(Double.isNaN(methyPrior)){
					continue;
				}
				double strand =splitLines[5].equalsIgnoreCase("-") ? 1.0 : 0.0;
				//filter out those regions with too much variation to the mean?
				//double baseQ = 1-Math.pow(10, (-Double.parseDouble(splitLines[8])/10.0));
				double DistToCenter = Double.parseDouble(splitLines[4])/2-Double.parseDouble(splitLines[10])+0.5;
				if(covOutlier > 0 && ((Double.parseDouble(splitLines[7])-stats[1].getMean())/stats[1].getStandardDeviation() > covOutlier || 
						(Double.parseDouble(splitLines[4])-stats[0].getMean())/stats[0].getStandardDeviation() > covOutlier ||
						(DistToCenter-stats[2].getMean())/stats[2].getStandardDeviation() > covOutlier)){
					continue;
				}
				double[] value;

				if(lowCoverage){
					value = new double[]{
							(Double.parseDouble(splitLines[4])-stats[0].getMean())/stats[0].getStandardDeviation(),
							(DistToCenter-stats[2].getMean())/stats[2].getStandardDeviation(),
					};
					
				}else{
					value = new double[]{
							(Double.parseDouble(splitLines[4])-stats[0].getMean())/stats[0].getStandardDeviation(),
							//Double.parseDouble(splitLines[4]),
							(Double.parseDouble(splitLines[7])-stats[1].getMean())/stats[1].getStandardDeviation(),
							//Double.parseDouble(splitLines[7]),
							(DistToCenter-stats[2].getMean())/stats[2].getStandardDeviation(),
							//Double.parseDouble(splitLines[10]),
							//(Double.parseDouble(splitLines[12])*Double.parseDouble(splitLines[4])-stats[3].getMean())/stats[3].getStandardDeviation(),
							//Double.parseDouble(splitLines[12])*Double.parseDouble(splitLines[4]),
							
					};
				}
				
				ObservationVector vector = new ObservationVector(value);
			////	if(!fragLen.containsKey(readName)){
			//		fragLen.put(readName, Integer.parseInt(splitLines[4]));
			//	}
				
				points++;
				if(matrixProcess.containsKey(readName)){
					TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>> readStat = matrixProcess.get(readName);
					if(readStat.containsKey(offset)){
						//throw new FileFormatException("Contain duplicate reads with the same read name and same offsets:" + readName + "\tOffset:" + offset);
						//TODO: find out why??
					}else{
						readStat.put(offset, Triple.of(splitLines[6], vector, new Pair<String, Double>(loc, methyPrior)));
					}
					matrixProcess.put(readName, readStat);
					
				}else{
					TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>> readStat =  new TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>>();
					readStat.put(offset, Triple.of(splitLines[6], vector, new Pair<String, Double>(loc, methyPrior)));
					matrixProcess.put(readName, readStat);
				}
				
	
				
			}
			if(matrixFile.endsWith(".gz")){
				gzipInputStream.close();
			}
			br.close();
			log.info("Number of point in total is loaded : " + points);
			//for(int i = 0; i < features; i++){
			//	log.info(stats.get("22")[i]);
				
			//}

			ArrayList<ObservationVector> matrixU = new ArrayList<ObservationVector>();
			ArrayList<ObservationVector> matrixM = new ArrayList<ObservationVector>();
			ArrayList<ArrayList<Integer>> matrixObserved = new ArrayList<ArrayList<Integer>>();
			TreeMap<Integer, Long[]> pi = new TreeMap<Integer, Long[]>();
			TreeMap<Integer, Long[]> aij = new TreeMap<Integer, Long[]>();
			
			for(String key : matrixProcess.keySet()){
				TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>> readStat = matrixProcess.get(key);
				
				
				ArrayList<ObservationVector> matrixRow = new ArrayList<ObservationVector>();
				HashMap<Integer, Pair<Integer, Double>> cpgDistRow = new HashMap<Integer, Pair<Integer, Double>>();
				ArrayList<String> locRow = new ArrayList<String>();
				Integer[] offsets = readStat.keySet().toArray(new Integer[readStat.keySet().size()]);
				boolean omitRead = false;
				for(int i = 0; i < offsets.length; i++){
					
					if(i == 0){
						if(offsets[i] < 0){
							omitRead = true;
							continue;
						}
						if(offsets[i] > maxCpgDist){
							omitRead = true;
							break;
						}
						cpgDistRow.put(i, new Pair<Integer, Double>((int)(offsets[i]/bin), readStat.get(offsets[i]).getRight().getSecond()));
					}else{
						int cpgDist = (int)((offsets[i] - offsets[i-1])/bin);
						if(cpgDist<0){
							omitRead = true;
							break;
						}
						if(cpgDist*bin > maxCpgDist){
							omitRead = true;
							break;
						}
						cpgDistRow.put(i, new Pair<Integer, Double>(cpgDist, readStat.get(offsets[i]).getRight().getSecond()));
					}
					matrixRow.add(readStat.get(offsets[i]).getMiddle());
					locRow.add(readStat.get(offsets[i]).getRight().getFirst());
				}
				if(omitRead){
					continue;
				}
				if(matrixRow.size() >= miniDataPoints && matrixRow.size() <= maxCpgs){
					matrix.add(Triple.of(cpgDistRow, matrixRow, locRow));
					ArrayList<Integer> matrixRowObserved = new ArrayList<Integer>();
					for(int i = 0; i < offsets.length; i++){
						if(i==0 && offsets[0] < 0){
							continue;
						}
						if(i>0 && (int)((offsets[i]-offsets[i-1]))<0 ){
							continue;
						}
						if(readStat.get(offsets[i]).getLeft().equalsIgnoreCase("u")){
							matrixU.add(readStat.get(offsets[i]).getMiddle());
							if(i==0){
								if(offsets[0] > maxCpgDist){
									continue;
								}
								
								if(pi.containsKey(offsets[0])){
									Long[] piTmp = pi.get(offsets[0]);
									piTmp[0]++;
									pi.put((int)offsets[0]/bin, piTmp);
								}else{
									Long[] piTmp = new Long[]{0L,0L};
									piTmp[0]++;
									pi.put((int)offsets[0]/bin, piTmp);
								}
									
								cpgDistFreq.add(new Pair<Integer, Double>(offsets[0]/bin,readStat.get(offsets[0]).getRight().getSecond()));
								//if(cpgDistFreq.containsKey(offsets[0])){
								//	cpgDistFreq.put(offsets[0], cpgDistFreq.get(offsets[0])+1);
								//}else{
								//	cpgDistFreq.put(offsets[0], 1L);
								//}
							}else{
								int cpgDist = (int)((offsets[i]-offsets[i-1])/bin);
								if(cpgDist*bin > maxCpgDist){
									continue;
								}
								if(readStat.get(offsets[i-1]).getLeft().equalsIgnoreCase("u")){
									
									if(aij.containsKey(cpgDist)){
										Long[] aijTmp = aij.get(cpgDist);
										aijTmp[0]++;
										aij.put(cpgDist, aijTmp);
									}else{
										Long[] aijTmp = new Long[]{0L,0L, 0L, 0L};
										aijTmp[0]++;
										aij.put(cpgDist, aijTmp);
									}
									
								}else if(readStat.get(offsets[i-1]).getLeft().equalsIgnoreCase("m")){
									if(aij.containsKey(cpgDist)){
										Long[] aijTmp = aij.get(cpgDist);
										aijTmp[2]++;
										aij.put(cpgDist, aijTmp);
									}else{
										Long[] aijTmp = new Long[]{0L,0L, 0L, 0L};
										aijTmp[2]++;
										aij.put(cpgDist, aijTmp);
									}
								}
								cpgDistFreq.add(new Pair<Integer, Double>(cpgDist,readStat.get(offsets[i]).getRight().getSecond()));
								//if(cpgDistFreq.containsKey(cpgDist)){
								//	cpgDistFreq.put(cpgDist, cpgDistFreq.get(cpgDist)+1);
								//}else{
								//	cpgDistFreq.put(cpgDist, 1L);
								//}
							}

							matrixRowObserved.add(0);
						}else if(readStat.get(offsets[i]).getLeft().equalsIgnoreCase("m")){
							matrixM.add(readStat.get(offsets[i]).getMiddle());
							if(i==0){
								if(offsets[0] > maxCpgDist){
									continue;
								}
								if(pi.containsKey(offsets[0])){
									Long[] piTmp = pi.get(offsets[0]);
									piTmp[1]++;
									pi.put((int)offsets[0]/bin, piTmp);
								}else{
									Long[] piTmp = new Long[]{0L,0L};
									piTmp[1]++;
									pi.put((int)offsets[0]/bin, piTmp);
								}
								cpgDistFreq.add(new Pair<Integer, Double>(offsets[0]/bin,readStat.get(offsets[0]).getRight().getSecond()));
								//if(cpgDistFreq.containsKey(offsets[0])){
								//	cpgDistFreq.put(offsets[0], cpgDistFreq.get(offsets[0])+1);
								//}else{
								//	cpgDistFreq.put(offsets[0], 1L);
								//}
							}else{
								int cpgDist = (int)((offsets[i]-offsets[i-1])/bin);
								if(cpgDist*bin > maxCpgDist){
									continue;
								}
								if(readStat.get(offsets[i-1]).getLeft().equalsIgnoreCase("u")){
									if(aij.containsKey(cpgDist)){
										Long[] aijTmp = aij.get(cpgDist);
										aijTmp[1]++;
										aij.put(cpgDist, aijTmp);
									}else{
										Long[] aijTmp = new Long[]{0L,0L, 0L, 0L};
										aijTmp[1]++;
										aij.put(cpgDist, aijTmp);
									}
									
								}else if(readStat.get(offsets[i-1]).getLeft().equalsIgnoreCase("m")){
									if(aij.containsKey(cpgDist)){
										Long[] aijTmp = aij.get(cpgDist);
										aijTmp[3]++;
										aij.put(cpgDist, aijTmp);
									}else{
										Long[] aijTmp = new Long[]{0L,0L, 0L, 0L};
										aijTmp[3]++;
										aij.put(cpgDist, aijTmp);
									}
								}
								cpgDistFreq.add(new Pair<Integer, Double>(cpgDist,readStat.get(offsets[i]).getRight().getSecond()));
								//if(cpgDistFreq.containsKey(cpgDist)){
								//	cpgDistFreq.put(cpgDist, cpgDistFreq.get(cpgDist)+1);
								//}else{
								//	cpgDistFreq.put(cpgDist, 1L);
								//}
							}
							
							matrixRowObserved.add(1);
						}
					}
					matrixObserved.add(matrixRowObserved);
					//double cpgDense = (double)matrixRow.size()/(double)fragLen.get(key);
					
					
					
					//if(cpgDense<this.minCpgNum){
					//	this.minCpgNum = cpgDense;
					//}
					//
				}
				if(cpgNumClip < 0){
					double cpgDense = (double)matrixRow.size();
					if(cpgDense>this.maxCpgNum){
						this.maxCpgNum = cpgDense;
					}
					
				}
			}
			TreeMap<Integer, double[]> piScale = new TreeMap<Integer, double[]>();
			for(Integer cpgDist : pi.keySet()){
				Long[] piTmp = pi.get(cpgDist);
				double[] piScaleTmp = new double[]{(double)piTmp[0]/(double)(piTmp[0] + piTmp[1]), 
						(double)piTmp[1]/(double)(piTmp[0] + piTmp[1])};
				piScale.put(cpgDist, piScaleTmp);
				
			}
			
			
			TreeMap<Integer, double[][]> aijScale = new TreeMap<Integer, double[][]>();
			for(Integer cpgDist : aij.keySet()){
				Long[] aijTmp = aij.get(cpgDist);
				double[][] aiScaleTmp = new double[][]{{(double)aijTmp[0]/(double)(aijTmp[0] + aijTmp[1]),
					(double)aijTmp[1]/(double)(aijTmp[0] + aijTmp[1])},
					{(double)aijTmp[2]/(double)(aijTmp[2] + aijTmp[3]),
						(double)aijTmp[3]/(double)(aijTmp[2] + aijTmp[3])}};
				aijScale.put(cpgDist, aiScaleTmp);
				
			}
		//	System.out.println(this.maxCpgNum + "\t" + this.minCpgNum);
			//output cpg distance distribution to see what is the best model...
			//System.out.println(matrix.size());
		
			log.info("maximum number of cpg in fragment is: " + this.maxCpgNum);
		log.info("The number of fragments used for training: " + matrix.size());
		return new MatrixObj(matrix, matrixU, matrixM, piScale, aijScale, matrixObserved, cpgDistFreq);
	}
	
	//initiate HMM && training HMM
	private void trainHmm(MatrixObj matrixObj, String modelFile) throws IOException, CloneNotSupportedException{
		System.out.println("HMM training new ....");
		List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> matrix = new ArrayList<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>>();
		
		for(Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>> row : matrixObj.matrix){
			matrix.add(new Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>(row.getLeft(), row.getMiddle()));
			
		}
		BayesianNhmmV5<ObservationVector> hmm = wgbs ? buildInitNhmm(matrixObj, true) : (gmm ? buildInitNhmmByGMM(matrixObj, matrix) :buildInitNhmmRandom(matrixObj, true));
		//System.err.println(matrixObj.cpgDistFreq.size() + matrixObj.cpgDistFreq.toArray(new Double[matrixObj.cpgDistFreq.size()])[0]);
		//Regression reg = new Regression(ArrayUtils.toPrimitive(matrixObj.cpgDistFreq.toArray(new Double[matrixObj.cpgDistFreq.size()])),1.0);
		//reg.poisson();
		//double lamda = reg.getBestEstimates()[0];
		//double lamda = matrixObj.lambda; //the mean of 
		//double lamda = 35.0;
		//System.err.println(lamda + "\t" + reg.getBestEstimates().length);
		hmm.setBayesianFactor(bayesianFactor);
		hmm.setMethyState(this.methylatedState);
		hmm.setMaxCpgNum(cpgNumClip < 0 ? maxCpgNum : cpgNumClip);
		hmm.setMinCpgNum(1);
		
		BaumWelchBayesianNhmmV5ScaledLearner bwl = new BaumWelchBayesianNhmmV5ScaledLearner();

		BayesianNhmmV5<ObservationVector> prevHmm = null;
		try {
			prevHmm = hmm.clone();
		} catch (CloneNotSupportedException e1) {
			e1.printStackTrace();
		}
		
		// This object measures the distance between two HMMs
		KullbackLeiblerDistanceBayesianNhmmV5Calculator klc = 
			new KullbackLeiblerDistanceBayesianNhmmV5Calculator(matrix);
		
		double distance = Double.MAX_VALUE;
		double distancePre = 0.01;
		int a = 0;
		// Incrementally improve the solution
		//while(Math.abs(distance) >= tol && a < iteration){
			while(Math.abs((Math.abs(distance)-Math.abs(distancePre))/distancePre) >= decayRate && Math.abs(distance) >= tol && a < iteration){
			a++;
			
			try {
				prevHmm = hmm.clone();
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println("HMM pre:\n" + prevHmm);
			//double[][] covU = ((OpdfMultiGaussian)prevHmm.getOpdf(0)).covariance();
			//String s = "";
			//for(int x = 0; x < covU.length; x++){
			//	for(int y = 0; y < covU[x].length; y++){
			//		s = s + "\t" + covU[x][y];
			//	}
			//	s += "\n";
			//}
			//System.out.println("State 0:\n" + s);
			//double[][] covM = ((OpdfMultiGaussian)prevHmm.getOpdf(1)).covariance();
			//s = "";
			//for(int x = 0; x < covM.length; x++){
			//	for(int y = 0; y < covM[x].length; y++){
			//		s = s + "\t" + covM[x][y];
			//	}
			//	s += "\n";
			//}
			//System.out.println("State 1:\n" + s);
			//hmm = bwl.iterate(hmm, result.value);
			//for (int r = 0; r <= hmm.nbCpgDistState(); r++) {
			//	for (int i = 0; i < hmm.nbStates(); i++) {
					
			//		for (int j = 0; j < hmm.nbStates(); j++) {
						
			//			System.err.println("Arji: " + hmm.getArij(r, i, j) + "\t" + r + "\t" + i + "\t" + j);
			//		}
			//	}
			
			//}
			hmm = bwl.iterate(hmm, matrix);
			distancePre = distance;
			distance = klc.distance(prevHmm, hmm, true);
			System.out.println("Distance at iteration " + a + ": " +
					distance + "\t" + Math.abs((Math.abs(distance)-Math.abs(distancePre))/distancePre));
			if(Double.isNaN(distance)){
				System.out.println("Random initiaton this time does not work. Restart at the new random point...");
				hmm =  wgbs ? buildInitNhmm(matrixObj, true) : (gmm ? buildInitNhmmByGMM(matrixObj, matrix) :buildInitNhmmRandom(matrixObj, true));
				distance = Double.MAX_VALUE;
				distancePre = 0.01;
			}
		}  
			/*
			BayesianNhmmV5<ObservationVector> hmmMethy1 = hmm.clone();
			
			methylatedState = 1-methylatedState;
			distance = Double.MAX_VALUE;
			distancePre = 0.01;
			a = 0;
			
				while(Math.abs((Math.abs(distance)-Math.abs(distancePre))/distancePre) >= decayRate && Math.abs(distance) >= tol && a < iteration){
				a++;
				
				try {
					prevHmm = hmm.clone();
				} catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				System.out.println("HMM pre:\n" + prevHmm);
				
				hmm = bwl.iterate(hmm, matrix);
				distancePre = distance;
				distance = klc.distance(prevHmm, hmm);
				System.out.println("Distance at iteration " + a + ": " +
						distance + "\t" + Math.abs((Math.abs(distance)-Math.abs(distancePre))/distancePre));
				if(Double.isNaN(distance)){
					System.out.println("Random initiaton this time does not work. Restart at the new random point...");
					hmm =  wgbs ? buildInitNhmm(matrixObj, true) : (kmeans ? buildInitNhmmByKmeans(matrixObj, matrix) :buildInitNhmmRandom(matrixObj, true));
					distance = Double.MAX_VALUE;
					distancePre = 0.01;
				}
			} 
				BayesianNhmmV5<ObservationVector> hmmMethy0 = hmm.clone();
			
				KullbackLeiblerDistanceBayesianNhmmV5Calculator klc2 = 
						new KullbackLeiblerDistanceBayesianNhmmV5Calculator(matrix);
				System.out.println("Distance between methy0 and methy 1 is:\n" + klc2.distance(hmmMethy0, hmmMethy1));
			*/
		System.out.println("Resulting NHMM:\n" + hmm);
		this.methylatedState = hmm.getMethyState(lowCoverage);
		//hmm.setMethyState(this.methylatedState);
		hmm.setMaxCpgNum(cpgNumClip < 0 ? maxCpgNum : cpgNumClip);
		hmm.setMinCpgNum(1);
		ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(modelFile));
		oos.writeObject(hmm);
		oos.close();
		//HmmWriter.write(hmmWriter, new OpdfMultiGaussianWriter(), hmm);
		
	}
	
	
	//strategy 1: when met state 3, choose between state 1 and state 2 that who give the highest probability of the sequences.
	//strategy 2: use sliding window to scan the region, and determin who is the closest state to state 3.
	/*
	private void decodeHmmAssignMethyState(MatrixObj matrixObj, String hmmFile, String outputFile) throws Exception{
		System.out.println("\nDecoding ...\n");	
		List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> matrix = new ArrayList<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>>();
		List<ArrayList<String>> matrixLoc = new ArrayList<ArrayList<String>>();
		for(Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>> row : matrixObj.matrix){
			matrix.add(new Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>(row.getLeft(), row.getMiddle()));
			matrixLoc.add(row.getRight());
		}
		
		ArrayList<ArrayList<Integer>> matrixObserved = matrixObj.matrixObserved;

	    ObjectInputStream objectinputstream = new ObjectInputStream(new FileInputStream(hmmFile));
	    BayesianNhmmV5<ObservationVector> hmm = (BayesianNhmmV5<ObservationVector>) objectinputstream.readObject();
		objectinputstream.close();
		long count = 0;
		long countCorrect = 0;
		long countMethy = 0;
		long countMethyCorrect = 0;
		long countUnmethy = 0;
		long countUnmethyCorrect = 0;
		
		
		HashMap<String,Pair<Pair<Integer,Integer>, Pair<Integer,Integer>>> methySummary = new HashMap<String,Pair<Pair<Integer,Integer>, Pair<Integer,Integer>>>(); //predict, observed
		for(int j=0; j < matrix.size(); j++){

			int[] hiddenState = (new ViterbiBayesianNhmmV5Calculator(matrix.get(j), hmm, decodeP)).stateSequence();
			Integer[] observedState = matrixObserved.get(j).toArray(new Integer[matrixObserved.get(j).size()]);
			
			ArrayList<String> locRow = matrixLoc.get(j);
			
			
			if(hiddenState.length != observedState.length){
				throw new IllegalArgumentException("HiddenState Length does not match with observed state length");
			}
			
			
			
			
			
			for(int i = 0; i < hiddenState.length; i++){
				String loc = locRow.get(i);
				
				int methyPredict = 0;
				int unmethyPredict = 0;
				if(hiddenState[i] % 2 == 1 ){
					methyPredict++;
				}else if(hiddenState[i]  % 2 == 0 ){
					unmethyPredict++;
				}
				
				if(methySummary.containsKey(loc)){
					Pair<Pair<Integer,Integer>, Pair<Integer,Integer>> tmp = methySummary.get(loc);
					tmp = new Pair<Pair<Integer,Integer>, Pair<Integer,Integer>>(new Pair<Integer,Integer>(tmp.getFirst().getFirst()+methyPredict,tmp.getFirst().getSecond()+(methyPredict + unmethyPredict)), 
																				new Pair<Integer,Integer>(tmp.getSecond().getFirst()+observedState[i],tmp.getSecond().getSecond()+1));
					methySummary.put(loc, tmp);
				}else{
					methySummary.put(loc, new Pair<Pair<Integer,Integer>, Pair<Integer,Integer>>(new Pair<Integer,Integer>(methyPredict,(methyPredict + unmethyPredict)), new Pair<Integer,Integer>(observedState[i],1)));
				}

				if(observedState[i]==0){
					countUnmethy++;
					if(hiddenState[i] % 2 == 0 ){
						countCorrect++;
						countUnmethyCorrect++;
					}
				}else{
					countMethy++;
					if(hiddenState[i] % 2 == 1){
						countCorrect++;
						countMethyCorrect++;
					}
				}
				count++;
				
			}		
			//ForwardBackwardScaledCalculator fbsc = new ForwardBackwardScaledCalculator(matrix.get(j), hmm);
		}
		
		FileOutputStream output = new FileOutputStream(outputFile);
		OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
		
		writer.write("#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict\ttotal_count_predict\tmethy_perc_obs\tmethy_count_obs\ttotal_count_obs\n");
		
		
		
		double[][] predData = new double[methySummary.size()][2];
		int i = 0;
		for(String loc : methySummary.keySet()){
			int methyPred = methySummary.get(loc).getFirst().getFirst();
			int totalPred = methySummary.get(loc).getFirst().getSecond();
			int methyObs = methySummary.get(loc).getSecond().getFirst();
			int totalObs = methySummary.get(loc).getSecond().getSecond();
			String[] locTmp = loc.split(":");
			String chr = locTmp[0];
			int start = Integer.parseInt(locTmp[1]);
			int end = Integer.parseInt(locTmp[2]);
			double pred = 100*(double)methyPred/(double)totalPred;
			double obs = 100*(double)methyObs/(double)totalObs;
			writer.write(chr + "\t" + start + "\t" + end + "\t" + pred + "\t" + methyPred + "\t" + totalPred + 
					"\t" + obs + "\t" + methyObs + "\t" + totalObs + "\n");
			predData[i][0]=pred;
			predData[i][1]=obs;
			i++;
		}
		writer.close();
		output.close();
		PearsonsCorrelation pearson =  new PearsonsCorrelation(predData);
		System.out.println(pearson.getCorrelationMatrix());
		System.out.println(pearson.getCorrelationPValues());
		System.out.println("counted point in total: " + count + "\tCorrect predicted:" + countCorrect + "\tPerc:" + 100*(double)countCorrect/(double)count + "%");
		System.out.println("counted point in methy: " + countMethy + "\tCorrect predicted:" + countMethyCorrect + "\tPerc:" + 100*(double)countMethyCorrect/(double)countMethy + "%");
		System.out.println("counted point in unmethy: " + countUnmethy + "\tCorrect predicted:" + countUnmethyCorrect + "\tPerc:" + 100*(double)countUnmethyCorrect/(double)countUnmethy + "%");
		
		
		
		
		//for(int f=0;f<=features;f++){
		//	System.err.println( unmethyCorrectSummary.get(f).getMean() + "\t"  + unmethyCorrectSummary.get(f).getStandardDeviation() + "\t" + unmethyCorrectSummary.get(f).getPercentile(50));
		//	System.err.println( unmethyWrongSummary.get(f).getMean() + "\t"  + unmethyWrongSummary.get(f).getStandardDeviation() + "\t" + unmethyWrongSummary.get(f).getPercentile(50));
		//	System.err.println( methyCorrectSummary.get(f).getMean() + "\t"  + methyCorrectSummary.get(f).getStandardDeviation() + "\t" + methyCorrectSummary.get(f).getPercentile(50));
		//	System.err.println( methyWrongSummary.get(f).getMean() + "\t"  + methyWrongSummary.get(f).getStandardDeviation() + "\t" + methyWrongSummary.get(f).getPercentile(50));
		//}
		
	}
	*/
	
	//decoding HMM 
	private double decodeHmm(MatrixObj matrixObj, String hmmFile, String outputFile, String inputFile, boolean reestimate) throws Exception{
		System.out.println("\nDecoding ...\n");
		//System.out.println("\nMethylation state is:" + methylatedState);
		List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> matrix = new ArrayList<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>>();
		List<ArrayList<String>> matrixLoc = new ArrayList<ArrayList<String>>();
		for(Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>> row : matrixObj.matrix){
			matrix.add(new Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>(row.getLeft(), row.getMiddle()));
			matrixLoc.add(row.getRight());
		}
		
		ArrayList<ArrayList<Integer>> matrixObserved = matrixObj.matrixObserved;

	    ObjectInputStream objectinputstream = new ObjectInputStream(new FileInputStream(hmmFile));
	    BayesianNhmmV5<ObservationVector> hmm = (BayesianNhmmV5<ObservationVector>) objectinputstream.readObject();
		objectinputstream.close();
		hmm.setBayesianFactor(bayesianFactor);
		hmm.setMethyState(this.methylatedState);
		hmm.setMaxCpgNum(cpgNumClip < 0 ? maxCpgNum : cpgNumClip);
		hmm.setMinCpgNum(1);
		long count = 0;
		long countCorrect = 0;
		long countMethy = 0;
		long countMethyCorrect = 0;
		long countUnmethy = 0;
		long countUnmethyCorrect = 0;
		
		
		HashMap<String,Pair<Pair<Integer,Integer>, Pair<Integer,Integer>>> methySummary = new HashMap<String,Pair<Pair<Integer,Integer>, Pair<Integer,Integer>>>(); //predict, observed
		
		double likelihood = 0;
		double likelihoodWithMethy = 0;
		for(int j=0; j < matrix.size(); j++){
			
			ViterbiBayesianNhmmV5Calculator vb = new ViterbiBayesianNhmmV5Calculator(matrix.get(j), hmm, methylatedState, decodeP);
			int[] hiddenState = vb.stateSequence();
			//System.err.println(vb.lnProbability(true) + "\t" + vb.lnProbability() + "\t" + hiddenState.length);
			//double prob = vb.lnProbability(true);
			double prob = vb.lnProbability();
			if(Double.isNaN(prob) || Double.isInfinite(prob)){
				
			}else{
				likelihood += prob/hiddenState.length;
			}
			double probWithMethy = vb.lnProbability(true);
			if(Double.isNaN(probWithMethy) || Double.isInfinite(probWithMethy)){
				
			}else{
				likelihoodWithMethy += probWithMethy/hiddenState.length;
			}
			
			Integer[] observedState = matrixObserved.get(j).toArray(new Integer[matrixObserved.get(j).size()]);
			
			ArrayList<String> locRow = matrixLoc.get(j);
			
			
			if(hiddenState.length != observedState.length){
				throw new IllegalArgumentException("HiddenState Length does not match with observed state length");
			}
			for(int i = 0; i < hiddenState.length; i++){
				
				String loc = locRow.get(i);
				
				int methyPredict = 0;
				int unmethyPredict = 0;
				if(hiddenState[i] == methylatedState ){
					methyPredict++;
				}else if(hiddenState[i]  == (1-methylatedState) ){
					unmethyPredict++;
				}
				
				
				if(methySummary.containsKey(loc)){
					Pair<Pair<Integer,Integer>, Pair<Integer,Integer>> tmp = methySummary.get(loc);
					tmp = new Pair<Pair<Integer,Integer>, Pair<Integer,Integer>>(new Pair<Integer,Integer>(tmp.getFirst().getFirst()+methyPredict,tmp.getFirst().getSecond()+(methyPredict + unmethyPredict)), 
																				new Pair<Integer,Integer>(tmp.getSecond().getFirst()+observedState[i],tmp.getSecond().getSecond()+1));
					methySummary.put(loc, tmp);
				}else{
					methySummary.put(loc, new Pair<Pair<Integer,Integer>, Pair<Integer,Integer>>(new Pair<Integer,Integer>(methyPredict,(methyPredict + unmethyPredict)), new Pair<Integer,Integer>(observedState[i],1)));
				}
				
				
				if(observedState[i]==0){
					countUnmethy++;
					if(hiddenState[i] == (1-methylatedState) ){
						countCorrect++;
						countUnmethyCorrect++;
					}
				}else{
					countMethy++;
					if(hiddenState[i] == methylatedState){
						countCorrect++;
						countMethyCorrect++;
					}
				}
				count++;
			}		
			//ForwardBackwardScaledCalculator fbsc = new ForwardBackwardScaledCalculator(matrix.get(j), hmm);
		}
		
		//if(!notAutomateIdentifyMethyState && !reestimate){
			//if((double)countMethyCorrect/(double)countMethy < 0.4){
			//System.out.println("methyState " + methylatedState + "\tLikelihood is:" + likelihood);	
		//	methylatedState = 1- methylatedState;
		//		System.out.println("\nre-estimate the reversed methylation states\n");
				//re-train
		//		int miniDataPointsPre = miniDataPoints;
		//		if(miniDataPoints==1){
		//			miniDataPoints = 2;
		//			matrixObj = processMatrixFile(inputFile);
		//		}
		//		trainHmm(matrixObj, hmmFile);
				
				//re-decoding
		//		if(miniDataPoints != miniDataPointsPre){
		//			miniDataPoints = miniDataPointsPre;
		//			matrixObj = processMatrixFile(inputFile);
		//		}
		//		double likelihoodCompared = decodeHmm(matrixObj, hmmFile, outputFile, inputFile, true);
				//System.out.println("methyState " + methylatedState + "\tLikelihood is:" + likelihoodCompared);	
				//if(likelihoodCompared > likelihood)
				//	return likelihoodCompared;
				
			//	methylatedState = 1- methylatedState;

			//}
		//}
		
		FileOutputStream output = new FileOutputStream(outputFile);
		OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
		
		writer.write("#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict\ttotal_count_predict\tmethy_perc_obs\tmethy_count_obs\ttotal_count_obs\n");
		
		
		double[][] predData = new double[methySummary.size()][2];
		int i = 0;
		for(String loc : methySummary.keySet()){
			int methyPred = methySummary.get(loc).getFirst().getFirst();
			int totalPred = methySummary.get(loc).getFirst().getSecond();
			int methyObs = methySummary.get(loc).getSecond().getFirst();
			int totalObs = methySummary.get(loc).getSecond().getSecond();
			String[] locTmp = loc.split(":");
			String chr = locTmp[0];
			int start = Integer.parseInt(locTmp[1]);
			int end = Integer.parseInt(locTmp[2]);
			double pred = 100*(double)methyPred/(double)totalPred;
			double obs = 100*(double)methyObs/(double)totalObs;
			writer.write(chr + "\t" + start + "\t" + end + "\t" + pred + "\t" + methyPred + "\t" + totalPred + 
					"\t" + obs + "\t" + methyObs + "\t" + totalObs + "\n");
			predData[i][0]=pred;
			predData[i][1]=obs;
			i++;
		}
		writer.close();
		output.close();
		PearsonsCorrelation pearson =  new PearsonsCorrelation(predData);
		System.out.println(pearson.getCorrelationMatrix());
		System.out.println(pearson.getCorrelationPValues());
		System.out.println("counted point in total: " + count + "\tCorrect predicted:" + countCorrect + "\tPerc:" + 100*(double)countCorrect/(double)count + "%");
		System.out.println("counted point in methy: " + countMethy + "\tCorrect predicted:" + countMethyCorrect + "\tPerc:" + 100*(double)countMethyCorrect/(double)countMethy + "%");
		System.out.println("counted point in unmethy: " + countUnmethy + "\tCorrect predicted:" + countUnmethyCorrect + "\tPerc:" + 100*(double)countUnmethyCorrect/(double)countUnmethy + "%");
		System.out.println("methyState " + methylatedState + "\tLikelihood is:" + likelihood);	
		System.out.println("methyState " + methylatedState + "\tLikelihoodWithMethy is:" + likelihoodWithMethy);	
		//return (double)countMethyCorrect/(double)countMethy;
		return likelihood;
		//for(int f=0;f<=features;f++){
		//	System.err.println( unmethyCorrectSummary.get(f).getMean() + "\t"  + unmethyCorrectSummary.get(f).getStandardDeviation() + "\t" + unmethyCorrectSummary.get(f).getPercentile(50));
		//	System.err.println( unmethyWrongSummary.get(f).getMean() + "\t"  + unmethyWrongSummary.get(f).getStandardDeviation() + "\t" + unmethyWrongSummary.get(f).getPercentile(50));
		//	System.err.println( methyCorrectSummary.get(f).getMean() + "\t"  + methyCorrectSummary.get(f).getStandardDeviation() + "\t" + methyCorrectSummary.get(f).getPercentile(50));
		//	System.err.println( methyWrongSummary.get(f).getMean() + "\t"  + methyWrongSummary.get(f).getStandardDeviation() + "\t" + methyWrongSummary.get(f).getPercentile(50));
		//}
		
	}
	
	
	//decoding HMM 
	private void aucMode(MatrixObj matrixObj, String hmmFile, String outputFile) throws Exception{
		System.out.println("\nROC curve ...\n");	
		List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> matrix = new ArrayList<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>>();
		
		for(Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>> row : matrixObj.matrix){
			matrix.add(new Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>(row.getLeft(), row.getMiddle()));
		
		}
		ArrayList<ArrayList<Integer>> matrixObserved = matrixObj.matrixObserved;

	    ObjectInputStream objectinputstream = new ObjectInputStream(new FileInputStream(hmmFile));
	    BayesianNhmmV5<ObservationVector> hmm = (BayesianNhmmV5<ObservationVector>) objectinputstream.readObject();
		objectinputstream.close();
		hmm.setBayesianFactor(bayesianFactor);
		hmm.setMethyState(this.methylatedState);
		hmm.setMaxCpgNum(cpgNumClip < 0 ? maxCpgNum : cpgNumClip);
		hmm.setMinCpgNum(1);
		boolean outFlag = false;
		double[] ps = new double[]{-1,-0.999, -0.99, -0.95, -0.9,-0.8,-0.7,-0.5,-0.1,0.0, 0.1,0.5,0.7,0.8,0.9,0.95,0.99,0.999,1.0};
		//for(double p=-1; p < 1; p+=0.1){
		for(double p : ps){
			long countMethy = 0;
			long countMethyCorrect = 0;
			long countUnmethy = 0;
			long countUnmethyCorrect = 0;

			
			
			
			
			
			for(int j=0; j < matrix.size(); j++){

				int[] hiddenState = (new ViterbiBayesianNhmmV5Calculator(matrix.get(j), hmm, methylatedState, p)).stateSequence();

				Integer[] observedState = matrixObserved.get(j).toArray(new Integer[matrixObserved.get(j).size()]);
				if(hiddenState.length != observedState.length){
					throw new IllegalArgumentException("HiddenState Length does not match with observed state length");
				}
				if(randomPerm){
					HashMap<Integer, Pair<Integer, Double>> cpgDistState = matrix.get(j).getFirst();
					for(int i = 0; i < observedState.length; i++){
						Double methyPrior = cpgDistState.get(i).getSecond();
						double rand = randomEngine.nextDouble();
						if(rand < methyPrior+p){
							hiddenState[i]=1;
						}else{
							hiddenState[i]=0;
						}
					}
				}
				
				
				

				for(int i = 0; i < hiddenState.length; i++){
					
					
					if(hiddenState[i] % 2 == observedState[i]){
						
						if(observedState[i] == 1){
							countMethyCorrect++;
							
						}else{
							countUnmethyCorrect++;
							
						}
					}
					
					if(observedState[i] == 1){
						countMethy++;
					}else{
						countUnmethy++;
					}
				}		
				
			}
			if(!outFlag){
				System.out.println(countMethy + "\t" + countUnmethy);
				outFlag=true;
			}
			double fpr = (double)(countUnmethy-countUnmethyCorrect)/(double)countUnmethy;
			double tpr = (double)(countMethyCorrect)/(double)countMethy;
			System.out.println(fpr + "\t" + tpr);
		}
		
		

		
		//for(int f=0;f<=features;f++){
		//	System.err.println( unmethyCorrectSummary.get(f).getMean() + "\t"  + unmethyCorrectSummary.get(f).getStandardDeviation() + "\t" + unmethyCorrectSummary.get(f).getPercentile(50));
		//	System.err.println( unmethyWrongSummary.get(f).getMean() + "\t"  + unmethyWrongSummary.get(f).getStandardDeviation() + "\t" + unmethyWrongSummary.get(f).getPercentile(50));
		//	System.err.println( methyCorrectSummary.get(f).getMean() + "\t"  + methyCorrectSummary.get(f).getStandardDeviation() + "\t" + methyCorrectSummary.get(f).getPercentile(50));
		//	System.err.println( methyWrongSummary.get(f).getMean() + "\t"  + methyWrongSummary.get(f).getStandardDeviation() + "\t" + methyWrongSummary.get(f).getPercentile(50));
		//}
		
	}

	
	protected BayesianNhmmV5<ObservationVector> buildInitNhmm(MatrixObj matrixObj, boolean addRandomFluct){	
		
		OpdfMultiMixtureGaussian omgU = new OpdfMultiMixtureGaussian(features, mixNumberInFeature);
		omgU.fit(matrixObj.matrixU);
		OpdfMultiMixtureGaussian omgM = new OpdfMultiMixtureGaussian(features, mixNumberInFeature);
		omgM.fit(matrixObj.matrixM);
		
		BayesianNhmmV5<ObservationVector> hmm = new BayesianNhmmV5<ObservationVector>(states, maxCpgDist/bin, new OpdfMultiMixtureGaussianFactory(features, mixNumberInFeature), bayesianFactor);
		for(int i = 0; i < states; i++){
			if(i % 2 == 0){
				hmm.setOpdf(i, omgU); //state 0 is unmethylated state
			}else{
				hmm.setOpdf(i, omgM);
			}
		}
		
	
		for(int z = 0; z <= maxCpgDist/bin; z++){
			if(matrixObj.pi.containsKey(z)){
				hmm.setPri (z, 0 ,matrixObj.pi.get(z)[0]);
				hmm.setPri (z, 1 ,matrixObj.pi.get(z)[1]);
			}else{
				double randomPi = randomEngine.nextDouble() * 0.5;	
				for(int i = 0; i < states; i++){
					if(i==0){
						hmm.setPri (z, i ,randomPi);
					}else{
						hmm.setPri (z, i ,(1-randomPi)/(states-1));
					}
				}
			}
			if(matrixObj.a.containsKey(z)){
				if(Double.isNaN(matrixObj.a.get(z)[0][0]) || Double.isNaN(matrixObj.a.get(z)[0][1]) || Double.isNaN(matrixObj.a.get(z)[1][0]) || Double.isNaN(matrixObj.a.get(z)[1][1])){
					for(int i = 0; i < states; i++){
						double randomAij = randomEngine.nextDouble()*0.5;
					
						for(int j = 0; j < states; j++){
							if(i==j){
								hmm.setArij (z, i, j, (1-randomAij)/(states-1));
							}else{
								hmm.setArij (z, i, j, randomAij);
							}
						}
					}
				}else{
					hmm.setArij (z, 0, 0, matrixObj.a.get(z)[0][0]);
					hmm.setArij (z, 0, 1, matrixObj.a.get(z)[0][1]);
					hmm.setArij (z, 1, 0, matrixObj.a.get(z)[1][0]);
					hmm.setArij (z, 1, 1, matrixObj.a.get(z)[1][1]);
				}
				
			}else{
				for(int i = 0; i < states; i++){
					double randomAij = randomEngine.nextDouble()*0.5;
				
					for(int j = 0; j < states; j++){
						if(i==j){
							hmm.setArij (z, i, j, (1-randomAij)/(states-1));
						}else{
							hmm.setArij (z, i, j, randomAij);
						}
					}
				}
			}
			
			
		}

		return hmm;
		
	}
	
	
	protected BayesianNhmmV5<ObservationVector> buildInitNhmmByGMM(MatrixObj matrixObj, List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> matrix){	
		System.out.println("GMMLearner...");
		//KMeansPlusLearner kl = new KMeansPlusLearner(states, new OpdfMultiMixtureGaussianFactory(features, mixNumberInFeature),matrix,maxCpgDist/bin,features, mixNumberInFeature,bayesianFactor, randomEngine, tolKmeans,decayKmeans, cpgNumClip, 1, lowCoverage);
		GMMLearner kl = new GMMLearner(states, new OpdfMultiMixtureGaussianFactory(features, mixNumberInFeature),matrix,maxCpgDist/bin,features, mixNumberInFeature,bayesianFactor, randomEngine, tolKmeans,decayKmeans, cpgNumClip, 1, lowCoverage);
		BayesianNhmmV5<ObservationVector> hmm = kl.learn();
		methylatedState = kl.returnMethyState();
		System.out.println("Methylated state is : " + methylatedState);
		return hmm;
	}
	
	protected BayesianNhmmV5<ObservationVector> buildInitNhmmRandom(MatrixObj matrixObj, boolean addRandomFluct){	
		
		ArrayList<ObservationVector> matrixAll = new ArrayList<ObservationVector>();
		matrixAll.addAll(matrixObj.matrixU);
		matrixAll.addAll(matrixObj.matrixM);
		Collections.shuffle(matrixAll);
		int size = matrixAll.size();
		OpdfMultiMixtureGaussian omgU = new OpdfMultiMixtureGaussian(features, mixNumberInFeature);
		omgU.fit(matrixAll.subList(0, (int)(size/2.0)));
		OpdfMultiMixtureGaussian omgM = new OpdfMultiMixtureGaussian(features, mixNumberInFeature);
		omgM.fit(matrixAll.subList((int)(size/2.0)+1, size));
		
		BayesianNhmmV5<ObservationVector> hmm = new BayesianNhmmV5<ObservationVector>(states,maxCpgDist/bin, new OpdfMultiMixtureGaussianFactory(features, mixNumberInFeature), bayesianFactor);
		
		for(int i = 0; i < states; i++){
			if(i % 2 == 0){
				hmm.setOpdf(i, omgU); //state 0 is unmethylated state
			}else{
				hmm.setOpdf(i, omgM);
			}
		}
		

	
		
			
			for(int r = 0; r <= maxCpgDist/bin; r++){
				double randomPi = randomEngine.nextDouble() * 0.5;	
				for(int i = 0; i < states; i++){
					if(i==0){
						hmm.setPri (r, i ,randomPi);
					}else{
						hmm.setPri (r, i ,(1-randomPi)/(states-1));
				}
				double randomAij = randomEngine.nextDouble()*0.5;
				for(int j = 0; j < states; j++){
						if(i==j){
							hmm.setArij (r, i, j, (1-randomAij)/(states-1));
						}else{
							hmm.setArij (r, i, j, randomAij);
						}
				}
				}
			}
			
			


		return hmm;
		
	}
	
	/*
	protected BayesianNhmm<ObservationVector> buildInitNhmmRandom(MatrixObj matrixObj, boolean addRandomFluct){	
		
		OpdfMultiGaussian omgU = new OpdfMultiGaussian(features);
		omgU.fit(matrixObj.matrixU);
		OpdfMultiGaussian omgM = new OpdfMultiGaussian(features);
		omgM.fit(matrixObj.matrixM);
		
		BayesianNhmm<ObservationVector> hmm = new BayesianNhmm<ObservationVector>(states, maxFragLen, new OpdfMultiGaussianFactory(features));
		hmm.setOpdf(0, omgU); //state 0 is unmethylated state
		hmm.setOpdf(1, omgM);

		for(int i = 0; i <= maxFragLen; i++){
			double randomPi = randomEngine.nextDouble()-0.5;
			
				hmm.setPri (i, 0 , 0.2);
				hmm.setPri (i, 1 , 0.8);
			
			
			double randomAij1 = randomEngine.nextDouble()-0.5;
			double randomAij2 = randomEngine.nextDouble()-0.5;
			
				hmm.setArij (i, 0, 0, 0.9);
				hmm.setArij (i, 0, 1, 0.1);
				hmm.setArij (i, 1, 0, 0.1);
				hmm.setArij (i, 1, 1, 0.9);
			
			
		}

		return hmm;
		
	}
	*/
	
	private void initiate() throws NumberFormatException, IOException{
		startTime = System.currentTimeMillis();
		if(seed < 0){
			randomEngine = new MersenneTwister();
		}else{
			randomEngine = new MersenneTwister(seed);
		}
		if(lowCoverage){
			features = 2;
		}
		
		if(mixNumberInFeature == null || mixNumberInFeature.isEmpty()){
			mixNumberInFeature = new ArrayList<Integer>();
			for(int i = 0; i < features; i++){
				mixNumberInFeature.add(1);
			}
		}else if(mixNumberInFeature.size() != features){
			throw new IllegalArgumentException("Wrong number of mixNumberInFeature" + mixNumberInFeature);
		}
		
		
	}

	private void finish(){
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;
		
		System.out.println("Counted " + points + " data points in total");
		System.out.println("CcBayesianNhmm's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
	
	public class MatrixObj{
		ArrayList<Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>>> matrix;
		ArrayList<ArrayList<Integer>> matrixObserved;
		ArrayList<ObservationVector> matrixU;
		ArrayList<ObservationVector> matrixM;
		TreeMap<Integer, double[]> pi;
		TreeMap<Integer, double[][]> a;
		ArrayList<Pair<Integer, Double>> cpgDistFreq; //number of instance at each CpG distance
		//HashMap<String, SummaryStatistics[]> stats;
		//double[] mean;
		//double[] variance;
		
		public MatrixObj(ArrayList<Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>>> matrix, ArrayList<ObservationVector> matrixU,
		ArrayList<ObservationVector> matrixM, TreeMap<Integer, double[]> pi, TreeMap<Integer, double[][]> a){
			this.matrix = matrix;
			this.matrixU = matrixU;
			this.matrixM = matrixM;
			this.pi = pi;
			this.a = a;
			
		}
		
		public MatrixObj(ArrayList<Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>>> matrix, ArrayList<ObservationVector> matrixU,
				ArrayList<ObservationVector> matrixM, TreeMap<Integer, double[]> pi, TreeMap<Integer, double[][]> a, ArrayList<ArrayList<Integer>> matrixObserved){
					this.matrix = matrix;
					this.matrixU = matrixU;
					this.matrixM = matrixM;
					this.pi = pi;
					this.a = a;
					this.matrixObserved = matrixObserved;
		}
		
		public MatrixObj(ArrayList<Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>>> matrix, ArrayList<ObservationVector> matrixU,
				ArrayList<ObservationVector> matrixM, TreeMap<Integer, double[]> pi, TreeMap<Integer, double[][]> a, 
				ArrayList<ArrayList<Integer>> matrixObserved, ArrayList<Pair<Integer, Double>> cpgDistFreq){
					this.matrix = matrix;
					this.matrixU = matrixU;
					this.matrixM = matrixM;
					this.pi = pi;
					this.a = a;
					this.matrixObserved = matrixObserved;
					this.cpgDistFreq = cpgDistFreq;
				//	this.stats = stats;
				//	mean = new double[3];
				//	variance = new double[3];
				//	for(int i = 0; i < 3; i++){
				//		mean[i] = stats[i].getMean();
				//		variance[i] = stats[i].getVariance();
				//	}
		}
		
		public List<Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>>> getSampled(int n){
			ArrayList<Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>>> matrixTmp = (ArrayList<Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>>>) matrix.clone();
			Collections.shuffle(matrixTmp);
			if(matrixTmp.size()<n){
				return matrixTmp;
			}else{
				return matrixTmp.subList(0, n);
			}
		}
	}

}
