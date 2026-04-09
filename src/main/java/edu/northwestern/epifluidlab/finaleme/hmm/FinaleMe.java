/**
 * CcNhmm.java
 * Apr 27, 2016
 * 10:25:30 AM
 * yaping    lyping1986@gmail.com
 */
package edu.northwestern.epifluidlab.finaleme.hmm;


import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixFormat;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.ObjectStreamClass;
import java.io.OutputStreamWriter;
import java.io.SequenceInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.Callable;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.nio.charset.StandardCharsets;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.jetbrains.bio.big.BedGraphParser;
import org.jetbrains.bio.big.BigWigFile;
import org.jetbrains.bio.big.WigSection;
import edu.northwestern.epifluidlab.finaleme.utils.ObservationVector;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;


import be.ac.ulg.montefiore.run.jahmm.io.FileFormatException;


/**
 * use cpg distance dependant transition probability matrix for the non homogenous HMM
 */
public class FinaleMe {

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

	@Option(name="-t",usage="number of threads for parallel training/decode. Use >0 to set explicitly; default uses all available cores.")
	public int threads = -1;

	@Option(name="-patOutput", usage="output UXM-compatible .pat.gz and .beta files for deconvolution. Requires -cpgIndexFile. default: false")
	public boolean patOutput = false;

	@Option(name="-cpgIndexFile", usage="CG_motif bedgraph/bed file listing all CpG positions genome-wide, for building CpG index. Required with -patOutput.")
	public String cpgIndexFile = null;

	@Option(name="-bwOutput", usage="output decode summary as bigWig files (.methy.bw/.cov.bw/.methy_count.bw). Requires -chromSizeFile. default: false")
	public boolean bwOutput = false;

	@Option(name="-chromSizeFile", usage="chromosome sizes file required for bigWig output when -bwOutput is enabled.")
	public String chromSizeFile = null;

	@Option(name="-bedGraphToBigWig", usage="path to UCSC bedGraphToBigWig executable. If unavailable, FinaleMe falls back to built-in Java bigWig writer when dependency is present. default: bedGraphToBigWig")
	public String bedGraphToBigWig = "bedGraphToBigWig";

	@Option(name="-bwStripChrPrefix", usage="strip leading 'chr' from chromosome names before bigWig conversion. default: false")
	public boolean bwStripChrPrefix = false;

	@Option(name="-bwConvertChrMToMT", usage="convert chrM/M chromosome name to MT before bigWig conversion. default: false")
	public boolean bwConvertChrMToMT = false;

	
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "FinaleMe [opts] model input_matrix.txt[.gz] prediction.txt.gz";
	
	private static final Logger log = LoggerFactory.getLogger(FinaleMe.class);

	private static long startTime = -1;
	private static long points = 0;
	private MersenneTwister randomEngine;
	private double maxCpgNum = Double.NEGATIVE_INFINITY;
	private static final String LEGACY_PACKAGE_PREFIX = "org.cchmc.epifluidlab.finaleme.";
	private static final String CURRENT_PACKAGE_PREFIX = "edu.northwestern.epifluidlab.finaleme.";

	private static String remapLegacyClassName(String className) {
		if (className == null || !className.contains(LEGACY_PACKAGE_PREFIX)) {
			return className;
		}
		return className.replace(LEGACY_PACKAGE_PREFIX, CURRENT_PACKAGE_PREFIX);
	}

	private static final class LegacyPackageObjectInputStream extends ObjectInputStream {
		LegacyPackageObjectInputStream(InputStream in) throws IOException {
			super(in);
		}

		@Override
		protected Class<?> resolveClass(ObjectStreamClass desc) throws IOException, ClassNotFoundException {
			String originalName = desc.getName();
			String remappedName = remapLegacyClassName(originalName);
			if (!originalName.equals(remappedName)) {
				ClassLoader loader = Thread.currentThread().getContextClassLoader();
				if (loader == null) {
					loader = FinaleMe.class.getClassLoader();
				}
				try {
					return Class.forName(remappedName, false, loader);
				} catch (ClassNotFoundException e) {
					// Fall back to default resolution if remap target cannot be loaded.
				}
			}
			return super.resolveClass(desc);
		}
	}

	@SuppressWarnings("unchecked")
	private BayesianNhmmV5<ObservationVector> loadHmmModel(String hmmFile) throws IOException, ClassNotFoundException {
		try (LegacyPackageObjectInputStream objectinputstream =
				new LegacyPackageObjectInputStream(new FileInputStream(hmmFile))) {
			return (BayesianNhmmV5<ObservationVector>) objectinputstream.readObject();
		}
	}

	private int resolveThreadCount()
	{
		int resolved = threads > 0 ? threads : Runtime.getRuntime().availableProcessors();
		return Math.max(1, resolved);
	}



	private HashMap<String, IntervalTree<Integer>> loadIntervalFile(String path) throws IOException {
		if (path == null) return null;
		log.info("Loading interval regions from " + path + " ...");
		HashMap<String, IntervalTree<Integer>> intervals = new HashMap<String, IntervalTree<Integer>>();
		GZIPInputStream gzipInputStream = null;
		BufferedReader br;
		if (path.endsWith(".gz")) {
			gzipInputStream = new GZIPInputStream(new FileInputStream(path));
			br = new BufferedReader(new InputStreamReader(gzipInputStream));
		} else {
			br = new BufferedReader(new FileReader(path));
		}
		String line;
		while ((line = br.readLine()) != null) {
			if (line.startsWith("#")) continue;
			String[] splitLines = line.split("\t");
			if (splitLines.length < 3) continue;
			String chr = splitLines[0];
			int start = Integer.parseInt(splitLines[1]);
			int end = Integer.parseInt(splitLines[2]);
			IntervalTree<Integer> tree;
			if (intervals.containsKey(chr)) {
				tree = intervals.get(chr);
			} else {
				tree = new IntervalTree<Integer>();
			}
			tree.put(start, end, 1);
			intervals.put(chr, tree);
		}
		if (path.endsWith(".gz")) {
			gzipInputStream.close();
		}
		br.close();
		return intervals;
	}

	private ParsedRow parseLine(String line,
								HashMap<String, IntervalTree<Integer>> overlapLoc,
								HashMap<String, IntervalTree<Integer>> excludeLoc) {
		if (line.startsWith("#")) return null;
		String[] splitLines = line.split("\t");
		if (splitLines.length < (features + 4) || splitLines[1].equalsIgnoreCase("start")
				|| Integer.parseInt(splitLines[4]) >= maxFragLen
				|| Integer.parseInt(splitLines[4]) <= minFragLen
				|| Double.parseDouble(splitLines[8]) <= 5) {
			return null;
		}
		String chr = splitLines[0];
		int start = Integer.parseInt(splitLines[1]);
		int end = Integer.parseInt(splitLines[2]);

		if (overlapLoc != null) {
			if (overlapLoc.containsKey(chr)) {
				if (overlapLoc.get(chr).minOverlapper(start, end) == null) {
					return null;
				}
			} else {
				return null;
			}
		}

		if (excludeLoc != null) {
			if (excludeLoc.containsKey(chr)) {
				if (excludeLoc.get(chr).minOverlapper(start, end) != null) {
					return null;
				}
			}
		}

		int offset = Integer.parseInt(splitLines[9]);
		if (offset < 0) return null;

		double methyPrior = Double.parseDouble(splitLines[11]);
		if (Double.isNaN(methyPrior)) return null;

		double fragLen = Double.parseDouble(splitLines[4]);
		double coverage = Double.parseDouble(splitLines[7]);
		double distToCenter = fragLen / 2 - Double.parseDouble(splitLines[10]) + 0.5;

		if (Double.compare(methyPrior, 100.0) == 0) {
			methyPrior -= 0.01;
		} else if (Double.compare(methyPrior, 0.0) == 0) {
			methyPrior += 0.01;
		}
		methyPrior /= 100;
		if (Double.isNaN(methyPrior)) return null;

		String loc = chr + ":" + start + ":" + end;
		String readName = splitLines[3];
		return new ParsedRow(readName, splitLines[6], loc, offset, methyPrior, fragLen, coverage, distToCenter);
	}

	private AssembledFragment assembleFragment(
			TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>> readStat) {
		ArrayList<ObservationVector> matrixRow = new ArrayList<ObservationVector>();
		HashMap<Integer, Pair<Integer, Double>> cpgDistRow = new HashMap<Integer, Pair<Integer, Double>>();
		ArrayList<String> locRow = new ArrayList<String>();
		Integer[] offsets = readStat.keySet().toArray(new Integer[readStat.keySet().size()]);
		boolean omitRead = false;

		for (int i = 0; i < offsets.length; i++) {
			if (i == 0) {
				if (offsets[i] < 0) {
					omitRead = true;
					continue;
				}
				if (offsets[i] > maxCpgDist) {
					omitRead = true;
					break;
				}
				cpgDistRow.put(i, new Pair<Integer, Double>((int)(offsets[i] / bin), readStat.get(offsets[i]).getRight().getSecond()));
			} else {
				int cpgDist = (int)((offsets[i] - offsets[i - 1]) / bin);
				if (cpgDist < 0) {
					omitRead = true;
					break;
				}
				if (cpgDist * bin > maxCpgDist) {
					omitRead = true;
					break;
				}
				cpgDistRow.put(i, new Pair<Integer, Double>(cpgDist, readStat.get(offsets[i]).getRight().getSecond()));
			}
			matrixRow.add(readStat.get(offsets[i]).getMiddle());
			locRow.add(readStat.get(offsets[i]).getRight().getFirst());
		}

		if (omitRead) return null;
		if (matrixRow.size() < miniDataPoints || matrixRow.size() > maxCpgs) return null;

		ArrayList<Integer> observedRow = new ArrayList<Integer>();
		for (int i = 0; i < offsets.length; i++) {
			if (i == 0 && offsets[0] < 0) continue;
			if (i > 0 && (int)((offsets[i] - offsets[i - 1])) < 0) continue;
			if (readStat.get(offsets[i]).getLeft().equalsIgnoreCase("u")) {
				observedRow.add(0);
			} else if (readStat.get(offsets[i]).getLeft().equalsIgnoreCase("m")) {
				observedRow.add(1);
			}
		}

		return new AssembledFragment(cpgDistRow, matrixRow, locRow, observedRow);
	}

	private SummaryStatistics[] collectStats(String matrixFile,
											 HashMap<String, IntervalTree<Integer>> overlapLoc,
											 HashMap<String, IntervalTree<Integer>> excludeLoc) throws IOException {
		SummaryStatistics[] stats = new SummaryStatistics[3];
		for (int i = 0; i < 3; i++) {
			stats[i] = new SummaryStatistics();
		}

		GZIPInputStream gzipInputStream = null;
		BufferedReader br;
		if (matrixFile.endsWith(".gz")) {
			gzipInputStream = new GZIPInputStream(new FileInputStream(matrixFile));
			br = new BufferedReader(new InputStreamReader(gzipInputStream));
		} else {
			br = new BufferedReader(new FileReader(matrixFile));
		}

		String line;
		while ((line = br.readLine()) != null) {
			ParsedRow row = parseLine(line, overlapLoc, excludeLoc);
			if (row == null) continue;
			stats[0].addValue(row.fragLen);
			stats[1].addValue(row.coverage);
			stats[2].addValue(row.distToCenter);
		}

		if (matrixFile.endsWith(".gz")) {
			gzipInputStream.close();
		}
		br.close();

		return stats;
	}

	private void logFeatureStats(SummaryStatistics[] stats) {
		final String[] featureNames = new String[]{"FragLen", "Norm_Frag_cov", "DistToCenter"};
		for (int i = 0; i < stats.length; i++) {
			String featureName = i < featureNames.length ? featureNames[i] : ("Feature" + i);
			logFeatureStat(i, featureName, stats[i]);
		}
	}

	private void logFeatureStat(int index, String featureName, SummaryStatistics stat) {
		long n = stat.getN();
		if (n == 0) {
			log.info("Feature " + index + " (" + featureName + "): n=0 (no usable points after filtering)");
			return;
		}
		String sd = n > 1 ? Double.toString(stat.getStandardDeviation()) : "NaN";
		log.info("Feature " + index + " (" + featureName + "): n=" + n
				+ ", min=" + stat.getMin()
				+ ", max=" + stat.getMax()
				+ ", mean=" + stat.getMean()
				+ ", sd=" + sd);
	}

	/**
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		FinaleMe cnh = new FinaleMe();
		cnh.doMain(args);
	}
	

	public void doMain(String[] args)
			throws Exception {

					CmdLineParser parser = new CmdLineParser(this);
					try
						{
							if(help || args.length < 3) throw new CmdLineException(parser, USAGE, new Throwable());
							parser.parseArgument(args);
							if (patOutput && cpgIndexFile == null) throw new CmdLineException(parser, "-patOutput requires -cpgIndexFile");
							if (bwOutput && chromSizeFile == null) throw new CmdLineException(parser, "-bwOutput requires -chromSizeFile");
							if (aucMode && bwOutput) throw new CmdLineException(parser, "-bwOutput is not supported with -aucMode");
							
						
						}
					catch (CmdLineException e)
					{
						System.err.println(e.getMessage());
						// print the list of available options
						parser.printUsage(System.err);
						System.err.println();
						return;
					}

					String modelFile = arguments.get(0);
					String inputFile = arguments.get(1);
					String outputFile = arguments.get(2);
						log.info("Using " + resolveThreadCount() + " threads for FinaleMe parallel sections ...");
					initiate();

					if (decodeModeOnly && !aucMode) {
						// Streaming decode path: bounded memory
						CpgIndex cpgIndex = null;
						if (patOutput) {
							cpgIndex = loadCpgIndex(cpgIndexFile);
						}
						LinkedHashMap<String, Integer> chromOrder = null;
						if (bwOutput) {
							chromOrder = loadChromOrder(chromSizeFile);
						}
						decodeOnlyStreaming(inputFile, modelFile, outputFile, cpgIndex, chromOrder);
					} else {
						// Original path: training + decode, aucMode
						MatrixObj matrixObj = processMatrixFile(inputFile);
						// cpgDistFreq is never read after construction; free it immediately
						matrixObj.cpgDistFreq = null;
						CpgIndex cpgIndex = null;
						if (patOutput) {
							cpgIndex = loadCpgIndex(cpgIndexFile);
						}
						LinkedHashMap<String, Integer> chromOrder = null;
						if (bwOutput) {
							chromOrder = loadChromOrder(chromSizeFile);
						}
						if(aucMode){
							aucMode(matrixObj, modelFile, outputFile);
						}else{
							if(!decodeModeOnly){
								int miniDataPointsPre = miniDataPoints;
								if(miniDataPoints < 2){
									miniDataPoints = 2;
									MatrixObj matrixObj2 = processMatrixFile(inputFile);
									matrixObj2.cpgDistFreq = null;
									trainHmm(matrixObj2, modelFile);
									matrixObj2 = null; // allow GC of second copy
								}else{
									trainHmm(matrixObj, modelFile);
								}
								miniDataPoints = miniDataPointsPre;
							}
							// Free training-only fields not needed for decode
							matrixObj.matrixU = null;
							matrixObj.matrixM = null;
							matrixObj.pi = null;
							matrixObj.a = null;
							decodeHmm(matrixObj, modelFile, outputFile, inputFile, false, cpgIndex, chromOrder);
						}
					}

					finish();
					

	}
	
	private MatrixObj processMatrixFile(String matrixFile) throws FileNotFoundException, IOException, FileFormatException{
		ArrayList<Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>>> matrix = new ArrayList<Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>>>();
		HashMap<String, TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>>> matrixProcess = new HashMap<String, TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>>>();
		ArrayList<Pair<Integer, Double>> cpgDistFreq = new ArrayList<Pair<Integer, Double>>();
		HashMap<String,IntervalTree<Integer>> overlapLoc = null;
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
		
		
		// Single-pass: read matrix file once, accumulate stats and store raw parsed rows
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
			for(int i = 0; i < 3; i++){
				stats[i] = new SummaryStatistics();
			}

			// Store raw parsed row data for second-phase processing
			// Each entry: [readName, methyStat, loc, offset, methyPrior, fragLen, coverage, distToCenter]
			ArrayList<Object[]> rawRows = new ArrayList<Object[]>();

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
				double fragLen = Double.parseDouble(splitLines[4]);
				double coverage = Double.parseDouble(splitLines[7]);
				double DistToCenter = fragLen/2-Double.parseDouble(splitLines[10])+0.5;
				stats[0].addValue(fragLen);
				stats[1].addValue(coverage);
				stats[2].addValue(DistToCenter);

				// Adjust methyPrior for boundary values
				if(Double.compare(methyPrior, 100.0)==0){
					methyPrior -= 0.01;
				}else if(Double.compare(methyPrior, 0.0)==0){
					methyPrior += 0.01;
				}
				methyPrior /= 100;
				if(Double.isNaN(methyPrior)){
					continue;
				}

				String loc = chr + ":" + start + ":" + end;
				String readName = splitLines[3];
				rawRows.add(new Object[]{readName, splitLines[6], loc, offset, methyPrior, fragLen, coverage, DistToCenter});
			}
			if(matrixFile.endsWith(".gz")){
				gzipInputStream.close();
			}
			br.close();

		logFeatureStats(stats);

		// Second phase: iterate in-memory raw rows, apply z-score normalization
		for(Object[] row : rawRows){
			String readName = (String) row[0];
			String methyStat = (String) row[1];
			String loc = (String) row[2];
			Integer offset = (Integer) row[3];
			Double methyPrior = (Double) row[4];
			double fragLen = (Double) row[5];
			double coverage = (Double) row[6];
			double DistToCenter = (Double) row[7];

			if(covOutlier > 0 && ((coverage-stats[1].getMean())/stats[1].getStandardDeviation() > covOutlier ||
					(fragLen-stats[0].getMean())/stats[0].getStandardDeviation() > covOutlier ||
					(DistToCenter-stats[2].getMean())/stats[2].getStandardDeviation() > covOutlier)){
				continue;
			}
			double[] value;

			if(lowCoverage){
				value = new double[]{
						(fragLen-stats[0].getMean())/stats[0].getStandardDeviation(),
						(DistToCenter-stats[2].getMean())/stats[2].getStandardDeviation(),
				};
			}else{
				value = new double[]{
						(fragLen-stats[0].getMean())/stats[0].getStandardDeviation(),
						(coverage-stats[1].getMean())/stats[1].getStandardDeviation(),
						(DistToCenter-stats[2].getMean())/stats[2].getStandardDeviation(),
				};
			}

			ObservationVector vector = new ObservationVector(value);

			points++;
			if(matrixProcess.containsKey(readName)){
				TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>> readStat = matrixProcess.get(readName);
				if(!readStat.containsKey(offset)){
					readStat.put(offset, Triple.of(methyStat, vector, new Pair<String, Double>(loc, methyPrior)));
				}
				matrixProcess.put(readName, readStat);
			}else{
				TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>> readStat =  new TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>>();
				readStat.put(offset, Triple.of(methyStat, vector, new Pair<String, Double>(loc, methyPrior)));
				matrixProcess.put(readName, readStat);
			}
		}
		rawRows = null; // allow GC
			log.info("Number of point in total is loaded : " + points);

			ArrayList<ObservationVector> matrixU = new ArrayList<ObservationVector>();
			ArrayList<ObservationVector> matrixM = new ArrayList<ObservationVector>();
			ArrayList<ArrayList<Integer>> matrixObserved = new ArrayList<ArrayList<Integer>>();
			TreeMap<Integer, Long[]> pi = new TreeMap<Integer, Long[]>();
			TreeMap<Integer, Long[]> aij = new TreeMap<Integer, Long[]>();
			
			double observedMaxCpgNum = 0.0;
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
							}
							
							matrixRowObserved.add(1);
						}
					}
					matrixObserved.add(matrixRowObserved);
				}
				double cpgDense = (double)matrixRow.size();
				if(cpgDense > observedMaxCpgNum){
					observedMaxCpgNum = cpgDense;
				}
				if(cpgNumClip < 0 && cpgDense > this.maxCpgNum){
					this.maxCpgNum = cpgDense;
				}
			}
			matrixProcess.clear();
			matrixProcess = null; // allow GC of the large per-read HashMap
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
			log.info("maximum number of cpg in fragment is: " + observedMaxCpgNum);
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
		hmm.setBayesianFactor(bayesianFactor);
		hmm.setMethyState(this.methylatedState);
		hmm.setMaxCpgNum(cpgNumClip < 0 ? maxCpgNum : cpgNumClip);
		hmm.setMinCpgNum(1);
		
		int nThreads = resolveThreadCount();
		BaumWelchBayesianNhmmV5ScaledLearner bwl = new BaumWelchBayesianNhmmV5ScaledLearner(nThreads);

		BayesianNhmmV5<ObservationVector> prevHmm = null;
		try {
			prevHmm = hmm.clone();
		} catch (CloneNotSupportedException e1) {
			e1.printStackTrace();
		}
		
		// This object measures the distance between two HMMs
		KullbackLeiblerDistanceBayesianNhmmV5Calculator klc = 
			new KullbackLeiblerDistanceBayesianNhmmV5Calculator(matrix, nThreads);
		
		double distance = Double.MAX_VALUE;
		double distancePre = 0.01;
		int a = 0;
		// Incrementally improve the solution
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
		System.out.println("Resulting NHMM:\n" + hmm);
		this.methylatedState = hmm.getMethyState(lowCoverage);
		hmm.setMaxCpgNum(cpgNumClip < 0 ? maxCpgNum : cpgNumClip);
		hmm.setMinCpgNum(1);
		ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(modelFile));
		oos.writeObject(hmm);
		oos.close();

	}
	

	//decoding HMM
	private double decodeHmm(MatrixObj matrixObj, String hmmFile, String outputFile, String inputFile, boolean reestimate, CpgIndex cpgIndex, LinkedHashMap<String, Integer> chromOrder) throws Exception{
		System.out.println("\nDecoding ...\n");
		List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> matrix = new ArrayList<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>>();
		List<ArrayList<String>> matrixLoc = new ArrayList<ArrayList<String>>();
		for(Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>> row : matrixObj.matrix){
			matrix.add(new Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>(row.getLeft(), row.getMiddle()));
			matrixLoc.add(row.getRight());
		}
		
		ArrayList<ArrayList<Integer>> matrixObserved = matrixObj.matrixObserved;

	    BayesianNhmmV5<ObservationVector> hmm = loadHmmModel(hmmFile);
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
		HashMap<String, Integer> patRecords = null;
		long skippedFragments = 0;
		if (patOutput && cpgIndex != null) {
			patRecords = new HashMap<String, Integer>();
		}

		double likelihood = 0;
		double likelihoodWithMethy = 0;

		// Parallelize Viterbi decoding across fragments
		int nThreads = resolveThreadCount();
		ExecutorService executor = Executors.newFixedThreadPool(nThreads);

		// Each task returns: [hiddenState, lnProb, lnProbWithMethy, fragmentIndex]
		List<Future<Object[]>> futures = new ArrayList<Future<Object[]>>(matrix.size());
		for(int j=0; j < matrix.size(); j++){
			final int idx = j;
			final Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>> seq = matrix.get(j);
			futures.add(executor.submit(new Callable<Object[]>() {
				@Override
				public Object[] call() {
					ViterbiBayesianNhmmV5Calculator vb = new ViterbiBayesianNhmmV5Calculator(seq, hmm, methylatedState, decodeP);
					int[] hiddenState = vb.stateSequence();
					double prob = vb.lnProbability();
					double probWithMethy = vb.lnProbability(true);
					return new Object[]{hiddenState, prob, probWithMethy, idx};
				}
			}));
		}

		// Aggregate results sequentially to maintain deterministic output
		for(Future<Object[]> future : futures){
			Object[] result = future.get();
			int[] hiddenState = (int[]) result[0];
			double prob = (Double) result[1];
			double probWithMethy = (Double) result[2];
			int j = (Integer) result[3];

			if(!Double.isNaN(prob) && !Double.isInfinite(prob)){
				likelihood += prob/hiddenState.length;
			}
			if(!Double.isNaN(probWithMethy) && !Double.isInfinite(probWithMethy)){
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

			if (patRecords != null) {
				if (locRow.size() != hiddenState.length) {
					skippedFragments++;
					continue;
				}
				String fragmentChr = null;
				int startCpgIndex = -1;
				int previousCpgIndex = -1;
				boolean contiguous = true;
				StringBuilder pattern = new StringBuilder(hiddenState.length);
				for (int i = 0; i < hiddenState.length; i++) {
					LocTuple locTuple = parseLoc(locRow.get(i));
					if (locTuple == null) {
						contiguous = false;
						break;
					}
					int globalIndex = cpgIndex.getGlobalIndex(locTuple.chr, locTuple.start);
					if (globalIndex < 0) {
						contiguous = false;
						break;
					}
					if (i == 0) {
						fragmentChr = locTuple.chr;
						startCpgIndex = globalIndex;
					} else {
						if (!fragmentChr.equals(locTuple.chr) || globalIndex != previousCpgIndex + 1) {
							contiguous = false;
							break;
						}
					}
					pattern.append(hiddenState[i] == methylatedState ? 'C' : 'T');
					previousCpgIndex = globalIndex;
				}
				if (contiguous && pattern.length() > 0) {
					String key = fragmentChr + "\t" + startCpgIndex + "\t" + pattern.toString();
					Integer countKey = patRecords.get(key);
					patRecords.put(key, countKey == null ? 1 : countKey + 1);
				} else {
					skippedFragments++;
				}
			}
		}
		executor.shutdown();

		double[][] predData = writeDecodePredictionOutputWithTabix(methySummary, outputFile, chromOrder);

		if (patRecords != null) {
			String patFile = derivePatFile(outputFile);
			ArrayList<PatRecord> sortedPatRecords = toSortedPatRecords(patRecords);
			writePatOutput(patFile, sortedPatRecords);
			writeBetaOutput(deriveBetaFile(patFile), cpgIndex, sortedPatRecords);
			if (skippedFragments > 0) {
				log.info("Skipped " + skippedFragments + " fragments for .pat/.beta output due to missing/non-consecutive CpG index mapping.");
			}
		}
		if (bwOutput && chromOrder != null) {
			writeDecodeBigWigOutputs(methySummary, outputFile, chromOrder);
		}

		PearsonsCorrelation pearson =  new PearsonsCorrelation(predData);
		System.out.println(pearson.getCorrelationMatrix());
		System.out.println(pearson.getCorrelationPValues());
		System.out.println("counted point in total: " + count + "\tCorrect predicted:" + countCorrect + "\tPerc:" + 100*(double)countCorrect/(double)count + "%");
		System.out.println("counted point in methy: " + countMethy + "\tCorrect predicted:" + countMethyCorrect + "\tPerc:" + 100*(double)countMethyCorrect/(double)countMethy + "%");
		System.out.println("counted point in unmethy: " + countUnmethy + "\tCorrect predicted:" + countUnmethyCorrect + "\tPerc:" + 100*(double)countUnmethyCorrect/(double)countUnmethy + "%");
		System.out.println("methyState " + methylatedState + "\tLikelihood is:" + likelihood);	
		System.out.println("methyState " + methylatedState + "\tLikelihoodWithMethy is:" + likelihoodWithMethy);	
		return likelihood;

	}

	/**
	 * Streaming decode-only path: reads the input file twice (once for stats, once for
	 * streaming decode) so that memory stays bounded regardless of coverage depth.
	 */
	private void decodeOnlyStreaming(String inputFile, String modelFile, String outputFile,
									CpgIndex cpgIndex, LinkedHashMap<String, Integer> chromOrder) throws Exception {
		System.out.println("\nStreaming decode-only mode ...\n");

		// Load region/exclude intervals
		HashMap<String, IntervalTree<Integer>> overlapLoc = loadIntervalFile(region);
		HashMap<String, IntervalTree<Integer>> excludeLoc = loadIntervalFile(exclude);

		// Phase 1: Collect stats for z-score normalization
		log.info("Phase 1: Collecting feature statistics ...");
		SummaryStatistics[] stats = collectStats(inputFile, overlapLoc, excludeLoc);
		logFeatureStats(stats);

		// Load HMM model
		BayesianNhmmV5<ObservationVector> hmm = loadHmmModel(modelFile);
		hmm.setBayesianFactor(bayesianFactor);
		hmm.setMethyState(this.methylatedState);
		// When cpgNumClip < 0 in decode-only mode, use maxCpgs as safe upper bound
		hmm.setMaxCpgNum(cpgNumClip < 0 ? maxCpgs : cpgNumClip);
		hmm.setMinCpgNum(1);

		// Phase 2: Streaming decode
		log.info("Phase 2: Streaming decode ...");

		// Accumulators (persist across batches)
		HashMap<String, Pair<Pair<Integer,Integer>, Pair<Integer,Integer>>> methySummary =
			new HashMap<String, Pair<Pair<Integer,Integer>, Pair<Integer,Integer>>>();
		HashMap<String, Integer> patRecords = null;
		long skippedFragments = 0;
		if (patOutput && cpgIndex != null) {
			patRecords = new HashMap<String, Integer>();
		}

		long count = 0;
		long countCorrect = 0;
		long countMethy = 0;
		long countMethyCorrect = 0;
		long countUnmethy = 0;
		long countUnmethyCorrect = 0;
		double likelihood = 0;
		double likelihoodWithMethy = 0;

		// Thread pool with bounded queue for backpressure
		int nThreads = resolveThreadCount();
		ExecutorService executor = new ThreadPoolExecutor(
			nThreads, nThreads, 0L, TimeUnit.MILLISECONDS,
			new ArrayBlockingQueue<Runnable>(nThreads * 2),
			new ThreadPoolExecutor.CallerRunsPolicy());

		// Batch processing state
		final int BATCH_SIZE = 500000;
		HashMap<String, TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>>> currentBatch =
			new HashMap<String, TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>>>();
		int batchReadCount = 0;
		String lastReadName = null;

		// Open input for second pass
		GZIPInputStream gzipInputStream = null;
		BufferedReader br;
		if (inputFile.endsWith(".gz")) {
			gzipInputStream = new GZIPInputStream(new FileInputStream(inputFile));
			br = new BufferedReader(new InputStreamReader(gzipInputStream));
		} else {
			br = new BufferedReader(new FileReader(inputFile));
		}

		String line;
		while ((line = br.readLine()) != null) {
			ParsedRow row = parseLine(line, overlapLoc, excludeLoc);
			if (row == null) continue;

			// Apply covOutlier filter
			if (covOutlier > 0 && ((row.coverage - stats[1].getMean()) / stats[1].getStandardDeviation() > covOutlier ||
					(row.fragLen - stats[0].getMean()) / stats[0].getStandardDeviation() > covOutlier ||
					(row.distToCenter - stats[2].getMean()) / stats[2].getStandardDeviation() > covOutlier)) {
				continue;
			}

			// Z-score normalize
			double[] value;
			if (lowCoverage) {
				value = new double[]{
					(row.fragLen - stats[0].getMean()) / stats[0].getStandardDeviation(),
					(row.distToCenter - stats[2].getMean()) / stats[2].getStandardDeviation(),
				};
			} else {
				value = new double[]{
					(row.fragLen - stats[0].getMean()) / stats[0].getStandardDeviation(),
					(row.coverage - stats[1].getMean()) / stats[1].getStandardDeviation(),
					(row.distToCenter - stats[2].getMean()) / stats[2].getStandardDeviation(),
				};
			}
			ObservationVector vector = new ObservationVector(value);
			points++;

			// Track new read names for batch counting
			if (!row.readName.equals(lastReadName) && lastReadName != null) {
				batchReadCount++;
			}
			lastReadName = row.readName;

			// Group by read name
			TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>> readStat = currentBatch.get(row.readName);
			if (readStat == null) {
				readStat = new TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>>();
				currentBatch.put(row.readName, readStat);
			}
			if (!readStat.containsKey(row.offset)) {
				readStat.put(row.offset, Triple.of(row.methyStat, vector, new Pair<String, Double>(row.loc, row.methyPrior)));
			}

			// Flush batch when threshold reached
			if (batchReadCount >= BATCH_SIZE) {
				// Keep the last read (may have more lines coming)
				TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>> lastReadStat = currentBatch.remove(lastReadName);

				// Process and decode all complete reads in this batch
				long[] batchCounters = decodeBatch(currentBatch, hmm, executor, cpgIndex, methySummary, patRecords);
				count += batchCounters[0];
				countCorrect += batchCounters[1];
				countMethy += batchCounters[2];
				countMethyCorrect += batchCounters[3];
				countUnmethy += batchCounters[4];
				countUnmethyCorrect += batchCounters[5];
				likelihood += Double.longBitsToDouble(batchCounters[6]);
				likelihoodWithMethy += Double.longBitsToDouble(batchCounters[7]);
				skippedFragments += batchCounters[8];

				currentBatch.clear();
				if (lastReadStat != null) {
					currentBatch.put(lastReadName, lastReadStat);
				}
				batchReadCount = 0;
			}
		}

		// Process final batch
		if (!currentBatch.isEmpty()) {
			long[] batchCounters = decodeBatch(currentBatch, hmm, executor, cpgIndex, methySummary, patRecords);
			count += batchCounters[0];
			countCorrect += batchCounters[1];
			countMethy += batchCounters[2];
			countMethyCorrect += batchCounters[3];
			countUnmethy += batchCounters[4];
			countUnmethyCorrect += batchCounters[5];
			likelihood += Double.longBitsToDouble(batchCounters[6]);
			likelihoodWithMethy += Double.longBitsToDouble(batchCounters[7]);
			skippedFragments += batchCounters[8];
		}

		if (inputFile.endsWith(".gz")) {
			gzipInputStream.close();
		}
		br.close();

		executor.shutdown();
		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS);

		log.info("Number of point in total is loaded : " + points);

		// Phase 3: Write outputs
		double[][] predData = writeDecodePredictionOutputWithTabix(methySummary, outputFile, chromOrder);

		if (patRecords != null) {
			String patFile = derivePatFile(outputFile);
			ArrayList<PatRecord> sortedPatRecords = toSortedPatRecords(patRecords);
			writePatOutput(patFile, sortedPatRecords);
			writeBetaOutput(deriveBetaFile(patFile), cpgIndex, sortedPatRecords);
			if (skippedFragments > 0) {
				log.info("Skipped " + skippedFragments + " fragments for .pat/.beta output due to missing/non-consecutive CpG index mapping.");
			}
		}
		if (bwOutput && chromOrder != null) {
			writeDecodeBigWigOutputs(methySummary, outputFile, chromOrder);
		}

		PearsonsCorrelation pearson = new PearsonsCorrelation(predData);
		System.out.println(pearson.getCorrelationMatrix());
		System.out.println(pearson.getCorrelationPValues());
		System.out.println("counted point in total: " + count + "\tCorrect predicted:" + countCorrect + "\tPerc:" + 100 * (double) countCorrect / (double) count + "%");
		System.out.println("counted point in methy: " + countMethy + "\tCorrect predicted:" + countMethyCorrect + "\tPerc:" + 100 * (double) countMethyCorrect / (double) countMethy + "%");
		System.out.println("counted point in unmethy: " + countUnmethy + "\tCorrect predicted:" + countUnmethyCorrect + "\tPerc:" + 100 * (double) countUnmethyCorrect / (double) countUnmethy + "%");
		System.out.println("methyState " + methylatedState + "\tLikelihood is:" + likelihood);
		System.out.println("methyState " + methylatedState + "\tLikelihoodWithMethy is:" + likelihoodWithMethy);
	}

	/**
	 * Decode one batch of reads: assemble fragments, run Viterbi in parallel,
	 * accumulate results into methySummary and patRecords.
	 * Returns counters: [count, countCorrect, countMethy, countMethyCorrect,
	 *                    countUnmethy, countUnmethyCorrect,
	 *                    likelihoodBits, likelihoodWithMethyBits, skippedFragments]
	 */
	private long[] decodeBatch(
			HashMap<String, TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>>> batch,
			BayesianNhmmV5<ObservationVector> hmm,
			ExecutorService executor,
			CpgIndex cpgIndex,
			HashMap<String, Pair<Pair<Integer,Integer>, Pair<Integer,Integer>>> methySummary,
			HashMap<String, Integer> patRecords) throws Exception {

		// Assemble fragments
		ArrayList<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> fragments =
			new ArrayList<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>>();
		ArrayList<ArrayList<String>> fragmentLocs = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<Integer>> fragmentObserved = new ArrayList<ArrayList<Integer>>();

		for (TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>> readStat : batch.values()) {
			AssembledFragment frag = assembleFragment(readStat);
			if (frag == null) continue;
			fragments.add(new Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>(
				frag.cpgDistRow, (List<ObservationVector>) (List<?>) frag.matrixRow));
			fragmentLocs.add(frag.locRow);
			fragmentObserved.add(frag.observedRow);
		}

		// Submit Viterbi tasks
		List<Future<Object[]>> futures = new ArrayList<Future<Object[]>>(fragments.size());
		for (int j = 0; j < fragments.size(); j++) {
			final int idx = j;
			final Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>> seq = fragments.get(j);
			futures.add(executor.submit(new Callable<Object[]>() {
				@Override
				public Object[] call() {
					ViterbiBayesianNhmmV5Calculator vb = new ViterbiBayesianNhmmV5Calculator(seq, hmm, methylatedState, decodeP);
					int[] hiddenState = vb.stateSequence();
					double prob = vb.lnProbability();
					double probWithMethy = vb.lnProbability(true);
					return new Object[]{hiddenState, prob, probWithMethy, idx};
				}
			}));
		}

		// Aggregate results
		long batchCount = 0, batchCorrect = 0, batchMethy = 0, batchMethyCorrect = 0;
		long batchUnmethy = 0, batchUnmethyCorrect = 0, batchSkipped = 0;
		double batchLikelihood = 0, batchLikelihoodWithMethy = 0;

		for (Future<Object[]> future : futures) {
			Object[] result = future.get();
			int[] hiddenState = (int[]) result[0];
			double prob = (Double) result[1];
			double probWithMethy = (Double) result[2];
			int j = (Integer) result[3];

			if (!Double.isNaN(prob) && !Double.isInfinite(prob)) {
				batchLikelihood += prob / hiddenState.length;
			}
			if (!Double.isNaN(probWithMethy) && !Double.isInfinite(probWithMethy)) {
				batchLikelihoodWithMethy += probWithMethy / hiddenState.length;
			}

			Integer[] observedState = fragmentObserved.get(j).toArray(new Integer[fragmentObserved.get(j).size()]);
			ArrayList<String> locRow = fragmentLocs.get(j);

			if (hiddenState.length != observedState.length) {
				throw new IllegalArgumentException("HiddenState Length does not match with observed state length");
			}
			for (int i = 0; i < hiddenState.length; i++) {
				String loc = locRow.get(i);

				int methyPredict = 0;
				int unmethyPredict = 0;
				if (hiddenState[i] == methylatedState) {
					methyPredict++;
				} else if (hiddenState[i] == (1 - methylatedState)) {
					unmethyPredict++;
				}

				if (methySummary.containsKey(loc)) {
					Pair<Pair<Integer,Integer>, Pair<Integer,Integer>> tmp = methySummary.get(loc);
					tmp = new Pair<Pair<Integer,Integer>, Pair<Integer,Integer>>(
						new Pair<Integer,Integer>(tmp.getFirst().getFirst() + methyPredict, tmp.getFirst().getSecond() + (methyPredict + unmethyPredict)),
						new Pair<Integer,Integer>(tmp.getSecond().getFirst() + observedState[i], tmp.getSecond().getSecond() + 1));
					methySummary.put(loc, tmp);
				} else {
					methySummary.put(loc, new Pair<Pair<Integer,Integer>, Pair<Integer,Integer>>(
						new Pair<Integer,Integer>(methyPredict, (methyPredict + unmethyPredict)),
						new Pair<Integer,Integer>(observedState[i], 1)));
				}

				if (observedState[i] == 0) {
					batchUnmethy++;
					if (hiddenState[i] == (1 - methylatedState)) {
						batchCorrect++;
						batchUnmethyCorrect++;
					}
				} else {
					batchMethy++;
					if (hiddenState[i] == methylatedState) {
						batchCorrect++;
						batchMethyCorrect++;
					}
				}
				batchCount++;
			}

			if (patRecords != null) {
				if (locRow.size() != hiddenState.length) {
					batchSkipped++;
					continue;
				}
				String fragmentChr = null;
				int startCpgIndex = -1;
				int previousCpgIndex = -1;
				boolean contiguous = true;
				StringBuilder pattern = new StringBuilder(hiddenState.length);
				for (int i = 0; i < hiddenState.length; i++) {
					LocTuple locTuple = parseLoc(locRow.get(i));
					if (locTuple == null) {
						contiguous = false;
						break;
					}
					int globalIndex = cpgIndex.getGlobalIndex(locTuple.chr, locTuple.start);
					if (globalIndex < 0) {
						contiguous = false;
						break;
					}
					if (i == 0) {
						fragmentChr = locTuple.chr;
						startCpgIndex = globalIndex;
					} else {
						if (!fragmentChr.equals(locTuple.chr) || globalIndex != previousCpgIndex + 1) {
							contiguous = false;
							break;
						}
					}
					pattern.append(hiddenState[i] == methylatedState ? 'C' : 'T');
					previousCpgIndex = globalIndex;
				}
				if (contiguous && pattern.length() > 0) {
					String key = fragmentChr + "\t" + startCpgIndex + "\t" + pattern.toString();
					Integer countKey = patRecords.get(key);
					patRecords.put(key, countKey == null ? 1 : countKey + 1);
				} else {
					batchSkipped++;
				}
			}
		}

		return new long[]{
			batchCount, batchCorrect, batchMethy, batchMethyCorrect,
			batchUnmethy, batchUnmethyCorrect,
			Double.doubleToLongBits(batchLikelihood),
			Double.doubleToLongBits(batchLikelihoodWithMethy),
			batchSkipped
		};
	}

	private CpgIndex loadCpgIndex(String cpgIndexFile) throws IOException {
		log.info("Loading CpG index file " + cpgIndexFile + " ...");
		LinkedHashMap<String, ArrayList<Integer>> chrToPositions = new LinkedHashMap<String, ArrayList<Integer>>();
		InputStream input = null;
		BufferedReader br = null;
		try {
			if (cpgIndexFile.endsWith(".gz")) {
				input = new GZIPInputStream(new FileInputStream(cpgIndexFile));
			} else {
				input = new FileInputStream(cpgIndexFile);
			}
			br = new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8));
			String line;
			while ((line = br.readLine()) != null) {
				if (line.isEmpty() || line.startsWith("#")) {
					continue;
				}
				String[] splitLines = line.split("\t");
				if (splitLines.length < 2) {
					continue;
				}
				String chr = splitLines[0];
				int start;
				try {
					start = Integer.parseInt(splitLines[1]);
				} catch (NumberFormatException e) {
					continue;
				}
				ArrayList<Integer> positions = chrToPositions.get(chr);
				if (positions == null) {
					positions = new ArrayList<Integer>();
					chrToPositions.put(chr, positions);
				}
				positions.add(start);
			}
		} finally {
			if (br != null) {
				br.close();
			}
		}

		LinkedHashMap<String, int[]> chrPositions = new LinkedHashMap<String, int[]>();
		HashMap<String, Integer> chrOffsets = new HashMap<String, Integer>();
		int totalSites = 0;
		for (Map.Entry<String, ArrayList<Integer>> entry : chrToPositions.entrySet()) {
			ArrayList<Integer> positions = entry.getValue();
			Collections.sort(positions);
			int uniqueCount = 0;
			int prev = Integer.MIN_VALUE;
			for (int pos : positions) {
				if (uniqueCount == 0 || pos != prev) {
					uniqueCount++;
					prev = pos;
				}
			}
			int[] posArray = new int[uniqueCount];
			int idx = 0;
			prev = Integer.MIN_VALUE;
			for (int pos : positions) {
				if (idx == 0 || pos != prev) {
					posArray[idx++] = pos;
					prev = pos;
				}
			}
			chrOffsets.put(entry.getKey(), totalSites);
			chrPositions.put(entry.getKey(), posArray);
			totalSites += posArray.length;
		}
		log.info("Loaded CpG index with " + totalSites + " sites across " + chrPositions.size() + " chromosomes.");
		return new CpgIndex(chrPositions, chrOffsets, totalSites);
	}

	private String derivePatFile(String outputFile) {
		if (outputFile.endsWith(".bed.gz")) {
			return outputFile.substring(0, outputFile.length() - ".bed.gz".length()) + ".pat.gz";
		}
		return outputFile + ".pat.gz";
	}

	private String deriveBetaFile(String patFile) {
		if (patFile.endsWith(".pat.gz")) {
			return patFile.substring(0, patFile.length() - ".pat.gz".length()) + ".beta";
		}
		return patFile + ".beta";
	}

	private ArrayList<PatRecord> toSortedPatRecords(HashMap<String, Integer> patRecords) {
		ArrayList<PatRecord> records = new ArrayList<PatRecord>(patRecords.size());
		for (Map.Entry<String, Integer> entry : patRecords.entrySet()) {
			String[] fields = entry.getKey().split("\t", 3);
			if (fields.length < 3) {
				continue;
			}
			int startCpgIndex;
			try {
				startCpgIndex = Integer.parseInt(fields[1]);
			} catch (NumberFormatException e) {
				continue;
			}
			records.add(new PatRecord(fields[0], startCpgIndex, fields[2], entry.getValue()));
		}
		Collections.sort(records, (a, b) -> {
			int cmp = Integer.compare(a.startCpgIndex, b.startCpgIndex);
			if (cmp != 0) {
				return cmp;
			}
			cmp = a.chr.compareTo(b.chr);
			if (cmp != 0) {
				return cmp;
			}
			return a.pattern.compareTo(b.pattern);
		});
		return records;
	}

	private void writePatOutput(String patFile, List<PatRecord> records) throws IOException {
		File patPath = new File(patFile);
		log.info("Writing UXM .pat.gz output: " + patPath.getPath());
		try (BlockCompressedOutputStream patOut = new BlockCompressedOutputStream(patPath)) {
			for (PatRecord record : records) {
				String line = record.chr + "\t" + record.startCpgIndex + "\t" + record.pattern + "\t" + record.count + "\n";
				patOut.write(line.getBytes(StandardCharsets.UTF_8));
			}
		}
	}

	private void writeBetaOutput(String betaFile, CpgIndex cpgIndex, List<PatRecord> records) throws IOException {
		log.info("Writing UXM .beta output: " + betaFile);
		int[] methyCounts = new int[cpgIndex.totalSites];
		int[] totalCounts = new int[cpgIndex.totalSites];
		for (PatRecord record : records) {
			for (int i = 0; i < record.pattern.length(); i++) {
				int arrayIndex = record.startCpgIndex + i - 1;
				if (arrayIndex < 0 || arrayIndex >= cpgIndex.totalSites) {
					continue;
				}
				if (record.pattern.charAt(i) == 'C') {
					methyCounts[arrayIndex] += record.count;
				}
				totalCounts[arrayIndex] += record.count;
			}
		}

		try (BufferedOutputStream betaOut = new BufferedOutputStream(new FileOutputStream(betaFile))) {
			for (int i = 0; i < cpgIndex.totalSites; i++) {
				int methy = methyCounts[i];
				int total = totalCounts[i];
				if (total > 255) {
					methy = (int) Math.round((double) methy * 255.0 / (double) total);
					total = 255;
				}
				if (methy > total) {
					methy = total;
				}
				betaOut.write((byte) (methy & 0xFF));
				betaOut.write((byte) (total & 0xFF));
			}
		}
	}

	private LinkedHashMap<String, Integer> loadChromOrder(String chromSizeFile) throws IOException {
		log.info("Loading chromosome sizes file " + chromSizeFile + " ...");
		LinkedHashMap<String, Integer> chromOrder = new LinkedHashMap<String, Integer>();
		InputStream input = null;
		BufferedReader br = null;
		int index = 0;
		try {
			if (chromSizeFile.endsWith(".gz")) {
				input = new GZIPInputStream(new FileInputStream(chromSizeFile));
			} else {
				input = new FileInputStream(chromSizeFile);
			}
			br = new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8));
			String line;
				while ((line = br.readLine()) != null) {
					if (line.isEmpty() || line.startsWith("#")) {
						continue;
					}
					String[] splitLines = line.split("\\s+");
					if (splitLines.length < 2) {
						continue;
					}
					String chr = splitLines[0];
					try {
						int chrSize = Integer.parseInt(splitLines[1]);
						if (chrSize <= 0) {
							continue;
						}
					} catch (NumberFormatException e) {
						continue;
					}
					if (!chromOrder.containsKey(chr)) {
						chromOrder.put(chr, index++);
					}
				}
		} finally {
			if (br != null) {
				br.close();
			}
		}
		if (chromOrder.isEmpty()) {
			throw new IllegalArgumentException("No chromosome entries found in -chromSizeFile: " + chromSizeFile);
		}
		log.info("Loaded " + chromOrder.size() + " chromosome entries for bigWig sort order.");
		return chromOrder;
	}

	private ArrayList<kotlin.Pair<String, Integer>> loadChromSizesForJavaBigWigWriter(String chromSizeFile) throws IOException {
		LinkedHashMap<String, Integer> chromSizeMap = new LinkedHashMap<String, Integer>();
		InputStream input = null;
		BufferedReader br = null;
		try {
			if (chromSizeFile.endsWith(".gz")) {
				input = new GZIPInputStream(new FileInputStream(chromSizeFile));
			} else {
				input = new FileInputStream(chromSizeFile);
			}
			br = new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8));
			String line;
			while ((line = br.readLine()) != null) {
				if (line.isEmpty() || line.startsWith("#")) {
					continue;
				}
				String[] splitLines = line.split("\\s+");
				if (splitLines.length < 2) {
					continue;
				}
				String chr = splitLines[0];
				int chrSize;
				try {
					chrSize = Integer.parseInt(splitLines[1]);
				} catch (NumberFormatException e) {
					continue;
				}
				if (chrSize <= 0) {
					continue;
				}
				if (!chromSizeMap.containsKey(chr)) {
					chromSizeMap.put(chr, chrSize);
				}
			}
		} finally {
			if (br != null) {
				br.close();
			}
		}
		if (chromSizeMap.isEmpty()) {
			throw new IllegalArgumentException("No chromosome size entries found in -chromSizeFile: " + chromSizeFile);
		}
		ArrayList<kotlin.Pair<String, Integer>> chromSizes = new ArrayList<kotlin.Pair<String, Integer>>(chromSizeMap.size());
		for (Map.Entry<String, Integer> entry : chromSizeMap.entrySet()) {
			chromSizes.add(new kotlin.Pair<String, Integer>(entry.getKey(), entry.getValue()));
		}
		return chromSizes;
	}

	private double[][] writeDecodePredictionOutputWithTabix(
			HashMap<String, Pair<Pair<Integer, Integer>, Pair<Integer, Integer>>> methySummary,
			String outputFile,
			LinkedHashMap<String, Integer> chromOrder) throws IOException {
		ArrayList<DecodeSummaryRow> rows = toSortedDecodeRowsForPrediction(methySummary, chromOrder);
		File bgzfOutput = new File(outputFile);
		log.info("Writing decode prediction output as BGZF: " + bgzfOutput.getPath());
		try (BlockCompressedOutputStream output = new BlockCompressedOutputStream(bgzfOutput);
				OutputStreamWriter writer = new OutputStreamWriter(output, StandardCharsets.UTF_8)) {
			writer.write("#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict\ttotal_count_predict\tmethy_perc_obs\tmethy_count_obs\ttotal_count_obs\n");
			for (DecodeSummaryRow row : rows) {
				writer.write(row.chr);
				writer.write('\t');
				writer.write(Integer.toString(row.start));
				writer.write('\t');
				writer.write(Integer.toString(row.end));
				writer.write('\t');
				writer.write(Double.toString(row.predictMethyPercent()));
				writer.write('\t');
				writer.write(Integer.toString(row.methyPred));
				writer.write('\t');
				writer.write(Integer.toString(row.totalPred));
				writer.write('\t');
				writer.write(Double.toString(row.observedMethyPercent()));
				writer.write('\t');
				writer.write(Integer.toString(row.methyObs));
				writer.write('\t');
				writer.write(Integer.toString(row.totalObs));
				writer.write('\n');
			}
		}

		if (!rows.isEmpty()) {
			createTabixIndexForDecodeOutput(bgzfOutput);
		} else {
			log.warn("Decode summary has no valid rows. Skipping tabix index generation for " + bgzfOutput.getPath());
		}

		double[][] predData = new double[rows.size()][2];
		for (int i = 0; i < rows.size(); i++) {
			predData[i][0] = rows.get(i).predictMethyPercent();
			predData[i][1] = rows.get(i).observedMethyPercent();
		}
		return predData;
	}

	private void createTabixIndexForDecodeOutput(File bgzfOutput) {
		try {
			IndexFactory.createTabixIndex(bgzfOutput, new BEDCodec(), TabixFormat.BED, null)
					.write(Tribble.tabixIndexFile(bgzfOutput).toPath());
			log.info("Wrote tabix index for decode output: " + Tribble.tabixIndexFile(bgzfOutput).getPath());
		} catch (Exception e) {
			log.warn("Failed to create tabix index for decode output " + bgzfOutput.getPath()
					+ ". You can run 'tabix -p bed " + bgzfOutput.getPath() + "' manually.", e);
		}
	}

	private void writeDecodeBigWigOutputs(HashMap<String, Pair<Pair<Integer, Integer>, Pair<Integer, Integer>>> methySummary, String outputFile, LinkedHashMap<String, Integer> chromOrder) throws Exception {
		if (methySummary == null || methySummary.isEmpty()) {
			log.info("Skip bigWig output because decode summary is empty.");
			return;
		}
		ArrayList<DecodeSummaryRow> rows = toSortedDecodeRowsForBigWig(methySummary, chromOrder);
		if (rows.isEmpty()) {
			log.info("Skip bigWig output because no valid decode intervals were parsed.");
			return;
		}

		String unknownChrom = null;
		for (DecodeSummaryRow row : rows) {
			if (!chromOrder.containsKey(row.chr)) {
				unknownChrom = row.chr;
				break;
			}
		}
		if (unknownChrom != null) {
			throw new IllegalArgumentException("Chromosome '" + unknownChrom + "' from decode output is not present in -chromSizeFile. " +
					"Consider -bwStripChrPrefix and/or -bwConvertChrMToMT if naming conventions differ.");
		}

		String base = deriveBigWigBase(outputFile);
		File methyBedGraph = new File(base + ".methy.bedgraph");
		File covBedGraph = new File(base + ".cov.bedgraph");
		File methyCountBedGraph = new File(base + ".methy_count.bedgraph");
		File methyBw = new File(base + ".methy.bw");
		File covBw = new File(base + ".cov.bw");
		File methyCountBw = new File(base + ".methy_count.bw");

		log.info("Writing temporary bedGraph files for bigWig conversion ...");
		writeDecodeBedGraph(methyBedGraph, rows, 0);
		writeDecodeBedGraph(covBedGraph, rows, 1);
		writeDecodeBedGraph(methyCountBedGraph, rows, 2);

		boolean success = false;
		try {
			try {
				runBedGraphToBigWig(methyBedGraph, methyBw);
				runBedGraphToBigWig(covBedGraph, covBw);
				runBedGraphToBigWig(methyCountBedGraph, methyCountBw);
				log.info("Wrote decode bigWig outputs with bedGraphToBigWig: " + methyBw.getPath() + ", " + covBw.getPath() + ", " + methyCountBw.getPath());
			} catch (IOException externalConvertError) {
				if (!isBedGraphToBigWigUnavailable(externalConvertError)) {
					throw externalConvertError;
				}
				log.warn("bedGraphToBigWig is unavailable ({}). Falling back to built-in Java BigWig writer (org.jetbrains.bio:big).",
						externalConvertError.getMessage());
				ArrayList<kotlin.Pair<String, Integer>> chromSizesForJavaWriter = loadChromSizesForJavaBigWigWriter(chromSizeFile);
				runBedGraphToBigWigWithJetBrainsBig(methyBedGraph, methyBw, chromSizesForJavaWriter);
				runBedGraphToBigWigWithJetBrainsBig(covBedGraph, covBw, chromSizesForJavaWriter);
				runBedGraphToBigWigWithJetBrainsBig(methyCountBedGraph, methyCountBw, chromSizesForJavaWriter);
				log.info("Wrote decode bigWig outputs with built-in Java writer: " + methyBw.getPath() + ", " + covBw.getPath() + ", " + methyCountBw.getPath());
			}
			success = true;
		} finally {
			if (success) {
				if (!methyBedGraph.delete()) log.warn("Could not delete temporary file: " + methyBedGraph.getPath());
				if (!covBedGraph.delete()) log.warn("Could not delete temporary file: " + covBedGraph.getPath());
				if (!methyCountBedGraph.delete()) log.warn("Could not delete temporary file: " + methyCountBedGraph.getPath());
			} else {
				log.warn("Keeping temporary bedGraph files after bigWig conversion failure for troubleshooting.");
			}
		}
	}

	private ArrayList<DecodeSummaryRow> toSortedDecodeRowsForPrediction(
			HashMap<String, Pair<Pair<Integer, Integer>, Pair<Integer, Integer>>> methySummary,
			LinkedHashMap<String, Integer> chromOrder) {
		return toSortedDecodeRows(methySummary, chromOrder, false);
	}

	private ArrayList<DecodeSummaryRow> toSortedDecodeRowsForBigWig(HashMap<String, Pair<Pair<Integer, Integer>, Pair<Integer, Integer>>> methySummary, LinkedHashMap<String, Integer> chromOrder) {
		return toSortedDecodeRows(methySummary, chromOrder, true);
	}

	private ArrayList<DecodeSummaryRow> toSortedDecodeRows(
			HashMap<String, Pair<Pair<Integer, Integer>, Pair<Integer, Integer>>> methySummary,
			LinkedHashMap<String, Integer> chromOrder,
			boolean normalizeChr) {
		ArrayList<DecodeSummaryRow> rows = new ArrayList<DecodeSummaryRow>(methySummary.size());
		for (Map.Entry<String, Pair<Pair<Integer, Integer>, Pair<Integer, Integer>>> entry : methySummary.entrySet()) {
			String[] locTmp = entry.getKey().split(":");
			if (locTmp.length < 3) {
				continue;
			}
			int start;
			int end;
			try {
				start = Integer.parseInt(locTmp[1]);
				end = Integer.parseInt(locTmp[2]);
			} catch (NumberFormatException e) {
				continue;
			}
			Pair<Pair<Integer, Integer>, Pair<Integer, Integer>> counts = entry.getValue();
			int methyPred = counts.getFirst().getFirst();
			int totalPred = counts.getFirst().getSecond();
			int methyObs = counts.getSecond().getFirst();
			int totalObs = counts.getSecond().getSecond();
			if (totalPred <= 0) {
				continue;
			}
			String chr = normalizeChr ? normalizeChrForBigWig(locTmp[0]) : locTmp[0];
			rows.add(new DecodeSummaryRow(chr, start, end, methyPred, totalPred, methyObs, totalObs));
		}

		Collections.sort(rows, (a, b) -> compareDecodeRows(a, b, chromOrder));
		return rows;
	}

	private int compareDecodeRows(DecodeSummaryRow a, DecodeSummaryRow b, LinkedHashMap<String, Integer> chromOrder) {
		if (chromOrder != null) {
			Integer aOrder = chromOrder.get(a.chr);
			Integer bOrder = chromOrder.get(b.chr);
			if (aOrder != null && bOrder != null) {
				int cmp = Integer.compare(aOrder, bOrder);
				if (cmp != 0) {
					return cmp;
				}
			} else if (aOrder != null) {
				return -1;
			} else if (bOrder != null) {
				return 1;
			}
		}
		int cmp = a.chr.compareTo(b.chr);
		if (cmp != 0) {
			return cmp;
		}
		cmp = Integer.compare(a.start, b.start);
		if (cmp != 0) {
			return cmp;
		}
		return Integer.compare(a.end, b.end);
	}

	private void writeDecodeBedGraph(File bedGraphFile, List<DecodeSummaryRow> rows, int metricMode) throws IOException {
		try (BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(bedGraphFile), StandardCharsets.UTF_8))) {
			for (DecodeSummaryRow row : rows) {
				String value;
				if (metricMode == 0) {
					value = Double.toString(row.predictMethyPercent());
				} else if (metricMode == 1) {
					value = Integer.toString(row.totalPred);
				} else {
					value = Integer.toString(row.methyPred);
				}
				writer.write(row.chr);
				writer.write('\t');
				writer.write(Integer.toString(row.start));
				writer.write('\t');
				writer.write(Integer.toString(row.end));
				writer.write('\t');
				writer.write(value);
				writer.write('\n');
			}
		}
	}

	private void runBedGraphToBigWig(File bedGraphFile, File bigWigFile) throws Exception {
		log.info("Converting " + bedGraphFile.getPath() + " -> " + bigWigFile.getPath());
		ProcessBuilder pb = new ProcessBuilder(bedGraphToBigWig, bedGraphFile.getPath(), chromSizeFile, bigWigFile.getPath());
		pb.redirectErrorStream(true);
		Process process = pb.start();
		StringBuilder commandOutput = new StringBuilder();
		try (BufferedReader processReader = new BufferedReader(new InputStreamReader(process.getInputStream(), StandardCharsets.UTF_8))) {
			String line;
			while ((line = processReader.readLine()) != null) {
				commandOutput.append(line).append('\n');
			}
		}
		int exitCode = process.waitFor();
		if (exitCode != 0) {
			throw new IOException("bedGraphToBigWig failed (exit " + exitCode + ") for " + bedGraphFile.getPath() + ":\n" + commandOutput.toString());
		}
	}

	private void runBedGraphToBigWigWithJetBrainsBig(File bedGraphFile, File bigWigFile, Iterable<kotlin.Pair<String, Integer>> chromSizesForJavaWriter) throws Exception {
		log.info("Converting " + bedGraphFile.getPath() + " -> " + bigWigFile.getPath() + " using built-in Java writer ...");
		byte[] bedGraphTrackHeader = "track type=bedGraph\n".getBytes(StandardCharsets.UTF_8);
		try (SequenceInputStream combinedInput = new SequenceInputStream(
				new ByteArrayInputStream(bedGraphTrackHeader),
				new FileInputStream(bedGraphFile));
				BufferedReader reader = new BufferedReader(new InputStreamReader(combinedInput, StandardCharsets.UTF_8));
				BedGraphParser parser = new BedGraphParser(reader)) {
			BigWigFile.write((Iterable<? extends WigSection>) parser, chromSizesForJavaWriter, bigWigFile.toPath());
		} catch (NoClassDefFoundError | ExceptionInInitializerError e) {
			throw new IOException("JetBrains BigWig writer classes are unavailable. Ensure dependency org.jetbrains.bio:big is on classpath, "
					+ "or provide -bedGraphToBigWig executable.", e);
		} catch (IllegalStateException e) {
			throw new IOException("Built-in Java BigWig writer could not parse temporary bedGraph file " + bedGraphFile.getPath()
					+ ". Ensure decode output is sorted and valid.", e);
		}
	}

	private boolean isBedGraphToBigWigUnavailable(IOException e) {
		if (e == null || e.getMessage() == null) {
			return false;
		}
		String messageLower = e.getMessage().toLowerCase(Locale.ROOT);
		return messageLower.contains("cannot run program")
				&& (messageLower.contains("error=2")
				|| messageLower.contains("no such file")
				|| messageLower.contains("not found"));
	}

	private String deriveBigWigBase(String outputFile) {
		if (outputFile.endsWith(".bed.gz")) {
			return outputFile.substring(0, outputFile.length() - ".bed.gz".length());
		}
		if (outputFile.endsWith(".gz")) {
			return outputFile.substring(0, outputFile.length() - ".gz".length());
		}
		if (outputFile.endsWith(".bed")) {
			return outputFile.substring(0, outputFile.length() - ".bed".length());
		}
		return outputFile;
	}

	private String normalizeChrForBigWig(String chr) {
		String normalized = chr;
		if (bwConvertChrMToMT && (normalized.equals("chrM") || normalized.equals("M") || normalized.equals("chrMT"))) {
			normalized = "MT";
		}
		if (bwStripChrPrefix && normalized.startsWith("chr")) {
			normalized = normalized.substring(3);
		}
		if (bwConvertChrMToMT && normalized.equals("M")) {
			normalized = "MT";
		}
		return normalized;
	}

	private LocTuple parseLoc(String loc) {
		if (loc == null) {
			return null;
		}
		int firstColon = loc.indexOf(':');
		int secondColon = loc.indexOf(':', firstColon + 1);
		if (firstColon < 0 || secondColon < 0) {
			return null;
		}
		String chr = loc.substring(0, firstColon);
		int start;
		try {
			start = Integer.parseInt(loc.substring(firstColon + 1, secondColon));
		} catch (NumberFormatException e) {
			return null;
		}
		return new LocTuple(chr, start);
	}
	
	
	//decoding HMM 
	private void aucMode(MatrixObj matrixObj, String hmmFile, String outputFile) throws Exception{
		System.out.println("\nROC curve ...\n");	
		List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> matrix = new ArrayList<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>>();
		
		for(Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>> row : matrixObj.matrix){
			matrix.add(new Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>(row.getLeft(), row.getMiddle()));
		
		}
		ArrayList<ArrayList<Integer>> matrixObserved = matrixObj.matrixObserved;

	    BayesianNhmmV5<ObservationVector> hmm = loadHmmModel(hmmFile);
		hmm.setBayesianFactor(bayesianFactor);
		hmm.setMethyState(this.methylatedState);
		hmm.setMaxCpgNum(cpgNumClip < 0 ? maxCpgNum : cpgNumClip);
		hmm.setMinCpgNum(1);
		boolean outFlag = false;
		double[] ps = new double[]{-1,-0.999, -0.99, -0.95, -0.9,-0.8,-0.7,-0.5,-0.1,0.0, 0.1,0.5,0.7,0.8,0.9,0.95,0.99,0.999,1.0};
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

	private static class LocTuple {
		final String chr;
		final int start;

		LocTuple(String chr, int start) {
			this.chr = chr;
			this.start = start;
		}
	}

	private static class PatRecord {
		final String chr;
		final int startCpgIndex;
		final String pattern;
		final int count;

		PatRecord(String chr, int startCpgIndex, String pattern, int count) {
			this.chr = chr;
			this.startCpgIndex = startCpgIndex;
			this.pattern = pattern;
			this.count = count;
		}
	}

	private static class DecodeSummaryRow {
		final String chr;
		final int start;
		final int end;
		final int methyPred;
		final int totalPred;
		final int methyObs;
		final int totalObs;

		DecodeSummaryRow(String chr, int start, int end, int methyPred, int totalPred, int methyObs, int totalObs) {
			this.chr = chr;
			this.start = start;
			this.end = end;
			this.methyPred = methyPred;
			this.totalPred = totalPred;
			this.methyObs = methyObs;
			this.totalObs = totalObs;
		}

		double predictMethyPercent() {
			return 100.0 * (double) methyPred / (double) totalPred;
		}

		double observedMethyPercent() {
			if (totalObs <= 0) {
				return Double.NaN;
			}
			return 100.0 * (double) methyObs / (double) totalObs;
		}
	}

	private static class CpgIndex {
		final LinkedHashMap<String, int[]> chrPositions;
		final HashMap<String, Integer> chrOffsets;
		final int totalSites;

		CpgIndex(LinkedHashMap<String, int[]> chrPositions, HashMap<String, Integer> chrOffsets, int totalSites) {
			this.chrPositions = chrPositions;
			this.chrOffsets = chrOffsets;
			this.totalSites = totalSites;
		}

		int getGlobalIndex(String chr, int start) {
			int[] positions = chrPositions.get(chr);
			Integer offset = chrOffsets.get(chr);
			if (positions == null || offset == null) {
				return -1;
			}
			// Try exact match first (same coordinate system)
			int localIndex = Arrays.binarySearch(positions, start);
			if (localIndex < 0) {
				// Try start+1: FinaleMe uses 0-based BED start (e.g. 10468),
				// wgbstools CpG.bed.gz uses 0-based CpG position (e.g. 10469, the C in CG)
				localIndex = Arrays.binarySearch(positions, start + 1);
			}
			if (localIndex < 0) {
				return -1;
			}
			return offset + localIndex + 1;
		}
	}

	private static class ParsedRow {
		final String readName;
		final String methyStat;
		final String loc;
		final int offset;
		final double methyPrior;
		final double fragLen;
		final double coverage;
		final double distToCenter;

		ParsedRow(String readName, String methyStat, String loc, int offset,
				  double methyPrior, double fragLen, double coverage, double distToCenter) {
			this.readName = readName;
			this.methyStat = methyStat;
			this.loc = loc;
			this.offset = offset;
			this.methyPrior = methyPrior;
			this.fragLen = fragLen;
			this.coverage = coverage;
			this.distToCenter = distToCenter;
		}
	}

	private static class AssembledFragment {
		final HashMap<Integer, Pair<Integer, Double>> cpgDistRow;
		final ArrayList<ObservationVector> matrixRow;
		final ArrayList<String> locRow;
		final ArrayList<Integer> observedRow;

		AssembledFragment(HashMap<Integer, Pair<Integer, Double>> cpgDistRow,
						  ArrayList<ObservationVector> matrixRow,
						  ArrayList<String> locRow,
						  ArrayList<Integer> observedRow) {
			this.cpgDistRow = cpgDistRow;
			this.matrixRow = matrixRow;
			this.locRow = locRow;
			this.observedRow = observedRow;
		}
	}

}
