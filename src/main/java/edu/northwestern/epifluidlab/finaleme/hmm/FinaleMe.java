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

	
	@Option(name="-useEndMotif",usage="include 5' end motif score as a third feature in lowCoverage mode (3D: fragLen, distToCenter, motifScore). Default: false")
	public boolean useEndMotif = false;

	@Option(name="-adaptEmissionOnly",usage="constrained Baum-Welch: freeze transitions/initiation, adapt emissions only. Requires -decodeModeOnly. Default: false")
	public boolean adaptEmissionOnly = false;

	@Option(name="-adaptLambda",usage="shrinkage regularization toward reference model (0=no regularization, 1=no adaptation). Default: 0.5")
	public double adaptLambda = 0.5;

	@Option(name="-adaptMaxIter",usage="max Baum-Welch iterations during emission adaptation. Default: 5")
	public int adaptMaxIter = 5;

	@Option(name="-adaptMinFragments",usage="minimum fragments with >= miniDataPoints CpGs to attempt adaptation; below this, use reference model directly. Default: 1000")
	public int adaptMinFragments = 1000;

	@Option(name="-saveAdaptedModel",usage="save the adapted model to this path. If not set, the adapted model is used for decoding but not saved, and the reference model is not modified.")
	public String saveAdaptedModel = null;

	@Option(name="-saveNormStats",usage="save per-feature normalization statistics (mean/sd) to TSV file during training. Required for proper adaptation re-centering.")
	public String saveNormStats = null;

	@Option(name="-loadNormStats",usage="load reference normalization statistics from TSV file. Used by -adaptEmissionOnly to re-center reference GMM into the target sample's z-score space.")
	public String loadNormStats = null;

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "FinaleMe [opts] model input_matrix.txt[.gz] prediction.txt.gz";
	
	private static final Logger log = LoggerFactory.getLogger(FinaleMe.class);

	// Last computed feature normalization stats from processMatrixFile/collectStats.
	// Used by adaptAndDecodeStreaming to access target sample stats for re-centering.
	private SummaryStatistics[] lastComputedStats = null;

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
			gzipInputStream = new GZIPInputStream(new FileInputStream(path), GZIP_BUFFER_SIZE);
			br = new BufferedReader(new InputStreamReader(gzipInputStream), READER_BUFFER_SIZE);
		} else {
			br = new BufferedReader(new FileReader(path), READER_BUFFER_SIZE);
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
		int minCols = features + 4 + (useEndMotif ? 1 : 0);
		if (splitLines.length < minCols || splitLines[1].equalsIgnoreCase("start")
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

		// When -useEndMotif is set, motif_score is at col [11] and methyPrior shifts to [12]
		double motifScore = Double.NaN;
		int methyPriorCol = 11;
		if (useEndMotif) {
			motifScore = Double.parseDouble(splitLines[11]);
			methyPriorCol = 12;
		}

		double methyPrior = Double.parseDouble(splitLines[methyPriorCol]);
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
		return new ParsedRow(readName, splitLines[6], loc, offset, methyPrior, fragLen, coverage, distToCenter, motifScore);
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

	// 64KB buffer for GZIPInputStream (default is 512 bytes — 128x larger reduces
	// inflate syscall overhead dramatically on large files like 356M-row matrices).
	private static final int GZIP_BUFFER_SIZE = 65536;
	private static final int READER_BUFFER_SIZE = 131072;

	/**
	 * Detect whether a gzipped file is actually BGZF (bgzipped). Bgzip files
	 * have a specific extra field 'BC' in the gzip header at bytes 14-15.
	 * Returns false if not bgzipped or on I/O error.
	 */
	private static boolean isBgzippedFile(String path) {
		try (FileInputStream fis = new FileInputStream(path)) {
			byte[] hdr = new byte[18];
			int n = fis.read(hdr);
			if (n < 18) return false;
			// gzip magic: 1f 8b, FLG.FEXTRA = 0x04
			if ((hdr[0] & 0xff) != 0x1f || (hdr[1] & 0xff) != 0x8b) return false;
			if ((hdr[3] & 0x04) == 0) return false;
			// BGZF subfield: SI1='B'=0x42, SI2='C'=0x43 at offset 12-13
			return hdr[12] == 0x42 && hdr[13] == 0x43;
		} catch (IOException e) {
			return false;
		}
	}

	/**
	 * Open a gzipped file for reading, using parallel bgzip decompression
	 * when the file is bgzipped and threads > 1 and `bgzip` is on PATH.
	 * Falls back to single-threaded GZIPInputStream otherwise.
	 *
	 * The returned BufferedReader must be closed by the caller. If a
	 * subprocess was spawned, closing the reader also closes its stream;
	 * the subprocess terminates when stdin EOF is reached.
	 */
	private BufferedReader openGzipReader(String path) throws IOException {
		if (!path.endsWith(".gz")) {
			return new BufferedReader(new FileReader(path), READER_BUFFER_SIZE);
		}
		int nThreads = resolveThreadCount();
		if (nThreads > 1 && isBgzippedFile(path)) {
			try {
				ProcessBuilder pb = new ProcessBuilder(
					"bgzip", "-d", "-c", "-@", String.valueOf(nThreads), path);
				pb.redirectErrorStream(false);
				Process proc = pb.start();
				log.info("Reading bgzipped file with parallel decompression: bgzip -d -@ " + nThreads);
				return new BufferedReader(new InputStreamReader(proc.getInputStream()), READER_BUFFER_SIZE);
			} catch (IOException e) {
				log.info("bgzip not available on PATH, falling back to single-threaded gzip read: " + e.getMessage());
			}
		}
		return new BufferedReader(
			new InputStreamReader(new GZIPInputStream(new FileInputStream(path), GZIP_BUFFER_SIZE)),
			READER_BUFFER_SIZE);
	}

	private SummaryStatistics[] collectStats(String matrixFile,
											 HashMap<String, IntervalTree<Integer>> overlapLoc,
											 HashMap<String, IntervalTree<Integer>> excludeLoc) throws IOException {
		int numStats = (lowCoverage && useEndMotif) ? 4 : 3;
		SummaryStatistics[] stats = new SummaryStatistics[numStats];
		for (int i = 0; i < numStats; i++) {
			stats[i] = new SummaryStatistics();
		}

		BufferedReader br = openGzipReader(matrixFile);

		String line;
		while ((line = br.readLine()) != null) {
			ParsedRow row = parseLine(line, overlapLoc, excludeLoc);
			if (row == null) continue;
			stats[0].addValue(row.fragLen);
			stats[1].addValue(row.coverage);
			stats[2].addValue(row.distToCenter);
			if (lowCoverage && useEndMotif && !Double.isNaN(row.motifScore)) {
				stats[3].addValue(row.motifScore);
			}
		}
		br.close();

		return stats;
	}

	private void logFeatureStats(SummaryStatistics[] stats) {
		final String[] featureNames = new String[]{"FragLen", "Norm_Frag_cov", "DistToCenter", "MotifScore"};
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

					if (adaptEmissionOnly && decodeModeOnly) {
						// Adaptation + decode path
						CpgIndex cpgIndex = null;
						if (patOutput) {
							cpgIndex = loadCpgIndex(cpgIndexFile);
						}
						LinkedHashMap<String, Integer> chromOrder = null;
						if (bwOutput) {
							chromOrder = loadChromOrder(chromSizeFile);
						}
						adaptAndDecodeStreaming(inputFile, modelFile, outputFile, cpgIndex, chromOrder);
					} else if (decodeModeOnly && !aucMode) {
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
				gzipInputStream = new GZIPInputStream(new FileInputStream(region), GZIP_BUFFER_SIZE);
				br = new BufferedReader(new InputStreamReader(gzipInputStream), READER_BUFFER_SIZE);
				
			}else{
				br = new BufferedReader(new FileReader(region), READER_BUFFER_SIZE);
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
				gzipInputStream = new GZIPInputStream(new FileInputStream(exclude), GZIP_BUFFER_SIZE);
				br = new BufferedReader(new InputStreamReader(gzipInputStream), READER_BUFFER_SIZE);
				
			}else{
				br = new BufferedReader(new FileReader(exclude), READER_BUFFER_SIZE);
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
		BufferedReader br = openGzipReader(matrixFile);

			// Two-pass file-based approach to avoid storing intermediate data.
			// Pass 1: collect feature statistics AND count CpGs per readName.
			//          The per-read CpG count lets Pass 2 skip reads that won't
			//          meet miniDataPoints — typically ~90% of all rows, saving
			//          ~60GB of matrixProcess memory.
			// Pass 2: re-read file, skip reads below threshold, z-score normalize,
			//          build matrixProcess directly.

			// === Pass 1: Stats + per-read CpG count ===
			String line;
			int numStats = (lowCoverage && useEndMotif) ? 4 : 3;
			SummaryStatistics[] stats = new SummaryStatistics[numStats];
			for(int i = 0; i < numStats; i++){
				stats[i] = new SummaryStatistics();
			}
			long statsRows = 0;

			// Count CpGs per readName so Pass 2 can skip reads with too few CpGs.
			// At 30X WGS: ~30M unique reads × ~100 bytes/entry ≈ 3-5GB.
			HashMap<String, int[]> readCpgCount = new HashMap<>();

			while( (line = br.readLine()) != null){
				ParsedRow row = parseLine(line, overlapLoc, excludeLoc);
				if(row == null) continue;
				stats[0].addValue(row.fragLen);
				stats[1].addValue(row.coverage);
				stats[2].addValue(row.distToCenter);
				if(lowCoverage && useEndMotif && !Double.isNaN(row.motifScore)){
					stats[3].addValue(row.motifScore);
				}
				statsRows++;
				// Count CpGs per read (use int[1] to avoid Integer boxing)
				int[] cnt = readCpgCount.get(row.readName);
				if(cnt == null){
					readCpgCount.put(row.readName, new int[]{1});
				}else{
					cnt[0]++;
				}
			}
			br.close();

		logFeatureStats(stats);
		lastComputedStats = stats; // expose for adaptAndDecodeStreaming

		// Save normalization stats for future adaptation re-centering
		if (saveNormStats != null) {
			saveNormalizationStats(stats, saveNormStats);
		}

		long totalReads = readCpgCount.size();
		long qualifyingReads = 0;
		for(int[] cnt : readCpgCount.values()){
			if(cnt[0] >= miniDataPoints) qualifyingReads++;
		}
		log.info("Pass 1 done: " + statsRows + " rows, " + totalReads +
				 " unique reads, " + qualifyingReads + " with >= " + miniDataPoints + " CpGs");

		// === Pass 2: Re-read file, skip non-qualifying reads, build matrixProcess ===
		log.info("Pass 2: Building observation vectors (skipping " +
				 (totalReads - qualifyingReads) + " reads with < " + miniDataPoints + " CpGs) ...");
		BufferedReader br2 = openGzipReader(matrixFile);
		long skippedRows = 0;

		while( (line = br2.readLine()) != null){
			ParsedRow row = parseLine(line, overlapLoc, excludeLoc);
			if(row == null) continue;

			// Early skip: if this read has fewer CpGs than miniDataPoints,
			// it will be discarded by assembleFragment anyway. Skip it now
			// to avoid allocating ObservationVector/Triple/Pair/loc String.
			int[] cnt = readCpgCount.get(row.readName);
			if(cnt == null || cnt[0] < miniDataPoints){
				skippedRows++;
				continue;
			}

			if(covOutlier > 0 && ((row.coverage-stats[1].getMean())/stats[1].getStandardDeviation() > covOutlier ||
					(row.fragLen-stats[0].getMean())/stats[0].getStandardDeviation() > covOutlier ||
					(row.distToCenter-stats[2].getMean())/stats[2].getStandardDeviation() > covOutlier)){
				continue;
			}
			double[] value;

			if(lowCoverage && useEndMotif){
				value = new double[]{
						(row.fragLen-stats[0].getMean())/stats[0].getStandardDeviation(),
						(row.distToCenter-stats[2].getMean())/stats[2].getStandardDeviation(),
						(row.motifScore-stats[3].getMean())/stats[3].getStandardDeviation(),
				};
			}else if(lowCoverage){
				value = new double[]{
						(row.fragLen-stats[0].getMean())/stats[0].getStandardDeviation(),
						(row.distToCenter-stats[2].getMean())/stats[2].getStandardDeviation(),
				};
			}else{
				value = new double[]{
						(row.fragLen-stats[0].getMean())/stats[0].getStandardDeviation(),
						(row.coverage-stats[1].getMean())/stats[1].getStandardDeviation(),
						(row.distToCenter-stats[2].getMean())/stats[2].getStandardDeviation(),
				};
			}

			ObservationVector vector = new ObservationVector(value);

			points++;
			if(matrixProcess.containsKey(row.readName)){
				TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>> readStat = matrixProcess.get(row.readName);
				if(!readStat.containsKey(row.offset)){
					readStat.put(row.offset, Triple.of(row.methyStat, vector, new Pair<String, Double>(row.loc, row.methyPrior)));
				}
				matrixProcess.put(row.readName, readStat);
			}else{
				TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>> readStat =  new TreeMap<Integer, Triple<String, ObservationVector, Pair<String, Double>>>();
				readStat.put(row.offset, Triple.of(row.methyStat, vector, new Pair<String, Double>(row.loc, row.methyPrior)));
				matrixProcess.put(row.readName, readStat);
			}
		}
		br2.close();
		readCpgCount = null; // free the count map
		log.info("Pass 2 done: " + points + " points loaded, " + skippedRows + " rows skipped");
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

		// Open input for second pass (uses parallel bgzip decompression if applicable)
		BufferedReader br = openGzipReader(inputFile);

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
			if (lowCoverage && useEndMotif) {
				value = new double[]{
					(row.fragLen - stats[0].getMean()) / stats[0].getStandardDeviation(),
					(row.distToCenter - stats[2].getMean()) / stats[2].getStandardDeviation(),
					(row.motifScore - stats[3].getMean()) / stats[3].getStandardDeviation(),
				};
			} else if (lowCoverage) {
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
				input = new GZIPInputStream(new FileInputStream(cpgIndexFile), GZIP_BUFFER_SIZE);
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
				input = new GZIPInputStream(new FileInputStream(chromSizeFile), GZIP_BUFFER_SIZE);
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
				input = new GZIPInputStream(new FileInputStream(chromSizeFile), GZIP_BUFFER_SIZE);
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
		try (BlockCompressedOutputStream output = new BlockCompressedOutputStream(bgzfOutput)) {
			OutputStreamWriter writer = new OutputStreamWriter(output, StandardCharsets.UTF_8);
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
			// Flush writer bytes into BGZF stream. We intentionally do not close
			// writer here because closing it would close `output`, which is already
			// managed by the outer try-with-resources.
			writer.flush();
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
			features = useEndMotif ? 3 : 2;
		}

		if(mixNumberInFeature == null || mixNumberInFeature.isEmpty()){
			mixNumberInFeature = new ArrayList<Integer>();
			for(int i = 0; i < features; i++){
				mixNumberInFeature.add(1);
			}
		}else if(mixNumberInFeature.size() != features){
			throw new IllegalArgumentException("Wrong number of mixNumberInFeature" + mixNumberInFeature);
		}

		// Validate adaptation flags
		if(adaptEmissionOnly && !decodeModeOnly){
			throw new IllegalArgumentException("-adaptEmissionOnly requires -decodeModeOnly");
		}
		if(saveAdaptedModel != null && !adaptEmissionOnly){
			throw new IllegalArgumentException("-saveAdaptedModel requires -adaptEmissionOnly");
		}
		if(saveAdaptedModel != null && !decodeModeOnly){
			throw new IllegalArgumentException("-saveAdaptedModel requires -decodeModeOnly");
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
		final double motifScore;

		ParsedRow(String readName, String methyStat, String loc, int offset,
				  double methyPrior, double fragLen, double coverage, double distToCenter,
				  double motifScore) {
			this.readName = readName;
			this.methyStat = methyStat;
			this.loc = loc;
			this.offset = offset;
			this.methyPrior = methyPrior;
			this.fragLen = fragLen;
			this.coverage = coverage;
			this.distToCenter = distToCenter;
			this.motifScore = motifScore;
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

	// ======================== Emission Adaptation Methods ========================

	/**
	 * Regularize MLE emission parameters toward reference model parameters.
	 * new = (1 - lambda) * mle + lambda * ref, applied to mu, sigma, proportions.
	 */
	@SuppressWarnings("unchecked")
	private OpdfMultiMixtureGaussian regularizeGmm(OpdfMultiMixtureGaussian mle,
												   OpdfMultiMixtureGaussian ref,
												   double lambda) {
		MultiMixtureGaussianDistribution mleDist = mle.getDistribution();
		MultiMixtureGaussianDistribution refDist = ref.getDistribution();

		int dim = mleDist.getDimension();
		ArrayList<Integer> mixNumbers = mleDist.getMixNumberInFeature();

		// Interpolate per-component parameters
		ArrayList<ArrayList<Double>> newMeans = new ArrayList<>();
		ArrayList<ArrayList<Double>> newVars = new ArrayList<>();
		ArrayList<ArrayList<Double>> newProps = new ArrayList<>();

		for (int d = 0; d < dim; d++) {
			ArrayList<Double> mleMeansD = mleDist.getMeanInEachGaussian().get(d);
			ArrayList<Double> refMeansD = refDist.getMeanInEachGaussian().get(d);
			ArrayList<Double> mleVarsD = mleDist.getVarianceInEachGaussian().get(d);
			ArrayList<Double> refVarsD = refDist.getVarianceInEachGaussian().get(d);
			ArrayList<Double> mlePropsD = mleDist.getPropInEachGaussian().get(d);
			ArrayList<Double> refPropsD = refDist.getPropInEachGaussian().get(d);

			int nComponents = mixNumbers.get(d);
			ArrayList<Double> interpMeans = new ArrayList<>();
			ArrayList<Double> interpVars = new ArrayList<>();
			ArrayList<Double> interpProps = new ArrayList<>();

			for (int k = 0; k < nComponents; k++) {
				interpMeans.add((1 - lambda) * mleMeansD.get(k) + lambda * refMeansD.get(k));
				interpVars.add((1 - lambda) * mleVarsD.get(k) + lambda * refVarsD.get(k));
				interpProps.add((1 - lambda) * mlePropsD.get(k) + lambda * refPropsD.get(k));
			}
			newMeans.add(interpMeans);
			newVars.add(interpVars);
			newProps.add(interpProps);
		}

		// Interpolate top-level mean and covariance
		double[] mleMean = mle.mean();
		double[] refMean = ref.mean();
		double[][] mleCov = mle.covariance();
		double[][] refCov = ref.covariance();

		double[] newMean = new double[dim];
		double[][] newCov = new double[dim][dim];
		for (int i = 0; i < dim; i++) {
			newMean[i] = (1 - lambda) * mleMean[i] + lambda * refMean[i];
			for (int j = 0; j < dim; j++) {
				newCov[i][j] = (1 - lambda) * mleCov[i][j] + lambda * refCov[i][j];
			}
		}

		return new OpdfMultiMixtureGaussian(newMean, newCov, mixNumbers, newMeans, newVars, newProps);
	}

	/**
	 * Perform one constrained Baum-Welch iteration:
	 * 1. Run standard BW iterate (updates everything)
	 * 2. Restore frozen transitions and pi from reference
	 * 3. Regularize emissions toward reference
	 */
	private BayesianNhmmV5<ObservationVector> adaptIteration(
			BayesianNhmmV5<ObservationVector> currentHmm,
			BayesianNhmmV5<ObservationVector> refHmm,
			List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> sequences,
			double lambda) {

		int nThreads = resolveThreadCount();
		BaumWelchBayesianNhmmV5ScaledLearner bwl = new BaumWelchBayesianNhmmV5ScaledLearner(nThreads);

		// Full BW iteration (updates everything: transitions, pi, emissions)
		BayesianNhmmV5<ObservationVector> updatedHmm = bwl.iterate(currentHmm, sequences);

		// Restore frozen transitions and pi from reference model
		for (int r = 0; r <= refHmm.nbCpgDistState(); r++) {
			for (int i = 0; i < refHmm.nbStates(); i++) {
				updatedHmm.setPri(r, i, refHmm.getPri(r, i));
				for (int j = 0; j < refHmm.nbStates(); j++) {
					updatedHmm.setArij(r, i, j, refHmm.getArij(r, i, j));
				}
			}
		}

		// Regularize emissions toward reference
		for (int s = 0; s < updatedHmm.nbStates(); s++) {
			OpdfMultiMixtureGaussian mleOpdf = (OpdfMultiMixtureGaussian) updatedHmm.getOpdf(s);
			OpdfMultiMixtureGaussian refOpdf = (OpdfMultiMixtureGaussian) refHmm.getOpdf(s);
			OpdfMultiMixtureGaussian regularized = regularizeGmm(mleOpdf, refOpdf, lambda);
			updatedHmm.setOpdf(s, regularized);
		}

		return updatedHmm;
	}

	/**
	 * Adaptation + decode streaming mode:
	 * 1. Load reference model
	 * 2. Collect stats, load data
	 * 3. Constrained BW adaptation (emissions only)
	 * 4. Viterbi decode all fragments using adapted model
	 *
	 * Note: We do not re-center reference GMM means to target z-score space.
	 * Since both reference and target are independently z-score normalized,
	 * adaptation iterations correct any distribution shift. This can be
	 * revisited if adaptation proves insufficient for large distribution shifts.
	 */
	private void adaptAndDecodeStreaming(String inputFile, String modelFile, String outputFile,
										 CpgIndex cpgIndex, LinkedHashMap<String, Integer> chromOrder) throws Exception {
		System.out.println("\nAdaptation + decode mode ...\n");

		// Load reference model
		BayesianNhmmV5<ObservationVector> refHmm = loadHmmModel(modelFile);
		refHmm.setBayesianFactor(bayesianFactor);
		refHmm.setMethyState(this.methylatedState);
		refHmm.setMaxCpgNum(cpgNumClip < 0 ? maxCpgs : cpgNumClip);
		refHmm.setMinCpgNum(1);

		// Load data for adaptation via processMatrixFile (handles its own
		// two-pass stats + CpG-count pre-filter internally; no separate
		// collectStats call needed).
		// Use the user's -miniDataPoints (e.g. 7) for training fragment selection.
		// Baum-Welch requires >= 2 CpGs per fragment; enforce as a floor.
		int adaptMiniDataPoints = Math.max(miniDataPoints, 2);
		int savedMiniDataPoints = miniDataPoints;
		miniDataPoints = adaptMiniDataPoints;
		log.info("Loading data for adaptation (miniDataPoints=" + adaptMiniDataPoints + ") ...");
		MatrixObj matrixObj = processMatrixFile(inputFile);
		miniDataPoints = savedMiniDataPoints; // restore original
		matrixObj.cpgDistFreq = null;

		List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> matrix =
			new ArrayList<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>>();
		for (org.apache.commons.lang3.tuple.Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>> row : matrixObj.matrix) {
			if (row.getMiddle().size() >= 2) {
				matrix.add(new Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>(row.getLeft(), row.getMiddle()));
			}
		}

		log.info("Loaded " + matrix.size() + " fragments for adaptation (>= " + adaptMiniDataPoints + " CpGs each)");

		if (matrix.size() < adaptMinFragments) {
			log.warn("Only " + matrix.size() + " qualifying fragments (< " + adaptMinFragments +
					 "); skipping adaptation, decoding with reference model directly.");
			matrixObj = null;
			matrix = null;
			miniDataPoints = 1;
			decodeOnlyStreaming(inputFile, modelFile, outputFile, cpgIndex, chromOrder);
			miniDataPoints = savedMiniDataPoints;
			return;
		}

		// Phase 3: Constrained Baum-Welch adaptation
		log.info("Phase 3: Constrained Baum-Welch emission adaptation (lambda=" + adaptLambda +
				 ", maxIter=" + adaptMaxIter + ") ...");

		logEmissionParams("REFERENCE MODEL (before adaptation)", refHmm);

		// Re-center reference GMM means from reference z-score space into
		// the target sample's z-score space (design doc §3.3.3):
		//   μ_adjusted = (μ_ref × σ_ref + mean_ref − mean_target) / σ_target
		// This corrects for the fact that z-score normalization uses each
		// sample's own mean/sd, so the same raw value gets different z-scores
		// in different samples.
		BayesianNhmmV5<ObservationVector> refHmmForAdapt;
		try {
			refHmmForAdapt = refHmm.clone();
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException("Failed to clone reference HMM", e);
		}
		if (loadNormStats != null) {
			double[][] refNormStats = loadNormalizationStats(loadNormStats);
			SummaryStatistics[] targetStats = lastComputedStats; // from processMatrixFile above
			if (targetStats != null && refNormStats != null) {
				log.info("Re-centering reference GMM from ref z-score space to target z-score space ...");
				recenterEmissions(refHmmForAdapt, refNormStats, targetStats);
				logEmissionParams("REFERENCE MODEL (after re-centering)", refHmmForAdapt);
			} else {
				log.warn("Cannot re-center: " +
						(targetStats == null ? "target stats not available" : "ref norm stats not loaded"));
			}
		} else {
			log.info("No -loadNormStats; skipping re-centering (ref and target z-score spaces may differ)");
		}

		BayesianNhmmV5<ObservationVector> adaptedHmm;
		try {
			adaptedHmm = refHmmForAdapt.clone();
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException("Failed to clone reference HMM", e);
		}

		KullbackLeiblerDistanceBayesianNhmmV5Calculator klc =
			new KullbackLeiblerDistanceBayesianNhmmV5Calculator(matrix, resolveThreadCount());

		double distance = Double.MAX_VALUE;
		for (int iter = 0; iter < adaptMaxIter; iter++) {
			BayesianNhmmV5<ObservationVector> prevHmm;
			try {
				prevHmm = adaptedHmm.clone();
			} catch (CloneNotSupportedException e) {
				throw new RuntimeException("Failed to clone HMM at iteration " + iter, e);
			}

			adaptedHmm = adaptIteration(adaptedHmm, refHmmForAdapt, matrix, adaptLambda);

			distance = klc.distance(prevHmm, adaptedHmm, true);
			log.info("Adaptation iteration " + (iter + 1) + ": KL distance = " + distance);
			logEmissionParams("After iteration " + (iter + 1), adaptedHmm);

			if (Double.isNaN(distance)) {
				log.warn("KL distance is NaN at iteration " + (iter + 1) +
						 "; falling back to reference model.");
				adaptedHmm = refHmmForAdapt;
				break;
			}
			if (Math.abs(distance) < tol) {
				log.info("Adaptation converged at iteration " + (iter + 1));
				break;
			}
		}

		logEmissionDelta("EMISSION DELTA (adapted - re-centered reference)", refHmmForAdapt, adaptedHmm);

		this.methylatedState = adaptedHmm.getMethyState(lowCoverage);
		adaptedHmm.setMaxCpgNum(cpgNumClip < 0 ? maxCpgs : cpgNumClip);
		adaptedHmm.setMinCpgNum(1);
		System.out.println("Adapted HMM:\n" + adaptedHmm);

		// Write adapted model to a temp file for decodeHmm (which reloads from disk).
		// The original reference model file is NEVER modified.
		File adaptedModelTmp = File.createTempFile("adapted_model_", ".ser");
		adaptedModelTmp.deleteOnExit();
		ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(adaptedModelTmp));
		oos.writeObject(adaptedHmm);
		oos.close();
		log.info("Adapted model written to temp file: " + adaptedModelTmp.getAbsolutePath());

		// Optionally save adapted model to a user-specified path
		if (saveAdaptedModel != null) {
			java.nio.file.Files.copy(adaptedModelTmp.toPath(),
				new File(saveAdaptedModel).toPath(),
				java.nio.file.StandardCopyOption.REPLACE_EXISTING);
			log.info("Adapted model saved to: " + saveAdaptedModel);
		}

		// Free adaptation data before decode to reduce memory
		matrix = null;
		matrixObj = null;

		// Phase 4: Viterbi decode ALL fragments (miniDataPoints=1) using the
		// adapted model via streaming decode. This outputs every fragment
		// regardless of the -miniDataPoints used for adaptation training.
		log.info("Phase 4: Streaming decode with adapted model (all fragments) ...");
		int decodeMiniDataPoints = miniDataPoints;
		miniDataPoints = 1;  // decode everything
		decodeOnlyStreaming(inputFile, adaptedModelTmp.getAbsolutePath(), outputFile, cpgIndex, chromOrder);
		miniDataPoints = decodeMiniDataPoints;  // restore
		adaptedModelTmp.delete();
	}

	/**
	 * Log per-state emission GMM parameters (mean, variance, mixture proportions)
	 * for adaptation diagnostics.
	 */
	private void logEmissionParams(String label, BayesianNhmmV5<ObservationVector> hmm) {
		final String[] featureNames = lowCoverage && useEndMotif
			? new String[]{"FragLen", "DistToCenter", "MotifScore"}
			: (lowCoverage ? new String[]{"FragLen", "DistToCenter"}
			: new String[]{"FragLen", "Coverage", "DistToCenter"});

		log.info("=== " + label + " ===");
		for (int s = 0; s < hmm.nbStates(); s++) {
			OpdfMultiMixtureGaussian opdf = (OpdfMultiMixtureGaussian) hmm.getOpdf(s);
			MultiMixtureGaussianDistribution dist = opdf.getDistribution();
			double[] mean = opdf.mean();
			double[][] cov = opdf.covariance();
			StringBuilder sb = new StringBuilder();
			sb.append("  State ").append(s).append(": mean=[");
			for (int d = 0; d < mean.length; d++) {
				if (d > 0) sb.append(", ");
				sb.append(String.format("%.6f", mean[d]));
			}
			sb.append("] var=[");
			for (int d = 0; d < cov.length; d++) {
				if (d > 0) sb.append(", ");
				sb.append(String.format("%.6f", cov[d][d]));
			}
			sb.append("]");
			log.info(sb.toString());

			// Per-feature per-component details
			for (int d = 0; d < dist.getDimension(); d++) {
				String fname = d < featureNames.length ? featureNames[d] : ("Feat" + d);
				ArrayList<Double> means = dist.getMeanInEachGaussian().get(d);
				ArrayList<Double> vars = dist.getVarianceInEachGaussian().get(d);
				ArrayList<Double> props = dist.getPropInEachGaussian().get(d);
				StringBuilder csb = new StringBuilder();
				csb.append("    ").append(fname).append(": ");
				for (int k = 0; k < means.size(); k++) {
					if (k > 0) csb.append(" | ");
					csb.append(String.format("mu=%.6f var=%.6f w=%.4f",
						means.get(k), vars.get(k), props.get(k)));
				}
				log.info(csb.toString());
			}
		}
	}

	/**
	 * Log the difference between reference and adapted emission parameters.
	 */
	/**
	 * Save per-feature normalization statistics (mean, sd) to a TSV file.
	 * Format: feature_index \t mean \t sd
	 */
	private void saveNormalizationStats(SummaryStatistics[] stats, String path) throws IOException {
		try (OutputStreamWriter w = new OutputStreamWriter(new FileOutputStream(path), StandardCharsets.UTF_8)) {
			w.write("feature\tmean\tsd\n");
			for (int i = 0; i < stats.length; i++) {
				w.write(i + "\t" + stats[i].getMean() + "\t" +
						(stats[i].getN() > 1 ? stats[i].getStandardDeviation() : 1.0) + "\n");
			}
		}
		log.info("Normalization stats saved to " + path + " (" + stats.length + " features)");
	}

	/**
	 * Load per-feature normalization statistics from a TSV file.
	 * Returns double[nFeatures][2] where [i][0] = mean, [i][1] = sd.
	 */
	private double[][] loadNormalizationStats(String path) throws IOException {
		ArrayList<double[]> rows = new ArrayList<>();
		try (BufferedReader br = new BufferedReader(new FileReader(path), READER_BUFFER_SIZE)) {
			String line;
			boolean header = true;
			while ((line = br.readLine()) != null) {
				if (header) { header = false; continue; }
				String[] parts = line.split("\t");
				if (parts.length >= 3) {
					rows.add(new double[]{Double.parseDouble(parts[1]), Double.parseDouble(parts[2])});
				}
			}
		}
		log.info("Loaded normalization stats from " + path + " (" + rows.size() + " features)");
		return rows.toArray(new double[0][]);
	}

	/**
	 * Re-center GMM emission means from reference z-score space to target z-score space.
	 *
	 * For each feature dimension d:
	 *   μ_adjusted = (μ_ref × σ_ref + mean_ref − mean_target) / σ_target
	 *   σ²_adjusted = σ²_ref × (σ_ref / σ_target)²
	 *
	 * This corrects for the fact that each sample's features are z-score
	 * normalized using its own mean/sd, so the reference GMM parameters
	 * are in a different coordinate system than the target data.
	 */
	private void recenterEmissions(BayesianNhmmV5<ObservationVector> hmm,
								   double[][] refNormStats,
								   SummaryStatistics[] targetStats) {
		for (int s = 0; s < hmm.nbStates(); s++) {
			OpdfMultiMixtureGaussian opdf = (OpdfMultiMixtureGaussian) hmm.getOpdf(s);
			MultiMixtureGaussianDistribution dist = opdf.getDistribution();
			int dim = dist.getDimension();

			// Build re-centered parameters
			double[] newMean = opdf.mean().clone();
			double[][] newCov = opdf.covariance().clone();
			for (int i = 0; i < dim; i++) {
				newCov[i] = newCov[i].clone();
			}
			ArrayList<ArrayList<Double>> newMeans = new ArrayList<>();
			ArrayList<ArrayList<Double>> newVars = new ArrayList<>();
			ArrayList<ArrayList<Double>> newProps = new ArrayList<>();

			for (int d = 0; d < dim; d++) {
				// Map feature dimension to stats index:
				// lowCoverage+endMotif: value[0]=fragLen(stats[0]), value[1]=distToCenter(stats[2]), value[2]=motifScore(stats[3])
				// lowCoverage: value[0]=fragLen(stats[0]), value[1]=distToCenter(stats[2])
				// normal: value[0]=fragLen(stats[0]), value[1]=coverage(stats[1]), value[2]=distToCenter(stats[2])
				int statsIdx;
				if (lowCoverage && useEndMotif) {
					statsIdx = (d == 0) ? 0 : (d == 1) ? 2 : 3;
				} else if (lowCoverage) {
					statsIdx = (d == 0) ? 0 : 2;
				} else {
					statsIdx = d; // 0→fragLen, 1→coverage, 2→distToCenter
				}

				double meanRef = refNormStats[statsIdx][0];
				double sdRef = refNormStats[statsIdx][1];
				double meanTarget = targetStats[statsIdx].getMean();
				double sdTarget = targetStats[statsIdx].getN() > 1 ? targetStats[statsIdx].getStandardDeviation() : 1.0;
				if (sdTarget == 0) sdTarget = 1.0;
				double scale = sdRef / sdTarget;

				// Re-center top-level mean and variance
				newMean[d] = (newMean[d] * sdRef + meanRef - meanTarget) / sdTarget;
				newCov[d][d] = newCov[d][d] * scale * scale;

				// Re-center per-component means and variances
				ArrayList<Double> compMeans = new ArrayList<>(dist.getMeanInEachGaussian().get(d));
				ArrayList<Double> compVars = new ArrayList<>(dist.getVarianceInEachGaussian().get(d));
				ArrayList<Double> compProps = new ArrayList<>(dist.getPropInEachGaussian().get(d));
				for (int k = 0; k < compMeans.size(); k++) {
					compMeans.set(k, (compMeans.get(k) * sdRef + meanRef - meanTarget) / sdTarget);
					compVars.set(k, compVars.get(k) * scale * scale);
				}
				newMeans.add(compMeans);
				newVars.add(compVars);
				newProps.add(compProps);
			}

			OpdfMultiMixtureGaussian recentered = new OpdfMultiMixtureGaussian(
				newMean, newCov, dist.getMixNumberInFeature(), newMeans, newVars, newProps);
			hmm.setOpdf(s, recentered);
		}
	}

	private void logEmissionDelta(String label, BayesianNhmmV5<ObservationVector> refHmm,
								  BayesianNhmmV5<ObservationVector> adaptedHmm) {
		final String[] featureNames = lowCoverage && useEndMotif
			? new String[]{"FragLen", "DistToCenter", "MotifScore"}
			: (lowCoverage ? new String[]{"FragLen", "DistToCenter"}
			: new String[]{"FragLen", "Coverage", "DistToCenter"});

		log.info("=== " + label + " ===");
		for (int s = 0; s < refHmm.nbStates(); s++) {
			OpdfMultiMixtureGaussian refOpdf = (OpdfMultiMixtureGaussian) refHmm.getOpdf(s);
			OpdfMultiMixtureGaussian adpOpdf = (OpdfMultiMixtureGaussian) adaptedHmm.getOpdf(s);
			double[] refMean = refOpdf.mean();
			double[] adpMean = adpOpdf.mean();
			double[][] refCov = refOpdf.covariance();
			double[][] adpCov = adpOpdf.covariance();

			StringBuilder sb = new StringBuilder();
			sb.append("  State ").append(s).append(": delta_mean=[");
			for (int d = 0; d < refMean.length; d++) {
				if (d > 0) sb.append(", ");
				sb.append(String.format("%+.6f", adpMean[d] - refMean[d]));
			}
			sb.append("] delta_var=[");
			for (int d = 0; d < refCov.length; d++) {
				if (d > 0) sb.append(", ");
				sb.append(String.format("%+.6f", adpCov[d][d] - refCov[d][d]));
			}
			sb.append("]");
			log.info(sb.toString());

			// Per-feature per-component deltas
			MultiMixtureGaussianDistribution refDist = refOpdf.getDistribution();
			MultiMixtureGaussianDistribution adpDist = adpOpdf.getDistribution();
			for (int d = 0; d < refDist.getDimension(); d++) {
				String fname = d < featureNames.length ? featureNames[d] : ("Feat" + d);
				ArrayList<Double> refMeans = refDist.getMeanInEachGaussian().get(d);
				ArrayList<Double> adpMeans = adpDist.getMeanInEachGaussian().get(d);
				ArrayList<Double> refVars = refDist.getVarianceInEachGaussian().get(d);
				ArrayList<Double> adpVars = adpDist.getVarianceInEachGaussian().get(d);
				ArrayList<Double> refProps = refDist.getPropInEachGaussian().get(d);
				ArrayList<Double> adpProps = adpDist.getPropInEachGaussian().get(d);
				StringBuilder csb = new StringBuilder();
				csb.append("    ").append(fname).append(": ");
				for (int k = 0; k < refMeans.size(); k++) {
					if (k > 0) csb.append(" | ");
					csb.append(String.format("Δmu=%+.6f Δvar=%+.6f Δw=%+.4f",
						adpMeans.get(k) - refMeans.get(k),
						adpVars.get(k) - refVars.get(k),
						adpProps.get(k) - refProps.get(k)));
				}
				log.info(csb.toString());
			}
		}
	}

}
