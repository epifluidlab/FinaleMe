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
import java.io.FileWriter;
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
import java.util.HashSet;
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

	@Option(name="-features",usage="dimension of the HMM observation vector. Must match the active feature flags: 3 + (-useEndMotif?1:0) + (-useBaseQ?1:0) - (-lowCoverage?1:0). Mismatch is caught at GMM-init time with a clear error, but it's better to set this correctly up front. Default: 3")
	public int features = 3;
	
	@Option(name="-tol",usage="tolerence level for the converge, default: 1e-5")
	public double tol = 1e-5;
	
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

	@Option(name="-aucMode",usage="evaluate the trained HMM against the methy_stat truth column in -aucMode -- sweeps a range of thresholds, prints (FPR, TPR) per threshold to stdout (legacy behavior), and now also computes scalar AUROC + AUPRC and (optionally) writes the full per-threshold curve to TSV / a 2-panel ROC+PR PDF. Default: false")
	public boolean aucMode = false;

	@Option(name="-aucCurveTsv",usage="write the per-threshold ROC/PR curve data (threshold, TP, FN, FP, TN, FPR, TPR, Precision, Recall) to this TSV path. AUROC and AUPRC are written as commented header lines. Implies -aucMode. Default: null")
	public String aucCurveTsv = null;

	@Option(name="-aucCurvePdf",usage="write a 2-panel ROC + Precision-Recall plot to this PDF path. Implies -aucMode and requires -aucCurveTsv (the data source) plus python3 with matplotlib on PATH (the plot is rendered by scripts/plot_roc_curve.py). Default: null")
	public String aucCurvePdf = null;

	@Option(name="-aucNThresholds",usage="number of evenly-spaced thresholds in [-1, 1] to sweep in -aucMode for the AUC integration. 0 = use the legacy 19-point non-uniform list (denser at extremes; fastest, lowest curve resolution). Higher values give smoother curves and tighter trapezoidal AUC estimates at the cost of N Viterbi passes over the input. Default: 0")
	public int aucNThresholds = 0;
	
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

	@Option(name="-useBaseQ",usage="add per-CpG base-quality score (Phred) as an additional feature dimension in the HMM observation vector. Driven by Volkov et al. 2026 (biorxiv 10.64898/2026.03.08.710357), which shows that per-base quality scores at fragment positions correlate with cancer-associated fragmentation features. Complements -useEndMotif (which aggregates baseQ-driven errors at the FRAGMENT 5' end into a per-motif score) by exposing the per-CpG-SITE quality variation along the fragment directly. Resulting feature dimension is 3 + (useEndMotif?1:0) + (useBaseQ?1:0) - (lowCoverage?1:0); set -features accordingly. Note: tabix-fragment input has constant baseQ (= -fragBaseQ); a Pass-1 degeneracy guard will throw if baseQ has SD=0. Default: false")
	public boolean useBaseQ = false;

	@Option(name="-adaptEmissionOnly",usage="constrained Baum-Welch: freeze transitions/initiation, adapt emissions only. Requires -decodeModeOnly. Default: false")
	public boolean adaptEmissionOnly = false;

	@Option(name="-adaptLambda",usage="shrinkage regularization toward reference model (0=no regularization, 1=no adaptation). Default: 0.5. Ignored when -autoAdaptLambda is set.")
	public double adaptLambda = 0.5;

	@Option(name="-autoAdaptLambda",usage="auto-tune -adaptLambda based on qualifying fragment count via a generalized (Hill) shrinkage: lambda = 1 / (1 + (N/N0)^beta). Low coverage -> more regularization toward reference. Default: false")
	public boolean autoAdaptLambda = false;

	@Option(name="-autoAdaptLambdaN0",usage="characteristic fragment count for -autoAdaptLambda: at N=N0, lambda=0.5. Default: 25000. With default beta=0.7, this gives lambda ~ 0.75 at 0.1X (5K), ~ 0.37 at 1X (53K), ~ 0.10 at 10X (530K), ~ 0.018 at 30X (6.3M).")
	public long autoAdaptLambdaN0 = 25_000L;

	@Option(name="-autoAdaptLambdaBeta",usage="curve steepness for -autoAdaptLambda: lambda = 1 / (1 + (N/N0)^beta). beta=1.0 reproduces the simple Bayesian shrinkage (linear in N). beta=0.7 (default) gives a gentle transition at low coverage (meaningful adaptation at 1X with ~53K fragments) while still strongly regularizing at 0.1X (~5K fragments, the effective noise floor for a 2-state 3-feature GMM). Larger beta -> sharper transition, smaller beta -> gentler.")
	public double autoAdaptLambdaBeta = 0.7;

	@Option(name="-autoTuneProbeLambda",usage="fixed lambda used for the 'signature probe' BW run that measures sample-vs-reference emission shift for -autoTuneBayesianFactor. Default: 0.0 (unregularized MLE -- maximum sample signature sensitivity). Decouples the shift signal from -autoAdaptLambda's decoding regularization. Set to a negative value to disable the probe and reuse the decoding-BW shift directly.")
	public double autoTuneProbeLambda = 0.0;

	@Option(name="-autoTuneProbeMaxIter",usage="max iterations for the signature probe BW. Default: 20 (needs more at lambda=0 to fully converge the MLE). Usually converges in 5-10 iter.")
	public int autoTuneProbeMaxIter = 20;

	@Option(name="-adaptMaxIter",usage="max Baum-Welch iterations during emission adaptation. Default: 50")
	public int adaptMaxIter = 50;

	@Option(name="-adaptMinFragments",usage="minimum fragments with >= miniDataPoints CpGs to attempt adaptation; below this, use reference model directly. Default: 1000")
	public int adaptMinFragments = 1000;

	@Option(name="-saveAdaptedModel",usage="save the adapted model to this path. If not set, the adapted model is used for decoding but not saved, and the reference model is not modified.")
	public String saveAdaptedModel = null;

	@Option(name="-saveNormStats",usage="save per-feature normalization statistics (mean/sd) to TSV file during training. Required for proper adaptation re-centering.")
	public String saveNormStats = null;

	@Option(name="-loadNormStats",usage="load reference normalization statistics from TSV file. Used by -adaptEmissionOnly to re-center reference GMM into the target sample's z-score space.")
	public String loadNormStats = null;

	@Option(name="-adaptReinitGmm",usage="re-initialize emission GMM on the target data before constrained Baum-Welch. Escapes the reference-model local optimum so emissions reflect the target distribution. Combine with -adaptLambda 0 for full emission retraining. Default: false")
	public boolean adaptReinitGmm = false;

	@Option(name="-adaptTransitions",usage="also adapt transitions and pi (initial state probabilities) during Baum-Welch. Without this flag, transitions/pi are frozen to the reference model. The same -adaptLambda regularization is applied. Default: false")
	public boolean adaptTransitions = false;

	@Option(name="-autoTuneBayesianFactor",usage="auto-tune bayesianFactor from sample-specific prior-vs-emission concordance. Low concordance (prior disagrees with HMM) -> low bayesianFactor; high concordance -> high bayesianFactor. Overrides -bayesianFactor after adaptation. Default: false")
	public boolean autoTuneBayesianFactor = false;

	@Option(name="-autoTuneSampleFragments",usage="[deprecated: kept for compatibility] number of random fragments to sample. No longer used. Default: 10000")
	public int autoTuneSampleFragments = 10000;

	@Option(name="-autoTuneScale",usage="[deprecated; use -autoTuneMidpoint/-autoTuneSteepness] legacy linear scale factor. Kept only for backward compatibility; ignored when -autoTuneMidpoint is non-zero. Default: 1.5")
	public double autoTuneScale = 1.5;

	@Option(name="-autoTuneMidpoint",usage="emission shift (in z-score sd) at which bayesianFactor is halfway between 0.05 and 0.9. The sigmoid transitions sharply around this value. Default: 0.60 (calibrated from observed healthy vs disease shifts using unregularized signature probe)")
	public double autoTuneMidpoint = 0.60;

	@Option(name="-autoTuneSteepness",usage="sigmoid steepness: higher = sharper transition between trust-prior and distrust-prior regimes. Default: 20.0")
	public double autoTuneSteepness = 20.0;

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

	// ----------------------------------------------------------------------
	// Distance-column compatibility check.
	//
	// CpgFeatureMatrixBuilder writes one of two distance columns:
	//   - "Dist_frag_end"     legacy unsigned distance-to-nearest-end
	//   - "Signed_Dist_Center" -useSignedDistCenter, strand-aware signed
	//                          distance from fragment center
	// The HMM reads features by COLUMN INDEX, not by name, so a model
	// trained on one form will silently mis-predict if decoded against a
	// BED of the other form. We persist the training-time column name in a
	// sidecar file next to the model (`<model>.meta.tsv`) and validate
	// matching at every decode entry point.
	// ----------------------------------------------------------------------

	private static final String DIST_COL_LEGACY = "Dist_frag_end";
	private static final String DIST_COL_SIGNED = "Signed_Dist_Center";

	/**
	 * Read the distance-column name from a CpG features BED's header line.
	 * Returns the canonical column name found in the header, or
	 * DIST_COL_LEGACY if no recognizable header is present (older BED files
	 * predating the option, which were always unsigned).
	 */
	private static String readDistColumnNameFromHeader(String bedFile) throws IOException {
		BufferedReader br;
		GZIPInputStream gz = null;
		if (bedFile.endsWith(".gz")) {
			gz = new GZIPInputStream(new FileInputStream(bedFile), 1 << 16);
			br = new BufferedReader(new InputStreamReader(gz), 1 << 16);
		} else {
			br = new BufferedReader(new FileReader(bedFile), 1 << 16);
		}
		try {
			String line;
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) {
					String stripped = line.substring(1);
					String[] cols = stripped.split("\t");
					for (String c : cols) {
						String t = c.trim();
						if (DIST_COL_SIGNED.equals(t) || DIST_COL_LEGACY.equals(t)) {
							return t;
						}
					}
					// Header line found but neither canonical name present.
					// Fall back to legacy: this matches older BED files that
					// used different conventions for the column header.
					return DIST_COL_LEGACY;
				}
				if (!line.isEmpty()) {
					// Hit data without finding a header; legacy BED.
					return DIST_COL_LEGACY;
				}
			}
		} finally {
			br.close();
			if (gz != null) gz.close();
		}
		return DIST_COL_LEGACY;
	}

	private static String modelMetaPath(String modelFile) {
		return modelFile + ".meta.tsv";
	}

	/**
	 * Write the model's training-time metadata sidecar. Called by trainHmm
	 * (and any other path that produces a fresh model file). Single-arg
	 * legacy form kept as a thin wrapper that records ONLY the
	 * distance-column name (dist_column key) — used by callers that don't
	 * have flag context. Newer callers should use the multi-arg form
	 * below to record the full feature-flag set, which the decode-side
	 * validator (validateFeatureFlagsAgainstModel) consumes to refuse
	 * mismatched -useEndMotif / -useBaseQ / -lowCoverage invocations.
	 */
	private static void writeModelMeta(String modelFile, String distColName) {
		writeModelMeta(modelFile, distColName, null, null, null, null);
	}

	private static void writeModelMeta(String modelFile, String distColName,
			Boolean useEndMotif, Boolean useBaseQ, Boolean lowCoverage, Integer features) {
		try (BufferedWriter bw = new BufferedWriter(new FileWriter(modelMetaPath(modelFile)))) {
			bw.write("key\tvalue\n");
			bw.write("dist_column\t" + distColName + "\n");
			if (useEndMotif != null)  bw.write("use_motif\t"    + useEndMotif  + "\n");
			if (useBaseQ != null)     bw.write("use_baseq\t"    + useBaseQ     + "\n");
			if (lowCoverage != null)  bw.write("low_coverage\t" + lowCoverage  + "\n");
			if (features != null)     bw.write("features\t"     + features     + "\n");
		} catch (IOException e) {
			// Sidecar write failure shouldn't kill training; log and continue.
			log.warn("Could not write model metadata sidecar " +
					modelMetaPath(modelFile) + ": " + e.getMessage());
		}
	}

	/**
	 * Read the entire model metadata sidecar as a key→value map. Returns
	 * an empty map when the sidecar is missing (older model trained before
	 * the metadata format existed).
	 */
	private static Map<String, String> readModelMeta(String modelFile) {
		Map<String, String> out = new HashMap<>();
		File meta = new File(modelMetaPath(modelFile));
		if (!meta.exists()) return out;
		try (BufferedReader br = new BufferedReader(new FileReader(meta))) {
			String line;
			while ((line = br.readLine()) != null) {
				if (line.startsWith("key\t")) continue;
				String[] kv = line.split("\t", 2);
				if (kv.length == 2) out.put(kv[0].trim(), kv[1].trim());
			}
		} catch (IOException e) {
			log.warn("Could not read model metadata sidecar " +
					modelMetaPath(modelFile) + ": " + e.getMessage());
		}
		return out;
	}

	/**
	 * Read the distance-column name recorded in the model's sidecar.
	 * Returns null if the sidecar is missing (older model trained before
	 * this validation existed).
	 */
	private static String readModelDistColumn(String modelFile) {
		File meta = new File(modelMetaPath(modelFile));
		if (!meta.exists()) return null;
		try (BufferedReader br = new BufferedReader(new FileReader(meta))) {
			String line;
			while ((line = br.readLine()) != null) {
				if (line.startsWith("key\t")) continue;
				String[] kv = line.split("\t", 2);
				if (kv.length == 2 && "dist_column".equals(kv[0].trim())) {
					return kv[1].trim();
				}
			}
		} catch (IOException e) {
			log.warn("Could not read model metadata sidecar " +
					modelMetaPath(modelFile) + ": " + e.getMessage());
		}
		return null;
	}

	// ----------------------------------------------------------------------
	// BED column-layout resolution.
	//
	// The CpG features BED has a fixed prefix (chr..Dist_frag_end at cols
	// 0-10) plus a variable trailing block whose layout depends on which
	// CpgFeatureMatrixBuilder flags were used: -useEndMotif adds a
	// motif_score column, -includeCpgDist adds dist_nearest_CpG, and
	// -valueWigs methyPrior:0:... adds a methyPrior column. parseLine used
	// to hard-code methyPrior at col 11 (or 12 with -useEndMotif on FinaleMe
	// itself), but this collapses three independent CpgFeatureMatrixBuilder
	// flags into one FinaleMe flag. Concretely: a BED built with
	// -useEndMotif -valueWigs methyPrior has motif_score at col 11 and
	// methyPrior at col 12 -- if FinaleMe is then called WITHOUT
	// -useEndMotif it reads motif_score as methyPrior, then divides it by
	// 100 (as if it were a percent), giving a near-zero "prior" that pulls
	// every prediction to U and tanks AUC to ~0.5.
	//
	// Fix: parse the BED's '#'-prefixed header row up front, find the
	// motif_score and methyPrior columns by NAME, and use those indices
	// throughout parseLine. Backward compat: if the BED has no header (or
	// uses the legacy un-prefixed "chr\tstart\tend\t..." form), fall back
	// to the legacy hardcoded indices.
	// ----------------------------------------------------------------------

	private int bedColMotifScore = -1;   // -1 = column not present in the BED
	private int bedColMethyPrior = -1;   // -1 = absent; will fall back below
	private int bedColDist = 10;         // legacy default: Dist_frag_end at col 10
	private boolean bedDistIsSigned = false;  // true when column is Signed_Dist_Center
	private boolean bedColumnsResolved = false;

	/**
	 * Open the BED, find the first header line (either '#'-prefixed or the
	 * legacy bare "chr start end ..." row), and resolve the motif_score and
	 * methyPrior column indices by name. If no header is present, leave the
	 * indices at -1 -- parseLine then falls back to legacy hardcoded indices
	 * driven by -useEndMotif.
	 */
	private void resolveBedColumnsFromHeader(String inputFile) throws IOException {
		if (bedColumnsResolved) return;
		BufferedReader br;
		GZIPInputStream gz = null;
		if (inputFile.endsWith(".gz")) {
			gz = new GZIPInputStream(new FileInputStream(inputFile), GZIP_BUFFER_SIZE);
			br = new BufferedReader(new InputStreamReader(gz), READER_BUFFER_SIZE);
		} else {
			br = new BufferedReader(new FileReader(inputFile), READER_BUFFER_SIZE);
		}
		try {
			String line;
			while ((line = br.readLine()) != null) {
				if (line.isEmpty()) continue;

				String[] cols;
				boolean isHeader = false;
				if (line.startsWith("#")) {
					cols = line.substring(1).split("\t");
					isHeader = true;
				} else {
					cols = line.split("\t");
					if (cols.length >= 2 && cols[1].equalsIgnoreCase("start")) {
						isHeader = true;
					}
				}

				if (isHeader) {
					Map<String, Integer> idx = new HashMap<>();
					for (int i = 0; i < cols.length; i++) {
						idx.put(cols[i].trim(), i);
					}
					Integer motifIdx = idx.get("motif_score");
					bedColMotifScore = motifIdx != null ? motifIdx : -1;
					Integer methyIdx = idx.get("methyPrior");
					bedColMethyPrior = methyIdx != null ? methyIdx : -1;
					// Resolve the distance column position AND its semantics:
					// Signed_Dist_Center  -> value is signed offset from center
					//                        (cpgOffset_from_5'end - fragLen/2),
					//                        use directly.
					// Dist_frag_end       -> value is unsigned distance to nearest
					//                        end (min(off, fragLen-1-off)),
					//                        convert via fragLen/2 - val + 0.5.
					Integer signedIdx = idx.get("Signed_Dist_Center");
					Integer unsignedIdx = idx.get("Dist_frag_end");
					if (signedIdx != null) {
						bedColDist = signedIdx;
						bedDistIsSigned = true;
					} else if (unsignedIdx != null) {
						bedColDist = unsignedIdx;
						bedDistIsSigned = false;
					} // else leave the legacy defaults (col 10, unsigned)
					log.info("BED column layout resolved from header: " +
							"motif_score=" + (bedColMotifScore >= 0 ? bedColMotifScore : "absent") +
							", methyPrior=" + (bedColMethyPrior >= 0 ? bedColMethyPrior : "absent") +
							", dist=col" + bedColDist +
							(bedDistIsSigned ? " (Signed_Dist_Center)" : " (Dist_frag_end)"));
					bedColumnsResolved = true;
					return;
				}
				// First non-empty line is data: legacy un-headered BED.
				log.info("BED has no header row; parseLine will use legacy column " +
						"indices (motif_score and methyPrior position controlled by " +
						"-useEndMotif).");
				bedColumnsResolved = true;
				return;
			}
		} finally {
			br.close();
			if (gz != null) gz.close();
		}
		bedColumnsResolved = true;
	}

	/**
	 * Validate that -useEndMotif on FinaleMe matches what the BED actually
	 * carries: -useEndMotif requires a motif_score column. Without this
	 * check, parseLine would either NaN-out the motif feature (silently bad
	 * model) or throw an ArrayIndexOutOfBoundsException at a random row.
	 */
	private void validateBedHasMotifIfRequested(String inputFile) {
		if (!useEndMotif) return;
		if (!bedColumnsResolved || bedColMotifScore < 0) {
			throw new IllegalArgumentException(
					"-useEndMotif requires the input BED " + inputFile +
					" to have a 'motif_score' column. The BED's header does not " +
					"include one, which means CpgFeatureMatrixBuilder was run " +
					"without -useEndMotif. Either drop -useEndMotif from this " +
					"FinaleMe invocation, or rebuild the feature BED with " +
					"-useEndMotif (and a -loadMotifLookup or -saveMotifLookup) in " +
					"CpgFeatureMatrixBuilder.");
		}
	}

	/**
	 * Confirm the input BED's distance-column form matches the model's
	 * training-time form. Throws on mismatch with an actionable message.
	 * If the model has no metadata sidecar (older model), warns and assumes
	 * legacy form.
	 */
	private static void validateDistColumnAgainstModel(String inputFile, String modelFile)
			throws IOException {
		String inputCol = readDistColumnNameFromHeader(inputFile);
		String modelCol = readModelDistColumn(modelFile);
		if (modelCol == null) {
			log.warn("No metadata sidecar found at " + modelMetaPath(modelFile) +
					" -- assuming this model was trained with the legacy " +
					DIST_COL_LEGACY + " distance form. The current input BED " +
					"reports " + inputCol + ". If those don't match, decode " +
					"results will be silently wrong; rebuild features with " +
					"the matching -useSignedDistCenter setting or re-train " +
					"the model on the matching feature form.");
			return;
		}
		if (!modelCol.equals(inputCol)) {
			throw new IllegalArgumentException(
					"Distance-column mismatch: model " + modelFile +
					" was trained with '" + modelCol + "' (per " +
					modelMetaPath(modelFile) + "), but input BED " + inputFile +
					" carries '" + inputCol + "'. The HMM reads features by " +
					"column index, so this mismatch would silently produce " +
					"wrong predictions. Rebuild the input BED with " +
					(DIST_COL_SIGNED.equals(modelCol)
							? "-useSignedDistCenter"
							: "the legacy unsigned distance (drop -useSignedDistCenter)") +
					" in CpgFeatureMatrixBuilder, OR re-train the model on " +
					"a BED matching the input's '" + inputCol + "' form.");
		}
	}

	/**
	 * Validate that the current invocation's feature-flag set matches what
	 * the model was trained with. Sidecar keys: use_motif, use_baseq,
	 * low_coverage. The HMM reads features by column INDEX, so toggling a
	 * flag between train and decode silently re-aligns observation vector
	 * dimensions to the wrong stats and produces nonsense predictions.
	 *
	 * Backward compat: if the sidecar is missing or doesn't have a
	 * particular key (older models trained pre-v0.64), log a warning and
	 * proceed assuming the current invocation is correct.
	 */
	private void validateFeatureFlagsAgainstModel(String modelFile) {
		Map<String, String> meta = readModelMeta(modelFile);
		if (meta.isEmpty()) {
			// Already covered by the dist-column validator's warning;
			// don't double-log.
			return;
		}
		checkFlagAgainstModel(meta, "use_motif",    "useEndMotif",  useEndMotif,  modelFile);
		checkFlagAgainstModel(meta, "use_baseq",    "useBaseQ",     useBaseQ,     modelFile);
		checkFlagAgainstModel(meta, "low_coverage", "lowCoverage",  lowCoverage,  modelFile);
	}

	private static void checkFlagAgainstModel(Map<String, String> meta,
			String key, String flagName, boolean current, String modelFile) {
		String trained = meta.get(key);
		if (trained == null) return;  // older sidecar without this key; tolerate.
		boolean trainedBool = Boolean.parseBoolean(trained);
		if (trainedBool != current) {
			throw new IllegalArgumentException(
					"Feature-flag mismatch: model " + modelFile + " was trained " +
					"with -" + flagName + "=" + trainedBool + " (per " +
					modelMetaPath(modelFile) + "), but the current invocation " +
					"has -" + flagName + "=" + current + ". The HMM reads " +
					"features by column index, so this mismatch would silently " +
					"produce wrong predictions. Toggle -" + flagName + " to " +
					"match the trained model, or retrain on a BED matching the " +
					"current invocation's flags.");
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
		// Header lines start with '#' (current format) or have "start" in
		// the second column (legacy un-prefixed form). Either way, skip.
		if (line.startsWith("#") || splitLines[1].equalsIgnoreCase("start")) {
			return null;
		}
		// Resolve which columns to read. The header was parsed up-front by
		// resolveBedColumnsFromHeader and stored in instance state; if that
		// returned absent (legacy un-headered BED), fall back to the legacy
		// hardcoded indices keyed off -useEndMotif so existing pipelines
		// keep working.
		int motifScoreCol = bedColMotifScore;
		int methyPriorCol = bedColMethyPrior;
		if (methyPriorCol < 0) {
			// Legacy fallback when the BED has no header at all.
			methyPriorCol = useEndMotif ? 12 : 11;
			if (motifScoreCol < 0 && useEndMotif) {
				motifScoreCol = 11;
			}
		}
		int minCols = Math.max(methyPriorCol + 1,
				motifScoreCol >= 0 ? motifScoreCol + 1 : 0);
		if (splitLines.length < minCols
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

		// motif_score: read only when both (a) the BED carries the column
		// AND (b) the user requested -useEndMotif so the HMM consumes it.
		// Reading without (a) is what produced the silent prior-shift bug
		// that this whole header-resolution mechanism exists to prevent.
		double motifScore = Double.NaN;
		if (useEndMotif && motifScoreCol >= 0) {
			motifScore = Double.parseDouble(splitLines[motifScoreCol]);
		}

		double methyPrior = Double.parseDouble(splitLines[methyPriorCol]);
		if (Double.isNaN(methyPrior)) return null;

		double fragLen = Double.parseDouble(splitLines[4]);
		double coverage = Double.parseDouble(splitLines[7]);
		// Distance column has two distinct value semantics, indicated by
		// the BED's column header (resolved up-front in
		// resolveBedColumnsFromHeader):
		//   - Dist_frag_end (legacy): unsigned distance to nearest end,
		//     value range [0, fragLen/2]. Convert with fragLen/2 - val + 0.5
		//     (matches the original FinaleMe paper's formula -- gives
		//     unsigned distance from center, range [0.5, fragLen/2+0.5]).
		//   - Signed_Dist_Center: already a signed distance from center,
		//     value range [-fragLen/2, +fragLen/2]. Use as-is.
		// Hard-coding splitLines[10] + the legacy formula is now wrong for
		// signed BEDs and would silently flip / scale every value.
		double rawDist = Double.parseDouble(splitLines[bedColDist]);
		double distToCenter = bedDistIsSigned
				? rawDist
				: (fragLen / 2 - rawDist + 0.5);

		if (Double.compare(methyPrior, 100.0) == 0) {
			methyPrior -= 0.01;
		} else if (Double.compare(methyPrior, 0.0) == 0) {
			methyPrior += 0.01;
		}
		methyPrior /= 100;
		if (Double.isNaN(methyPrior)) return null;

		// Per-CpG baseQ (Phred). Already read above for the >5 filter; capture
		// it explicitly into ParsedRow so downstream paths can route it into
		// the HMM observation vector when -useBaseQ is set.
		double baseQ = Double.parseDouble(splitLines[8]);

		String loc = chr + ":" + start + ":" + end;
		String readName = splitLines[3];
		return new ParsedRow(readName, splitLines[6], loc, offset, methyPrior, fragLen, coverage, distToCenter, motifScore, baseQ);
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

	/**
	 * Index of the baseQ slot in stats[]. -1 when -useBaseQ is off; equals 3
	 * when motif is off, 4 when motif is on. baseQ goes after motif so it
	 * doesn't disturb the existing motif slot index.
	 */
	private int baseQStatsIndex() {
		if (!useBaseQ) return -1;
		return useEndMotif ? 4 : 3;
	}

	/**
	 * Map an observation-vector dimension index (in the order produced by
	 * buildObservationVector) to its slot in the stats[] array. Used by
	 * adaptation re-centering to look up the right normalization stats per
	 * feature. Layout (with flag legend in parentheses):
	 *
	 *   d=0: FragLen        → stats[0]
	 *   d=1: Coverage       → stats[1]   (only when !lowCoverage)
	 *        else d=1: DistToCenter → stats[2]
	 *   d=2: DistToCenter   → stats[2]   (when !lowCoverage)
	 *        else d=2: MotifScore → stats[3]    (when lowCoverage && useEndMotif)
	 *        or  d=2: baseQ        → bqIdx     (when lowCoverage && useBaseQ only)
	 *   d=3: MotifScore     → stats[3]   (when !lowCoverage && useEndMotif)
	 *        or  d=3: baseQ → bqIdx     (when !lowCoverage && !useEndMotif && useBaseQ)
	 *   d=4: baseQ          → bqIdx     (when !lowCoverage && useEndMotif && useBaseQ)
	 *
	 * Equivalent to walking the same conditional blocks as
	 * buildObservationVector and incrementing through the stats indices.
	 */
	private int obsDimToStatsIdx(int dim) {
		int d = 0;
		// FragLen at d=0 -> stats[0]
		if (dim == d) return 0;
		d++;
		if (!lowCoverage) {
			if (dim == d) return 1; // Coverage
			d++;
		}
		if (dim == d) return 2; // DistToCenter
		d++;
		if (useEndMotif) {
			if (dim == d) return 3; // MotifScore
			d++;
		}
		if (useBaseQ) {
			if (dim == d) return baseQStatsIndex(); // baseQ
			d++;
		}
		// Defensive: shouldn't happen if the caller respects observation-vector dimensions.
		throw new IllegalStateException("obsDimToStatsIdx: dimension " + dim +
				" out of range for current flag configuration");
	}

	/**
	 * Build the z-scored observation vector from a parsed row according to
	 * the active feature flags. Handles all four legacy cases plus the two
	 * new cases that opt-in -useBaseQ. baseQ is appended at the END of the
	 * vector to keep existing slot indices stable for motif-bearing models.
	 *
	 * Resulting dimension equals 3 + (useEndMotif ? 1 : 0) + (useBaseQ ? 1 : 0)
	 * - (lowCoverage ? 1 : 0).
	 */
	private double[] buildObservationVector(ParsedRow row, SummaryStatistics[] stats, int bqIdx) {
		java.util.ArrayList<Double> dims = new java.util.ArrayList<>(stats.length);
		dims.add((row.fragLen - stats[0].getMean()) / stats[0].getStandardDeviation());
		if (!lowCoverage) {
			dims.add((row.coverage - stats[1].getMean()) / stats[1].getStandardDeviation());
		}
		dims.add((row.distToCenter - stats[2].getMean()) / stats[2].getStandardDeviation());
		if (useEndMotif) {
			dims.add((row.motifScore - stats[3].getMean()) / stats[3].getStandardDeviation());
		}
		if (useBaseQ && bqIdx >= 0) {
			double sd = stats[bqIdx].getStandardDeviation();
			double zscored = sd > 0
					? (row.baseQ - stats[bqIdx].getMean()) / sd
					: 0.0;  // SD=0 already caught by the degeneracy guard; defensive fallback
			dims.add(zscored);
		}
		double[] out = new double[dims.size()];
		for (int i = 0; i < out.length; i++) out[i] = dims.get(i);
		return out;
	}

	private SummaryStatistics[] collectStats(String matrixFile,
											 HashMap<String, IntervalTree<Integer>> overlapLoc,
											 HashMap<String, IntervalTree<Integer>> excludeLoc) throws IOException {
		int numStats = 3 + (useEndMotif ? 1 : 0) + (useBaseQ ? 1 : 0);
		SummaryStatistics[] stats = new SummaryStatistics[numStats];
		for (int i = 0; i < numStats; i++) {
			stats[i] = new SummaryStatistics();
		}
		int bqIdx = baseQStatsIndex();

		BufferedReader br = openGzipReader(matrixFile);

		String line;
		while ((line = br.readLine()) != null) {
			ParsedRow row = parseLine(line, overlapLoc, excludeLoc);
			if (row == null) continue;
			stats[0].addValue(row.fragLen);
			stats[1].addValue(row.coverage);
			stats[2].addValue(row.distToCenter);
			if (useEndMotif && !Double.isNaN(row.motifScore)) {
				stats[3].addValue(row.motifScore);
			}
			if (useBaseQ && bqIdx >= 0) {
				stats[bqIdx].addValue(row.baseQ);
			}
		}
		br.close();

		return stats;
	}

	private void logFeatureStats(SummaryStatistics[] stats) {
		// Build dynamic feature-name list based on which flags are active.
		// Slots 0..2 are always FragLen / Norm_Frag_cov / DistToCenter.
		// Slot 3 = MotifScore (if -useEndMotif) OR baseQ (if -useBaseQ alone)
		// Slot 4 = baseQ (only when both -useEndMotif AND -useBaseQ)
		String[] featureNames = new String[stats.length];
		featureNames[0] = "FragLen";
		if (stats.length > 1) featureNames[1] = "Norm_Frag_cov";
		if (stats.length > 2) featureNames[2] = "DistToCenter";
		int nextSlot = 3;
		if (useEndMotif && nextSlot < stats.length) {
			featureNames[nextSlot++] = "MotifScore";
		}
		if (useBaseQ && nextSlot < stats.length) {
			featureNames[nextSlot++] = "baseQ";
		}
		for (int i = 0; i < stats.length; i++) {
			String featureName = featureNames[i] != null ? featureNames[i] : ("Feature" + i);
			logFeatureStat(i, featureName, stats[i]);
		}

		// Degeneracy guard for -useBaseQ: tabix-fragment input writes a
		// CONSTANT baseQ (= -fragBaseQ default), and feeding a constant
		// feature into the GMM gives a singular covariance matrix at
		// init. Fail fast with an actionable message rather than letting
		// EigenDecomposition.getInverse() throw a cryptic stack trace.
		int bqIdx = baseQStatsIndex();
		if (bqIdx >= 0 && bqIdx < stats.length) {
			SummaryStatistics bqStat = stats[bqIdx];
			if (bqStat.getN() > 1 && bqStat.getStandardDeviation() == 0.0) {
				throw new IllegalStateException(
						"-useBaseQ feature has SD=0 across all CpGs (constant " +
						"baseQ = " + bqStat.getMean() + " over " + bqStat.getN() + " " +
						"observations). This typically means the BED came from " +
						"tabix-fragment input where baseQ is a constant from " +
						"-fragBaseQ. Either drop -useBaseQ, or rebuild the BED " +
						"from a BAM input where baseQ varies per CpG site.");
			}
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

					// Resolve where motif_score / methyPrior live in the BED by
					// reading the '#'-prefixed header. parseLine uses the resolved
					// indices instead of guessing from the -useEndMotif flag, which
					// silently mis-parses when the two flag-sets diverge between
					// CpgFeatureMatrixBuilder (which determines what's in the BED)
					// and FinaleMe (which determines what gets used).
					resolveBedColumnsFromHeader(inputFile);
					validateBedHasMotifIfRequested(inputFile);

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
									trainHmm(matrixObj2, modelFile, inputFile);
									matrixObj2 = null; // allow GC of second copy
								}else{
									trainHmm(matrixObj, modelFile, inputFile);
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
		return processMatrixFile(matrixFile, false);
	}

	/**
	 * Process matrix file for adaptation/training.
	 *
	 * @param matrixFile Input feature matrix (possibly bgzipped).
	 * @param minimalForAdaptation When true, skip building the auxiliary structures
	 *   (matrixU, matrixM, pi, a, matrixObserved, cpgDistFreq) that are only used
	 *   by the normal training's random/GMM initialization paths. The adaptation
	 *   path uses only the assembled `matrix` list, so skipping these saves a
	 *   second full copy of observation vector references (~12-20 GB at 30X WGS).
	 */
	private MatrixObj processMatrixFile(String matrixFile, boolean minimalForAdaptation) throws FileNotFoundException, IOException, FileFormatException{
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
			int numStats = 3 + (useEndMotif ? 1 : 0) + (useBaseQ ? 1 : 0);
			SummaryStatistics[] stats = new SummaryStatistics[numStats];
			for(int i = 0; i < numStats; i++){
				stats[i] = new SummaryStatistics();
			}
			int bqIdx = baseQStatsIndex();
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
				if(useEndMotif && !Double.isNaN(row.motifScore)){
					stats[3].addValue(row.motifScore);
				}
				if(useBaseQ && bqIdx >= 0){
					stats[bqIdx].addValue(row.baseQ);
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
		// Convert the count map into a compact HashSet of qualifying reads.
		// Saves ~75% memory vs keeping full HashMap<String,int[]> through Pass 2:
		// at 150M reads, the HashMap is ~22GB vs the HashSet ~600MB (after filter).
		HashSet<String> qualifyingReadSet = new HashSet<>(Math.max(16, (int)(totalReads / 20)));
		for(java.util.Map.Entry<String, int[]> e : readCpgCount.entrySet()){
			if(e.getValue()[0] >= miniDataPoints){
				qualifyingReadSet.add(e.getKey());
			}
		}
		long qualifyingReads = qualifyingReadSet.size();
		readCpgCount = null; // drop the full count map — frees ~22GB for 150M reads
		log.info("Pass 1 done: " + statsRows + " rows, " + totalReads +
				 " unique reads, " + qualifyingReads + " with >= " + miniDataPoints + " CpGs");
		log.info("Dropped readCpgCount map; retaining " + qualifyingReads + "-entry qualifying-read set");

		// === Pass 2: Re-read file, skip non-qualifying reads, build matrixProcess ===
		log.info("Pass 2: Building observation vectors (skipping " +
				 (totalReads - qualifyingReads) + " reads with < " + miniDataPoints + " CpGs) ...");
		BufferedReader br2 = openGzipReader(matrixFile);
		long skippedRows = 0;

		while( (line = br2.readLine()) != null){
			ParsedRow row = parseLine(line, overlapLoc, excludeLoc);
			if(row == null) continue;

			// Early skip: if this read is not in the qualifying set (< miniDataPoints CpGs),
			// skip it now to avoid allocating ObservationVector/Triple/Pair/loc String.
			if(!qualifyingReadSet.contains(row.readName)){
				skippedRows++;
				continue;
			}

			if(covOutlier > 0 && ((row.coverage-stats[1].getMean())/stats[1].getStandardDeviation() > covOutlier ||
					(row.fragLen-stats[0].getMean())/stats[0].getStandardDeviation() > covOutlier ||
					(row.distToCenter-stats[2].getMean())/stats[2].getStandardDeviation() > covOutlier ||
					(useBaseQ && bqIdx >= 0 && stats[bqIdx].getStandardDeviation() > 0
							&& (row.baseQ - stats[bqIdx].getMean()) / stats[bqIdx].getStandardDeviation() > covOutlier))){
				continue;
			}
			double[] value = buildObservationVector(row, stats, bqIdx);

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
		qualifyingReadSet = null; // free the qualifying-read set
		log.info("Pass 2 done: " + points + " points loaded, " + skippedRows + " rows skipped");
			log.info("Number of point in total is loaded : " + points);

			// Auxiliary structures populated below; they are only consumed by
			// the normal training's random/GMM initialization and the final
			// decode bookkeeping. When minimalForAdaptation is true we keep
			// them empty to avoid a second copy of observation vector refs
			// (~12-20 GB at 30X WGS).
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
					// Skip the expensive per-CpG bookkeeping when only the assembled
					// matrix is needed (adaptation path).
					if (minimalForAdaptation) {
						double cpgDenseMinimal = (double) matrixRow.size();
						if (cpgDenseMinimal > observedMaxCpgNum) {
							observedMaxCpgNum = cpgDenseMinimal;
						}
						if (cpgNumClip < 0 && cpgDenseMinimal > this.maxCpgNum) {
							this.maxCpgNum = cpgDenseMinimal;
						}
						continue; // skip matrixU/M/pi/aij/observed population
					}
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
		trainHmm(matrixObj, modelFile, null);
	}

	/**
	 * Train + save the HMM model AND record the input BED's distance-column
	 * form in a sidecar metadata file. trainingInputFile is used solely to
	 * read the BED header for that purpose; if null, the sidecar is skipped
	 * (legacy callers).
	 */
	private void trainHmm(MatrixObj matrixObj, String modelFile, String trainingInputFile) throws IOException, CloneNotSupportedException{
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

		// Record the training-time distance-column form AND feature-flag
		// set so future decode runs can refuse mismatched BEDs / flags. See
		// validateDistColumnAgainstModel and validateFeatureFlagsAgainstModel.
		if (trainingInputFile != null) {
			String distCol = readDistColumnNameFromHeader(trainingInputFile);
			log.info("Recording training distance column '" + distCol +
					"' and feature flags (use_motif=" + useEndMotif +
					", use_baseq=" + useBaseQ + ", low_coverage=" + lowCoverage +
					", features=" + features + ") to " + modelMetaPath(modelFile));
			writeModelMeta(modelFile, distCol, useEndMotif, useBaseQ, lowCoverage, features);
		}

	}


	//decoding HMM
	private double decodeHmm(MatrixObj matrixObj, String hmmFile, String outputFile, String inputFile, boolean reestimate, CpgIndex cpgIndex, LinkedHashMap<String, Integer> chromOrder) throws Exception{
		// Reject mismatched distance-column form. In the train+decode path
		// trainHmm has just written the metadata sidecar matching this same
		// input, so this passes by construction; the check matters when
		// decodeHmm is invoked with a pre-existing model.
		validateDistColumnAgainstModel(inputFile, hmmFile);
		validateFeatureFlagsAgainstModel(hmmFile);
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
		decodeOnlyStreaming(inputFile, modelFile, outputFile, cpgIndex, chromOrder, null);
	}

	/**
	 * Streaming decode with optional pre-computed stats.
	 * When {@code precomputedStats != null}, Pass 1 (collectStats) is skipped —
	 * saving one full pass over the (large) bgzipped input when the caller
	 * already computed the stats (e.g. adaptAndDecodeStreaming).
	 */
	private void decodeOnlyStreaming(String inputFile, String modelFile, String outputFile,
									CpgIndex cpgIndex, LinkedHashMap<String, Integer> chromOrder,
									SummaryStatistics[] precomputedStats) throws Exception {
		System.out.println("\nStreaming decode-only mode ...\n");

		// Reject input BEDs whose distance-column form differs from the
		// model's training-time form (silent mis-prediction otherwise).
		validateDistColumnAgainstModel(inputFile, modelFile);
		validateFeatureFlagsAgainstModel(modelFile);

		// Load region/exclude intervals
		HashMap<String, IntervalTree<Integer>> overlapLoc = loadIntervalFile(region);
		HashMap<String, IntervalTree<Integer>> excludeLoc = loadIntervalFile(exclude);

		SummaryStatistics[] stats;
		if (precomputedStats != null) {
			log.info("Phase 1: Using pre-computed feature statistics (skipping file scan) ...");
			stats = precomputedStats;
			logFeatureStats(stats);
		} else {
			// Phase 1: Collect stats for z-score normalization
			log.info("Phase 1: Collecting feature statistics ...");
			stats = collectStats(inputFile, overlapLoc, excludeLoc);
			logFeatureStats(stats);
		}

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
		int bqIdxStream = baseQStatsIndex();
		while ((line = br.readLine()) != null) {
			ParsedRow row = parseLine(line, overlapLoc, excludeLoc);
			if (row == null) continue;

			// Apply covOutlier filter
			if (covOutlier > 0 && ((row.coverage - stats[1].getMean()) / stats[1].getStandardDeviation() > covOutlier ||
					(row.fragLen - stats[0].getMean()) / stats[0].getStandardDeviation() > covOutlier ||
					(row.distToCenter - stats[2].getMean()) / stats[2].getStandardDeviation() > covOutlier ||
					(useBaseQ && bqIdxStream >= 0 && stats[bqIdxStream].getStandardDeviation() > 0
							&& (row.baseQ - stats[bqIdxStream].getMean()) / stats[bqIdxStream].getStandardDeviation() > covOutlier))) {
				continue;
			}

			// Z-score normalize via the shared helper so the streaming-decode
			// path stays in sync with the training path's vector layout.
			double[] value = buildObservationVector(row, stats, bqIdxStream);
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
		// Validate feature-flag consistency with the trained model. The AUC
		// path shares parseLine + processMatrixFile with the training path,
		// so the same flag-mismatch silent-corruption applies if a user
		// trains with one flag set and AUC-evaluates with another.
		validateFeatureFlagsAgainstModel(hmmFile);
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

		// Threshold sweep. -aucNThresholds 0 (default) keeps the legacy
		// 19-point non-uniform list (denser at the extremes); >0 generates
		// that many evenly-spaced points in [-1, 1] for smoother curves.
		double[] ps;
		if(aucNThresholds > 0){
			ps = new double[aucNThresholds];
			for(int i = 0; i < aucNThresholds; i++){
				ps[i] = -1.0 + 2.0 * i / (aucNThresholds - 1);
			}
		}else{
			ps = new double[]{-1, -0.999, -0.99, -0.95, -0.9, -0.8, -0.7, -0.5, -0.1, 0.0,
					0.1, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999, 1.0};
		}

		// Per-threshold confusion-matrix entries. We retain everything so
		// AUROC / AUPRC can be computed by trapezoidal integration after
		// the sweep, and the curve data can be written out as TSV / PDF.
		int n = ps.length;
		long[] tps = new long[n], fns = new long[n], fps = new long[n], tns = new long[n];

		// Parallelize the per-fragment Viterbi inner loop. The HMM object
		// is read-only after loadHmmModel + setters, so each task can
		// instantiate its own ViterbiBayesianNhmmV5Calculator over a
		// disjoint chunk of fragments and reduce four long counters into
		// a shared array. The threshold sweep stays sequential (each
		// threshold writes one stdout line, and the legacy stdout order
		// matters); the heavy work (Viterbi) is what gets parallelized.
		final int nThreads = resolveThreadCount();
		final int matrixSize = matrix.size();
		final int chunkSize = Math.max(1, (matrixSize + nThreads - 1) / nThreads);
		log.info("AUC sweep: " + n + " thresholds x " + matrixSize +
				" fragments using " + nThreads + " threads");
		ExecutorService aucExecutor = Executors.newFixedThreadPool(nThreads);
		final List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> matrixFinal = matrix;
		final ArrayList<ArrayList<Integer>> matrixObservedFinal = matrixObserved;

		boolean outFlag = false;
		try {
			for(int t = 0; t < n; t++){
				final double p = ps[t];
				final int thresholdIdx = t;
				List<Future<long[]>> futures = new ArrayList<>(nThreads);
				for(int chunkStart = 0; chunkStart < matrixSize; chunkStart += chunkSize){
					final int cs = chunkStart;
					final int ce = Math.min(cs + chunkSize, matrixSize);
					futures.add(aucExecutor.submit(() -> {
						long lMethy = 0, lMethyCorrect = 0, lUnmethy = 0, lUnmethyCorrect = 0;
						// Per-task RNG keeps -randomPerm reproducible regardless
						// of thread count: seed = global-seed XOR threshold-index
						// XOR chunk-start. MersenneTwister is not thread-safe,
						// so we cannot share `randomEngine` across tasks.
						MersenneTwister localRng = randomPerm
								? new MersenneTwister(0x9E3779B97F4A7C15L
										^ ((long) thresholdIdx << 32) ^ (long) cs)
								: null;
						for(int j = cs; j < ce; j++){
							int[] hiddenState = (new ViterbiBayesianNhmmV5Calculator(
									matrixFinal.get(j), hmm, methylatedState, p)).stateSequence();
							Integer[] observedState = matrixObservedFinal.get(j)
									.toArray(new Integer[matrixObservedFinal.get(j).size()]);
							if(hiddenState.length != observedState.length){
								throw new IllegalArgumentException(
										"HiddenState Length does not match with observed state length");
							}
							if(randomPerm){
								HashMap<Integer, Pair<Integer, Double>> cpgDistState =
										matrixFinal.get(j).getFirst();
								for(int i = 0; i < observedState.length; i++){
									Double methyPrior = cpgDistState.get(i).getSecond();
									double rand = localRng.nextDouble();
									if(rand < methyPrior + p){
										hiddenState[i] = 1;
									}else{
										hiddenState[i] = 0;
									}
								}
							}
							for(int i = 0; i < hiddenState.length; i++){
								if(hiddenState[i] % 2 == observedState[i]){
									if(observedState[i] == 1){
										lMethyCorrect++;
									}else{
										lUnmethyCorrect++;
									}
								}
								if(observedState[i] == 1){
									lMethy++;
								}else{
									lUnmethy++;
								}
							}
						}
						return new long[]{lMethy, lMethyCorrect, lUnmethy, lUnmethyCorrect};
					}));
				}
				long countMethy = 0, countMethyCorrect = 0, countUnmethy = 0, countUnmethyCorrect = 0;
				for(Future<long[]> f : futures){
					long[] c = f.get();
					countMethy        += c[0];
					countMethyCorrect += c[1];
					countUnmethy      += c[2];
					countUnmethyCorrect += c[3];
				}

				if(!outFlag){
					System.out.println(countMethy + "\t" + countUnmethy);
					outFlag = true;
				}
				double fpr = (double)(countUnmethy - countUnmethyCorrect) / (double)countUnmethy;
				double tpr = (double)(countMethyCorrect) / (double)countMethy;
				System.out.println(fpr + "\t" + tpr); // legacy stdout, preserved

				tps[t] = countMethyCorrect;
				fns[t] = countMethy - countMethyCorrect;
				fps[t] = countUnmethy - countUnmethyCorrect;
				tns[t] = countUnmethyCorrect;
			}
		} finally {
			aucExecutor.shutdown();
			aucExecutor.awaitTermination(10, TimeUnit.SECONDS);
		}

		// Compute AUROC and AUPRC by trapezoidal integration. The PR curve
		// gets a (recall=0, precision=1) virtual endpoint to anchor the
		// left side -- standard convention when the threshold sweep doesn't
		// reach 0% recall by itself.
		double[] fprArr = new double[n];
		double[] tprArr = new double[n];
		double[] precArr = new double[n];
		double[] recArr = new double[n];
		for(int t = 0; t < n; t++){
			double tp = tps[t], fn = fns[t], fp = fps[t], tn = tns[t];
			fprArr[t] = (fp + tn) > 0 ? fp / (fp + tn) : 0.0;
			tprArr[t] = (tp + fn) > 0 ? tp / (tp + fn) : 0.0;
			precArr[t] = (tp + fp) > 0 ? tp / (tp + fp) : 1.0; // undefined -> 1.0 (no positives called)
			recArr[t] = tprArr[t]; // recall == TPR
		}
		double auroc = trapezoidAUC(fprArr, tprArr);
		double[] recForPR = new double[n + 1];
		double[] precForPR = new double[n + 1];
		recForPR[0] = 0.0; precForPR[0] = 1.0; // anchor
		System.arraycopy(recArr, 0, recForPR, 1, n);
		System.arraycopy(precArr, 0, precForPR, 1, n);
		double auprc = trapezoidAUC(recForPR, precForPR);

		System.out.println();
		System.out.println("AUROC: " + String.format("%.6f", auroc));
		System.out.println("AUPRC: " + String.format("%.6f", auprc));
		log.info(String.format("AUROC = %.6f  AUPRC = %.6f  (n_methylated=%d, n_unmethylated=%d, %d thresholds)",
				auroc, auprc, tps[0] + fns[0], fps[0] + tns[0], n));

		// TSV with per-threshold curve data. Forced-on if -aucCurvePdf is
		// requested (since the Python plotter consumes the TSV).
		String tsvPath = aucCurveTsv;
		if(tsvPath == null && aucCurvePdf != null){
			tsvPath = aucCurvePdf + ".curve.tsv";
			log.info("-aucCurvePdf set without -aucCurveTsv; writing curve data to " + tsvPath);
		}
		if(tsvPath != null){
			writeAucCurveTsv(tsvPath, ps, tps, fns, fps, tns, fprArr, tprArr, precArr, recArr,
					auroc, auprc);
			log.info("Wrote ROC/PR curve data: " + tsvPath);
		}
		if(aucCurvePdf != null){
			renderAucCurvePdf(tsvPath, aucCurvePdf, auroc, auprc);
		}
	}

	/**
	 * Trapezoidal AUC. Sorts the input pairs by x ascending, then sums
	 * (x_i+1 - x_i) * (y_i + y_i+1) / 2 over consecutive pairs. NaN /
	 * Infinity inputs are skipped.
	 */
	private static double trapezoidAUC(double[] xs, double[] ys){
		if(xs.length != ys.length || xs.length < 2) return 0.0;
		int n = xs.length;
		double[][] pairs = new double[n][2];
		for(int i = 0; i < n; i++){
			pairs[i][0] = xs[i];
			pairs[i][1] = ys[i];
		}
		java.util.Arrays.sort(pairs, (a, b) -> Double.compare(a[0], b[0]));
		double area = 0.0;
		for(int i = 1; i < n; i++){
			double x0 = pairs[i - 1][0], y0 = pairs[i - 1][1];
			double x1 = pairs[i][0], y1 = pairs[i][1];
			if(Double.isNaN(x0) || Double.isNaN(x1) || Double.isNaN(y0) || Double.isNaN(y1)) continue;
			area += (x1 - x0) * (y0 + y1) / 2.0;
		}
		return area;
	}

	private static void writeAucCurveTsv(String path, double[] ps, long[] tps, long[] fns,
			long[] fps, long[] tns, double[] fpr, double[] tpr, double[] prec, double[] rec,
			double auroc, double auprc) throws IOException{
		try(BufferedWriter bw = new BufferedWriter(new FileWriter(path))){
			bw.write("# AUROC\t" + String.format("%.6f", auroc) + "\n");
			bw.write("# AUPRC\t" + String.format("%.6f", auprc) + "\n");
			bw.write("# n_methylated\t" + (tps[0] + fns[0]) + "\n");
			bw.write("# n_unmethylated\t" + (fps[0] + tns[0]) + "\n");
			bw.write("threshold\tTP\tFN\tFP\tTN\tFPR\tTPR\tPrecision\tRecall\n");
			for(int t = 0; t < ps.length; t++){
				bw.write(String.format("%.6f\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f%n",
						ps[t], tps[t], fns[t], fps[t], tns[t],
						fpr[t], tpr[t], prec[t], rec[t]));
			}
		}
	}

	/**
	 * Render the ROC + PR curves to PDF by invoking scripts/plot_roc_curve.py.
	 * Tries CWD-relative path first, then JAR-relative ../scripts. The Python
	 * helper requires matplotlib; if it's missing or python3 isn't on PATH,
	 * we log a clear error pointing at the TSV (which already has all the
	 * data the user needs to plot externally).
	 */
	private static void renderAucCurvePdf(String tsvPath, String pdfPath, double auroc, double auprc) {
		File script = null;
		File cwdScript = new File("scripts/plot_roc_curve.py");
		if(cwdScript.exists()){
			script = cwdScript;
		}else{
			try {
				File jarLoc = new File(FinaleMe.class.getProtectionDomain()
						.getCodeSource().getLocation().toURI());
				File rel = new File(jarLoc.getParentFile().getParentFile(),
						"scripts/plot_roc_curve.py");
				if(rel.exists()) script = rel;
			} catch (Exception e) { /* ignore */ }
		}
		if(script == null){
			log.warn("Could not locate scripts/plot_roc_curve.py to render " + pdfPath +
					"; the curve data is in " + tsvPath +
					" -- plot it manually with: python3 scripts/plot_roc_curve.py " +
					tsvPath + " -o " + pdfPath);
			return;
		}
		try {
			ProcessBuilder pb = new ProcessBuilder("python3", script.getAbsolutePath(),
					tsvPath, "-o", pdfPath);
			pb.redirectErrorStream(true);
			Process proc = pb.start();
			java.io.BufferedReader reader = new java.io.BufferedReader(
					new java.io.InputStreamReader(proc.getInputStream()));
			StringBuilder out = new StringBuilder();
			String line;
			while((line = reader.readLine()) != null){
				out.append(line).append("\n");
			}
			int rc = proc.waitFor();
			if(rc != 0){
				log.warn("plot_roc_curve.py exited " + rc + "; PDF not produced. " +
						"Output:\n" + out + "\nThe curve data is in " + tsvPath +
						" -- you can plot it manually.");
			}else{
				log.info("Wrote ROC/PR PDF: " + pdfPath);
			}
		} catch (Exception e) {
			log.warn("Could not invoke python3 scripts/plot_roc_curve.py: " +
					e.getMessage() + ". The curve data is in " + tsvPath +
					" -- plot it manually if needed.");
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
		// Sanity check: observation vector dimension must equal features.
		// Mismatch typically means -features N was inconsistent with the file's
		// feature columns (e.g. -useEndMotif mismatch with file schema).
		if (!matrix.isEmpty() && !matrix.get(0).getSecond().isEmpty()) {
			int obsDim = matrix.get(0).getSecond().get(0).dimension();
			if (obsDim != features) {
				int expected = 3 + (useEndMotif ? 1 : 0) + (useBaseQ ? 1 : 0) - (lowCoverage ? 1 : 0);
				throw new IllegalStateException(
					"Observation vector dimension (" + obsDim + ") does not match -features (" + features + "). " +
					"For the current flag set (-useEndMotif=" + useEndMotif + ", -useBaseQ=" + useBaseQ +
					", -lowCoverage=" + lowCoverage + ") the expected -features value is " + expected + ". " +
					"Set -features=" + expected + " or toggle the relevant flags so the HMM dimension matches " +
					"the input observation vector.");
			}
		}
		// Clip z-scored feature values to +/-OBS_CLIP_SD to prevent extreme
		// outliers (e.g. 0.5 default motif score when target motif mean=0.99
		// sd=0.001 produces z-score ~-420) from making the GMM covariance
		// singular and returning NaN emissions. Replace NaN/Inf with 0.
		final double OBS_CLIP_SD = 5.0;
		List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> clippedMatrix =
			new ArrayList<>(matrix.size());
		for (Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>> p : matrix) {
			ArrayList<ObservationVector> clippedObs = new ArrayList<>(p.getSecond().size());
			for (ObservationVector ov : p.getSecond()) {
				double[] vals = ov.values();
				double[] clipped = new double[vals.length];
				for (int d = 0; d < vals.length; d++) {
					double v = vals[d];
					if (Double.isNaN(v) || Double.isInfinite(v)) v = 0.0;
					clipped[d] = Math.max(-OBS_CLIP_SD, Math.min(OBS_CLIP_SD, v));
				}
				clippedObs.add(new ObservationVector(clipped));
			}
			clippedMatrix.add(new Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>(
				p.getFirst(), clippedObs));
		}
		GMMLearner kl = new GMMLearner(states, new OpdfMultiMixtureGaussianFactory(features, mixNumberInFeature),clippedMatrix,maxCpgDist/bin,features, mixNumberInFeature,bayesianFactor, randomEngine, tolKmeans,decayKmeans, cpgNumClip, 1, lowCoverage);
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
		if(adaptReinitGmm && !adaptEmissionOnly){
			throw new IllegalArgumentException("-adaptReinitGmm requires -adaptEmissionOnly");
		}
		if(adaptTransitions && !adaptEmissionOnly){
			throw new IllegalArgumentException("-adaptTransitions requires -adaptEmissionOnly");
		}
		if(autoTuneBayesianFactor && !adaptEmissionOnly){
			throw new IllegalArgumentException("-autoTuneBayesianFactor requires -adaptEmissionOnly");
		}
		if(autoAdaptLambda && !adaptEmissionOnly){
			throw new IllegalArgumentException("-autoAdaptLambda requires -adaptEmissionOnly");
		}
		if(autoAdaptLambda && autoAdaptLambdaN0 <= 0){
			throw new IllegalArgumentException("-autoAdaptLambdaN0 must be positive");
		}
		if(autoAdaptLambda && autoAdaptLambdaBeta <= 0.0){
			throw new IllegalArgumentException("-autoAdaptLambdaBeta must be positive");
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

		/** Build a MatrixObj sharing this object's matrixU/matrixM/pi/a/matrixObserved/cpgDistFreq
		 *  but with a different matrix reference. Used by adaptReinitGmm to feed a subsampled
		 *  matrix to GMMLearner without copying unrelated fields. */
		public MatrixObj cloneShallowWithMatrix(ArrayList<Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>>> newMatrix) {
			MatrixObj clone = new MatrixObj(newMatrix, this.matrixU, this.matrixM, this.pi, this.a, this.matrixObserved, this.cpgDistFreq);
			return clone;
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
		final double baseQ;       // per-CpG Phred score from BED col 8; consumed by HMM iff -useBaseQ

		ParsedRow(String readName, String methyStat, String loc, int offset,
				  double methyPrior, double fragLen, double coverage, double distToCenter,
				  double motifScore, double baseQ) {
			this.readName = readName;
			this.methyStat = methyStat;
			this.loc = loc;
			this.offset = offset;
			this.methyPrior = methyPrior;
			this.fragLen = fragLen;
			this.coverage = coverage;
			this.distToCenter = distToCenter;
			this.motifScore = motifScore;
			this.baseQ = baseQ;
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

		if (adaptTransitions) {
			// Regularize transitions and pi toward reference:
			// new = (1 - lambda) * MLE + lambda * ref
			// Both pi and Arij are probability distributions, so the interpolation
			// preserves the sum-to-1 property (both inputs sum to 1).
			for (int r = 0; r <= refHmm.nbCpgDistState(); r++) {
				for (int i = 0; i < refHmm.nbStates(); i++) {
					double mlePi = updatedHmm.getPri(r, i);
					double refPi = refHmm.getPri(r, i);
					updatedHmm.setPri(r, i, (1 - lambda) * mlePi + lambda * refPi);
					for (int j = 0; j < refHmm.nbStates(); j++) {
						double mleA = updatedHmm.getArij(r, i, j);
						double refA = refHmm.getArij(r, i, j);
						updatedHmm.setArij(r, i, j, (1 - lambda) * mleA + lambda * refA);
					}
				}
			}
		} else {
			// Freeze transitions and pi: restore reference values
			for (int r = 0; r <= refHmm.nbCpgDistState(); r++) {
				for (int i = 0; i < refHmm.nbStates(); i++) {
					updatedHmm.setPri(r, i, refHmm.getPri(r, i));
					for (int j = 0; j < refHmm.nbStates(); j++) {
						updatedHmm.setArij(r, i, j, refHmm.getArij(r, i, j));
					}
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

		// Reject input BEDs whose distance-column form differs from the
		// reference model's training-time form.
		validateDistColumnAgainstModel(inputFile, modelFile);
		validateFeatureFlagsAgainstModel(modelFile);

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
		log.info("Loading data for adaptation (miniDataPoints=" + adaptMiniDataPoints + ", minimal mode) ...");
		MatrixObj matrixObj = processMatrixFile(inputFile, /*minimalForAdaptation=*/ true);
		miniDataPoints = savedMiniDataPoints; // restore original
		// matrixU/M/pi/a/matrixObserved/cpgDistFreq are not populated in minimal mode
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
			SummaryStatistics[] fallbackStats = lastComputedStats;
			matrixObj = null;
			matrix = null;
			miniDataPoints = 1;
			decodeOnlyStreaming(inputFile, modelFile, outputFile, cpgIndex, chromOrder, fallbackStats);
			miniDataPoints = savedMiniDataPoints;
			return;
		}

		// Auto-tune -adaptLambda via generalized (Hill) Bayesian shrinkage based on fragment count.
		// Fewer fragments -> noisier MLE -> more regularization toward reference.
		//   lambda = 1 / (1 + (N/N0)^beta)
		// At N=N0: lambda=0.5. At N>>N0: lambda->0. At N<<N0: lambda->1.
		// beta=1.0 recovers the simple Bayesian shrinkage 1 - N/(N+N0).
		// Defaults (beta=0.7, N0=25K) produce:
		//   0.1X (5K)    -> lambda ~ 0.75 (strong reg, noise floor for 2-state 3-feature GMM)
		//   1X   (53K)   -> lambda ~ 0.37 (meaningful adaptation; prior does not dominate)
		//   10X  (530K)  -> lambda ~ 0.10 (light reg; sample signature dominates)
		//   30X  (6.3M)  -> lambda ~ 0.018 (near-unregularized MLE)
		if (autoAdaptLambda) {
			long n = matrix.size();
			double ratio = (double) n / (double) autoAdaptLambdaN0;
			double tunedLambda = 1.0 / (1.0 + Math.pow(ratio, autoAdaptLambdaBeta));
			// Clip to [0.0, 0.95] so there's always at least some adaptation
			if (tunedLambda < 0.0) tunedLambda = 0.0;
			if (tunedLambda > 0.95) tunedLambda = 0.95;
			log.info("Auto-adapt lambda: N=" + n + " qualifying fragments, N0=" + autoAdaptLambdaN0 +
					 ", beta=" + String.format("%.3f", autoAdaptLambdaBeta) +
					 " -> lambda = 1 / (1 + (" + n + "/" + autoAdaptLambdaN0 + ")^" +
					 String.format("%.3f", autoAdaptLambdaBeta) + ") = " +
					 String.format("%.4f", tunedLambda) + " (was " + adaptLambda + ")");
			this.adaptLambda = tunedLambda;
		}

		// Phase 3: Constrained Baum-Welch adaptation
		log.info("Phase 3: Baum-Welch adaptation (lambda=" + String.format("%.4f", adaptLambda) +
				 ", maxIter=" + adaptMaxIter +
				 ", transitions=" + (adaptTransitions ? "ADAPTED" : "FROZEN") + ") ...");

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

		// Keep a pristine clone of the re-centered reference BEFORE any GMM
		// re-init, so the delta/auto-tune metrics compare against the true
		// reference rather than the GMM-init version.
		BayesianNhmmV5<ObservationVector> refHmmPristine;
		try {
			refHmmPristine = refHmmForAdapt.clone();
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException("Failed to clone re-centered reference HMM", e);
		}

		// Optional: re-initialize emissions via GMM on the target data.
		// This escapes the reference-model's local optimum in EM. Transitions
		// and pi stay frozen to the reference model; only the emission GMMs
		// are replaced with target-data-derived initial estimates.
		if (adaptReinitGmm) {
			// Subsample the matrix before GMM init: ClustersGMM allocates a
			// per-ObservationVector Hashtable entry (~100 bytes each), which
			// blows up at ~60M observations. A random sample of ~100K fragments
			// (~1M observations) is plenty to capture the emission distribution
			// for a 2-state 2-3-feature GMM.
			final int GMM_MAX_FRAGMENTS = 100_000;
			List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> gmmMatrix;
			MatrixObj gmmMatrixObj;
			if (matrix.size() > GMM_MAX_FRAGMENTS) {
				log.info("Subsampling " + GMM_MAX_FRAGMENTS + " of " + matrix.size() +
						 " fragments for GMM initialization ...");
				java.util.Random rng = seed >= 0 ? new java.util.Random(seed) : new java.util.Random();
				// Reservoir sampling: pick GMM_MAX_FRAGMENTS uniformly at random.
				// Outlier clipping is handled inside buildInitNhmmByGMM.
				List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>>> sampled =
					new ArrayList<>(GMM_MAX_FRAGMENTS);
				for (int i = 0; i < matrix.size(); i++) {
					if (i < GMM_MAX_FRAGMENTS) {
						sampled.add(matrix.get(i));
					} else {
						int j = rng.nextInt(i + 1);
						if (j < GMM_MAX_FRAGMENTS) {
							sampled.set(j, matrix.get(i));
						}
					}
				}
				gmmMatrix = sampled;
				// Build a shallow matrixObj wrapping the sampled matrix for GMMLearner
				ArrayList<org.apache.commons.lang3.tuple.Triple<HashMap<Integer, Pair<Integer, Double>>, ArrayList<ObservationVector>, ArrayList<String>>> tripleMatrix =
					new ArrayList<>(sampled.size());
				for (Pair<HashMap<Integer, Pair<Integer, Double>>, List<ObservationVector>> p : sampled) {
					tripleMatrix.add(org.apache.commons.lang3.tuple.Triple.of(
						p.getFirst(), new ArrayList<>(p.getSecond()), new ArrayList<String>()));
				}
				gmmMatrixObj = matrixObj.cloneShallowWithMatrix(tripleMatrix);
			} else {
				gmmMatrix = matrix;
				gmmMatrixObj = matrixObj;
			}
			log.info("Re-initializing emission GMM on target data (transitions/pi stay frozen) ...");
			BayesianNhmmV5<ObservationVector> gmmInitHmm = buildInitNhmmByGMM(gmmMatrixObj, gmmMatrix);

			// Validate GMM output — fall back to re-centered reference if NaN.
			boolean gmmHasNaN = false;
			for (int s = 0; s < gmmInitHmm.nbStates() && !gmmHasNaN; s++) {
				double[] m = ((OpdfMultiMixtureGaussian) gmmInitHmm.getOpdf(s)).mean();
				for (double v : m) {
					if (Double.isNaN(v) || Double.isInfinite(v)) { gmmHasNaN = true; break; }
				}
			}
			if (gmmHasNaN) {
				log.warn("GMM re-init produced NaN emissions (likely from extreme feature outliers); " +
						 "keeping re-centered reference emissions as initialization");
			} else {
				// Copy target-data GMM emissions into refHmmForAdapt (replacing
				// reference emissions). Transitions/pi in refHmmForAdapt stay intact.
				for (int s = 0; s < refHmmForAdapt.nbStates(); s++) {
					refHmmForAdapt.setOpdf(s, gmmInitHmm.getOpdf(s));
				}
				logEmissionParams("REFERENCE MODEL (after GMM re-init on target)", refHmmForAdapt);
			}
			// Free subsampled structures
			gmmMatrix = null;
			gmmMatrixObj = null;
			gmmInitHmm = null;
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

			// Regularize toward the PRISTINE re-centered reference, NOT
			// refHmmForAdapt (which has its emissions replaced with GMM-init
			// when -adaptReinitGmm is set). Using GMM-init as the reg target
			// defeats lambda's purpose: BW would just stay near GMM-init
			// regardless of lambda since GMM-init is already target-data-fit.
			adaptedHmm = adaptIteration(adaptedHmm, refHmmPristine, matrix, adaptLambda);

			distance = klc.distance(prevHmm, adaptedHmm, true);
			log.info("Adaptation iteration " + (iter + 1) + ": KL distance = " + distance);
			logEmissionParams("After iteration " + (iter + 1), adaptedHmm);

			if (Double.isNaN(distance)) {
				log.warn("KL distance is NaN at iteration " + (iter + 1) +
						 "; falling back to pristine re-centered reference model.");
				adaptedHmm = refHmmPristine;
				break;
			}
			if (Math.abs(distance) < tol) {
				log.info("Adaptation converged at iteration " + (iter + 1));
				break;
			}
		}

		// Use the pristine re-centered reference (not the GMM-init-overwritten
		// refHmmForAdapt) so the delta reflects "how far from reference did
		// BW end up" rather than "how far did BW move from the GMM init".
		logEmissionDelta("EMISSION DELTA (adapted - re-centered reference)", refHmmPristine, adaptedHmm);

		this.methylatedState = adaptedHmm.getMethyState(lowCoverage);
		adaptedHmm.setMaxCpgNum(cpgNumClip < 0 ? maxCpgs : cpgNumClip);
		adaptedHmm.setMinCpgNum(1);

		// Auto-tune bayesianFactor based on HOW MUCH the HMM emissions moved
		// from the (re-centered) reference.
		//
		// When -autoAdaptLambda is ON and lambda is high (low-cov), the
		// decoding BW heavily regularizes toward the reference, so the
		// adapted emissions look ~= reference regardless of sample type.
		// This artificially hides the sample signature from the shift metric.
		//
		// Solution: run a short "signature probe" BW with a FIXED low
		// lambda (autoTuneProbeLambda, default 0.0 = unregularized MLE)
		// to preserve the sample signature. Measure shift from the probe,
		// not from the heavily-regularized decoding model.
		if (autoTuneBayesianFactor) {
			BayesianNhmmV5<ObservationVector> shiftSourceHmm;
			// Run unregularized probe BW whenever decoding BW is regularized
			// (adaptLambda > 0) to measure the true sample signature. The
			// probe's lambda (default 0) must be strictly less than adaptLambda
			// to be informative. Set -autoTuneProbeLambda to a negative value
			// to disable and reuse the adapted model's shift directly.
			if (autoTuneProbeLambda >= 0 && autoTuneProbeLambda < adaptLambda) {
				log.info("Running signature probe BW (lambda=" + autoTuneProbeLambda +
						 ", maxIter=" + autoTuneProbeMaxIter +
						 ") to decouple shift signal from decoding regularization (lambda=" +
						 String.format("%.4f", adaptLambda) + ") ...");
				BayesianNhmmV5<ObservationVector> probeHmm;
				try {
					probeHmm = refHmmForAdapt.clone();
				} catch (CloneNotSupportedException e) {
					throw new RuntimeException("Failed to clone for probe BW", e);
				}
				double probeDistance = Double.MAX_VALUE;
				for (int iter = 0; iter < autoTuneProbeMaxIter; iter++) {
					BayesianNhmmV5<ObservationVector> prev;
					try {
						prev = probeHmm.clone();
					} catch (CloneNotSupportedException e) {
						throw new RuntimeException("Failed to clone probe HMM", e);
					}
					probeHmm = adaptIteration(probeHmm, refHmmPristine, matrix, autoTuneProbeLambda);
					probeDistance = klc.distance(prev, probeHmm, true);
					if (Double.isNaN(probeDistance)) {
						log.warn("Probe KL is NaN; aborting probe, will use adapted-model shift instead.");
						probeHmm = adaptedHmm;
						break;
					}
					if (Math.abs(probeDistance) < tol) {
						log.info("Probe BW converged at iteration " + (iter + 1));
						break;
					}
				}
				logEmissionParams("SIGNATURE PROBE (lambda=" + autoTuneProbeLambda + ")", probeHmm);
				shiftSourceHmm = probeHmm;
			} else {
				shiftSourceHmm = adaptedHmm;
			}
			// Compare against the pristine re-centered reference.
			double tunedFactor = computeAutoTunedBayesianFactor(refHmmPristine, shiftSourceHmm, klc, autoTuneMidpoint, autoTuneSteepness);
			log.info("Auto-tuned bayesianFactor: " + String.format("%.4f", tunedFactor) +
					 " (was " + bayesianFactor + ")");
			this.bayesianFactor = tunedFactor;
			adaptedHmm.setBayesianFactor(tunedFactor);
		}

		System.out.println("Adapted HMM:\n" + adaptedHmm);

		// Write adapted model to a temp file for decodeHmm (which reloads from disk).
		// The original reference model file is NEVER modified.
		File adaptedModelTmp = File.createTempFile("adapted_model_", ".ser");
		adaptedModelTmp.deleteOnExit();
		ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(adaptedModelTmp));
		oos.writeObject(adaptedHmm);
		oos.close();
		log.info("Adapted model written to temp file: " + adaptedModelTmp.getAbsolutePath());

		// Propagate the reference model's distance-column metadata + feature
		// flags to the adapted temp model so the downstream decodeOnlyStreaming
		// call (which validates) doesn't reject it. The adapted model uses
		// the same feature schema as the reference; reuse the reference's
		// recorded flag values (or fall back to the current invocation's
		// flags if the reference predates the v0.64 sidecar format).
		Map<String, String> refMeta = readModelMeta(modelFile);
		String refDistCol = refMeta.get("dist_column");
		String adaptedDistCol = (refDistCol != null)
				? refDistCol
				: readDistColumnNameFromHeader(inputFile);
		Boolean adaptedUseMotif = refMeta.containsKey("use_motif")
				? Boolean.parseBoolean(refMeta.get("use_motif")) : useEndMotif;
		Boolean adaptedUseBaseQ = refMeta.containsKey("use_baseq")
				? Boolean.parseBoolean(refMeta.get("use_baseq")) : useBaseQ;
		Boolean adaptedLowCov   = refMeta.containsKey("low_coverage")
				? Boolean.parseBoolean(refMeta.get("low_coverage")) : lowCoverage;
		Integer adaptedFeatures = refMeta.containsKey("features")
				? Integer.parseInt(refMeta.get("features")) : features;
		writeModelMeta(adaptedModelTmp.getAbsolutePath(), adaptedDistCol,
				adaptedUseMotif, adaptedUseBaseQ, adaptedLowCov, adaptedFeatures);
		File adaptedTmpMeta = new File(modelMetaPath(adaptedModelTmp.getAbsolutePath()));
		adaptedTmpMeta.deleteOnExit();

		// Optionally save adapted model to a user-specified path
		if (saveAdaptedModel != null) {
			java.nio.file.Files.copy(adaptedModelTmp.toPath(),
				new File(saveAdaptedModel).toPath(),
				java.nio.file.StandardCopyOption.REPLACE_EXISTING);
			writeModelMeta(saveAdaptedModel, adaptedDistCol,
					adaptedUseMotif, adaptedUseBaseQ, adaptedLowCov, adaptedFeatures);
			log.info("Adapted model saved to: " + saveAdaptedModel +
					" (with metadata sidecar)");
		}

		// Free adaptation data before decode to reduce memory
		matrix = null;
		matrixObj = null;

		// Phase 4: Viterbi decode ALL fragments (miniDataPoints=1) using the
		// adapted model via streaming decode. This outputs every fragment
		// regardless of the -miniDataPoints used for adaptation training.
		// Pass the already-computed stats to skip a redundant file scan.
		log.info("Phase 4: Streaming decode with adapted model (all fragments, reusing cached stats) ...");
		int decodeMiniDataPoints = miniDataPoints;
		miniDataPoints = 1;  // decode everything
		SummaryStatistics[] cachedStats = lastComputedStats;
		lastComputedStats = null; // free for GC
		decodeOnlyStreaming(inputFile, adaptedModelTmp.getAbsolutePath(), outputFile, cpgIndex, chromOrder, cachedStats);
		miniDataPoints = decodeMiniDataPoints;  // restore
		adaptedModelTmp.delete();
	}

	/**
	 * Log per-state emission GMM parameters (mean, variance, mixture proportions)
	 * for adaptation diagnostics.
	 */
	private void logEmissionParams(String label, BayesianNhmmV5<ObservationVector> hmm) {
		// Build feature-name list dynamically based on the active flags so
		// that -useBaseQ (and any future flag) propagates to the diagnostic
		// logs without manual table edits.
		java.util.ArrayList<String> names = new java.util.ArrayList<>(5);
		names.add("FragLen");
		if (!lowCoverage) names.add("Coverage");
		names.add("DistToCenter");
		if (useEndMotif) names.add("MotifScore");
		if (useBaseQ) names.add("baseQ");
		final String[] featureNames = names.toArray(new String[0]);

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
	/**
	 * Auto-tune bayesianFactor using a sigmoid on the emission shift between
	 * the (re-centered) reference HMM and the adapted HMM. Shift is the max
	 * over states of the L2 distance between the per-state emission means
	 * in z-score feature space.
	 *
	 * Mapping:
	 *   sigmoidInput = (shift - midpoint) * steepness
	 *   sigmoidOutput = 1 / (1 + exp(-sigmoidInput))     -- in [0, 1]
	 *   bayesianFactor = clip(0.9 - 0.85 * sigmoidOutput, 0.05, 0.9)
	 *
	 * Saturates near 0.9 for shift << midpoint (trust prior) and near 0.05
	 * for shift >> midpoint (distrust prior), with a sharp transition
	 * centered on `midpoint`. Steepness controls transition width.
	 *
	 * Defaults (midpoint=0.55, steepness=15) calibrated from observed data:
	 *   healthy ~ shift 0.37  -> bayesianFactor ~ 0.85
	 *   disease ~ shift 0.72  -> bayesianFactor ~ 0.11
	 *
	 * Secondary diagnostic: KL divergence between reference and adapted
	 * HMMs. Logged for interpretation; not used in tuning.
	 */
	private double computeAutoTunedBayesianFactor(
			BayesianNhmmV5<ObservationVector> refHmm,
			BayesianNhmmV5<ObservationVector> adaptedHmm,
			KullbackLeiblerDistanceBayesianNhmmV5Calculator klc,
			double midpoint, double steepness) {
		int nStates = refHmm.nbStates();
		double maxShift = 0.0;
		double[] perStateShift = new double[nStates];
		for (int s = 0; s < nStates; s++) {
			double[] refMean = ((OpdfMultiMixtureGaussian) refHmm.getOpdf(s)).mean();
			double[] adpMean = ((OpdfMultiMixtureGaussian) adaptedHmm.getOpdf(s)).mean();
			double sumSq = 0;
			int d = Math.min(refMean.length, adpMean.length);
			for (int k = 0; k < d; k++) {
				double diff = adpMean[k] - refMean[k];
				sumSq += diff * diff;
			}
			double shift = Math.sqrt(sumSq);
			perStateShift[s] = shift;
			if (shift > maxShift) maxShift = shift;
		}

		// Also compute KL divergence as a supplementary diagnostic
		double klDist = Double.NaN;
		try {
			klDist = klc.distance(refHmm, adaptedHmm, true);
		} catch (Exception e) {
			log.info("KL divergence computation failed (non-fatal): " + e.getMessage());
		}

		// Sigmoid mapping: sharp transition around midpoint.
		double sigmoidInput = (maxShift - midpoint) * steepness;
		double sigmoidOutput;
		if (sigmoidInput > 20) sigmoidOutput = 1.0;
		else if (sigmoidInput < -20) sigmoidOutput = 0.0;
		else sigmoidOutput = 1.0 / (1.0 + Math.exp(-sigmoidInput));
		double tuned = 0.9 - 0.85 * sigmoidOutput;
		if (tuned < 0.05) tuned = 0.05;
		if (tuned > 0.9) tuned = 0.9;

		StringBuilder shiftStr = new StringBuilder();
		shiftStr.append("per-state L2 shifts: [");
		for (int s = 0; s < nStates; s++) {
			if (s > 0) shiftStr.append(", ");
			shiftStr.append(String.format("state%d=%.4f", s, perStateShift[s]));
		}
		shiftStr.append("]");
		log.info("Auto-tune emission shift (adapted vs re-centered reference): " + shiftStr.toString());
		log.info("Auto-tune max shift = " + String.format("%.4f", maxShift) +
				 " sd, KL divergence = " + String.format("%.6f", klDist) +
				 ", midpoint = " + midpoint + ", steepness = " + steepness +
				 " -> sigmoid((" + String.format("%.4f", maxShift) + " - " + midpoint + ") * " + steepness +
				 ") = " + String.format("%.4f", sigmoidOutput) +
				 " -> bayesianFactor = " + String.format("%.4f", tuned));

		if (maxShift < midpoint - 0.15) {
			log.info("Adapted emissions close to reference (shift < midpoint - 0.15): sample appears " +
					 "SIMILAR to reference tissue type. bayesianFactor saturates near 0.9 (high prior trust).");
		} else if (maxShift > midpoint + 0.15) {
			log.warn("Adapted emissions DIVERGED substantially from reference (shift > midpoint + 0.15): " +
					 "sample differs from reference tissue type. bayesianFactor saturates near 0.05 " +
					 "(low prior trust). Consider using a sample-appropriate -valueWigs track.");
		}
		return tuned;
	}

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
				// Map feature-dimension index in the observation vector to
				// its slot in the stats[] array. Uses the same flag matrix
				// as buildObservationVector(); see that helper for layout
				// reasoning. Order in the obs vector:
				//   FragLen, [Coverage], DistToCenter, [MotifScore], [baseQ]
				// where bracketed dims are present only when the
				// corresponding flag is on. stats[] is always
				//   [0]=FragLen, [1]=Coverage, [2]=DistToCenter,
				//   [3]=MotifScore (if motif), [bqIdx]=baseQ.
				int statsIdx = obsDimToStatsIdx(d);

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
		// Build feature-name list dynamically (same logic as logEmissionParams).
		java.util.ArrayList<String> names = new java.util.ArrayList<>(5);
		names.add("FragLen");
		if (!lowCoverage) names.add("Coverage");
		names.add("DistToCenter");
		if (useEndMotif) names.add("MotifScore");
		if (useBaseQ) names.add("baseQ");
		final String[] featureNames = names.toArray(new String[0]);

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
