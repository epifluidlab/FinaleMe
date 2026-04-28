/**
 * CpgFeatureMatrixBuilder.java
 * Feb 27, 2016
 * 5:28:37 PM
 * yaping    lyping1986@gmail.com
 */
package edu.northwestern.epifluidlab.finaleme.utils;


import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.IntervalTree.Node;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.TabixReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.biojava.nbio.genome.parsers.twobit.TwoBitParser;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.unc.genomics.io.BigWigFileReader;

/**
 *
 */
public class CpgFeatureMatrixBuilder extends AbstractCpgMultiMetricsStats {

	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;

	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-maxFragLen",usage="maximum fragment length allowed to check. Default: 500")
	public int maxFragLen = 500;

	@Option(name="-maxDistToFragEnd",usage="maximum distant to the end of the fragment allowed to check. in order to be copnsistent with training model. Default: 250")
	public int maxDistToFragEnd = 250;

	@Option(name="-useSignedDistCenter",usage="write a STRAND-AWARE SIGNED distance from fragment center in the Dist_frag_end column instead of the legacy unsigned distance-to-nearest-end. Value = (cpgOffset_from_5'end) - fragLen/2.0, sign-corrected for strand: negative = CpG is in the 5' half of its fragment, positive = 3' half, magnitude = bp from center. The legacy unsigned form is symmetric and discards 5'/3' position information; the signed form preserves it and is strictly more informative for the HMM. Header column name becomes 'Signed_Dist_Center' so downstream tools can detect which form is in use. NOTE: a model trained on one form CANNOT be reused with the other (the feature distribution differs); set this consistently across train and decode runs. Default: false")
	public boolean useSignedDistCenter = false;

	@Option(name="-totalReadsInBam",usage="total number of reads (or fragments, depending on -coverageReadLevel) used to normalize coverage column. default estimate from bam file by program. Default: -1")
	public long totalReadsInBam = -1;

	@Option(name="-coverageReadLevel",usage="legacy mode: count Norm_Frag_cov as raw_reads_at_CpG / total_filtered_reads (each PE read counted separately on both numerator and denominator). The new default (false) uses fragments for both: numerator = unique fragments at CpG (deduplicated by readName), denominator = total filtered fragments. Fragment-level matches tabix fragment mode and gives a Norm_Frag_cov scale that's independent of paired-end vs single-end input. Set to true ONLY when re-using a model trained before v0.64 (which was trained at the read-level scale). Default: false")
	public boolean coverageReadLevel = false;
	
	@Option(name="-maxCov",usage="maximum coverage allowed to check. Default: 250")
	public int maxCov = 250;
	
	@Option(name="-kmerLen",usage="the K-mer length to check. Default: 0")
	public int kmerLen = 0;

	@Option(name="-kmerString",usage="the fiel contain selected K-mer to check. Otherwise, use -kmerLen to automately generate all the k-mer. Default: null")
	public String kmerString = null;

	@Option(name="-kmerExt",usage="the +/- region in reference genome to check the k-mer frequency. default is +/- 100bp around CpGs. Default: 100")
	public int kmerExt = 100;
	
	@Option(name="-excludeRegions",usage="bed files indicated excluded regions. -excludeRegions trackFileName. Default: null")
	public ArrayList<String> excludeRegions = null;

	@Option(name="-overlapRegions",usage="bed files to check if regions are overlapped. -overlapRegions trasckName:trackFileName. Default: null")
	public ArrayList<String> overlapRegions = null;

	@Option(name="-distantRegions",usage="bed files to check the distance to these regions. like the distance to TSS. -distantRegions trasckName:trackFileName. Default: null")
	public ArrayList<String> distantRegions = null;

	@Option(name="-valueWigs",usage="bigwig files to check the value in these regions. like the recombination rate of some region.extRegion indicate the +/- regionBp from CpG, -valueWigs trasckName:extRegion:trackFileName. Default: null")
	public ArrayList<String> valueWigs = null;

	@Option(name="-valueBeds",usage="tabixed bed.gz files to check the value in these regions. like the recombination rate of some region. -valueBeds trasckName:extRegion:trackFileName. Default: null")
	public ArrayList<String> valueBeds = null;

	@Option(name="-wgsMode",usage="used for WGS but not bisulfite space. . Default: false")
	public boolean wgsMode = false;

	@Option(name="-skipSecondEnd",usage="skip the 2nd end for the statistics. Default: false")
	public boolean skipSecondEnd = false;

	@Option(name="-includeCpgDist",usage="include the distance between CpGs. Default: false")
	public boolean includeCpgDist = false;

	@Option(name="-stringentPaired",usage="Only use paired end reads that faced to each other. Default: false")
	public boolean stringentPaired = false;
	
	@Option(name="-useFragBaseKmer",usage="add k-mer in fragment. Default: false")
	public boolean useFragBaseKmer = false;
	
	@Option(name="-useStrandSpecificFragBase",usage="use k-mer generated by aware the strand of fragment status. Default: false")
	public boolean useStrandSpecificFragBase = false;

	
	@Option(name="-useNoChrPrefixBam",usage="use bam file with GRch37 instead of hg19 coordinate. Default: false")
	public boolean useNoChrPrefixBam = false;

	@Option(name="-fragmentInputTabix",usage="interpret the 4th argument as bgzipped/tabix-indexed fragment BED/TSV file instead of BAM/CRAM. default: auto-detect by extension")
	public boolean fragmentInputTabix = false;

	@Option(name="-fragStrandColumn",usage="1-based column index for fragment strand in tabix BED/TSV input. 0 means auto-detect from col4/col6. default: 0")
	public int fragStrandColumn = 0;

	@Option(name="-fragNameColumn",usage="1-based column index for fragment name in tabix BED/TSV input. 0 = synthesize a per-row coordinate-based name (chr_start_end_strand_line; guaranteed unique per fragment). Set this to a positive index ONLY when the input has genuine unique read names (e.g., a BAM-derived fragment file with original SAM read names) -- the column value must NOT be a BED placeholder ('.', '*', empty) or every fragment collapses to one 'read' and breaks HMM training with a singular covariance matrix. default: 0")
	public int fragNameColumn = 0;

	@Option(name="-fragMethyColumn",usage="1-based column index for methylation state (m/u) in tabix BED/TSV input. 0 means infer from valueWig/default. default: 0")
	public int fragMethyColumn = 0;

	@Option(name="-fragBaseQ",usage="synthetic baseQ used for tabix fragment input when no per-base quality exists. default: 60")
	public int fragBaseQ = 60;

	@Option(name="-defaultMethyStat",usage="default methylation state ('m' or 'u') for the methy_stat output column in tabix fragment mode when -fragMethyColumn is unset and no inline 'm'/'u' token is detected in extra columns. methy_stat is NOT consumed by HMM training or decoding -- the HMM is unsupervised over the feature vector -- so this only affects the output column used for AUC/QC reporting. The default 'm' matches the behavior of BAM input under -wgsMode (where every covered CpG is labeled 'm' because the read base is the unconverted reference). default: m")
	public String defaultMethyStat = "m";

	@Option(name="-useEndMotif",usage="add 5' end 4-mer motif score as a feature column in the output. Default: false")
	public boolean useEndMotif = false;

	@Option(name="-noCoverage",usage="write NaN for the Norm_Frag_cov column (coverage-free mode). Default: false")
	public boolean noCoverage = false;

	@Option(name="-saveMotifLookup",usage="save 5' end motif score lookup table (256 4-mer -> methylation rate) to TSV file during training. Default: null")
	public String saveMotifLookup = null;

	@Option(name="-loadMotifLookup",usage="load 5' end motif score lookup table from TSV file for decode mode. Default: null")
	public String loadMotifLookup = null;

	@Option(name="-t",usage="number of threads for parallel 5Mb bin processing. Use >0 to set explicitly; default uses all available cores.")
	public int threads = -1;

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "CpgFeatureMatrixBuilder [opts] hg19.2bit cpg_list.bed all_cpg.bed wgs.bam|fragments.tsv.gz cpg_detail.txt.gz";
	
	private static final Logger log = LoggerFactory.getLogger(CpgFeatureMatrixBuilder.class);

	private static long startTime = -1;
	private static long points = 0;

	// 5' end motif scoring: 256 4-mer -> methylation rate mapping
	private java.util.concurrent.ConcurrentHashMap<String, java.util.concurrent.atomic.AtomicLongArray> motifCounts = null;
	private HashMap<String, Double> motifScoreMap = null;


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		CpgFeatureMatrixBuilder cmms = new CpgFeatureMatrixBuilder();
		//BasicConfigurator.configure();
		cmms.doMain(args);
	}
	
	@SuppressWarnings("resource")
	public void doMain(String[] args)
			throws Exception {

					CmdLineParser parser = new CmdLineParser(this);
					//parser.setUsageWidth(80);
					try
					{
						if(help || args.length < 5) throw new CmdLineException(parser, USAGE, new Throwable());
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

					// Validate motif lookup options
					if (saveMotifLookup != null && !useEndMotif) {
						throw new CmdLineException(parser, "-saveMotifLookup requires -useEndMotif", new Throwable());
					}
					if (loadMotifLookup != null && !useEndMotif) {
						throw new CmdLineException(parser, "-loadMotifLookup requires -useEndMotif", new Throwable());
					}
					if (saveMotifLookup != null && loadMotifLookup != null) {
						throw new CmdLineException(parser, "-saveMotifLookup and -loadMotifLookup are mutually exclusive", new Throwable());
					}
					// -saveMotifLookup needs ground-truth methylation per CpG to compute
					// methylated/total per 4-mer. The only mode that supplies this is
					// bisulfite-converted BAM input (default BAM, -wgsMode=false). Under
					// -wgsMode every covered CpG is labeled 'm' (the read base is the
					// unconverted reference C/G), so -saveMotifLookup would silently
					// produce a degenerate lookup with score=1.0 for every 4-mer.
					// Error out early instead of producing a useless file. The tabix
					// equivalent of this check happens after useTabixFragmentInput is
					// resolved below (line ~268).
					if (saveMotifLookup != null && wgsMode) {
						throw new CmdLineException(parser,
								"-saveMotifLookup requires bisulfite-converted BAM input " +
								"(-wgsMode=false). Under -wgsMode every covered CpG is " +
								"labeled methylated, so the lookup table would be " +
								"degenerate (score=1.0 for every 4-mer). Train the motif " +
								"lookup on a paired bisulfite/WGBS BAM, then re-use it " +
								"in -wgsMode/tabix runs via -loadMotifLookup.",
								new Throwable());
					}

					// Initialize motif data structures
					if (saveMotifLookup != null) {
						motifCounts = new java.util.concurrent.ConcurrentHashMap<>();
					}
					if (loadMotifLookup != null) {
						motifScoreMap = loadMotifLookupFile(loadMotifLookup);
						log.info("Loaded motif lookup with " + motifScoreMap.size() + " entries from " + loadMotifLookup);
					}

					//read input bed file, for each row,
					//String intervalFile = arguments.get(0);
					String refFile = arguments.get(0);
					String cpgListFile = arguments.get(1);
					String allCpgFile = arguments.get(2);
					String wgsBamFile = arguments.get(3);
					String detailFile = arguments.get(4);
					boolean useTabixFragmentInput = isTabixFragmentInput(wgsBamFile);
					log.info("Input fragment source mode: " + (useTabixFragmentInput ? "tabix BED/TSV" : "BAM/CRAM"));
					if(useTabixFragmentInput){
						validateTabixIndexExists(wgsBamFile);
					}
					// Companion check to the -saveMotifLookup + -wgsMode validation
					// above: tabix fragment input also lacks per-CpG bisulfite
					// information, so methy_stat is constant ('m' under
					// -defaultMethyStat) and the resulting lookup would be
					// degenerate.
					if (saveMotifLookup != null && useTabixFragmentInput) {
						throw new CmdLineException(parser,
								"-saveMotifLookup is not supported with tabix fragment " +
								"input. Tabix fragment files lack per-CpG bisulfite " +
								"information, so methy_stat is constant and the resulting " +
								"motif lookup would be degenerate (score=1.0 for every " +
								"4-mer). Train the motif lookup on a paired WGBS BAM, " +
								"then re-use it on tabix/WGS runs via -loadMotifLookup.",
								new Throwable());
					}

					// In tabix fragment mode we cannot train a motif lookup
					// (the check above), so -useEndMotif requires a pre-trained
					// lookup. Without -loadMotifLookup, getMotifScore() returns
					// 0.5 for every fragment -> the motif_score column is a
					// constant -> downstream FinaleMe HMM training hits a
					// singular covariance matrix (sd=0 feature). Catch this
					// at startup with a clear message rather than producing
					// a feature matrix that silently breaks training.
					if (useEndMotif && useTabixFragmentInput && loadMotifLookup == null) {
						throw new CmdLineException(parser,
								"-useEndMotif on tabix fragment input requires " +
								"-loadMotifLookup <pretrained.tsv>. Without it the " +
								"motif_score column is constant (0.5 for every " +
								"fragment) because the motif lookup cannot be " +
								"trained in tabix mode (methy_stat is degenerate). " +
								"That makes the motif_score feature degenerate and " +
								"causes 'matrix is singular' errors in FinaleMe HMM " +
								"training. Either drop -useEndMotif, or train the " +
								"lookup on a paired bisulfite/WGBS BAM (default BAM " +
								"mode + -saveMotifLookup) and pass the resulting " +
								"TSV here via -loadMotifLookup.",
								new Throwable());
					}

					initiate();			
					
					//initiate different kinds of reader
					//reference genome
					TwoBitParser refParser = new TwoBitParser(new File(refFile));
					
					//String[] names = p.getSequenceNames();
					//for(int i=0;i<names.length;i++) {
					//  p.setCurrentSequence(names[i]);
					//  p.printFastaSequence();
					//  p.close();
					//}
					//loading exlusion interval file

					//load interval files
					log.info("Processing interval file ... ");
					SamReader wgsReader = null;
					if(!useTabixFragmentInput){
						wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(wgsBamFile));
					}
					//SAMSequenceDictionary dictSeq = SAMSequenceDictionaryExtractor.extractDictionary(new File(wgsBamFile));
					//GenomeLocParser glpSeq = new GenomeLocParser(dictSeq);
					
						HashMap<String,IntervalTree<Integer>> ignoreLocCollections = loadExcludeRegions(excludeRegions, log);
					
					
							HashMap<String,IntervalTree<String>> allCpgLocCollections = new HashMap<String,IntervalTree<String>>();
							if(includeCpgDist){
								allCpgLocCollections = loadStrandedIntervals(allCpgFile, log, "Loading all CpG intervals ... ");
							}
						
						
						
						
					HashMap<String, HashMap<String,IntervalTree<Integer>>> overlapLocStringCollections = null;
					LinkedHashSet<String> overlapLocString = new LinkedHashSet<String>();
					if(overlapRegions!=null && !overlapRegions.isEmpty()){
						log.info("Overalpped intervals ... ");
						overlapLocStringCollections = new HashMap<String, HashMap<String,IntervalTree<Integer>>>();
						for(String overlapRegionString : overlapRegions){
							
							String[] splitStrings = overlapRegionString.split(":");
							if(splitStrings.length < 2){
								throw new IllegalArgumentException("need to provide trackname:trackFile for overlapRegions");
							}
							HashMap<String,IntervalTree<Integer>> overlapLocCollections =  new HashMap<String,IntervalTree<Integer>>();
							String overlapRegionName = splitStrings[0];
							String overlapRegion = splitStrings[1];
							overlapLocString.add(overlapRegionName);
							BufferedReader br = new BufferedReader(new FileReader(overlapRegion));
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
								if(overlapLocCollections.containsKey(chr)){
									tree = overlapLocCollections.get(chr);
								}else{
									tree = new IntervalTree<Integer>();
								}
								tree.put(start, end, 1);
								overlapLocCollections.put(chr, tree);
								
							}
							br.close();
							overlapLocStringCollections.put(overlapRegionName, overlapLocCollections);
						
						}
					}
					
					HashMap<String, HashMap<String,IntervalTree<String>>> distantLocStringCollections = null;
					LinkedHashSet<String> distantLocString = new LinkedHashSet<String>();
					if(distantRegions!=null && !distantRegions.isEmpty()){
						log.info(" Intervals used to calculate distances... ");
						distantLocStringCollections = new HashMap<String, HashMap<String,IntervalTree<String>>>();
						for(String distantRegionString : distantRegions){
							
							String[] splitStrings = distantRegionString.split(":");
							if(splitStrings.length < 2){
								throw new IllegalArgumentException("need to provide trackname:trackFile for distantRegions");
							}
							HashMap<String,IntervalTree<String>> distantLocCollections =  new HashMap<String,IntervalTree<String>>();
							String distantRegionName = splitStrings[0];
							String distantRegion = splitStrings[1];
							distantLocString.add(distantRegionName);
							
							BufferedReader br = new BufferedReader(new FileReader(distantRegion));
							String line;
							while( (line = br.readLine()) != null){
								if(line.startsWith("#"))
									continue;
								String[] splitLines = line.split("\t");
								if(splitLines.length<3){
									continue;
								}else{
									String chr = splitLines[0];
									int start = Integer.parseInt(splitLines[1]);
									int end = Integer.parseInt(splitLines[2]);
									IntervalTree<String> tree;
									if(distantLocCollections.containsKey(chr)){
										tree = distantLocCollections.get(chr);
									}else{
										tree = new IntervalTree<String>();
									}
									String strand = ".";
									if(splitLines.length >= 6){
										if(splitLines[5].equalsIgnoreCase("-")){
											strand = "-";
										}else if(splitLines[5].equalsIgnoreCase("+")){
											strand = "+";
										}
										//strand = splitLines[5].equalsIgnoreCase("-") ? "-" : "+";
									}
									tree.put(start, end, strand);
									distantLocCollections.put(chr, tree);
								}
								
								
								
								
							}
							br.close();
							distantLocStringCollections.put(distantRegionName, distantLocCollections);
						
						}
					}
					
					
					HashMap<String, Pair<Integer, TabixFeatureReader<BEDFeature, ?>>> valueBedReaders = null;
					LinkedHashSet<String> valueBedLocString = new LinkedHashSet<String>();
					if(valueBeds != null){
						log.info("Loading value interval bed file ... ");
						valueBedReaders =  new HashMap<String, Pair<Integer, TabixFeatureReader<BEDFeature, ?>>>();
						for(String valueBedString : valueBeds){
							String[] splitStrings = valueBedString.split(":");
							if(splitStrings.length < 3){
								throw new IllegalArgumentException("need to provide trackname:extRegion:trackFile for valueBeds");
							}
							String valueBedName = splitStrings[0];
							int valueBedExt = Integer.parseInt(splitStrings[1]);
							String valueRegion = splitStrings[2];
							valueBedLocString.add(valueBedName);
							valueBedReaders.put(valueBedName,new Pair<Integer, TabixFeatureReader<BEDFeature, ?>>(valueBedExt,new TabixFeatureReader(valueRegion, new BEDCodec())));
						}
						
					}
					
					HashMap<String, Pair<Integer, BigWigFileReader>> valueWigReaders = null;
					LinkedHashSet<String> valueWigLocString = new LinkedHashSet<String>();
					if(valueWigs != null){
						log.info("Loading value interval big wig file ... ");
						valueWigReaders =  new HashMap<String, Pair<Integer,BigWigFileReader>>();
						for(String valueWigString : valueWigs){
							String[] splitStrings = valueWigString.split(":");
							if(splitStrings.length < 3){
								throw new IllegalArgumentException("need to provide trackname:extRegion:trackFile for valueWigs");
							}
							String valueWigName = splitStrings[0];
							int valueWigExt = Integer.parseInt(splitStrings[1]);
							String valueRegion = splitStrings[2];
							valueWigLocString.add(valueWigName);
							valueWigReaders.put(valueWigName,new Pair<Integer, BigWigFileReader>(valueWigExt,new BigWigFileReader(new File(valueRegion).toPath())));
						}
						
					}
					
					LinkedHashSet<String> kmerCollections = new LinkedHashSet<String>();
					if(kmerString != null){
						log.info("Loading selected K-mer file ... ");	
						BufferedReader br = new BufferedReader(new FileReader(kmerString));
						String line;
						while( (line = br.readLine()) != null){
							if(line.startsWith("#"))
								continue;
							kmerCollections.add(line);
							
							
						}
						br.close();
					}else{
						log.info("Automate generate all k-mer until length " + kmerLen);
						for(int i = 2; i <=  kmerLen; i++){
							for(byte[] kmer : SequenceUtil.generateAllKmers(i)){
								kmerCollections.add(new String(kmer));
								
							}
						}
						
					}
					
					String header = "";
					if(overlapLocString.size()>0){
						for(String key : overlapLocString){
							header = header + "\t" + key;
						}
					}
					if(distantLocString.size()>0){
						for(String key : distantLocString){
							header = header + "\t" + key;
						}
					}
					if(valueBedLocString.size()>0){
						for(String key : valueBedLocString){
							header = header + "\t" + key;
						}
					}
					if(valueWigLocString.size()>0){
						for(String key : valueWigLocString){
							header = header + "\t" + key;
						}
					}
					if(kmerCollections.size()>0){
						if(!useFragBaseKmer){
							for(String key : kmerCollections){
								header = header + "\t" + key;
							}
						}else{
							for(String key : kmerCollections){
								header = header + "\t" + key + "_frag";
							}
						}
						
					}
					
					
						HashMap<String,IntervalTree<String>> cpgCollections =
								loadStrandedIntervals(cpgListFile, log, "Loading CpG interval file ... ");

					
					double readsNumTotal = 0;
					if (noCoverage) {
						// Coverage column is NaN in -noCoverage mode; skip the expensive
						// total-reads count entirely (saves minutes on large BAMs).
						log.info("Skipping total reads count (-noCoverage mode; coverage column will be NaN)");
						readsNumTotal = 1;  // avoid division by zero downstream
					} else if(totalReadsInBam > 0){
						log.info("Get total reads number used for scaling from input option -totalReadsInBam ... ");
						readsNumTotal = totalReadsInBam;
					}else if(useTabixFragmentInput){
						log.info("Get total fragments number used for scaling from tabix fragment file (BGZF sampling) ... ");
						readsNumTotal = estimateTotalFragmentsFromTabixInput(wgsBamFile);
						log.info("Estimated total fragments: " + (long)readsNumTotal +
								" (sampled from BGZF; pass -totalReadsInBam <N> to skip estimation)");
					}else{
						// Strategy 1 (fast, ~2 seconds): Get raw total from BAM index or
						// samtools idxstats, then sample ~1M reads to measure the pass
						// rate under failFlagFilter + stringentPaired, and scale.
						// This gives an accurate filtered count without a full scan.
						log.info("Estimating filtered read count (index + sampling) ...");
						boolean usedIndex = false;

						// Step A: Get raw (unfiltered) total from index
						long rawIndexTotal = 0;
						if(wgsReader.indexing() != null && wgsReader.indexing().hasBrowseableIndex()){
							try {
								htsjdk.samtools.BAMIndex bamIndex = wgsReader.indexing().getIndex();
								int nRefs = wgsReader.getFileHeader().getSequenceDictionary().getSequences().size();
								for(int refIdx = 0; refIdx < nRefs; refIdx++){
									try {
										htsjdk.samtools.BAMIndexMetaData meta = bamIndex.getMetaData(refIdx);
										if(meta != null){
											rawIndexTotal += meta.getAlignedRecordCount();
										}
									} catch(Exception refEx){
										// skip bad references
									}
								}
								if (rawIndexTotal > 0) {
									log.info("BAM index: " + rawIndexTotal + " raw aligned reads");
								}
							} catch(Exception e){
								log.info("BAM index metadata failed: " + e.getMessage());
							}
						}
						// If BAM index metadata returned 0 (absent metadata bin), try samtools idxstats
						if(rawIndexTotal == 0){
							try {
								ProcessBuilder pb = new ProcessBuilder("samtools", "idxstats", wgsBamFile);
								pb.redirectErrorStream(true);
								Process proc = pb.start();
								BufferedReader idxReader = new BufferedReader(
									new InputStreamReader(proc.getInputStream()));
								String idxLine;
								while((idxLine = idxReader.readLine()) != null){
									String[] parts = idxLine.split("\t");
									if(parts.length >= 3){
										try {
											rawIndexTotal += Long.parseLong(parts[2]);
										} catch(NumberFormatException nfe){
											// skip malformed lines
										}
									}
								}
								proc.waitFor();
								if (rawIndexTotal > 0) {
									log.info("samtools idxstats: " + rawIndexTotal + " raw aligned reads");
								}
							} catch(Exception e){
								log.info("samtools idxstats not available: " + e.getMessage());
							}
						}

						// Step B: If raw total < 5M, go straight to full scan (fast enough).
						// Otherwise, sample ~50K reads from the midpoint of each major
						// chromosome via the BAM index to measure the filter pass rate.
						final long FULL_SCAN_THRESHOLD = 5_000_000L;
						final int TARGET_PER_CHROM = 50_000;
						final int REGION_SIZE = 1_000_000; // 1Mb window around each midpoint

						if(rawIndexTotal > 0 && rawIndexTotal < FULL_SCAN_THRESHOLD){
							log.info("Raw total " + rawIndexTotal + " < 5M; using full scan for exact count");
							// fall through to full scan below
						} else if(rawIndexTotal >= FULL_SCAN_THRESHOLD){
							long sampled = 0;
							long passed = 0;
							// Track unique read names among PASSING reads so we can
							// compute a fragments-per-passing-read ratio. For typical
							// PE data this is ~0.5 (each fragment has two reads sharing
							// the same readName); for SE data it's 1.0; mixed PE+SE
							// gives somewhere in between. The HashSet caps at ~ TARGET_PER_CHROM
							// × ~22 chromosomes ≈ 1M entries, ~80 MB at typical readName lengths.
							java.util.HashSet<String> sampleFragmentNames = new java.util.HashSet<>();
							SamReader sampleReader = SamReaderFactory.makeDefault()
								.validationStringency(ValidationStringency.SILENT)
								.open(new File(wgsBamFile));
							List<htsjdk.samtools.SAMSequenceRecord> sequences =
								sampleReader.getFileHeader().getSequenceDictionary().getSequences();
							int chromsWithReads = 0;
							for(htsjdk.samtools.SAMSequenceRecord seq : sequences){
								int seqLen = seq.getSequenceLength();
								if(seqLen < REGION_SIZE) continue;
								String seqName = seq.getSequenceName();
								int midpoint = seqLen / 2;
								int queryStart = Math.max(1, midpoint - REGION_SIZE / 2);
								int queryEnd = Math.min(seqLen, midpoint + REGION_SIZE / 2);
								int chromSampled = 0;
								try {
									SAMRecordIterator regionIt = sampleReader.queryOverlapping(
										seqName, queryStart, queryEnd);
									while(regionIt.hasNext() && chromSampled < TARGET_PER_CHROM){
										SAMRecord r = regionIt.next();
										sampled++;
										chromSampled++;
										if(failFlagFilter(r)){
											continue;
										}
										if(stringentPaired && !CcInferenceUtils.passReadPairOrientation(r)){
											continue;
										}
										passed++;
										sampleFragmentNames.add(r.getReadName());
									}
									regionIt.close();
								} catch(Exception regionEx){
									// Some references may not be queryable; skip
								}
								if(chromSampled > 0) chromsWithReads++;
							}
							sampleReader.close();

							if(sampled > 0 && passed > 0){
								double passRate = (double) passed / sampled;
								long filteredReadsTotal = (long)(rawIndexTotal * passRate);
								if(coverageReadLevel){
									// Legacy: denominator = filtered raw reads.
									readsNumTotal = filteredReadsTotal;
									usedIndex = true;
									log.info("Sampled " + sampled + " reads across " + chromsWithReads +
											 " chromosomes: " + passed + " passed filters (" +
											 String.format("%.1f%%", passRate * 100) + "); estimated " +
											 filteredReadsTotal + " filtered reads " +
											 "(read-level coverage scale; legacy)");
								} else {
									// New default: denominator = filtered fragments.
									// fragmentRatio = unique fragment names / passing reads.
									// For PE: ~0.5; for SE: 1.0; mixed: in-between.
									double fragmentRatio = (double) sampleFragmentNames.size() / (double) passed;
									long filteredFragmentsTotal = (long)(filteredReadsTotal * fragmentRatio);
									readsNumTotal = filteredFragmentsTotal;
									usedIndex = true;
									log.info("Sampled " + sampled + " reads across " + chromsWithReads +
											 " chromosomes: " + passed + " passed filters (" +
											 String.format("%.1f%%", passRate * 100) + "); " +
											 sampleFragmentNames.size() + " unique fragments " +
											 String.format("(fragment ratio %.3f)", fragmentRatio) +
											 "; estimated " + filteredReadsTotal + " filtered reads -> " +
											 filteredFragmentsTotal + " filtered fragments " +
											 "(fragment-level coverage scale)");
								}
							} else if(sampled > 0){
								log.warn("All " + sampled + " sampled reads were filtered; using raw index total");
								readsNumTotal = rawIndexTotal;
								usedIndex = true;
							}
						}

						// Strategy 2: Full scan (last resort, when no index available)
						if(!usedIndex){
							log.info("No BAM index available. Full scan for read count ...");
							SAMRecordIterator wgsIt = wgsReader.iterator();
							java.util.HashSet<String> fragmentNames =
									coverageReadLevel ? null : new java.util.HashSet<>();
							long readsPassed = 0;
							while(wgsIt.hasNext()){
								SAMRecord r = wgsIt.next();
								if(failFlagFilter(r)){
									continue;
								}else{
									if(stringentPaired && !CcInferenceUtils.passReadPairOrientation(r)){
										continue;
									}
								}
								readsPassed++;
								if(fragmentNames != null){
									fragmentNames.add(r.getReadName());
								}
							}
							wgsIt.close();
							if(coverageReadLevel){
								readsNumTotal = readsPassed;
								log.info("Full scan: " + readsPassed + " filtered reads " +
										"(read-level coverage scale; legacy)");
							} else {
								readsNumTotal = (fragmentNames == null ? readsPassed : fragmentNames.size());
								log.info("Full scan: " + readsPassed + " filtered reads -> " +
										readsNumTotal + " unique fragments " +
										"(fragment-level coverage scale)");
							}
						}
					}

					if(coverageReadLevel){
						log.info((long)readsNumTotal + " reads in total (read-level coverage scale) ...");
					} else {
						log.info((long)readsNumTotal + " fragments in total (fragment-level coverage scale) ...");
					}
					readsNumTotal = readsNumTotal/1000000;
					log.info("Output value for each CpG in each DNA fragment ... ");

					// Build list of 5Mb genomic bins from CpG collections
					final int BIN_SIZE = 5_000_000;
					List<String[]> genomicBins = new ArrayList<String[]>(); // [chr, binStart, binEnd]
					for(String chr : cpgCollections.keySet()){
						if(chr.equalsIgnoreCase("chrM")){
							continue;
						}
						IntervalTree<String> cpgChrCollections = cpgCollections.get(chr);
						// Find min/max coordinates for this chromosome's CpGs
						int minCoord = Integer.MAX_VALUE;
						int maxCoord = Integer.MIN_VALUE;
						Iterator<Node<String>> it = cpgChrCollections.iterator();
						while(it.hasNext()){
							Node<String> node = it.next();
							if(node.getStart() < minCoord) minCoord = node.getStart();
							if(node.getEnd() > maxCoord) maxCoord = node.getEnd();
						}
						// Create bins
						for(int binStart = (minCoord / BIN_SIZE) * BIN_SIZE; binStart <= maxCoord; binStart += BIN_SIZE){
							int binEnd = binStart + BIN_SIZE;
							// Check if any CpGs exist in this bin
							Iterator<Node<String>> binIt = cpgChrCollections.overlappers(binStart, binEnd);
							if(binIt.hasNext()){
								genomicBins.add(new String[]{chr, String.valueOf(binStart), String.valueOf(binEnd)});
							}
						}
					}

						log.info("Processing " + genomicBins.size() + " genomic bins (" + BIN_SIZE/1_000_000 + "Mb each) in parallel ...");
						long totalCpgTargets = 0;
						for(String[] bin : genomicBins){
							String binChr = bin[0];
							int binStart = Integer.parseInt(bin[1]);
							int binEnd = Integer.parseInt(bin[2]);
							IntervalTree<String> cpgChrCollections = cpgCollections.get(binChr);
							Iterator<Node<String>> binCpgIt = cpgChrCollections.overlappers(binStart, binEnd);
							while(binCpgIt.hasNext()){
								binCpgIt.next();
								totalCpgTargets++;
							}
						}
						log.info("Total CpGs scheduled for processing: " + totalCpgTargets);

					// Write header to output file using BGZF (bgzip) format.
					// Bgzipped files are valid gzip but split into ~64KB independent
					// blocks, enabling parallel decompression (via `bgzip -d -@ N`)
					// and random access via tabix index. File extension stays .gz for
					// compatibility with standard gzip readers.
					FileOutputStream output = new FileOutputStream(detailFile);
					BlockCompressedOutputStream bgzipOut = new BlockCompressedOutputStream(output, (java.io.File) null);
					OutputStreamWriter writer = new OutputStreamWriter(bgzipOut, "UTF-8");
					// Header prefixed with '#' so tabix skips it automatically.
					String distColumnHeader = useSignedDistCenter ? "Signed_Dist_Center" : "Dist_frag_end";
					writer.write("#chr\tstart\tend\treadName\tFragLen\tFrag_strand\tmethy_stat\tNorm_Frag_cov\tbaseQ\tOffset_frag\t" + distColumnHeader);
					if (useEndMotif) {
						writer.write("\tmotif_score");
					}
					if(includeCpgDist){
						writer.write("\tdist_nearest_CpG");
					}
					writer.write(header + "\n");
					writer.flush();

					// Capture final variables for lambda/anonymous class access
					final double finalReadsNumTotal = readsNumTotal;
					final String finalRefFile = refFile;
					final String finalWgsBamFile = wgsBamFile;
					final boolean finalUseTabixFragmentInput = useTabixFragmentInput;
					final HashMap<String,IntervalTree<String>> finalCpgCollections = cpgCollections;
					final HashMap<String,IntervalTree<String>> finalAllCpgLocCollections = allCpgLocCollections;
					final HashMap<String, HashMap<String,IntervalTree<Integer>>> finalOverlapLocStringCollections = overlapLocStringCollections;
					final HashMap<String, HashMap<String,IntervalTree<String>>> finalDistantLocStringCollections = distantLocStringCollections;
					final LinkedHashSet<String> finalOverlapLocString = overlapLocString;
						final LinkedHashSet<String> finalDistantLocString = distantLocString;
						final LinkedHashSet<String> finalValueBedLocString = valueBedLocString;
						final LinkedHashSet<String> finalValueWigLocString = valueWigLocString;
						final LinkedHashSet<String> finalKmerCollections = kmerCollections;
						final long finalTotalCpgTargets = Math.max(1L, totalCpgTargets);

					// Parallel processing of genomic bins
					int nThreads = threads > 0 ? threads : Runtime.getRuntime().availableProcessors();
					nThreads = Math.max(1, nThreads);
					log.info("Using " + nThreads + " threads for parallel genomic-bin processing ...");
					ExecutorService executor = Executors.newFixedThreadPool(nThreads);
					AtomicLong globalCpgCount = new AtomicLong(0);
					AtomicLong globalPoints = new AtomicLong(0);

					// Create temp files for each bin, submit tasks
					List<File> tempFiles = new ArrayList<File>();
					List<Future<Long>> futures = new ArrayList<Future<Long>>();

					for(int binIdx = 0; binIdx < genomicBins.size(); binIdx++){
						final String[] bin = genomicBins.get(binIdx);
						final String binChr = bin[0];
						final int binStart = Integer.parseInt(bin[1]);
						final int binEnd = Integer.parseInt(bin[2]);
						final File tempFile = File.createTempFile("cpg_bin_" + binIdx + "_", ".tmp");
						tempFile.deleteOnExit();
						tempFiles.add(tempFile);

						futures.add(executor.submit(new Callable<Long>() {
							@Override
							public Long call() throws Exception {
								long binPoints = 0;

								// Per-thread I/O resources
								SamReader binReader = null;
								TabixReader binFragmentReader = null;
								if(finalUseTabixFragmentInput){
									binFragmentReader = new TabixReader(finalWgsBamFile);
								}else{
									binReader = SamReaderFactory.makeDefault()
										.validationStringency(ValidationStringency.SILENT)
										.open(new File(finalWgsBamFile));
								}
								TwoBitParser binRefParser = new TwoBitParser(new File(finalRefFile));
								binRefParser.setCurrentSequence(binChr);

								// Per-thread TabixFeatureReaders
								HashMap<String, Pair<Integer, TabixFeatureReader<BEDFeature, ?>>> binValueBedReaders = null;
								if(valueBeds != null){
									binValueBedReaders = new HashMap<String, Pair<Integer, TabixFeatureReader<BEDFeature, ?>>>();
									for(String valueBedString : valueBeds){
										String[] splitStrings = valueBedString.split(":");
										String vbName = splitStrings[0];
										int vbExt = Integer.parseInt(splitStrings[1]);
										String vbFile = splitStrings[2];
										binValueBedReaders.put(vbName, new Pair<Integer, TabixFeatureReader<BEDFeature, ?>>(vbExt, new TabixFeatureReader(vbFile, new BEDCodec())));
									}
								}

								// Per-thread BigWigFileReaders
								HashMap<String, Pair<Integer, BigWigFileReader>> binValueWigReaders = null;
								if(valueWigs != null){
									binValueWigReaders = new HashMap<String, Pair<Integer, BigWigFileReader>>();
									for(String valueWigString : valueWigs){
										String[] splitStrings = valueWigString.split(":");
										String vwName = splitStrings[0];
										int vwExt = Integer.parseInt(splitStrings[1]);
										String vwFile = splitStrings[2];
										binValueWigReaders.put(vwName, new Pair<Integer, BigWigFileReader>(vwExt, new BigWigFileReader(new File(vwFile).toPath())));
									}
								}

								OutputStreamWriter binWriter = new OutputStreamWriter(new FileOutputStream(tempFile), "UTF-8");

								String bamChr = binChr;
								if(useNoChrPrefixBam){
									Pattern replace = Pattern.compile("^chr");
									Matcher matcher1 = replace.matcher(bamChr);
									bamChr = matcher1.replaceAll("");
								}

									IntervalTree<String> cpgChrCollections = finalCpgCollections.get(binChr);
									Iterator<Node<String>> cpgIterator = cpgChrCollections.overlappers(binStart, binEnd);
									List<Node<String>> cpgNodes = new ArrayList<Node<String>>();
									while(cpgIterator.hasNext()){
										cpgNodes.add(cpgIterator.next());
									}
									Collections.sort(cpgNodes, new Comparator<Node<String>>() {
										@Override
										public int compare(Node<String> a, Node<String> b) {
											int cmp = Integer.compare(a.getStart(), b.getStart());
											if(cmp != 0){
												return cmp;
											}
											return Integer.compare(a.getEnd(), b.getEnd());
										}
									});

									if(!finalUseTabixFragmentInput){
										// Streaming read window: single pass over BAM iterator for this bin
										SAMRecordIterator binReadIt = binReader.queryOverlapping(bamChr, binStart + 1, binEnd);
										List<SAMRecord> activeReads = new ArrayList<SAMRecord>();
										SAMRecord nextCandidateRead = null;
										long localCpgCount = 0;
										// Per-fragment motif cache: avoids redundant .2bit lookups when
										// the same read covers multiple CpGs. Cleared per CpG.
										HashMap<String, Double> motifCache = useEndMotif ? new HashMap<>() : null;
										HashMap<String, String> motifStringCache = (useEndMotif && saveMotifLookup != null) ? new HashMap<>() : null;

										for(Node<String> cpg : cpgNodes){
											int start = cpg.getStart();
											int end = cpg.getEnd();
											int fragMostLeft = start+1;
											int fragMostRight = end;

											// Add reads whose alignment start is now in range for current CpG
											while(true){
												if(nextCandidateRead == null){
													nextCandidateRead = nextValidBinRead(binReadIt);
												}
												if(nextCandidateRead == null){
													break;
												}
												if(nextCandidateRead.getAlignmentStart() <= end){
													activeReads.add(nextCandidateRead);
													nextCandidateRead = null;
												}else{
													break;
												}
											}

											// Expire reads that cannot overlap this CpG or any downstream CpG
											for(Iterator<SAMRecord> activeIt = activeReads.iterator(); activeIt.hasNext();){
												SAMRecord activeRead = activeIt.next();
												if(activeRead.getAlignmentEnd() < (start + 1)){
													activeIt.remove();
												}
											}

											HashMap<String, SAMRecord> countedReads = new HashMap<String, SAMRecord>();
											int readNumber = 0;
											for(SAMRecord r : activeReads){
												// We stream by alignment windows; keep exact overlap semantics per CpG.
												if(r.getAlignmentStart() > end || r.getAlignmentEnd() < (start + 1)){
													continue;
												}
												readNumber++;
												boolean negStrand = r.getReadNegativeStrandFlag();
												boolean secondEnd = r.getReadPairedFlag() && r.getSecondOfPairFlag();
												if(secondEnd){
													negStrand = !negStrand;
												}

												int bisulfitePos = 0;
												if(!wgsMode){
													if(r.getTransientAttribute("BS") != null){
														bisulfitePos = Integer.parseInt((String) r.getTransientAttribute("BS"));
													}else{
														bisulfitePos = CcInferenceUtils.bisulfiteIncompleteReads(r);
														r.setTransientAttribute("BS", bisulfitePos);
													}
												}

												int offSet = r.getReadPositionAtReferencePosition(end)-1;
												if(bisulfitePos < 0){
													continue;
												}else if(bisulfitePos > 0){
													if((!negStrand && offSet < bisulfitePos) || (negStrand && offSet >= bisulfitePos)){
														continue;
													}
												}
												if(offSet<0){
													continue;
												}

												if(r.getAlignmentStart() < fragMostLeft){
													fragMostLeft = r.getAlignmentStart();
												}
												if(r.getMateAlignmentStart() < fragMostLeft){
													fragMostLeft = r.getMateAlignmentStart();
												}

												if(r.getAlignmentEnd() > fragMostRight){
													fragMostRight = r.getAlignmentEnd();
												}
												int mateEnd = CcInferenceUtils.getMateAlignmentEndByMateCigar(r);
												if(mateEnd > fragMostRight){
													fragMostRight = mateEnd;
												}

												String readName = r.getReadName();
												if(countedReads.containsKey(readName)){
													SAMRecord prev = countedReads.get(readName);
													int offSetPrev = prev.getReadPositionAtReferencePosition(end)-1;
													if(offSet < r.getBaseQualities().length && offSetPrev < prev.getBaseQualities().length){
														byte baseQ = r.getBaseQualities()[offSet];
														byte base = CcInferenceUtils.toUpperCase(r.getReadBases()[offSet]);
														byte baseQPrev = prev.getBaseQualities()[offSetPrev];
														byte basePrev = CcInferenceUtils.toUpperCase(prev.getReadBases()[offSetPrev]);

														if(!BaseUtilsMore.basesAreEqual(base, basePrev)){
															if(baseQ > baseQPrev){
																countedReads.put(readName, r);
															}else if(baseQ == baseQPrev && !secondEnd){
																countedReads.put(readName, r);
															}
														}
													}
												}else{
													countedReads.put(readName, r);
												}
											}

											if(readNumber >= maxCov || countedReads.size()==0){
												localCpgCount++;
												if(localCpgCount % 1000 == 0){
													long total = globalCpgCount.addAndGet(1000);
													logCpgProgress(total, finalTotalCpgTargets);
												}
												continue;
											}

										// SCALE: by default v0.64+, BAM coverage is fragment-level
										// to match tabix mode: numerator = countedReads.size()
										// (deduplicated by readName) and denominator = total
										// filtered fragments. Use -coverageReadLevel for legacy
										// pre-v0.64 behavior (raw reads on both numerator and
										// denominator). See tutorial.md §4.3.
										double normalizedFragCov = coverageReadLevel
												? (double) readNumber / finalReadsNumTotal
												: (double) countedReads.size() / finalReadsNumTotal;

										byte[] refBasesExt = CcInferenceUtils.toUpperCase(binRefParser.loadFragment(end-1-kmerExt, kmerExt*2+1).getBytes());
										byte refBase = refBasesExt[kmerExt];

										HashMap<String, Double> kmerMapsRef = new HashMap<String, Double>();
										if(!useFragBaseKmer){
											for(int j = 2; j <= kmerLen; j++){
												kmerMapsRef.putAll(CcInferenceUtils.kmerFreqSearch(refBasesExt, j));
											}
										}

										// nearest CpG distance
										double nearestCpg = Double.NaN;
										if(includeCpgDist){
											IntervalTree<String> cpgLocCollections = finalAllCpgLocCollections.get(binChr);
											Iterator<Node<String>> upstreamCpgIt = null;
											Iterator<Node<String>> downstreamCpgIt = null;
											if(BaseUtilsMore.basesAreEqual(refBase, BaseUtilsMore.C)){
												upstreamCpgIt = cpgLocCollections.reverseIterator(start-1, end-1);
												downstreamCpgIt = cpgLocCollections.iterator(start+2, end+2);
											}else if(BaseUtilsMore.basesAreEqual(refBase, BaseUtilsMore.G)){
												upstreamCpgIt = cpgLocCollections.reverseIterator(start-2, end-2);
												downstreamCpgIt = cpgLocCollections.iterator(start+1, end+1);
												}else{
													localCpgCount++;
													if(localCpgCount % 1000 == 0){
														long total = globalCpgCount.addAndGet(1000);
														logCpgProgress(total, finalTotalCpgTargets);
													}
													continue;
												}

											if(upstreamCpgIt== null || !upstreamCpgIt.hasNext()){
												IntervalTree.Node<String> downstream = downstreamCpgIt.next();
												nearestCpg = CcInferenceUtils.intervalDistance(downstream,cpg);
											}else if(downstreamCpgIt==null || !downstreamCpgIt.hasNext()){
												IntervalTree.Node<String> upstream = upstreamCpgIt.next();
												nearestCpg = CcInferenceUtils.intervalDistance(upstream,cpg);
											}else{
												IntervalTree.Node<String> upstream = upstreamCpgIt.next();
												IntervalTree.Node<String> downstream = downstreamCpgIt.next();
												int dist1 = CcInferenceUtils.intervalDistance(upstream, cpg);
												int dist2 = CcInferenceUtils.intervalDistance(downstream, cpg);
												nearestCpg = Math.abs(dist1) < Math.abs(dist2) ? dist1 : dist2;
											}
										}

										// overlap with feature in reference genome
										HashMap<String, Integer> overlapStatCollections = new HashMap<String, Integer>();
										if(overlapRegions!=null && !overlapRegions.isEmpty()){
											for(String key : finalOverlapLocStringCollections.keySet()){
												HashMap<String,IntervalTree<Integer>> tmp = finalOverlapLocStringCollections.get(key);
												if(tmp.containsKey(binChr)){
													overlapStatCollections.put(key, tmp.get(binChr).minOverlapper(start, end)==null ? 0 : 1);
												}else{
													overlapStatCollections.put(key, 0);
												}
											}
										}

										// distance with feature in reference genome
										HashMap<String, Integer> distStatCollections = new HashMap<String, Integer>();
										if(distantRegions!=null && !distantRegions.isEmpty()){
											for(String key : finalDistantLocStringCollections.keySet()){
												IntervalTree<String> locCollections = finalDistantLocStringCollections.get(key).get(binChr);
												int distanceNearest = Integer.MAX_VALUE;
												if(locCollections!=null && locCollections.size()>0){
													Iterator<Node<String>> upstreamIt = locCollections.reverseIterator(start, end);
													Iterator<Node<String>> downstreamIt = locCollections.iterator(start, end);
													if(!upstreamIt.hasNext()){
														IntervalTree.Node<String> downstream = locCollections.min(start, end);
														distanceNearest = CcInferenceUtils.intervalDistance(downstream,cpg);
													}else if(!downstreamIt.hasNext()){
														IntervalTree.Node<String> upstream = locCollections.max(start, end);
														distanceNearest = CcInferenceUtils.intervalDistance(upstream,cpg);
													}else{
														IntervalTree.Node<String> upstream = locCollections.max(start, end);
														IntervalTree.Node<String> downstream = locCollections.min(start, end);
														int dist1 = CcInferenceUtils.intervalDistance(upstream, cpg);
														int dist2 = CcInferenceUtils.intervalDistance(downstream, cpg);
														distanceNearest = Math.abs(dist1) < Math.abs(dist2) ? dist1 : dist2;
													}
												}
												distStatCollections.put(key, distanceNearest);
											}
										}

										// value in bed file
										HashMap<String, Double> valBedStatCollections = new HashMap<String, Double>();
										if(binValueBedReaders != null){
											for(String key : binValueBedReaders.keySet()){
												int range = binValueBedReaders.get(key).getFirst();
												boolean mean0 = range < 0;
												if(mean0) range = -range;
												TabixFeatureReader<BEDFeature, ?> bedReader = binValueBedReaders.get(key).getSecond();
												CloseableTribbleIterator<BEDFeature> featureIt = bedReader.query(binChr, (start-range < 0 ? 1 : start-range+1), end+range);
												DescriptiveStatistics statFeature = new DescriptiveStatistics();
												while(featureIt.hasNext()){
													BEDFeature term = featureIt.next();
													if(!Double.isNaN(term.getScore())) statFeature.addValue(term.getScore());
												}
												featureIt.close();
												if(statFeature.getN()>0){
													valBedStatCollections.put(key, mean0 ? (double)statFeature.getSum()/(double)(range*2+1) : statFeature.getMean());
												}else{
													valBedStatCollections.put(key, Double.NaN);
												}
											}
										}

										// value in wig file
										HashMap<String, Double> valWigStatCollections = new HashMap<String, Double>();
										if(binValueWigReaders != null){
											for(String key : binValueWigReaders.keySet()){
												int range = binValueWigReaders.get(key).getFirst();
												BigWigFileReader wigReader = binValueWigReaders.get(key).getSecond();
												if(range < 0){
													range = -range;
													SummaryStatistics statFeature = wigReader.queryStats(binChr, (start-range < 0 ? 1 : start-range), end+range);
													valWigStatCollections.put(key, statFeature.getN()>0 ? (double)statFeature.getSum()/(double)(range*2+1) : Double.NaN);
												}else{
													SummaryStatistics statFeature = wigReader.queryStats(binChr, start-range, end+range);
													valWigStatCollections.put(key, statFeature.getN()>0 ? statFeature.getMean() : Double.NaN);
												}
											}
										}

										for(String readName : countedReads.keySet()){
											SAMRecord r = countedReads.get(readName);
											boolean negStrand = r.getReadNegativeStrandFlag();
											boolean secondEnd = r.getReadPairedFlag() && r.getSecondOfPairFlag();
											if(secondEnd){
												negStrand = !negStrand;
											}
											int offSet = r.getReadPositionAtReferencePosition(end)-1;
											if(offSet<0) continue;
											char methyStat = '.';
											byte[] bases = r.getReadBases();
											byte base = bases[offSet];
											byte[] baseQs = r.getBaseQualities();
											byte baseQ = baseQs[offSet];
											if(baseQ <= minBaseQ) continue;

											if(negStrand){
												if(BaseUtilsMore.basesAreEqual(refBase, BaseUtilsMore.G)){
													if(BaseUtilsMore.basesAreEqual(base, BaseUtilsMore.G)){
														methyStat = 'm';
													}else if(BaseUtilsMore.basesAreEqual(base, BaseUtilsMore.A)){
														methyStat = 'u';
													}else{
														continue;
													}
												}else{
													continue;
												}
											}else{
												if(BaseUtilsMore.basesAreEqual(refBase, BaseUtilsMore.C)){
													if(BaseUtilsMore.basesAreEqual(base, BaseUtilsMore.C)){
														methyStat = 'm';
													}else if(BaseUtilsMore.basesAreEqual(base, BaseUtilsMore.T)){
														methyStat = 'u';
													}else{
														continue;
													}
												}else{
													continue;
												}
											}

											int fragLen = Math.abs(r.getInferredInsertSize());
											if(fragLen > maxFragLen) continue;
											int cpgOffset = CcInferenceUtils.getFragOffsetFromReadsOffset(r, offSet);
											if(cpgOffset < 0 || cpgOffset >= fragLen) continue;
											int distToFragEnd = CcInferenceUtils.getDistFragEndFromReadsOffset(r, offSet);
											if(distToFragEnd < 0) continue;
											if(distToFragEnd > maxDistToFragEnd) continue;

											char fragStrand = negStrand ? '-' : '+';

											int fragStart = Math.min(r.getAlignmentStart(), r.getAlignmentStart()+r.getInferredInsertSize());
											int fragEnd = Math.max(r.getAlignmentStart(), r.getAlignmentStart()+r.getInferredInsertSize());
											if(r.getInferredInsertSize() == 0) continue;

											// Kmer reference sequence: only load when kmer features are requested.
											// Previously this large loadFragment() was always called even when
											// unused, wasting ~500-1000 bytes per read per CpG.
											HashMap<String, Double> kmerMapsFrag = new HashMap<String, Double>();
											if (useFragBaseKmer || useStrandSpecificFragBase) {
												try {
													byte[] refBasesFrag = binRefParser.loadFragment(fragMostLeft, fragMostRight-fragMostLeft+1).getBytes();
													if(negStrand && useStrandSpecificFragBase){
														SequenceUtil.reverseComplement(refBasesFrag);
													}
													if(useFragBaseKmer){
														for(int j = 2; j <= kmerLen; j++){
															kmerMapsFrag.putAll(CcInferenceUtils.kmerFreqSearch(refBasesFrag, j));
														}
													}
												} catch (Exception e) {
													// TwoBitParser can fail at N-block boundaries; skip kmer for this CpG
												}
											}

											// 5' end motif extraction and scoring (BAM path).
											// The motif is a per-fragment property (same for all CpGs within
											// the same fragment). Cache by read name within this CpG bin to
											// avoid redundant .2bit lookups when the same read covers multiple CpGs.
											String bamMotif = null;
											double bamMotifScore = Double.NaN;
											if (useEndMotif) {
												String motifCacheKey = readName;
												Double cachedScore = motifCache.get(motifCacheKey);
												if (cachedScore != null) {
													bamMotifScore = cachedScore;
													// Lookup motif string only if training mode needs it
													if (saveMotifLookup != null) {
														bamMotif = motifStringCache.get(motifCacheKey);
													}
												} else {
													bamMotif = extractFivePrimeMotif(binRefParser, fragStart - 1, fragEnd, negStrand);
													bamMotifScore = getMotifScore(bamMotif);
													motifCache.put(motifCacheKey, bamMotifScore);
													if (saveMotifLookup != null) {
														accumulateMotifCount(bamMotif, methyStat);
														motifStringCache.put(motifCacheKey, bamMotif);
													}
												}
											}

											double covValue = noCoverage ? Double.NaN : normalizedFragCov;
											// Distance feature: legacy unsigned distance-to-nearest-end,
											// or signed strand-aware distance-from-center when -useSignedDistCenter.
											// The unsigned distToFragEnd computed above is still used for the
											// -maxDistToFragEnd filter (so the filter semantics don't change).
											String distColValue = useSignedDistCenter
													? String.format("%.1f", signedDistFromCenter(cpgOffset, fragLen, negStrand))
													: Integer.toString(distToFragEnd);
											binWriter.write(binChr + "\t" + start + "\t" + end + "\t" + readName + "\t" + fragLen + "\t" + fragStrand + "\t" + methyStat + "\t" + String.format("%.6f",covValue)
													 + "\t" + (int)baseQ + "\t" + cpgOffset + "\t" + distColValue);
											if (useEndMotif) {
												if (saveMotifLookup != null && loadMotifLookup == null) {
													// Training mode: write 4-mer string as placeholder.
													// The concatenation phase replaces it with the real score
													// after all bins have contributed to motifCounts.
													binWriter.write("\t" + (bamMotif != null ? bamMotif : "NNNN"));
												} else {
													binWriter.write("\t" + (Double.isNaN(bamMotifScore) ? "NaN" : String.format("%.6f", bamMotifScore)));
												}
											}
											if(includeCpgDist) binWriter.write("\t" + nearestCpg);
											if(overlapStatCollections.size()>0){
												for(String key : finalOverlapLocString) binWriter.write("\t" + overlapStatCollections.get(key));
											}
											if(distStatCollections.size() > 0){
												for(String key : finalDistantLocString) binWriter.write("\t" + distStatCollections.get(key));
											}
											if(valBedStatCollections.size()>0){
												for(String key : finalValueBedLocString) binWriter.write("\t" + String.format("%.3f",valBedStatCollections.get(key)));
											}
											if(valWigStatCollections.size()>0){
												for(String key : finalValueWigLocString) binWriter.write("\t" + String.format("%.3f",valWigStatCollections.get(key)));
											}
											if(kmerMapsRef.size()>0){
												for(String key : finalKmerCollections) binWriter.write("\t" + String.format("%.3f",kmerMapsRef.get(key)));
											}
											if(kmerMapsFrag.size()>0){
												for(String key : finalKmerCollections) binWriter.write("\t" + String.format("%.3f",kmerMapsFrag.get(key)));
											}
											binWriter.write("\n");
											binPoints++;
										}
											localCpgCount++;
											if(localCpgCount % 1000 == 0){
												long total = globalCpgCount.addAndGet(1000);
												logCpgProgress(total, finalTotalCpgTargets);
													binWriter.flush();
												}
											}
										binReadIt.close();

											// Flush remaining CpG count
											if(localCpgCount % 1000 != 0){
												long total = globalCpgCount.addAndGet(localCpgCount % 1000);
												logCpgProgress(total, finalTotalCpgTargets);
											}
									}else{
										TabixReader.Iterator fragIt = queryTabixIterator(binFragmentReader, bamChr, binStart, binEnd);
										ArrayList<FragmentRecord> binFragments = new ArrayList<FragmentRecord>();
										long fragLineNumber = 0;
										if(fragIt != null){
											String fragLine;
											while((fragLine = fragIt.next()) != null){
												fragLineNumber++;
												FragmentRecord fragment = parseFragmentRecord(fragLine, fragLineNumber);
												if(fragment == null){
													continue;
												}
												if(fragment.end <= binStart || fragment.start >= binEnd){
													continue;
												}
												if(fragment.end <= fragment.start){
													continue;
												}
												binFragments.add(fragment);
											}
										}

										Collections.sort(binFragments, new Comparator<FragmentRecord>() {
											@Override
											public int compare(FragmentRecord a, FragmentRecord b) {
												int cmp = Integer.compare(a.start, b.start);
												if(cmp != 0){
													return cmp;
												}
												return Integer.compare(a.end, b.end);
											}
										});

										int fragmentCursor = 0;
										ArrayList<FragmentRecord> activeFragments = new ArrayList<FragmentRecord>();
										long localCpgCount = 0;
										// Per-fragment motif cache (Tabix path)
										HashMap<String, Double> tabixMotifCache = useEndMotif ? new HashMap<>() : null;
										HashMap<String, String> tabixMotifStringCache = (useEndMotif && saveMotifLookup != null) ? new HashMap<>() : null;

										for(Node<String> cpg : cpgNodes){
											int start = cpg.getStart();
											int end = cpg.getEnd();

											while(fragmentCursor < binFragments.size() && binFragments.get(fragmentCursor).start < end){
												activeFragments.add(binFragments.get(fragmentCursor));
												fragmentCursor++;
											}
											for(Iterator<FragmentRecord> activeIt = activeFragments.iterator(); activeIt.hasNext();){
												FragmentRecord activeFragment = activeIt.next();
												if(activeFragment.end <= start){
													activeIt.remove();
												}
											}

											HashMap<String, FragmentRecord> countedReads = new HashMap<String, FragmentRecord>();
											int readNumber = 0;
											for(FragmentRecord fragment : activeFragments){
												if(fragment.start >= end || fragment.end <= start){
													continue;
												}
												readNumber++;
												if(!countedReads.containsKey(fragment.readName)){
													countedReads.put(fragment.readName, fragment);
												}
											}

											if(readNumber >= maxCov || countedReads.size() == 0){
												localCpgCount++;
												if(localCpgCount % 1000 == 0){
													long total = globalCpgCount.addAndGet(1000);
													logCpgProgress(total, finalTotalCpgTargets);
												}
												continue;
											}

											double normalizedFragCov = (double)readNumber/finalReadsNumTotal;

											byte[] refBasesExt = CcInferenceUtils.toUpperCase(binRefParser.loadFragment(end-1-kmerExt, kmerExt*2+1).getBytes());
											byte refBase = refBasesExt[kmerExt];

											HashMap<String, Double> kmerMapsRef = new HashMap<String, Double>();
											if(!useFragBaseKmer){
												for(int j = 2; j <= kmerLen; j++){
													kmerMapsRef.putAll(CcInferenceUtils.kmerFreqSearch(refBasesExt, j));
												}
											}

											// nearest CpG distance
											double nearestCpg = Double.NaN;
											if(includeCpgDist){
												IntervalTree<String> cpgLocCollections = finalAllCpgLocCollections.get(binChr);
												Iterator<Node<String>> upstreamCpgIt = null;
												Iterator<Node<String>> downstreamCpgIt = null;
												if(BaseUtilsMore.basesAreEqual(refBase, BaseUtilsMore.C)){
													upstreamCpgIt = cpgLocCollections.reverseIterator(start-1, end-1);
													downstreamCpgIt = cpgLocCollections.iterator(start+2, end+2);
												}else if(BaseUtilsMore.basesAreEqual(refBase, BaseUtilsMore.G)){
													upstreamCpgIt = cpgLocCollections.reverseIterator(start-2, end-2);
													downstreamCpgIt = cpgLocCollections.iterator(start+1, end+1);
												}else{
													localCpgCount++;
													if(localCpgCount % 1000 == 0){
														long total = globalCpgCount.addAndGet(1000);
														logCpgProgress(total, finalTotalCpgTargets);
													}
													continue;
												}

												if(upstreamCpgIt== null || !upstreamCpgIt.hasNext()){
													IntervalTree.Node<String> downstream = downstreamCpgIt.next();
													nearestCpg = CcInferenceUtils.intervalDistance(downstream,cpg);
												}else if(downstreamCpgIt==null || !downstreamCpgIt.hasNext()){
													IntervalTree.Node<String> upstream = upstreamCpgIt.next();
													nearestCpg = CcInferenceUtils.intervalDistance(upstream,cpg);
												}else{
													IntervalTree.Node<String> upstream = upstreamCpgIt.next();
													IntervalTree.Node<String> downstream = downstreamCpgIt.next();
													int dist1 = CcInferenceUtils.intervalDistance(upstream, cpg);
													int dist2 = CcInferenceUtils.intervalDistance(downstream, cpg);
													nearestCpg = Math.abs(dist1) < Math.abs(dist2) ? dist1 : dist2;
												}
											}

											// overlap with feature in reference genome
											HashMap<String, Integer> overlapStatCollections = new HashMap<String, Integer>();
											if(overlapRegions!=null && !overlapRegions.isEmpty()){
												for(String key : finalOverlapLocStringCollections.keySet()){
													HashMap<String,IntervalTree<Integer>> tmp = finalOverlapLocStringCollections.get(key);
													if(tmp.containsKey(binChr)){
														overlapStatCollections.put(key, tmp.get(binChr).minOverlapper(start, end)==null ? 0 : 1);
													}else{
														overlapStatCollections.put(key, 0);
													}
												}
											}

											// distance with feature in reference genome
											HashMap<String, Integer> distStatCollections = new HashMap<String, Integer>();
											if(distantRegions!=null && !distantRegions.isEmpty()){
												for(String key : finalDistantLocStringCollections.keySet()){
													IntervalTree<String> locCollections = finalDistantLocStringCollections.get(key).get(binChr);
													int distanceNearest = Integer.MAX_VALUE;
													if(locCollections!=null && locCollections.size()>0){
														Iterator<Node<String>> upstreamIt = locCollections.reverseIterator(start, end);
														Iterator<Node<String>> downstreamIt = locCollections.iterator(start, end);
														if(!upstreamIt.hasNext()){
															IntervalTree.Node<String> downstream = locCollections.min(start, end);
															distanceNearest = CcInferenceUtils.intervalDistance(downstream,cpg);
														}else if(!downstreamIt.hasNext()){
															IntervalTree.Node<String> upstream = locCollections.max(start, end);
															distanceNearest = CcInferenceUtils.intervalDistance(upstream,cpg);
														}else{
															IntervalTree.Node<String> upstream = locCollections.max(start, end);
															IntervalTree.Node<String> downstream = locCollections.min(start, end);
															int dist1 = CcInferenceUtils.intervalDistance(upstream, cpg);
															int dist2 = CcInferenceUtils.intervalDistance(downstream, cpg);
															distanceNearest = Math.abs(dist1) < Math.abs(dist2) ? dist1 : dist2;
														}
													}
													distStatCollections.put(key, distanceNearest);
												}
											}

											// value in bed file
											HashMap<String, Double> valBedStatCollections = new HashMap<String, Double>();
											if(binValueBedReaders != null){
												for(String key : binValueBedReaders.keySet()){
													int range = binValueBedReaders.get(key).getFirst();
													boolean mean0 = range < 0;
													if(mean0) range = -range;
													TabixFeatureReader<BEDFeature, ?> bedReader = binValueBedReaders.get(key).getSecond();
													CloseableTribbleIterator<BEDFeature> featureIt = bedReader.query(binChr, (start-range < 0 ? 1 : start-range+1), end+range);
													DescriptiveStatistics statFeature = new DescriptiveStatistics();
													while(featureIt.hasNext()){
														BEDFeature term = featureIt.next();
														if(!Double.isNaN(term.getScore())) statFeature.addValue(term.getScore());
													}
													featureIt.close();
													if(statFeature.getN()>0){
														valBedStatCollections.put(key, mean0 ? (double)statFeature.getSum()/(double)(range*2+1) : statFeature.getMean());
													}else{
														valBedStatCollections.put(key, Double.NaN);
													}
												}
											}

											// value in wig file
											HashMap<String, Double> valWigStatCollections = new HashMap<String, Double>();
											if(binValueWigReaders != null){
												for(String key : binValueWigReaders.keySet()){
													int range = binValueWigReaders.get(key).getFirst();
													BigWigFileReader wigReader = binValueWigReaders.get(key).getSecond();
													if(range < 0){
														range = -range;
														SummaryStatistics statFeature = wigReader.queryStats(binChr, (start-range < 0 ? 1 : start-range), end+range);
														valWigStatCollections.put(key, statFeature.getN()>0 ? (double)statFeature.getSum()/(double)(range*2+1) : Double.NaN);
													}else{
														SummaryStatistics statFeature = wigReader.queryStats(binChr, start-range, end+range);
														valWigStatCollections.put(key, statFeature.getN()>0 ? statFeature.getMean() : Double.NaN);
													}
												}
											}

											int cpgPos = end - 1;
											for(String readName : countedReads.keySet()){
												FragmentRecord fragment = countedReads.get(readName);
												boolean negStrand = fragment.strand == '-';
												if(negStrand){
													if(!BaseUtilsMore.basesAreEqual(refBase, BaseUtilsMore.G)){
														continue;
													}
												}else{
													if(!BaseUtilsMore.basesAreEqual(refBase, BaseUtilsMore.C)){
														continue;
													}
												}

												char methyStat = resolveFragmentMethyStat(fragment);
												if(methyStat != 'm' && methyStat != 'u'){
													continue;
												}

												int fragLen = fragment.end - fragment.start;
												if(fragLen <= 0 || fragLen > maxFragLen){
													continue;
												}
												int cpgOffset = cpgPos - fragment.start;
												if(cpgOffset < 0 || cpgOffset >= fragLen){
													continue;
												}
												int distToFragEnd = Math.min(cpgOffset, (fragLen - 1) - cpgOffset);
												if(distToFragEnd > maxDistToFragEnd){
													continue;
												}

												int baseQ = fragment.baseQ;
												if(baseQ <= minBaseQ){
													continue;
												}
												char fragStrand = fragment.strand;

												HashMap<String, Double> kmerMapsFrag = new HashMap<String, Double>();
												if(useFragBaseKmer){
													byte[] refBasesFrag = binRefParser.loadFragment(fragment.start, fragLen).getBytes();
													if(negStrand && useStrandSpecificFragBase){
														SequenceUtil.reverseComplement(refBasesFrag);
													}
													for(int j = 2; j <= kmerLen; j++){
														kmerMapsFrag.putAll(CcInferenceUtils.kmerFreqSearch(refBasesFrag, j));
													}
												}

												// 5' end motif extraction and scoring (Tabix path)
												// Cache by readName: same fragment at multiple CpGs gets same motif.
												String tabixMotif = null;
												double tabixMotifScore = Double.NaN;
												if (useEndMotif) {
													Double cachedScore = tabixMotifCache.get(readName);
													if (cachedScore != null) {
														tabixMotifScore = cachedScore;
														if (saveMotifLookup != null) {
															tabixMotif = tabixMotifStringCache.get(readName);
														}
													} else {
														tabixMotif = extractFivePrimeMotif(binRefParser, fragment.start, fragment.end, negStrand);
														tabixMotifScore = getMotifScore(tabixMotif);
														tabixMotifCache.put(readName, tabixMotifScore);
														if (saveMotifLookup != null) {
															accumulateMotifCount(tabixMotif, methyStat);
															tabixMotifStringCache.put(readName, tabixMotif);
														}
													}
												}

												double tabixCovValue = noCoverage ? Double.NaN : normalizedFragCov;
												// Distance feature: legacy unsigned distance-to-nearest-end,
												// or signed strand-aware distance-from-center when -useSignedDistCenter.
												// distToFragEnd above is still used for the -maxDistToFragEnd filter
												// so the filter semantics are unchanged across the two modes.
												String distColValue = useSignedDistCenter
														? String.format("%.1f", signedDistFromCenter(cpgOffset, fragLen, negStrand))
														: Integer.toString(distToFragEnd);
												binWriter.write(binChr + "\t" + start + "\t" + end + "\t" + readName + "\t" + fragLen + "\t" + fragStrand + "\t" + methyStat + "\t" + String.format("%.6f",tabixCovValue)
														 + "\t" + baseQ + "\t" + cpgOffset + "\t" + distColValue);
												if (useEndMotif) {
													if (saveMotifLookup != null && loadMotifLookup == null) {
														// Training mode: write 4-mer string as placeholder
														binWriter.write("\t" + (tabixMotif != null ? tabixMotif : "NNNN"));
													} else {
														binWriter.write("\t" + (Double.isNaN(tabixMotifScore) ? "NaN" : String.format("%.6f", tabixMotifScore)));
													}
												}
												if(includeCpgDist) binWriter.write("\t" + nearestCpg);
												if(overlapStatCollections.size()>0){
													for(String key : finalOverlapLocString) binWriter.write("\t" + overlapStatCollections.get(key));
												}
												if(distStatCollections.size() > 0){
													for(String key : finalDistantLocString) binWriter.write("\t" + distStatCollections.get(key));
												}
												if(valBedStatCollections.size()>0){
													for(String key : finalValueBedLocString) binWriter.write("\t" + String.format("%.3f",valBedStatCollections.get(key)));
												}
												if(valWigStatCollections.size()>0){
													for(String key : finalValueWigLocString) binWriter.write("\t" + String.format("%.3f",valWigStatCollections.get(key)));
												}
												if(kmerMapsRef.size()>0){
													for(String key : finalKmerCollections) binWriter.write("\t" + String.format("%.3f",kmerMapsRef.get(key)));
												}
												if(kmerMapsFrag.size()>0){
													for(String key : finalKmerCollections) binWriter.write("\t" + String.format("%.3f",kmerMapsFrag.get(key)));
												}
												binWriter.write("\n");
												binPoints++;
											}

											localCpgCount++;
											if(localCpgCount % 1000 == 0){
												long total = globalCpgCount.addAndGet(1000);
												logCpgProgress(total, finalTotalCpgTargets);
												binWriter.flush();
											}
										}

										if(localCpgCount % 1000 != 0){
											long total = globalCpgCount.addAndGet(localCpgCount % 1000);
											logCpgProgress(total, finalTotalCpgTargets);
										}
									}

								binWriter.close();
								if(binReader != null){
									binReader.close();
								}
								if(binFragmentReader != null){
									binFragmentReader.close();
								}
								binRefParser.closeParser();
								if(binValueBedReaders != null){
									for(String key : binValueBedReaders.keySet()) binValueBedReaders.get(key).getSecond().close();
								}
								if(binValueWigReaders != null){
									for(String key : binValueWigReaders.keySet()) binValueWigReaders.get(key).getSecond().close();
								}

								globalPoints.addAndGet(binPoints);
								return binPoints;
							}
						}));
					}

					// Wait for all bin tasks to complete
					try {
						for(Future<Long> future : futures){
							future.get();
						}
						} catch (Exception e) {
							throw new RuntimeException("Error in parallel feature extraction", e);
						}
						executor.shutdown();
						logCpgProgress(globalCpgCount.get(), finalTotalCpgTargets);

					// If training mode (-saveMotifLookup without -loadMotifLookup),
					// compute motifScoreMap from accumulated counts BEFORE concatenation
					// so that the output contains real motif scores, not placeholders.
					if (saveMotifLookup != null && motifCounts != null && motifScoreMap == null) {
						log.info("Computing motif scores from " + motifCounts.size() + " accumulated 4-mers ...");
						motifScoreMap = new HashMap<>();
						for (java.util.Map.Entry<String, java.util.concurrent.atomic.AtomicLongArray> entry : motifCounts.entrySet()) {
							long methylated = entry.getValue().get(0);
							long total = entry.getValue().get(1);
							motifScoreMap.put(entry.getKey(), (double)(methylated + 1) / (total + 2));
						}
						log.info("Motif scores computed. Replacing placeholder values in output ...");
					}

					// Concatenate temp files into the output.
					// When motifScoreMap was just computed (training mode), replace the
					// NaN motif_score column (col 11) with actual values looked up from
					// each line's fragment coordinates in the .2bit reference.
					boolean needMotifRewrite = (useEndMotif && saveMotifLookup != null && loadMotifLookup == null);
					// motif_score is column index 11 (0-based) when -useEndMotif is set
					final int MOTIF_COL = 11;

					for(File tempFile : tempFiles){
						if (!needMotifRewrite) {
							// Fast path: raw byte copy when no rewriting needed
							byte[] buffer = new byte[8192];
							FileInputStream fis = new FileInputStream(tempFile);
							int bytesRead;
							while((bytesRead = fis.read(buffer)) != -1){
								writer.write(new String(buffer, 0, bytesRead, "UTF-8"));
							}
							fis.close();
						} else {
							// Line-by-line rewrite: replace motif_score NaN with real value.
							// Each line has readName at col[3], FragLen at col[4], strand at col[5],
							// fragStart can be derived from col[1](cpg_start) and col[9](offset).
							// But it's simpler to look up the motif from the readName using the
							// per-read motif cache that was populated during processing... except
							// caches are per-bin and gone. Instead, extract the motif score from
							// the methy_stat (col[6]) and the motif lookup computed above.
							// Since we don't have the 4-mer stored in the output, we need to
							// recompute it. But that requires .2bit access again.
							//
							// Simpler approach: store the 4-mer in a hidden column during
							// processing, or accept that we must re-derive it.
							//
							// SIMPLEST: during parallel processing, when saveMotifLookup mode,
							// write the 4-mer string in the motif_score column (not NaN or 0.5).
							// Then here, replace the 4-mer with its score.
							BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(tempFile), "UTF-8"));
							String line;
							while ((line = br.readLine()) != null) {
								String[] cols = line.split("\t", MOTIF_COL + 2);
								if (cols.length > MOTIF_COL) {
									String motifKey = cols[MOTIF_COL].trim();
									double score = motifScoreMap != null ? motifScoreMap.getOrDefault(motifKey, 0.5) : 0.5;
									cols[MOTIF_COL] = String.format("%.6f", score);
									// Rejoin — but we split with limit, so cols[MOTIF_COL+1] has the rest
									StringBuilder sb = new StringBuilder();
									for (int ci = 0; ci < cols.length; ci++) {
										if (ci > 0) sb.append('\t');
										sb.append(cols[ci]);
									}
									writer.write(sb.toString());
								} else {
									writer.write(line);
								}
								writer.write("\n");
							}
							br.close();
						}
						tempFile.delete();
					}

					writer.close();
					output.close();

					// Build tabix index (.tbi) for the bgzipped output.
					// The output columns start with chr/start/end, so we use
					// the standard BED format with 0-based half-open semantics.
					// The index enables parallel block-level reads downstream.
					try {
						log.info("Building tabix index for " + detailFile + " ...");
						buildTabixIndex(detailFile);
						log.info("Tabix index written to " + detailFile + ".tbi");
					} catch (Exception tbxEx) {
						log.warn("Failed to build tabix index (output is still bgzipped and readable): " + tbxEx.getMessage());
					}

					// Save motif lookup table
					if (saveMotifLookup != null && motifCounts != null) {
						saveMotifLookupFile(saveMotifLookup);
					}

					if(wgsReader != null){
						wgsReader.close();
					}
					refParser.closeParser();

					if(valueBedReaders != null){
						for(String key : valueBedReaders.keySet()){
							valueBedReaders.get(key).getSecond().close();
						}
					}
					if(valueWigReaders != null){
						for(String key : valueWigReaders.keySet()){
							valueWigReaders.get(key).getSecond().close();
						}
					}

					points = globalPoints.get();
					finish();

		}
		
		private boolean isTabixFragmentInput(String fragmentInputFile){
			String lower = fragmentInputFile.toLowerCase(Locale.US);
			if(fragmentInputTabix){
				return true;
			}
			return !(lower.endsWith(".bam") || lower.endsWith(".cram"));
		}

		private void validateTabixIndexExists(String tabixFile){
			File tbi = new File(tabixFile + ".tbi");
			if(!tbi.exists()){
				throw new IllegalArgumentException("Tabix index not found for fragment input: " + tabixFile + " (expected " + tabixFile + ".tbi)");
			}
		}

		/**
		 * Estimate the total number of fragments in a bgzipped tabix-indexed
		 * fragment file using a fast sampling-based approach analogous to the
		 * BAM-index estimator. Order of strategies, fastest first:
		 *
		 *   1. BGZF sampling (default for .bed.gz / .tsv.gz with bgzip header):
		 *        Read first {@code TABIX_COUNT_SAMPLE_LINES} lines via
		 *        BlockCompressedInputStream, track the compressed-byte
		 *        offset (from the BGZF virtual file pointer's high 48 bits),
		 *        compute compressed_bytes_per_line, then extrapolate over
		 *        the total compressed file size. Typical cost: tens of
		 *        milliseconds on a multi-GB file.
		 *
		 *   2. Full scan (fallback): the legacy line-by-line counter, used
		 *        when the input is not bgzip (so virtual file pointers are
		 *        unavailable), when sampling fails for any reason, or when
		 *        the file is small enough that scanning is cheaper than
		 *        opening BCI.
		 *
		 * Returns an exact count for path 2 and an approximate count for
		 * path 1. If users need an exact count, they can pass it directly
		 * via -totalReadsInBam.
		 */
		private long estimateTotalFragmentsFromTabixInput(String fragmentInputFile) throws IOException{
			Long sampled = estimateTotalFragmentsFromTabixSampling(fragmentInputFile);
			if(sampled != null){
				return sampled;
			}
			log.info("BGZF sampling not available for " + fragmentInputFile + "; falling back to full scan");
			return estimateTotalFragmentsFromTabixFullScan(fragmentInputFile);
		}

		/**
		 * Cap on the number of records read during the BGZF sample. Large
		 * enough to span many chromosomes naturally (a sorted tabix file
		 * crosses chromosome boundaries every few hundred to few thousand
		 * records on whole-genome data), small enough to stay under ~1 second.
		 */
		private static final int TABIX_COUNT_SAMPLE_LINES = 200_000;

		/**
		 * Minimum chromosomes the sample must span before we accept the
		 * estimate. If we haven't crossed at least this many chromosome
		 * boundaries, the bytes-per-record ratio may be biased toward whatever
		 * chromosome happens to be at the start of the file. Below this we
		 * keep reading until we either reach the line cap or run out of file.
		 */
		private static final int TABIX_COUNT_MIN_CHROMOSOMES = 5;

		/**
		 * Minimum compressed bytes consumed before we trust the
		 * bytes-per-record ratio. If the sample fits inside the first BGZF
		 * block (no block boundary crossed), the consumed-byte counter
		 * rounds to 0 and the estimate would be infinite.
		 */
		private static final long TABIX_COUNT_MIN_COMPRESSED_BYTES = 1L << 16; // 64 KB

		/**
		 * BGZF stratified-by-chromosome sampler.
		 *
		 * Reads from the start of the bgzipped fragment file and tracks both
		 * the BGZF virtual file pointer (for compressed-bytes-per-record) and
		 * the set of distinct chromosomes seen so far. Stops when EITHER:
		 *
		 *   * we have read {@link #TABIX_COUNT_SAMPLE_LINES} records, OR
		 *   * we have seen records from at least
		 *     {@link #TABIX_COUNT_MIN_CHROMOSOMES} distinct chromosomes AND
		 *     have crossed at least 64 KB of compressed bytes.
		 *
		 * Tabix files are sorted by chromosome, position, so iterating forward
		 * naturally pulls samples from successive chromosomes — this is the
		 * closest practical analog to BAM-style per-chromosome sampling
		 * (BAM index supplies exact per-chromosome record counts for free,
		 * but tabix indexes only store byte offsets and BGZF block-boundary
		 * seeks are fragile to off-block compressed offsets).
		 *
		 * Falls back to a full scan when sampling fails or when the whole
		 * file fits in fewer than the line cap.
		 */
		private Long estimateTotalFragmentsFromTabixSampling(String fragmentInputFile) {
			File file = new File(fragmentInputFile);
			long fileSize = file.length();
			if (fileSize <= 0) {
				return null;
			}
			BlockCompressedInputStream bgzf = null;
			BufferedReader br = null;
			try {
				bgzf = new BlockCompressedInputStream(file);
				br = new BufferedReader(new InputStreamReader(bgzf));
				long sampledLines = 0;
				long compressedBytesConsumed = 0;
				java.util.LinkedHashSet<String> chromosomesSeen = new java.util.LinkedHashSet<>();
				String line;
				boolean exhausted = false;
				while ((line = br.readLine()) != null) {
					if (line.startsWith("#") || line.trim().isEmpty()) continue;
					int firstTab = line.indexOf('\t');
					if (firstTab <= 0) continue;
					int secondTab = line.indexOf('\t', firstTab + 1);
					if (secondTab < 0) continue; // <3 columns
					sampledLines++;
					if (chromosomesSeen.size() < TABIX_COUNT_MIN_CHROMOSOMES + 1) {
						chromosomesSeen.add(line.substring(0, firstTab));
					}
					compressedBytesConsumed = bgzf.getFilePointer() >> 16;

					boolean haveEnoughChroms =
							chromosomesSeen.size() >= TABIX_COUNT_MIN_CHROMOSOMES;
					boolean haveEnoughBytes =
							compressedBytesConsumed >= TABIX_COUNT_MIN_COMPRESSED_BYTES;
					if (haveEnoughChroms && haveEnoughBytes &&
							sampledLines >= TABIX_COUNT_SAMPLE_LINES / 4) {
						// Stop early: we've already crossed multiple chromosomes
						// and have enough compressed-byte advance for a stable
						// bytes-per-record estimate. Reading further would only
						// shrink the relative uncertainty marginally.
						break;
					}
					if (sampledLines >= TABIX_COUNT_SAMPLE_LINES) {
						break;
					}
				}
				if (br.readLine() == null && sampledLines > 0 &&
						sampledLines < TABIX_COUNT_SAMPLE_LINES) {
					exhausted = true;
				}

				if (sampledLines == 0) {
					return 0L;
				}
				if (exhausted) {
					// Whole file fits in the sample → exact count.
					log.info(String.format(
							"Tabix fragment count (whole file fits within sample): " +
									"%,d fragments across %d chromosome(s)",
							sampledLines, chromosomesSeen.size()));
					return sampledLines;
				}
				if (compressedBytesConsumed < TABIX_COUNT_MIN_COMPRESSED_BYTES) {
					// Sampling didn't cross enough compressed bytes for a
					// reliable bytes-per-record ratio. Caller falls back.
					return null;
				}
				double bytesPerLine =
						(double) compressedBytesConsumed / (double) sampledLines;
				double estimated = (double) fileSize / bytesPerLine;
				long estimatedLong = Math.max(sampledLines, Math.round(estimated));
				log.info(String.format(
						"Tabix fragment count (BGZF sampling, %d chromosomes spanned): " +
								"%,d sampled lines used %,d compressed bytes -> %.2f " +
								"bytes/line; total compressed size %,d bytes -> ~%,d fragments",
						chromosomesSeen.size(), sampledLines, compressedBytesConsumed,
						bytesPerLine, fileSize, estimatedLong));
				return estimatedLong;
			} catch (IOException e) {
				log.warn("BGZF sampling failed for " + fragmentInputFile + ": " + e.getMessage());
				return null;
			} finally {
				try { if (br != null) br.close(); } catch (IOException ignored) {}
				try { if (bgzf != null) bgzf.close(); } catch (IOException ignored) {}
			}
		}

		private long estimateTotalFragmentsFromTabixFullScan(String fragmentInputFile) throws IOException{
			long total = 0;
			BufferedReader br = openMaybeGzipReader(fragmentInputFile);
			String line;
			while((line = br.readLine()) != null){
				if(line.startsWith("#") || line.trim().isEmpty()){
					continue;
				}
				String[] splitLines = line.split("\t");
				if(splitLines.length < 3){
					continue;
				}
				total++;
			}
			br.close();
			return total;
		}

		private TabixReader.Iterator queryTabixIterator(TabixReader tabixReader, String chr, int start, int end) throws IOException{
			TabixReader.Iterator iterator = tabixReader.query(chr, Math.max(0, start), end);
			if(iterator == null){
				if(chr.startsWith("chr")){
					iterator = tabixReader.query(chr.substring(3), Math.max(0, start), end);
				}else{
					iterator = tabixReader.query("chr" + chr, Math.max(0, start), end);
				}
			}
			return iterator;
		}

		private FragmentRecord parseFragmentRecord(String line, long lineNumber){
			if(line == null || line.startsWith("#") || line.trim().isEmpty()){
				return null;
			}
			String[] splitLines = line.split("\t");
			if(splitLines.length < 3){
				return null;
			}

			String chr = splitLines[0];
			int start;
			int end;
			try{
				start = Integer.parseInt(splitLines[1]);
				end = Integer.parseInt(splitLines[2]);
			}catch(NumberFormatException e){
				return null;
			}

			char strand = '+';
			if(fragStrandColumn > 0){
				int idx = fragStrandColumn - 1;
				if(idx < splitLines.length){
					Character parsed = parseStrandToken(splitLines[idx]);
					if(parsed != null){
						strand = parsed;
					}
				}
			}else{
				if(splitLines.length >= 6){
					Character parsed = parseStrandToken(splitLines[5]);
					if(parsed != null){
						strand = parsed;
					}
				}
				if(splitLines.length >= 4){
					Character parsed = parseStrandToken(splitLines[3]);
					if(parsed != null){
						strand = parsed;
					}
				}
			}

			// Resolve readName. The HMM uses readName to GROUP CpG observations
			// belonging to the same fragment (the sequential observations along
			// one fragment); collapsing all rows to a single readName collapses
			// the training set to one giant pseudo-fragment of millions of
			// observations and produces a singular GMM covariance matrix.
			//
			// In tabix mode each row IS one fragment, so we synthesize a
			// coordinate+strand+line-number readName by default - guaranteed
			// unique per row. The user can override via -fragNameColumn when
			// the input has genuine unique read names (e.g., a BAM-derived
			// fragment file with original SAM read names). Even then we reject
			// BED placeholder values (".", "*", empty) that otherwise collapse
			// every row to one "read".
			//
			// We deliberately do NOT auto-detect a name from col4/col5 (older
			// behavior): that consumed BED placeholder values and produced
			// "1 unique read" pathologies that broke training silently.
			String readName = null;
			if(fragNameColumn > 0){
				int idx = fragNameColumn - 1;
				if(idx < splitLines.length){
					String candidate = splitLines[idx].trim();
					if(!candidate.isEmpty() && !candidate.equals(".") && !candidate.equals("*")){
						readName = candidate;
					}
				}
			}
			if(readName == null || readName.isEmpty()){
				readName = "frag_" + lineNumber + "_" + chr + "_" + start + "_" + end + "_" + strand;
			}

			char methyStat = '.';
			if(fragMethyColumn > 0){
				int idx = fragMethyColumn - 1;
				if(idx < splitLines.length){
					Character parsedMethy = parseMethyToken(splitLines[idx]);
					if(parsedMethy != null){
						methyStat = parsedMethy;
					}
				}
			}else{
				for(int i = 3; i < splitLines.length; i++){
					Character parsedMethy = parseMethyToken(splitLines[i]);
					if(parsedMethy != null){
						methyStat = parsedMethy;
						break;
					}
				}
			}

			return new FragmentRecord(chr, start, end, strand, readName, methyStat, fragBaseQ);
		}

		/**
		 * Resolve the methy_stat output column for a tabix-fragment row.
		 *
		 * Order of precedence:
		 *   1. If the input row had an explicit 'm'/'u' token (from
		 *      -fragMethyColumn or auto-detected in extra columns), use it.
		 *   2. Otherwise, return -defaultMethyStat (default 'm').
		 *
		 * methy_stat is NOT consumed by HMM training or decoding -- the HMM
		 * is unsupervised over the feature vector and predicts the methylation
		 * label as the hidden state. methy_stat is written to the output for
		 * AUC/QC reporting only. Under -wgsMode for BAM input, every covered
		 * CpG is also labeled 'm' for the same reason (the read base is the
		 * unconverted reference C/G), so 'm' is the natural default for tabix
		 * fragment input as well.
		 */
		private char resolveFragmentMethyStat(FragmentRecord fragment){
			if(fragment.methyStat == 'm' || fragment.methyStat == 'u'){
				return fragment.methyStat;
			}
			return defaultMethyStat.toLowerCase(Locale.US).startsWith("m") ? 'm' : 'u';
		}

		private boolean isStrandToken(String token){
			return parseStrandToken(token) != null;
		}

		private Character parseStrandToken(String token){
			if(token == null){
				return null;
			}
			String value = token.trim();
			if(value.equals("+") || value.equalsIgnoreCase("f") || value.equals("1")){
				return '+';
			}
			if(value.equals("-") || value.equalsIgnoreCase("r") || value.equals("-1")){
				return '-';
			}
			return null;
		}

		private Character parseMethyToken(String token){
			if(token == null){
				return null;
			}
			String value = token.trim().toLowerCase(Locale.US);
			if(value.equals("m") || value.equals("meth") || value.equals("methylated")){
				return 'm';
			}
			if(value.equals("u") || value.equals("unmeth") || value.equals("unmethylated")){
				return 'u';
			}
			return null;
		}
		
		private static class FragmentRecord {
			final String chr;
			final int start;
			final int end;
			final char strand;
			final String readName;
			final char methyStat;
			final int baseQ;

			FragmentRecord(String chr, int start, int end, char strand, String readName, char methyStat, int baseQ){
				this.chr = chr;
				this.start = start;
				this.end = end;
				this.strand = strand;
				this.readName = readName;
				this.methyStat = methyStat;
				this.baseQ = baseQ;
			}
		}
		
		
		private SAMRecord nextValidBinRead(SAMRecordIterator iterator){
			while(iterator.hasNext()){
				SAMRecord r = iterator.next();
				if(failFlagFilter(r)){
					continue;
				}
				if(stringentPaired && !CcInferenceUtils.passReadPairOrientation(r)){
					continue;
				}
				return r;
			}
			return null;
		}

		private void logCpgProgress(long processedCpgs, long totalCpgs){
			long cappedProcessed = Math.min(processedCpgs, totalCpgs);
			double percent = totalCpgs <= 0 ? 100.0 : (100.0 * cappedProcessed / (double) totalCpgs);
			log.info(String.format(Locale.US, "CpG progress: %,d/%,d (%.2f%%)", cappedProcessed, totalCpgs, percent));
		}
		
		
		private boolean failFlagFilter(SAMRecord r){
			return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
					|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag() || !r.getProperPairFlag()
					|| (skipSecondEnd && r.getReadPairedFlag() && r.getSecondOfPairFlag());
		}
	

	
	
	private void initiate(){
		startTime = System.currentTimeMillis();
	}

	private void finish(){
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;

		log.info("Counted " + points + " data points in total");
		log.info("CpgFeatureMatrixBuilder's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}

	/**
	 * Extract the 4-mer at the 5' end of a fragment from the .2bit reference genome.
	 *
	 * For + strand fragments: 5' end is at fragStart, extract [fragStart, fragStart+4).
	 * For - strand fragments: 5' end is at fragEnd, extract reverse complement of [fragEnd-4, fragEnd).
	 *
	 * Coordinate convention (0-based half-open, matching BED/TwoBitParser):
	 *   fragStart = 0-based start position
	 *   fragEnd   = 0-based exclusive end position
	 *
	 * For BAM input: fragStart = SAMRecord.getAlignmentStart() - 1,
	 *                fragEnd   = getAlignmentStart() + getInferredInsertSize()
	 *                (already 0-based exclusive since TLEN includes the +1).
	 * For Tabix input: coords are already 0-based half-open.
	 *
	 * @param binRefParser TwoBitParser set to the current chromosome
	 * @param fragStart    0-based fragment start (genomic)
	 * @param fragEnd      0-based exclusive fragment end (genomic)
	 * @param negStrand    whether fragment is on negative strand
	 * @return uppercase 4-mer string, or null if lookup fails or contains N
	 */
	private String extractFivePrimeMotif(TwoBitParser binRefParser, int fragStart, int fragEnd, boolean negStrand) {
		try {
			if (fragStart < 0 || fragEnd <= fragStart) return null;
			String seq;
			if (!negStrand) {
				// + strand: 5' is at fragStart, extract [fragStart, fragStart+4)
				if (fragStart + 4 > fragEnd) return null;  // fragment too short
				seq = binRefParser.loadFragment(fragStart, 4);
			} else {
				// - strand: 5' is at fragEnd-1, extract RC of [fragEnd-4, fragEnd)
				int motifStart = fragEnd - 4;
				if (motifStart < 0 || motifStart < fragStart) return null;
				seq = binRefParser.loadFragment(motifStart, 4);
				byte[] seqBytes = seq.getBytes();
				SequenceUtil.reverseComplement(seqBytes);
				seq = new String(seqBytes);
			}
			seq = seq.toUpperCase();
			if (seq.length() != 4 || seq.contains("N")) return null;
			return seq;
		} catch (Exception e) {
			return null;
		}
	}

	/**
	 * Signed CpG-to-fragment-center distance, strand-corrected so that
	 * negative = CpG is in the 5' half of its fragment, positive = 3' half.
	 *
	 * Both BAM and tabix code paths compute cpgOffset as the offset from
	 * the LEFTMOST coordinate of the fragment (i.e., from fragment.start).
	 * On + strand fragments, the leftmost coordinate IS the 5' end, so
	 * cpgOffset already measures from the 5'. On - strand fragments, the
	 * leftmost coordinate is the 3' end, so we flip to (fragLen-1)-cpgOffset
	 * to get the offset from the 5'. Then subtract fragLen/2 to center.
	 */
	private static double signedDistFromCenter(int cpgOffset, int fragLen, boolean negStrand) {
		int offsetFrom5Prime = negStrand ? (fragLen - 1 - cpgOffset) : cpgOffset;
		return offsetFrom5Prime - fragLen / 2.0;
	}

	/**
	 * Look up motif score from the loaded motif score map.
	 * Returns 0.5 (neutral) for null/absent motifs.
	 */
	private double getMotifScore(String motif) {
		if (motif == null || motifScoreMap == null) return 0.5;
		return motifScoreMap.getOrDefault(motif, 0.5);
	}

	/**
	 * Accumulate motif counts for training mode (thread-safe).
	 */
	private void accumulateMotifCount(String motif, char methyStat) {
		if (motif == null || motifCounts == null) return;
		java.util.concurrent.atomic.AtomicLongArray counts = motifCounts.computeIfAbsent(motif,
				k -> new java.util.concurrent.atomic.AtomicLongArray(2));
		counts.incrementAndGet(1); // total
		if (methyStat == 'm') {
			counts.incrementAndGet(0); // methylated
		}
	}

	/**
	 * Load motif score lookup table from TSV file.
	 * Format: header line, then rows of "motif\tscore".
	 */
	private HashMap<String, Double> loadMotifLookupFile(String path) throws IOException {
		HashMap<String, Double> map = new HashMap<>();
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line;
		boolean header = true;
		while ((line = br.readLine()) != null) {
			if (header) { header = false; continue; }
			String[] parts = line.split("\t");
			if (parts.length >= 2) {
				map.put(parts[0].trim(), Double.parseDouble(parts[1].trim()));
			}
		}
		br.close();
		return map;
	}

	/**
	 * Build a tabix index (.tbi) for a bgzipped TSV file whose first three
	 * columns are chr, start, end. Uses htsjdk's TabixIndexCreator with
	 * BED-like semantics (0-based start, half-open end).
	 */
	private void buildTabixIndex(String bgzippedFile) throws IOException {
		// Feature input: read bgzipped file, use a BED-style codec to emit
		// a minimal Feature per line (chr, start, end). The tabix index
		// only needs chr/start/end, so we build a lightweight one.
		htsjdk.tribble.index.tabix.TabixFormat fmt = new htsjdk.tribble.index.tabix.TabixFormat(
			htsjdk.tribble.index.tabix.TabixFormat.GENERIC_FLAGS,
			1,      // seq column (1-based: chr)
			2,      // begin column (1-based: start)
			3,      // end column (1-based: end)
			'#',    // meta (comment) char
			0       // number of header lines to skip (header starts with # or "chr"/"start")
		);
		htsjdk.tribble.index.tabix.TabixIndexCreator creator =
			new htsjdk.tribble.index.tabix.TabixIndexCreator(fmt);

		htsjdk.samtools.util.BlockCompressedInputStream bcis =
			new htsjdk.samtools.util.BlockCompressedInputStream(new java.io.File(bgzippedFile));
		try {
			long virtualOffset = bcis.getFilePointer();
			String line;
			while ((line = bcis.readLine()) != null) {
				if (line.length() == 0 || line.charAt(0) == '#') {
					virtualOffset = bcis.getFilePointer();
					continue;
				}
				String[] cols = line.split("\t", 4);
				if (cols.length < 3) {
					virtualOffset = bcis.getFilePointer();
					continue;
				}
				try {
					String chr = cols[0];
					final int start = Integer.parseInt(cols[1]);
					final int end = Integer.parseInt(cols[2]);
					final String finalChr = chr;
					htsjdk.tribble.Feature feat = new htsjdk.tribble.Feature() {
						@Override public String getContig() { return finalChr; }
						@Override public int getStart() { return start + 1; } // tabix uses 1-based
						@Override public int getEnd() { return end; }
					};
					creator.addFeature(feat, virtualOffset);
				} catch (NumberFormatException nfe) {
					// skip malformed lines
				}
				virtualOffset = bcis.getFilePointer();
			}
		} finally {
			bcis.close();
		}
		htsjdk.tribble.index.Index idx = creator.finalizeIndex(0L);
		idx.writeBasedOnFeatureFile(new java.io.File(bgzippedFile));
	}

	/**
	 * Save motif score lookup table to TSV file.
	 * Uses Laplace smoothing: score = (methylated + 1) / (total + 2).
	 */
	private void saveMotifLookupFile(String path) throws IOException {
		OutputStreamWriter mw = new OutputStreamWriter(new FileOutputStream(path), "UTF-8");
		mw.write("motif\tscore\n");
		for (java.util.Map.Entry<String, java.util.concurrent.atomic.AtomicLongArray> entry : motifCounts.entrySet()) {
			long methylated = entry.getValue().get(0);
			long total = entry.getValue().get(1);
			double score = (double)(methylated + 1) / (total + 2);  // Laplace smoothing
			mw.write(entry.getKey() + "\t" + String.format("%.6f", score) + "\n");
		}
		mw.close();
		log.info("Motif lookup saved to " + path + " with " + motifCounts.size() + " motifs");
	}

}
