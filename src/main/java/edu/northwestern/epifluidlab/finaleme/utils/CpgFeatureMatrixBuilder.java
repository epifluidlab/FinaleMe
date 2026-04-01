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
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.IntervalTree.Node;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
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

	@Option(name="-totalReadsInBam",usage="total number of reads used to normalize coverage column. default estimate from bam file by program. Default: -1")
	public long totalReadsInBam = -1;
	
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

	@Option(name="-t",usage="number of threads for parallel 5Mb bin processing. Use >0 to set explicitly; default uses all available cores.")
	public int threads = -1;
	
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "CpgFeatureMatrixBuilder [opts] hg19.2bit cpg_list.bed all_cpg.bed wgs.bam cpg_detail.txt.gz";
	
	private static final Logger log = LoggerFactory.getLogger(CpgFeatureMatrixBuilder.class);

	private static long startTime = -1;
	private static long points = 0;


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

					//read input bed file, for each row,
					//String intervalFile = arguments.get(0);
					String refFile = arguments.get(0);
					String cpgListFile = arguments.get(1);
					String allCpgFile = arguments.get(2);
					String wgsBamFile = arguments.get(3);
					String detailFile = arguments.get(4);

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
					SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(wgsBamFile));
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
					if(totalReadsInBam > 0){
						log.info("Get total reads number used for scaling from input option -totalReadsInBam ... ");
						readsNumTotal = totalReadsInBam;
					}else{
						// Use BAM index statistics for fast approximate read count
						log.info("Get total reads number used for scaling from bam index... ");
						boolean usedIndex = false;
						if(wgsReader.indexing() != null && wgsReader.indexing().hasBrowseableIndex()){
							try {
								htsjdk.samtools.BAMIndexMetaData[] indexStats = null;
								htsjdk.samtools.BAMIndex bamIndex = wgsReader.indexing().getIndex();
								int nRefs = wgsReader.getFileHeader().getSequenceDictionary().getSequences().size();
								for(int refIdx = 0; refIdx < nRefs; refIdx++){
									htsjdk.samtools.BAMIndexMetaData meta = bamIndex.getMetaData(refIdx);
									if(meta != null){
										readsNumTotal += meta.getAlignedRecordCount();
									}
								}
								usedIndex = true;
								log.info("Estimated " + (long)readsNumTotal + " aligned reads from BAM index (fast path)");
							} catch(Exception e){
								log.info("Failed to get counts from BAM index, falling back to full scan...");
								usedIndex = false;
								readsNumTotal = 0;
							}
						}
						if(!usedIndex){
							log.info("Get total reads number used for scaling from bam file (full scan)... ");
							SAMRecordIterator wgsIt = wgsReader.iterator();
							while(wgsIt.hasNext()){
								SAMRecord r = wgsIt.next();
								if(failFlagFilter(r)){
									continue;
								}else{
									if(stringentPaired && !CcInferenceUtils.passReadPairOrientation(r)){
										continue;
									}
								}
								readsNumTotal++;
							}
							wgsIt.close();
						}
					}
					
					log.info((long)readsNumTotal + " reads in total ...");
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

					// Write header to output file
					FileOutputStream output = new FileOutputStream(detailFile);
					OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
					writer.write("chr\tstart\tend\treadName\tFragLen\tFrag_strand\tmethy_stat\tNorm_Frag_cov\tbaseQ\tOffset_frag\tDist_frag_end");
					if(includeCpgDist){
						writer.write("\tdist_nearest_CpG");
					}
					writer.write(header + "\n");
					writer.flush();

					// Capture final variables for lambda/anonymous class access
					final double finalReadsNumTotal = readsNumTotal;
					final String finalRefFile = refFile;
					final String finalWgsBamFile = wgsBamFile;
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
								SamReader binReader = SamReaderFactory.makeDefault()
									.validationStringency(ValidationStringency.SILENT)
									.open(new File(finalWgsBamFile));
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

									// Streaming read window: single pass over BAM iterator for this bin
									SAMRecordIterator binReadIt = binReader.queryOverlapping(bamChr, binStart + 1, binEnd);
									List<SAMRecord> activeReads = new ArrayList<SAMRecord>();
									SAMRecord nextCandidateRead = null;
									long localCpgCount = 0;

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
										int distToFragEnd = CcInferenceUtils.getDistFragEndFromReadsOffset(r, offSet);
										if(distToFragEnd > maxDistToFragEnd) continue;

										char fragStrand = negStrand ? '-' : '+';

										int fragStart = Math.min(r.getAlignmentStart(), r.getAlignmentStart()+r.getInferredInsertSize());
										int fragEnd = Math.max(r.getAlignmentStart(), r.getAlignmentStart()+r.getInferredInsertSize());
										if(r.getInferredInsertSize() == 0) continue;

										byte[] refBasesFrag = binRefParser.loadFragment(fragMostLeft, fragMostRight-fragMostLeft+1).getBytes();
										if(negStrand && useStrandSpecificFragBase){
											SequenceUtil.reverseComplement(refBasesFrag);
										}
										HashMap<String, Double> kmerMapsFrag = new HashMap<String, Double>();
										if(useFragBaseKmer){
											for(int j = 2; j <= kmerLen; j++){
												kmerMapsFrag.putAll(CcInferenceUtils.kmerFreqSearch(refBasesFrag, j));
											}
										}

										binWriter.write(binChr + "\t" + start + "\t" + end + "\t" + readName + "\t" + fragLen + "\t" + fragStrand + "\t" + methyStat + "\t" + String.format("%.6f",normalizedFragCov)
												 + "\t" + (int)baseQ + "\t" + cpgOffset + "\t" + distToFragEnd);
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

								binWriter.close();
								binReader.close();
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

					// Concatenate temp files into the output
					byte[] buffer = new byte[8192];
					for(File tempFile : tempFiles){
						FileInputStream fis = new FileInputStream(tempFile);
						int bytesRead;
						while((bytesRead = fis.read(buffer)) != -1){
							writer.write(new String(buffer, 0, bytesRead, "UTF-8"));
						}
						fis.close();
						tempFile.delete();
					}

					writer.close();
					output.close();

					wgsReader.close();
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
	

}
