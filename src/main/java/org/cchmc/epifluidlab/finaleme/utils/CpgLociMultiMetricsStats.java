/**
 * CpgMultiMetricsStats.java
 * Feb 27, 2016
 * 5:28:37 PM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.finaleme.utils;


import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.IntervalTree.Node;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.Logger;
import org.biojava.nbio.genome.parsers.twobit.TwoBitParser;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.unc.genomics.io.BigWigFileReader;

/**
 *
 */
public class CpgLociMultiMetricsStats {

	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;

	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-maxFragLen",usage="maximum fragment length allowed to check. Default: 1000")
	public int maxFragLen = 1000;

	@Option(name="-maxCov",usage="maximum coverage allowed to check. Default: 250")
	public int maxCov = 250;
	
	@Option(name="-kmerLen",usage="the K-mer length to check. Default: 4")
	public int kmerLen = 4;

	@Option(name="-kmerString",usage="the fiel contain selected K-mer to check. Otherwise, use -kmerLen to automately generate all the k-mer. Default: null")
	public String kmerString = null;

	@Option(name="-kmerExt",usage="the +/- region in reference genome to check the k-mer frequency. default is +/- 200bp around CpGs. Default: 200")
	public int kmerExt = 200;
	
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
		
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "CpgLociMultiMetricsStats [opts] hg19.2bit cpg_list.bed wgs.bam cpg_loc_detail.txt.gz";
	
	private static Logger log = Logger.getLogger(CpgLociMultiMetricsStats.class);

	private static long startTime = -1;
	private static long points = 0;


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		CpgLociMultiMetricsStats clmms = new CpgLociMultiMetricsStats();
		//BasicConfigurator.configure();
		clmms.doMain(args);
	}
	
	@SuppressWarnings("resource")
	public void doMain(String[] args)
			throws Exception {

					CmdLineParser parser = new CmdLineParser(this);
					//parser.setUsageWidth(80);
					try
					{
						if(help || args.length < 4) throw new CmdLineException(parser, USAGE, new Throwable());
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
					String wgsBamFile = arguments.get(2);
					String detailFile = arguments.get(3);

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
					
					HashMap<String,IntervalTree<Integer>> ignoreLocCollections = null;
					if(excludeRegions!=null && !excludeRegions.isEmpty()){
						log.info("Excluding intervals ... ");
						ignoreLocCollections = new HashMap<String,IntervalTree<Integer>>();
						
						for(String excludeRegion : excludeRegions){
							BufferedReader br = new BufferedReader(new FileReader(excludeRegion));
							String line;
							
							while( (line = br.readLine()) != null){
								if(line.startsWith("#"))
									continue;
								String[] splitLines = line.split("\t");
								if(splitLines.length<3){
									continue;
								}
								String chr = splitLines[0];
								int start = Integer.parseInt(splitLines[1])+1;
								int end = Integer.parseInt(splitLines[2]);
								IntervalTree<Integer> tree;
								if(ignoreLocCollections.containsKey(chr)){
									tree = ignoreLocCollections.get(chr);
								}else{
									tree = new IntervalTree<Integer>();
								}
								tree.put(start, end, 1);
								ignoreLocCollections.put(chr, tree);
							}
							br.close();
						
						}
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
								int start = Integer.parseInt(splitLines[1])+1;
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
									int start = Integer.parseInt(splitLines[1])+1;
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
						
							for(String key : kmerCollections){
								header = header + "\t" + key;
							}
						
						
					}
					
					
					log.info("Loading CpG interval file ... ");	 //CpG file need to be merged two strand file (so the distance to next CpG is not always the other strand one...)
					IntervalList cpgCollectionsUnsorted = new IntervalList(wgsReader.getFileHeader());
					BufferedReader br = new BufferedReader(new FileReader(cpgListFile));
					String line;
					while( (line = br.readLine()) != null){
						if(line.startsWith("#"))
							continue;
						String[] splitLines = line.split("\t");
						if(splitLines.length<3){
							continue;
						}else{
							String chr = splitLines[0];
							int start = Integer.parseInt(splitLines[1])+1;
							int end = Integer.parseInt(splitLines[2]);
							if(splitLines.length < 6){
								cpgCollectionsUnsorted.add(new Interval(chr, start, end));
							}else{
								boolean negStrand = splitLines[5].equalsIgnoreCase("-") ? true : false;
								cpgCollectionsUnsorted.add(new Interval(chr, start, end, negStrand, splitLines[3]));
							}
						}
						
					}
					br.close();
					IntervalList cpgCollections = cpgCollectionsUnsorted.sorted();
					
					log.info("Get total reads number used for scaling ... ");
					SAMRecordIterator wgsIt = wgsReader.iterator();
					double readsNumTotal = 0;
					while(wgsIt.hasNext()){
						wgsIt.next();
						readsNumTotal++;
					}
					wgsIt.close();
					log.info((long)readsNumTotal + " reads in total ...");
					readsNumTotal = readsNumTotal/1000000;
					log.info("Output value for each CpG loci ... ");
					FileOutputStream output = new FileOutputStream(detailFile);
					OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
					
					//basic
					writer.write("chr\tstart\tend\tmethyLevel\tnorm_frag_cov\tcpg_strand\tmethylationReads\tfragSB\tdist_nearest_CpG");
					//record 10 percentile FragLen and Offset in Fragment
					writer.write("\tfragLen_mean\tfragLen_sd");
					for(int i=0;i<10;i++){
						writer.write("\tfragLen_" + i);
						
					}
					writer.write("\tfragOffset_mean\tfragOffset_sd");
					for(int i=0;i<10;i++){
						writer.write("\tfragOffset_" + i);
						
					}
					writer.write("\tfragOffsetRatio_mean\tfragOffsetRatio_sd");
					for(int i=0;i<10;i++){
						writer.write("\tfragOffsetRatio_" + i);
						
					}
					writer.write("\tbaseQ_mean\tbaseQ_sd");
					for(int i=0;i<10;i++){
						writer.write("\tbaseQ_" + i);
						
					}
					writer.write(header + "\n");
					
					
					String prevChr = "";
					for( int i = 0; i < cpgCollections.size(); i++){
						//log.info("start" + i);
						Interval cpg = cpgCollections.getIntervals().get(i);
						String chr = cpg.getContig();
						int start = cpg.getStart();
						int end = cpg.getEnd();
						int fragMostLeft = start;
						int fragMostRight = end;
						wgsIt = wgsReader.queryOverlapping(chr,start,end);
						HashMap<String, SAMRecord> countedReads = new HashMap<String, SAMRecord>();
						
						//log.info("testincpg" + i + "\t" + cpgCollections.size() + "\t" + wgsIt.hasNext() + "\t" + chr + "\t" + start + "\t" + end);
						int readNumber =0;
						while(wgsIt.hasNext()){
							SAMRecord r = wgsIt.next();
							//log.info(r.getReadName() + "\t" + failFlagFilter(r));
							if(failFlagFilter(r)){
								continue;
							}else{
								readNumber++;
								boolean negStrand = r.getReadNegativeStrandFlag();
								boolean secondEnd = r.getReadPairedFlag() && r.getSecondOfPairFlag();
								if(secondEnd){
									negStrand = !negStrand;
								}
								
								int bisulfitePos = 0;
								if(!wgsMode){
									if(r.getTransientAttribute("BS") != null){ //if the reads had been processed before
										bisulfitePos = Integer.parseInt((String) r.getTransientAttribute("BS"));
									}else{
										bisulfitePos = CcInferenceUtils.bisulfiteIncompleteReads(r);
										r.setTransientAttribute("BS", bisulfitePos);
										
									}
								}
								
								//log.info("testin" + readsNumTotal + "\t" + bisulfitePos);
								int offSet = r.getReadPositionAtReferencePosition(end)-1;
								if(bisulfitePos < 0){
									continue;
								}else if(bisulfitePos > 0){
									if((!negStrand && offSet < bisulfitePos) || (negStrand && offSet >= bisulfitePos)){
										continue;
									}
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
								if(r.getAlignmentStart() + r.getInferredInsertSize() > fragMostRight){
									fragMostRight = r.getAlignmentStart() + r.getInferredInsertSize();
								}
								
								String readName = r.getReadName();
								//int fragLen = Math.abs(r.getInferredInsertSize());
								if(countedReads.containsKey(readName)){//to filter overlapped fragments, which affect a lot in cfDNA
									byte baseQ = r.getBaseQualities()[offSet];
									byte base = r.getReadBases()[offSet];
									
									SAMRecord prev = countedReads.get(readName);
									int offSetPrev = prev.getReadPositionAtReferencePosition(end)-1;
									byte baseQPrev = prev.getBaseQualities()[offSetPrev];
									byte basePrev = prev.getReadBases()[offSetPrev];
									//int fragLenPrev = Math.abs(prev.getInferredInsertSize());
									
									if(!BaseUtils.basesAreEqual(base, basePrev)){
										if(baseQ > baseQPrev){
											countedReads.put(readName, r);
											
										}else if(baseQ < baseQPrev){
											
										}else{
											if(!secondEnd){
												countedReads.put(readName, r);
												
											}
										}
									}
									
								}else{
									countedReads.put(readName, r);
									
								}
							}
							
						}
						wgsIt.close();
						
						if(readNumber >= maxCov || countedReads.size()==0){
							continue;
						}
						//log.info("test" + readsNumTotal + "\t" + countedReads.size());
						double normalizedFragCov = (double)countedReads.size()/readsNumTotal;
						//System.err.println(normalizedFragCov + "\t" + countedReads.size() + "\t" + readsNumTotal);
						if(!chr.equalsIgnoreCase(prevChr)){
							refParser.setCurrentSequence(chr);
							prevChr = chr;
						}
						byte[] refBasesExt = refParser.loadFragment(end - kmerExt+1, kmerExt*2+1).getBytes();
						byte refBase = refBasesExt[kmerExt];
						
						
						
						
						//nearest cpg's distance in reference genome
						boolean cpgNegStrand = false;
						double nearestCpg = Double.NaN;
						if(BaseUtils.basesAreEqual(refBase, BaseUtilsMore.C)){
							if(i < cpgCollections.size()-2){
								if(i>=1){
									Interval downstreamCpg = cpgCollections.getIntervals().get(i+2);
									Interval upstreamCpg = cpgCollections.getIntervals().get(i-1);
									nearestCpg = Math.min(Math.abs(CcInferenceUtils.intervalDistance(downstreamCpg, cpg)),Math.abs(CcInferenceUtils.intervalDistance(upstreamCpg, cpg)));
									
								}else{
									Interval downstreamCpg = cpgCollections.getIntervals().get(i+2);
									nearestCpg = Math.abs(CcInferenceUtils.intervalDistance(downstreamCpg, cpg));
								}
							}else{
								if(i>=1){
									Interval upstreamCpg = cpgCollections.getIntervals().get(i-1);
									nearestCpg = Math.abs(CcInferenceUtils.intervalDistance(upstreamCpg, cpg));
								}else{
									
								}
							}
						}else if(BaseUtils.basesAreEqual(refBase, BaseUtilsMore.G)){
							cpgNegStrand = true;
							if(i >= 2){
								if(i < cpgCollections.size()-1){
									Interval downstreamCpg = cpgCollections.getIntervals().get(i+1);
									Interval upstreamCpg = cpgCollections.getIntervals().get(i-2);
									nearestCpg = Math.min(Math.abs(CcInferenceUtils.intervalDistance(downstreamCpg, cpg)),Math.abs(CcInferenceUtils.intervalDistance(upstreamCpg, cpg)));
									
								}else{
									Interval upstreamCpg = cpgCollections.getIntervals().get(i-2);
									nearestCpg = Math.abs(CcInferenceUtils.intervalDistance(upstreamCpg, cpg));
								}
							}else{
								if(i < cpgCollections.size()-1){
									Interval downstreamCpg = cpgCollections.getIntervals().get(i+1);
									nearestCpg = Math.abs(CcInferenceUtils.intervalDistance(downstreamCpg, cpg));
								}else{
									
								}
							}
						}
						if(cpgNegStrand){
							refBasesExt = BaseUtils.simpleReverseComplement(refBasesExt);
						}
						HashMap<String, Double> kmerMapsRef = new HashMap<String, Double>();
						for(int j = 2; j <= kmerLen; j++){
								kmerMapsRef.putAll(CcInferenceUtils.kmerFreqSearch(refBasesExt, j));
								
						}
						//overlap with feature in reference genome
						
						HashMap<String, Integer> overlapStatCollections = new HashMap<String, Integer>();
						for(String key : overlapLocStringCollections.keySet()){
							HashMap<String,IntervalTree<Integer>> tmp = overlapLocStringCollections.get(key);
							if(tmp.containsKey(cpg.getContig())){
								if(tmp.get(cpg.getContig()).minOverlapper(start, end)==null){
									overlapStatCollections.put(key, 0);
								}else{
									overlapStatCollections.put(key, 1);
								}
								
							}else{
								overlapStatCollections.put(key,0);
							}
							
						}
						//distance with feature in reference genome
						HashMap<String, Integer> distStatCollections = new HashMap<String, Integer>();
						for(String key : distantLocStringCollections.keySet()){
							IntervalTree<String> locCollections = distantLocStringCollections.get(key).get(cpg.getContig());
							
							//IntervalTree.Node<String> upstream = locCollections.max(start, end);
							//IntervalTree.Node<String> downstream = locCollections.min(start, end);
							int distanceNearest = Integer.MAX_VALUE;
							if(locCollections!=null && locCollections.size()>0){
								Iterator<Node<String>> upstreamIt = locCollections.reverseIterator(start, end);
								Iterator<Node<String>> downstreamIt = locCollections.iterator(start, end);
								if(!upstreamIt.hasNext()){
									IntervalTree.Node<String> downstream = locCollections.min(start, end);
									distanceNearest = CcInferenceUtils.intervalDistance(downstream,cpg);
									//System.err.println(downstream.toString());
								}else if(!downstreamIt.hasNext()){
									IntervalTree.Node<String> upstream = locCollections.max(start, end);
									distanceNearest = CcInferenceUtils.intervalDistance(upstream,cpg);
									//System.err.println(upstream.toString());
								}else{
									IntervalTree.Node<String> upstream = locCollections.max(start, end);
									IntervalTree.Node<String> downstream = locCollections.min(start, end);
									//System.err.println(upstream.toString());
									//System.err.println(downstream.toString());
									
									int dist1 = CcInferenceUtils.intervalDistance(upstream, cpg);
									int dist2 = CcInferenceUtils.intervalDistance(downstream, cpg);
									if(Math.abs(dist1) < Math.abs(dist2)){
										distanceNearest = dist1;
									}else{
										distanceNearest = dist2;
									}
								}
							}
							
							
							distStatCollections.put(key, distanceNearest);
						}
						//value in bed file
						HashMap<String, Double> valBedStatCollections = new HashMap<String, Double>();
						if(valueBedReaders != null){
							for(String key : valueBedReaders.keySet()){
								int range = valueBedReaders.get(key).getFirst();
								boolean mean0 = range < 0 ? true : false;
								if(mean0){
									range = 0-range;
								}
								TabixFeatureReader<BEDFeature, ?> bedReader = valueBedReaders.get(key).getSecond();
								CloseableTribbleIterator<BEDFeature> featureIt = bedReader.query(chr, (start-range <= 0 ? 1 : start-range), end+range);
								DescriptiveStatistics statFeature = new DescriptiveStatistics();
								while(featureIt.hasNext()){
									BEDFeature term = featureIt.next();
									if(!Double.isNaN(term.getScore())){
										statFeature.addValue(term.getScore());
									}
									
								}
								featureIt.close();
								if(statFeature.getN()>0){
									if(mean0){
										valBedStatCollections.put(key, (double)statFeature.getSum()/(double)(range*2+1));
									}else{
										valBedStatCollections.put(key, statFeature.getMean());
									}
									
								}else{
									valBedStatCollections.put(key, Double.NaN);
								}
								
							}
						}
						
						
						//value in wig file
						HashMap<String, Double> valWigStatCollections = new HashMap<String, Double>();
						if(valueWigReaders != null){
							for(String key : valueWigReaders.keySet()){
								int range = valueWigReaders.get(key).getFirst();
								if(range < 0){
									BigWigFileReader wigReader = valueWigReaders.get(key).getSecond();
									range = 0-range;
									SummaryStatistics statFeature = wigReader.queryStats(chr, (start-range <= 0 ? 1 : start-range), end+range);
									
									if(statFeature.getN()>0){
										valWigStatCollections.put(key, (double)statFeature.getSum()/(double)(range*2+1));
									}else{
										valWigStatCollections.put(key, Double.NaN);
									}
									
								}else{
									BigWigFileReader wigReader = valueWigReaders.get(key).getSecond();
									SummaryStatistics statFeature = wigReader.queryStats(chr, start-range, end+range);
									if(statFeature.getN()>0){
										valWigStatCollections.put(key, statFeature.getMean());
									}else{
										valWigStatCollections.put(key, Double.NaN);
									}
								}
								
								
								
							}
						}
						
						DescriptiveStatistics offSetStat = new DescriptiveStatistics();
						DescriptiveStatistics offSetRatioStat = new DescriptiveStatistics();
						DescriptiveStatistics fragLenStat = new DescriptiveStatistics();
						DescriptiveStatistics baseQStat = new DescriptiveStatistics();
						int methyReads = 0;
						int unmethyReads = 0;
						char cpgStrand = cpgNegStrand ? '-' : '+';
						int negStrandNum = 0;
						for(String readName : countedReads.keySet()){
							SAMRecord r = countedReads.get(readName);
							boolean negStrand = r.getReadNegativeStrandFlag();
							boolean secondEnd = r.getReadPairedFlag() && r.getSecondOfPairFlag();
							if(secondEnd){
								negStrand = !negStrand;
							}
							int offSet = r.getReadPositionAtReferencePosition(end)-1;
							
							byte[] bases = r.getReadBases();
							byte base = bases[offSet];
							byte[] baseQs = r.getBaseQualities();
							byte baseQ = baseQs[offSet];
							if(baseQ <= minBaseQ){
								continue;
							}
							
							int fragLen = Math.abs(r.getInferredInsertSize());
							
							if(negStrand){
								offSet = fragLen-offSet;
							}
							
							offSetStat.addValue(offSet);
							fragLenStat.addValue(fragLen);
							baseQStat.addValue(baseQ);
							offSetRatioStat.addValue((double)offSet/(double)fragLen);
							
							if(negStrand){
								negStrandNum++;
								if(BaseUtils.basesAreEqual(refBase, BaseUtilsMore.G)){
									if(BaseUtils.basesAreEqual(base, BaseUtilsMore.G)){
										methyReads++;
									}else if(BaseUtils.basesAreEqual(base, BaseUtilsMore.A)){
										unmethyReads++;
									}else{
										continue;
									}
								}else{
									continue;
								}
							}else{
								if(BaseUtils.basesAreEqual(refBase, BaseUtilsMore.C)){
									if(BaseUtils.basesAreEqual(base, BaseUtilsMore.C)){
										methyReads++;
									}else if(BaseUtils.basesAreEqual(base, BaseUtilsMore.T)){
										unmethyReads++;
									}else{
										continue;
									}
								}else{
									continue;
								}
							}
							
							
							
							
							points++;
							if(points % 1000000 == 0){
								log.info("Processing data points " + points + " ...");
							}
							
						}
						if((methyReads + unmethyReads) == 0){
							continue;
						}
						
						writer.write(chr + "\t" + (start-1) + "\t" + end + "\t" + String.format("%.3f",100*(double)methyReads/(double)(methyReads + unmethyReads)) + "\t" + String.format("%.3f",normalizedFragCov)
								+ "\t" + cpgStrand + "\t" + (methyReads+unmethyReads)  + "\t" + String.format("%.3f",100*(double)negStrandNum/(double)(countedReads.size())) 
								+ "\t" + nearestCpg);
						writer.write("\t" + String.format("%.2f", fragLenStat.getMean()) + "\t" + String.format("%.2f", fragLenStat.getStandardDeviation()));
						for(int j=10;j<=100;j+=10){
							writer.write("\t" + String.format("%.2f", fragLenStat.getPercentile(j)));
							
						}
						writer.write("\t" + String.format("%.2f", offSetStat.getMean()) + "\t" + String.format("%.2f", offSetStat.getStandardDeviation()));
						for(int j=10;j<=100;j+=10){
							writer.write("\t" + String.format("%.2f", offSetStat.getPercentile(j)));
							
						}
						writer.write("\t" + String.format("%.2f", offSetRatioStat.getMean()) + "\t" + String.format("%.2f", offSetRatioStat.getStandardDeviation()));
						for(int j=10;j<=100;j+=10){
							writer.write("\t" + String.format("%.2f", offSetRatioStat.getPercentile(j)));
							
						}
						writer.write("\t" + String.format("%.2f", baseQStat.getMean()) + "\t" + String.format("%.2f", baseQStat.getStandardDeviation()));
						for(int j=10;j<=100;j+=10){
							writer.write("\t" + String.format("%.2f", baseQStat.getPercentile(j)));
							
						}
						
						//overlap regions
						if(overlapStatCollections.size()>0){
							for(String key : overlapLocString){
								writer.write("\t" + overlapStatCollections.get(key));
							}
						}
					
						//distant regions
						if(distStatCollections.size() > 0){
							for(String key : distantLocString){
								writer.write("\t" + distStatCollections.get(key));
							}
						}
					
						//valBed regions
						if(valBedStatCollections.size()>0){
							for(String key : valueBedLocString){
								writer.write("\t" + String.format("%.3f",valBedStatCollections.get(key)));
							}
						}
					
						//valWig regions
						if(valWigStatCollections.size()>0){
							for(String key : valueWigLocString){
								writer.write("\t" + String.format("%.3f",valWigStatCollections.get(key)));
							}
						}
						//k-mer in reference genome
						if(kmerMapsRef.size()>0 ){
							for(String key : kmerCollections){
								writer.write("\t" + String.format("%.3f",kmerMapsRef.get(key)));
							}
						}

						writer.write("\n");
						
					}
					writer.close();
					output.close();
					
					
					wgsReader.close();
					refParser.close();

					
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
					
					finish();

	}
	
	
	private boolean failFlagFilter(SAMRecord r){
		return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
				|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag() || !r.getProperPairFlag();
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
		
		
		log.info("CpgLociMultiMetricsStats's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
	

}
