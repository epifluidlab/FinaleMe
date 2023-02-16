/**
 * EstimateMaxPrediction.java
 * Mar 1, 2017
 * 9:49:47 AM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.methhmm.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.IntervalTree.Node;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;



/**
 *
 */
public class EstimateMaxPrediction {

	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;

	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-maxFragLen",usage="maximum fragment length allowed to check. Default: 500")
	public int maxFragLen = 500;
	
	@Option(name="-minDist",usage="minimum distance between two fragment. Default: 1")
	public int minDist = 1;
	
	@Option(name="-useNoChrPrefixBam",usage="use bam file with GRch37 instead of hg19 coordinate. Default: false")
	public boolean useNoChrPrefixBam = false;

	/**
	 * 
	 */
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "EstimateMaxPrediction [opts] cpg_list.bed[.gz] wgbs/wgs.bam";
	
	private static Logger log = Logger.getLogger(EstimateMaxPrediction.class);

	private static long startTime = -1;
	private static long pointTotal = 0;
	private static long pointConcord = 0;
	private static long pointMethyConcord = 0;
	private static long pointUnmethyConcord = 0;


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		EstimateMaxPrediction emp = new EstimateMaxPrediction();
		BasicConfigurator.configure();
		emp.doMain(args);
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
					String cpgListFile = arguments.get(0);
					String wgsBamFile = arguments.get(1);

					initiate();
					
					

					//load interval files
					log.info("Processing interval file ... ");
					SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(wgsBamFile));
						log.info("Processing interval file ... ");
						HashMap<String,IntervalTree<String>> cpgList = new HashMap<String,IntervalTree<String>>();

							GZIPInputStream gzipInputStream = null;
							BufferedReader br1;
							if(cpgListFile.endsWith(".gz")){
								gzipInputStream = new GZIPInputStream(new FileInputStream(cpgListFile));
								br1 = new BufferedReader(new InputStreamReader(gzipInputStream));
								
							}else{
								br1 = new BufferedReader(new FileReader(cpgListFile));
							}
								
								String line1;
								
								while( (line1 = br1.readLine()) != null){
									if(line1.startsWith("#"))
										continue;
									String[] splitLines = line1.split("\t");
									if(splitLines.length<3){
										continue;
									}
									String chr = splitLines[0];
									int start = Integer.parseInt(splitLines[1]);
									int end = Integer.parseInt(splitLines[2]);
									
									
									IntervalTree<String> tree;
									
									if(cpgList.containsKey(chr)){
										tree = cpgList.get(chr);
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
									cpgList.put(chr, tree);
								}
								if(cpgListFile.endsWith(".gz")){
									gzipInputStream.close();
								}
						br1.close();
						
						
						log.info("Processing bam file ... ");//too much time and not good when 
						
						
						
						
						
						for(String chrKey : cpgList.keySet()){
							if(chrKey.equalsIgnoreCase("chrM")){
								continue;
								
							}
							
							String bamChr = chrKey;
							if(useNoChrPrefixBam){
								
								Pattern replace = Pattern.compile("^chr");
				                Matcher matcher1 = replace.matcher(bamChr);
				                bamChr=matcher1.replaceAll("");
							}
							long prevPos = 0;
							
							
							IntervalTree<String> cpgChrCollections = cpgList.get(chrKey);
							Iterator<Node<String>> cpgIterator = cpgChrCollections.iterator();
							while(cpgIterator.hasNext()){
								Node<String> cpg = cpgIterator.next();
								int start = cpg.getStart();
								int end = cpg.getEnd();
								
								ArrayList<SAMRecord> countedReads = new ArrayList<SAMRecord>();
								//System.err.println(bamChr + "\t" + start + "\t" + end);
								SAMRecordIterator wgsIt = wgsReader.queryOverlapping(bamChr,start+1,end);

								while(wgsIt.hasNext()){
									SAMRecord r = wgsIt.next();
									if(failFlagFilter(r)){
										continue;
									}else{
										countedReads.add(r);
									}
								}
								wgsIt.close();
								if(countedReads.size() <= 2)
									continue;
								//System.err.println(start + "\t"+ end + "\t" + countedReads.size());
								for(int i = 0; i < countedReads.size(); i++){
									SAMRecord r1 = countedReads.get(i);
									int readOffset1 = r1.getReadPositionAtReferencePosition(end)-1;
									byte[] refBases1 = CcInferenceUtils.toUpperCase(SequenceUtil.makeReferenceFromAlignment(r1, false));
									//System.err.println(readOffset1 + "\t" + refBases1.length);
									if(readOffset1<0 || readOffset1 >= refBases1.length){
										continue;
									}
									byte refBase1 = refBases1[readOffset1];
									byte readBasesQ1 = r1.getBaseQualities()[readOffset1];
									byte readBases1 = r1.getReadBases()[readOffset1];
									if(readBasesQ1 < minBaseQ){
										continue;
									}
									boolean negStrand1 = r1.getReadNegativeStrandFlag();
									boolean secondPair1 = r1.getReadPairedFlag() && r1.getSecondOfPairFlag();
									if(secondPair1){
										negStrand1=!negStrand1;
									}
									int fragOffset1 = CcInferenceUtils.getFragOffsetFromReadsOffset(r1, readOffset1);
									int distFragEnd1 = CcInferenceUtils.getDistFragEndFromReadsOffset(r1, readOffset1);
									
									
									
									for(int j = i+1; j < countedReads.size(); j++){
										SAMRecord r2 = countedReads.get(j);
										boolean negStrand2 = r2.getReadNegativeStrandFlag();
										boolean secondPair2 = r2.getReadPairedFlag() && r2.getSecondOfPairFlag();
										if(secondPair2){
											negStrand2=!negStrand2;
										}
										//System.err.println(r1.getReadName() + "\t" + r2.getReadName() + "\t" + negStrand1 + '\t' + negStrand2);
										if(r2.getReadName().equalsIgnoreCase(r1.getReadName()) || negStrand1!=negStrand2){
											continue;
										}
										int readOffset2 = r2.getReadPositionAtReferencePosition(end)-1;
										//System.err.println(readOffset2 + "\t" + r2.getReadBases().length);
										if(readOffset2<0 || readOffset2 >= r2.getReadBases().length){
											continue;
										}
										
										int fragOffset2 = CcInferenceUtils.getFragOffsetFromReadsOffset(r2, readOffset2);
										int distFragEnd2 = CcInferenceUtils.getDistFragEndFromReadsOffset(r2, readOffset2);
										
										byte readBasesQ2 = r2.getBaseQualities()[readOffset2];
										byte readBases2 = r2.getReadBases()[readOffset2];
										if(readBasesQ2 < minBaseQ){
											continue;
										}
										
										//System.err.println(end + "\t" + r1.getAlignmentStart() + "\t" + r1.getAlignmentEnd() + "\t"  + r1.getMateAlignmentStart() + "\t" + readOffset1 + "\t" + fragOffset1 + "\t" + distFragEnd1 + "\t" 
										//		+ new String(refBases1) + "\t" + (char)(refBase1) + "\t" + new String(r1.getReadBases()) + "\t" + (char)readBases1 + "\t" + negStrand1 + "\t" + secondPair1);
										//System.err.println(r2.getAlignmentStart() + "\t" + r2.getAlignmentEnd() + "\t"  + r2.getMateAlignmentStart() + "\t" + readOffset2 + "\t" + fragOffset2 + "\t" + distFragEnd2 + "\t" 
										//		+ new String(r2.getReadBases()) + "\t" + (char)readBases2 + "\t" + negStrand2 + "\t" + secondPair2);
										
										//System.err.println(refBase1 + "\t"+ readBases1 + "\t" + readBases2 + "\t" + fragOffset1 + "\t" + fragOffset2 + "\t" + distFragEnd1 + "\t" + distFragEnd2 + "\t" + r1.getReadName() + "\t" + r2.getReadName());
										
										if(Math.abs(fragOffset1-fragOffset2) == minDist && Math.abs(distFragEnd1-distFragEnd2) == minDist && !(Math.abs(fragOffset1-fragOffset2) == 0 && Math.abs(distFragEnd1-distFragEnd2) == 0)){
											//System.err.println(end + "\t" + r1.getReadName() + "\t" + r1.getAlignmentStart() + "\t" + r1.getAlignmentEnd() + "\t"  + r1.getMateAlignmentStart() + "\t" + readOffset1 + "\t" + fragOffset1 + "\t" + distFragEnd1 + "\t" 
											//		+ new String(refBases1) + "\t" + (char)(refBase1) + "\t" + new String(r1.getReadBases()) + "\t" + (char)readBases1 + "\t" + negStrand1 + "\t" + secondPair1);
											//System.err.println(r2.getReadName() + "\t" + r2.getAlignmentStart() + "\t" + r2.getAlignmentEnd() + "\t"  + r2.getMateAlignmentStart() + "\t" + readOffset2 + "\t" + fragOffset2 + "\t" + distFragEnd2 + "\t" 
											//		+ new String(r2.getReadBases()) + "\t" + (char)readBases2 + "\t" + negStrand2 + "\t" + secondPair2);
											if(SequenceUtil.basesEqual(refBase1, SequenceUtil.C)){
												if(!negStrand1){
													pointTotal++;
													if(SequenceUtil.basesEqual(readBases1, readBases2)){
														pointConcord++;
														if(SequenceUtil.basesEqual(readBases1, SequenceUtil.C)){
															pointMethyConcord++;
														}else if(SequenceUtil.basesEqual(readBases1, SequenceUtil.T)){
															pointUnmethyConcord++;
														}
													}
													if(pointTotal % 10000 == 0 ){
														System.err.println("pointConcord: " + pointConcord + "\tpointMethyConcord: " + pointMethyConcord + "\tpointUnmethyConcord: " + pointUnmethyConcord + "\tpointTotal: " + pointTotal);
													}
													
												}
											}else if(SequenceUtil.basesEqual(refBase1, SequenceUtil.G)){
												if(negStrand1){
													pointTotal++;
													if(SequenceUtil.basesEqual(readBases1, readBases2)){
														pointConcord++;
														if(SequenceUtil.basesEqual(readBases1, SequenceUtil.G)){
															pointMethyConcord++;
														}else if(SequenceUtil.basesEqual(readBases1, SequenceUtil.A)){
															pointUnmethyConcord++;
														}
													}
													if(pointTotal % 10000 == 0 ){
														System.err.println("pointConcord: " + pointConcord + "\tpointMethyConcord: " + pointMethyConcord + "\tpointUnmethyConcord: " + pointUnmethyConcord + "\tpointTotal: " + pointTotal);
													}
													
												}
											}
										}
										//if(pointTotal % 100 == 0 ){
											
										//}
									}
								}
								
							}
						}

						wgsReader.close();
					
					finish();
	}
	
	private boolean failFlagFilter(SAMRecord r){
		return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
				|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag() || !r.getProperPairFlag() || !CcInferenceUtils.passReadPairOrientation(r) || r.getStringAttribute("MD") == null;
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
		
		log.info("Counted  " + pointTotal + " pairs in total");
		log.info("Counted concordant " + pointConcord + " pairs in total");
		log.info("CpgMultiMetricsToHdf5's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
	


}
