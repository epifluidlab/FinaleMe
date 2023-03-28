/**
 * BamToFragBasesPair.java
 * Nov 30, 2016
 * 10:28:25 AM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.finaleme.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.SequenceUtil;
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
import java.util.HashMap;
import java.util.List;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.log4j.Logger;
import org.biojava.nbio.genome.parsers.twobit.TwoBitParser;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;



/**
 * convert paired end bam file to pair of sequence context of the fragment and its relative cpg methylation status (e.g. ACTGCGTT,  00001000) 
 * (currently, only use reference genome, filtered out CpG methylation status in dbSNP position by indexed vcf.gz/bed.gz file provided); so A,T,G,C,N 5 different characters.
 * output also chr, start(frag), end(frag), mappingQ, frag_name, +/- strand, fragmentLen, seqBases, methyBinaryVector
 * CURRENT:: each CpG point as a single data point. so multiple CpG in the same fragment will be split into multiple training point. .
 */
public class BamToFragBasesPair {

	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;

	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-maxFragLen",usage="maximum fragment length allowed to check. Default: 200")
	public int maxFragLen = 200;

	@Option(name="-minCpg",usage="minimum number of CpG in the fragment. Default: 1")
	public int minCpg = 1;
	
	@Option(name="-snpBed",usage="tabixed bed.gz files to check if the region is in dbSNP/SNP position. -valueBeds trackFileName. Default: null")
	public String snpBed = null;

	@Option(name="-wgsMode",usage="used for WGS, not bisulfite space. . Default: false")
	public boolean wgsMode = false;
	
	@Option(name="-splitCpgMode",usage="each CpG point as a single data point. so multiple CpG in the same fragment will be split into multiple training point.so (e.g. ACGGCGTT,  01001000) will be (ACGGCGTT,  01000000) and (ACGGCGTT,  00001000). Default: false")
	public boolean splitCpgMode = false;

	@Option(name="-padSeq",usage="pad sequence to the maxFragLen. Default: false")
	public boolean padSeq = false;

	
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "BamToFragBasesPair [opts] hg19.2bit cpg_list.bed[.gz] wgbs/wgs.bam fragBasePair.bed.gz";
	
	private static Logger log = Logger.getLogger(BamToFragBasesPair.class);

	private static long startTime = -1;
	private static long reads = 0;


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		BamToFragBasesPair btfbp = new BamToFragBasesPair();
		//BasicConfigurator.configure();
		btfbp.doMain(args);
	}
	

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
					TwoBitParser refParser = new TwoBitParser(new File(refFile));
					
					//load interval files
					//log.info("Processing interval file ... ");
					SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(wgsBamFile));
					
					TabixFeatureReader<BEDFeature, ?> snpReader = null;
					if(snpBed != null){
						snpReader = new TabixFeatureReader(snpBed, new BEDCodec());
					}
					
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
								if(snpReader != null){
									CloseableTribbleIterator<BEDFeature> it = snpReader.query(chr, start, end);
									if(it.hasNext()){
										//System.err.println(it.next().getEnd());
										continue;
									}
								}
								
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
					if(snpReader != null){
						snpReader.close();
					}
					
					
					
					log.info("Processing bam file ... ");
					FileOutputStream output = new FileOutputStream(detailFile);
					OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
					
					//header file
					if(splitCpgMode){
						writer.write("#chr\tstart\tend\tfragLen\tfragName\tfragStrand\tseqBases\tcpgVector\tmethyStatus\n");
					}else{
						writer.write("#chr\tstart\tend\tfragLen\tfragName\tfragStrand\tseqBases\tmethyVector\n");
					}
					
					
						SAMRecordIterator wgsIt = wgsReader.iterator();
						//HashMap<String, SAMRecord> countedReads = new HashMap<String, SAMRecord>();
						
						long readNumber =0;
						String prevChr = "";
						while(wgsIt.hasNext()){
							SAMRecord r1 = wgsIt.next();
							SAMRecord r2 = wgsIt.next();
							String chr1 = r1.getContig();
							int start1 = r1.getAlignmentStart();
							int end1 = r1.getAlignmentEnd();
							int start2 = r2.getAlignmentStart();
							int end2 = r2.getAlignmentEnd();
							int fragMostLeft = start1 < start2 ? start1 : start2;
							int fragMostRight = end1 < end2 ? end2 : end1;
							//another way to do the same thing..
							//Iterator<Node<String>> it1 = null;
							//Iterator<Node<String>> it2 = null;
							
							if(cpgList.containsKey(chr1)){
								IntervalTree<String> cg = cpgList.get(chr1);
								
								if(cg.minOverlapper(start1, end1)==null && cg.minOverlapper(start2, end2)==null){
									continue;
								}else{
									//if(cg.minOverlapper(start1, end1)!=null){
									//	it1 = cg.overlappers(start1, end1);
									//}
									//if(cg.minOverlapper(start2, end2)!=null){
									//	it2 = cg.overlappers(start2, end2);
									//}
								}
							}else{
								continue;
							}
							
							readNumber+=2;
							if(readNumber % 1000000 == 0){
								log.info("Processing reads " + readNumber + " ...");
								writer.flush();
							}
							int fragLen1 = Math.abs(r1.getInferredInsertSize());
							int fragLen2 = Math.abs(r1.getInferredInsertSize());
							if(fragLen1 != fragLen2 || (fragMostRight-fragMostLeft+1) > maxFragLen){
								continue;
							}
							if(!r1.getReadName().equalsIgnoreCase(r2.getReadName())){
								throw new IllegalArgumentException("bam file need to be sorted by read name: " + r1.getReadName() + " is not equal to " + r2.getReadName());
							}
							//log.info(r.getReadName() + "\t" + failFlagFilter(r));
							if(failFlagFilter(r1) || failFlagFilter(r2)){
								continue;
							}else{
								boolean negStrand = r1.getReadNegativeStrandFlag();
								boolean secondEnd = r1.getReadPairedFlag() && r1.getSecondOfPairFlag();
								if(secondEnd){
									negStrand = !negStrand;
								}
								if(!chr1.equalsIgnoreCase(prevChr)){
									refParser.close();
									refParser.setCurrentSequence(chr1);
									prevChr = chr1;
								}
								//get reference genome in the fragment
								//System.err.println(fragMostLeft + "\t" + fragMostRight + "\t" + (fragMostRight-fragMostLeft+1) + "\t" + r1.getReadName());
								String refBasesFrag = refParser.loadFragment(fragMostLeft+1, fragMostRight-fragMostLeft+1);
								if(negStrand){
									refBasesFrag = SequenceUtil.reverseComplement(refBasesFrag);
								}
								
								//get methylation vector
								//HashMap<Integer,Integer> cpgStat = new HashMap<Integer,Integer>();
								//if(it1 != null){
								//	while(it1.hasNext()){
								//		Node<String> cg = it1.next();
								//		
								//	}
								//}
								
								
								byte[] refBases1 = CcInferenceUtils.toUpperCase(CcInferenceUtils.modifyRefSeqByCigar(CcInferenceUtils.refStrFromMd(r1), r1.getCigarString()));
								byte[] basesQ1 = CcInferenceUtils.getClippedReadsBaseQuality(r1);
								byte[] bases1 = CcInferenceUtils.toUpperCase(CcInferenceUtils.getClippedReadsBase(r1));
								//if(r1.getReadNegativeStrandFlag()){
									//bases1 = CcInferenceUtils.complementArray(bases1);
								//}
								byte[] refBases2 = CcInferenceUtils.toUpperCase(CcInferenceUtils.modifyRefSeqByCigar(CcInferenceUtils.refStrFromMd(r2), r2.getCigarString()));
								byte[] basesQ2 = CcInferenceUtils.getClippedReadsBaseQuality(r2);
								byte[] bases2 = CcInferenceUtils.toUpperCase(CcInferenceUtils.getClippedReadsBase(r2));
								//if(r2.getReadNegativeStrandFlag()){
								//	bases2 = CcInferenceUtils.complementArray(bases2);
								//}
								int bisulfiteStartPos1st = 0;
								int bisulfiteStartPos2nd = 0;
								if(!wgsMode){
									bisulfiteStartPos1st = CcInferenceUtils.bisulfiteIncompleteReads(bases1, refBases1, r1.getReadNegativeStrandFlag(), r1.getReadPairedFlag() & r1.getSecondOfPairFlag()); //cell free dna in 2nd end, too many mismatches, so not use mismatches filter..
									bisulfiteStartPos2nd = CcInferenceUtils.bisulfiteIncompleteReads(bases2, refBases2, r2.getReadNegativeStrandFlag(), r2.getReadPairedFlag() & r2.getSecondOfPairFlag());
									//log.debug(bisulfiteStartPos1st + "\t" + bisulfiteStartPos2nd);
									//System.err.println(bisulfiteStartPos1st + "\t" + bisulfiteStartPos2nd);
									if(bisulfiteStartPos1st < 0 || bisulfiteStartPos2nd <0){
										continue;
									}
								}
								
								
								
								TreeMap<Integer, Integer> cpgMethy1 = CcInferenceUtils.readsMethyAndLocWithFivePrimeFilterNoAdjustNegStrand(refBases1, bases1, basesQ1, 
										r1.getReadNegativeStrandFlag(), r1.getReadPairedFlag() & r1.getSecondOfPairFlag(), minBaseQ, bisulfiteStartPos1st, maxFragLen);
								TreeMap<Integer, Integer> cpgMethy2 = CcInferenceUtils.readsMethyAndLocWithFivePrimeFilterNoAdjustNegStrand(refBases2, bases2, basesQ2, 
										r2.getReadNegativeStrandFlag(), r2.getReadPairedFlag() & r2.getSecondOfPairFlag(), minBaseQ, bisulfiteStartPos2nd, maxFragLen);
								

								
								TreeMap<Integer, Integer> cpgMethy =  new TreeMap<Integer, Integer>();
								
								if(r1.getReadPairedFlag() & r1.getSecondOfPairFlag()){
									int pos1 = r1.getAlignmentStart();
									int pos2 = r2.getAlignmentStart();
									for(Integer offset : cpgMethy2.keySet()){
										cpgMethy.put(pos2+offset, cpgMethy2.get(offset));
										
									}
									for(Integer offset : cpgMethy1.keySet()){
										if(!cpgMethy.containsKey(pos1+offset)){
											cpgMethy.put(pos1+offset, cpgMethy1.get(offset));
											
										}
										
									}
								}else{
									int pos1 = r1.getAlignmentStart();
									int pos2 = r2.getAlignmentStart();
									for(Integer offset : cpgMethy1.keySet()){
										cpgMethy.put(pos1+offset, cpgMethy1.get(offset));
										
									}
									for(Integer offset : cpgMethy2.keySet()){
										if(!cpgMethy.containsKey(pos2+offset)){
											cpgMethy.put(pos2+offset, cpgMethy2.get(offset));
											
										}
										
									}
								}
								
								if(cpgMethy.size() < minCpg){
									continue;
								}
								
								
								
								 
								if(splitCpgMode){// output also chr, start(cpg), end(cpg), fragmentLen, frag_name, +/- strand, seqBases, methyBinaryVector, cpg_methy_status
									for(int cgCor : cpgMethy.keySet()){
										int methyStatus = cpgMethy.get(cgCor);
										String methyVector = "";
										if(negStrand){
											int j=0;
											for(int i = fragMostRight; i >= fragMostLeft; i--,j++){
												if(cgCor == (i+2)){
													methyVector += "2";
												}else{
													methyVector += "1";
												}
											}
											if(padSeq){
												for(;j<maxFragLen; j++){
													methyVector += "0";
												}
											}
											
										}else{
											int j=0;
											for(int i = fragMostLeft; i <= fragMostRight; i++,j++){
												if(cgCor == (i+2)){
													methyVector += "2";
												}else{
													methyVector += "1";
												}
											}
											if(padSeq){
												for(;j<maxFragLen; j++){
													methyVector += "0";
												}
											}
											
										}
										writer.write(chr1 + "\t" + (cgCor-1) + "\t" + cgCor + "\t" + fragLen1 + "\t" + r1.getReadName() + "\t" + (negStrand ? "-" : "+") + "\t" + (padSeq ? CcInferenceUtils.charToVector(refBasesFrag, maxFragLen) : CcInferenceUtils.charToVector(refBasesFrag)) + "\t" + methyVector + "\t" + methyStatus + "\n");
									}
								}else{// output also chr, start(frag), end(frag), fragmentLen, frag_name, +/- strand, seqBases, methyBinaryVector
									String methyVector = "";
									if(negStrand){
										int j=0;
										for(int i = fragMostRight; i >= fragMostLeft; i--,j++){
											if(cpgMethy.containsKey(i+2)){
												methyVector += cpgMethy.get(i+2)+1;
											}else{
												methyVector += 1;
											}
										}
										if(padSeq){
											for(;j<maxFragLen; j++){
												methyVector += "0";
											}
										}
										
									}else{
										int j=0;
										for(int i = fragMostLeft; i <= fragMostRight; i++,j++){
											if(cpgMethy.containsKey(i+2)){
												methyVector += cpgMethy.get(i+2)+1;
											}else{
												methyVector += 1;
											}
										}
										if(padSeq){
											for(;j<maxFragLen; j++){
												methyVector += "0";
											}
										}
										
									}
									writer.write(chr1 + "\t" + fragMostLeft + "\t" + fragMostRight + "\t" + fragLen1 + "\t" + r1.getReadName() + "\t" + (negStrand ? "-" : "+") + "\t" + (padSeq ? CcInferenceUtils.charToVector(refBasesFrag, maxFragLen) : CcInferenceUtils.charToVector(refBasesFrag)) + "\t" + methyVector + "\n");
								}
								
								
									reads++;
								
							}
							
							
						}
						wgsIt.close();
						refParser.close();
					
					writer.close();
					output.close();
					wgsReader.close();

					finish();

	}
	
	
	
	
	private boolean failFlagFilter(SAMRecord r){
		return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
				|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag() || !r.getProperPairFlag() || !CcInferenceUtils.passReadPairOrientation(r);
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
		
		log.info("Counted " + reads + " informative reads in total");

		log.info("BamToFragBasesPair's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
		

}
