/**
 * FragMethyPlusRefMethy.java
 * Aug 3, 2016
 * 2:46:09 PM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.methhmm.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.TreeMap;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.unc.genomics.io.BigWigFileReader;
import edu.unc.genomics.io.WigFileException;

/**
 *
 */
public class FragMethyPlusRefMethy {

	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;

	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-maxFragLen",usage="maximum fragment length allowed to check. Default: 500")
	public int maxFragLen = 500;

	@Option(name="-minCpg",usage="minimum number of CpG in the fragment. Default: 1")
	public int minCpg = 1;

	@Option(name="-onlyReadsNoNA",usage="only use reads without any missing value in any of the samples. Default: false")
	public boolean onlyReadsNoNA = false;
	
	
	@Option(name="-valueWigs",usage="bigwig files to check the value in these regions, -valueWigs trasckName:trackFileName. Default: null")
	public ArrayList<String> valueWigs = null;

	@Option(name="-covWigs",usage="bigwig files to check the coverage value in these regions. -valueWigs trasckName:trackFileName. Default: null")
	public ArrayList<String> covWigs = null;
	
	@Option(name="-getRidOffChrPrefix",usage="get rid of chr prefix for bam file mapped in hg19, while provided big wig file is in GRch37. Default: false")
	public boolean getRidOffChrPrefix = false;
		
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "FragMethyPlusRefMethy [opts] wgs.bam reads_tissue_of_origin_probability.txt.gz";
	
	private static Logger log = Logger.getLogger(FragMethyPlusRefMethy.class);

	private static long startTime = -1;
	private static long reads = 0;


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		FragMethyPlusRefMethy rtoo = new FragMethyPlusRefMethy();
		//BasicConfigurator.configure();
		rtoo.doMain(args);
	}
	

	public void doMain(String[] args)
			throws Exception {

					CmdLineParser parser = new CmdLineParser(this);
					//parser.setUsageWidth(80);
					try
					{
						if(help || args.length < 3) throw new CmdLineException(parser, USAGE, new Throwable());
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

					//String cpgListFile = arguments.get(0);
					String wgsBamFile = arguments.get(0);
					String detailFile = arguments.get(1);

					initiate();			
					
					//load interval files
					//log.info("Processing interval file ... ");
					SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(wgsBamFile));
					
					HashMap<String, BigWigFileReader> valueWigReaders = null;
					LinkedHashSet<String> valueWigLocString = new LinkedHashSet<String>();
					if(valueWigs != null){
						log.info("Loading value interval big wig file ... ");
						valueWigReaders =  new HashMap<String, BigWigFileReader>();
						for(String valueWigString : valueWigs){
							String[] splitStrings = valueWigString.split(":");
							if(splitStrings.length < 2){
								throw new IllegalArgumentException("need to provide trackname:trackFile for valueWigs");
							}
							String valueWigName = splitStrings[0];
							String valueRegion = splitStrings[1];
							valueWigLocString.add(valueWigName);
							valueWigReaders.put(valueWigName,new BigWigFileReader(new File(valueRegion).toPath()));
						}
						
					}
					
					HashMap<String, BigWigFileReader> covWigReaders = null;
					LinkedHashSet<String> covWigLocString = new LinkedHashSet<String>();
					if(covWigs != null){
						log.info("Loading coverage interval big wig file ... ");
						covWigReaders =  new HashMap<String, BigWigFileReader>();
						for(String covWigString : covWigs){
							String[] splitStrings = covWigString.split(":");
							if(splitStrings.length < 2){
								throw new IllegalArgumentException("need to provide trackname:trackFile for covWigs");
							}
							String covWigName = splitStrings[0];
							String covRegion = splitStrings[1];
							covWigLocString.add(covWigName);
							covWigReaders.put(covWigName,new BigWigFileReader(new File(covRegion).toPath()));
						}
						
					}
					
					log.info("Processing bam file ... ");
					FileOutputStream output = new FileOutputStream(detailFile);
					OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
					
					//header file
					writer.write("#readName\tnumCg\treadMethy");
					for(String key : valueWigReaders.keySet()){
						writer.write("\t" + key + "\t" + key + "_cov");
					}
					
					writer.write("\n");
					
						SAMRecordIterator wgsIt = wgsReader.iterator();
						//HashMap<String, SAMRecord> countedReads = new HashMap<String, SAMRecord>();
						
						long readNumber =0;
						while(wgsIt.hasNext()){
							SAMRecord r1 = wgsIt.next();
							SAMRecord r2 = wgsIt.next();
							readNumber+=2;
							if(readNumber % 10000 == 0){
								log.info("Processing reads " + readNumber + " ...");
								writer.flush();
							}
							if(!r1.getReadName().equalsIgnoreCase(r2.getReadName())){
								throw new IllegalArgumentException("bam file need to be sorted by read name: " + r1.getReadName() + " is not equal to " + r2.getReadName());
							}
							//log.info(r.getReadName() + "\t" + failFlagFilter(r));
							if(failFlagFilter(r1) || failFlagFilter(r2)){
								continue;
							}else{
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
								int bisulfiteStartPos1st = CcInferenceUtils.bisulfiteIncompleteReads(bases1, refBases1, r1.getReadNegativeStrandFlag(), r1.getReadPairedFlag() & r1.getSecondOfPairFlag()); //cell free dna in 2nd end, too many mismatches, so not use mismatches filter..
								int bisulfiteStartPos2nd = CcInferenceUtils.bisulfiteIncompleteReads(bases2, refBases2, r2.getReadNegativeStrandFlag(), r2.getReadPairedFlag() & r2.getSecondOfPairFlag());
								//log.debug(bisulfiteStartPos1st + "\t" + bisulfiteStartPos2nd);
								//System.err.println(bisulfiteStartPos1st + "\t" + bisulfiteStartPos2nd);
								if(bisulfiteStartPos1st < 0 || bisulfiteStartPos2nd <0){
									continue;
								}
								
								TreeMap<Integer, Integer> cpgMethy1 = CcInferenceUtils.readsMethyAndLocWithFivePrimeFilter(refBases1, bases1, basesQ1, 
										r1.getReadNegativeStrandFlag(), r1.getReadPairedFlag() & r1.getSecondOfPairFlag(), minBaseQ, bisulfiteStartPos1st, maxFragLen);
								TreeMap<Integer, Integer> cpgMethy2 = CcInferenceUtils.readsMethyAndLocWithFivePrimeFilter(refBases2, bases2, basesQ2, 
										r2.getReadNegativeStrandFlag(), r2.getReadPairedFlag() & r2.getSecondOfPairFlag(), minBaseQ, bisulfiteStartPos2nd, maxFragLen);
								

								
								TreeMap<Integer, Integer> cpgMethy =  new TreeMap<Integer, Integer>();
								int methyCount = 0;
								if(r1.getReadPairedFlag() & r1.getSecondOfPairFlag()){
									int pos1 = r1.getAlignmentStart();
									int pos2 = r2.getAlignmentStart();
									for(Integer offset : cpgMethy2.keySet()){
										cpgMethy.put(pos2+offset, cpgMethy2.get(offset));
										methyCount+=cpgMethy2.get(offset);
									}
									for(Integer offset : cpgMethy1.keySet()){
										if(!cpgMethy.containsKey(pos1+offset)){
											cpgMethy.put(pos1+offset, cpgMethy1.get(offset));
											methyCount+=cpgMethy1.get(offset);
										}
										
									}
								}else{
									int pos1 = r1.getAlignmentStart();
									int pos2 = r2.getAlignmentStart();
									for(Integer offset : cpgMethy1.keySet()){
										cpgMethy.put(pos1+offset, cpgMethy1.get(offset));
										methyCount+=cpgMethy1.get(offset);
									}
									for(Integer offset : cpgMethy2.keySet()){
										if(!cpgMethy.containsKey(pos2+offset)){
											cpgMethy.put(pos2+offset, cpgMethy2.get(offset));
											methyCount+=cpgMethy2.get(offset);
										}
										
									}
								}
								
								if(cpgMethy.size() < minCpg){
									continue;
								}
								
								String chr = r1.getReferenceName();
								
								
								//calculate the probaility belong to multiple tissues ...
								if(valueWigReaders != null){
									HashMap<String, Pair<Double, Integer>> probabilities = new HashMap<String, Pair<Double, Integer>>();
									boolean probNAflag = false;
									
									
											
									for(String key : valueWigReaders.keySet()){
										Pair<Double, Integer> refMethySummary = fragRefMethy(valueWigReaders.get(key), covWigReaders.get(key), cpgMethy, chr);
										//Pair<Double, Integer> refMethySummary = new Pair(1,1);
										if(refMethySummary.getSecond()<1){
											probNAflag=true;
										}
										probabilities.put(key, refMethySummary);
									}
									
									if(probNAflag && onlyReadsNoNA){
										continue;
									}
									writer.write(r1.getReadName() + "\t" + cpgMethy.size() + "\t" + (double)methyCount/(double)cpgMethy.size());
									//System.err.print(r1.getReadName() + "\t" + cpgMethy.size() + "\t" + (double)methyCount/(double)cpgMethy.size());
									for(String key : valueWigReaders.keySet()){
										Pair<Double, Integer> refMethySummary = probabilities.get(key);
										writer.write("\t" + refMethySummary.getFirst()  + "\t" + refMethySummary.getSecond());
										//System.err.print("\t" + refMethySummary.getFirst()  + "\t" + refMethySummary.getSecond());
										//System.err.print("\t" + key);
									}
									writer.write("\n");
									//System.err.println();
									reads++;
									
								}
								
							}
							
							
						}
						wgsIt.close();
						
						if(valueWigs != null){
							for(String key : valueWigReaders.keySet()){
								valueWigReaders.get(key).close();
							}
						}
						
					
						if(covWigs != null){
							for(String key : covWigReaders.keySet()){
								covWigReaders.get(key).close();
							}
						}	
					
					writer.close();
					output.close();
					wgsReader.close();

					finish();

	}
	
	
	private boolean failFlagFilter(SAMRecord r){
		return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
				|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag() || !r.getProperPairFlag();
	}
	
	private Pair<Double, Integer> fragRefMethy(BigWigFileReader reader, BigWigFileReader covReader, TreeMap<Integer, Integer> cpgMethy, String chr) throws WigFileException, IOException{

		if(getRidOffChrPrefix){
			chr = chr.replace("chr", "");
		}

		int methyCountSumReads = 0;
		int totalCountSumReads = 0;
		
		for(Integer pos : cpgMethy.keySet()){
			SummaryStatistics stat = reader.queryStats(chr, pos, pos);
			SummaryStatistics statCov = covReader.queryStats(chr, pos, pos);
			if(stat.getN()>0){
				//double methyMean = (double)stat.getMean();
				int methyCov = (int)statCov.getMean();
				//int methylatedCount = (int)(methyCov*methyMean);
				int methylatedCount = (int)stat.getMean();
				methyCountSumReads +=methylatedCount;
				totalCountSumReads +=methyCov;
				
			}	
			
				
		}	

		//double methyReads = (double)methyCountSumReads/(100*(double)totalCountSumReads);
		double methyReads = (double)methyCountSumReads/(double)totalCountSumReads;
		return new Pair<Double,Integer>(methyReads, totalCountSumReads);
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

		log.info("FragMethyPlusRefMethy's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
					
	
}
