/**
 * FragComplexMetricsStats.java
 * Mar 2, 2016
 * 11:11:39 AM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.finaleme.utils;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.SequenceUtil;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.lang3.tuple.Triple;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 *
 */
public class FragComplexMetricsStats {

	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;

	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-maxFragLen",usage="maximum fragment length allowed to check. Default: 10000")
	public int maxFragLen = 200;

	@Option(name="-kmerLen",usage="the K-mer length to check. Default: 4")
	public int kmerLen = 4;
	
	@Option(name="-excludeFragNoCG",usage="Exclude fragment without CG in the summary file. Default: false")
	public boolean excludeFragNoCG = false;

	@Option(name="-excludeFragCHmethy",usage="Exclude fragment with CH methylation more than 0.1. methylation scale 0-1. Default: 0.1")
	public double excludeFragCHmethy = 0.1;
	
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "FragComplexMetricsStats [opts] wgs.bam(sorted by name) fragment_detail.txt.gz";
	
	private static Logger log = Logger.getLogger(FragComplexMetricsStats.class);

	private static long startTime = -1;
	private static long readPairs = 0;


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		FragComplexMetricsStats fcms = new FragComplexMetricsStats();
		BasicConfigurator.configure();
		fcms.doMain(args);
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
					String wgsBamFile = arguments.get(0);
					String detailFile = arguments.get(1);

					initiate();			
					
					//load interval files
					log.info("Processing bam file ... ");
					SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(wgsBamFile));
					FileOutputStream output = new FileOutputStream(detailFile);
					OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
					writer.write("#chr\tstart\tend\treadName\tFragLen\tFrag_strand\tfrac_methy_CG\tmethy_density_CG\tmethy_CG\tunmethy_CG\tGC_perc");
					for(int i=2; i <= kmerLen; i++){
						for(byte[] kmer : SequenceUtil.generateAllKmers(i)){
							writer.write("\t" + new String(kmer));
						}
					}
					
					writer.write("\n");
					SAMRecordIterator wgsIt = wgsReader.iterator();
					
					while(wgsIt.hasNext()){
						SAMRecord first = wgsIt.next();
						String readNameFirst = first.getReadName();
						if(!wgsIt.hasNext()){
							throw new SAMException("Read " + readNameFirst + "does not have paired read!");
						}
						SAMRecord second = wgsIt.next();
						String readNameSecond = second.getReadName();
						if(!readNameFirst.equalsIgnoreCase(readNameSecond)){
							throw new SAMFormatException("Read " + readNameFirst + "is not sorted by name!");
						}
						
						if(!processReadPair(first, second, writer)){
							continue;
						}
							
							

						
						readPairs++;
						if(readPairs % 1000000 == 0){
							log.info("Processing read pairs " + readPairs + " ...");
						}
					}
					wgsIt.close();
					wgsReader.close();
					
					//
					

					
					writer.close();
					output.close();

					
					finish();

	}
	
	private boolean processReadPair(SAMRecord first, SAMRecord second, OutputStreamWriter writer) throws Exception{
		//filter by flag, mapping quality, too many mismatches, bisulfite-incomplete conversion
		
		if(failFlagFilter(first) || failFlagFilter(second)){
			return false;
		}
		if(first.getAlignmentStart() >= second.getAlignmentStart()){ //first read is always at the upstream in reference genome
			SAMRecord r = first;
			first = second;
			second = r;
		}
		//filter by orientation, only consider when two end face with each other
		if(!CcInferenceUtils.passReadPairOrientation(first, second)){
			//log.debug(first);
			return false;
		}
		
		//filter by insertion size
		int fragLen = Math.abs(first.getInferredInsertSize());
		if(fragLen > maxFragLen || fragLen != Math.abs(second.getInferredInsertSize())){
			return false;
		}
		
		//summarize methylation information
		int bisulfiteStartPos1st = CcInferenceUtils.bisulfiteIncompleteReads(first); //cell free dna in 2nd end, too many mismatches, so not use mismatches filter..
		int bisulfiteStartPos2nd = CcInferenceUtils.bisulfiteIncompleteReads(second);
		//log.debug(bisulfiteStartPos1st + "\t" + bisulfiteStartPos2nd);
		if(bisulfiteStartPos1st < 0 || bisulfiteStartPos2nd <0){
			return false;
		}
		ArrayList<byte[]> fragmentNonOverlapInfo = CcInferenceUtils.constructNonOverlapFragment(first, second);
		// summarize fragment length info considering overlap information
		Triple<Integer, Integer, Integer[]> methySummary = CcInferenceUtils.readsMethyAndLocSummaryWithFivePrimeFilter(CcInferenceUtils.toUpperCase(fragmentNonOverlapInfo.get(0)),fragmentNonOverlapInfo.get(1),fragmentNonOverlapInfo.get(2),
															first.getReadNegativeStrandFlag(), first.getSecondOfPairFlag(), minBaseQ, bisulfiteStartPos1st, maxFragLen);
		//summarize frag_len, (strand info may affect information, +/- of 1st end), CpG number in each 10bin, G+C% in each 10bp bin, all 4mer number, mean,sd,median CpG methylation of this fragment
		//two end is from the same strand DNA fragment, so it is determined by 1st end's strand status
		//log.debug(methySummary.getFirst() + "\t" + methySummary.getSecond() + "\t" + (double)methySummary.getFirst()/(double)(methySummary.getFirst() + methySummary.getSecond()));
		int methy=methySummary.getLeft();
		int unmethy=methySummary.getMiddle();
		double methyDensity = 100*(double)methy/(double)fragLen;
		if(fragLen>(first.getReadLength() + second.getReadLength())){
			methyDensity = 100*(double)methy/(double)(first.getReadLength() + second.getReadLength());
		}
		double methyFrac = 100*(double)methy/(double)(methy + unmethy);
		Integer[] cpgLocs = methySummary.getRight();
		byte[] refBase = CcInferenceUtils.toUpperCase(fragmentNonOverlapInfo.get(0));
		char strand = '+';
		if((first.getFirstOfPairFlag() && first.getReadNegativeStrandFlag()) || (first.getSecondOfPairFlag() && !first.getReadNegativeStrandFlag())){
			refBase = BaseUtilsMore.simpleReverse(CcInferenceUtils.toUpperCase(fragmentNonOverlapInfo.get(0)));
			strand = '-';
		}
		double gcFreq = CcInferenceUtils.patternFreqSearch(refBase, "G") + CcInferenceUtils.patternFreqSearch(refBase, "C");
		//double cpgFreq = CcInferenceUtils.patternFreqSearch(refBase, "CG");
		writer.write(first.getContig() + "\t" + first.getAlignmentStart() + "\t" + second.getAlignmentEnd() + "\t" + first.getReadName()+ "\t" + fragLen + "\t" + strand  + "\t" + String.format("%.4f", methyFrac)  + "\t" + String.format("%.4f", methyDensity) + "\t" + 
				methy + "\t" + unmethy + "\t" + String.format("%.3f", 100*gcFreq));
		//add K-mer frequency (maximum 4-mer now)
		for(int i=2; i <= kmerLen; i++){
			for(byte[] kmer : SequenceUtil.generateAllKmers(i)){
				double freq = CcInferenceUtils.patternFreqSearch(refBase, new String(kmer));
				writer.write("\t" + String.format("%.3f", 100*freq));
			}
		}
		
		//for(Integer cpgLoc : cpgLocs){
		//	writer.write("\t" + cpgLoc);
		//}
		writer.write("\n");
		return true;
	}
	
	private boolean failFlagFilter(SAMRecord r){
		return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
				|| r.getReadFailsVendorQualityCheckFlag() && r.getDuplicateReadFlag() || (r.getReadPairedFlag() && !r.getProperPairFlag());
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
		
		log.info("Counted " + readPairs + " read pairs in total");
		log.info("FragComplexMetricsStats's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}

//	public class FragComplexMetrics{
		
//		public FragComplexMetrics(){
			
//		}
//	}
	

}
