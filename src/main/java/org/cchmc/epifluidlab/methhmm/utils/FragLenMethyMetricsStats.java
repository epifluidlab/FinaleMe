/**
 * CovFragLenMetricsStats.java
 * Feb 18, 2016
 * 11:09:54 AM
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
import java.util.List;




import java.util.TreeMap;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.MathIllegalStateException;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 *
 */
public class FragLenMethyMetricsStats {


//	@Option(name="-minCov",usage="minimum coverage for methylation calculation. default: 1")
//	public int minCov = 1;
	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;

	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-maxFragLen",usage="maximum fragment length allowed to check. Default: 10000")
	public int maxFragLen = 10000;
	
	@Option(name="-excludeFragNoCG",usage="Exclude fragment without CG in the summary file. Default: false")
	public boolean excludeFragNoCG = false;

	@Option(name="-excludeFragCHmethy",usage="Exclude fragment with CH methylation more than 0.1. methylation scale 0-1. Default: 0.1")
	public double excludeFragCHmethy = 0.1;

	@Option(name="-skipSecondEnd",usage="skip the 2nd end for the statistics. Default: false")
	public boolean skipSecondEnd = false;
	
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "FragLenMethyMetricsStats [opts] wgs.bam fragment_detail.txt.gz fragLenMethySummary.txt.gz fragMethyLenSummary.txt.gz";
	
	private static Logger log = Logger.getLogger(FragLenMethyMetricsStats.class);

	private static long startTime = -1;


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		FragLenMethyMetricsStats flmms = new FragLenMethyMetricsStats();
		BasicConfigurator.configure();
		flmms.doMain(args);
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
					String wgsBamFile = arguments.get(0);
					String detailFile = arguments.get(1);
					String  fragLenMethySummaryFile = arguments.get(2);
					String  fragMethyLenSummaryFile = arguments.get(3);
					initiate();
					
					
					
					//load interval files
					log.info("Processing bam file ... ");
					SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(wgsBamFile));
					
					SAMRecordIterator wgsIt = wgsReader.iterator();
					HashMap<String,FragMethyDatum> readsStat = new HashMap<String,FragMethyDatum>();
					long num = 0;
					
					while(wgsIt.hasNext()){
						SAMRecord r = wgsIt.next();
						if(skipSecondEnd && r.getReadPairedFlag() && r.getSecondOfPairFlag())
							continue;
						FragMethyDatum fragMethyDatum = processRead(r, wgsReader);
						String readName = r.getReadName();
						if(readsStat.containsKey(readName)){ //second pair of read, only add methylation, not fragment length
								readsStat.put(readName, readsStat.get(readName).add(fragMethyDatum));
						}else{
								readsStat.put(readName, fragMethyDatum);
						}

						
						num++;
						if(num % 1000000 == 0){
							log.info("Processing reads " + num + " ...");
						}
					}
					wgsIt.close();
					wgsReader.close();
					
					//
					FileOutputStream output = new FileOutputStream(detailFile);
					OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
					writer.write("#readName\tFragLen\tfrac_methy_CG\tmethy_CG\tunmethy_CG\tfrac_methy_CH\tmethy_CH\tunmethy_CH\n");
					
					log.info("Summarize fragment informaiton ... ");
					TreeMap<Integer,DescriptiveStatistics> fragLenMethySummary = new TreeMap<Integer,DescriptiveStatistics>();
					TreeMap<Integer,DescriptiveStatistics> fragMethyLenSummary = new TreeMap<Integer,DescriptiveStatistics>();
					for(String key : readsStat.keySet()){
						FragMethyDatum fm = readsStat.get(key);
						if(excludeFragNoCG && (fm.methy_CG + fm.unmethy_CG)==0){
							continue;
						}
						if(fm.fragLen > maxFragLen || (double)fm.methy_CH/(double)(fm.methy_CH + fm.unmethy_CH) > excludeFragCHmethy){
							continue;
						}
						writer.write(key + "\t" + fm.fragLen+ "\t" + String.format("%.3f",100*(double)fm.methy_CG/(double)(fm.methy_CG + fm.unmethy_CG))  + "\t" + fm.methy_CG + "\t" + fm.unmethy_CG + "\t" 
								 + String.format("%.3f",100*(double)fm.methy_CH/(double)(fm.methy_CH + fm.unmethy_CH))  + "\t" + fm.methy_CH + "\t" + fm.unmethy_CH + "\n");
						if(fragLenMethySummary.containsKey(fm.fragLen)){
							DescriptiveStatistics stats = fragLenMethySummary.get(fm.fragLen);
							if(!Double.isNaN((double)fm.methy_CG/(double)(fm.methy_CG + fm.unmethy_CG))){
								stats.addValue(100*(double)fm.methy_CG/(double)(fm.methy_CG + fm.unmethy_CG));
								fragLenMethySummary.put(fm.fragLen, stats);
							}
							
						}else{
							DescriptiveStatistics stats = new DescriptiveStatistics();
							if(!Double.isNaN((double)fm.methy_CG/(double)(fm.methy_CG + fm.unmethy_CG))){
								stats.addValue(100*(double)fm.methy_CG/(double)(fm.methy_CG + fm.unmethy_CG));
								fragLenMethySummary.put(fm.fragLen, stats);
							}
						}
						if(!Double.isNaN((double)fm.methy_CG/(double)(fm.methy_CG + fm.unmethy_CG))){
							int methy = (int)(10*(double)fm.methy_CG/(double)(fm.methy_CG + fm.unmethy_CG));
							if(fragMethyLenSummary.containsKey(methy)){
								DescriptiveStatistics stats = fragMethyLenSummary.get(methy);
								stats.addValue(fm.fragLen);
								fragMethyLenSummary.put(methy, stats);
							}else{
								DescriptiveStatistics stats = new DescriptiveStatistics();
								stats.addValue(fm.fragLen);
								fragMethyLenSummary.put(methy, stats);
							}
						}
						
						
					}
					
					writer.close();
					output.close();
					
					
					log.info("Output fragment summary informaiton ... ");
					FileOutputStream outputFragLenMethySummary = new FileOutputStream(fragLenMethySummaryFile);
					OutputStreamWriter writerFragLenMethySummary = new OutputStreamWriter(new GZIPOutputStream(outputFragLenMethySummary), "UTF-8");
					writerFragLenMethySummary.write("#FragLen\tFracMethy_mean\tFracMethy_sd\tFracMethy_median\n");
					for(Integer key : fragLenMethySummary.keySet()){
						DescriptiveStatistics stat = fragLenMethySummary.get(key);
						writerFragLenMethySummary.write(key + "\t" + String.format("%.3f", stat.getMean()) + "\t" + String.format("%.3f", stat.getStandardDeviation())
								+ "\t" + String.format("%.3f", stat.getPercentile(50)) + "\t" + String.format("%d", stat.getN()) + "\n");
					}
					writerFragLenMethySummary.close();
					outputFragLenMethySummary.close();
					
					FileOutputStream outputFragMethyLenSummary = new FileOutputStream(fragMethyLenSummaryFile);
					OutputStreamWriter writerFragMethyLenSummary = new OutputStreamWriter(new GZIPOutputStream(outputFragMethyLenSummary), "UTF-8");
					writerFragMethyLenSummary.write("#FragLen\tFracLen_mean\tFracLen_sd\tFracLen_median\n");
					for(Integer key : fragMethyLenSummary.keySet()){
						DescriptiveStatistics stat = fragMethyLenSummary.get(key);
						writerFragMethyLenSummary.write(key*10 + "\t" + String.format("%.3f", stat.getMean()) + "\t" + String.format("%.3f", stat.getStandardDeviation())
								+ "\t" + String.format("%.3f", stat.getPercentile(50)) + "\t" + String.format("%d", stat.getN()) + "\n");
					}
					writerFragMethyLenSummary.close();
					outputFragMethyLenSummary.close();
					
					finish();

	}
	
	//summary coverage and fragment length distribution
	private FragMethyDatum processRead(SAMRecord r, SamReader reader) throws MathIllegalStateException, MathIllegalArgumentException, IOException{
		
		Pair<Integer, Integer> readsMethyCG = CcInferenceUtils.readsMethySummary(r, minMapQ, minBaseQ, reader, false);
		Pair<Integer, Integer> readsMethyCH = CcInferenceUtils.readsMethySummaryNonCG(r, minMapQ, minBaseQ, reader, false);
		FragMethyDatum fragMethyDatum = new FragMethyDatum(readsMethyCG.getFirst(), readsMethyCG.getSecond(), readsMethyCH.getFirst(), readsMethyCH.getSecond(), Math.abs(r.getInferredInsertSize()));
		return fragMethyDatum;
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
		
		
		log.info("FragLenMethyMetricsStats's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
	
	public class FragMethyDatum{
		public long methy_CG = 0L;
		public long unmethy_CG = 0L;
		public long methy_CH = 0L;
		public long unmethy_CH = 0L;
		public int fragLen = 0;
		
		public FragMethyDatum(){
			
		}
		
		public FragMethyDatum(long methy_CG, long unmethy_CG, long methy_CH, long unmethy_CH, int fragLen){
			this.methy_CG = methy_CG;
			this.unmethy_CG = unmethy_CG ;
			this.methy_CH = methy_CH;
			this.unmethy_CH = unmethy_CH;
			this.fragLen = fragLen;
		}
		
		public FragMethyDatum add(FragMethyDatum other){ //only merge methylation, but fragment length not.
			this.methy_CG = this.methy_CG + other.methy_CG;
			this.unmethy_CG = this.unmethy_CG + other.unmethy_CG;
			this.methy_CH = this.methy_CH + other.methy_CH;
			this.unmethy_CH = this.unmethy_CH + other.unmethy_CH;
			
			return this;
		}
		
	}

}
