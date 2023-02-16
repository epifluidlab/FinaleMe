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
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;




import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.MathIllegalStateException;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 *
 */
public class CovFragLenWindowMetricsStats {

	@Option(name="-values",usage="files provide methylation/GC_percent profiles in the same window. should be gziped and tabix indexed. or bigwig file. default: null")
	public ArrayList<String> values;
	@Option(name="-scale",usage="scale the coverage. Usually the total number of used reads. default option just use the raw reads count. default: -1")
	public long scale = -1;

//	@Option(name="-minCov",usage="minimum coverage for methylation calculation. default: 1")
//	public int minCov = 1;

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "CovFragLenMetricsStats [opts] interval.bed wgs.bam summary.txt.gz";
	
	private static Logger log = Logger.getLogger(CovFragLenWindowMetricsStats.class);

	private static long startTime = -1;


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		CovFragLenWindowMetricsStats cfms = new CovFragLenWindowMetricsStats();
		BasicConfigurator.configure();
		cfms.doMain(args);
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
					String intervalFile = arguments.get(0);
					String wgsBamFile = arguments.get(1);
					String summaryFile = arguments.get(2);
					initiate();
					FileOutputStream output = new FileOutputStream(summaryFile);
					OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
					
					writer.write("#chr\tstart\tend\tnum_frag\tmean_frag_len\tsd_frag_len\tmedian_frag_len");
					ArrayList<TabixFeatureReader<BEDFeature, ?>> bedReaders = new ArrayList<TabixFeatureReader<BEDFeature, ?>>();
					if(values != null){
						for(int i=1, j=0; i<= values.size(); i++, j++){
							writer.write("\tValues" + i +"_data_points\tValues" + i +"_mean\tValues" + i +"_sd\tValues" + i +"_median");
							bedReaders.add(new TabixFeatureReader(values.get(j), new BEDCodec()));
						}
					}
					
					writer.write("\n");
					
					//load interval files
					log.info("Processing interval file ... ");
					SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(wgsBamFile));
					
					BufferedReader br = new BufferedReader(new FileReader(intervalFile));
					String line;
					long lineNum = 0;
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
						
						processBam(writer, chr, start, end, wgsReader);
						
						if(values != null){
							for(int i=0; i< values.size(); i++){
								processBed(writer, chr, start, end, bedReaders.get(i));
							}
							
						}
						writer.write("\n");
						
						lineNum++;
						if(lineNum % 10000 == 0){
							log.info("Processing " + lineNum + " lines ... "); 
							writer.flush();
						}
					}
					br.close();
					wgsReader.close();
					if(values != null){
						for(int i=0; i< values.size(); i++){
							bedReaders.get(i).close();
						}
					}
					//
					
					
					
					writer.close();
					output.close();
					finish();

	}
	
	//summary coverage and fragment length distribution
	private void processBam(OutputStreamWriter writer, String chr, int start, int end, SamReader wgsReader) throws MathIllegalStateException, MathIllegalArgumentException, IOException{
		SAMRecordIterator wgsIt = wgsReader.queryOverlapping(chr,start,end);
		HashSet<String> countedReads = new HashSet<String>();
		DescriptiveStatistics statFrag = new DescriptiveStatistics();
		while(wgsIt.hasNext()){
			SAMRecord r = wgsIt.next();
			if( r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() 
						|| r.getReadFailsVendorQualityCheckFlag() && r.getDuplicateReadFlag() || (r.getReadPairedFlag() && !r.getProperPairFlag())){
				continue;
			}else{
				if(countedReads.contains(r.getReadName())){
					continue;
				}
				statFrag.addValue(Math.abs(r.getInferredInsertSize()));
				
				//if(r.getReadPairedFlag()){
				//	countedReads.add(r.g);
				//}
				countedReads.add(r.getReadName());
			}
			
		}
		wgsIt.close();
		double coverage = statFrag.getN();
		if(scale > 0){
			coverage =(double)coverage/(double)scale;
		}
		writer.write(chr + "\t" + (start-1) + "\t" + end + "\t" +  coverage);
		if(statFrag.getN()>0){
			writer.write("\t" + String.format("%.3f", statFrag.getMean()) + "\t" + String.format("%.3f", statFrag.getStandardDeviation()) + "\t" + String.format("%.3f", statFrag.getPercentile(50)));
		}else{
			writer.write("\tNA\tNA\tNA");
		}
		
	}
	
	//summary coverage and fragment length distribution
	private void processBed(OutputStreamWriter writer, String chr, int start, int end, TabixFeatureReader<BEDFeature, ?> bedReader) throws IOException{
		CloseableTribbleIterator<BEDFeature> featureIt = bedReader.query(chr, start, end);
		DescriptiveStatistics statFeature = new DescriptiveStatistics();
		while(featureIt.hasNext()){
			BEDFeature term = featureIt.next();
			statFeature.addValue(term.getScore());
		}
		featureIt.close();
		
		if(statFeature.getN()>0){
			writer.write("\t"+ statFeature.getN() + "\t" + String.format("%.3f", statFeature.getMean()) + "\t" + String.format("%.3f", statFeature.getStandardDeviation()) + "\t" + String.format("%.3f", statFeature.getPercentile(50)));
		}else{
			writer.write("\t0\tNA\tNA\tNA");
		}
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
		
		
		log.info("CovFragLenMetricsStats's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}

}
