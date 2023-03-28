/**
 * CcinferenceMetricsFromBed.java
 * Jun 24, 2017
 * 8:51:06 AM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.finaleme.utils;

import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.tribble.annotation.Strand;
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
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.unc.genomics.io.BigWigFileReader;

/**
 *
 */
public class CcinferenceMetricsFromBed {

	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-maxFragLen",usage="maximum fragment length allowed to check. Default: 500")
	public int maxFragLen = 500;

	@Option(name="-totalFragments",usage="total number of fragments used to normalize coverage column. default estimate from wgs/wgbs bed file by program. Default: -1")
	public long totalFragments = -1;
	
	@Option(name="-maxCov",usage="maximum coverage allowed to check. Default: 250")
	public int maxCov = 250;
	
	@Option(name="-valueWig",usage="bigwig files to annotate the regions.need to be consistent with cpg_list.bed's coordinate system. Default: null")
	public String valueWig = null;
	
	@Option(name="-useNoChrPrefixBam",usage="use wgs/wgbs.bed.gz file with GRch37 instead of hg19 coordinate. Default: false")
	public boolean useNoChrPrefixBam = false;
	
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "CcinferenceMetricsFromBed [opts] cpg_list.hg19.bed[.gz] wgs/wgbs.hg19.bed.gz cpg_detail.txt.gz";
	
	private static Logger log = Logger.getLogger(CcinferenceMetricsFromBed.class);

	private static long startTime = -1;
	private static long points = 0;


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		CcinferenceMetricsFromBed cmfb = new CcinferenceMetricsFromBed();
		//BasicConfigurator.configure();
		cmfb.doMain(args);
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
					String cpgListFile = arguments.get(0);
					String wgsBedFile = arguments.get(1);
					String detailFile = arguments.get(2);

					initiate();		
					TabixFeatureReader<BEDFeature, ?> wgsReader = new TabixFeatureReader(wgsBedFile, new BEDCodec());
					
					double readsNumTotal = 0;
					if(totalFragments > 0){
						log.info("Get total reads number used for scaling from input option -totalReadsInBam ... ");
						readsNumTotal = totalFragments;
					}else{
						log.info("Get total reads number used for scaling from bam file... ");
						CloseableTribbleIterator<BEDFeature> wgsIt = wgsReader.iterator();
						
						while(wgsIt.hasNext()){
							BEDFeature r = wgsIt.next();
							if(r.getScore() < minMapQ || r.getEnd() - r.getStart() > maxFragLen){
								continue;
							}
							readsNumTotal++;
						}
						wgsIt.close();
					}
					
					log.info((long)readsNumTotal + " fragments in total ...");
					readsNumTotal = readsNumTotal/1000000;
					
					log.info("Output value for each CpG in each DNA fragment ... ");
					
					FileOutputStream output = new FileOutputStream(detailFile);
					OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
					
					//basic
					writer.write("chr\tstart\tend\treadName\tFragLen\tFrag_strand\tmethy_stat\tNorm_Frag_cov\tbaseQ\tOffset_frag\tDist_frag_end");
					
					BigWigFileReader wigReader = null;
					if(valueWig != null){
						wigReader = new BigWigFileReader(new File(valueWig).toPath());
						writer.write("\tcgmethy");
					}
					writer.write("\n");
					
					GZIPInputStream gzipInputStream = null;
					BufferedReader br1;
					if(cpgListFile.endsWith(".gz")){
						gzipInputStream = new GZIPInputStream(new FileInputStream(cpgListFile));
						br1 = new BufferedReader(new InputStreamReader(gzipInputStream));
						
					}else{
						br1 = new BufferedReader(new FileReader(cpgListFile));
					}
						
						String line1;
						long i = 0;
						while( (line1 = br1.readLine()) != null){
							if(line1.startsWith("#"))
								continue;
							String[] splitLines = line1.split("\t");
							if(splitLines.length<6){
								continue;
							}
							String chr = splitLines[0];
							int start = Integer.parseInt(splitLines[1]);
							int end = Integer.parseInt(splitLines[2]);
							String strand = splitLines[5];
							
							double methyPrior = Double.NaN;
							if(wigReader != null){
								SummaryStatistics statFeature = wigReader.queryStats(chr, start, end);
								if(statFeature.getN()>0){
									methyPrior = statFeature.getMean();
								}else{
									methyPrior = Double.NaN;
								}
							}
							//System.err.println(chr + "\t" + start + "\t" + end + "\t" + methyPrior);
							String bamChr = chr;
							if(useNoChrPrefixBam){
								
								Pattern replace = Pattern.compile("^chr");
								Matcher matcher1 = replace.matcher(bamChr);
								bamChr=matcher1.replaceAll("");
							}
							CloseableTribbleIterator<BEDFeature> fragIt = wgsReader.query(bamChr, start+1, end);
							
							ArrayList<BEDFeature> fragments = new ArrayList<BEDFeature>();
							while(fragIt.hasNext()){
								BEDFeature r = fragIt.next();
								if(r.getScore() < minMapQ || r.getEnd() - r.getStart() > maxFragLen){
									continue;
								}
								fragments.add(r);
							}
							fragIt.close();
							if(fragments.size() > maxCov || fragments.size() == 0){
								continue;
							}
							char methyStat = 'm';
							int baseQ = 40;
							double normalizedFragCov = (double)fragments.size()/(double)readsNumTotal;
							//System.err.println(normalizedFragCov);
							for(BEDFeature fragment : fragments){
								String fragStrand = "+";
								if(fragment.getStrand().equals(Strand.NEGATIVE)){
									fragStrand = "-";
								}
								
								//System.err.println(fragStrand + "\t" + strand);
								if(strand.equalsIgnoreCase(fragStrand)){
									int fragLen = fragment.getEnd() - fragment.getStart() + 1;
									String readName = fragment.getName();
									int offSet = end - fragment.getStart() + 1;
									if(fragStrand.equalsIgnoreCase("-")){
										offSet = fragment.getEnd() - end + 1;
									}
									int distToFragEnd = end - fragment.getStart() + 1;
									int theOtherDist = fragment.getEnd() - end + 1;
									if(theOtherDist < distToFragEnd){
										distToFragEnd = theOtherDist;
									}
									writer.write(chr + "\t" + start + "\t" + end + "\t" + readName + "\t" + fragLen + "\t" + fragStrand + "\t" + methyStat + "\t" + normalizedFragCov
											 + "\t" + (int)baseQ + "\t" + offSet + "\t" + distToFragEnd + "\t" + methyPrior + "\n");
									//System.err.println(chr + "\t" + start + "\t" + end + "\t" + readName + "\t" + fragLen + "\t" + fragStrand + "\t" + methyStat + "\t" + normalizedFragCov
									//		 + "\t" + (int)baseQ + "\t" + offSet + "\t" + distToFragEnd + "\t" + methyPrior);
									points++;
								}
								
								
							}
							i++;
							if(i % 1000 == 0){
								log.info("Processing Cpg " + i + " ...");
								writer.flush();
								//if(i % 10000 == 0){
								//	wgsReader.close();
								//	wgsReader = new TabixFeatureReader(wgsBedFile, new BEDCodec());
								//}
							}
						}
						if(cpgListFile.endsWith(".gz")){
							gzipInputStream.close();
						}
						br1.close();
						if(valueWig != null){
							wigReader.close();
						}
						wgsReader.close();
						writer.close();
						output.close();
					finish();
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
		log.info("CcinferenceMetricsFromBed's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
	

}
