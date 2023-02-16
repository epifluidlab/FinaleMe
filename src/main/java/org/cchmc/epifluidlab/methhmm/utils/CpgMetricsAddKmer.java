/**
 * CpgMetricsAddKmer.java
 * Jul 17, 2017
 * 5:07:03 PM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.methhmm.utils;

import htsjdk.samtools.util.SequenceUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.log4j.Logger;
import org.biojava.nbio.genome.parsers.twobit.TwoBitParser;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 *
 */
public class CpgMetricsAddKmer {


	@Option(name="-kmerLen",usage="the K-mer length to check. Default: 4")
	public int kmerLen = 4;
	
	@Option(name="-addChr",usage="if the ref file are hg system and different from input.txt.gz file. Default: false")
	public boolean addChr = false;
	
	@Option(name="-h",usage="show option information")
	public boolean help = false;
	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "CpgMetricsAddKmer [opts] hg19.2bit input.txt[.gz] output.txt.gz";
	
	private static Logger log = Logger.getLogger(CpgMetricsAddKmer.class);

	private static long startTime = -1;
	private static long points = 0;


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		CpgMetricsAddKmer cmak = new CpgMetricsAddKmer();
		//BasicConfigurator.configure();
		cmak.doMain(args);
	}
	
	@SuppressWarnings("resource")
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

					String refFile = arguments.get(0);
					String inputFile = arguments.get(1);
					String outputFile = arguments.get(2);

					initiate();			
					
					//initiate different kinds of reader
					//reference genome
					TwoBitParser refParser = new TwoBitParser(new File(refFile));
					OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile)), "UTF-8");
					
					//basic
					
					log.info("Loading file ... ");
					
						GZIPInputStream gzipInputStream1 = null;
						BufferedReader br;
						if(inputFile.endsWith(".gz")){
							gzipInputStream1 = new GZIPInputStream(new FileInputStream(inputFile));
							br = new BufferedReader(new InputStreamReader(gzipInputStream1));
							
						}else{
							br = new BufferedReader(new FileReader(inputFile));
						}
							
							String line;
							String prevChr = "";
							while( (line = br.readLine()) != null){
								String[] splitLines = line.split("\t");
								if(line.startsWith("#") || splitLines[1].equalsIgnoreCase("start") || splitLines.length<3)
									continue;
								
								String chr = splitLines[0];
								if(addChr){
									chr = "chr" + chr;
								}
								int start = Integer.parseInt(splitLines[1]);
								int end = Integer.parseInt(splitLines[2]);
								String strand = splitLines[5];
								if(!chr.equalsIgnoreCase(prevChr)){
									refParser.close();
									refParser.setCurrentSequence(chr);
									prevChr = chr;
								}
								byte[] refBasesExt = CcInferenceUtils.toUpperCase(refParser.loadFragment(end-1-kmerLen, kmerLen*2+1).getBytes());
								
								if(strand.equalsIgnoreCase("-")){
									SequenceUtil.reverseComplement(refBasesExt);
								}
								writer.write(line + "\t" + new String(refBasesExt) + "\n");
								points++;
							}
							if(inputFile.endsWith(".gz")){
								gzipInputStream1.close();
							}
							br.close();
							writer.close();

							refParser.closeParser();
					
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
		log.info("CpgMetricsAddKmer's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
	
	

}
