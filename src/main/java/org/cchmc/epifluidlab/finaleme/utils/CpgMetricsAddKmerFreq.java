/**
 * CpgMetricsAddKmer.java
 * Jul 17, 2017
 * 5:07:03 PM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.finaleme.utils;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.util.Pair;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * only count CpG now...
 */
public class CpgMetricsAddKmerFreq {
	
	@Option(name="-h",usage="show option information")
	public boolean help = false;
	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "CpgMetricsAddKmerFreq [opts] input.txt[.gz] output.txt.gz";
	
	private static Logger log = Logger.getLogger(CpgMetricsAddKmerFreq.class);

	private static long startTime = -1;
	private static long points = 0;


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		CpgMetricsAddKmerFreq cmak = new CpgMetricsAddKmerFreq();
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

					String inputFile = arguments.get(0);
					String outputFile = arguments.get(1);

					initiate();			
					
					//initiate different kinds of reader
					//reference genome
					
					
					HashMap<String, Pair<Integer,Integer>> matrixProcess = new HashMap<String, Pair<Integer,Integer>>();
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

							while( (line = br.readLine()) != null){
								String[] splitLines = line.split("\t");
								if(line.startsWith("#") || splitLines[1].equalsIgnoreCase("start") || splitLines.length<3)
									continue;
								
								String key = splitLines[3];
								int fragLen = Integer.parseInt(splitLines[4]);
								if(matrixProcess.containsKey(key)){
									Pair<Integer,Integer> value = matrixProcess.get(key);
									value = new Pair<Integer,Integer>(value.getFirst()+1,fragLen);
									matrixProcess.put(key, value);
								}else{
									Pair<Integer,Integer> value = new Pair<Integer,Integer>(1,fragLen);
									matrixProcess.put(key, value);
								}

							}
							if(inputFile.endsWith(".gz")){
								gzipInputStream1.close();
							}
							br.close();
							log.info("Writing file ... ");

							OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile)), "UTF-8");
							gzipInputStream1 = null;
							if(inputFile.endsWith(".gz")){
								gzipInputStream1 = new GZIPInputStream(new FileInputStream(inputFile));
								br = new BufferedReader(new InputStreamReader(gzipInputStream1));
								
							}else{
								br = new BufferedReader(new FileReader(inputFile));
							}
								

								while( (line = br.readLine()) != null){
									String[] splitLines = line.split("\t");
									if(line.startsWith("#") || splitLines[1].equalsIgnoreCase("start") || splitLines.length<3)
										continue;
									
									String key = splitLines[3];
									if(matrixProcess.containsKey(key)){
										Pair<Integer,Integer> value = matrixProcess.get(key);
										double freq = (double)value.getFirst()/(double)value.getSecond();
										writer.write(line + "\t" + freq + "\n");
									}

									
									points++;
								}
								if(inputFile.endsWith(".gz")){
									gzipInputStream1.close();
								}
								br.close();
							
							
							writer.close();
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
		log.info("CpgMetricsAddKmerFreq's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
	
	

}
