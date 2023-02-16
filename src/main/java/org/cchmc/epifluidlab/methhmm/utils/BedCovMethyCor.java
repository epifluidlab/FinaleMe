/**
 * BedCovMethyCor.java
 * Sep 27, 2016
 * 7:19:24 AM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.methhmm.utils;



import htsjdk.samtools.util.IntervalTree;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * in each window sizes, look at the correlation between methylation level and cfDNA's coverage 
 * output: window name: chr:start-end, cor, p-value
 */
public class BedCovMethyCor {

	
	@Option(name="-window",usage="the window size. Default: 1000")
	public int window = 1000;

	@Option(name="-step",usage="the step size. Default: 1000")
	public int step = 1000;

	@Option(name="-skipProcessFirstRow",usage="don't process first row, default: not enabled")
	public boolean skipProcessFirstRow = false;
	

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "BedCovMethyCor [opts] hg19.chrom.size input.6plus2.bed[.gz] output.methy_cov_cor.txt.gz";
	
	private static Logger log = Logger.getLogger(BedCovMethyCor.class);

	private static long startTime = -1;
	private static long lineNum=0;

	private OutputStreamWriter writer = null; 


	/**
	 * @param args
	 */
	public static void main(String[] args)
			throws Exception {
			BedCovMethyCor bcmc = new BedCovMethyCor();
				BasicConfigurator.configure();
				bcmc.doMain(args);
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
							String refSizeFile = arguments.get(0);
							String inputFile = arguments.get(1);
							String outputFile = arguments.get(2);

							initiate(outputFile);
							log.info("Parsing input bed file ...");
							
							GZIPInputStream gzipInputStream = null;
							BufferedReader br;
							if(inputFile.endsWith(".gz")){
								gzipInputStream = new GZIPInputStream(new FileInputStream(inputFile));
								br = new BufferedReader(new InputStreamReader(gzipInputStream));
								
							}else{
								br = new BufferedReader(new FileReader(inputFile));
							}
							String line;

							HashMap<String, IntervalTree<String>> cpgCollections = new HashMap<String, IntervalTree<String>>();
							while( (line = br.readLine()) != null){
								if((skipProcessFirstRow && lineNum==0) || line.startsWith("#")){
									continue;
								}else{
									String[] splitin = line.split("\t");
									String chr = splitin[0];
									int start = Integer.parseInt(splitin[1]);
									int end = Integer.parseInt(splitin[2]);
									IntervalTree<String> tree = null;
									if(cpgCollections.containsKey(chr)){
										tree = cpgCollections.get(chr);
									}else{
										tree = new IntervalTree<String>();
									}
									tree.put(start,  end,  line);
									cpgCollections.put(chr, tree);
								}

							}
							if(inputFile.endsWith(".gz")){
								gzipInputStream.close();
							}
							br.close();
							
							log.info("Processing each window ...");
							BufferedReader br1 = new BufferedReader(new FileReader(refSizeFile));
							while( (line = br1.readLine()) != null){
								String[] splitin = line.split("\t");
								String chr = splitin[0];
								int len = Integer.parseInt(splitin[1]);
								if(cpgCollections.containsKey(chr)){
									IntervalTree<String> tree = cpgCollections.get(chr);
									for(int start = 0; start < len - window; start += step){
										int end = start + window;
										ArrayList<Double> methyList = new ArrayList<Double>();
										ArrayList<Double> covList = new ArrayList<Double>();
										//ArrayList<Integer> numCList = new ArrayList<Integer>();
										//ArrayList<Integer> numTList = new ArrayList<Integer>();
										Iterator<IntervalTree.Node<String>> it = tree.overlappers(start, end);
										while(it.hasNext()){
											String[] splitit = it.next().getValue().split("\t");
											double meth = Double.parseDouble(splitit[4]);
											double cov = Double.parseDouble(splitit[6]);
											//methyList.add(meth);
											//covList.add(Double.parseDouble(splitit[6]));
											int numCT = Integer.parseInt(splitit[5]);
											int numC = (int)(meth*numCT);
											int numT = numCT - numC;
											for(int i=1; i <= numC; i++){
												methyList.add(1.0);
												covList.add(cov);
											}
											for(int i=1; i <= numT; i++){
												methyList.add(0.0);
												covList.add(cov);
											}
											//System.err.println(chr + "\t" + start + "\t" + end + "\t" + splitit[0] + "\t" + splitit[1] + "\t" + splitit[2]);
										}
										//System.err.println(methyList.size());
										//System.err.println(covList.size());
										if(methyList.size()>=3){
											double[][] data = new double[methyList.size()][2];
											for(int i=0; i < methyList.size(); i++){
												data[i][0] = methyList.get(i);
												data[i][1] = covList.get(i);
											}
											PearsonsCorrelation pc = new PearsonsCorrelation(data);
											double cor = pc.getCorrelationMatrix().getEntry(0, 1);
											if(!Double.isNaN(cor)){
												//System.err.println(methyList.size());
												//System.err.println(covList.size());
												//System.err.println(s1.getStandardDeviation());
												//System.err.println(s2.getStandardDeviation());
												//System.err.println(pc.getCorrelationMatrix());
												//System.err.println(pc.getCorrelationMatrix().getEntry(0, 1));
												//System.err.println(pc.getCorrelationPValues());
												//System.err.println(pc.getCorrelationPValues().getEntry(0, 1));
												double pvalue = pc.getCorrelationPValues().getEntry(0, 1);
												String name = chr + ":" + start + "-" + end; 
												writer.write(name + "\t" + cor + "\t" + pvalue + "\t" + methyList.size() + "\n");
											}
											
										}
										
									}
								}
							}
							br1.close();
							
							
							
							finish();
	}
	
	private void initiate(String outputFile) throws Exception{
		startTime = System.currentTimeMillis();
		writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile)), "UTF-8");

	}

	private void finish() throws Exception{
		writer.close();

		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;
		
		log.info("BedCovMethyCor's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
	

}
