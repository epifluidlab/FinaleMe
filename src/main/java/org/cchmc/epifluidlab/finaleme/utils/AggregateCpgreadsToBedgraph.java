/**
 * AggregateCpgreadsToBedgraph.java
 * Mar 11, 2016
 * 1:50:54 PM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.finaleme.utils;

import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.unc.genomics.BedEntry;



/**
 *
 */
public class AggregateCpgreadsToBedgraph {

	
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	final private static String USAGE = "AggregateCpgreadsToBedgraph [opts] output.bedgraph.gz cpg.bedgraph[.gz] value1.bedgraph.gz";

	@Argument
	private List<String> arguments = new ArrayList<String>();


	private static Logger log = Logger.getLogger(AggregateCpgreadsToBedgraph.class);
	private static long startTime = -1;


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		AggregateCpgreadsToBedgraph actb = new AggregateCpgreadsToBedgraph();
		
		actb.doMain(args);

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
					
					initiate();
					String outputFile = arguments.get(0);
					String cpgFile = arguments.get(1);
					String valueFile = arguments.get(2);
					//parse interval file
					log.info("Parsing input interval's feature file, duplicated lines will be automately removed in the final results ...");
					HashSet<BedEntry> inputIntervalSet = new HashSet<BedEntry>();
					BufferedReader br;
					GZIPInputStream gzipInputStream = null;
					if(cpgFile.endsWith(".gz")){
						gzipInputStream = new GZIPInputStream(new FileInputStream(cpgFile));
						br = new BufferedReader(new InputStreamReader(gzipInputStream));
						
					}else{
						br = new BufferedReader(new FileReader(cpgFile));
					}
					
					String line;

					while( (line = br.readLine()) != null){
						if(line.startsWith("#"))
							continue;
						BedEntry interval = BedEntry.parse(line);
						inputIntervalSet.add(interval);
						
					}
					if(cpgFile.endsWith(".gz")){
						gzipInputStream.close();
					}
					br.close();
					//IntervalFileReader<? extends Interval> reader = autodetectFileType(new File(intervalFile).toPath());
					//inputIntervalSet.addAll(reader.loadAll());
					
					log.info(inputIntervalSet.size() + " unique line in total at feature file ..."); 
					
					
						TabixFeatureReader<BEDFeature, ?> bf = new TabixFeatureReader(valueFile, new BEDCodec());

						
						FileOutputStream output = new FileOutputStream(outputFile);
						OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
						
						long lineNum=0;
						for(BedEntry entry : inputIntervalSet){
							CloseableTribbleIterator<BEDFeature> featureItBg = bf.query(entry.getChr(), entry.low(), entry.high());
							int numFeaturesMethy = 0;
							int numFeaturesUnmethy = 0;
							while(featureItBg.hasNext()){
								BEDFeature t = featureItBg.next();
								if(t.getName().equalsIgnoreCase("m")){
									numFeaturesMethy++;
								}else if(t.getName().equalsIgnoreCase("u")){
									numFeaturesUnmethy++;
								}
								
							}
							
							if((numFeaturesMethy + numFeaturesUnmethy)==0){
								//writer.write(entry.getChr() + "\t" + (entry.low()-1) + "\t" + entry.high() + "\t" + Double.NaN);
								continue;
							}else{
								writer.write(entry.getChr() + "\t" + (entry.low()-1) + "\t" + entry.high() + "\t" + (double)numFeaturesMethy/(double)(numFeaturesMethy+numFeaturesUnmethy) + 
										"\t" + (numFeaturesMethy+numFeaturesUnmethy));
							}
							writer.write("\n");
							//divided by background file's density
							lineNum++;
							if(lineNum % 1000000 == 0){
								log.info("Processing " + lineNum + " lines in file " + valueFile + " ... "); 
							}
						}
						bf.close();
						writer.close();
						output.close();
					
					finish();
	}
	
	

	
	
	private void initiate() throws IOException{
		startTime = System.currentTimeMillis();
		
	}

	private void finish() throws IOException{
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;
		
		log.info("AlignValueToFeature's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}

}
