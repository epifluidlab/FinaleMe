/**
 * BasesMetricsWithMethy.java
 * Jun 6, 2017
 * 1:45:29 PM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.finaleme.utils;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.apache.commons.math3.util.Pair;
import org.apache.log4j.Logger;
import org.biojava.nbio.genome.parsers.twobit.TwoBitParser;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import ch.systemsx.cisd.base.mdarray.MDDoubleArray;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.HDF5FloatStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5IntStorageFeatures;
import ch.systemsx.cisd.hdf5.IHDF5Writer;

/**
 * convert bam files to most naive 8 bases metrics
 * samtools view -u -F 1804 -f 2 WGBS_ctDNA_27y_nonPreg_female.merged.hg19.chr22.100kGoodReads.sortByFrag.bam | bamToBed -bedpe -mate1 -i stdin | perl -ne 'chomp;@f=split "\t";use List::Util qw[min max];$c=$f[0];$s=min($f[1],$f[2],$f[4],$f[5]);$e=max($f[1],$f[2],$f[4],$f[5]);print "$f[0]\t$s\t$e\t.\t.\t$f[8]\n";' | bgzip -c > WGBS_ctDNA_27y_nonPreg_female.merged.hg19.chr22.100kGoodReads.sortByFrag.bedpe.gz
 * perl ~/compbio/code/mytools/perl/vcf2bed6plus2.strand.pl WGBS_ctDNA_27y_nonPreg_female.merged.hg19.reorder.mdups.calmd.cpg.filtered.sort.vcf CG
 */
public class BasesMetricsWithMethy {

	@Option(name="-maxCov",usage="maximum coverage used. Default: 64")
	public int maxCov = 64;
	
	@Option(name="-maxFragLen",usage="maximum fragment length used. Default: 500")
	public int maxFragLen = 500;
	
	@Option(name="-useNoChrPrefixBam",usage="use bam file with GRch37 instead of hg19 coordinate. Default: false")
	public boolean useNoChrPrefixBam = false;
	
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "BasesMetricsWithMethy [opts] hg19.2bit interval.bed[.gz] cpg_methy.strand.6plus2.bed[.gz] wgbs/wgs.bedpe[.gz] metrics.h5";
	
	private static Logger log = Logger.getLogger(CpgMultiMetricsStats.class);

	private static long startTime = -1;
	private static long points = 0;
	private IHDF5Writer writer = null; 


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		BasesMetricsWithMethy bmwm = new BasesMetricsWithMethy();
		//BasicConfigurator.configure();
		bmwm.doMain(args);
	}
	
	@SuppressWarnings("resource")
	public void doMain(String[] args)
			throws Exception {

					CmdLineParser parser = new CmdLineParser(this);
					//parser.setUsageWidth(80);
					try
					{
						if(help || args.length < 5) throw new CmdLineException(parser, USAGE, new Throwable());
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
					String intervalFile = arguments.get(1);
					String cpgMethyFile = arguments.get(2);
					String wgsBedFile = arguments.get(3);
					String outputFile = arguments.get(4);

					initiate(outputFile);	
					
					TwoBitParser refParser = new TwoBitParser(new File(refFile));
					log.info("Loading CpG interval file ... ");
					HashMap<String,IntervalTree<Pair<Integer, Double>>> cpgCollections = new HashMap<String,IntervalTree<Pair<Integer, Double>>>();
						GZIPInputStream gzipInputStream1 = null;
						BufferedReader br;
						if(cpgMethyFile.endsWith(".gz")){
							gzipInputStream1 = new GZIPInputStream(new FileInputStream(cpgMethyFile));
							br = new BufferedReader(new InputStreamReader(gzipInputStream1));
							
						}else{
							br = new BufferedReader(new FileReader(cpgMethyFile));
						}
							
							String line;
							
							while( (line = br.readLine()) != null){
								if(line.startsWith("#") || line.startsWith("track"))
									continue;
								String[] splitLines = line.split("\t");
								if(splitLines.length<8){
									continue;
								}
								String chr = splitLines[0];
								int start = Integer.parseInt(splitLines[1]);
								int end = Integer.parseInt(splitLines[2]);
								IntervalTree<Pair<Integer, Double>> tree;
								
								if(cpgCollections.containsKey(chr)){
									tree = cpgCollections.get(chr);
								}else{
									tree = new IntervalTree<Pair<Integer, Double>>();
								}
								int strand = 0;
								double methy = Double.parseDouble(splitLines[6]);
									
								if(splitLines[5].equalsIgnoreCase("-")){
										strand = 1;
								}else if(splitLines[5].equalsIgnoreCase("+")){
										strand = 0;
								}
								
								Pair<Integer, Double> value = new Pair<Integer, Double>(strand, methy);	
									
								tree.put(end, end, value);
								cpgCollections.put(chr, tree);
							}
							if(cpgMethyFile.endsWith(".gz")){
								gzipInputStream1.close();
							}
							br.close();
					
						log.info("Loading wgs/wgbs bedpe file ... ");
						HashMap<String,IntervalTree<Integer>> wgsBedCollections = new HashMap<String,IntervalTree<Integer>>();
							GZIPInputStream gzipInputStream2 = null;
							BufferedReader br2;
							if(wgsBedFile.endsWith(".gz")){
								gzipInputStream2 = new GZIPInputStream(new FileInputStream(wgsBedFile));
								br2 = new BufferedReader(new InputStreamReader(gzipInputStream2));
								
							}else{
								br2 = new BufferedReader(new FileReader(wgsBedFile));
							}
							
								
							while( (line = br2.readLine()) != null){
									if(line.startsWith("#") || line.startsWith("track"))
										continue;
									String[] splitLines = line.split("\t");
									if(splitLines.length<6){
										continue;
									}
									String chr = splitLines[0];
									int start = Integer.parseInt(splitLines[1]);
									int end = Integer.parseInt(splitLines[2]);
									if(end - start > maxFragLen){
										continue;
									}
									IntervalTree<Integer> tree;
									
									if(wgsBedCollections.containsKey(chr)){
										tree = wgsBedCollections.get(chr);
									}else{
										tree = new IntervalTree<Integer>();
									}
									int strand = 0;
										
									if(splitLines[5].equalsIgnoreCase("-")){
											strand = 1;
									}else if(splitLines[5].equalsIgnoreCase("+")){
											strand = 0;
									}	
									tree.put(start+1, end, strand);
									wgsBedCollections.put(chr, tree);
								}
								if(wgsBedFile.endsWith(".gz")){
									gzipInputStream2.close();
								}
								br2.close();
					
					log.info("Output hdf5 file by interval regions... ");
					String methyPath = "/methyInfo";
					String basePath = "/baseInfo";
					String windowCordinatePath = "/windowCordinateInfo";
					writer.float64().createMDArray(basePath, new long[]{1,maxCov, 4}, new int[]{1,maxCov,4}, HDF5FloatStorageFeatures.FLOAT_DEFLATE_DELETE);  //
					writer.float64().createArray(methyPath, 1, 1, HDF5FloatStorageFeatures.FLOAT_DEFLATE_DELETE);  //methylation level
					writer.int32().createMatrix(windowCordinatePath, 1,3, 1, 3, HDF5IntStorageFeatures.INT_DEFLATE_DELETE); //chr, start, end 
					
					GZIPInputStream gzipInputStream3 = null;
					BufferedReader br3;
					if(intervalFile.endsWith(".gz")){
						gzipInputStream3 = new GZIPInputStream(new FileInputStream(intervalFile));
						br3 = new BufferedReader(new InputStreamReader(gzipInputStream3));
						
					}else{
						br3 = new BufferedReader(new FileReader(intervalFile));
					}
					
					long lineNum=0;
					long intNum = 0;
					String pprevChr = "";
					while( (line = br3.readLine()) != null){
							if(line.startsWith("#") || line.startsWith("track"))
								continue;
							String[] splitLines = line.split("\t");
							if(splitLines.length<3){
								continue;
							}
							String chr = splitLines[0];
							int start = Integer.parseInt(splitLines[1]);
							int end = Integer.parseInt(splitLines[2]);
							int contigIndex =CcInferenceUtils.contigNameToId(chr);
							
							if(!chr.equalsIgnoreCase(pprevChr)){
								refParser.close();
								refParser.setCurrentSequence(chr);
								pprevChr = chr;
							}
							
							for(int pos = start+1; pos <= end; pos++){
								
								if(wgsBedCollections.containsKey(chr) && cpgCollections.containsKey(chr)){
									Iterator<Node<Integer>> wgsIt = wgsBedCollections.get(chr).overlappers(pos, pos);
									Iterator<Node<Pair<Integer, Double>>> cgIt = cpgCollections.get(chr).overlappers(pos, pos);
									double methy = 0.0;
									int strand = 0;
									if(cgIt != null && cgIt.hasNext()){
										Pair<Integer, Double> value = cgIt.next().getValue();
										strand = value.getFirst();
										methy = value.getSecond();
										//System.err.println(methy + "\t" + pos + "\t" + strand);
									}
									ArrayList<Integer> baseCollection = new ArrayList<Integer>();
									if(wgsIt != null ){
										while(wgsIt.hasNext()){
											baseCollection.add(wgsIt.next().getValue());
										}
									}
									//System.err.println(baseCollection.size() + "\t" + pos);
									byte refBase = refParser.loadFragment(pos-1, 1).getBytes()[0];
									MDDoubleArray mdaRefBases = new MDDoubleArray(new int[]{1,maxCov, 4});
									int i = 0;
									for(; i < baseCollection.size() & i < maxCov; i++){
										int fragStrand = baseCollection.get(i);
										double[] refVector = charToVector(refBase, fragStrand);
										for(int j = 0; j < refVector.length; j++){
											mdaRefBases.set(refVector[j], 0,i,j);
										}
										
									}
									for(; i < maxCov; i++){
										double[] refVector = new double[]{0,0,0,0};
										for(int j = 0; j < refVector.length; j++){
											mdaRefBases.set(refVector[j], 0,i,j);
										}
									}
									writer.float64().writeArrayBlockWithOffset(methyPath, new double[]{methy}, 1, lineNum);
									writer.float64().writeMDArrayBlockWithOffset(basePath, mdaRefBases, new long[]{lineNum, 0,0});
									
									if(lineNum % 1000000 == 0){
										log.info("Processing line: " + lineNum);
										
									}
									lineNum++;
									
								}else{
									continue; //if no reads for the entire chromosome.
								}
								
							}
							
							writer.int32().writeMatrixBlockWithOffset(windowCordinatePath, new int[][]{{contigIndex,start, end}}, intNum, 0);
							intNum++;
					}
					if(intervalFile.endsWith(".gz")){
						gzipInputStream3.close();
					}
					br3.close();
					refParser.closeParser();
					finish();
	}
	
	private double[] charToVector(byte base, int strand) {
		
			
		base = CcInferenceUtils.toUpperCase(base);
		
			if(strand == 0){
				switch (base) {
				case 'A':
					return new double[]{1.0,0.0,0.0,0.0};
				case 'T':
					return new double[]{0.0,1.0,0.0,0.0};
				case 'C':
					return new double[]{0.0,0.0,1.0,0.0};
				case 'G':
					return new double[]{0.0,0.0,0.0,1.0};
				default:
					return new double[]{0.25,0.25,0.25,0.25};
				}
			}else{
				switch (base) {
				case 'A':
					return new double[]{0.0,1.0,0.0,0.0};
				case 'T':
					return new double[]{1.0,0.0,0.0,0.0};
				case 'C':
					return new double[]{0.0,0.0,0.0,1.0};
				case 'G':
					return new double[]{0.0,0.0,1.0,0.0};
				default:
					return new double[]{0.25,0.25,0.25,0.25};
				}
				
			}
			
	}
	
	private void initiate(String outputFile){
		startTime = System.currentTimeMillis();
		writer = HDF5Factory.configure(outputFile).overwrite().writer();
	}

	private void finish(){
		writer.close();
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;
		
		log.info("Counted " + points + " data points in total");
		log.info("BasesMetricsWithMethy's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}

}
