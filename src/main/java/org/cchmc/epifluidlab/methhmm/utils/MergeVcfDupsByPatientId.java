/**
 * MergeVcfDupsByPatientId.java
 * Mar 9, 2017
 * 11:14:48 AM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.methhmm.utils;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.exception.UnsortedFileException;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;



/**
 * when VCF contain genotype results at multiple bam files from the same individual, merge these genotype call.
 * When they agree with each other on GT, add up AD, DP, GQ, PL; when not agree with each other on GT, add up the agreed ones first, then choose the ones with largest GQ
 * sample name sheet format: 1st column: patient id;  2nd column: library id
 */
public class MergeVcfDupsByPatientId {

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "MergeVcfDupsByPatientId [opts] sampleSheet.txt input.vcf.gz outputMerged.vcf.gz";
	
	private static Logger log = Logger.getLogger(MergeVcfDupsByPatientId.class);
	private static long startTime = -1;
	private static long number = 0;
	private static VariantContextWriter writer;

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		MergeVcfDupsByPatientId mvdbp = new MergeVcfDupsByPatientId();
		BasicConfigurator.configure();
		mvdbp.doMain(args);
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
					String sampleFile = arguments.get(0);
					String inputFile = arguments.get(1);
					String outputFile = arguments.get(2);

					VCFFileReader vcfFileReader = initiate(inputFile, outputFile);	
					VCFHeader vcfHeader = vcfFileReader.getFileHeader();
					HashSet<String> sampleNameSet = new HashSet<String>(vcfHeader.getGenotypeSamples());
					log.info("Processing sample name file ... ");
					TreeMap<String,HashSet<String>> sampleNameList = new TreeMap<String,HashSet<String>>();//each patient id have sets of library name contained
					TreeMap<String,String> sampleNamePatientIdList = new TreeMap<String,String>();//each library name have its belonged patient id
					
					
						GZIPInputStream gzipInputStream = null;
						BufferedReader br1;
						if(sampleFile.endsWith(".gz")){
							gzipInputStream = new GZIPInputStream(new FileInputStream(sampleFile));
							br1 = new BufferedReader(new InputStreamReader(gzipInputStream));
							
						}else{
							br1 = new BufferedReader(new FileReader(sampleFile));
						}
							
							String line1;
							
							while( (line1 = br1.readLine()) != null){
								if(line1.startsWith("#"))
									continue;
								String[] splitLines = line1.split("\t");
								if(splitLines.length<2){
									continue;
								}
								
								String patientId = splitLines[0];
								String libraryId = splitLines[1];
								if(libraryId.equalsIgnoreCase("NA") || !sampleNameSet.contains(libraryId)){
									continue;
								}
								
								if(sampleNameList.containsKey(patientId)){
									HashSet<String> tmp = sampleNameList.get(patientId);
									tmp.add(libraryId);
									sampleNameList.put(patientId, tmp);
									sampleNamePatientIdList.put(libraryId, patientId);
								}else{
									HashSet<String> tmp = new HashSet<String>();
									tmp.add(libraryId);
									sampleNameList.put(patientId, tmp);
									sampleNamePatientIdList.put(libraryId, patientId);
								}
							}
							if(sampleFile.endsWith(".gz")){
								gzipInputStream.close();
							}
						br1.close();

						
					log.info("Processing vcf file ... ");
					Set<VCFHeaderLine> headMetaData = vcfHeader.getMetaDataInInputOrder();
					VCFHeader newHeader = new VCFHeader(headMetaData, sampleNameList.keySet());
					writer.writeHeader(newHeader);
					CloseableIterator<VariantContext> itSnp = vcfFileReader.iterator();
					
					
					while(itSnp.hasNext()){
						VariantContext vc = itSnp.next();	
						if(vc.isFiltered())
							continue;
						
						TreeMap<String, Genotype> genotypeListByPatient = new TreeMap<String, Genotype>();
						
						for(String patientId : sampleNameList.keySet()){
							HashSet<String> libraryIds = sampleNameList.get(patientId);
							
							
							Genotype gtPatientDefault = null;
							HashMap<String, Genotype> genotypeListWithinSamePatient = new HashMap<String, Genotype>();
							//maybe should hashmap to store the genotype list with the same genotype...
							for(String libraryId : libraryIds){
								if(vc.hasGenotype(libraryId)){
									Genotype gt = vc.getGenotype(libraryId);
									if(gtPatientDefault == null){
										gtPatientDefault = new GenotypeBuilder(gt).name(patientId).make();
									}
									
									if(gt.isCalled()){
										String genotypeString = gt.getGenotypeString();
										if(genotypeListWithinSamePatient.containsKey(genotypeString)){
											List<Allele> allele = gt.getAlleles();
											int[] gtAD = gt.getAD();
											int gtDP = gt.getDP();
											int[] gtPL = gt.getPL();
											int gtGQ = gt.getGQ();
											Genotype gtOld = genotypeListWithinSamePatient.get(genotypeString);
											int[] gtOldAD = gtOld.getAD();
											int gtOldDP = gtOld.getDP();
											int[] gtOldPL = gtOld.getPL();
											int gtOldGQ = gtOld.getGQ();
											for(int i = 0; i < gtAD.length; i++){
												gtAD[i]+=gtOldAD[i];
											}
											for(int i = 0; i < gtPL.length; i++){
												gtPL[i]+=gtOldPL[i];
											}
											gtDP+=gtOldDP;
											gtGQ+=gtOldGQ;
											
											GenotypeBuilder builder = new GenotypeBuilder();
											Genotype gtNew = builder.name(patientId).alleles(allele).AD(gtAD).DP(gtDP).GQ(gtGQ).PL(gtPL).make();
											genotypeListWithinSamePatient.put(genotypeString, gtNew);
											
										}else{
											genotypeListWithinSamePatient.put(genotypeString, new GenotypeBuilder(gt).name(patientId).make());
										}
									}
								}
							}
							if(gtPatientDefault == null){
								continue;
							}
							
							//select the best genotype for it:
							if(genotypeListWithinSamePatient.isEmpty()){
								genotypeListByPatient.put(patientId, gtPatientDefault);
							}else{
								LinkedList<Genotype> bestGenotype = new LinkedList<Genotype>();
								bestGenotype.add(gtPatientDefault);
								for(String genotypeString : genotypeListWithinSamePatient.keySet()){
									Genotype gtOld = genotypeListWithinSamePatient.get(genotypeString);
									Genotype gt = bestGenotype.getFirst();
										int gtDP = gt.getDP();
										int gtGQ = gt.getGQ();
										int gtOldDP = gtOld.getDP();
										int gtOldGQ = gtOld.getGQ();
										if(gtOldGQ > gtGQ){
											bestGenotype.addFirst(new GenotypeBuilder(gtOld).name(patientId).make());
										}else if(gtOldGQ == gtGQ){
											if(gtOldDP<gtDP){
												bestGenotype.addFirst(new GenotypeBuilder(gtOld).name(patientId).make());
											}
										}
									
								}
								genotypeListByPatient.put(patientId, bestGenotype.getFirst());
							}
							
							
						}
						
						if(!genotypeListByPatient.isEmpty()){
							ArrayList<Genotype> gcList = new ArrayList<Genotype>();
							Map<String, Integer> sampleNameToOffset = new HashMap<String, Integer>();
							List<String> sampleNames = new ArrayList<String>(sampleNameList.keySet());
							int j = 0;
							for(String key : sampleNameList.keySet()){
								gcList.add(genotypeListByPatient.get(key));
								sampleNameToOffset.put(key,j);
								j++;
								//if(vc.getStart() == 10247){
								//	System.err.println(genotypeListByPatient.get(key));
								//}
								//System.err.println(genotypeListByPatient.get(key));
							}
							VariantContext vcNew = new VariantContextBuilder(vc).genotypes(GenotypesContext.create(gcList, sampleNameToOffset, sampleNames)).make();
							//if(vc.getStart() == 10247){
								//System.err.println(vcNew);
							//}
							writer.add(vcNew);
							number++;
							if(number % 1000000 == 0){
								log.info("Processing VCF lines: " + number);
							}
						}
						
						
					}
					itSnp.close();
					vcfFileReader.close();
					
					finish();
	}
	
	private VCFFileReader initiate(String vcfFileName, String outputFile){
		startTime = System.currentTimeMillis();
		File vcfFile = new File(vcfFileName);
		
		File indexFile = Tribble.tabixIndexFile(vcfFile);
		if(!indexFile.exists() || !indexFile.canRead()){
			throw new UnsortedFileException(vcfFileName + " file's index " + indexFile.getName() + " does not exist, please use tabix to index it ...");
		}
		VCFFileReader reader = new VCFFileReader(vcfFile, true);
		
		VariantContextWriterBuilder vcfBuilder = new VariantContextWriterBuilder().setOutputFile(outputFile).setIndexCreator(new TabixIndexCreator(TabixFormat.VCF))
		.setReferenceDictionary(VCFFileReader.getSequenceDictionary(vcfFile)).setOption(Options.INDEX_ON_THE_FLY);
		if(outputFile.endsWith(".vcf.gz")){
			vcfBuilder.setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF);
		}else if(outputFile.endsWith(".vcf")){
			vcfBuilder.setOutputFileType(VariantContextWriterBuilder.OutputType.VCF);
		}else if(outputFile.endsWith(".bcf")){
			vcfBuilder.setOutputFileType(VariantContextWriterBuilder.OutputType.BCF);
		}else{
			throw new IllegalArgumentException("Wrong file suffix for vcf format output!");
		}
		writer = vcfBuilder.build();
		
		return reader;
	}

	private void finish(){
		writer.close();
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;
		log.info("Total VCF lines: " + number);
		log.info("MergeVcfDupsByPatientId's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}


}
