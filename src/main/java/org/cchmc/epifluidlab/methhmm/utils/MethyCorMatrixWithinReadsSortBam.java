/**
 * MethyCorMatrixWithinReads.java
 * Apr 9, 2016
 * 9:11:45 AM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.methhmm.utils;

import htsjdk.samtools.*;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;


import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;

/**
 *
 */
public class MethyCorMatrixWithinReadsSortBam {

	/**
	 * @param args
	 */

	
	@Option(name="-extend",usage="summary plot region, so the correlation matrix will be calculated from 1100+500 to 1300+1500, for DKO1 paper F2 means -100(-100 to +400) to +100 (+100 to +500) of TSS. default:500")
	public int extend = 500;
	
	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;

	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-stringentPaired",usage="Only use paired end reads that faced to each other. Default: false")
	public boolean stringentPaired = false;




	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "MethyCorMatrixWithinReadsSortBam [opts] wgbs.sortByName.bam cor.detail.txt.gz";
	
	private static Logger log = Logger.getLogger(MethyCorMatrixWithinReadsSortBam.class);

	private static long startTime = -1;
	private static long points = 1;


	
	public static void main(String[] args)
	throws Exception
	{
		MethyCorMatrixWithinReadsSortBam mcmwr = new MethyCorMatrixWithinReadsSortBam();
		mcmwr.doMain(args);
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
		

		String wgsBamFile = arguments.get(0);
		String corDetailFile = arguments.get(1);
		//String corSummaryFile = arguments.get(2);
		initiate();
		
		SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(wgsBamFile));
		SAMRecordIterator wgsIt = wgsReader.iterator();
		/*

		HashMap<String, Pair<SAMRecord,SAMRecord>> fragments = new HashMap<String, Pair<SAMRecord,SAMRecord>>();
		log.info("Reading DNA fragments ... ");
		long numberOfFrag = 0;
		while(wgsIt.hasNext()){
			SAMRecord r = wgsIt.next();
			//log.info(r.getReadName() + "\t" + failFlagFilter(r));
			if(failFlagFilter(r)){
				continue;
			}else{
				if(stringentPaired && !CcInferenceUtils.passReadPairOrientation(r)){
					continue;
				}

				String readName = r.getReadName();
				if(fragments.containsKey(readName)){
					Pair<SAMRecord,SAMRecord> readPair = fragments.get(readName);
					if(r.getFirstOfPairFlag() && ((SAMRecord)readPair.getLeft()).getSecondOfPairFlag()){
						readPair = new Pair<SAMRecord,SAMRecord>(r,readPair.getLeft());
						fragments.put(readName,readPair);
					}else if(r.getSecondOfPairFlag() && ((SAMRecord)readPair.getLeft()).getFirstOfPairFlag()){
						readPair = new Pair<SAMRecord,SAMRecord>(readPair.getLeft(),r);
						fragments.put(readName,readPair);
					}else{
						fragments.remove(readName);
					}

				}else{
					Pair<SAMRecord,SAMRecord> readPair = new Pair<SAMRecord,SAMRecord>(r,r);
					fragments.put(readName,readPair);
				}

			}
			numberOfFrag++;
			if(numberOfFrag % 10000000 == 0){
				log.info("Importing reads: " + numberOfFrag + " ...");
		}
		wgsIt.close();
		wgsReader.close();
*/
		log.info("Output value for each Locus in each DNA fragment ... ");
		FileOutputStream outputCor = new FileOutputStream(corDetailFile);
		OutputStreamWriter writerCor = new OutputStreamWriter(new GZIPOutputStream(outputCor), "UTF-8");
		//FileOutputStream output = new FileOutputStream(corSummaryFile);
		//OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
		//ArrayList<Double[]> dataMatrix = new ArrayList<Double[]>();
		while(wgsIt.hasNext()){
			SAMRecord first = wgsIt.next();
			if(wgsIt.hasNext()){
				SAMRecord second = wgsIt.next();
				while(!first.getReadName().equalsIgnoreCase(second.getReadName())){
					if(wgsIt.hasNext()){
						first = second;
						second = wgsIt.next();
					}else{
						break;
					}
				}
				if(failFlagFilter(first) || failFlagFilter(second)){
					continue;
				}else{
					if(stringentPaired && !CcInferenceUtils.passReadPairOrientation(first)){
						//System.err.println(!CcInferenceUtils.passReadPairOrientation(first));
						continue;
					}
					if(first.getAlignmentStart() >= second.getAlignmentStart()){ //first read is always at the upstream in reference genome
						SAMRecord r = first;
						first = second;
						second = r;
					}
					ArrayList<byte[]> fragmentNonOverlapInfo = CcInferenceUtils.constructNonOverlapFragment(first, second);
					boolean negStrand = first.getReadNegativeStrandFlag();
					boolean secondFlag = first.getSecondOfPairFlag();

					//int bisulfiteStartPos1st = CcInferenceUtils.bisulfiteIncompleteReads(first); //cell free dna in 2nd end, too many mismatches, so not use mismatches filter..
					//int bisulfiteStartPos2nd = CcInferenceUtils.bisulfiteIncompleteReads(second);

					//if(bisulfiteStartPos1st < 0 || bisulfiteStartPos2nd < 0){
						//System.err.println(bisulfiteStartPos1st + "\t" + bisulfiteStartPos2nd);
					//	continue;
					//}

					if(secondFlag){
						negStrand = second.getReadNegativeStrandFlag();
					}
					//System.err.println(CcInferenceUtils.toUpperCase(fragmentNonOverlapInfo.get(0)));
					//System.err.println(fragmentNonOverlapInfo.get(1));
					//System.err.println(fragmentNonOverlapInfo.get(2));
					Double[] methy = CcInferenceUtils.fragMethyAtEachPos(CcInferenceUtils.toUpperCase(fragmentNonOverlapInfo.get(0)),fragmentNonOverlapInfo.get(1),fragmentNonOverlapInfo.get(2),
							negStrand, minBaseQ, 0, extend);

					writerCor.write(first.getReadName());
					//System.err.print(first.getReadName());
					for(int i = 0; i < methy.length; i++){
						//for(int j = i+1; j < pearsonScoreMatrix[0].length; j++){
						writerCor.write("\t" + methy[i]);
						//System.err.print("\t" + methy[i]);
						//}

					}
					writerCor.write("\n");
					//System.err.print("\n");
					//dataMatrix.add(methy);
					points++;
					if(points % 1000000 == 0){
						log.info("Processing data point " + points + " ...");
						writerCor.flush();
						//wgsReader.close();
						//wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(wgsBamFile));

					}
				}
			}
		}
		wgsIt.close();
		wgsReader.close();
		//writerCor.write("\n");
		writerCor.close();
		outputCor.close();
/*
		log.info("Calculating correlation ... ");
		DefaultDenseDoubleMatrix2D doubleMatrix = listToMatrix(dataMatrix);
		double[][] pearsonScoreMatrix = doubleMatrix.corrcoef(org.ujmp.core.calculation.Calculation.Ret.NEW, true, false).toDoubleArray();
		//System.err.println(pearsonScoreMatrix.length + "\t" + pearsonScoreMatrix[0].length);	
		for(int i = 0; i < pearsonScoreMatrix.length; i++){
			//for(int j = i+1; j < pearsonScoreMatrix[0].length; j++){
				writer.write(pearsonScoreMatrix[0][i] + "\t");
			//}
			
		}


		writer.close();
		output.close();
	*/
		
		
		finish();
	}

		
		private boolean failFlagFilter(SAMRecord r){
			return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
					|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag()
					|| !r.getProperPairFlag();
		}
/*
	private DefaultDenseDoubleMatrix2D listToMatrix(ArrayList<Double[]> dataMatrix){
		DefaultDenseDoubleMatrix2D doubleMatrix = new DefaultDenseDoubleMatrix2D(dataMatrix.size(), extend);
		for(int i = 0; i < dataMatrix.size(); i++){
			Double[] t = dataMatrix.get(i);
			for(int j = 0; j < t.length; j++){
				doubleMatrix.setDouble(t[j], i, j);
				
			}
				
			
		}
		
		return doubleMatrix;
	}
*/

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
		log.info("MethyCorMatrixWithinReadsSortBam's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
	
}
