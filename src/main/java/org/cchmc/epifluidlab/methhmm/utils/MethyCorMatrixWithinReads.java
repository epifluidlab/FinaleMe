/**
 * MethyCorMatrixWithinReads.java
 * Apr 9, 2016
 * 9:11:45 AM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.methhmm.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
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
import org.ujmp.core.doublematrix.impl.DefaultDenseDoubleMatrix2D;


/**
 *
 */
public class MethyCorMatrixWithinReads {

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

	@Option(name="-skipSecondEnd",usage="skip the 2nd end for the statistics. Default: false")
	public boolean skipSecondEnd = false;



	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "MethyCorMatrixWithinReads [opts] locus.bed wgbs.bam cor.detail.txt.gz methy.detail.txt.gz";
	
	private static Logger log = Logger.getLogger(MethyCorMatrixWithinReads.class);

	private static long startTime = -1;
	private static long points = 1;


	
	public static void main(String[] args)
	throws Exception
	{
		MethyCorMatrixWithinReads mcmwr = new MethyCorMatrixWithinReads();
		mcmwr.doMain(args);
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
		

		String cpgListFile = arguments.get(0);
		String wgsBamFile = arguments.get(1);
		String corDetailFile = arguments.get(2);
		String methyDetailFile = arguments.get(3);
		initiate();
		
		SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(wgsBamFile));
		
		log.info("Loading locus interval file ... ");	 //CpG file need to be merged two strand file (so the distance to next CpG is not always the other strand one...)
		IntervalList cpgCollectionsUnsorted = new IntervalList(wgsReader.getFileHeader());
		BufferedReader br;
		GZIPInputStream gzipInputStream = null;
		if(cpgListFile.endsWith(".gz")){
			gzipInputStream = new GZIPInputStream(new FileInputStream(cpgListFile));
			br = new BufferedReader(new InputStreamReader(gzipInputStream));
			
		}else{
			br = new BufferedReader(new FileReader(cpgListFile));
		}
		String line;
		while( (line = br.readLine()) != null){
			if(line.startsWith("#"))
				continue;
			String[] splitLines = line.split("\t");
			if(splitLines.length<3){
				continue;
			}else{
				String chr = splitLines[0];
				int start = Integer.parseInt(splitLines[1])+1;
				int end = Integer.parseInt(splitLines[2]);
				if(splitLines.length < 6){
					cpgCollectionsUnsorted.add(new Interval(chr, start, end));
				}else{
					boolean negStrand = splitLines[5].equalsIgnoreCase("-") ? true : false;
					cpgCollectionsUnsorted.add(new Interval(chr, start, end, negStrand, splitLines[3]));
				}
			}
			
		}
		if(cpgListFile.endsWith(".gz")){
			gzipInputStream.close();
		}
		br.close();
		IntervalList cpgCollections = cpgCollectionsUnsorted.sorted();
		

		
		log.info("Output value for each Locus in each DNA fragment ... ");
		FileOutputStream outputCor = new FileOutputStream(corDetailFile);
		OutputStreamWriter writerCor = new OutputStreamWriter(new GZIPOutputStream(outputCor), "UTF-8");
		FileOutputStream outputMethy = new FileOutputStream(methyDetailFile);
		OutputStreamWriter writerMethy = new OutputStreamWriter(new GZIPOutputStream(outputMethy), "UTF-8");
		ArrayList<Double[]> dataMatrix = new ArrayList<Double[]>();

		for( int i = 0; i < cpgCollections.size(); i++){
			//log.info("start" + i);
			Interval cpg = cpgCollections.getIntervals().get(i);
			String chr = cpg.getContig();
			int start = cpg.getStart();
			int end = cpg.getEnd();
			SAMRecordIterator wgsIt = wgsReader.queryOverlapping(chr,start,end);
			ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();

			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				//log.info(r.getReadName() + "\t" + failFlagFilter(r));
				if(failFlagFilter(r)){
					continue;
				}else{
					if(stringentPaired && !CcInferenceUtils.passReadPairOrientation(r)){
						continue;
					}
					reads.add(r);
				}
			}
			wgsIt.close();

			HashSet<String> countedReads = new HashSet<String>();
			for(SAMRecord first : reads){
				if(!countedReads.contains(first.getReadName())){
					SAMRecord second =  wgsReader.queryMate(first);
					if(first.getAlignmentStart() >= second.getAlignmentStart()){ //first read is always at the upstream in reference genome
						SAMRecord r = first;
						first = second;
						second = r;
					}
					ArrayList<byte[]> fragmentNonOverlapInfo = CcInferenceUtils.constructNonOverlapFragment(first, second);
					boolean negStrand = first.getReadNegativeStrandFlag();
					boolean secondFlag = first.getSecondOfPairFlag();
					
					int bisulfiteStartPos1st = CcInferenceUtils.bisulfiteIncompleteReads(first); //cell free dna in 2nd end, too many mismatches, so not use mismatches filter..
					int bisulfiteStartPos2nd = CcInferenceUtils.bisulfiteIncompleteReads(second);
					
					if(bisulfiteStartPos1st < 0 || bisulfiteStartPos2nd < 0){
						continue;
					}
					
					if(secondFlag){
						negStrand = second.getReadNegativeStrandFlag();
						secondFlag = false;
					}
					
					int locOffsetStart = first.getReadPositionAtReferencePosition(end)-1;
					if(locOffsetStart < 0){
						locOffsetStart = end-first.getAlignmentStart();
					}
					
					Double[] methy = CcInferenceUtils.fragMethyAtEachPos(CcInferenceUtils.toUpperCase(fragmentNonOverlapInfo.get(0)),fragmentNonOverlapInfo.get(1),fragmentNonOverlapInfo.get(2),
							negStrand, minBaseQ, locOffsetStart, extend);
					dataMatrix.add(methy);
					writerMethy.write(chr + "\t" + start + "\t" + end + "\t" + (negStrand ? '-' : '+'));
					for(Double m : methy){
						writerMethy.write("\t" + m);
					}
					writerMethy.write("\n");
					countedReads.add(first.getReadName());
					points++;
					if(points % 100000 == 0){
						log.info("Processing data point " + points + " ...");
						writerMethy.flush();
						wgsReader.close();
						wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(wgsBamFile));
						
					}
				}
			}
			
			
		}
		wgsReader.close();
		writerMethy.close();
		outputMethy.close();
		
		log.info("Calculating correlation ... ");
		DefaultDenseDoubleMatrix2D doubleMatrix = listToMatrix(dataMatrix);
		double[][] pearsonScoreMatrix = doubleMatrix.corrcoef(org.ujmp.core.calculation.Calculation.Ret.NEW, true, false).toDoubleArray();
		//System.err.println(pearsonScoreMatrix.length + "\t" + pearsonScoreMatrix[0].length);	
		for(int i = 0; i < pearsonScoreMatrix.length; i++){
			//for(int j = i+1; j < pearsonScoreMatrix[0].length; j++){
				writerCor.write(pearsonScoreMatrix[0][i] + "\t");
			//}
			
		}
		writerCor.write("\n");	
		writerCor.close();
		outputCor.close();
		
		
		
		finish();
	}

		
		private boolean failFlagFilter(SAMRecord r){
			return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
					|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag()
					|| !r.getProperPairFlag()|| (skipSecondEnd && r.getReadPairedFlag() && r.getSecondOfPairFlag());
		}
		
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
		log.info("MethyCorMatrixWithinReads's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
	
}
