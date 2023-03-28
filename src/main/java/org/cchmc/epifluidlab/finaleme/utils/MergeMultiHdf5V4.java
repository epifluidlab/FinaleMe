/**
 * MergeMultiHdf5V4.java
 * Mar 21, 2017
 * 3:35:33 PM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.finaleme.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.HDF5FloatStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5IntStorageFeatures;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;

/**
 *
 */
public class MergeMultiHdf5V4 {

	/**
	 * 
	 */
	@Option(name="-maxFragLen",usage="maximum fragment length allowed to check. Default: 500")
	public int maxFragLen = 500;

	@Option(name="-maxCgPerFrag",usage="maximum number of CpG in each fragment for the outputy MDarray. currently 24 is the maximum for training dataset...Default: 30")
	public int maxCgPerFrag = 30;
	
	@Option(name="-blockSize",usage="specify the blockSize (row number) of alignment for HDF5 reader/writer access. default: 10000")
	public int blockSize = 10000;

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "MergeMultiHdf5V4 [opts] output.h5 input1.h5 input2.h5 ... ";
	
	private static Logger log = Logger.getLogger(MergeMultiHdf5V4.class);

	private static long startTime = -1;
	private static long pointTotal = 0;
	private IHDF5Writer writer = null; 


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		MergeMultiHdf5V4 mmh = new MergeMultiHdf5V4();
		BasicConfigurator.configure();
		mmh.doMain(args);
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

					//read input bed file, for each row,
					//String intervalFile = arguments.get(0);
					String outputFile = arguments.remove(0);
					List<String> inputFiles = arguments;

					initiate(outputFile);
					boolean writeHeader = false;
					String cpgMatrixStatsPath = "/cpgMatrixStats"; //store the methylated, unmethylated data points; then min and max of each feature... 
					String headerInfoPath = "/cpgMatrixHeaderInfo";
					String cpgMatrixPath = "/cpgMatrixInfo";
					String fragBasePath = "/fragBaseInfo";
					String refBasePath = "/refBaseInfo";
					String responsePath = "/responseInfo";
					String cpgCordinatePath = "/cpgCordinatePath"; //write down the cpg coordinate and read name belongs to ...
					int headerSize = 0;
					
					long rowTotal=0L;
					double[][] cpgMatrixStat = null;
					double pointsTotal = 0;
					//need to sort by cpg coordinate and then read name by each of the files...
					for(String inputFile : inputFiles){
						IHDF5Reader reader = HDF5Factory.openForReading(inputFile);
						long[] dims = reader.object().getDimensions(responsePath);
						long[] cpgCordinateDim = reader.object().getDimensions(cpgCordinatePath);
						long[] cpgMatrixDim = reader.object().getDimensions(cpgMatrixPath);
						long[] fragBaseDim = reader.object().getDimensions(fragBasePath);
						long[] refBaseDim = reader.object().getDimensions(refBasePath);
						long[] cpgMatrixStatsDim = reader.object().getDimensions(cpgMatrixStatsPath);
						
						if(dims[0] != cpgMatrixDim[0] || dims[0] != fragBaseDim[0] || dims[0] != refBaseDim[0] || dims[0] != cpgCordinateDim[0]){
							throw new IllegalArgumentException("Number of elements is not equal between daatsets");
						}
						
						if(!writeHeader){
							reader.object().copy(headerInfoPath, writer, headerInfoPath);
							headerSize = (int)reader.object().getDimensions(headerInfoPath)[0];
							writer.int32().createArray(responsePath, 1,1, HDF5IntStorageFeatures.INT_DEFLATE_DELETE); 
							writer.int32().createMatrix(cpgCordinatePath, 1,3, 1, 3, HDF5IntStorageFeatures.INT_DEFLATE_DELETE); 

							writer.float64().createMDArray(cpgMatrixPath, new long[]{1,maxCgPerFrag,headerSize}, new int[]{1,maxCgPerFrag,headerSize}, HDF5FloatStorageFeatures.FLOAT_DEFLATE_DELETE);
							writer.float64().createMDArray(fragBasePath, new long[]{1,maxFragLen,4}, new int[]{1,maxFragLen,4}, HDF5FloatStorageFeatures.FLOAT_DEFLATE_DELETE);
							writer.float64().createMDArray(refBasePath, new long[]{1,maxFragLen*2+1,4}, new int[]{1,maxFragLen*2+1,4}, HDF5FloatStorageFeatures.FLOAT_DEFLATE_DELETE);
							//write response result array
							writer.float64().createMatrix(cpgMatrixStatsPath, 1,4, 1, 4, HDF5FloatStorageFeatures.FLOAT_DEFLATE_DELETE); //store the methylated, unmethylated  data points, TotalFragment , 0; then min,max, mean, sd of each feature... 
							
							cpgMatrixStat = new double[(int)cpgMatrixStatsDim[0]][(int)cpgMatrixStatsDim[1]];
							for(int i = 0; i < (int)cpgMatrixStatsDim[0]; i++){
								for(int j = 0; j < (int)cpgMatrixStatsDim[1]; j++){
									if(i==0){
										cpgMatrixStat[i][j] = 0.;
									}else{
										if(j==0){
											cpgMatrixStat[i][j] = Double.MAX_VALUE;
										}else if(j==1){
											cpgMatrixStat[i][j] = Double.MIN_VALUE;
										}else{
											cpgMatrixStat[i][j] = 0.;
										}
									}
									
								}
							}
							writeHeader = true;
						}
						//need to sort it here for each files..
						//long[] sortIndex = sortGenomeLocRead(reader.int32().readMatrix(cpgCordinatePath));
						
						long rowStart = 0L;
						//for(int i = 0; i < sortIndex.length; i++){
						while(rowStart<dims[0]){
							long size = blockSize;
							if((rowStart+blockSize)>=dims[0]){
									size = dims[0]-rowStart;
							}
						//	long rowStart = sortIndex[i];
						//	long size = 1L;
							writer.int32().writeArrayBlockWithOffset(responsePath,reader.int32().readArrayBlockWithOffset(responsePath, (int)size, rowStart),(int)size,rowTotal);
							writer.int32().writeMatrixBlockWithOffset(cpgCordinatePath,reader.int32().readMatrixBlockWithOffset(cpgCordinatePath, (int)size, (int)cpgCordinateDim[1], rowStart, 0L),(int)size,(int)cpgCordinateDim[1], rowTotal, 0L);

							writer.float64().writeMDArrayBlockWithOffset(cpgMatrixPath,reader.float64().readMDArrayBlockWithOffset(cpgMatrixPath, new int[]{(int)size,(int)cpgMatrixDim[1],(int)cpgMatrixDim[2]}, new long[]{rowStart,0,0}),new long[]{rowTotal,0,0});
							writer.float64().writeMDArrayBlockWithOffset(fragBasePath,reader.float64().readMDArrayBlockWithOffset(fragBasePath, new int[]{(int)size,(int)fragBaseDim[1],(int)fragBaseDim[2]}, new long[]{rowStart,0,0}),new long[]{rowTotal,0,0});
							writer.float64().writeMDArrayBlockWithOffset(refBasePath,reader.float64().readMDArrayBlockWithOffset(refBasePath, new int[]{(int)size,(int)refBaseDim[1],(int)refBaseDim[2]}, new long[]{rowStart,0,0}),new long[]{rowTotal,0,0});
							
							if(rowTotal % 100000 == 0){
								log.info("Processing rows ... " + rowTotal);
							}
							rowStart+=blockSize;
							rowTotal += size;
						}
						double[][] cpgMatrixStatsTmp = reader.float64().readMatrix(cpgMatrixStatsPath);
						double pointSum = cpgMatrixStatsTmp[0][2] + cpgMatrixStatsTmp[0][3];
						pointsTotal += pointSum;
						for(int i = 0; i < (int)cpgMatrixStatsDim[0]; i++){
							if(i==0){
								for(int j = 0; j < (int)cpgMatrixStatsDim[1]; j++){
									cpgMatrixStat[i][j] += cpgMatrixStatsTmp[i][j];
								}
							}else{
								cpgMatrixStat[i][0] = cpgMatrixStatsTmp[i][0] < cpgMatrixStat[i][0] ? cpgMatrixStatsTmp[i][0] : cpgMatrixStat[i][0];
								cpgMatrixStat[i][1] = cpgMatrixStatsTmp[i][1] > cpgMatrixStat[i][1] ? cpgMatrixStatsTmp[i][1] : cpgMatrixStat[i][1];
								cpgMatrixStat[i][2] += cpgMatrixStatsTmp[i][2]*pointSum;
								cpgMatrixStat[i][3] += cpgMatrixStatsTmp[i][3]*cpgMatrixStatsTmp[i][3]*pointSum;
							}
						}

						reader.close();
					}
					
					writer.float64().writeMatrixBlockWithOffset(cpgMatrixStatsPath,new double[][]{cpgMatrixStat[0]},1,4, 0L, 0L);
					
					for(int i = 1; i < cpgMatrixStat.length; i++){ //min, max, mean, sd
						writer.float64().writeMatrixBlockWithOffset(cpgMatrixStatsPath,new double[][]{{cpgMatrixStat[i][0], cpgMatrixStat[i][1], cpgMatrixStat[i][2]/pointsTotal, Math.sqrt(cpgMatrixStat[i][3]/pointsTotal)}},1,4, i, 0L);
						
					}
					
					finish();
	}
	
	private long[] sortGenomeLocRead(int[][] matrix){
		TreeMap<GenomeLocRead, Long> sortMatrix = new TreeMap<GenomeLocRead, Long>();
		for(long i = 0L; i < matrix.length; i++){
			GenomeLocRead tmp = new GenomeLocRead(matrix[(int)i][0], matrix[(int)i][1], matrix[(int)i][2]);
			//System.err.println(matrix[(int)i][0] + "\t" + matrix[(int)i][1] + "\t" +  matrix[(int)i][2] + "\t" + i + tmp.hashCode());
			sortMatrix.put(tmp, i);
		}
		System.err.println(sortMatrix.size() + "\t" + matrix.length);
		long[] sortedIndexs = new long[sortMatrix.size()];
		int j = 0;
		for(GenomeLocRead g : sortMatrix.keySet()){
			if(sortMatrix.containsKey(g)){
				sortedIndexs[j] = sortMatrix.get(g);
				j++;
			}
			
		}
		System.err.println(j);
		return sortedIndexs;
	}
	
	class GenomeLocRead implements Comparable<GenomeLocRead>{
		
		public int chr;
		public int pos;
		public int readNames;
		
		public GenomeLocRead(int chr, int pos, int readNames){
			this.chr = chr;
			this.pos = pos;
			this.readNames = readNames;
			
		}


		
		@Override
		public int compareTo(GenomeLocRead o) {
			if(this.chr == o.chr){
				if(this.pos == o.pos){
					return this.readNames - o.readNames;
					
				}else{
					return this.pos - o.pos;
				}
				
			}else{
				return this.chr - o.chr;
			}

		}
		
		@Override
		public int hashCode(){
			int hashcode = 17;
			hashcode = 37 * hashcode + this.chr ;
			hashcode = 37 * hashcode + this.pos ;
			hashcode = 37 * hashcode + this.readNames;

			return hashcode;
			
		}
		
		@Override
		public boolean equals(Object o){
			if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            
			if(this.chr == ((GenomeLocRead)o).chr && this.pos == ((GenomeLocRead)o).pos && this.readNames == ((GenomeLocRead)o).readNames){
				return true;
			}else{
				return false;
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
		
		log.info("Merged" + pointTotal + " points in total");
		log.info("MergeMultiHdf5V4's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
	

}
