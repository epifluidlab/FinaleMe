
package org.cchmc.epifluidlab.finaleme.utils;


import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;






















import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.gatk.utils.BaseUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CigarUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.SequenceUtil;



/**
 *
 */
public class CcInferenceUtils implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 6787330726376812239L;

	
	static public String iupacCodeComplement(String bases) {
		byte[] tmp = bases.getBytes();
		byte[] revTmp = new byte[tmp.length];
		for(int i = 0; i < tmp.length; i++){
			revTmp[i] = iupacCodeComplement(tmp[i]);
		}
		return new String(revTmp);
	}

	static public byte[] toUpperCase(byte[] b) {
		byte[] upperCases = new byte[b.length];
		for(int i=0;i<b.length; i++){
			upperCases[i] = toUpperCase(b[i]);
		}
		return upperCases;
	}
	
	static public byte toUpperCase(byte b) {
		return (byte) Character.toUpperCase((char) b);
	}
	
	// convert IUPAC code pattern to its complement code (in reverse strand)
	static public byte iupacCodeComplement(byte base) {
		base = toUpperCase(base);
		switch (base) {
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		case 'R':
			return 'Y';
		case 'Y':
			return 'R';
		case 'S':
			return 'S';
		case 'W':
			return 'W';
		case 'K':
			return 'M';
		case 'M':
			return 'K';
		case 'B':
			return 'V';
		case 'H':
			return 'D';
		case 'D':
			return 'H';
		case 'V':
			return 'B';
		default:
			return base;
		}
	}
	
	static public int byte2num(byte base) {
		base = toUpperCase(base);
		switch (base) {
		case 'A':
			return 1;
		case 'T':
			return 2;
		case 'C':
			return 3;
		case 'G':
			return 4;
		default:
			return 5;
		}
	}
	
	static public String charToVector(String bases) {
		byte[] tmp = bases.getBytes();
		String vectors = "";
		for(int i = 0; i < tmp.length; i++){
			vectors += byte2num(tmp[i]);
		}
		return vectors;
	}

	static public String charToVector(String bases, int maxLen) {
		byte[] tmp = bases.getBytes();
		String vectors = "";
		int j=0;
		for(int i = 0; i < tmp.length; i++, j++){
			vectors += byte2num(tmp[i]);
		}
		for(;j<maxLen; j++){
			vectors += 0;
		}
		return vectors;
	}

	
	/**
	 * It is a mismatch in Bisulfite space or not
	 * 
	 * @param refBase
	 * @param base
	 * @param negStrand
	 * @param secondPair
	 * @return matchOrNotInBisulfiteSpace
	 */
	static public boolean baseMatchInBisulfiteSpace(byte refBase, byte base, boolean negStrand, boolean secondPair) {
		if (!negStrand) {
			if (secondPair) {
				if (SequenceUtil.basesEqual(refBase, SequenceUtil.G) && SequenceUtil.basesEqual(base, SequenceUtil.A)) {
					return true;
				} else {
					return SequenceUtil.basesEqual(refBase, base);
				}
			} else {
				if (SequenceUtil.basesEqual(refBase, SequenceUtil.C) && SequenceUtil.basesEqual(base, SequenceUtil.T)) {
					return true;
				} else {
					return SequenceUtil.basesEqual(refBase, base);
				}
			}

		} else {
			if (secondPair) {
				if (SequenceUtil.basesEqual(refBase, SequenceUtil.C) && SequenceUtil.basesEqual(base, SequenceUtil.T)) {
					return true;
				} else {
					return SequenceUtil.basesEqual(refBase, base);
				}
			} else {
				if (SequenceUtil.basesEqual(refBase, SequenceUtil.G) && SequenceUtil.basesEqual(base, SequenceUtil.A)) {
					return true;
				} else {
					return SequenceUtil.basesEqual(refBase, base);
				}
			}

		}
	}
	
	static public Pair<Integer, Integer> readsMethySummary(byte[] refBases, byte[] bases, boolean negStrand, boolean secondPair){
		
		int numC = 0;
		int numT = 0;
		if((!secondPair && !negStrand) || (secondPair && negStrand)){
			for(int i = 0; i < refBases.length -1; i ++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.C) && SequenceUtil.basesEqual(refBases[i+1], SequenceUtil.G)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.C)){
						numC++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.T)){
						numT++;
					}
				}
			}
		}else{
			for(int i = 1; i < refBases.length; i ++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.G) && SequenceUtil.basesEqual(refBases[i-1], SequenceUtil.C)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.G)){
						numC++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.A)){
						numT++;
					}
				}
			}
		}
		return new Pair<Integer, Integer>(numC,numT); //first is C, second is T
	}
	
	static public Pair<Integer, Integer> readsMethySummary(byte[] refBases, byte[] bases, byte[] basesQ, boolean negStrand, boolean secondPair, int minBaseQ){
		
		int numC = 0;
		int numT = 0;
		if((!secondPair && !negStrand) || (secondPair && negStrand)){
			for(int i = 0; i < refBases.length -1; i ++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.C) && SequenceUtil.basesEqual(refBases[i+1], SequenceUtil.G)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.C) && basesQ[i] >= minBaseQ){
						numC++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.T) && basesQ[i] >= minBaseQ){
						numT++;
					}
				}
			}
		}else{
			for(int i = 1; i < refBases.length; i ++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.G) && SequenceUtil.basesEqual(refBases[i-1], SequenceUtil.C)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.G) && basesQ[i] >= minBaseQ){
						numC++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.A) && basesQ[i] >= minBaseQ){
						numT++;
					}
				}
			}
		}
		return new Pair<Integer, Integer>(numC,numT); //first is C, second is T
	}
	
	static public Pair<Integer, Integer> readsMethySummary(SAMRecord read, int minMapQ, int minBaseQ, SamReader reader, boolean useMate){
		boolean negStrand = read.getReadNegativeStrandFlag();
		boolean secondPair = read.getReadPairedFlag() && read.getSecondOfPairFlag();
		//filter if it is bad reads
		if(read.getReadUnmappedFlag() || read.getNotPrimaryAlignmentFlag() || read.getDuplicateReadFlag() 
				|| read.getReadFailsVendorQualityCheckFlag() || (read.getReadPairedFlag() && !read.getProperPairFlag()) 
				|| read.getMappingQuality() < minMapQ){
			return new Pair<Integer, Integer>(0,0);
		}
		
		byte[] refBases = toUpperCase(modifyRefSeqByCigar(SequenceUtil.makeReferenceFromAlignment(read, false), read.getCigarString()));
		byte[] readBasesQ = getClippedReadsBaseQuality(read);
		byte[] readBases = toUpperCase(getClippedReadsBase(read));
		if(negStrand){
			readBases = complementArray(readBases);
		}
		Pair<Integer, Integer> ctSummary = readsMethySummary(refBases, readBases, readBasesQ, negStrand, secondPair, minBaseQ);
		
		
		
		if(read.getReadPairedFlag() && useMate){
			SAMRecord m = reader.queryMate(read);
			byte[] refBasesMate = toUpperCase(modifyRefSeqByCigar(SequenceUtil.makeReferenceFromAlignment(m, false), m.getCigarString()));
			byte[] basesMate = toUpperCase(getClippedReadsBase(m));
			byte[] basesQMate = getClippedReadsBaseQuality(read);
			boolean negStrandMate = m.getReadNegativeStrandFlag();
			if(negStrandMate){
				basesMate = complementArray(basesMate);
			}
			boolean secondPairMate = m.getReadPairedFlag() && m.getSecondOfPairFlag();
			Pair<Integer, Integer> ctSummaryMate = readsMethySummary(refBasesMate, basesMate, basesQMate, negStrandMate, secondPairMate, minBaseQ);
			ctSummary = new Pair<Integer, Integer>(ctSummary.getFirst()+ctSummaryMate.getFirst(), ctSummary.getSecond()+ctSummaryMate.getSecond());
		}
		
		return ctSummary; //first is C, second is T
	}

	static public Pair<Integer, Integer> readsMethySummaryWithFivePrimeFilter(byte[] refBases, byte[] bases, byte[] basesQ, boolean negStrand, boolean secondPair, int minBaseQ, int bisulfiteConvertStart){
		
		int numC = 0;
		int numT = 0;
		if((!secondPair && !negStrand) || (secondPair && negStrand)){
			for(int i = bisulfiteConvertStart; i < refBases.length -1; i++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.C) && SequenceUtil.basesEqual(refBases[i+1], SequenceUtil.G)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.C) && basesQ[i] >= minBaseQ){
						numC++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.T) && basesQ[i] >= minBaseQ){
						numT++;
					}
				}
			}
		}else{
			for(int i = 1; i < refBases.length-bisulfiteConvertStart; i++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.G) && SequenceUtil.basesEqual(refBases[i-1], SequenceUtil.C)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.G) && basesQ[i] >= minBaseQ){
						numC++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.A) && basesQ[i] >= minBaseQ){
						numT++;
					}
				}
			}
		}
		return new Pair<Integer, Integer>(numC,numT); //first is C, second is T
	}
	
	static public Triple<Integer, Integer, Integer[]> readsMethyAndLocSummaryWithFivePrimeFilter(byte[] refBases, byte[] bases, byte[] basesQ, 
			boolean negStrand, boolean secondPair, int minBaseQ, int bisulfiteConvertStart, int maxFragLen){
		Integer[] cpgLoc = new Integer[maxFragLen];
		int numC = 0;
		int numT = 0;
		int middle = (int)refBases.length/2;
		HashSet<Integer> cpgLocationSets = new HashSet<Integer>();
		
		if((!secondPair && !negStrand) || (secondPair && negStrand)){
			for(int i = bisulfiteConvertStart; i < refBases.length -1; i++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.C) && SequenceUtil.basesEqual(refBases[i+1], SequenceUtil.G)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.C) && basesQ[i] >= minBaseQ){
						numC++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.T) && basesQ[i] >= minBaseQ){
						numT++;
					}
					cpgLocationSets.add(i-middle);
				}
			}
		}else{
			for(int i = 1; i < refBases.length-bisulfiteConvertStart; i++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.G) && SequenceUtil.basesEqual(refBases[i-1], SequenceUtil.C)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.G) && basesQ[i] >= minBaseQ){
						numC++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.A) && basesQ[i] >= minBaseQ){
						numT++;
					}
					cpgLocationSets.add(i-middle);
				}
				
			}
		}
		for(int j = 0-maxFragLen/2, z=0; z < maxFragLen; j++,z++){
			if(Math.abs(j)>=middle){
				cpgLoc[z]=-1;
			}else{
				if(cpgLocationSets.contains(j)){//cpg location
					cpgLoc[z]=1;
				}else{
					cpgLoc[z]=0;
				}
			}
		}
		//System.err.println(new String(refBases));
		//for(Integer t : cpgLocationSets){
		//	System.err.print("\t" + t);
		//}
		//System.err.println();
		//for(Integer t : cpgLoc){
		//	System.err.print(t);
		//}
		//System.err.println();
		return Triple.of(numC,numT, cpgLoc); //first is C, second is T, third is cpg location relative to the center of fragment
	}
	
	
	
	static public TreeMap<Integer, Integer> readsMethyAndLocWithFivePrimeFilter(byte[] refBases, byte[] bases, byte[] basesQ, 
			boolean negStrand, boolean secondPair, int minBaseQ, int bisulfiteConvertStart, int maxFragLen){
		//ArrayList<Interval> cpgLoc = new ArrayList<Interval>();
		//ArrayList<Integer> cpgMethy = new ArrayList<Integer>();
		//new Interval(chr, start, end);
		TreeMap<Integer, Integer> cpgMethy = new TreeMap<Integer, Integer>(); //key: relative location from beginning of fragment, value: 0:unmethy, 1: methy

		
		if((!secondPair && !negStrand) || (secondPair && negStrand)){
			for(int i = bisulfiteConvertStart; i < refBases.length -1; i++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.C) && SequenceUtil.basesEqual(refBases[i+1], SequenceUtil.G)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.C) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i, 1);
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.T) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i, 0);
					}
				}
			}
		}else{
			for(int i = 1; i < refBases.length-bisulfiteConvertStart; i++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.G) && SequenceUtil.basesEqual(refBases[i-1], SequenceUtil.C)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.G) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i-1, 1);
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.A) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i-1, 0);
					}
				}
				
			}
		}
		
		
		return cpgMethy; //key: relative location from beginning of fragment, value: 0:unmethy, 1: methy
	}
	
	static public TreeMap<Integer, Integer> readsMethyAndLocWithFivePrimeFilterNoAdjustNegStrand(byte[] refBases, byte[] bases, byte[] basesQ, 
			boolean negStrand, boolean secondPair, int minBaseQ, int bisulfiteConvertStart, int maxFragLen){
		//ArrayList<Interval> cpgLoc = new ArrayList<Interval>();
		//ArrayList<Integer> cpgMethy = new ArrayList<Integer>();
		//new Interval(chr, start, end);
		TreeMap<Integer, Integer> cpgMethy = new TreeMap<Integer, Integer>(); //key: relative location from beginning of fragment, value: 0:unmethy, 1: methy

		
		if((!secondPair && !negStrand) || (secondPair && negStrand)){
			for(int i = bisulfiteConvertStart; i < refBases.length -1; i++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.C) && SequenceUtil.basesEqual(refBases[i+1], SequenceUtil.G)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.C) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i, 1);
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.T) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i, 0);
					}
				}
			}
		}else{
			for(int i = 1; i < refBases.length-bisulfiteConvertStart; i++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.G) && SequenceUtil.basesEqual(refBases[i-1], SequenceUtil.C)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.G) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i, 1);
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.A) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i, 0);
					}
				}
				
			}
		}
		
		
		return cpgMethy; //key: relative location from beginning of fragment, value: 0:unmethy, 1: methy
	}
	
	static public TreeMap<Integer, Pair<Integer, Double>> readsMethyBaseqAndLocWithFivePrimeFilterNoAdjustNegStrand(byte[] refBases, byte[] bases, byte[] basesQ, 
			boolean negStrand, boolean secondPair, int minBaseQ, int bisulfiteConvertStart, int maxFragLen){

		TreeMap<Integer, Pair<Integer, Double>> cpgMethy = new TreeMap<Integer, Pair<Integer, Double>>(); //key: relative location from beginning of fragment, value: 0:unmethy, 1: methy; 2nd one is for baseQ (transform back pred scale to 0-1 scale).

		
		if((!secondPair && !negStrand) || (secondPair && negStrand)){
			for(int i = bisulfiteConvertStart; i < refBases.length -1; i++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.C) && SequenceUtil.basesEqual(refBases[i+1], SequenceUtil.G)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.C) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i, new Pair<Integer, Double>(1, phredScaleToDouble(basesQ[i])));
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.T) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i, new Pair<Integer, Double>(0, phredScaleToDouble(basesQ[i])));
					}
				}
			}
		}else{
			for(int i = 1; i < refBases.length-bisulfiteConvertStart; i++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.G) && SequenceUtil.basesEqual(refBases[i-1], SequenceUtil.C)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.G) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i, new Pair<Integer, Double>(1, phredScaleToDouble(basesQ[i])));
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.A) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i, new Pair<Integer, Double>(0, phredScaleToDouble(basesQ[i])));
					}
				}
				
			}
		}
		
		
		return cpgMethy; //key: relative location from beginning of fragment, value: 0:unmethy, 1: methy
	}
	
	static public double phredScaleToDouble(byte baseQ){
		return 1-Math.pow(10, (0-baseQ/10));
	}
	
	static public Double[] fragMethyAtEachPos(byte[] refBases, byte[] bases, byte[] basesQ, 
			boolean negStrand,  int minBaseQ, int locOffsetStart, int extend){
		Double[] methylation = new Double[extend];
		
		
		if(negStrand){
			for(int i = 0, offset = locOffsetStart; i < extend; i++, offset--){
				if(offset >= 1){
					if(SequenceUtil.basesEqual(refBases[offset], SequenceUtil.G) && SequenceUtil.basesEqual(refBases[offset-1], SequenceUtil.C)){
						if(SequenceUtil.basesEqual(bases[offset], SequenceUtil.G) && basesQ[offset] >= minBaseQ){
							methylation[i]=1.0;
						}else if(SequenceUtil.basesEqual(bases[offset], SequenceUtil.A) && basesQ[offset] >= minBaseQ){
							methylation[i]=0.0;
						}else{
							methylation[i]=Double.NaN;
						}
					}else{
						methylation[i]=Double.NaN;
					}
				}else{
					methylation[i]=Double.NaN;
				}
			}
		}else{
			for(int i = 0, offset = locOffsetStart; i < extend; i++, offset++){
				if(offset < refBases.length -1){
					if(SequenceUtil.basesEqual(refBases[offset], SequenceUtil.C) && SequenceUtil.basesEqual(refBases[offset+1], SequenceUtil.G)){
						if(SequenceUtil.basesEqual(bases[offset], SequenceUtil.C) && basesQ[offset] >= minBaseQ){
							methylation[i]=1.0;
						}else if(SequenceUtil.basesEqual(bases[offset], SequenceUtil.T) && basesQ[offset] >= minBaseQ){
							methylation[i]=0.0;
						}else{
							methylation[i]=Double.NaN;
						}
					}else{
						methylation[i]=Double.NaN;
					}
				}else{
					methylation[i]=Double.NaN;
				}
			}
		}
		
		return methylation; 
	}
	
	
	//reads CpG summary for two end together, reads should have already been filter out correctly by flag and by orientations
	static public ArrayList<byte[]> constructNonOverlapFragment(SAMRecord read1, SAMRecord read2){
		//boolean negStrandRead1 = read1.getReadNegativeStrandFlag();
		if(read1.getAlignmentStart() >= read2.getAlignmentStart()){ //first read is always at the upstream in reference genome
			SAMRecord r = read1;
			read1 = read2;
			read2 = r;
		}
		
		boolean secondPairRead1 = read1.getReadPairedFlag() && read1.getSecondOfPairFlag();
		
		
		//boolean negStrandRead2 = read2.getReadNegativeStrandFlag();
		//boolean secondPairRead2 = read2.getReadPairedFlag() && read2.getSecondOfPairFlag();
		
		
		
		byte[] refBasesRead1 = SequenceUtil.makeReferenceFromAlignment(read1, false);
		byte[] readBasesQRead1 = read1.getBaseQualities();
		byte[] readBasesRead1 = read1.getReadBases();
		
		byte[] refBasesRead2 = SequenceUtil.makeReferenceFromAlignment(read2, false);
		byte[] readBasesQRead2 = read2.getBaseQualities();
		byte[] readBasesRead2 = read2.getReadBases();
		
		ArrayList<Byte> refBases = new ArrayList<Byte>();
		ArrayList<Byte> readBases = new ArrayList<Byte>();
		ArrayList<Byte> readBaseQs = new ArrayList<Byte>();
		
		byte[] refBasesMerge;
		byte[] readBasesMerge;
		byte[] readBasesQMerge;
		
		if(!secondPairRead1){//read 1 is 1st end, if 1st and 2nd end is different, trust one with higher baseQ, if baseQ is the same, trust 1st end always.
			int read2Start=-1;
			for(int i = 0, j=1; i < refBasesRead1.length; i++, j++){
				int refPos = read1.getReferencePositionAtReadPosition(j);
				int read2Pos = read2.getReadPositionAtReferencePosition(refPos)-1;
				if(read2Pos >= 0){
					if(!BaseUtils.basesAreEqual(readBasesRead1[i], readBasesRead2[read2Pos])){
						if(readBasesQRead1[i] < readBasesRead2[read2Pos]){
							refBases.add(refBasesRead2[read2Pos]);
							readBases.add(readBasesRead2[read2Pos]);
							readBaseQs.add(readBasesQRead2[read2Pos]);
						}else{
							refBases.add(refBasesRead1[i]);
							readBases.add(readBasesRead1[i]);
							readBaseQs.add(readBasesQRead1[i]);
						}
					}else{
						refBases.add(refBasesRead1[i]);
						readBases.add(readBasesRead1[i]);
						readBaseQs.add(readBasesQRead1[i]);
					}
					read2Start = read2Pos;
				}else{
					refBases.add(refBasesRead1[i]);
					readBases.add(readBasesRead1[i]);
					readBaseQs.add(readBasesQRead1[i]);
				}
			}
			if(read2Start < 0){ //two seperate end, not overlapped
				for(int i = read1.getAlignmentEnd()-1; i < read2.getAlignmentStart(); i++){
					refBases.add(BaseUtilsMore.N);
					readBases.add(BaseUtilsMore.N);
					readBaseQs.add(BaseUtilsMore.N);
				}
			}
			
			read2Start++;
			for(int i = read2Start; i < refBasesRead2.length; i++){
				refBases.add(refBasesRead2[i]);
				readBases.add(readBasesRead2[i]);
				readBaseQs.add(readBasesQRead2[i]);
			}
			refBasesMerge =  ArrayUtils.toPrimitive(refBases.toArray(new Byte[refBases.size()]));
			readBasesMerge = ArrayUtils.toPrimitive(readBases.toArray(new Byte[readBases.size()]));
			readBasesQMerge = ArrayUtils.toPrimitive(readBaseQs.toArray(new Byte[readBases.size()]));
			//System.err.println(new String(refBasesRead1) + "\t" + new String(readBasesRead1) + "\t" + new String(readBasesQRead1)
			//				+ "\t" + new String(refBasesRead2) + "\t" + new String(readBasesRead2) + "\t" + new String(readBasesQRead2));
			//System.err.println(new String(refBasesMerge) + "\t" + new String(readBasesMerge) + "\t" + new String(readBasesQMerge));
			
		}else{
			
			int read1End=refBasesRead1.length;
			for(int i = refBasesRead2.length-1, j=i+1; i >=0 ; i--, j--){
				int refPos = read2.getReferencePositionAtReadPosition(j);
				int read1Pos = read1.getReadPositionAtReferencePosition(refPos)-1;
				if(read1Pos >= 0){
					if(!BaseUtils.basesAreEqual(readBasesRead2[i], readBasesRead1[read1Pos])){
						if(readBasesQRead2[i] < readBasesRead1[read1Pos]){
							refBases.add(refBasesRead1[read1Pos]);
							readBases.add(readBasesRead1[read1Pos]);
							readBaseQs.add(readBasesQRead1[read1Pos]);
						}else{
							refBases.add(refBasesRead2[i]);
							readBases.add(readBasesRead2[i]);
							readBaseQs.add(readBasesQRead2[i]);
						}
					}else{
						refBases.add(refBasesRead2[i]);
						readBases.add(readBasesRead2[i]);
						readBaseQs.add(readBasesQRead2[i]);
					}
					read1End = read1Pos;
				}else{
					refBases.add(refBasesRead2[i]);
					readBases.add(readBasesRead2[i]);
					readBaseQs.add(readBasesQRead2[i]);
				}
			}
			if(read1End == refBasesRead1.length){
				for(int i = read1.getAlignmentEnd()-1; i < read2.getAlignmentStart(); i++){
					refBases.add(BaseUtilsMore.N);
					readBases.add(BaseUtilsMore.N);
					readBaseQs.add(BaseUtilsMore.N);
				}
			}
			
			read1End--;
			for(int i = read1End; i >= 0; i--){
				refBases.add(refBasesRead1[i]);
				readBases.add(readBasesRead1[i]);
				readBaseQs.add(readBasesQRead1[i]);
			}
			refBasesMerge =  BaseUtilsMore.simpleReverse(ArrayUtils.toPrimitive(refBases.toArray(new Byte[refBases.size()])));
			readBasesMerge = BaseUtilsMore.simpleReverse(ArrayUtils.toPrimitive(readBases.toArray(new Byte[readBases.size()])));
			readBasesQMerge = ArrayUtils.toPrimitive(readBaseQs.toArray(new Byte[readBases.size()]));
			
			SequenceUtil.reverseQualities(readBasesQMerge);
			//System.err.println(new String(refBasesRead1) + "\t" + new String(readBasesRead1) + "\t" + new String(readBasesQRead1)
			//				+ "\t" + new String(refBasesRead2) + "\t" + new String(readBasesRead2) + "\t" + new String(readBasesQRead2));
			//System.err.println( "MERGED: " + new String(refBasesMerge) + "\t" + new String(readBasesMerge) + "\t" + new String(readBasesQMerge));
		}
		
		ArrayList<byte[]> result = new ArrayList<byte[]>();
		result.add(refBasesMerge);
		result.add(readBasesMerge);
		result.add(readBasesQMerge);
		return result;
	}
	
	
	
	//search for CH
	static public Pair<Integer, Integer> readsMethySummaryNonCG(byte[] refBases, byte[] bases, byte[] basesQ, boolean negStrand, boolean secondPair, int minBaseQ){
		
		int numC = 0;
		int numT = 0;
		if((!secondPair && !negStrand) || (secondPair && negStrand)){
			for(int i = 0; i < refBases.length -1; i ++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.C) && !SequenceUtil.basesEqual(refBases[i+1], SequenceUtil.G)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.C) && basesQ[i] >= minBaseQ){
						numC++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.T) && basesQ[i] >= minBaseQ){
						numT++;
					}
				}
			}
		}else{
			for(int i = 1; i < refBases.length; i ++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.G) && !SequenceUtil.basesEqual(refBases[i-1], SequenceUtil.C)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.G) && basesQ[i] >= minBaseQ){
						numC++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.A) && basesQ[i] >= minBaseQ){
						numT++;
					}
				}
			}
		}
		return new Pair<Integer, Integer>(numC,numT); //first is C, second is T
	}
	
	static public Pair<Integer, Integer> readsMethySummaryNonCG(SAMRecord read, int minMapQ, int minBaseQ, SamReader reader, boolean useMate){
		boolean negStrand = read.getReadNegativeStrandFlag();
		boolean secondPair = read.getReadPairedFlag() && read.getSecondOfPairFlag();
		//filter if it is bad reads
		if(read.getReadUnmappedFlag() || read.getNotPrimaryAlignmentFlag() || read.getDuplicateReadFlag() 
				|| read.getReadFailsVendorQualityCheckFlag() || (read.getReadPairedFlag() && !read.getProperPairFlag()) 
				|| read.getMappingQuality() < minMapQ){
			return new Pair<Integer, Integer>(0,0);
		}
		
		byte[] refBases = toUpperCase(modifyRefSeqByCigar(SequenceUtil.makeReferenceFromAlignment(read, false), read.getCigarString()));
		byte[] readBasesQ = getClippedReadsBaseQuality(read);
		byte[] readBases = toUpperCase(getClippedReadsBase(read));
		if(negStrand){
			readBases = complementArray(readBases);
		}
		
		Pair<Integer, Integer> ctSummary = readsMethySummaryNonCG(refBases, readBases, readBasesQ, negStrand, secondPair, minBaseQ);
		
		
		
		if(read.getReadPairedFlag() && useMate){
			SAMRecord m = reader.queryMate(read);
			byte[] refBasesMate = toUpperCase(modifyRefSeqByCigar(SequenceUtil.makeReferenceFromAlignment(m, false), m.getCigarString()));
			byte[] basesMate = toUpperCase(getClippedReadsBase(m));
			byte[] basesQMate = getClippedReadsBaseQuality(read);
			boolean negStrandMate = m.getReadNegativeStrandFlag();
			if(negStrandMate){
				basesMate = complementArray(basesMate);
			}
			boolean secondPairMate = m.getReadPairedFlag() && m.getSecondOfPairFlag();
			Pair<Integer, Integer> ctSummaryMate = readsMethySummaryNonCG(refBasesMate, basesMate, basesQMate, negStrandMate, secondPairMate, minBaseQ);
			ctSummary = new Pair<Integer, Integer>(ctSummary.getFirst()+ctSummaryMate.getFirst(), ctSummary.getSecond()+ctSummaryMate.getSecond());
		}
		
		return ctSummary; //first is C, second is T
	}

///also output CG number and CG locations Pair<Pair<Integer, Integer>, HashMap<Integer,Integer>>
	
	


	//return -1 means the whole reads is filtered out, otherwise, the >=0 result indicates the bisulfite conversion start position. 0 means the whole reads is passed, and there is no bisulfite conversion
	//start position
	static public int bisulfiteIncompleteReads(SAMRecord read) throws Exception{
		return bisulfiteIncompleteReads(read.getReadNegativeStrandFlag() ? complementArray(toUpperCase(getClippedReadsBase(read))) : toUpperCase(getClippedReadsBase(read)), toUpperCase(modifyRefSeqByCigar(refStrFromMd(read), read.getCigarString()))
				, read.getReadNegativeStrandFlag(), read.getReadPairedFlag() & read.getSecondOfPairFlag(), false, true, true, 0.1, 0.1, 0, "CH", 1);
	}

	//return -1 means the whole reads is filtered out, otherwise, the >=0 result indicates the bisulfite conversion start position. 0 means the whole reads is passed, and there is no bisulfite conversion
		//start position
	static public int bisulfiteIncompleteReads(byte[] bases, byte[] refBases, boolean negativeStrand, boolean secondEnd) throws Exception{
		return bisulfiteIncompleteReads(bases, refBases, negativeStrand, secondEnd, true, true, true, 0.1, 0.4, 0, "CH", 1);
	}
	
	//return -1 means the whole reads is filtered out, otherwise, the >=0 result indicates the bisulfite conversion start position. 0 means the whole reads is passed, and there is no bisulfite conversion
	//start position
	static public int bisulfiteIncompleteReads(byte[] bases, byte[] refBases, boolean negativeStrand, boolean secondEnd, boolean mismatchFilter, 
			boolean bisulfiteIncompletReads, boolean bisulfiteFivePrimeConvFilter, 
			double maxMismatches, double minPatConvRate, int minConv, String patConv5, int posCinPatConv5) throws Exception{
				
		int readLength = bases.length;
		boolean originalNegStrand = negativeStrand;
		if(secondEnd)
			negativeStrand = !negativeStrand;


		if(bisulfiteIncompletReads){
			if(negativeStrand){
				bases = BaseUtils.simpleReverseComplement(bases);
				refBases = BaseUtils.simpleReverseComplement(refBases);
				
			}
		}
		byte[] patterns = patConv5.getBytes();
		
		
		int numberOfMismatches = 0;
		int numberOfPattern = 0;
		int numberOfPatternInRef = 0; 
		int convertedCount = 0;
		
		for(int i = 0; i < readLength; i++){
			if(mismatchFilter){
				if( !BaseUtils.basesAreEqual(refBases[i], bases[i]) && BaseUtilsMore.isBisulfiteMismatch(refBases[i], bases[i],false, false)) //already get reverse complement for negative strand, so no need to provide flag again...
					numberOfMismatches++;
				if(numberOfMismatches > maxMismatches * readLength){
					//System.err.println(new String(refBases) + "\t" + new String(bases) + "\t" + numberOfMismatches + "\t" + negativeStrand + "\t" + secondEnd);
					return -1;
				}
					
			}
			
			if((bisulfiteIncompletReads || bisulfiteFivePrimeConvFilter) && i <= readLength-patterns.length){
				short numMatchesInRef=0;
				short numMatchesInReads=0;
				boolean conv = false;
				
				for(int j = i, index = 0; index < patterns.length; index++, j++){
					
					if(numberOfPatternInRef != -1){ //mean reads is mapped
						//System.err.println(index + "\t" + j + "\t" + readLength + "\t" + patterns.length + "\t" + i + "\t" + refBases.length + '\t' + originalNegStrand);
						if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(patterns[index], refBases[j])){
							numMatchesInRef++;
							
						}
						
					}
					if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(patterns[index], bases[j])){
						numMatchesInReads++;
					}
					if(index == (posCinPatConv5-1) && BaseUtils.basesAreEqual(bases[j], BaseUtilsMore.T)){
						conv=true;
					}
					
				}
				if(numMatchesInReads == patterns.length){
					numberOfPattern++;
				}
				if(numberOfPatternInRef != -1 && numMatchesInRef == patterns.length){
					numberOfPatternInRef++;
				}
				
				if(bisulfiteFivePrimeConvFilter && numMatchesInReads == patterns.length-1){
					if( numMatchesInRef == (patterns.length)){
						if(conv)
							convertedCount++;
						
						if(convertedCount >= minConv && minConv>0){
							if(negativeStrand){
								return readLength-i;
							}else{
								return i;
							}
							
						}
					}
				}
			}
			
			
			
			
		}
		//if(new String(refBases).equalsIgnoreCase("tatctgcgttacactcgtcg")){
		//	System.err.println("res: " + numberOfPattern + "\t" + numberOfPatternInRef);
		//}
		if(bisulfiteIncompletReads){
			//System.err.println("res: " + numberOfPattern + "\t" + numberOfPatternInRef);
			if(((numberOfPatternInRef == -1 || numberOfPatternInRef == 0) && numberOfPattern == 0) || (numberOfPatternInRef != -1 && (double)numberOfPattern/(double)numberOfPatternInRef < minPatConvRate)){

			}else{
				
				return -1;
			}
		}

		return 0;
	}
	
	public static int patternCountSearch(byte[] strings, String pattern){
		//byte[] strings = string.getBytes();
		byte[] patterns = pattern.getBytes();
		int count = 0;
		for(int i = 0; i < strings.length - patterns.length + 1; i++){
			for(int j=i, index=0; index < patterns.length; j++,index++){
				if(BaseUtils.basesAreEqual(strings[j], patterns[index])){
					count++;
				}
			}
		}
		return count;
	}

	public static double patternFreqSearch(byte[] strings, String pattern){
		//byte[] strings = string.getBytes();
		byte[] patterns = pattern.getBytes();
		int count = 0;
		int totalCount=0;
		for(int i = 0; i < strings.length - patterns.length + 1; i++, totalCount++){
			int matches = 0;
			for(int j=i, index = 0; index < patterns.length; j++,index++){
				if(BaseUtils.basesAreEqual(strings[j], patterns[index])){
					matches++;
				}
			}
			if(matches == patterns.length){
				count++;
			}
		}
		//System.err.println(pattern + "\t" + new String(strings) + "\t" + count + "\t" + totalCount);
		return (double)count/(double)totalCount;
	}	
	
	
	public static LinkedHashMap<String, Double> kmerFreqSearch(byte[] strings, int kmer){
		//byte[] strings = string.getBytes();
		HashMap<String, Integer> kmerCount = new HashMap<String, Integer>();
		

		int totalCount=0;
		for(int i = 0; i < strings.length - kmer + 1; i++, totalCount++){
			String bases = "";
			for(int j=i, index = 0; index < kmer; j++,index++){
				bases = bases + Character.toUpperCase((char)strings[j]);
			}
			if(kmerCount.containsKey(bases)){
				kmerCount.put(bases, kmerCount.get(bases)+1);
			}else{
				kmerCount.put(bases, 1);
			}
		}
		//System.err.println(new String(strings) + "\t" + kmer + "\t" + kmerCount.size() + "\t" + totalCount);
		//for(String key : kmerCount.keySet()){
		//	System.err.print(key + "\t" + kmerCount.get(key) + "\t");
		//}
		//System.err.println();
		LinkedHashMap<String, Double> kmerFreq = new LinkedHashMap<String, Double>();
		for(byte[] k : SequenceUtil.generateAllKmers(kmer)){
			String p = new String(k);
			if(kmerCount.containsKey(p)){
				kmerFreq.put(p, 100*(double)kmerCount.get(p)/(double)totalCount);
			//	System.err.print("\t" + kmerCount.get(p));
			}else{
				kmerFreq.put(p, 0.0);
			}
		}
		//System.err.println();
		//System.err.println(pattern + "\t" + new String(strings) + "\t" + count + "\t" + totalCount);
		return kmerFreq;
	}	
	
	
	//From ben's picardUtils, need to test it in paired end space
	public static byte[] refStrFromMd(String seq, String md, Cigar cigar, String readName)
			throws Exception
			{
				if (seq == null) throw new Exception("Can not run refStrFromMd with a null seq variable");
				if (md == null) throw new Exception("Can not run refStrFromMd with a null MD variable. Reads is: " + readName);
				
				//Pattern mdPat = Pattern.compile("\\G(?:([0-9]+)|([ACTGNactgn])|(\\^[ACTGNactgn]+))");
				Pattern mdPat = Pattern.compile("\\G(?:([0-9]+)|([ABCDGHKMNRSTVWXYabcdghkmnrstvwxy])|(\\^[ABCDGHKMNRSTVWXYabcdghkmnrstvwxy]+))");
				
				// Use sb as the reference output string
				StringBuilder sb = new StringBuilder(500*2+1);

				Matcher match = mdPat.matcher(md);
				int curSeqPos = 0;
				//int curMdPos = 0; // Not the same as seq pos when you have indels

				int savedBases = 0;
				for (final CigarElement cigEl : cigar.getCigarElements()) 
				{
					int cigElLen = cigEl.getLength();
					CigarOperator cigElOp = cigEl.getOperator();
//					System.err.printf("\tCigar El: len=%d, op=%s, consumesRead=%b, consumesRef=%b\n",
//							cigElLen,cigElOp,cigElOp.consumesReadBases(), cigElOp.consumesReferenceBases());
					
					
					// If it consumes reference bases, it's either a match or a deletion in the sequence
					// read.  Either way, we're going to need to parse throught the MD.
					if (cigElOp.consumesReferenceBases())
					{
						// We have a match region, go through the MD
						int basesMatched = 0;
						
						// Do we have any saved matched bases?
						while ((savedBases>0) && (basesMatched < cigElLen))
						{
							sb.append(seq.charAt(curSeqPos++));
							savedBases--;
							basesMatched++;
//							System.err.printf("\t\tDepleting saved bases, saved=%d, curSeqPos=%d, basesMatched=%d\n",savedBases,curSeqPos,basesMatched); 
						}

						while (basesMatched < cigElLen)
						{
							boolean matched = match.find();
							if (matched)
							{
//								System.err.println("Matched , basesMatched=" + basesMatched + ", match=" + match.group() + "," + match.group(1) + "," + match.group(2) + ", start=" + match.start());
								String mg;
								if ( ((mg = match.group(1)) !=null) && (mg.length() > 0) )
								{
									// It's a number , meaning a series of matches
									int num = Integer.parseInt(mg);
									for (int i = 0; i < num; i++)
									{
										if (basesMatched<cigElLen)
										{
											sb.append(seq.charAt(curSeqPos));
											curSeqPos++;
										}
										else
										{
											savedBases++;
										}
										basesMatched++;
									}
								}

								else if ( ((mg = match.group(2)) !=null) && (mg.length() > 0) )
								{
									// It's a single nucleotide, meaning a mismatch
									if (basesMatched<cigElLen)
									{
										sb.append(mg.charAt(0));
										curSeqPos++;
									}
									else
									{
										savedBases++;
									}
									basesMatched++;
								}
								else if ( ((mg = match.group(3)) !=null) && (mg.length() > 0) )
								{
									// It's a deletion, starting with a caret
									// don't include caret
									for (int i = 1; i < mg.length(); i++)
									{
										// Since this function is actually just meant to make a reference that lines up nucleotide 
										//  for nucleotide with the sequence read, we don't actually add the insertion to the reference.
										//sb.append(mg.charAt(i));
										basesMatched++;
									}
									
									// Check just to make sure.
									if (basesMatched != cigElLen)
									{
										throw new Exception("Got a deletion in CIGAR (" + cigar + ", deletion " + cigElLen + 
												" length) with an unequal ref insertion in MD (" + md + ", md " + basesMatched + " length" + "in Read: " + readName);
									}
									if (cigElOp != CigarOperator.DELETION)
									{
										throw new Exception ("Got an insertion in MD ("+md+") without a corresponding deletion in cigar ("+cigar+")" + "in Read: " + readName);
									}
									
								}
								else
								{
									matched = false;
									
								}
							}

							if (!matched)
							{
								throw new Exception("Illegal MD pattern: " + md + "in Read: " + readName);
							}

//							System.err.println("SavedBasesMatched=" + savedBases);
						}

					}
					else if (cigElOp.consumesReadBases())
					{
						// We have an insertion in read
						for (int i = 0; i < cigElLen; i++)
						{
							char c = (cigElOp == CigarOperator.SOFT_CLIP) ? '0' : '-';
							sb.append( c );
							curSeqPos++;
						}
					}
					else
					{
						// It's an op that consumes neither read nor reference bases.  Do we just ignore??
					}

				}
				
				return sb.toString().getBytes();
			}
	

	public static byte[] refStrFromMd(SAMRecord read)
			throws Exception
			{
				return refStrFromMd(read.getReadString(), read.getStringAttribute("MD"),read.getCigar(), read.getReadName());
				
		
			}

	
	/**
	 * 
	 * @param contig
	 *            give contig string
	 * @return contigIndex return a number represent this contig
	 */
	static public int contigNameToId(String contig) {
		Pattern replace = Pattern.compile("^chr");
        Matcher matcher1 = replace.matcher(contig);
        contig=matcher1.replaceAll("");
		int contigs = 0;
		if(contig.equalsIgnoreCase("MT")){
			contigs = 25;
		}else if(contig.equalsIgnoreCase("X")){
			contigs = 23;
		}else if(contig.equalsIgnoreCase("Y")){
			contigs = 24;
		}else{
			contigs = Integer.parseInt(contig);
		}
		return contigs;
	}


	public static boolean passReadPairOrientation(SAMRecord first, SAMRecord second){
		boolean firstNegStrand = first.getReadNegativeStrandFlag();
		boolean secondNegStrand = second.getReadNegativeStrandFlag();
		boolean firstSecondEnd = first.getSecondOfPairFlag();
		boolean secondSecondEnd = second.getSecondOfPairFlag();
		
		
		int firstStartPos = first.getAlignmentStart();
		int firstEndPos = first.getAlignmentEnd();
		int secondStartPos = second.getAlignmentStart();
		int secondEndPos = second.getAlignmentEnd();
		boolean corFlag = (firstEndPos >= firstStartPos && firstEndPos <= secondEndPos) && (secondStartPos >= firstStartPos && secondStartPos <= secondEndPos);
		//need to test if alignment start always < alignmentEnd, yes, that is true;
		if(!firstSecondEnd && !firstNegStrand && secondSecondEnd && secondNegStrand){
			return corFlag;
		}else if(firstSecondEnd && firstNegStrand && !secondSecondEnd && !secondNegStrand){
			return corFlag;
		}else if(!firstSecondEnd && firstNegStrand && secondSecondEnd && !secondNegStrand){
			return corFlag;
		}else if(firstSecondEnd && !firstNegStrand && !secondSecondEnd && secondNegStrand){
			return corFlag;
		}else{
			return false;
		}
		
	}
	
	//when can't pass 2nd end reads, just use first end reads length and strand to calculate the possible mate reads end cor. 
	public static boolean passReadPairOrientation(SAMRecord first){
		boolean firstNegStrand = first.getReadNegativeStrandFlag();
		boolean secondNegStrand = first.getMateNegativeStrandFlag();
		boolean firstSecondEnd = first.getSecondOfPairFlag();
		boolean secondSecondEnd = !firstSecondEnd;
		
		
		int firstStartPos = first.getAlignmentStart();
		int firstEndPos = first.getAlignmentEnd();
		int secondStartPos = first.getMateAlignmentStart();
		int secondEndPos = secondStartPos + first.getReadLength();
		if(firstStartPos >= secondStartPos){
			firstStartPos = first.getMateAlignmentStart();
			firstEndPos = firstStartPos + first.getReadLength();
			secondStartPos = first.getAlignmentStart();
			secondEndPos = first.getAlignmentEnd();
			
			firstNegStrand = first.getMateNegativeStrandFlag();
			secondNegStrand = first.getReadNegativeStrandFlag();
			secondSecondEnd = first.getSecondOfPairFlag();
			firstSecondEnd = !secondSecondEnd;
			
		}
		
		boolean corFlag = (firstEndPos >= firstStartPos && firstEndPos <= secondEndPos) && (secondStartPos >= firstStartPos && secondStartPos <= secondEndPos);
		//need to test if alignment start always < alignmentEnd, yes, that is true;
		if(!firstSecondEnd && !firstNegStrand && secondSecondEnd && secondNegStrand){
			return corFlag;
		//}else if(firstSecondEnd && firstNegStrand && !secondSecondEnd && !secondNegStrand){
		//	return corFlag;
		//}else if(!firstSecondEnd && firstNegStrand && secondSecondEnd && !secondNegStrand){
		//	return corFlag;
		}else if(firstSecondEnd && !firstNegStrand && !secondSecondEnd && secondNegStrand){
			return corFlag;
		}else{
			return false;
		}
		
	}
	
	public static int getFragOffsetFromReadsOffset(SAMRecord r, int readOffset){
		int rStartPos = r.getAlignmentStart();
		int mateStartPos = r.getMateAlignmentStart();

		int fragLen = Math.abs(r.getInferredInsertSize());
		boolean leftEndInFragment = false;
		if(rStartPos < mateStartPos){
			leftEndInFragment = true;
		}

		if(leftEndInFragment){
			return readOffset;
		}else{
			//System.err.println(leftEndInFragment + "\t" + (fragLen-(r.getReadLength()-readOffset)));
			return fragLen-(r.getReadLength()-readOffset);
			
		}
		
	}
	
	public static int getDistFragEndFromReadsOffset(SAMRecord r, int readOffset){
		int rStartPos = r.getAlignmentStart();
		int mateStartPos = r.getMateAlignmentStart();

		boolean leftEndInFragment = false;
		if(rStartPos < mateStartPos){
			leftEndInFragment = true;
		}

		if(leftEndInFragment){
			return readOffset;
		}else{
			//System.err.println(leftEndInFragment + "\t" + (r.getReadLength()-readOffset));
			return r.getReadLength()-readOffset;
		}
		
	}
	
	public static boolean trimmedBase(int offset, int readLen, boolean trim2nd, int trim5, int trim3, boolean negativeStrand, boolean secondEnd){

		if(trim2nd){
			if(secondEnd){
				if(negativeStrand){
					return offset >= trim3 && offset <= readLen - trim5;
					
				}
				else{
					return offset >= trim5 && offset <= readLen - trim3;
						
				}
			}
			else{
				return true;
			}
		}
		else{
			if(negativeStrand){
				return offset >= trim3 && offset <= readLen - trim5;
					
			}
			else{
				return offset >= trim5 && offset <= readLen - trim3;

			}
		}
		
	}
	
	
	//return from interval a's boundary to interval b's boundary, if a has strand information, based on a's strand orientation for the direction
	public static int intervalDistance(Interval a, Interval b){
		if(a == null || b == null){
			throw new RuntimeException("Interval is null");
		}
		if(!a.getContig().equalsIgnoreCase(b.getContig())){
			return Integer.MAX_VALUE;
		}else{
			if(a.abuts(b)){
				return 0;
			}else{
				int[] distances = new int[]{
					b.getStart() - a.getStart(),
					b.getStart() - a.getEnd(),
					b.getEnd() - a.getStart(),
					b.getEnd() - a.getEnd()
				};
				int minAbs = Math.abs(distances[0]);
				int minDis = distances[0];
				for(int i = 1; i < distances.length; i++){
					if(Math.abs(distances[i]) < minAbs){
						minAbs = Math.abs(distances[i]);
						minDis = distances[i];
					}
				}
				
				if(a.isNegativeStrand()){
					return 0-minDis;
				}else{
					return minDis;
				}
				
			}
		}
		
	}
		
	//if no "+" or "-" strand information, it will return abolute distance information
	public static int intervalDistance(IntervalTree.Node<String> a, Interval b){
		
		
			//overlap
		if((b.getStart() >= a.getStart() && b.getStart() <= a.getEnd())||(b.getEnd() >= a.getStart() && b.getEnd() <= a.getEnd())){
			return 0;
		}
		
				int[] distances = new int[]{
					b.getStart() - a.getStart(),
					b.getStart() - a.getEnd(),
					b.getEnd() - a.getStart(),
					b.getEnd() - a.getEnd()
				};
				int minAbs = Math.abs(distances[0]);
				int minDis = distances[0];
				for(int i = 1; i < distances.length; i++){
					if(Math.abs(distances[i]) < minAbs){
						minAbs = Math.abs(distances[i]);
						minDis = distances[i];
					}
				}
				
				if(a.getValue().equalsIgnoreCase("-")){
					return 0-minDis;
				}else if(a.getValue().equalsIgnoreCase("+")){
					return minDis;
				}else{
					return minAbs;
				}
				
			
		
	}
	
	public static int intervalDistance(IntervalTree.Node<String> a, IntervalTree.Node<String> b){
		
		
		//overlap
	if((b.getStart() >= a.getStart() && b.getStart() <= a.getEnd())||(b.getEnd() >= a.getStart() && b.getEnd() <= a.getEnd())){
		return 0;
	}
	
			int[] distances = new int[]{
				b.getStart() - a.getStart(),
				b.getStart() - a.getEnd(),
				b.getEnd() - a.getStart(),
				b.getEnd() - a.getEnd()
			};
			int minAbs = Math.abs(distances[0]);
			int minDis = distances[0];
			for(int i = 1; i < distances.length; i++){
				if(Math.abs(distances[i]) < minAbs){
					minAbs = Math.abs(distances[i]);
					minDis = distances[i];
				}
			}
			
			if(a.getValue().equalsIgnoreCase("-")){
				return 0-minDis;
			}else if(a.getValue().equalsIgnoreCase("+")){
				return minDis;
			}else{
				return minAbs;
			}
			
		
	
}
	
	//for HMM utils 
	public static <T> List<T> flat(List<Pair<HashMap<Integer, Integer>, List<T>>> lists)
	{	
		List<T> v = new ArrayList<T>();
		
		for (Pair<HashMap<Integer, Integer>, List<T>> list : lists)
			v.addAll(list.getSecond());
		
		return v;
	}
	
	public static <T> List<T> flatPair(List<Pair<HashMap<Integer, Pair<Integer, Double>>, List<T>>> lists)
	{	
		List<T> v = new ArrayList<T>();
		
		for (Pair<HashMap<Integer, Pair<Integer, Double>>, List<T>> list : lists)
			v.addAll(list.getSecond());
		
		return v;
	}
	
	public static <T> List<T> flatTripleString(List<Pair<HashMap<Integer, Triple<Integer, Double, String>>, List<T>>> lists)
	{	
		List<T> v = new ArrayList<T>();
		
		for (Pair<HashMap<Integer, Triple<Integer, Double, String>>, List<T>> list : lists)
			v.addAll(list.getSecond());
		
		return v;
	}
	
	public static <T> List<T> flatTriple(List<Pair<HashMap<Integer, Triple<Integer, Double, Integer>>, List<T>>> lists)
	{	
		List<T> v = new ArrayList<T>();
		
		for (Pair<HashMap<Integer, Triple<Integer, Double, Integer>>, List<T>> list : lists)
			v.addAll(list.getSecond());
		
		return v;
	}
	
	public static <T> List<T> flatTriple(List<Pair<HashMap<Integer, Triple<Integer, Double, Integer>>, List<T>>> lists, int nbGenicStates)
	{	
		List<T> v = new ArrayList<T>();
		
		for (Pair<HashMap<Integer, Triple<Integer, Double, Integer>>, List<T>> list : lists){
			List<T> listContent = list.getSecond();
			HashMap<Integer, Triple<Integer, Double, Integer>> listState = list.getFirst();
			for(int i = 0; i < listContent.size(); i++){
				if(listState.get(i).getRight() == nbGenicStates){
					v.add(listContent.get(i));
				}
			}
			
		}
			
		
		return v;
	}
	
	static public byte[] modifyRefSeqByCigar(byte[] seqByte, String cigarString){ //mainly for "I" and "D" CIGAR string
		char[] cigarList = CigarUtil.cigarArrayFromString(cigarString);

		ArrayList<Byte> seqsNew = new ArrayList<Byte>();
		//System.err.println(new String(seqByte) + "\t" + seqByte.length + "\t" + cigarString);
		int offSet = 0;
		for(char cigar : cigarList){
			if(cigar == 'M' || cigar == 'I'){
				seqsNew.add(seqByte[offSet]);
				offSet++;
			}else if(cigar == 'D'){
				//seqsNew.add(SequenceUtil.N);
				
			}else if(cigar == 'S'){
				offSet++;
			}
		}
		return ArrayUtils.toPrimitive(seqsNew.toArray(new Byte[seqsNew.size()]));
	}

	static public byte[] modifyRefBasesByCigar(byte[] seqByte, String cigarString){ //mainly for "I" and "D" CIGAR string
		char[] cigarList = CigarUtil.cigarArrayFromString(cigarString);

		ArrayList<Byte> seqsNew = new ArrayList<Byte>();
		//System.err.println(new String(seqByte) + "\t" + seqByte.length + "\t" + cigarString);
		int offSet = 0;
		for(char cigar : cigarList){
			if(cigar == 'M'){
				seqsNew.add(seqByte[offSet]);
				offSet++;
			}else if(cigar == 'I'){
				if(offSet>0){
					seqsNew.add((byte)('-'));
				}
				
				
				
			}else if(cigar == 'D'){
				//seqsNew.add(SequenceUtil.N);
				offSet++;
				
			}else if(cigar == 'S'){
				
			}
		}
		return ArrayUtils.toPrimitive(seqsNew.toArray(new Byte[seqsNew.size()]));
	}
	
	static public int adjPosByCigar(int pos, SAMRecord r){ //mainly for "I" and "D" CIGAR string
		char[] cigarList = CigarUtil.cigarArrayFromString(r.getCigarString());
		int start = r.getAlignmentStart();
		int offset = 0;
		int posAdj = 0;
		for(char cigar : cigarList){
			if(start+offset >= pos){
				break;
			}
			
			if(cigar == 'M'){
				offset++;
			}else if(cigar == 'I'){
				if(offset>0){
					posAdj--;
				}
				
			}else if(cigar == 'D'){
				posAdj++;
				offset++;
				
			}else if(cigar == 'S'){
				
			}
		}
		return pos+posAdj;
	}

	
	static public int getMateAlignmentEndByMateCigar(SAMRecord r){ //mainly for "I" and "D" CIGAR string
		
		if(r.getStringAttribute("MC") == null){
			return r.getReadLength() + r.getMateAlignmentStart() - 1;
		}
		char[] cigarList = CigarUtil.cigarArrayFromString(r.getStringAttribute("MC"));

		int len = 0;
		for(char cigar : cigarList){
			if(cigar == 'M' || cigar == 'I'){
				
				len++;
			}else if(cigar == 'D'){
				//seqsNew.add(SequenceUtil.N);
				
			}else if(cigar == 'S'){
				
			}
		}
		return len + r.getMateAlignmentStart() - 1;
	}
	
	static public byte[] getClippedReadsBaseQuality(SAMRecord r){
		int truncateStart = r.getReadPositionAtReferencePosition(r.getAlignmentStart())-1;
		int truncateEnd = r.getReadPositionAtReferencePosition(r.getAlignmentEnd());
		byte[] truncatedBaseQs = new byte[truncateEnd-truncateStart];
		byte[] baseQs = r.getBaseQualities();
		if(truncateStart<0){
			truncateStart=0;
		}
		if(truncateEnd>baseQs.length){
			truncateEnd = baseQs.length;
		}
		//if(r.getReadNegativeStrandFlag()){//TODO: remove the whole if, after using the correct bam file
		//	ArrayUtils.reverse(baseQs);
		//}
		for(int i = truncateStart, index=0; i < truncateEnd; i++, index++){
			truncatedBaseQs[index] = baseQs[i];
		}
		
		return truncatedBaseQs;
	}
	
	static public byte[] getClippedReadsBase(SAMRecord r){
		int truncateStart = r.getReadPositionAtReferencePosition(r.getAlignmentStart())-1;
		int truncateEnd = r.getReadPositionAtReferencePosition(r.getAlignmentEnd());
		byte[] truncatedBases = new byte[truncateEnd-truncateStart];
		byte[] bases = r.getReadBases();
		//if(r.getReadNegativeStrandFlag()){ //TODO: remove the whole if, after using the correct bam file
		//	ArrayUtils.reverse(bases);
		//}
		if(truncateStart<0){
			truncateStart=0;
		}
		if(truncateEnd>bases.length){
			truncateEnd = bases.length;
		}
		for(int i = truncateStart, index=0; i < truncateEnd; i++, index++){
			truncatedBases[index] = bases[i];
		}
		
		return truncatedBases;
	}
	
	static public byte[] complementArray(byte[] a){
		byte[] b = new byte[a.length];
		for(int i = 0; i < a.length; i++){
			b[i]=SequenceUtil.complement(a[i]);
		}
		return b;
	}

}
