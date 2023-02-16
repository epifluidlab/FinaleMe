/**
 * 
 */
package org.cchmc.epifluidlab.methhmm.utils;

import htsjdk.samtools.SAMRecord;

import org.broadinstitute.gatk.engine.filters.ReadFilter;





/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Dec 11, 2012 11:22:51 AM
 * 
 */
public class InvertedDupsReadFilter extends ReadFilter {

	/* (non-Javadoc)
	 * @see net.sf.picard.filter.SamRecordFilter#filterOut(net.sf.samtools.SAMRecord)
	 */
	@Override
	public boolean filterOut(SAMRecord samRecord) {
		if (samRecord.getReadPairedFlag() && (samRecord.getAlignmentStart() == samRecord.getMateAlignmentStart() && samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag()))
		{
			if (samRecord.getSecondOfPairFlag()) return true;
		}
		return false;
	}

	/* (non-Javadoc)
	 * @see net.sf.picard.filter.SamRecordFilter#filterOut(net.sf.samtools.SAMRecord, net.sf.samtools.SAMRecord)
	 */
	@Override
	public boolean filterOut(SAMRecord arg0, SAMRecord arg1) {
		// TODO Auto-generated method stub
		return false;
	}

}
