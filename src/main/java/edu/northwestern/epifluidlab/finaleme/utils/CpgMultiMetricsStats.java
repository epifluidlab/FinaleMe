package edu.northwestern.epifluidlab.finaleme.utils;

/**
 * @deprecated Use {@link CpgFeatureMatrixBuilder} instead.
 *             This alias is kept for backward-compatible script entrypoints.
 */
@Deprecated
public class CpgMultiMetricsStats extends CpgFeatureMatrixBuilder {

	public static void main(String[] args) throws Exception {
		CpgMultiMetricsStats alias = new CpgMultiMetricsStats();
		alias.doMain(args);
	}
}

