/**
 * Shared utilities for CpgMultiMetricsStats variants.
 */
package edu.northwestern.epifluidlab.finaleme.utils;

import htsjdk.samtools.util.IntervalTree;
import org.slf4j.Logger;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

abstract class AbstractCpgMultiMetricsStats {

	protected BufferedReader openMaybeGzipReader(String path) throws IOException {
		if(path.endsWith(".gz")){
			return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path))));
		}
		return new BufferedReader(new FileReader(path));
	}

	protected HashMap<String, IntervalTree<Integer>> loadExcludeRegions(ArrayList<String> excludeRegions, Logger log)
			throws IOException {
		if(excludeRegions == null || excludeRegions.isEmpty()){
			return null;
		}

		log.info("Excluding intervals ... ");
		HashMap<String,IntervalTree<Integer>> ignoreLocCollections = new HashMap<String,IntervalTree<Integer>>();
		for(String excludeRegion : excludeRegions){
			BufferedReader br = openMaybeGzipReader(excludeRegion);
			String line;
			while((line = br.readLine()) != null){
				if(line.startsWith("#")){
					continue;
				}
				String[] splitLines = line.split("\t");
				if(splitLines.length < 3){
					continue;
				}
				String chr = splitLines[0];
				int start = Integer.parseInt(splitLines[1]);
				int end = Integer.parseInt(splitLines[2]);
				IntervalTree<Integer> tree = ignoreLocCollections.containsKey(chr)
						? ignoreLocCollections.get(chr)
						: new IntervalTree<Integer>();
				tree.put(start, end, 1);
				ignoreLocCollections.put(chr, tree);
			}
			br.close();
		}
		return ignoreLocCollections;
	}

	protected HashMap<String, IntervalTree<String>> loadStrandedIntervals(String inputFile, Logger log, String logLabel)
			throws IOException {
		log.info(logLabel);
		HashMap<String,IntervalTree<String>> collections = new HashMap<String,IntervalTree<String>>();
		BufferedReader br = openMaybeGzipReader(inputFile);
		String line;
		while((line = br.readLine()) != null){
			if(line.startsWith("#")){
				continue;
			}
			String[] splitLines = line.split("\t");
			if(splitLines.length < 3){
				continue;
			}
			String chr = splitLines[0];
			int start = Integer.parseInt(splitLines[1]);
			int end = Integer.parseInt(splitLines[2]);

			IntervalTree<String> tree = collections.containsKey(chr)
					? collections.get(chr)
					: new IntervalTree<String>();

			String strand = ".";
			if(splitLines.length >= 6){
				if(splitLines[5].equalsIgnoreCase("-")){
					strand = "-";
				}else if(splitLines[5].equalsIgnoreCase("+")){
					strand = "+";
				}
			}
			tree.put(start, end, strand);
			collections.put(chr, tree);
		}
		br.close();
		return collections;
	}
}
