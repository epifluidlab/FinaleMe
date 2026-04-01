# FinaleMe Optimization Plan

## Context

FinaleMe is a Java bioinformatics tool (~9,600 lines) for predicting DNA methylation from cell-free DNA WGS data using a Non-Homogenous Bayesian HMM. The entire pipeline was **single-threaded** with several I/O and algorithmic inefficiencies. The tool processes BAM files (often 30x WGS, tens of GB), generates per-CpG feature matrices, trains an HMM via Baum-Welch, and decodes via Viterbi. Typical runs require `-Xmx100G` and can take many hours. The goal is to dramatically reduce runtime while preserving numerical correctness.

---

## Implementation Status

| Item | Status |
|------|--------|
| 1A. Parallelize Viterbi Decoding | Done |
| 1B. Parallelize Baum-Welch Forward-Backward | Done |
| 1C. Parallelize Feature Extraction by 5Mb Genomic Bins | Done |
| 1D. Parallelize KL-Divergence Computation | Done |
| 2A. Single-Pass Matrix File Processing | Done |
| 2B. Eliminate Full BAM Scan for Read Counting | Done |
| 2C. Cache flatPair in Baum-Welch | Done |
| 3A. Replace Boxed Types with Primitives in HMM Model | Done |
| 3B. Reduce allGamma Memory Footprint | Done |
| 3C. Streaming Read Processing in Feature Extraction | Not yet (high risk) |
| 4A. Consolidate CpgMultiMetricsStats Variants | Not yet |
| 4B. Dependency and Java Version Upgrades | Done |
| 4C. Add Unit Tests for HMM Core | Done |
| 4D. Remove Dead Code | Done |

---

## TIER 1: Parallelization (Highest Impact)

### 1A. Parallelize Viterbi Decoding (Done)
- **File**: [FinaleMe.java](src/main/java/org/cchmc/epifluidlab/finaleme/hmm/FinaleMe.java) `decodeHmm()` (~line 997)
- **Previous**: Sequential `for(int j=0; j < matrix.size(); j++)` loop, each fragment independently decoded
- **Change**: Use `ExecutorService` with `ForkJoinPool` to decode fragments in parallel. The HMM model is read-only during decode. Collect results into a pre-sized array, then aggregate `methySummary` sequentially afterward.
- **Expected speedup**: ~Nx where N = core count (8-16x typical)

### 1B. Parallelize Baum-Welch Forward-Backward (Done)
- **File**: [BaumWelchBayesianNhmmV5ScaledLearner.java](src/main/java/org/cchmc/epifluidlab/finaleme/hmm/BaumWelchBayesianNhmmV5ScaledLearner.java) `iterate()`
- **Previous**: Sequential loop over all sequences computing forward-backward, xi, gamma, then accumulating into shared `aijNum/aijDen/arij`
- **Change**: Map-reduce pattern -- parallel forward-backward/xi/gamma computation across sequences, then sequential accumulation of results into transition matrices
- **Expected speedup**: ~Nx per Baum-Welch iteration (up to 50 iterations)

### 1C. Parallelize Feature Extraction by 5Mb Genomic Bins (Done)
- **File**: [CpgMultiMetricsStats.java](src/main/java/org/cchmc/epifluidlab/finaleme/utils/CpgMultiMetricsStats.java)
- **Previous**: Sequential `for(String chr : cpgCollections.keySet())` with per-CpG BAM index queries
- **Change**: Divide genome into 5Mb non-overlapping bins. Each bin processed in parallel with its own `SamReader`, `TwoBitParser`, `TabixFeatureReader`, and `BigWigFileReader` instances. Each bin writes to a temp file, then results are concatenated. Thread pool sized to `Runtime.getRuntime().availableProcessors()`. Shared interval trees (exclude/overlap/distance regions) are read-only, safe to share across threads.
- **Expected speedup**: ~Nx where N = core count, most effective on SSD/NVMe storage

### 1D. Parallelize KL-Divergence Computation (Done)
- **File**: [KullbackLeiblerDistanceBayesianNhmmV5Calculator.java](src/main/java/org/cchmc/epifluidlab/finaleme/hmm/KullbackLeiblerDistanceBayesianNhmmV5Calculator.java)
- **Previous**: Sequential forward-backward on all sequences for distance computation
- **Change**: Parallel execution with `ExecutorService`, local `Double` futures accumulated after completion

---

## TIER 2: I/O and Algorithm Fixes (High Impact, Low Risk)

### 2A. Single-Pass Matrix File Processing (Done)
- **File**: [FinaleMe.java](src/main/java/org/cchmc/epifluidlab/finaleme/hmm/FinaleMe.java) `processMatrixFile()`
- **Previous**: Read the (often multi-GB gzipped) matrix file TWICE -- first for summary statistics, then for data loading
- **Change**: Single pass -- store raw parsed values in memory during stats computation (they end up in memory anyway via `matrixProcess`), then compute z-scores from the in-memory data
- **Achieved speedup**: 2x for this step (eliminates one full gzip decompression + I/O pass)

### 2B. Eliminate Full BAM Scan for Read Counting (Done)
- **File**: [CpgMultiMetricsStats.java](src/main/java/org/cchmc/epifluidlab/finaleme/utils/CpgMultiMetricsStats.java)
- **Previous**: When `-totalReadsInBam` not provided, iterated through the ENTIRE BAM to count reads
- **Change**: Use BAM index statistics (`BAMIndexMetaData.getAlignedRecordCount()`) for fast approximate count with full-scan fallback
- **Achieved speedup**: Saves 30-60 min for a typical 30x WGS BAM

### 2C. Cache `flatPair` in Baum-Welch (Done)
- **File**: [BaumWelchBayesianNhmmV5ScaledLearner.java](src/main/java/org/cchmc/epifluidlab/finaleme/hmm/BaumWelchBayesianNhmmV5ScaledLearner.java)
- **Previous**: `CcInferenceUtils.flatPair(sequences)` called inside `for (int i = 0; i < hmm.nbStates(); i++)` loop -- created a new list of ALL observations per state per iteration
- **Change**: Compute once before the loop, reuse the list

---

## TIER 3: Memory and Data Structure Optimizations

### 3A. Replace Boxed Types with Primitives in HMM Model (Done)
- **File**: [BayesianNhmmV5.java](src/main/java/org/cchmc/epifluidlab/finaleme/hmm/BayesianNhmmV5.java)
- **Previous**: `TreeMap<Integer, Double[]>` for pi, `TreeMap<Integer, Double[][]>` for transitions
- **Change**: Direct array indexing: `double[][] pi = new double[nbCpgDistStates+1][nbStates]`, `double[][][] a = new double[nbCpgDistStates+1][nbStates][nbStates]`. Eliminates TreeMap boxing/unboxing overhead, autoboxing, and hash lookups. Updated serialVersionUID for incompatible format change.

### 3B. Reduce `allGamma` Memory Footprint (Done)
- **File**: [BaumWelchBayesianNhmmV5ScaledLearner.java](src/main/java/org/cchmc/epifluidlab/finaleme/hmm/BaumWelchBayesianNhmmV5ScaledLearner.java)
- **Previous**: Stored gamma arrays for ALL sequences simultaneously in `allGamma[sequences.size()][][]`
- **Change**: Eliminated `allGamma` array entirely. During sequential accumulation of xi results, gamma is recomputed on-the-fly from xi for each sequence. Only `firstGamma[seq][state]` is stored for pi computation, and `pdfWeights[state][totalObs]` accumulates PDF weights incrementally. This avoids storing both allXi and allGamma simultaneously.

### 3C. Streaming Read Processing in Feature Extraction (Not Yet Implemented)
- **File**: [CpgMultiMetricsStats.java](src/main/java/org/cchmc/epifluidlab/finaleme/utils/CpgMultiMetricsStats.java)
- **Current**: Per-CpG BAM index query with new `HashMap<String, SAMRecord>` at each position
- **Planned Change**: Sliding window approach -- maintain active reads as CpG position advances, converting O(N*M) index queries to O(R) total read operations
- **Risk**: High -- fundamental algorithm change requiring careful validation

---

## TIER 4: Architecture and Maintenance

### 4A. Consolidate CpgMultiMetricsStats Variants (Not Yet Implemented)
- Four near-identical classes (~3,559 lines total): `CpgMultiMetricsStats`, `CpgMultiMetricsStatsInMemory`, `CpgMultiMetricsStatsNoBam`, `CpgMultiMetricsStatsV2`
- Extract abstract base class with shared logic; each variant overrides reference/BAM access strategy

### 4B. Dependency and Java Version Upgrades (Done)
- **log4j 1.2.17**: Replaced with SLF4J 2.0.9 + Logback 1.4.14 across all 23+ Java files
- **Java target**: Upgraded from Java 1.8 to Java 21 (maven-compiler-plugin 3.11.0, source/target 21)
- **GATK 3.3**: Removed entirely. `BaseUtils` replaced with inline `basesAreEqual()` and `simpleReverseComplement()` in `BaseUtilsMore.java`. Deleted unused `NotProperPairedReadFilter` and `InvertedDupsReadFilter`.
- **htsjdk 1.141**: Kept as-is (upgrade requires API migration across all BAM processing code)
- **JUnit 5**: Added junit-jupiter 5.10.1 + maven-surefire-plugin 3.2.3

### 4C. Add Unit Tests for HMM Core (Done)
- **File**: [BayesianNhmmV5Test.java](src/test/java/org/cchmc/epifluidlab/finaleme/hmm/BayesianNhmmV5Test.java)
- 13 tests covering: model properties, pi/transition get/set, clone independence, toString, opdf initialization, forward-backward probability, Viterbi decoding, Baum-Welch single iteration (validates probability constraints), KL-divergence (zero for identical models)

### 4D. Remove Dead Code (Done)
- Removed unused method `decodeHmmAssignMethyState` (~130 lines)
- Removed unused methods `estimateXiWithMethy` and 4-parameter `estimateGamma` (~100 lines)
- Removed ~320 lines of commented-out code across all files

---

## Implementation Order (for remaining work)

| Phase | Work | Effort |
|-------|------|--------|
| 1 | Streaming read processing (3C) | 1-2 weeks |
| 2 | Code consolidation (4A) | 1-2 weeks |
| 3 | htsjdk upgrade (4B partial) | 2-3 weeks |

---

## Verification

For every change:
1. Run on small test dataset (single chromosome subset)
2. Compare prediction output within floating-point tolerance
3. Compare trained model parameters (pi, arij, opdf) within tolerance
4. For parallelization: verify `nThreads=1` produces identical output to sequential version
5. Full dataset regression test comparing methylation prediction accuracy

---

## Critical Files

- [FinaleMe.java](src/main/java/org/cchmc/epifluidlab/finaleme/hmm/FinaleMe.java) -- matrix processing, training loop, decode loop
- [BaumWelchBayesianNhmmV5ScaledLearner.java](src/main/java/org/cchmc/epifluidlab/finaleme/hmm/BaumWelchBayesianNhmmV5ScaledLearner.java) -- Baum-Welch iteration (main hotspot)
- [CpgMultiMetricsStats.java](src/main/java/org/cchmc/epifluidlab/finaleme/utils/CpgMultiMetricsStats.java) -- feature extraction (BAM processing)
- [BayesianNhmmV5.java](src/main/java/org/cchmc/epifluidlab/finaleme/hmm/BayesianNhmmV5.java) -- HMM model data structures
- [KullbackLeiblerDistanceBayesianNhmmV5Calculator.java](src/main/java/org/cchmc/epifluidlab/finaleme/hmm/KullbackLeiblerDistanceBayesianNhmmV5Calculator.java) -- convergence check
- [ViterbiBayesianNhmmV5Calculator.java](src/main/java/org/cchmc/epifluidlab/finaleme/hmm/ViterbiBayesianNhmmV5Calculator.java) -- Viterbi decoding
