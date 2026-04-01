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
| 3C. Streaming Read Processing in Feature Extraction | Done |
| 4A. Consolidate CpG Feature-Matrix Builder Variants | Done |
| 4B. Dependency and Java Version Upgrades | Done |
| 4C. Add Unit Tests for HMM Core | Done |
| 4D. Remove Dead Code | Done |
| 4E. Update README for v0.60 | Done |
| 4F. Clean Packaging and Dependency Resolution | Done |

---

## TIER 1: Parallelization (Highest Impact)

### 1A. Parallelize Viterbi Decoding (Done)
- **File**: [FinaleMe.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/FinaleMe.java) `decodeHmm()` (~line 997)
- **Previous**: Sequential `for(int j=0; j < matrix.size(); j++)` loop, each fragment independently decoded
- **Change**: Use `ExecutorService` with `ForkJoinPool` to decode fragments in parallel. The HMM model is read-only during decode. Collect results into a pre-sized array, then aggregate `methySummary` sequentially afterward.
- **Expected speedup**: ~Nx where N = core count (8-16x typical)

### 1B. Parallelize Baum-Welch Forward-Backward (Done)
- **File**: [BaumWelchBayesianNhmmV5ScaledLearner.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/BaumWelchBayesianNhmmV5ScaledLearner.java) `iterate()`
- **Previous**: Sequential loop over all sequences computing forward-backward, xi, gamma, then accumulating into shared `aijNum/aijDen/arij`
- **Change**: Map-reduce pattern -- parallel forward-backward/xi/gamma computation across sequences, then sequential accumulation of results into transition matrices
- **Expected speedup**: ~Nx per Baum-Welch iteration (up to 50 iterations)

### 1C. Parallelize Feature Extraction by 5Mb Genomic Bins (Done)
- **File**: [CpgFeatureMatrixBuilder.java](src/main/java/edu/northwestern/epifluidlab/finaleme/utils/CpgFeatureMatrixBuilder.java)
- **Previous**: Sequential `for(String chr : cpgCollections.keySet())` with per-CpG BAM index queries
- **Change**: Divide genome into 5Mb non-overlapping bins. Each bin processed in parallel with its own `SamReader`, `TwoBitParser`, `TabixFeatureReader`, and `BigWigFileReader` instances. Each bin writes to a temp file, then results are concatenated. Thread pool sized to `Runtime.getRuntime().availableProcessors()`. Shared interval trees (exclude/overlap/distance regions) are read-only, safe to share across threads.
- **Expected speedup**: ~Nx where N = core count, most effective on SSD/NVMe storage

### 1D. Parallelize KL-Divergence Computation (Done)
- **File**: [KullbackLeiblerDistanceBayesianNhmmV5Calculator.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/KullbackLeiblerDistanceBayesianNhmmV5Calculator.java)
- **Previous**: Sequential forward-backward on all sequences for distance computation
- **Change**: Parallel execution with `ExecutorService`, local `Double` futures accumulated after completion

---

## TIER 2: I/O and Algorithm Fixes (High Impact, Low Risk)

### 2A. Single-Pass Matrix File Processing (Done)
- **File**: [FinaleMe.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/FinaleMe.java) `processMatrixFile()`
- **Previous**: Read the (often multi-GB gzipped) matrix file TWICE -- first for summary statistics, then for data loading
- **Change**: Single pass -- store raw parsed values in memory during stats computation (they end up in memory anyway via `matrixProcess`), then compute z-scores from the in-memory data
- **Achieved speedup**: 2x for this step (eliminates one full gzip decompression + I/O pass)

### 2B. Eliminate Full BAM Scan for Read Counting (Done)
- **File**: [CpgFeatureMatrixBuilder.java](src/main/java/edu/northwestern/epifluidlab/finaleme/utils/CpgFeatureMatrixBuilder.java)
- **Previous**: When `-totalReadsInBam` not provided, iterated through the ENTIRE BAM to count reads
- **Change**: Use BAM index statistics (`BAMIndexMetaData.getAlignedRecordCount()`) for fast approximate count with full-scan fallback
- **Achieved speedup**: Saves 30-60 min for a typical 30x WGS BAM

### 2C. Cache `flatPair` in Baum-Welch (Done)
- **File**: [BaumWelchBayesianNhmmV5ScaledLearner.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/BaumWelchBayesianNhmmV5ScaledLearner.java)
- **Previous**: `CcInferenceUtils.flatPair(sequences)` called inside `for (int i = 0; i < hmm.nbStates(); i++)` loop -- created a new list of ALL observations per state per iteration
- **Change**: Compute once before the loop, reuse the list

---

## TIER 3: Memory and Data Structure Optimizations

### 3A. Replace Boxed Types with Primitives in HMM Model (Done)
- **File**: [BayesianNhmmV5.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/BayesianNhmmV5.java)
- **Previous**: `TreeMap<Integer, Double[]>` for pi, `TreeMap<Integer, Double[][]>` for transitions
- **Change**: Direct array indexing: `double[][] pi = new double[nbCpgDistStates+1][nbStates]`, `double[][][] a = new double[nbCpgDistStates+1][nbStates][nbStates]`. Eliminates TreeMap boxing/unboxing overhead, autoboxing, and hash lookups. Updated serialVersionUID for incompatible format change.

### 3B. Reduce `allGamma` Memory Footprint (Done)
- **File**: [BaumWelchBayesianNhmmV5ScaledLearner.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/BaumWelchBayesianNhmmV5ScaledLearner.java)
- **Previous**: Stored gamma arrays for ALL sequences simultaneously in `allGamma[sequences.size()][][]`
- **Change**: Eliminated `allGamma` array entirely. During sequential accumulation of xi results, gamma is recomputed on-the-fly from xi for each sequence. Only `firstGamma[seq][state]` is stored for pi computation, and `pdfWeights[state][totalObs]` accumulates PDF weights incrementally. This avoids storing both allXi and allGamma simultaneously.

### 3C. Streaming Read Processing in Feature Extraction (Done)
- **File**: [CpgFeatureMatrixBuilder.java](src/main/java/edu/northwestern/epifluidlab/finaleme/utils/CpgFeatureMatrixBuilder.java)
- **Previous**: Per-CpG BAM index query (`queryOverlapping`) with new `HashMap<String, SAMRecord>` at each CpG position
- **Change**: Implemented per-bin streaming read processing. Each 5Mb worker now scans BAM once for the bin, maintains an active-read sliding window, and evaluates CpGs in coordinate order. This removes per-CpG BAM index seeks while preserving per-CpG deduplication and methylation-calling logic.
- **Result**: O(R) BAM iterator pass per bin (plus per-overlap scoring), with significantly reduced random-index overhead on large BAMs.

---

## TIER 4: Architecture and Maintenance

### 4A. Consolidate CpG Feature-Matrix Builder Variants (Done)
- Added shared abstract base class: [AbstractCpgMultiMetricsStats.java](src/main/java/edu/northwestern/epifluidlab/finaleme/utils/AbstractCpgMultiMetricsStats.java)
- Refactored `CpgFeatureMatrixBuilder`, `CpgMultiMetricsStats` (deprecated alias), and `CpgMultiMetricsStatsV2` to extend the shared base.
- Moved shared interval-loading logic (gz/plain BED readers, exclude-region loading, stranded CpG interval loading) into the base class and reused across all four variants.

### 4B. Dependency and Java Version Upgrades (Done)
- **log4j 1.2.17**: Replaced with SLF4J 2.0.9 + Logback 1.4.14 across all 23+ Java files
- **Java target**: Upgraded from Java 1.8 to Java 21 (maven-compiler-plugin 3.11.0, source/target 21)
- **GATK 3.3**: Removed entirely. `BaseUtils` replaced with inline `basesAreEqual()` and `simpleReverseComplement()` in `BaseUtilsMore.java`. Deleted unused `NotProperPairedReadFilter` and `InvertedDupsReadFilter`.
- **htsjdk**: Upgraded from 1.141 to 2.24.1 in `pom.xml` and removed legacy HTTP Maven repository configuration.
- **JUnit 5**: Added junit-jupiter 5.10.1 + maven-surefire-plugin 3.2.3

### 4C. Add Unit Tests for HMM Core (Done)
- **File**: [BayesianNhmmV5Test.java](src/test/java/edu/northwestern/epifluidlab/finaleme/hmm/BayesianNhmmV5Test.java)
- 13 tests covering: model properties, pi/transition get/set, clone independence, toString, opdf initialization, forward-backward probability, Viterbi decoding, Baum-Welch single iteration (validates probability constraints), KL-divergence (zero for identical models)

### 4D. Remove Dead Code (Done)
- Removed unused method `decodeHmmAssignMethyState` (~130 lines)
- Removed unused methods `estimateXiWithMethy` and 4-parameter `estimateGamma` (~100 lines)
- Removed ~320 lines of commented-out code across all files

### 4E. Update README for v0.60 (Done)
- Updated Java requirement from 1.8 to Java 21, Maven from 3.8.6 to 3.8+
- Added and maintained "Build from source" documentation for the current packaging model
- Added "Running tests" section (`mvn test`, 13 JUnit 5 tests)
- Added "What's new in v0.60" changelog with model incompatibility warning
- Updated all `java -cp` commands for v0.60 and single-fat-jar usage in core FinaleMe Steps 1-3.
- Added notes about parallelization in Steps 1-3

### 4F. Clean Packaging and Dependency Resolution (Done)
- Removed all `systemPath` dependencies from `pom.xml` and switched to normal Maven dependencies.
- Added vendored Maven repository support at `lib-repo/` (file repository) and sync script [scripts/sync-vendored-repo.sh](scripts/sync-vendored-repo.sh) to materialize Maven-layout artifacts from `lib/*.jar`.
- Enabled standard `mvn package` assembly execution to always produce `target/FinaleMe-0.60-jar-with-dependencies.jar`.
- Added project-local Maven cache configuration (`.mvn/maven.config`) for reproducible builds without relying on `~/.m2`.
- Updated runtime commands for Steps 1-3 to use only the single fat jar in classpath.

---

## Implementation Order (Completed)

All previously remaining items have been implemented:
- Streaming read processing (3C)
- CpG feature-matrix builder consolidation (4A)
- htsjdk upgrade (4B partial)
- clean packaging and dependency resolution (4F)

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

- [FinaleMe.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/FinaleMe.java) -- matrix processing, training loop, decode loop
- [BaumWelchBayesianNhmmV5ScaledLearner.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/BaumWelchBayesianNhmmV5ScaledLearner.java) -- Baum-Welch iteration (main hotspot)
- [CpgFeatureMatrixBuilder.java](src/main/java/edu/northwestern/epifluidlab/finaleme/utils/CpgFeatureMatrixBuilder.java) -- feature extraction (BAM processing)
- [BayesianNhmmV5.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/BayesianNhmmV5.java) -- HMM model data structures
- [KullbackLeiblerDistanceBayesianNhmmV5Calculator.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/KullbackLeiblerDistanceBayesianNhmmV5Calculator.java) -- convergence check
- [ViterbiBayesianNhmmV5Calculator.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/ViterbiBayesianNhmmV5Calculator.java) -- Viterbi decoding
