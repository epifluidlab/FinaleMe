# FinaleMe

FinaleMe (FragmentatIoN AnaLysis of cEll-free DNA Methylation) is a Java program to predict DNA methylation in deep and low-coverage cell-free DNA WGS data without other training data.

## Citation

Cite our paper:

Liu Y# et al. (2024) FinaleMe: Predicting DNA methylation by the fragmentation patterns of plasma cell-free DNA. Nature Communications doi: [https://doi.org/10.1038/s41467-024-47196-6](https://doi.org/10.1038/s41467-024-47196-6)


## Installation

### System requirements:

- Java 21 or later — download Oracle JDK 21 from https://www.oracle.com/java/technologies/downloads/#java21 (select your platform: macOS arm64/x64, Linux x64/aarch64). On HPC clusters without root access, use the tarball and set `JAVA_HOME` and `PATH`.
- Apache Maven 3.8+ (only if you need to compile the source code from scratch)
- Perl (tested in v5.26.3), bedGraphToBigWig from UCSC tools (tested in v4), bedtools (tested in v2.29.2). (only if you need to convert predicted methylation level to big wig files)
- R (tested in v4.2.1) and quadprog package. (only if you need to perform tissues-of-origin analysis)

### Build from source:

    ./scripts/sync-vendored-repo.sh
    mvn clean package

This produces `target/FinaleMe-0.60-jar-with-dependencies.jar`, which is self-contained for FinaleMe runtime usage (no extra `lib/*.jar` entries needed in `-cp` for Steps 1-3 below).

`./scripts/sync-vendored-repo.sh` materializes a Maven-layout local repository under `lib-repo/` from `lib/*.jar`, so the build uses normal Maven dependencies instead of `systemPath`.

Or use the precompiled `target/FinaleMe-0.60-jar-with-dependencies.jar` (from Releases/build) directly for FinaleMe Steps 1-3.

### Quick setup of reference data

Use the setup helper to build FinaleMe and download all required reference files for both hg19 and hg38:

```bash
./scripts/setup_references.sh
```

By default files are written to `data/` under this repo. To use a different location:

```bash
export FINALEME_DATA_DIR=/path/to/finaleme_data
./scripts/setup_references.sh
```

Run individual setup stages:

```bash
./scripts/setup_references.sh deps
./scripts/setup_references.sh build
./scripts/setup_references.sh genomes
./scripts/setup_references.sh chromsizes
./scripts/setup_references.sh cpg
./scripts/setup_references.sh darkregions
./scripts/setup_references.sh methylation
./scripts/setup_references.sh summary
```

Reference files downloaded by the setup script:

| File | Purpose | Used in |
|------|---------|---------|
| `hg19.2bit`, `hg38.2bit` | Reference genome in UCSC 2bit format | Step 1 |
| `hg19.chrom.sizes`, `hg38.chrom.sizes` | Chromosome size files for bigWig conversion | Step 3 (`-bwOutput`) |
| `CG_motif.hg19.common_chr.pos_only.bedgraph.gz`, `CG_motif_seqkit.pos_only.hg38.bedgraph.gz` | CpG motif coordinates | Step 1 |
| `CpG_index.hg19.bed.gz`, `CpG_index.hg38.bed.gz` (+ `.csi`) | CpG index for pat/beta outputs | Step 3 (`-patOutput`) |
| `wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed`, `hg38-blacklist.v2.bed` | Dark/blacklist regions | Step 1 filtering |
| `wgbs_buffyCoat_jensen2015GB.methy.hg19.bw`, `wgbs_buffyCoat_jensen2015GB.methy.hg38.bw` | Methylation prior tracks | Step 1 (`-valueWigs`) |

Java recommendation: Oracle JDK 21 from https://www.oracle.com/java/technologies/downloads/#java21. Choose the matching platform package. On clusters without root access, download the tarball and set `JAVA_HOME` and `PATH`.

### Running tests:

    mvn test

This runs 13 JUnit 5 unit tests covering the HMM core: model properties, forward-backward probability, Viterbi decoding, Baum-Welch training iteration, and KL-divergence convergence.

### Other required data:

> **Tip:** Run `./scripts/setup_references.sh` to download all files below automatically for both hg19 and hg38. See [Quick setup](#quick-setup-of-reference-data).

- methylation prior file in standard big wig format (use `wgbs_buffyCoat_jensen2015GB.methy.hg19.bw` and `wgbs_buffyCoat_jensen2015GB.methy.hg38.bw` from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19392525.svg)](https://doi.org/10.5281/zenodo.19392525)). Or generate your own from WGBS data in healthy buffy coat (our data source: Jensen et al. 2015 Genome Biology).
- bed files to mask dark regions in the genome (use `wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed` and `hg38-blacklist.v2.bed` from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19392525.svg)](https://doi.org/10.5281/zenodo.19392525)). Or provide your own dark region files for other references.
- Chromosome sizes: can be obtained from the FASTA file of the reference genome: `samtools faidx input.fa`. See this [example](https://github.com/epifluidlab/cragr/blob/3d419a49/inst/extdata/human_g1k_v37.chrom.sizes).
- CG_motif.bedgraph: bedgraph file with CpG coordinates in the reference genome (**forward strand only**, one entry per CpG dinucleotide). Each CpG dinucleotide (5'-CG-3') is listed once at the position of the C on the forward strand. For hg19, this yields ~28M positions. Both forward- and reverse-strand reads are captured at the same CpG position during feature extraction, so no fragments are lost. Use `CG_motif.hg19.common_chr.pos_only.bedgraph.gz` (hg19) and `CG_motif_seqkit.pos_only.hg38.bedgraph.gz` (hg38) from Zenodo.

  To generate this file from a reference FASTA:
  ```
  # Extract all CG dinucleotide positions (forward-strand C position, 0-based BED format)
  python3 -c "
  from pysam import FastaFile
  ref = FastaFile('hg19.fa')
  for chrom in ref.references:
      seq = ref.fetch(chrom).upper()
      for i in range(len(seq) - 1):
          if seq[i] == 'C' and seq[i+1] == 'G':
              print(f'{chrom}\t{i}\t{i+1}\t1')
  " > CG_motif.hg19.bedgraph
  ```
  Or download pre-built files from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19392525.svg)](https://doi.org/10.5281/zenodo.19392525).

- hg19.2bit / hg38.2bit: binary version of reference genome, which can be downloaded from [UCSC genome browser](http://hgdownload.soe.ucsc.edu/goldenPath/) or converted from .fasta files by [faToTwoBit](https://github.com/ENCODE-DCC/kentUtils)

### Small test input data
- bam files from chr22 in healthy individuals can be downloaded here [https://zenodo.org/records/6914806/files/BH01.chr22.bam?download=1](https://zenodo.org/records/6914806/files/BH01.chr22.bam?download=1)

## What's new in v0.60

- **Java 21 required** (upgraded from Java 1.8)
- **Parallelized feature extraction** by 5Mb genomic bins for multi-core speedup
- **Parallelized HMM training** (Baum-Welch forward-backward) and Viterbi decoding
- **Parallelized KL-divergence** convergence check
- **Single-pass matrix file processing** (was reading the file twice)
- **Fast BAM read counting** via index statistics (no more full BAM scan)
- **Streaming BAM feature extraction** within each 5Mb bin (active-read window; no per-CpG BAM index query)
- **Feature extraction progress reporting** with processed/total CpGs and percentage
- **TwoBitParser debug logs disabled by default** (via bundled Logback config)
- **Reduced memory usage** by replacing boxed types with primitive arrays and eliminating allGamma storage
- **Replaced log4j** with SLF4J + Logback (fixes CVE-2021-4104)
- **Removed GATK dependency** (replaced with inline utilities)
- **Upgraded htsjdk** from 1.141 to 2.24.1
- **Consolidated Cpg feature-matrix builders** using a shared abstract base and common interval-loading utilities
- **Clean packaging model**: normal Maven dependencies (no `systemPath`) + single self-contained runtime jar for Steps 1-3
- **Added unit tests** for HMM core (13 tests)

See [PLAN.md](PLAN.md) for full details on all optimizations.

> **Note:** Models trained with v0.59 or earlier are not compatible with v0.60 due to the internal data structure change from `TreeMap<Integer, Double[]>` to primitive arrays. You will need to retrain models with v0.60.

## Getting started

### Input
- coordinate-sorted and indexed bam file
- or bgzipped/tabix-indexed fragment BED/TSV file with at least `chr`, `start`, `end`, and `strand`

### Start analysis

The analysis consists of several steps:

1. Calculate the feature file for each CpG in each DNA fragment from an indexed bam file.
2. Train HMM model from the feature file (only utilize fragments with >=7 CpGs).
3. Predict the methylation status at each CpG.
4. Convert the prediction result into a bigwig file for the visualization and analysis.
5. Perform tissues-of-origin analysis by predicted DNA methylation level in cfDNA and methylome from reference panel.

#### Step 1: extract features from bam files for the training and decoding
```
java -Xmx20G -cp "target/FinaleMe-0.60-jar-with-dependencies.jar" \
  edu.northwestern.epifluidlab.finaleme.utils.CpgFeatureMatrixBuilder \
  hg19.2bit \
  CG_motif.hg19.common_chr.pos_only.bedgraph \
  CG_motif.hg19.common_chr.pos_only.bedgraph \
  input.bam \
  CpgFeatureMatrixBuilder.hg19.details.bed.gz \
  -stringentPaired \
  -excludeRegions wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed \
  -valueWigs methyPrior:0:wgbs_buffyCoat_jensen2015GB.methy.hg19.bw \
  -wgsMode
```

`edu.northwestern.epifluidlab.finaleme.utils.CpgMultiMetricsStats` is still available as a deprecated alias for backward-compatible scripts.

Alternative Step 1 input mode (tabix fragment BED/TSV):
```
java -Xmx20G -cp "target/FinaleMe-0.60-jar-with-dependencies.jar" \
  edu.northwestern.epifluidlab.finaleme.utils.CpgFeatureMatrixBuilder \
  hg19.2bit \
  CG_motif.hg19.common_chr.pos_only.bedgraph \
  CG_motif.hg19.common_chr.pos_only.bedgraph \
  fragments.bed.gz \
  CpgFeatureMatrixBuilder.hg19.details.bed.gz \
  -fragmentInputTabix \
  -fragStrandColumn 4 \
  -valueWigs methyPrior:0:wgbs_buffyCoat_jensen2015GB.methy.hg19.bw \
  -inferMethyFromValueWig
```

Notes for tabix fragment mode:
- If methylation state is already present in the fragment file, provide `-fragMethyColumn` (1-based column index, values `m/u`).
- If methylation state is absent, the builder can infer it from the first `-valueWigs` track (commonly methylation prior).
- If neither is available, `-defaultMethyStat` is used (`u` by default).

Feature extraction is now parallelized by 5Mb genomic bins and processed with a per-bin streaming BAM read window, automatically using all available CPU cores.
Runtime progress is logged as CpGs processed out of total with percent completion.
You can set thread count explicitly with `-t` (for example, `-t 8`).

* CG_motif.hg19.common_chr.pos_only.bedgraph is the bedgraph file with forward-strand CpG coordinates (see "Other required data" above for how to generate it). It can be downloaded here [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19392525.svg)](https://doi.org/10.5281/zenodo.19392525)
* hg19.2bit is the binary input of reference genome, which can be downloaded from [UCSC genome browser](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) or converted from .fasta files by [faToTwoBit](https://github.com/ENCODE-DCC/kentUtils)

#### Step 2: train the model 
```
java -Xmx40G -cp "target/FinaleMe-0.60-jar-with-dependencies.jar" \
  edu.northwestern.epifluidlab.finaleme.hmm.FinaleMe \
  test.FinaleMe.mincg7.model \
  CpgFeatureMatrixBuilder.hg19.details.bed.gz \
  test.FinaleMe.mincg7.prediction.bed.gz \
  -miniDataPoints 7 -gmm -covOutlier 3
```

Training uses parallelized Baum-Welch (forward-backward computed across sequences in parallel) and parallelized KL-divergence convergence checking.
Set FinaleMe thread count with `-t` (for example, `-t 8`).

#### Step 3: decode and make the prediction of CpG methylation level
```
java -Xmx20G -cp "target/FinaleMe-0.60-jar-with-dependencies.jar" \
  edu.northwestern.epifluidlab.finaleme.hmm.FinaleMe \
  test.FinaleMe.mincg7.model \
  CpgFeatureMatrixBuilder.hg19.details.bed.gz \
  test.FinaleMe.mincg7.prediction.bed.gz \
  -decodeModeOnly
```

Viterbi decoding is parallelized across fragments.
Set FinaleMe thread count with `-t` (for example, `-t 8`).

You can also generate bigWig files directly during decode (so Step 4 is optional):
```
  -bwOutput \
  -chromSizeFile hg19.chrom.sizes
```
This produces:
- `*.methy.bw` (predicted methylation percent)
- `*.cov.bw` (predicted total coverage)
- `*.methy_count.bw` (predicted methylated count)

If chromosome naming differs between decode output and your chrom size file, use:
```
  -bwStripChrPrefix \
  -bwConvertChrMToMT
```

To emit UXM-compatible per-fragment outputs for deconvolution (`.pat.gz` + `.beta`), add:
```
  -patOutput \
  -cpgIndexFile data/CpG_index.hg19.bed.gz
```
For hg38, use `-cpgIndexFile data/CpG_index.hg38.bed.gz`.
These pre-built CpG index files are downloaded by `./scripts/setup_references.sh cpg` and are equivalent to the `wgbs_tools/references/*/CpG.bed.gz` files generated by `wgbstools init_genome`.

This writes:
- `*.pat.gz`: bgzip-compressed, tab-separated `chr`, `global_cpg_index`, `C/T pattern`, `count` — compatible with [UXM_deconv](https://github.com/nloyfer/UXM_deconv) and [wgbstools](https://github.com/nloyfer/wgbs_tools)
- `*.beta`: binary uint8 pairs `(methylated_count, total_count)` for each CpG index row


Then run UXM deconvolution directly:
```
uxm deconv test.FinaleMe.mincg7.prediction.pat.gz \
  -o test.FinaleMe.mincg7.prediction.uxm_result.csv \
  -a /path/to/UXM_deconv/supplemental/Atlas.U25.l4.hg19.tsv
```

#### Step 4 (optional): convert predicted result to .bw file for the visualization in genome browser
```
perl src/perl/bedpredict2bw.b37.pl test test.FinaleMe.mincg7.prediction.bed.gz
```

#### Step 5: tissues-of-origin analysis

There are two approaches for tissues-of-origin deconvolution:

##### Option A: UXM deconvolution (recommended)

This approach uses the per-fragment methylation patterns from the `.pat.gz` file generated in Step 3 (with `-patOutput`). It leverages the [UXM_deconv](https://github.com/nloyfer/UXM_deconv) method, which classifies each fragment as Unmethylated (U), miXed (X), or Methylated (M) and deconvolves against a reference atlas.

**Prerequisites — install wgbstools and UXM_deconv:**

1. Install [wgbstools](https://github.com/nloyfer/wgbs_tools) (required by UXM_deconv for pat file processing):
   ```
   git clone https://github.com/nloyfer/wgbs_tools.git
   cd wgbs_tools
   python setup.py
   ```

2. Install [UXM_deconv](https://github.com/nloyfer/UXM_deconv):
   ```
   git clone https://github.com/nloyfer/UXM_deconv.git
   cd UXM_deconv
   pip install -r requirements.txt
   ```


**Run deconvolution:**
```
uxm deconv test.FinaleMe.mincg7.prediction.pat.gz \
  -o test.FinaleMe.mincg7.prediction.uxm_result.csv \
  -a /path/to/UXM_deconv/supplemental/Atlas.U25.l4.hg19.tsv
```

The output CSV contains estimated cell-type proportions for each sample. Reference atlases for hg19 and hg38 are included in the `UXM_deconv/supplemental/` directory. See the [UXM_deconv README](https://github.com/nloyfer/UXM_deconv#readme) for additional options (e.g., `--ignore`, `--include` for selecting specific cell types).

##### Option B: quadratic programming with methylation density (original method)

This approach uses the bigWig files from Step 4 to compute methylation density in genomic windows and deconvolves using quadratic programming against a reference methylome panel.

* Generate methylation density in 1kb window at autosomes across all available methylation prediction files:
```
ls *WGS.FinaleMe.mincg7.merged.cov.b37.bw | perl -ne 'chomp;$cov=$_;$m=$cov;$m=~s/cov/methy_count/;print " -bigWig $m -useMean0 0 -regionMode 0 -bigWig $cov -useMean0 0 -regionMode 0";' >> cfdna.methy_summary.cmd.txt

cat cfdna.methy_summary.cmd.txt | perl -ne 'chomp;@f=split " -useMean0 0 -regionMode 0";for($i=1,$j=0;$j<=$#f;$j<=$#f;$j+=2,$i++){$name=$f[$j];$name=~s/ -bigWig (\S+)\S+methy_count.b37.bw/$1/;print "$i\t$name\n";}' > cfdna.names_order.txt

perl -e '$cmd=`cat cfdna.methy_summary.cmd.txt`;chomp($cmd); `java -Xmx10G -cp "lib/dnaaseUtils-0.14-jar-with-dependencies.jar:lib/java-genomics-io.jar:lib/igv.jar" main.java.edu.mit.compbio.utils.AlignMultiWigInsideBed autosome_1kb_intervals.UCSC.cpgIsland_plus_shore.b37.bed output.add_value.methy.bed.gz $cmd`;'
```

* R script is available within src/R/TissueOfOriginExampleScript.R
* The reference methylomes used in the paper can be downloaded here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19392525.svg)](https://doi.org/10.5281/zenodo.19392525)

## License

For academic research, please refer to MIT license. For commerical usage, please contact the authors.
