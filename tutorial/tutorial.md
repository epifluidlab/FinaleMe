# FinaleMe Tutorial

This tutorial contains the full, detailed usage guide for FinaleMe v0.63.

If you want the shortest path to run the pipeline, start from [README.md](../README.md).

## 1. What FinaleMe does

FinaleMe predicts CpG methylation from cfDNA fragment features derived from BAM/CRAM or tabix-indexed fragment files, then deconvolves the predictions into tissue-of-origin (TOO) cell-type fractions. The standard workflow has six steps:

1. Build feature matrix (`CpgFeatureMatrixBuilder`)
2. Train HMM model (`FinaleMe`)
3. Decode methylation (`FinaleMe`)
4. Optional legacy conversion to bigWig (Perl helper)
5. Tissues-of-origin deconvolution per sample (`BetaValueDeconvolution` recommended; UXM compatibility optional)
6. Cohort-level differential cell-type analysis (`scripts/too_diff_analysis.py`)

## 2. Installation and reference setup

Default pipeline requirement summary:

- Steps 1-5 in this tutorial require only FinaleMe (Java + jar build) and the pre-built reference panels downloaded by `setup_references.sh`.
- Step 6 (differential analysis) requires Python 3.9+ with `pandas`, `numpy`, `statsmodels`.
- `wgbstools` and `UXM_deconv` are only needed for the optional custom atlas generation workflow (see [tutorial/tutorial_ref_maps.md](tutorial_ref_maps.md)).

## 2.1 Install source

```bash
git clone https://github.com/epifluidlab/FinaleMe.git
cd FinaleMe
```

## 2.2 One-command setup (recommended)

```bash
./scripts/setup_references.sh
```

This command:

- checks dependencies
- builds `target/FinaleMe-0.63-jar-with-dependencies.jar` (if missing)
- downloads hg19/hg38 reference data into `data/`, including the pre-built TOO reference panels for Step 5

Use a custom data directory with:

```bash
export FINALEME_DATA_DIR=/path/to/finaleme_data
./scripts/setup_references.sh
```

Useful subcommands:

```bash
./scripts/setup_references.sh deps
./scripts/setup_references.sh build
./scripts/setup_references.sh genomes
./scripts/setup_references.sh chromsizes
./scripts/setup_references.sh cpg
./scripts/setup_references.sh darkregions
./scripts/setup_references.sh methylation
./scripts/setup_references.sh atlas         # TOO deconvolution reference panels (hg19/hg38)
./scripts/setup_references.sh summary
```

## 2.3 Reference files used in examples

- `data/hg19.2bit`
- `data/hg19.chrom.sizes`
- `data/CG_motif.hg19.common_chr.pos_only.bedgraph.gz`
- `data/wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed`
- `data/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw`
- `data/CpG_index.hg19.bed.gz`
- `data/Atlas.CGI_shore.U250.l3.hg19.tsv` (TOO reference panel for Step 5; hg38 also provided)
- `data/Atlas.CGI_shore.U250.l3.hg38.tsv`
- `data/Atlas.pluse_microglia_astrocyte.CGI_shore.U250.l3.hg38.tsv` (hg38 panel + microglia/astrocyte)

## 2.4 Small test BAM

```bash
mkdir -p test results
curl -L "https://zenodo.org/records/6914806/files/BH01.chr22.bam?download=1" -o test/BH01.chr22.bam
curl -L "https://zenodo.org/records/6914806/files/BH01.chr22.bam.bai?download=1" -o test/BH01.chr22.bam.bai || samtools index test/BH01.chr22.bam
```

## 3. Inputs and expected formats

## 3.1 Step 1 positional arguments

`CpgFeatureMatrixBuilder` usage:

```bash
CpgFeatureMatrixBuilder [opts] hg19.2bit cpg_list.bed all_cpg.bed wgs.bam|fragments.tsv.gz cpg_detail.txt.gz
```

Meaning:

- argument 1: reference `.2bit`
- argument 2: target CpG list to process
- argument 3: all CpGs list (used when `-includeCpgDist` is enabled)
- argument 4: BAM/CRAM or bgzipped+tabix fragment BED/TSV
- argument 5: output feature matrix (`.gz` recommended)

## 3.2 Supported fragment input modes

- BAM/CRAM mode (default)
- tabix fragment mode (`-fragmentInputTabix` or auto-detected by file extension)

For tabix fragment BED/TSV, the file must include at least `chr`, `start`, `end`, and strand information.

## 4. Step 1: Feature extraction (detailed)

## 4.1 Standard BAM command

```bash
JAR="target/FinaleMe-0.63-jar-with-dependencies.jar"

java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.utils.CpgFeatureMatrixBuilder \
  data/hg19.2bit \
  data/CG_motif.hg19.common_chr.pos_only.bedgraph.gz \
  data/CG_motif.hg19.common_chr.pos_only.bedgraph.gz \
  test/BH01.chr22.bam \
  results/BH01.cpg_features.hg19.bed.gz \
  -stringentPaired \
  -excludeRegions data/wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed \
  -valueWigs methyPrior:0:data/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw \
  -useNoChrPrefixBam \
  -wgsMode \
  -t 4
```

## 4.2 Tabix fragment mode command

```bash
java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.utils.CpgFeatureMatrixBuilder \
  data/hg19.2bit \
  data/CG_motif.hg19.common_chr.pos_only.bedgraph.gz \
  data/CG_motif.hg19.common_chr.pos_only.bedgraph.gz \
  fragments.bed.gz \
  results/fragments.cpg_features.hg19.bed.gz \
  -fragmentInputTabix \
  -fragStrandColumn 4 \
  -valueWigs methyPrior:0:data/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw \
  -t 4
```

## 4.3 Step 1 options

### Core filters and behavior

| Option | Description |
|---|---|
| `-minBaseQ` | Minimum base quality threshold (default 5). |
| `-minMapQ` | Minimum mapping quality threshold (default 30). |
| `-maxFragLen` | Maximum fragment length considered (default 500). |
| `-maxDistToFragEnd` | Max allowed distance from CpG to fragment end (default 250). |
| `-maxCov` | Maximum coverage threshold per CpG (default 250). |
| `-totalReadsInBam` | Override auto-estimated total reads/fragments for normalization. |
| `-wgsMode` | Enable WGS mode (non-bisulfite-space behavior). |
| `-skipSecondEnd` | Ignore read2 in paired-end statistics. |
| `-stringentPaired` | Keep only properly oriented read pairs. |
| `-includeCpgDist` | Add nearest-CpG distance feature column. |
| `-excludeRegions` | BED file(s) of excluded intervals. |
| `-useNoChrPrefixBam` | Use BAM contig naming without `chr` prefix. |
| `-t` | Number of threads for parallel 5Mb bins. |

### Additional track-derived features

| Option | Description |
|---|---|
| `-overlapRegions track:file` | Add overlap flag(s) against BED track(s). |
| `-distantRegions track:file` | Add nearest-distance feature(s) to interval track(s). |
| `-valueWigs track:ext:file` | Add averaged value feature(s) from bigWig around CpGs. |
| `-valueBeds track:ext:file` | Add averaged value feature(s) from tabix BED around CpGs. |

### Sequence/k-mer features

| Option | Description |
|---|---|
| `-kmerLen` | Auto-generate all k-mers up to this length. |
| `-kmerString` | Provide explicit k-mer list file. |
| `-kmerExt` | +/- region for k-mer extraction around CpG (default 100). |
| `-useFragBaseKmer` | Compute k-mer from fragment sequence context. |
| `-useStrandSpecificFragBase` | Strand-aware fragment k-mer mode. |

### Tabix fragment mode specific

| Option | Description |
|---|---|
| `-fragmentInputTabix` | Force tabix fragment mode. |
| `-fragStrandColumn` | 1-based strand column index (0=auto). |
| `-fragNameColumn` | 1-based fragment-name column index (0=auto/synthetic). |
| `-fragMethyColumn` | 1-based methylation column index (`m`/`u`) (0=use `-defaultMethyStat`). Rare; use only if your tabix file carries per-fragment methylation calls. |
| `-fragBaseQ` | Synthetic base quality for fragment mode (default 60). |
| `-defaultMethyStat` | Value written to the `methy_stat` output column when `-fragMethyColumn` is unset and no inline `m`/`u` token is detected. Default: `m`. The HMM does NOT use `methy_stat` as a label or input — it is unsupervised over the feature vector — so this only affects the AUC/QC reporting columns. The default `m` matches BAM behavior under `-wgsMode`. |

#### How the total fragment count is computed in tabix mode

The `Norm_Frag_cov` column normalizes per-CpG coverage by the total record count. In tabix mode the count is determined as follows, in priority order:

1. `-totalReadsInBam <N>`: skip estimation entirely and use `N`.
2. **BGZF sampling** (default): open the bgzipped file with `BlockCompressedInputStream`, read the first 50,000 data rows, extract the compressed-byte offset from the BGZF virtual file pointer, compute compressed bytes per line, then extrapolate to the full compressed file size. Typical cost: tens of milliseconds even on multi-GB files.
3. **Full scan** (fallback): the legacy line-by-line counter. Used only when BGZF sampling fails (plain-gzip header, very small file, or the sample fits inside the first BGZF block).

Sampling accuracy on a synthetic 5M-row file: estimate 5,148,318 vs exact 5,000,000 (3% over) in 42 ms vs 1311 ms full-scan (31× speedup). For exact counts pass `-totalReadsInBam` directly.

#### `Norm_Frag_cov` is fragment-level by default in v0.63+

Both BAM and tabix modes now produce **fragment-level** coverage by default:

| Mode | Numerator at each CpG | Denominator |
|---|---|---|
| BAM (default v0.63+) | unique fragments overlapping the CpG (deduplicated by readName) | total filtered fragments (estimated from sample's unique-readName ratio) |
| tabix | fragments overlapping the CpG (1 row = 1 fragment) | total fragments in the tabix file |

Both ratios produce the same `Norm_Frag_cov` scale on the same biological sample, so models trained on BAM transfer cleanly to tabix-fragment input and vice versa.

#### Legacy read-level scale (`-coverageReadLevel`)

For models trained on BAM input **before v0.63**, the historical scale was read-level:

- Numerator = raw SAMRecords overlapping the CpG (PE: each end counted separately when both overlap)
- Denominator = filtered raw reads (≈ 2 × fragments for PE)

Pass `-coverageReadLevel` to recover this legacy behavior when re-using a pre-v0.63 BAM-trained model. With this flag:

- Pre-v0.63 BAM-trained models keep working unchanged on BAM input.
- The flag is **not honored in tabix fragment mode** (tabix files have no read-level information; the count is inherently fragment-level).

Cross-mode practice:

| Trained on | Decode input | Recommended flags |
|---|---|---|
| BAM, v0.63+ (default fragment-level) | BAM **or** tabix | none (default works) |
| BAM, pre-v0.63 (read-level) | BAM | `-coverageReadLevel` |
| BAM, pre-v0.63 (read-level) | tabix | retrain or rescale; the legacy read-level scale cannot be reproduced from tabix |
| Tabix, any version | BAM **or** tabix | none (default works) |

#### `-useEndMotif` and `-saveMotifLookup` in tabix mode

`-useEndMotif` is **fully supported** in tabix mode. The 5' end 4-mer is extracted from the reference at `(fragment.start, fragment.end)` (with strand awareness), cached per `readName` so the same fragment at multiple CpGs reuses the same motif, and either the motif string (training mode) or the lookup score (`-loadMotifLookup`) is written as the trailing output column.

`-saveMotifLookup` is **NOT supported** in tabix mode. The lookup table is built as `methylated_count / total_count` per 4-mer, which requires per-CpG bisulfite ground truth. Tabix fragment input has no such info — `methy_stat` is constant `'m'` from `-defaultMethyStat` — so the resulting lookup would be degenerate (score 1.0 for every 4-mer). The script therefore **errors out at startup** with a clear message:

```
-saveMotifLookup is not supported with tabix fragment input. Tabix fragment
files lack per-CpG bisulfite information, so methy_stat is constant and the
resulting motif lookup would be degenerate (score=1.0 for every 4-mer). Train
the motif lookup on a paired WGBS BAM, then re-use it on tabix/WGS runs via
-loadMotifLookup.
```

The same restriction applies to BAM input under `-wgsMode` (every covered CpG is labeled `'m'` because the read base is the unconverted reference). Standard recipe:

1. Train the motif lookup once on a paired WGBS/bisulfite-converted BAM:
   ```bash
   java ... CpgFeatureMatrixBuilder ... wgbs.bam wgbs.cpg_features.bed.gz \
     -useEndMotif -saveMotifLookup motif_lookup.tsv
   ```
2. Re-use it on tabix or WGS-mode runs:
   ```bash
   java ... CpgFeatureMatrixBuilder ... fragments.bed.gz cpg_features.bed.gz \
     -fragmentInputTabix -useEndMotif -loadMotifLookup motif_lookup.tsv
   ```

## 4.4 Step 1 output format (`*.cpg_features*.bed.gz`)

Header starts with:

```text
chr	start	end	readName	FragLen	Frag_strand	methy_stat	Norm_Frag_cov	baseQ	Offset_frag	Dist_frag_end
```

Optional columns can follow in this order:

1. `dist_nearest_CpG` (if `-includeCpgDist`)
2. one column per `-overlapRegions` track
3. one column per `-distantRegions` track
4. one column per `-valueBeds` track
5. one column per `-valueWigs` track
6. one column per k-mer feature

Key fields:

- `methy_stat`: observed methylation label (`m`/`u`) per CpG record
- `Norm_Frag_cov`: normalized fragment coverage feature
- `Offset_frag`: CpG offset index within fragment
- `Dist_frag_end`: minimum distance to fragment ends

## 5. Step 2: Train HMM (detailed)

## 5.1 Training command

```bash
java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.hmm.FinaleMe \
  results/BH01.FinaleMe.model \
  results/BH01.cpg_features.hg19.bed.gz \
  results/BH01.train.prediction.bed.gz \
  -miniDataPoints 7 \
  -gmm \
  -covOutlier 3 \
  -t 4
```

This writes a serialized model file (`.model`) and a prediction table.

## 5.2 Training-related options (`FinaleMe`)

| Option | Description |
|---|---|
| `-states` | Number of hidden states (even number expected). |
| `-features` | Number of features per observation vector. |
| `-miniDataPoints` | Minimum CpGs per fragment to include. |
| `-maxCpgs` | Maximum CpGs per fragment to include. |
| `-maxFragLen` | Maximum fragment-length state bound. |
| `-minFragLen` | Minimum fragment-length threshold. |
| `-maxCpgDist` | Max CpG distance used for transition bins. |
| `-bin` | Distance bin size for non-homogeneous priors/transitions. |
| `-covOutlier` | Outlier filter by z-score in feature loading. |
| `-gmm` | Initialize HMM using GMM. |
| `-wgbs` | WGBS-oriented initialization mode. |
| `-iteration` | Max Baum-Welch iterations. |
| `-tol` | Convergence tolerance. |
| `-decayRate` | Relative convergence threshold. |
| `-tolKmeans` | K-means tolerance used by initialization. |
| `-decayKmeans` | K-means decay criterion. |
| `-mixNumberInFeature` | Mixture count(s) for Gaussian emissions. |
| `-bayesianFactor` | Prior weighting factor in decoding/training. |
| `-cpgNumClip` | Clip for CpG-count scaling in HMM. |
| `-methylatedState` | Which state is interpreted as methylated (0/1). |
| `-seed` | Random seed (`<0` for non-deterministic). |
| `-t` | Parallel worker count for training/decoding internals. |

## 6. Step 3: Decode methylation (detailed)

## 6.1 Decode command

```bash
java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.hmm.FinaleMe \
  results/BH01.FinaleMe.model \
  results/BH01.cpg_features.hg19.bed.gz \
  results/BH01.decode.prediction.bed.gz \
  -decodeModeOnly \
  -t 4 \
  -bwOutput \
  -chromSizeFile data/hg19.chrom.sizes \
  -patOutput \
  -cpgIndexFile data/CpG_index.hg19.bed.gz
```

## 6.2 Decode/output options (`FinaleMe`)

| Option | Description |
|---|---|
| `-decodeModeOnly` | Skip training and decode directly with existing model. |
| `-decodeP` | Decision criterion used by Viterbi methylation labeling. |
| `-randomPerm` | Randomize labels from prior instead of trained HMM. |
| `-lowCoverage` | Low-coverage mode with alternate feature handling. |
| `-region` | Decode only within BED intervals. |
| `-exclude` | Exclude BED intervals from decode. |
| `-patOutput` | Write UXM-compatible `.pat.gz` and `.beta` outputs. |
| `-cpgIndexFile` | Required CpG index for `-patOutput` (use `data/CpG_index.*.bed.gz`). |
| `-bwOutput` | Write decode summary bigWig outputs. |
| `-chromSizeFile` | Required with `-bwOutput` (used by both UCSC converter and Java fallback writer). |
| `-bedGraphToBigWig` | Path to UCSC converter executable; if missing, FinaleMe auto-falls back to Java BigWig writer. |
| `-bwStripChrPrefix` | Remove `chr` prefix in bigWig conversion. |
| `-bwConvertChrMToMT` | Convert chrM/M naming to MT. |
| `-t` | Parallel decoding thread count. |

### AUC mode option

| Option | Description |
|---|---|
| `-aucMode` | Compute ROC/AUC-style summaries across decode thresholds. |

`-bwOutput` is not supported together with `-aucMode`.

## 6.3 Step 3 prediction output format (`*.prediction.bed.gz`)

Header:

```text
#chr	start	end	methy_perc_predict	methy_count_predict	total_count_predict	methy_perc_obs	methy_count_obs	total_count_obs
```

Column meaning:

- `methy_perc_predict`: predicted methylation percentage at locus
- `methy_count_predict`: predicted methylated count
- `total_count_predict`: predicted total count
- `methy_perc_obs`: observed methylation percentage from feature input labels
- `methy_count_obs`: observed methylated count
- `total_count_obs`: observed total count

## 6.4 Optional UXM output formats

### `.pat.gz`

Tab-separated rows:

```text
chr	start_cpg_index	CT_pattern	count
```

- `start_cpg_index`: global CpG index of first CpG in fragment pattern
- `CT_pattern`: per-CpG decoded pattern (`C` for methylated, `T` for unmethylated)
- `count`: multiplicity of identical fragment pattern records

### `.beta`

Binary file storing per-index `(methylated_count, total_count)` as uint8 pairs.

## 6.5 Optional bigWig outputs

When `-bwOutput` is enabled, FinaleMe writes:

- `*.methy.bw`: predicted methylation percentage track
- `*.cov.bw`: predicted total count track
- `*.methy_count.bw`: predicted methylated count track

## 7. Step 4: Legacy Perl bigWig workflow

If you prefer the old conversion utility:

```bash
perl src/perl/bedpredict2bw.b37.pl results/BH01 results/BH01.decode.prediction.bed.gz
```

This is optional when Step 3 already runs with `-bwOutput`.

## 8. Step 5: Tissues-of-origin deconvolution (per sample)

## 8.1 Recommended: `BetaValueDeconvolution` with the pre-built atlas

The simplest workflow uses the atlas downloaded by `setup_references.sh`:

```bash
JAR="target/FinaleMe-0.63-jar-with-dependencies.jar"

java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.utils.BetaValueDeconvolution \
  -refPanel data/Atlas.CGI_shore.U250.l3.hg19.tsv \
  -cpgIndex data/CpG_index.hg19.bed.gz \
  -output results/BH01.deconv.beta.tsv \
  results/BH01.decode.prediction.bed.gz
```

For hg38 inputs, use `data/Atlas.CGI_shore.U250.l3.hg38.tsv` (or the extended `Atlas.pluse_microglia_astrocyte.CGI_shore.U250.l3.hg38.tsv` which adds microglia and astrocyte cell types).

### What the defaults do

With the v0.63 production defaults, the command above automatically:

- Binarizes both reference and sample at `-binarizeThreshold 0.1`.
- Solves the deconvolution with `-solver NNLS` (Lawson-Hanson + sum=1 normalization; matches the original method).
- Runs `-bootstrap` (1000 stratified replicates by `target` cell type) for 95% percentile confidence intervals on each cell-type fraction.
- Runs `-permutationTest` (10000 column-shuffles per cell type, `-permutationMode celltype` with pooled null) for frequentist-calibrated p-values.
- Applies `-fdrAlpha 0.05` to flag significant cell types.
- Uses fixed seeds (`-bootstrapSeed 42 -permutationSeed 42`) and 10 threads (`-bootstrapThreads 10`) so reruns are bit-identical and reasonably fast.

Output (long format, one row per `(sample, cell_type)`):

```
sample   cell_type   proportion   CI_lower   CI_upper   p_value   q_value   significant   p_source      n_replicates
BH01     Blood-T     0.156        0.142      0.171      5.3e-05   5.3e-05   YES           permutation   10000
BH01     Adipocyte   0.012        0.000      0.025      0.213     0.213     NO            permutation   10000
...
```

### Multiple samples

Pass several query files in one invocation; they are processed in parallel and emit a single combined long-format TSV:

```bash
java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.utils.BetaValueDeconvolution \
  -refPanel data/Atlas.CGI_shore.U250.l3.hg19.tsv \
  -cpgIndex data/CpG_index.hg19.bed.gz \
  -output results/cohort.deconv.beta.tsv \
  results/SAMPLE_01.decode.prediction.bed.gz \
  results/SAMPLE_02.decode.prediction.bed.gz \
  results/SAMPLE_03.decode.prediction.bed.gz
```

Or run one job per sample (e.g. on a cluster) and concatenate later. Step 6 (`too_diff_analysis.py`) accepts either layout.

## 8.2 Step 5 options (`BetaValueDeconvolution`)

### Reference panel

| Option | Description |
|---|---|
| `-refPanel` | Pre-aggregated atlas TSV with embedded per-cell-type values (recommended; use `data/Atlas.CGI_shore.U250.l3.{hg19,hg38}.tsv` or `data/Atlas.pluse_microglia_astrocyte.CGI_shore.U250.l3.hg38.tsv`). |
| `-markerRegions` | Atlas/marker BED with `startCpG/endCpG` columns. Used when the atlas does not embed per-cell-type values, in combination with `-refBetas`/`-refGroups`. |
| `-refBetas` | Comma-separated reference `.beta` files, or a text file listing one path per line. Required only when `-refPanel` is not used. |
| `-refGroups` | CSV mapping sample names to cell-type groups. Required only with `-refBetas`. |
| `-cpgIndex` | Required CpG index BED for marker coordinate -> CpG mapping. Use `data/CpG_index.{hg19,hg38}.bed.gz`. |
| `-cgiShoreRegions` | Legacy mode: CGI+shore BED file used together with `-refBetas`/`-refGroups` to tile windows on the fly. Not needed when `-refPanel` or `-markerRegions` is provided. |
| `-windowSize` | Window size (bp) for legacy CGI+shore tiling. Default: 1000. |
| `-minCoverage` | Minimum total Cs+Ts per window across reference samples. Default: 10. |
| `-topVariablePercent` | Fraction of most-variable windows kept in legacy mode. Default: 0.01. |
| `-replicateMode` | How to combine replicates: `aggregate` (sum counts, default) or `average` (mean of ratios). |

### Solver and binarization

| Option | Description |
|---|---|
| `-solver` | `NNLS` (Lawson-Hanson + simplex renormalization, default) or `QP` (Goldfarb-Idnani quadratic programming). |
| `-binarizeThreshold` | Methylation density threshold for binarizing both reference and sample before NNLS. Default: 0.1. |

### Bootstrap (CI on the point estimates)

| Option | Description |
|---|---|
| `-bootstrap` | Enable bootstrap CI + per-cell-type stats. **Default: true.** |
| `-nBootstrap` | Bootstrap replicates. Default: 1000. Use >=1000 for tight small p-values. |
| `-ciLevel` | Two-sided CI level. Default: 0.95 (-> percentile CI [2.5%, 97.5%]). |
| `-bootstrapStratified` | Stratify resampling by the `target` cell-type column in the reference panel. Eliminates the inflated-CI bias for cell types with fewer markers. **Default: true.** |
| `-bootstrapThreads` | Parallel workers for bootstrap and permutation. Default: 10. |
| `-bootstrapSeed` | RNG seed for bootstrap. Default: 42. Use `-1` for fresh entropy per run. |

### Permutation (frequentist p-values, BH q-values)

| Option | Description |
|---|---|
| `-permutationTest` | Enable permutation test in addition to bootstrap. **Default: true.** |
| `-permutationMode` | `celltype` (default; at each marker row, shuffle the K cell-type values among columns -- symmetric null suitable for cell-type-specific atlases like U250) or `marker` (per cell type, shuffle the rows of that cell type's reference column). |
| `-nPermutations` | Permutation replicates. Default: 10000. |
| `-permutationSeed` | RNG seed for permutations. Default: 42. |
| `-permutationNullPooled` | Under `-permutationMode celltype`, pool all K x P permuted values into a single null shared across cell types. Achieves p-value resolution of `1/(K*P+1)`. **Default: true.** |
| `-permutationBHCorrect` | Apply BH FDR correction even when the pooled null is in use. Default: false (the pooled null is already calibrated across cell types). |
| `-fdrAlpha` | q-value threshold for the `significant` flag. Default: 0.05. |

### Output

| Option | Description |
|---|---|
| `-output` | Output TSV path. With bootstrap on (the default), long format: one row per `(sample, cell_type)`. With bootstrap off, wide format: one row per cell type, one column per sample. |

Two practical choices that tune sensitivity vs. cost:

- For routine cohort runs, the defaults (B=1000, P=10000, 10 threads) need ~5-30 minutes per sample on a typical 38-cell-type panel.
- To dial things back for a quick scan: `-nBootstrap 200 -nPermutations 1000`. To go more sensitive: `-nPermutations 50000` or `-permutationBHCorrect`.

## 8.3 Output schema

Long-format columns (when bootstrap or permutation is on, which is the default):

| Column | Meaning |
|---|---|
| `sample` | Query file basename. |
| `cell_type` | One of the K cell types in the reference panel. |
| `proportion` | Point estimate from the full-data NNLS solve, normalized to sum to 1.0 across cell types. |
| `CI_lower`, `CI_upper` | Percentile bootstrap CI at `-ciLevel`. |
| `p_value` | Permutation p (default) or bootstrap-tail p if `-permutationTest=false`. |
| `q_value` | BH-adjusted p (or = p when `-permutationNullPooled=true` and `-permutationBHCorrect=false`). |
| `significant` | `YES` if `q_value <= fdrAlpha`, else `NO`. |
| `p_source` | `permutation` or `bootstrap` -- which procedure produced the p-value. |
| `n_replicates` | B (bootstrap) or P (permutation) actually run. |

Wide format (legacy; emitted when `-bootstrap=false`): rows are cell types, columns are samples, values are point estimates.

## 8.4 Optional: build a custom atlas

The pre-built atlas in `data/` covers 38 standard cell types from the Loyfer et al. reference. To build a custom atlas (different cell types, different shore size, etc.), see [tutorial/tutorial_ref_maps.md](tutorial_ref_maps.md). This requires `wgbstools` + `UXM_deconv` and is an advanced workflow.

## 8.5 Optional legacy mode: UXM deconvolution

This requires running Step 3 with `-patOutput`.

```bash
uxm deconv results/BH01.decode.prediction.pat.gz \
  -o results/BH01.uxm_result.csv \
  -a /path/to/UXM_deconv/supplemental/Atlas.U25.l4.hg19.tsv
```

UXM does not produce CIs, p-values, or q-values; for those, use `BetaValueDeconvolution` instead.

## 9. Step 6: Differential TOO analysis (cohort comparison)

The per-sample `significant` flag from Step 5 tests `H0: w_c = 0` *in that sample* -- it is a QC diagnostic for individual samples, **not** a group comparison test. To find cell types that differ between groups (e.g., Disease vs Control), use `scripts/too_diff_analysis.py`, which:

1. Pivots per-sample point estimates into a samples x cell_types matrix.
2. Optionally refines NNLS-floored zeros using the bootstrap `CI_upper`.
3. Applies the centered log-ratio (CLR) transform to handle the simplex constraint.
4. Fits a linear model per cell type: `clr(w_c) ~ group + clinical_covariates + technical_covariates`.
5. Applies BH FDR correction across cell types.

## 9.1 Why not just count "YES" flags across samples?

| Pitfall of per-sample YES counting | Why it matters |
|---|---|
| Significance != presence | A 5% true cell-type fraction may fail q<0.05 in a low-coverage sample purely from power loss; counting "YES" loses this. |
| Loss of magnitude | "Prostate at 8%" and "Prostate at 35%" both register as YES; biologically very different. |
| Wrong null hypothesis | Per-sample q tests "is this fraction non-zero in this sample?", not "do groups differ?". |
| Threshold flicker | Borderline cell types switch YES/NO at q=0.05 due to noise; this generates spurious between-group differences. |

The CLR + linear model approach uses point-estimate **magnitudes** as the data, models between-sample variance correctly, and lets you adjust for confounders like batch and coverage tier.

## 9.2 Inputs

You can feed `too_diff_analysis.py` either:

- A single combined long-format TSV produced by running `BetaValueDeconvolution` on multiple query files in one invocation (`--deconv-tsv results/cohort.deconv.beta.tsv`), or
- Per-sample TSVs from independent invocations (`--deconv-files results/cohort/*.deconv.beta.tsv`); the script auto-concatenates them and infers the `sample` field from the filename if needed.

You also need a metadata TSV with at minimum a `sample` column matching the deconv outputs and a column with the group label:

```
sample      disease_status   age   sex   library_batch
14230_1     Disease          67    M     B1
14230_2     Disease          71    F     B2
HD_45       Control          62    M     B1
HD_46       Control          59    F     B2
...
```

## 9.3 Recommended command

```bash
python scripts/too_diff_analysis.py \
  --deconv-tsv results/cohort.deconv.beta.tsv \
  --metadata samples.tsv \
  --group-col disease_status \
  --reference-group Control \
  --covariates age,sex,library_batch \
  --test omnibus \
  --refine-zeros-with-ci \
  --output results/diff_celltypes.tsv \
  --verbose
```

For two groups (Disease vs Control), `omnibus` is equivalent to a t-test on the group factor. For 3+ groups, `omnibus` performs an F-test (any group differs); use `--test pairwise` with `--reference-group <name>` to get one row per non-reference group.

## 9.4 Step 6 options (`too_diff_analysis.py`)

### Input

| Option | Description |
|---|---|
| `--deconv-tsv` | Single combined long-format TSV from `BetaValueDeconvolution -bootstrap`. |
| `--deconv-files` | Per-sample TSVs (multiple paths or a glob pattern). Use this when each sample was a separate Step 5 job. |
| `--metadata` | Sample metadata TSV with at minimum `sample` and the `--group-col` column. |

### Group definition

| Option | Description |
|---|---|
| `--group-col` | Metadata column holding the group label. Default: `group`. |
| `--reference-group` | Group level to use as the reference for pairwise comparisons. If omitted, the alphabetically first level is used. Required for `--test pairwise`. |
| `--covariates` | Comma-separated metadata columns to adjust for, e.g. `age,sex,library_batch,coverage_tier`. Categorical columns are auto-detected by dtype. |
| `--test` | `omnibus` (default; F-test on the group factor, any group differs) or `pairwise` (each non-reference group vs the reference, one row per `(cell_type, level)`). |

### Compositional preprocessing

| Option | Description |
|---|---|
| `--exclude-cell-types` | Comma-separated cell types to skip (e.g. `Unknown` if the input came from `finaleme-too`; `BetaValueDeconvolution` does not estimate Unknown). |
| `--refine-zeros-with-ci` | Replace exact-zero point estimates with `CI_upper/2` where `CI_upper > --refine-zeros-ci-min`. NNLS sometimes drives small genuine fractions to 0 because of the non-negative constraint; the bootstrap CI can detect "below detection" rather than "absent". |
| `--refine-zeros-ci-min` | Threshold above which a zero point estimate is treated as below-detection. Default: 1e-4. |
| `--clr-eps` | Pseudocount added to zero proportions before log. Default: 1e-6. |
| `--weighted` | Inverse-variance weight samples by mean CI width across cell types. Down-weights low-quality samples. Requires CI columns. |

### Output

| Option | Description |
|---|---|
| `--fdr-alpha` | q-value threshold for the `significant` flag. Default: 0.05. |
| `--output` | Output TSV path. |
| `--verbose` | Print group sizes, drop counts, and a top-significant summary table. |

## 9.5 Output schema

Per-cell-type results, sorted by q-value ascending:

| Column | Meaning |
|---|---|
| `cell_type` | Cell type from the reference panel. |
| `n_samples` | Number of samples that contributed to this cell type's fit (after dropping rows with NaN in group/covariates). |
| `r_squared` | Fraction of CLR variance explained by `group + covariates`. |
| `F_stat` | (omnibus only) F statistic for the group factor. |
| `p_value` | omnibus: F-test p-value. pairwise: t-test p-value for the level vs reference. |
| `q_value` | Benjamini-Hochberg-adjusted p-value across cell types. |
| `significant` | `True` if `q_value <= --fdr-alpha`. |
| `max_effect_clr` | (omnibus only) Maximum-magnitude per-level coefficient in CLR space; positive = elevated in that group vs reference. |
| `max_effect_level` | (omnibus only) Group level corresponding to `max_effect_clr`. |
| `level`, `effect_clr`, `std_err` | (pairwise only) Per-level CLR-space effect and standard error. |

CLR-space effects are interpretable as multiplicative log-ratios: `effect_clr = +1.0` means the cell type's geometric-mean fraction in that group is `e^1 ≈ 2.7x` higher than in the reference, after accounting for covariates and the simplex constraint.

## 9.6 Worked example

```
$ python scripts/too_diff_analysis.py \
    --deconv-tsv results/cohort.deconv.beta.tsv \
    --metadata samples.tsv \
    --group-col disease_status \
    --reference-group Control \
    --covariates age,sex,library_batch \
    --output results/diff_celltypes.tsv --verbose

Samples in deconv: 30; in metadata: 30; intersection used: 30
Group sizes:
disease_status
Control    15
Disease    15
Cell types: 38
Wrote results/diff_celltypes.tsv: 38 row(s), 5 significant at q<=0.05.

Top significant cell types:
   cell_type  max_effect_clr max_effect_level   p_value   q_value  n_samples
    Prostate-Ep        +2.55          Disease  5.4e-15   2.7e-14         30
        Blood-T        -0.98          Disease  2.4e-11   6.0e-11         30
      Liver-Hep        -1.24          Disease  5.5e-06   9.1e-06         30
       Endothel        -0.50          Disease  9.3e-03   1.2e-02         30
   Lung-Ep-Bron        -0.31          Disease  3.5e-02   4.3e-02         30
```

In this example the disease cohort has a strong elevation of Prostate-Ep (+2.55 CLR -> ~12x geometric-mean increase) accompanied by drops in blood and liver fractions, all robust to age/sex/library_batch.

## 10. Troubleshooting

## 10.1 `ClassNotFoundException` on old model files

If you decode an old model trained before package migration, use the current v0.63 jar. Backward-compatible class-name remapping is implemented for legacy serialized model class names.

## 10.2 No `bedGraphToBigWig` in PATH

FinaleMe now auto-falls back to Java BigWig writing when this executable is missing.
Install UCSC `bedGraphToBigWig` if you still prefer/require the UCSC binary path.

## 10.3 Missing CpG index for `-patOutput` or `-cpgIndex`

Use setup-provided files:

- hg19: `data/CpG_index.hg19.bed.gz`
- hg38: `data/CpG_index.hg38.bed.gz`

## 10.4 Memory usage in high coverage WGS data

Try different `-Xmx` and appropriate `-t`. We tested with `-Xmx20G` and `-t 5` for HD_46 dataset (~16X depth), but may need `-Xmx80G` and `-t 5` for 14230_1 dataset (~39X depth) in the paper.

## 10.5 Chromosome naming mismatch in bigWig

Use:

- `-bwStripChrPrefix`
- `-bwConvertChrMToMT`

as needed for your chrom-size naming convention.

## 10.6 BetaValueDeconvolution defaults are too slow

The v0.62 defaults run 1000 bootstrap replicates and 10000 permutations per cell type. For a quick scan use:

```bash
-nBootstrap 200 -nPermutations 1000
```

For a one-off run without inference (legacy wide-format output, fastest):

```bash
-bootstrap=false -permutationTest=false
```

## 10.7 Step 6 differential analysis: "no overlap between deconv and metadata"

`too_diff_analysis.py` matches samples by exact string equality on the `sample` column. If your deconv files have basenames like `BH01.decode.prediction.bed.gz` but your metadata uses `BH01`, either:

- Strip the suffix from the deconv `sample` column, or
- Update the metadata `sample` column to match the deconv basenames exactly.

## 10.8 Step 6 reports "fit failed" for some cell types

This usually means too few non-zero observations to fit the linear model. Check the `n_samples` column. If a cell type is zero in nearly every sample (common for very rare tissues at low coverage), exclude it via `--exclude-cell-types CT1,CT2`.

## 11. Performance and reproducibility notes

- Step 1 is parallelized by 5Mb genomic bins (`-t` controls worker count).
- Training and decode are parallelized in FinaleMe (`-t` controls worker count).
- Step 5 bootstrap and permutation are parallelized via `-bootstrapThreads`.
- All randomness in Step 5 is seeded by default (`-bootstrapSeed 42 -permutationSeed 42`); reruns at any thread count produce bit-identical output. Set `-1` for fresh entropy.
- Use fixed `-seed` for reproducible randomized operations in Steps 2-3.

## 12. Backward-compatible command alias

`edu.northwestern.epifluidlab.finaleme.utils.CpgMultiMetricsStats` is kept as a deprecated alias to `CpgFeatureMatrixBuilder` for script compatibility.

## 13. References

- FinaleMe paper: https://doi.org/10.1038/s41467-024-47196-6
- Reference data (Zenodo): https://doi.org/10.5281/zenodo.19392525
- TOO reference panels (Zenodo): https://doi.org/10.5281/zenodo.19742408
- wgbstools: https://github.com/nloyfer/wgbs_tools
- UXM_deconv: https://github.com/nloyfer/UXM_deconv
