# FinaleMe Tutorial

This tutorial contains the full, detailed usage guide for FinaleMe v0.60.

If you want the shortest path to run the pipeline, start from [README.md](../README.md).

## 1. What FinaleMe does

FinaleMe predicts CpG methylation from cfDNA fragment features derived from BAM/CRAM or tabix-indexed fragment files. The standard workflow has five steps:

1. Build feature matrix (`CpgFeatureMatrixBuilder`)
2. Train HMM model (`FinaleMe`)
3. Decode methylation (`FinaleMe`)
4. Optional legacy conversion to bigWig (Perl helper)
5. Tissues-of-origin deconvolution (UXM or original methylation-density workflow)

## 2. Installation and reference setup

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
- builds `target/FinaleMe-0.60-jar-with-dependencies.jar` (if missing)
- downloads hg19/hg38 reference data into `data/`

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
./scripts/setup_references.sh summary
```

## 2.3 Reference files used in examples

- `data/hg19.2bit`
- `data/hg19.chrom.sizes`
- `data/CG_motif.hg19.common_chr.pos_only.bedgraph.gz`
- `data/wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed`
- `data/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw`
- `data/CpG_index.hg19.bed.gz`

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
JAR="target/FinaleMe-0.60-jar-with-dependencies.jar"

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
  -inferMethyFromValueWig \
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
| `-fragMethyColumn` | 1-based methylation column index (`m/u`) (0=infer/default). |
| `-fragBaseQ` | Synthetic base quality for fragment mode (default 60). |
| `-defaultMethyStat` | Default methylation state if not provided/inferred (`m` or `u`). |
| `-inferMethyFromValueWig` | Infer methylation from first `-valueWigs` track (`>=50 => m`). |

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
java -Xmx40G -cp "$JAR" \
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
| `-chromSizeFile` | Required with `-bwOutput`. |
| `-bedGraphToBigWig` | Path to UCSC converter executable. |
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

## 8. Step 5: Tissues-of-origin analysis

## 8.1 Option A (recommended): UXM deconvolution

Install tools:

```bash
git clone https://github.com/nloyfer/wgbs_tools.git
cd wgbs_tools
python setup.py

cd ..
git clone https://github.com/nloyfer/UXM_deconv.git
cd UXM_deconv
pip install -r requirements.txt
```

Run deconvolution:

```bash
uxm deconv results/BH01.decode.prediction.pat.gz \
  -o results/BH01.uxm_result.csv \
  -a /path/to/UXM_deconv/supplemental/Atlas.U25.l4.hg19.tsv
```

Notes:

- `data/CpG_index.hg19.bed.gz` and `data/CpG_index.hg38.bed.gz` from setup are equivalent to wgbstools generated CpG index files.
- `wgbstools init_genome` is only needed if you want to regenerate references manually.

## 8.2 Option B: Original methylation-density/QP workflow

This workflow uses bigWig aggregation and quadratic programming.

```bash
ls *WGS.FinaleMe.mincg7.merged.cov.b37.bw | perl -ne 'chomp;$cov=$_;$m=$cov;$m=~s/cov/methy_count/;print " -bigWig $m -useMean0 0 -regionMode 0 -bigWig $cov -useMean0 0 -regionMode 0";' >> cfdna.methy_summary.cmd.txt

cat cfdna.methy_summary.cmd.txt | perl -ne 'chomp;@f=split " -useMean0 0 -regionMode 0";for($i=1,$j=0;$j<=$#f;$j+=2,$i++){$name=$f[$j];$name=~s/ -bigWig (\S+)\S+methy_count.b37.bw/$1/;print "$i\t$name\n";}' > cfdna.names_order.txt

perl -e '$cmd=`cat cfdna.methy_summary.cmd.txt`;chomp($cmd); `java -Xmx10G -cp "lib/dnaaseUtils-0.14-jar-with-dependencies.jar:lib/java-genomics-io.jar:lib/igv.jar" main.java.edu.mit.compbio.utils.AlignMultiWigInsideBed autosome_1kb_intervals.UCSC.cpgIsland_plus_shore.b37.bed output.add_value.methy.bed.gz $cmd`;'
```

Reference R script: `src/R/TissueOfOriginExampleScript.R`

## 9. Troubleshooting

## 9.1 `ClassNotFoundException` on old model files

If you decode an old model trained before package migration, use the current v0.60 jar. Backward-compatible class-name remapping is implemented for legacy serialized model class names.

## 9.2 No `bedGraphToBigWig` in PATH

Install UCSC binary matching your platform, or pass explicit path with `-bedGraphToBigWig`.

## 9.3 Missing CpG index for `-patOutput`

Use setup-provided files:

- hg19: `data/CpG_index.hg19.bed.gz`
- hg38: `data/CpG_index.hg38.bed.gz`

## 9.4 High memory usage in decode

Use `-decodeModeOnly` streaming mode with moderate `-Xmx` and appropriate `-t`.

## 9.5 Chromosome naming mismatch in bigWig

Use:

- `-bwStripChrPrefix`
- `-bwConvertChrMToMT`

as needed for your chrom-size naming convention.

## 10. Performance and reproducibility notes

- Step 1 is parallelized by 5Mb genomic bins (`-t` controls worker count).
- Training and decode are parallelized in FinaleMe (`-t` controls worker count).
- Use fixed `-seed` for reproducible randomized operations.

## 11. Backward-compatible command alias

`edu.northwestern.epifluidlab.finaleme.utils.CpgMultiMetricsStats` is kept as a deprecated alias to `CpgFeatureMatrixBuilder` for script compatibility.

## 12. References

- FinaleMe paper: https://doi.org/10.1038/s41467-024-47196-6
- Reference data (Zenodo): https://doi.org/10.5281/zenodo.19392525
- wgbstools: https://github.com/nloyfer/wgbs_tools
- UXM_deconv: https://github.com/nloyfer/UXM_deconv
