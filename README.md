# FinaleMe

FinaleMe (FragmentatIoN AnaLysis of cEll-free DNA Methylation) predicts CpG methylation from cfDNA WGS fragment features using an HMM pipeline.

## Citation

Liu Y# et al. (2024) FinaleMe: Predicting DNA methylation by the fragmentation patterns of plasma cell-free DNA. Nature Communications. https://doi.org/10.1038/s41467-024-47196-6

## System requirements

- Java 21 or later (Oracle JDK 21 recommended): https://www.oracle.com/java/technologies/downloads/#java21
- For source build only: Apache Maven 3.8+
- Optional for bigWig conversion: `bedGraphToBigWig` (UCSC tools)
- Optional for Step 5 tissues-of-origin: `wgbstools` and `UXM_deconv`

## Quick install

```bash
git clone https://github.com/epifluidlab/FinaleMe.git
cd FinaleMe
```

## Quick setup

Run one command to build FinaleMe and download required hg19/hg38 reference files into `data/`:

```bash
./scripts/setup_references.sh
```

If you want a different reference directory:

```bash
export FINALEME_DATA_DIR=/path/to/finaleme_data
./scripts/setup_references.sh
```

## Small test dataset

Download a chr22 BAM test file:

```bash
mkdir -p test results
curl -L "https://zenodo.org/records/6914806/files/BH01.chr22.bam?download=1" -o test/BH01.chr22.bam
curl -L "https://zenodo.org/records/6914806/files/BH01.chr22.bam.bai?download=1" -o test/BH01.chr22.bam.bai || samtools index test/BH01.chr22.bam
```

## Getting started

Set a jar variable once:

```bash
JAR="target/FinaleMe-0.60-jar-with-dependencies.jar"
```

### Step 1: Build CpG feature matrix

```bash
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

Output: `results/BH01.cpg_features.hg19.bed.gz`

### Step 2: Train HMM model

```bash
java -Xmx40G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.hmm.FinaleMe \
  results/BH01.FinaleMe.model \
  results/BH01.cpg_features.hg19.bed.gz \
  results/BH01.train.prediction.bed.gz \
  -miniDataPoints 7 -gmm -covOutlier 3 -t 4
```

Outputs: model `results/BH01.FinaleMe.model` and training prediction file.

### Step 3: Decode CpG methylation

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

Outputs:
- `results/BH01.decode.prediction.bed.gz`
- `results/BH01.decode.prediction.methy.bw`
- `results/BH01.decode.prediction.cov.bw`
- `results/BH01.decode.prediction.methy_count.bw`
- `results/BH01.decode.prediction.pat.gz`
- `results/BH01.decode.prediction.beta`

### Step 4: Optional legacy bigWig conversion script

If you prefer the older Perl-based conversion workflow:

```bash
perl src/perl/bedpredict2bw.b37.pl results/BH01 results/BH01.decode.prediction.bed.gz
```

### Step 5: Tissues-of-origin analysis

Using UXM deconvolution:

```bash
uxm deconv results/BH01.decode.prediction.pat.gz \
  -o results/BH01.uxm_result.csv \
  -a /path/to/UXM_deconv/supplemental/Atlas.U25.l4.hg19.tsv
```

## Full tutorial

For full option-by-option documentation, file format details, advanced workflows (including tabix fragment input), and troubleshooting, see:

- [tutorial/tutorial.md](tutorial/tutorial.md)

## License

For academic research, please refer to MIT license. For commercial usage, please contact the authors.
