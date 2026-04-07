# FinaleMe

FinaleMe (FragmentatIoN AnaLysis of cEll-free DNA Methylation) predicts CpG methylation from cfDNA WGS fragment features using an HMM pipeline.

## Citation

Liu Y# et al. (2024) FinaleMe: Predicting DNA methylation by the fragmentation patterns of plasma cell-free DNA. Nature Communications. https://doi.org/10.1038/s41467-024-47196-6

## System requirements

- Java 21 or later (Oracle JDK 21 recommended): https://www.oracle.com/java/technologies/downloads/#java21
- Apache Maven 3.8+: 
https://maven.apache.org/install.html (or use the build .jar file from release)
- Optional for bigWig conversion speed/compatibility: `bedGraphToBigWig` (UCSC tools from here: https://hgdownload.soe.ucsc.edu/admin/exe/ and modify $PATH to allow the direct usage). If missing, FinaleMe falls back to built-in Java BigWig writing.
- Optional for custom Step 4 tissue reference-map generation: `wgbstools` (https://github.com/nloyfer/wgbs_tools) and `UXM_deconv` (https://github.com/nloyfer/UXM_deconv)

## Quick install

```bash
git clone https://github.com/epifluidlab/FinaleMe.git
cd FinaleMe
./scripts/sync-vendored-repo.sh
mvn clean package
```

## Quick setup

Run one command to build FinaleMe and download required hg19/hg38 reference files into `data/`:

```bash
./scripts/setup_references.sh
```

## Test dataset

Download a chr22 BAM test file (a ~100X cfDNA WGS data from Snyder et al. 2016 Cell paper):

```bash
mkdir -p test results
curl -L "https://zenodo.org/records/6914806/files/BH01.chr22.bam?download=1" -o test/BH01.chr22.bam
curl -L "https://zenodo.org/records/6914806/files/BH01.chr22.bam.bai?download=1" -o test/BH01.chr22.bam.bai || samtools index test/BH01.chr22.bam
```

## Getting started

Set a jar variable once:

```bash
JAR="target/FinaleMe-0.61-jar-with-dependencies.jar"
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
  -useNoChrPrefixBam \
  -wgsMode \
  -t 4
```

Output: `results/BH01.cpg_features.hg19.bed.gz`. 

It will cost ~25 min for the test dataset.

### Step 2: Train HMM model

```bash
java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.hmm.FinaleMe \
  results/BH01.FinaleMe.model \
  results/BH01.cpg_features.hg19.bed.gz \
  results/BH01.train.prediction.bed.gz \
  -miniDataPoints 7 -gmm -covOutlier 3 -t 4
```

Outputs: model `results/BH01.FinaleMe.model` and training prediction file.  

It will cost < 1 min for the test dataset.

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
```

Outputs:
- `results/BH01.decode.prediction.bed.gz`
- `results/BH01.decode.prediction.methy.bw`
- `results/BH01.decode.prediction.cov.bw`
- `results/BH01.decode.prediction.methy_count.bw`. 

It will cost < 1 min for the test dataset.

### Step 4: Tissues-of-origin analysis

Using FinaleMe's beta-value deconvolution (`BetaValueDeconvolution`, recommended):

```bash
java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.utils.BetaValueDeconvolution \
  -binarizeThreshold 0.1 \
  -markerRegions results/cgi_shore_atlas/Atlas.CGI_shore.U250.l3.hg19.tsv \
  -refBetas results/cgi_shore_atlas/reference_wgbs/betas/beta_list.txt \
  -refGroups results/cgi_shore_atlas/groups_fixed.csv \
  -cpgIndex data/CpG_index.hg19.bed.gz \
  -solver NNLS \
  -output results/BH01.deconv.beta.tsv \
  results/BH01.decode.prediction.bed.gz
```

Use [tutorial/tutorial_ref_maps.md](tutorial/tutorial_ref_maps.md) to generate the marker atlas (`-markerRegions`) and reference panel files (`-refBetas`, `-refGroups`).

## Full tutorial

For full option-by-option documentation, file format details, advanced workflows (including tabix fragment input files from FinaleDB), and troubleshooting, see:

- [tutorial/tutorial.md](tutorial/tutorial.md)

## License

For academic research, please refer to MIT license. For commercial usage, please contact the authors.
