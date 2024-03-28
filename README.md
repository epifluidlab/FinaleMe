# FinaleMe

FinaleMe (FragmentatIoN AnaLysis of cEll-free DNA Methylation) is a Java program to predict DNA methylation in deep and low-coverage cell-free DNA WGS data without other training data.

## Citation

Cite our paper:

Liu Y# et al. (2024) FinaleMe: Predicting DNA methylation by the fragmentation patterns of plasma cell-free DNA. Nature Communications doi: [https://doi.org/10.1038/s41467-024-47196-6](https://doi.org/10.1038/s41467-024-47196-6)


## Installation

### System requirements:

- Java (tested in openjdk-1.8.0 at Linux CentOS 8 (x64 platform) and Mac OSX 12.4 (aarch64 platform))
- Apache Maven (tested in v3.8.6, only if you need to compile the source code from scratch)
- Perl (tested in v5.26.3), bedGraphToBigWig from UCSC tools (tested in v4), bedtools (tested in v2.29.2). (only if you need to convert predicted methylation level to big wig files)
- R (tested in v4.2.1) and quadprog package. (only if you need to perform tissues-of-origin analysis)

### Quick installation:

    mvn compile assembly:single

Or use the precompiled .jar file from Releases and other .jars from lib/ directory without installation

### Other required data:

- methylation prior file in standard big wig format (can use wgbs_buffyCoat_jensen2015GB.methy.hg19.bw file directly from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7647046.svg)](https://doi.org/10.5281/zenodo.7647046)). Or you can generate your own by using WGBS data in buffycoat from healthy individuals (our data is from Jensen et al. 2015 Genome Biology paper)
- bed files to mask the Dark regions in the genome (can use wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed files directly from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7647046.svg)](https://doi.org/10.5281/zenodo.7647046)). Or you can download these dark region files for other reference genome.
- Chromosome sizes: can be obtained from the FASTA file of the reference genome: `samtools faidx input.fa`. See this [example](https://github.com/epifluidlab/cragr/blob/3d419a49/inst/extdata/human_g1k_v37.chrom.sizes).
- CG_motif.bedgraph: bedgraph file with CpG's coordinate in the reference genome
- hg19.bit: binary version of reference genome, which can be downloaded from [UCSC genome browser](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) or converted from .fastq files by [faToTwoBit](https://github.com/ENCODE-DCC/kentUtils)

### Small test input data
- bam files from chr22 in healthy individuals can be downloaded here [https://zenodo.org/records/6914806/files/BH01.chr22.bam?download=1](https://zenodo.org/records/6914806/files/BH01.chr22.bam?download=1)

## Getting started

### Input
- coordinate-sorted and indexed bam file

### Start analysis

The analysis consists of several steps:

1. Calculate the feature file for each CpG in each DNA fragment from an indexed bam file.
2. Train HMM model from the feature file (only utilize fragments with >=7 CpGs).
3. Predict the methylation status at each CpG.
4. Convert the prediction result into a bigwig file for the visualization and analysis.
5. Perform tissues-of-origin analysis by predicted DNA methylation level in cfDNA and methylome from reference panel.

#### Step 1: extract features from bam files for the training and decoding
```
java -Xmx20G -cp "target/FinaleMe-0.58-jar-with-dependencies.jar:lib/gatk-package-distribution-3.3.jar:lib/sis-jhdf5-batteries_included.jar:lib/java-genomics-io.jar:lib/igv.jar" org.cchmc.epifluidlab.finaleme.utils.CpgMultiMetricsStats hg19.2bit CG_motif.hg19.common_chr.pos_only.bedgraph CG_motif.hg19.common_chr.pos_only.bedgraph input.bam CpgMultiMetricsStats.hg19.details.bed.gz -stringentPaired -excludeRegions wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed -valueWigs methyPrior:0:wgbs_buffyCoat_jensen2015GB.methy.hg19.bw -wgsMode
```

* CG_motif.hg19.common_chr.pos_only.bedgraph is the bedgraph file with CpG's coordinate in the reference genome
* hg19.bit is the binary input of reference genome, which can be downloaded from [UCSC genome browser](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) or converted from .fastq files by [faToTwoBit](https://github.com/ENCODE-DCC/kentUtils)

#### Step 2: train the model 
```
java -Xmx100G -cp "target/FinaleMe-0.58-jar-with-dependencies.jar:lib/jahmm-0.6.2.jar" org.cchmc.epifluidlab.finaleme.hmm.FinaleMe test.FinaleMe.mincg7.model CpgMultiMetricsStats.hg19.details.bed.gz test.FinaleMe.mincg7.prediction.bed.gz -miniDataPoints 7 -gmm -covOutlier 3
```
#### Step 3: decode and make the prediction of CpG methylation level
```
java -Xmx100G -cp "target/FinaleMe-0.58-jar-with-dependencies.jar:lib/jahmm-0.6.2.jar" org.cchmc.epifluidlab.finaleme.hmm.FinaleMe test.FinaleMe.mincg7.model CpgMultiMetricsStats.hg19.details.bed.gz test.FinaleMe.mincg7.prediction.bed.gz -decodeModeOnly
```

#### Step 4: convert predicted result to .bw file for the visualization in genome browser
```
perl src/perl/bedpredict2bw.b37.pl test test.FinaleMe.mincg7.prediction.bed.gz
```

#### Step 5: tissues-of-origin analysis
* generate methylation density in 1kb window at autosomes across all available methylation prediction files
```
ls *WGS.FinaleMe.mincg7.merged.cov.b37.bw | perl -ne 'chomp;$cov=$_;$m=$cov;$m=~s/cov/methy_count/;print " -bigWig $m -useMean0 0 -regionMode 0 -bigWig $cov -useMean0 0 -regionMode 0";' >> cfdna.methy_summary.cmd.txt

cat cfdna.methy_summary.cmd.txt | perl -ne 'chomp;@f=split " -useMean0 0 -regionMode 0";for($i=1,$j=0;$j<=$#f;$j+=2,$i++){$name=$f[$j];$name=~s/ -bigWig (\S+)\S+methy_count.b37.bw/$1/;print "$i\t$name\n";}' > cfdna.names_order.txt

perl -e '$cmd=`cat cfdna.methy_summary.cmd.txt`;chomp($cmd); `java -Xmx10G -cp "lib/dnaaseUtils-0.14-jar-with-dependencies.jar:lib/java-genomics-io.jar:lib/igv.jar" main.java.edu.mit.compbio.utils.AlignMultiWigInsideBed autosome_1kb_intervals.UCSC.cpgIsland_plus_shore.b37.bed output.add_value.methy.bed.gz $cmd`;'

```

* R script is available within src/R/TissueOfOriginExampleScript.R


## License

For academic research, please refer to MIT license. For commerical usage, please contact the authors.
