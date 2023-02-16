# MethHMM

MethHMM is a Java program to predict DNA methylation in deep and low-coverage cell-free DNA WGS data without other training data.

## Citation

Cite our paper:

Liu Y# et al. (2023) MethHMM: Predicting DNA methylation by the fragmentation patterns of plasma cell-free DNA. Preprint


## Installation

### System requirements:

- Java (tested in openjdk-18.0.1.1 at x64 platform, CentOS Linux 8)
- Apache Maven (tested in v3.8.6, only if you need to compile the source code from scratch)
- Perl (tested in v5.26.3), bedGraphToBigWig from UCSC tools (tested in v4), bedtools (tested in v2.29.2). (only if you need to convert predicted methylation level to big wig files)
- R (tested in v4.2.1) and quadprog package. (only if you need to perform tissues-of-origin analysis)

### Quick installation:

    mvn compile assembly:single

Or use the precompiled .jar file from Releases and other .jars from lib/ directory without installation

### Other required data:

- methylation prior file in standard big wig format (can use wgbs_buffyCoat_jensen2015GB.methy.hg19.bw file directly from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7647046.svg)](https://doi.org/10.5281/zenodo.7647046))
- bed files to mask the Dark regions in the genome (can use wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed files directly from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7647046.svg)](https://doi.org/10.5281/zenodo.7647046))
- Chromosome sizes: can be obtained from the FASTA file of the reference genome: `samtools faidx input.fa`. See this [example](https://github.com/epifluidlab/cragr/blob/3d419a49/inst/extdata/human_g1k_v37.chrom.sizes).

## Getting started

### Input
- coordinate-sorted and indexed bam file

### Start analysis

The analysis consists of several stages:

1. Calculate the feature file for each CpG in each DNA fragment from an indexed bam file.
2. Train HMM model from the feature file (only utilize fragments with >=7 CpGs).
3. Predict the methylation status at each CpG in each DNA fragment.
4. Aggregate the methylation status across DNA fragments overlapped with the same CpG sites in the reference genome (0-100% scale).
5. Perform tissues-of-origin analysis by predicted DNA methylation level in cfDNA and methylome from reference panel.

#### Stage 1
TBD

#### Stage 2
TBD

#### Stage 3
TBD

#### Stage 4
TBD

#### Stage 5
TBD

## License

For academic research, please refer to MIT license. For commerical usage, please contact the authors.