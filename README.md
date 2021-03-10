![.](ONT_logo.png  "Oxford Nanopore Technologies")
******************
# pipeline-umi-amplicon

### Overview

`pipeline-umi-amplicon` is a pipeline for generating high accuracy single
molecule reads using unique molecular identifiers (UMIs) from amplicon data.
The pipeline accepts FASTQ-format sequence files as input and outputs both
aligned reads and QC stats.

### Features

The pipeline performs the following steps:
- Reads are mapped to reference genome using minimap2
- Separate into amplicons
- Extract UMI sequences for all reads
- Cluster UMI sequences per amplicon using vsearch and compute high accuracy consensus reads
- Align high accuracy conesnsus reads and perform simple variant calling (optional)

******************
# Getting Started

### Requirements
The following software packages must be installed prior to running:

-  [miniconda3](https://conda.io/miniconda.html) - please refer to installation [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Installation
After installing miniconda3, install the pipeline as follows:
```bash
# Get pipeline
$ git clone https://github.com/nanoporetech/pipeline-umi-amplicon.git 
# Change to directory
$ cd pipeline-umi-amplicon
# Create conda environment with all dependencies
$ conda env create -f environment.yml
# Activate environment
$ conda activate pipeline-umi-amplicon
# Install python packages provided by pipeline-umi-amplicon
$ cd lib && pip install . && cd ..

# To test if the installation was successful run
$ snakemake -j 1 -pr --configfile config.yml
# Deactivate environment
$ conda deactivate
```

### Input

To run the pipeline the following input files are required:

| Input | Description |
|-------|-------------|
| Reference genome | FASTA file containing the reference genome (e.g. GRCh38 for human) |
| Nanopore reads| Folder containing FASTQ files or a single concatenated FASTQ file. Reads should be **q-score filtered**|
| Targets / Amplicons | A BED file containing the chromosome, start and end coordinate and the name of all amplicons|

# BED format
Tab separated and needs a unique name:
```
chr1    107167322       107168239       target_a_chr1_107167756_T_C
```

### Output

 The main output files created by the pipeline are:

| Output | Description |
|--------|-------------|
| Aligned reads | Aligned reads in indexed and sorted BAM format |
| Variant calls (optional) | Called variants in VCF format |

After the a pipeline analysis has completed, the aligned reads can be found at `{output_folder}/{run_name}/align/{amplicon_name}_final.bam` e.g. `example_egfr_single_read_run/align/EGFR_917_final.bam`.

### Usage:

To run the pipeline with default settings invoke snakemake as follows.

```bash
$ snakemake -j 30 reads --configfile config.yml
```

`-j` specifies how many CPU cores will be used by the pipeline. `reads` is the default target (see Targets); this will run all steps required to produce aligned high accuracy consensus reads. Please see the example config files for the required parameters.

### Targets

|Name| Description |
|--|--|
| reads | Only prodcues high accuracy consensus and align them to the reference |
| variants | Same as `reads` + calls variants using `varscan2` |

### Options

The pipeline accepts several input parameters. They can either be changed in the `config.yml` file or specified when running snakemake.

For example:
```bash
snakemake -j 30 reads --config input_fastq=/data/pass/ reference_fasta=/ref/hg19.fasta
```

##### Required parameters

These parameters have to be specified to run the pipeline.

| Parameter | Allowed | Description |
|-----------|---------|-------------|
| sample_name | String | Name of the output folder |
| input_fastq | Absolute file path | FASTQ file or folder containing FASTQ files |
| reference_fasta | Absolute file path | FASTA file containing the reference genome |
| targets_bed | Absolute file path | BED file containing amplicon coordinates and names |

##### Optional parameters

See `config.yml`

******************
# Help
## Licence and Copyright

(c) 2020 Oxford Nanopore Technologies Ltd.

This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

### References and Supporting Information
If you use this pipeline please cite:

- Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191
- Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9 [19505943]
- Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584
- For optional variant calling: Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing Genome Research DOI: 10.1101/gr.129684.111

### Research Release

Research releases are provided as technology demonstrators to provide early
access to features or stimulate Community development of tools. Support for this
software will be minimal and is only provided directly by the developers.
Feature requests, improvements, and discussions are welcome and can be
implemented by forking and pull requests. However much as we would like to
rectify every issue and piece of feedback users may have, the developers may
have limited resource for support of this software. Research releases may be
unstable and subject to rapid iteration by Oxford Nanopore Technologies.
