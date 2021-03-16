![.](ONT_logo.png  "Oxford Nanopore Technologies")
******************
# pipeline-umi-amplicon

# This is a modified pipeline to work with a custom cDNA protocol. This pipeline will not work optimally with data produced according to the official UMI protocol. Please use the latest release or master branch on GitHub.

# How to run the cDNA version of the pipeline

```bash
$ conda env create -n umi-cluster ./environment.yml
$ conda activae umi-cluster
$ pip install lib/ # Only required when run for the first time
# Set correct paths in config file and run pipeline
$ snakemake --snakefile Snakefile --configfile config_high_acc_r941.yml -d results/ all --cores 30 -pr
```


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
