"""
UMI amplicion tools setup script

Copyright (c) 2020 by Oxford Nanopore Technologies Ltd.
"""
import os
from setuptools import setup

__version__ = '1.0.0'

setup(
    name='umi-amplicon-tools',
    version=__version__,
    author='Nanopore Technologies Ltd.',
    description='Toolset to work with ONT amplicon sequencing using UMIs',
    zip_safe=False,
    install_requires=[
        'tqdm',
        'pysam',
        'numpy',
        'pandas',
        'seaborn',
        'edlib',
        'biopython'
    ],
    packages=['umi_amplicon_tools'],
    package_data={
        'umi_amplicon_tools': ['data/*'],
    },
    entry_points={
        "console_scripts": [
            'umi_filter_reads = umi_amplicon_tools.filter_reads:main',
            'umi_extract = umi_amplicon_tools.extract_umis:main',
            'umi_reformat_consensus = umi_amplicon_tools.reformat_consensus:main',
            'umi_parse_clusters = umi_amplicon_tools.parse_clusters:main',
            'umi_bam_to_phred = umi_amplicon_tools.bam_to_phred:main',
            'umi_stats = umi_amplicon_tools.umi_stats:main'
        ]
    },
)

