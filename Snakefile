

########################
### FIXED PARAMETERS ###
########################
all_samples = config.get('input_fastqs')
if not all_samples:
    raise RuntimeError("No input FASTQ files found. Please specify 'input_fastqs' in config file")

#########################
## Optional parameters ##
#########################
allowed_umi_errors = config.get("umi_errors", 3)
subset_reads = config.get("downsample_to", 0)
min_reads_per_cluster = config.get("min_reads_per_cluster", 20)
max_reads_per_cluster = config.get("max_reads_per_cluster", 60)
balance_strands = config.get("balance_strands", True)
mm = config.get("medaka_model", "r941_min_high_g360")
########################
########################
########################

name = config.get("name")

min_length = 40
max_length = 60

balance_strands_param = "--balance_strands"
if not balance_strands:
    balance_strands_param = ""

sample_names = all_samples.keys()

########################
######### RULES ########
########################

rule all:
    input:
        expand("{name}_trimmed_umi_consensus_min{min}.fasta", name=sample_names, min=min_reads_per_cluster)


rule trim:
    input:
        lambda wildcards: all_samples[wildcards.name]
    output:
        "{name}_trimmed.fastq"
    params:
        read_number = subset_reads,
    threads: 10
    shell:
        """
        catfishq --max_n {params.read_number} {input} > {output}_tmp
        porechop -t {threads} --discard_middle -i {output}_tmp -o {output}
        rm -f {output}_tmp
        """

rule detect_umi_fasta:
    input:
        "{name}_trimmed.fastq"
    output:
        "tmp/{name}_trimmed_detected_umis.fasta"
    params:
        errors = allowed_umi_errors,
    shell:
        """
        umi_extract --max-error {params.errors} {input} -o {output} --tsv {output}.tsv
        """

rule detect_umi_consensus_fasta:
    input:
        "tmp/{name}_trimmed_umi_pre-consensus_min{min}.fasta"
    output:
        "tmp/{name}_trimmed_umi_pre-consensus_detected_umi_min{min}.fasta"
    params:
        errors = allowed_umi_errors,
    shell:
        """
        umi_extract --max-error {params.errors} {input} -o {output} --tsv {output}.tsv
        """

rule cluster:
    input: "tmp/{name}_trimmed_detected_umis.fasta"
    output:
        CENT = "tmp/{name}_clustering/clusters_centroid.fasta",
        CONS = "tmp/{name}_clustering/clusters_consensus.fasta",
        DIR = directory("tmp/{name}_clustering/vsearch_clusters")
    params:
        min_length = min_length,
        max_length = max_length
    threads: 10
    shell:
        "mkdir -p tmp/{wildcards.name}_clustering/vsearch_clusters && vsearch --clusterout_id --clusters tmp/{wildcards.name}_clustering/vsearch_clusters/test --centroids {output.CENT} --consout {output.CONS} --minseqlength {params.min_length} --maxseqlength {params.max_length} --qmask none --threads {threads} --cluster_fast {input} --clusterout_sort --gapopen 0E/5I --gapext 0E/2I --mismatch -8 --match 6 --iddef 0 --minwordmatches 0 --qmask none -id 0.85"


rule cluster_consensus:
    input: "tmp/{name}_trimmed_umi_pre-consensus_detected_umi_min{min}.fasta"
    output:
        CENT = "tmp/{name}_clustering_pre-consensus_min{min}/clusters_centroid.fasta",
        CONS = "tmp/{name}_clustering_pre-consensus_min{min}/clusters_consensus.fasta",
        DIR = directory("tmp/{name}_clustering_pre-consensus_min{min}/vsearch_clusters")
    params:
        min_length = min_length,
        max_length = max_length
    threads: 10
    shell:
        """
        if [ -s {input} ]
        then
            mkdir -p tmp/{wildcards.name}_clustering_pre-consensus_min{wildcards.min}/vsearch_clusters && vsearch --clusterout_id --clusters tmp/{wildcards.name}_clustering_pre-consensus_min{wildcards.min}/vsearch_clusters/test --centroids {output.CENT} --consout {output.CONS} --minseqlength {params.min_length} --maxseqlength {params.max_length} --qmask none --threads {threads} --cluster_fast {input} --clusterout_sort --gapopen 0E/5I --gapext 0E/2I --mismatch -8 --match 6 --iddef 0 --minwordmatches 0 --qmask none -id 0.85
        else
            touch {output.CENT} {output.CONS} 
            mkdir -p tmp/{wildcards.name}_clustering_pre-consensus_min{wildcards.min}/vsearch_clusters
        fi
        """


rule reformat_consensus_clusters:
    input:
        "tmp/{name}_clustering_pre-consensus_min{min}/clusters_consensus.fasta"
    output:
        "{name}_trimmed_umi_consensus_min{min}.fasta"
    shell:
        "cat {input} | umi_reformat_consensus > {output}"


rule reformat_filter_clusters:
    input:
        "tmp/{name}_clustering/clusters_consensus.fasta",
        "tmp/{name}_clustering/vsearch_clusters"
    params:
        min_reads_per_cluster = min_reads_per_cluster,
        max_reads_per_cluster = max_reads_per_cluster,
        balance_strands_param = balance_strands_param
    output:
        stats = "{name}_trimmed_umi_cluster_stats_min{min}.tsv",
        out_file = "tmp/{name}_clustering/smolecule_clusters_min{min}.fa"
    shell:
        "umi_parse_clusters -o {output.out_file} {params.balance_strands_param} --min_reads_per_clusters {params.min_reads_per_cluster} --max_reads_per_clusters {params.max_reads_per_cluster} --stats_out {output.stats} {input}"


rule polish_clusters:
    input:
        "tmp/{name}_clustering/smolecule_clusters_min{min}.fa"
    output:
        "tmp/{name}_trimmed_umi_pre-consensus_min{min}.fasta"
    params:
        medaka_model = mm,
        FOLDER = "tmp/{name}_consensus_tmp_min{min}"
    threads: 1
    shell:
        """
        rm -rf {params.FOLDER}
        if [ -s {input} ]
        then
            medaka smolecule --threads {threads} --length 50 --depth 1 --model {params.medaka_model} --method spoa {input} {params.FOLDER}
            cp {params.FOLDER}/consensus.fasta {output}
        else
            touch {output}
        fi
        """


