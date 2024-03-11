def read_bed_names(filename):
    names = []
    with open(filename) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                print("Warning: ignoring {}. No name found".format(line))
                continue
            names.append(cols[3])
    return names

########################
### FIXED PARAMETERS ###
########################
reference_fasta = config.get("reference_fasta")
if not reference_fasta:
    raise RuntimeError("No reference FASTA found. Please specify 'reference_fasta' in config file")

input_folder = config.get('input_fastq')
if not input_folder:
    raise RuntimeError("No input FASTQ files found. Please specify 'input_fastq' in config file")

target_bed = config.get('targets_bed')
if not target_bed:
    raise RuntimeError("No target BED file found. Please spcify 'targets_bed' in config file")

sample_name = config.get("sample_name", "umi_sample")

#########################
## Optional parameters ##
#########################
allowed_umi_errors = config.get("umi_errors", 3)
subset_reads = config.get("downsample_to", 0)
min_reads_per_cluster = config.get("min_reads_per_cluster", 20)
max_reads_per_cluster = config.get("max_reads_per_cluster", 60)
min_overlap = config.get("min_overlap", 0.9)
balance_strands = config.get("balance_strands", True)
mm = config.get("medaka_model", "r941_min_high_g360")
fwd_context = config.get("fwd_context", "GTATCGTGTAGAGACTGCGTAGG")
rev_context = config.get("rev_context", "AGTGATCGAGTCAGTGCGAGTG")
fwd_umi = config.get("fwd_umi", "TTTVVVVTTVVVVTTVVVVTTVVVVTTT")
rev_umi = config.get("rev_umi", "AAABBBBAABBBBAABBBBAABBBBAAA")
min_length = config.get("min_length", 40)
max_length = config.get("max_length", 60)
filter_reads = config.get("filter_reads", False)
min_read_len = config.get("min_read_len", 100)
min_mean_qual = config.get("min_mean_qual", 70)
varscan_params = config.get("varscan_params", '--variants 1 --output-vcf 1 --min-coverage 8 --min-avg-qual 0 --min-var-freq 0.01 --strand-filter 0 --p-value 1 --min-reads2 2')

########################
########################
########################

target = read_bed_names(target_bed)

balance_strands_param = "--balance_strands"
if not balance_strands:
    balance_strands_param = ""

minimap2_param = "-ax map-ont -k 13"

print("Targets: {}".format(" ".join(target)), file=sys.stderr)

########################
######### RULES ########
########################

rule reads:
    input:
        expand("{name}/targets.bed", name=sample_name),
        expand("{name}/align/{target}_final.bam.bai", name=sample_name, target=target),
        expand("{name}/stats/{target}_vsearch_cluster_stats.tsv", name=sample_name, target=target),
        expand("{name}/stats/{target}_consensus_size_vs_acc.tsv", name=sample_name, target=target)

rule variants:
    input:
        expand("{name}/variants/{target}_final.vcf.gz", name=sample_name, target=target)

rule all:
    input:
        expand("{name}/targets.bed", name=sample_name),
        expand("{name}/align/{target}_final.bam.bai", name=sample_name, target=target),
        expand("{name}/stats/{target}_vsearch_cluster_stats.tsv", name=sample_name, target=target), 
        expand("{name}/variants/{target}_final.vcf.gz", name=sample_name, target=target),
        expand("{name}/variants/{target}_consensus.vcf.gz", name=sample_name, target=target),
        expand("{name}/variants/{target}_d.vcf.gz", name=sample_name, target="1"),
        expand("{name}/stats/{target}_consensus_size_vs_acc.tsv", name=sample_name, target=target),

rule copy_bed:
    input:
        target_bed
    output:
        "{name}/targets.bed"
    shell:
        "cp {input} {output}"

rule filter_reads:
    input:
        FQ = input_folder
    params:
        min_read_len = min_read_len,
        min_mean_qual = min_mean_qual,
        filter_reads = filter_reads
    output:
        FQ = "{name}/read.filt.fastq.gz",
        STATS = "{name}/stats/reads_stats.txt"
    threads: 1
    shell:
        """
        printf 'Total reads in file pre filtering: ' 2>&1 | tee {output.STATS}
        if [[ {input.FQ} =~ \.gz$ ]]
        then
            zcat {input.FQ} | echo $((`wc -l`/4)) 2>&1 | tee {output.STATS}
        else
            cat {input.FQ} | echo $((`wc -l`/4)) 2>&1 | tee {output.STATS}
        fi
        
        if [[ {params.filter_reads} == "True" ]]
        then
            filtlong --min_length {params.min_read_len} --min_mean_q {params.min_mean_qual} {input.FQ} | gzip > {output.FQ}
            printf 'Total reads in file post filtering: ' 2>&1 | tee {output.STATS}
            zcat {output.FQ} | echo $((`wc -l`/4)) 2>&1 | tee {output.STATS}
        else
            cp {input.FQ} {output.FQ}
        fi
        """

rule map_1d:
    input:
        FQ = "{name}/read.filt.fastq.gz",
        REF = reference_fasta
    params:
        read_number = subset_reads,
        minimap2_param = minimap2_param
    output:
        BAM = "{name}/align/1_d.bam",
        BAI = "{name}/align/1_d.bam.bai"
    threads: 30
    shell:
        """
        catfishq --max_n {params.read_number} {input.FQ} | minimap2 {params.minimap2_param} -t {threads} {input.REF} - | samtools sort -@ 5 -o {output.BAM} - && samtools index -@ {threads} {output.BAM}
        rm {input.FQ} #because this is a copy now
        """

# Split reads by amplicons
rule split_reads:
    input:
        "{name}/align/1_d.bam"
    output:
        DIR = directory("{name}/fasta_filtered/"),
        STATS = "{name}/stats/umi_filter_reads_stats.txt"
    params:
        bed = target_bed,
        min_overlap = min_overlap
    shell:
        """
        mkdir -p {output.DIR}
        umi_filter_reads --min_overlap {params.min_overlap} -o {output.DIR} {params.bed} {input} 2>&1 | tee {output.STATS}
        """


# Map consensus reads after polishing
rule map_consensus:
    input:
        FA = "{name}/fasta/{target}_{type}.fasta",
        REF = reference_fasta
    params:
        minimap2_param = minimap2_param
    output:
        BAM = "{name}/align/{target}_{type}.bam",
        BAI = "{name}/align/{target}_{type}.bam.bai"
    threads: 3
    shell:
        "minimap2 {params.minimap2_param} -t {threads} {input.REF} {input.FA} | samtools sort -@ 5 -o {output.BAM} - && samtools index -@ {threads} {output.BAM}"


rule detect_umi_fasta:
    input:
        "{name}/fasta_filtered/"
    output:
        "{name}/fasta_umi/{target}_detected_umis.fasta"
    params:
        errors = allowed_umi_errors,
        fwd_context = fwd_context,
        rev_context = rev_context,
        fwd_umi = fwd_umi,
        rev_umi = rev_umi,
    shell:
        """
        umi_extract --fwd-context {params.fwd_context} --rev-context {params.rev_context} --fwd-umi {params.fwd_umi} --rev-umi {params.rev_umi} --max-error {params.errors} {input}/{wildcards.target}.fastq -o {output} --tsv {output}.tsv
        """

rule detect_umi_consensus_fasta:
    input:
        "{name}/fasta/{target}_consensus.fasta"
    output:
        "{name}/fasta_umi/{target}_detected_umis_final.fasta"
    params:
        errors = allowed_umi_errors,
        fwd_context = fwd_context,
        rev_context = rev_context,
        fwd_umi = fwd_umi,
        rev_umi = rev_umi,
    shell:
        """
        umi_extract --fwd-context {params.fwd_context} --rev-context {params.rev_context} --fwd-umi {params.fwd_umi} --rev-umi {params.rev_umi} --max-error {params.errors} {input} -o {output} --tsv {output}.tsv
        """

rule cluster:
    input: "{name}/fasta_umi/{target}_detected_umis.fasta"
    output:
        CENT = "{name}/clustering/{target}/clusters_centroid.fasta",
        CONS = "{name}/clustering/{target}/clusters_consensus.fasta",
        DIR = directory("{name}/clustering/{target}/vsearch_clusters")
    params:
        min_length = min_length,
        max_length = max_length
    threads: 10
    shell:
        "mkdir -p {wildcards.name}/clustering/{wildcards.target}/vsearch_clusters && vsearch --clusterout_id --clusters {wildcards.name}/clustering/{wildcards.target}/vsearch_clusters/test --centroids {output.CENT} --consout {output.CONS} --minseqlength {params.min_length} --maxseqlength {params.max_length} --qmask none --threads {threads} --cluster_fast {input} --clusterout_sort --gapopen 0E/5I --gapext 0E/2I --mismatch -8 --match 6 --iddef 0 --minwordmatches 0 -id 0.85"



rule cluster_consensus:
    input: "{name}/fasta_umi/{target}_detected_umis_final.fasta"
    output:
        CENT = "{name}/clustering_consensus/{target}/clusters_centroid.fasta",
        CONS = "{name}/clustering_consensus/{target}/clusters_consensus.fasta",
        DIR = directory("{name}/clustering_consensus/{target}/vsearch_clusters")
    params:
        min_length = min_length,
        max_length = max_length
    threads: 10
    shell:
        " mkdir -p {wildcards.name}/clustering_consensus/{wildcards.target}/vsearch_clusters && vsearch --clusterout_id --clusters {wildcards.name}/clustering_consensus/{wildcards.target}/vsearch_clusters/test --centroids {output.CENT} --consout {output.CONS} --minseqlength {params.min_length} --maxseqlength {params.max_length} --qmask none --threads {threads} --cluster_fast {input} --clusterout_sort --gapopen 0E/5I --gapext 0E/2I --mismatch -8 --match 6 --iddef 0 --minwordmatches 0 -id 0.85"

rule reformat_consensus_clusters:
    input:
        "{name}/clustering_consensus/{target}/clusters_consensus.fasta"
    output:
        "{name}/fasta/{target}_final.fasta"
    shell:
        "cat {input} | umi_reformat_consensus > {output}"


rule reformat_filter_clusters:
    input:
        "{name}/clustering/{target}/clusters_consensus.fasta",
        "{name}/clustering/{target}/vsearch_clusters"
    params:
        min_reads_per_cluster = min_reads_per_cluster,
        max_reads_per_cluster = max_reads_per_cluster,
        balance_strands_param = balance_strands_param
    output:
        out_dir = directory("{name}/clustering/{target}/clusters_fa/"),
        stats = "{name}/stats/{target}_vsearch_cluster_stats.tsv",
        out_file = "{name}/clustering/{target}/smolecule_clusters.fa"
    shell:
        "umi_parse_clusters --smolecule_out {output.out_file} {params.balance_strands_param} --min_reads_per_clusters {params.min_reads_per_cluster} --max_reads_per_clusters {params.max_reads_per_cluster} --stats_out {output.stats} -o {output.out_dir} {input}"


rule polish_clusters:
    input:
        I1 = "{name}/clustering/{target}/clusters_fa/",
        I2 = "{name}/clustering/{target}/smolecule_clusters.fa"
    output:
        FOLDER = directory("{name}/fasta/{target}_consensus_tmp"),
        BAM = "{name}/fasta/{target}_consensus.bam",
        F = "{name}/fasta/{target}_consensus.fasta"
    params:
        medaka_model = mm
    threads: 30
    shell:
        """
        rm -rf {output.FOLDER}
        medaka smolecule --threads {threads} --length 50 --depth 2 --model {params.medaka_model} --method spoa {output.FOLDER} {input.I2} 2> {output.BAM}_smolecule.log
        cp {output.FOLDER}/consensus.fasta {output.F}
        cp {output.FOLDER}/subreads_to_spoa.bam {output.BAM} && cp {output.FOLDER}/subreads_to_spoa.bam.bai {output.BAM}.bai
        """

rule call_variants:
    input:
        BAM = "{name}/align/{target}_{type}.bam",
        REF = reference_fasta
    output:
        "{name}/variants/{target}_{type}.vcf"
    params:
        varscan_params = varscan_params
    shell:
        "samtools mpileup -q 0 -Q 0 -B -d 10000000 -A -f {input.REF} {input.BAM} | varscan mpileup2cns {params.varscan_params} > {output}"

rule index_variants:
    input:
        "{name}/variants/{target}_{type}.vcf"
    output:
        "{name}/variants/{target}_{type}.vcf.gz"
    shell:
        "bedtools sort -header -i {input} | bgzip > {output} && tabix {output}"

rule seqkit_bam_acc_tsv:
    input:
        "{name}/align/{target}_{stage}.bam"
    output:
        "{name}/stats/{target}_{stage}_size_vs_acc.tsv"
    shell:
        """
        echo -e "Read\tCluster_size\tRef\tMapQual\tAcc\tReadLen\tRefLen\tRefAln\tRefCov\tReadAln\tReadCov\tStrand\tMeanQual\tLeftClip\tRightClip\tFlags\tIsSec\tIsSup" > {output} && seqkit bam {input} 2>&1 | sed 's/_/\t/' | tail -n +2 >> {output}
        """
