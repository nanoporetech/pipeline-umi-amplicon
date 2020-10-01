import argparse
import pysam
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args(argv):
    usage = "Annotate homopolymers in FASTA file"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-m", "--min", dest="MIN", help="Min read count for clusters")
    parser.add_argument("-M", "--max", dest="MAX", help="Max read count for clusters")
    parser.add_argument(
        "-n",
        "--name",
        dest="NAME",
        help="Amplicon name. Name used in the target BED file.)",
    )
    parser.add_argument(
        "-f", "--figures", dest="FIGURES", help="Prefix for output figures"
    )
    parser.add_argument(
        "FOLDERS", type=str, nargs="+", help="UMI pipeline output folder"
    )

    args = parser.parse_args(argv)

    return args


def get_readstast_bam(folder):
    bam_1d_file = "{}/align/1d.bam".format(folder)

    mapped = 0
    unmapped = 0
    with pysam.AlignmentFile(bam_1d_file, "rb") as bam_1d:
        for r in bam_1d:
            if r.is_unmapped:
                unmapped += 1
                continue
            if not r.is_supplementary and not r.is_secondary:
                mapped += 1
    sequenced = mapped + unmapped
    return sequenced, mapped, unmapped


def count_fastx_reads_with_umi(folder, amplicon_name):
    fastq_with_umi_file = "{}/fasta_umi/{}_detected_umis.fasta".format(
        folder, amplicon_name
    )
    return count_fastx_reads(fastq_with_umi_file)


def count_fastx_reads_on_target(folder, amplicon_name):
    fastq_ontarget_file = "{}/fasta_filtered/{}.fastq".format(folder, amplicon_name)
    return count_fastx_reads(fastq_ontarget_file)


def count_fastx_reads(filename):
    count = 0
    with pysam.FastxFile(filename) as fh:
        for entry in fh:
            count += 1
    return count


def get_median_acc(folder, amplicon_name):
    acc_stats_file = "{}/stats/{}_consensus_size_vs_acc.tsv".format(
        folder, amplicon_name
    )
    acc_stats = pd.read_csv(acc_stats_file, sep="\t")
    acc_stats_prim = acc_stats.query("Flags in [16, 0]")
    return acc_stats_prim["Acc"].median()


def get_reads_usage_stats(folder, amplicon_name, min_size=20, max_size=60):
    cluster_stats_file = "{}/stats/{}_vsearch_cluster_stats.tsv".format(
        folder, amplicon_name
    )
    cluster_stats = pd.read_csv(cluster_stats_file, sep="\t")

    min_fwd = min_size / 2.0
    min_rev = min_size / 2.0

    total_reads = cluster_stats["n"].sum()

    small_clusters = cluster_stats[
        (cluster_stats["n_fwd"] < min_fwd) | (cluster_stats["n_rev"] < min_rev)
    ]
    reads_in_small_clusters = small_clusters["n"].sum()

    ok_clusters = cluster_stats[
        (cluster_stats["n_fwd"] >= min_fwd) & (cluster_stats["n_rev"] >= min_rev)
    ]
    reads_in_large_clusters = ok_clusters["n"].sum()
    cluster_count = len(ok_clusters)

    excess_reads = sum(
        cluster_stats[(cluster_stats["n"] - max_size) > 0]["n"] - max_size
    )

    usable_reads = reads_in_large_clusters - excess_reads

    usable_perc = usable_reads * 100.0 / total_reads
    reads_in_small_clusters_perc = reads_in_small_clusters * 100.0 / total_reads
    excess_reads_perc = excess_reads * 100.0 / total_reads

    return (
        cluster_count,
        usable_reads,
        total_reads,
        usable_perc,
        reads_in_small_clusters_perc,
        excess_reads_perc,
    )


def get_stats(
    folder,
    amplicon_name,
    display_name,
    min_cluster_size,
    max_cluster_size,
    figure_prefix,
):
    if figure_prefix:
        cluster_read_hist(folder, amplicon_name, display_name, figure_prefix)
    sequenced, mapped, unmapped = get_readstast_bam(folder)
    on_target = count_fastx_reads_on_target(folder, amplicon_name)
    with_umi = count_fastx_reads_with_umi(folder, amplicon_name)
    median_acc = get_median_acc(folder, amplicon_name)
    (
        cluster_count,
        usable_reads,
        total_reads,
        usable_perc,
        reads_in_small_clusters_perc,
        excess_reads_perc,
    ) = get_reads_usage_stats(folder, amplicon_name, min_cluster_size, max_cluster_size)

    if total_reads != with_umi:
        print(
            "Warning: clustering numbers don't add up. {} reads missing. This can happen if vsearch excludes UMIs because of length.".format(
                abs(total_reads - with_umi)
            ),
            file=sys.stderr,
        )

    return {
        "name": display_name,
        "sequenced": int(sequenced),
        "mapped": int(mapped),
        "on_target": int(on_target),
        "with_umi": int(with_umi),
        "min_size": min_cluster_size,
        "max_size": max_cluster_size,
        "median_acc": float(median_acc),
        "cluster_count": int(cluster_count),
        "usable_reads": int(usable_reads),
        "usable_perc": round(usable_perc, 2),
        "reads_in_small_clusters_perc": round(reads_in_small_clusters_perc, 2),
        "excess_reads_perc": round(excess_reads_perc, 2),
    }


def cluster_read_hist(folder, amplicon_name, display_name, figure_prefix):
    filename = "{}/stats/{}_vsearch_cluster_stats.tsv".format(folder, amplicon_name)
    cluster_stats = pd.read_csv(filename, sep="\t")

    current_palette = sns.color_palette()

    fig = plt.figure(figsize=(15, 8))
    sns.set_style(
        "whitegrid",
        {
            "grid.linestyle": "--",
        },
    )
    sns.set_context("paper", font_scale=3)

    tmp = cluster_stats["n"].value_counts().reset_index()
    tmp["Read_count"] = tmp["n"] * tmp["index"]

    ax = sns.barplot(x="index", color=current_palette[0], y="Read_count", data=tmp)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set(
        xlabel="Cluster size",
        ylabel="Number of reads in clusters",
        title="Cluster size distribution - Sample: {}".format(display_name),
    )

    ax.axvline(20, ls="--", zorder=1, label="Min cluster size", color="red", alpha=0.8)
    ax.axvline(60, ls="--", zorder=1, label="Max cluster size", color="red", alpha=0.8)
    ax.set_xlim(0, 200)

    xmin, xmax = ax.get_xlim()
    custom_ticks = np.linspace(xmin, xmax, 5, dtype=int)
    ax.set_xticks(custom_ticks)
    ax.set_xticklabels(custom_ticks)

    sns.despine()
    plt.tight_layout()
    fig.savefig(str(figure_prefix) + "_cluster_size_distribution.pdf", dpi=300)
    fig.savefig(str(figure_prefix) + "_cluster_size_distribution.png", dpi=100)
    plt.show()


def main(argv=sys.argv[1:]):
    args = parse_args(argv=argv)

    folders = args.FOLDERS
    figures = args.FIGURES
    name = args.NAME
    c_min = int(args.MIN)
    c_max = int(args.MAX)

    df_stats = pd.DataFrame(
        columns=[
            "name",
            "sequenced",
            "mapped",
            "on_target",
            "with_umi",
            "min_size",
            "max_size",
            "median_acc",
            "cluster_count",
            "usable_reads",
            "usable_perc",
            "reads_in_small_clusters_perc",
            "excess_reads_perc",
        ]
    )
    i = 0
    for folder in folders:
        if ":" in folder:
            cols = folder.split(":")
            display_name = cols[0]
            folder = cols[1]
        else:
            display_name = name

        figures_ext = None
        if figures:
            figures_ext = "{}_{}_{}".format(figures, name, i)
        stats = get_stats(folder, name, display_name, c_min, c_max, figures_ext)
        df_stats = df_stats.append(stats, ignore_index=True)

        i += 1
    print(df_stats.to_csv(sep="\t", index=False, header=True), end="")


if __name__ == "__main__":
    main()
