import argparse
import logging
import os
import sys

import pysam
from tqdm import tqdm


def parse_args(argv):
    """
    Commandline parser

    :param argv: Command line arguments
    :type argv: List
    """
    usage = "Command line interface to telemap"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-l",
        "--log",
        dest="log",
        choices=[
            "DEBUG",
            "INFO",
            "WARNING",
            "ERROR",
            "CRITICAL",
            "debug",
            "info",
            "warning",
            "error",
            "critical",
        ],
        default="INFO",
        help="Print debug information",
    )

    parser.add_argument(
        "--max_clusters",
        dest="MAX_CLUSTERS",
        type=int,
        default=0,
        help="Stop after N clusters.",
    )

    parser.add_argument(
        "-b",
        "--balance_strands",
        dest="BAL_STRANDS",
        action="store_true",
        help="Balance strands in clusters",
    )

    parser.add_argument(
        "--min_reads_per_clusters",
        dest="MIN_CLUSTER_READS",
        type=int,
        default=20,
        help="Reads per cluster. Clusters with less reads will be discarded, clusters with more will be downsampled. 50%% must be forward and 50% reverse reads",
    )

    parser.add_argument(
        "--max_reads_per_clusters",
        dest="MAX_CLUSTER_READS",
        type=int,
        default=100,
        help="Reads per cluster. Clusters with less reads will be discarded, clusters with more will be downsampled. 50%% must be forward and 50% reverse reads",
    )

    parser.add_argument(
        "-o", "--output", dest="OUTPUT", default="clusters_fa/", help="Output folder"
    )

    parser.add_argument(
        "--stats_out",
        dest="STATS_OUT",
        default="/dev/null",
        help="Output stats file. Contains cluster size, etc. for each cluster found",
    )

    parser.add_argument(
        "--smolecule_out",
        dest="SMOLECULE_OUT",
        default="/dev/null",
        help="Input file for medaka smolecule",
    )

    parser.add_argument(
        "VSEARCH_CONSENSUS", type=str, nargs=1, help="VSearch consensus FASTA"
    )

    parser.add_argument(
        "VSEARCH_FOLDER", type=str, nargs=1, help="VSearch cluster folder"
    )

    args = parser.parse_args(argv)

    return args


def polish_cluster(
    id_cluster,
    n_cluster,
    cluster_folder,
    output_folder,
    min_reads=0,
    max_reads=10000000,
    stats_out=sys.stdout,
    balance_strands=True,
    smolecule_out=None,
    cons_umi=None,
):

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    reads_found = 0
    reads_written = 0
    reads_written_fwd = 0
    reads_written_rev = 0
    reads_written_clusters = 0

    n_fwd = 0
    n_rev = 0
    reads_skipped = 0
    cluster_written = 0

    cluster_filename = os.path.join(output_folder, "cluster{}.fasta".format(id_cluster))

    reads_fwd = {}
    reads_rev = {}

    required_fwd = max_reads
    required_rev = max_reads

    with pysam.FastxFile(
        os.path.join(cluster_folder, "test{}".format(id_cluster))
    ) as fh:
        for entry in fh:
            cols = entry.name.split(";")
            if len(cols) > 1:
                if "=" in cols[1]:
                    strand = cols[1].split("=")[1]
                else:
                    strand = "+"
            else:
                strand = "+"
            read_id = cols[0]

            # TODO: compute edit distance between cons_umi and umi. Ignore read if too large
            reads_found += 1

            if strand == "+":
                if required_fwd > 0:
                    reads_fwd[read_id] = entry
                    required_fwd -= 1
                n_fwd += 1
            else:
                if required_rev > 0:
                    reads_rev[read_id] = entry
                    required_rev -= 1
                n_rev += 1

    if balance_strands:
        min_fwd = int(min_reads / 2)
        min_rev = int(min_reads / 2)
        max_reads = min(n_fwd * 2, n_rev * 2, max_reads)
        max_fwd = int(max_reads / 2)
        max_rev = int(max_reads / 2)

        n_reads = reads_found
    else:
        min_fwd = 0
        min_rev = 0

        if n_fwd > n_rev:
            max_rev = min(n_rev, int(max_reads / 2))
            max_fwd = min(max_reads - max_rev, n_fwd)
        else:
            max_fwd = min(n_fwd, int(max_reads / 2))
            max_rev = min(max_reads - max_fwd, n_rev)

        n_reads = max_fwd + max_rev

    logging.debug(
        "Cluster: {} has {}/{} fwd and {}/{} rev reads ({} skipped)".format(
            cluster_filename, n_fwd, max_fwd, n_rev, max_rev, reads_skipped
        )
    )

    if n_fwd >= min_fwd and n_rev >= min_rev and n_reads >= min_reads:

        with open(cluster_filename, "w") as out:
            for entry in list(reads_fwd.values()) + list(reads_rev.values()):
                cols = entry.name.split(";")
                strand = "+"
                read_seq = entry.sequence
                # print(cols)
                if len(cols) > 6:
                    strand = cols[1].split("=")[1]
                    read_seq = cols[6].split("=")[1]
                read_id = cols[0]

                if strand == "+":
                    if max_fwd > 0:
                        print(">{}".format(read_id), file=out)
                        print("{}".format(read_seq), file=out)
                        if smolecule_out:
                            print(
                                ">{}_{}".format(id_cluster, reads_written),
                                file=smolecule_out,
                            )
                            print("{}".format(read_seq), file=smolecule_out)
                        reads_written += 1
                        max_fwd -= 1
                        reads_written_fwd += 1
                else:
                    if max_rev > 0:
                        print(">{}".format(read_id), file=out)
                        print("{}".format(read_seq), file=out)
                        if smolecule_out:
                            print(
                                ">{}_{}".format(id_cluster, reads_written),
                                file=smolecule_out,
                            )
                            print("{}".format(read_seq), file=smolecule_out)
                        reads_written += 1
                        max_rev -= 1
                        reads_written_rev += 1

            cluster_written = 1
            reads_written_clusters += reads_found
    else:
        logging.debug("Cluster skipped")

    logging.debug(
        "Cluster: {} has {} reads written: {} fwd - {} rev".format(
            cluster_filename, reads_written, reads_written_fwd, reads_written_rev
        )
    )

    print(
        "cluster{}".format(id_cluster),
        n_fwd,
        n_rev,
        int(max_reads / 2) - 1 - max_fwd,
        int(max_reads / 2) - 1 - max_rev,
        reads_found,
        reads_skipped,
        reads_written,
        cluster_written,
        sep="\t",
        file=stats_out,
    )
    return (
        cluster_written,
        reads_found,
        reads_skipped,
        reads_written,
        reads_written_clusters,
    )


def parse_clusters(args):
    cluster_filename = args.VSEARCH_CONSENSUS[0]
    min_read_per_cluster = args.MIN_CLUSTER_READS
    max_read_per_cluster = args.MAX_CLUSTER_READS
    stats_out_filename = args.STATS_OUT
    smolecule_filname = args.SMOLECULE_OUT
    max_clusters = args.MAX_CLUSTERS

    n_clusters = 0
    n_written = 0
    reads_found = 0
    reads_skipped = 0
    reads_written = 0
    reads_written_clusters = 0

    with pysam.FastxFile(cluster_filename) as fh:
        for entry in fh:
            n_clusters += 1

    with open(stats_out_filename, "w") as stats_out, open(
        smolecule_filname, "w"
    ) as smolecule_out:
        print(
            "id_cluster",
            "n_fwd",
            "n_rev",
            "written_fwd",
            "written_rev",
            "n",
            "skipped",
            "written",
            "cluster_written",
            sep="\t",
            file=stats_out,
        )
        with tqdm(total=n_clusters) as pbar:
            with pysam.FastxFile(cluster_filename) as fh:
                for entry in fh:
                    pbar.update(1)

                    cols = entry.name.split(";")
                    n_cluster = int(cols[-2].split("=")[1])
                    id_cluster = int(cols[-1].split("=")[1])
                    cons_umi = None  # entry.sequence.replace("T", "")

                    # if n_cluster < min_read_per_cluster:
                    #     continue

                    a, b, bb, c, d = polish_cluster(
                        id_cluster,
                        n_cluster,
                        args.VSEARCH_FOLDER[0],
                        args.OUTPUT,
                        min_read_per_cluster,
                        max_read_per_cluster,
                        stats_out,
                        args.BAL_STRANDS,
                        smolecule_out,
                        cons_umi=cons_umi,
                    )
                    n_written += a
                    reads_found += b
                    reads_skipped += bb
                    reads_written += c
                    reads_written_clusters += d
                    if max_clusters and n_written > max_clusters:
                        break

    if n_written == 0 or reads_found == 0:
        raise RuntimeError(
            "ERROR - did not find a single cluster passing the min_read threashold!"
        )

    logging.info(
        "Clusters: {}% written ({})".format(
            int(n_written * 100.0 / n_clusters), n_written
        )
    )
    logging.info("Reads: {} found".format(reads_found))
    logging.info(
        "Reads: {} removed ({}%)".format(
            reads_skipped, reads_skipped * 100.0 / (reads_skipped + reads_found)
        )
    )
    logging.info("Reads: {}% written".format(int(reads_written * 100.0 / reads_found)))
    logging.info(
        "Reads: {}% in written clusters".format(
            int(reads_written_clusters * 100.0 / reads_found)
        )
    )


def main(argv=sys.argv[1:]):
    """
    Basic command line interface to telemap.

    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.log.upper())
    logging.basicConfig(level=numeric_level, format="%(message)s")

    try:
        parse_clusters(args)
    except RuntimeError as e:
        logging.error(e)


if __name__ == "__main__":
    main()
