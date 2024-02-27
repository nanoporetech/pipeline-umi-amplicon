import argparse
import logging
import os
import sys

import pysam


def rev_comp(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement.get(base, base) for base in reversed(seq))


def str2bool(v):
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


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
        "-t", "--threads", dest="THREADS", type=int, default=1, help="Number of threads."
    )

    parser.add_argument(
        "--min_overlap",
        dest="MIN_OVERLAP",
        type=float,
        default=0.9,
        help="Min overlap with target region",
    )

    parser.add_argument(
        "--include_sec",
        dest="INCL_SEC",
        action="store_true",
        help="Include secondary alignments",
    )

    parser.add_argument(
        "-o", "--output", dest="OUT", type=str, required=False, help="Output folder"
    )

    parser.add_argument("BED", type=str, nargs=1, help="BED file")

    parser.add_argument(
        "BAM", type=str, nargs="?", default="/dev/stdin", help="BAM file"
    )

    args = parser.parse_args(argv)

    return args


def parse_bed(bed_regions):
    with open(bed_regions) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                logging.warning("Ignoring BED entry: {}".format(line))
                continue

            region = {
                "chr": cols[0],
                "start": int(cols[1]),
                "end": int(cols[2]),
                "name": cols[3],
            }
            yield region


def count_reads(bam_file):
    n_total = 0
    n_unmapped = 0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            n_total += 1
            if read.is_unmapped:
                n_unmapped += 1
    return n_total, n_unmapped


def filter_reads(args):
    bed_regions = args.BED[0]
    bam_file = args.BAM
    max_clipping = 250
    min_overlap = args.MIN_OVERLAP
    incl_sec = args.INCL_SEC
    output_folder = args.OUT

    n_concatamer = 0
    n_short = 0
    n_ontarget = 0
    n_total, n_unmapped = count_reads(bam_file)

    logging.info("Reads found: {}".format(n_total))
    unmapped_perc = 0
    if n_total:
        unmapped_perc = int(100.0 * n_unmapped / n_total)

    logging.info("Reads unmapped: {} ({}%)".format(n_unmapped, unmapped_perc))

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for region in parse_bed(bed_regions):
            logging.info(region["name"])
            n_reads_region = 0
            output_fastq = os.path.join(
                output_folder, "{}.fastq".format(region["name"])
            )
            with open(output_fastq, "w") as out:
                region_length = region["end"] - region["start"]
                for read in bam.fetch(
                    contig=region["chr"], start=region["start"], stop=region["end"]
                ):
                    if read.is_secondary and not incl_sec:
                        continue
                    if read.is_supplementary:
                        continue

                    n_ontarget += 1
                    if read.query_alignment_length < (
                        read.query_length - 2 * max_clipping
                    ):
                        n_concatamer += 1
                        continue

                    if read.reference_length < (region_length * min_overlap):
                        n_short += 1
                        continue
                    n_reads_region += 1

                    read_strand = "+"
                    if read.is_reverse:
                        read_strand = "-"
                    print(
                        ">{} strand={}".format(read.query_name, read_strand), file=out
                    )
                    if read.is_reverse:
                        print(rev_comp(read.query_sequence), file=out)
                    else:
                        print(read.query_sequence, file=out)

            logging.info("Reads found: {}".format(n_reads_region))

    ontarget_perc = 0
    if n_total:
        ontarget_perc = int(100.0 * n_ontarget / n_total)

    concatermer_perc = 0
    if n_ontarget:
        concatermer_perc = int(100.0 * n_concatamer / n_ontarget)
    if n_ontarget:
        short_perc = int(100.0 * n_short / n_ontarget)
    logging.info("On target: {} ({}%)".format(n_ontarget, ontarget_perc))
    logging.info("{} concatamers - {}%".format(n_concatamer, concatermer_perc))
    logging.info("{} short - {}%".format(n_short, short_perc))


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

    filter_reads(args)


if __name__ == "__main__":
    main()
