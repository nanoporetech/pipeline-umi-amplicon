import argparse
import logging
import sys

import pysam


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

    args = parser.parse_args(argv)

    return args


def parse_stdin(args):
    cluster_filename = "/dev/stdout"
    consensus_filename = "/dev/stdin"

    with open(cluster_filename, "w") as out:
        with pysam.FastxFile(consensus_filename) as fh:
            for entry in fh:
                cols = entry.name.split(";")
                read_id = cols[0]
                read_seq = cols[6].split("=")[1]
                print(">{}".format(read_id), file=out)
                print(read_seq, file=out)


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

    parse_stdin(args)


if __name__ == "__main__":
    main()
