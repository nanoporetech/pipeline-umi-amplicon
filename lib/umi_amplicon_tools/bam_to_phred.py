import math
import pysam
import glob
import parasail
import pysam
import os
import sys
import argparse

from tqdm import tqdm



def parse_args(argv):
    usage = "Annotate homopolymers in FASTA file"
    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("BAM",
                        type=str,
                        nargs=1,
                        help="BAM file")

    parser.add_argument("BED",
                        type=str,
                        nargs=1,
                        help="BED file")

    args = parser.parse_args(argv)

    return args

def _do_hit2tel(hit):
    """
    Converts minimap2 hit object to telemetry JSON entry.

    :param name: read name
    :param qry_len: int
    :param mean_qscore: float
    :param hit: minimap2 hit object
    :return: dict
    """
    qry_len = hit.query_length
    query_name = hit.query_name

    mean_qscore = 0 # _compute_mean_qscore(hit.query_qualities)

    summary = {}
    summary['qry_len'] = qry_len
    summary['qry_name'] = query_name
    summary['read_id'] = query_name
    summary['aligner'] = "minimap2"
    summary['mean_qscore'] = mean_qscore

    if hit.is_unmapped:
        return None

    reference_name = hit.reference_name

    el_count = [0] * 10
    for el in hit.cigartuples:
        el_count[el[0]] += el[1]

    nn = hit.get_tag("nn")
    nm = hit.get_tag("NM")
    dels = el_count[2]
    ins = el_count[1]

    mismatches = nm - dels - ins
    matches = el_count[0] - mismatches

    aln_score = 0

    accuracy = float(matches) / (matches + nm) * 100
    identity = float(matches) / (matches + mismatches) * 100
    query_alignment_start = hit.query_alignment_start
    query_alignment_end = hit.query_alignment_end
    query_alignment_length = hit.query_alignment_length

    ref_alignment_start = hit.reference_start
    ref_alignment_end = hit.reference_end
    ref_alignment_length = hit.reference_length

    qry_coverage = float(query_alignment_length) / qry_len * 100

    summary['accuracy'] = round(accuracy, 2)
    summary['identity'] = round(identity, 2)
    summary['matches'] = matches #+ mismatches
    summary['mismatches'] = mismatches
    summary['deletions'] = dels
    summary['insertions'] = ins
    summary['aln_len'] = matches + nm
    summary['qry_cov'] = round(qry_coverage, 3)
    summary['ref_name'] = reference_name
    summary['ref_start'] = ref_alignment_start
    summary['ref_end'] = ref_alignment_end
    summary['qry_start'] = query_alignment_start
    summary['qry_end'] = query_alignment_end
    summary['forward'] = not hit.is_reverse
    summary['aln_score'] = aln_score

    return summary


def main(argv=sys.argv[1:]):
    args = parse_args(argv=argv)
    bam = args.BAM[0]
    bed_file = args.BED[0]

    print("read_name", "aplicon_id", "cluster_size", "forward_reads", "reverse_reads", "accuracy", "phred", "phred_low", "read_len", "aln_len", "matches", "mismatches", 'deletions', 'insertions', 'strand', sep='\t')
    with open(bed_file) as bed:
        for line in bed:
            if not line:
                continue
            cols = line.strip().split("\t")
            chrom = cols[0]
            start = int(cols[1]) + 1
            end = int(cols[2]) + 1
            name = cols[3]
            with pysam.AlignmentFile(bam, "rb") as reads:
                for read in reads.fetch(region="{}:{}-{}".format(chrom, start, end)):
                    if read.is_supplementary or read.is_secondary:
                        continue

                    results = _do_hit2tel(read)
                    #if "_" in results['qry_name']:
                    #    qname, cluster_size = results['qry_name'].split('_')
                    #else:
                    qname = results['qry_name']
                    cluster_size = 0

                    fwd = int(cluster_size) / 2
                    p = (results['mismatches'] + results['insertions'] + results['deletions']) / results['aln_len']
                    if p:
                        qscore = -10 * math.log10(p)
                        qscore_low = qscore
                    else:
                        p = (1) / results['aln_len']
                        qscore = 'inf'
                        qscore_low = -10 * math.log10(p)
    #                         qscore = int(qscore + 0.5)

                    print(qname, name, cluster_size, fwd, fwd, results['accuracy'], qscore, qscore_low, results['qry_len'], results['aln_len'], results['matches'], results['mismatches'], results['deletions'], results['insertions'], results['forward'], sep='\t')


if __name__ == '__main__':
    main()

