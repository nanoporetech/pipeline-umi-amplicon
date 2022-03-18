import pysam
import glob
import parasail
import pysam
import os
import sys
import argparse

from tqdm import tqdm

#               'A','C','G','T','V','B','N'
_FWD_MATRIX = ([ 0, -1, -1, -1,  0, -1,  0], # A
               [-1,  0, -1, -1,  0,  0,  0], # C
               [-1, -1,  0, -1,  0,  0,  0], # G
               [-1, -1, -1,  0, -1,  0,  0], # T
               [ 0,  0,  0, -1,  0, -1,  0], # V A,G,C
               [-1,  0,  0,  0, -1,  0,  0], # B G,C,T
               [ 0,  0,  0,  0,  0,  0,  0]) # N

_FWD_SYMBOLS = ('A', 'C', 'G', 'T', 'V', 'B', 'N')

def create_matrix(matrix, symbols):
    parsail_matrix = parasail.matrix_create(''.join(symbols), 2, -1)
    for ii, line in enumerate(matrix):
        for jj, element in enumerate(line):
            parsail_matrix[ii, jj] = element

    return parsail_matrix

m = create_matrix(_FWD_MATRIX, _FWD_SYMBOLS)


def score(query, pattern):
    result = parasail.nw_scan_8(query, pattern, 1, 1, m)

    return abs(result.score)


def parse_args(argv):
    usage = "Annotate homopolymers in FASTA file"
    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("VSEARCH_DIR",
                        type=str,
                        nargs=1,
                        help="Directory containing vsearch clusters")

    args = parser.parse_args(argv)

    return args

def main(argv=sys.argv[1:]):
    args = parse_args(argv=argv)
    max_edit_dist = 3
    max_len_diff = 0.15
    max_umi_error = 3

    total = 0
    written = 0
    skipped_edit = 0
    skipped_len = 0
    skipped_umi = 0
    for folder in args.VSEARCH_DIR:
        for f in glob.glob("{}/*".format(folder)):
            c_name = os.path.basename(f)
            with pysam.FastxFile(f) as fh:
                reads_fwd = []
                reads_rev = []
                for rec in fh:
                    total += 1
                    split = rec.name.split(';')
                    read = {
                        'name': split[0],
                        'umi': rec.sequence
                    }
                    for el in split[1:]:
                        k, v = el.split('=')
                        read[k] = v

                    if int(read['umi_fwd']) > max_umi_error or int(read['umi_rev']) > max_umi_error:
                        skipped_umi += 1
                        continue

                    if read['strand'] == '+':
                        reads_fwd.append(read)
                    else:
                        reads_rev.append(read)

                for fwd, rev in zip(reads_fwd, reads_rev):
                    s = score(fwd['umi'], rev['umi'])
                    fwd_len = len(fwd['seq'])
                    rev_len = len(rev['seq'])
                    if s > max_edit_dist:
                        skipped_edit += 1
                        continue
                    if (abs(fwd_len - rev_len) / max(fwd_len, rev_len)) > max_len_diff:
                        skipped_len += 1
                        continue
                    print(fwd['name'], rev['name'], fwd_len, rev_len, s, fwd['umi'], rev['umi'], c_name, sep='\t')
#                     print(fwd['seq'] + rev['seq'], file=fq)
                    written += 1

    print("Total: {}".format(total), file=sys.stderr)
    print("Reads written: {}".format(written), file=sys.stderr)
    print("Skipped UMI: {}".format(skipped_umi), file=sys.stderr)
    print("Skipped edit: {}".format(skipped_edit), file=sys.stderr)
    print("Skipped len: {}".format(skipped_len), file=sys.stderr)

if __name__ == '__main__':
    main()


