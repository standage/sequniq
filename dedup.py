#!/usr/bin/env python

# Copyright (c) 2015, Daniel S. Standage and CONTRIBUTORS
#
# HymHub is distributed under the CC BY 4.0 License. See the
# 'LICENSE' file in the HymHub code distribution or online at
# https://github.com/BrendelGroup/HymHub/blob/master/LICENSE.

import argparse
import sys
from sequniq import parse
from sequniq import records


def get_args():
    desc = 'Remove duplicate sequences (such as PCR duplicates) from input'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--fasta', action='store_true',
                        help='data are in Fasta format; default is Fastq')
    parser.add_argument('-o', '--out', default=sys.stdout, metavar='OUT',
                        help='output file; default is terminal (stdout)')
    parser.add_argument('infile', type=argparse.FileType('r'), nargs='+',
                        help='Input filename(s)')
    return parser.parse_args()


def chained_file_contents(filehandles):
    for fh in filehandles:
        for line in fh:
            yield line


def main(args=get_args()):
    assert isinstance(args.infile, list)
    fastq = not args.fasta
    for record in records.uniqseqs(chained_file_contents(args.infile),
                                   fastq=fastq):
        parse.write(record, args.out)


if __name__ == '__main__':
    main()
