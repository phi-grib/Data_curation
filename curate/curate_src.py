"""
    CLI for curation code.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 10/02/2021, 17:32 PM
"""

import argparse
import os
import pathlib
import sys

def main():

    parser = argparse.ArgumentParser(
                                    description='Curation tool CLI for handling structure curation and data selection from input files.\n')

    parser.add_argument('-i', '--infile',
                        help='Input file.',
                        required=False)

    parser.add_argument('-o', '--outfile',
                        help='Output file name',
                        required=False)
    
    parser.add_argument('-f', '--format',
                        help='Format of the output file',
                        choices=['xlsx', 'csv', 'tsv', 'sdf'],
                        required=False)

    parser.add_argument('-c', '--command',
                        action='store',
                        choices=['curate', 'split'],
                        help='Action type: \'curate\' or \'split\'',
                        required=True)
    
    args = parser.parse_args()
    
    if args.infile is not None:
        if not os.path.isfile(args.infile):
            sys.stderr.write('Input file {} not found\n'.format(args.infile))
            return

    if args.command == 'curate':
        if (args.outfile is None) or (args.infile is None):
            sys.stderr.write('datacur curate : input and output file arguments are compulsory')
            return

if __name__ == '__main__':
    main()
