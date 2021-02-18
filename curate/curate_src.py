"""
    CLI for curation code.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 10/02/2021, 17:32 PM
"""

import argparse
import os
import sys

import curate.dataset_curation as datacur

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
    
    parser.add_argument('-s', '--smiles_col',
                        help='Column name where the SMILES string is found.',
                        required=False)
    
    parser.add_argument('-r', '--remove',
                        action='store_true',
                        help='Remove problematic structures after SMILES curation',
                        required=False)

    args = parser.parse_args()
    
    if args.infile is not None:
        if not os.path.isfile(args.infile):
            sys.stderr.write('Input file {} not found\n'.format(args.infile))
            return

    if args.command == 'curate':
        if (args.outfile is None) or (args.infile is None) or (args.format is None):
            sys.stderr.write('datacur curate : input, output and output format arguments are compulsory\n')
            return

        if args.smiles_col is None:
            smiles_ = 'structure'
        else:
            smiles_ = args.smiles_col

        curating = datacur.DataCuration(data_input=args.infile)
        curating.curate_data(structure_column=smiles_, remove_problematic=args.remove)
        curating.get_output_file(outfile_name=args.outfile, outfile_type=args.format)

if __name__ == '__main__':
    main()