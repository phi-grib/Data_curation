"""
    CLI for curation code.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 10/02/2021, 17:32 PM
"""

import argparse
import os
import sys

import curate.context as context

from curate.util import utils, config

def main():

    parser = argparse.ArgumentParser(
                                    description='Curation tool CLI for handling structure \
                                        curation and data selection from input files.\n')

    parser.add_argument('-i', '--infile',
                        help='Input file.',
                        required=False)

    parser.add_argument('-o', '--outfile',
                        help='Output file name.',
                        required=False)
    
    parser.add_argument('-f', '--format',
                        help='Format of the output file.',
                        choices=['xlsx', 'csv', 'tsv', 'sdf'],
                        required=False)

    parser.add_argument('-a', '--action',
                        action='store',
                        help='Manage action.',
                        choices=['silent'],
                        required=False)
    
    parser.add_argument('-c', '--command',
                        action='store',
                        choices=['curate', 'split', 'config'],
                        help='command type: \'curate\' or \'split\' or \'config\'.',
                        required=True)
    
    parser.add_argument('-d', '--directory',
                        help='Defines the root directory for the curation repository.',
                        required=False)

    parser.add_argument('-id', '--id_column',
                        help='Column name containing the molecular identifiers.\
                            CAS number is recommended, although other ids can be passed.',
                        required=False)

    parser.add_argument('-s', '--smiles_col',
                        help='Column name where the SMILES string is found.',
                        required=False)
    
    parser.add_argument('-sep', '--separator',
                        help='If added, takes this argument as the file separator.',
                        required=False)

    parser.add_argument('-r', '--remove',
                        action='store_true',
                        help='Remove problematic structures after SMILES curation.',
                        required=False)

    args = parser.parse_args()
    
    if args.infile is not None:
        if not os.path.isfile(args.infile):
            sys.stderr.write('Input file {} not found\n'.format(args.infile))
            return
    
    if args.command != 'config':
        utils.config_test()

    if args.command == 'curate':
        if (args.outfile is None) or (args.infile is None) or (args.format is None):
            sys.stderr.write("datacur curate : input, output and output format arguments are compulsory\n")
            return
        
        if args.id_column is None:
            id_ = 'name'
        else:
            id_ = args.id_column

        if args.smiles_col is None:
            smiles_ = 'structure'
        else:
            smiles_ = args.smiles_col

        if args.separator:
            sep = args.separator
        else:
            sep = None

        context.curation_cmd(data_input=args.infile,
                             molecule_identifier=id_,
                             structure_column=smiles_,
                             separator=sep,
                             remove_problematic=args.remove,
                             outfile_name=args.outfile,
                             outfile_type=args.format)
    
    elif args.command == 'config':
        success, results = config.configure(args.directory, (args.action == 'silent'))
        if not success:
            sys.stderr.write("{}, configuration unchanged\n".format(results))

if __name__ == '__main__':
    main()