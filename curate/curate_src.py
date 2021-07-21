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

    parser.add_argument('-e', '--endpoint',
                        help='Endpoint curation name.',
                        required=False)
    
    parser.add_argument('-f', '--format',
                        help='Format of the output file.',
                        choices=['xlsx', 'csv', 'tsv', 'sdf','json'],
                        required=False)

    parser.add_argument('-a', '--action',
                        action='store',
                        help='Manage action.',
                        choices=['silent','new','list','remove','chembl','export', 'download'],
                        required=False)
    
    parser.add_argument('-c', '--command',
                        action='store',
                        choices=['curate', 'split', 'config', 'manage'],
                        help='command type: \'curate\' or \'split\' or \'config\' or \'manage\'.',
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
    
    parser.add_argument('-m', '--metadata',
                        help='Selects the metadata columns to be processed. Optional.',
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
        if args.action == 'chembl':
            pass
        elif not os.path.isfile(args.infile):
            sys.stderr.write('Input file {} not found\n'.format(args.infile))
            return
    
    if args.command != 'config':
        utils.config_test()

    if args.command == 'curate':

        if args.action == 'chembl':
            from curate.chem import chembl_extraction
            input_file = chembl_extraction.get_dataframe_from_target(args.infile)
            args.id_column = 'molecule_chembl_id'
            args.smiles_col = 'canonical_smiles'
        else:
            input_file = args.infile

        if (input_file is None) or (args.endpoint is None):
            sys.stderr.write("datacur curate : input and endpoint arguments are compulsory\n")
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

        if args.metadata:
            meta_ = args.metadata
        else:
            meta_ = None
        
        command = {'data_input':args.infile,
                    'molecule_identifier':id_,
                    'structure_column':smiles_,
                    'metadata':meta_,
                    'separator':sep,
                    'remove_problematic':args.remove,
                    'endpoint':args.endpoint}

        context.curation_cmd(command)
    
    elif args.command == 'config':
        success, results = config.configure(args.directory, (args.action == 'silent'))
        if not success:
            sys.stderr.write("{}, configuration unchanged\n".format(results))
    
    elif args.command == 'manage':
        success, results = context.manage_cmd(args)
        if not success:
            sys.stderr.write(results)

if __name__ == '__main__':
    main()