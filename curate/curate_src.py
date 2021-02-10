"""
    CLI for curation code.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 10/02/2021, 17:32 PM
"""

import argparse

def main():

    parser = argparse.ArgumentParser(
                                    description='Curation tool CLI for handling structure curation and data selection from input files.\n')

    parser.add_argument('-i', '--infile',
                        help='Input file.',
                        required=False)

    parser.add_argument('-o', '--outfile',
                        help='Output file',
                        required=False)
    