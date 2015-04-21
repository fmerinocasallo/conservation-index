#! /usr/bin/python3.4
"""Generate a report based on the analysis of the given alignment."""
# Support for python2.7
from __future__ import print_function

import argparse
from errno import ENOENT, EIO
from math import ceil
from os.path import basename, splitext
import sys

from Bio import AlignIO

from _Conservation_Index import Conservation_Index, Report


__author__ = 'Francisco Merino'
__copyright__ = 'Copyright 2013, Final Project'
__credits__ = ['Francisco Merino', 'Jorge Alvarez', 'Elvira Mayordomo']
__license__ = 'GPL'
__version__ = '3.0'
__maintainer__ = 'Francisco Merino'
__email__ = 'fmerino@unizar.es'
__status__ = 'Development'

def read_arguments():
    """
    Read command line's arguments
    """
    parser = argparse.ArgumentParser(description=('Generate a report based on '
                                                  'the analysis of the given '
                                                  'alignment')
                                    )
    parser.add_argument('-if', '--input_filename', required=True,
                        help=('filename of the file storing the alignment to '
                              'be analyzed (FASTA format to be aligned)'))
    parser.add_argument('-st', '--sequence_type', required=True,
                        # choices=['dna', 'protein'],
                        choices=['dna'],
                        help='type of sequences in the alignment')
    parser.add_argument('-rt', '--report_type', required=True,
                        choices=['basic', 'detailed'],
                        help='type of report to be generated')
    parser.add_argument('-fm', '--frequencies_method', required=True,
                        choices=['unweighted', 'weighted'],
                        help='method to be used to estimated frequencies')
    parser.add_argument('-cm', '--conservation_method', required=True,
                        choices=['entropy', 'variance'],
                        help='method to be used to estimated conservation')
    parser.add_argument('-co', '--condition', required=True,
                        choices=['greater', 'less'],
                        help='condition to be met by columns')
    parser.add_argument('-th', '--threshold', required=True, type=float,
                        help='threshold associated to the given condition')
    parser.add_argument('-od', '--output_directory', default='',
                        help=('location of the report to be generated. '
                              'If it is not specified, it will be stored '
                              'in the same directory as this script'))
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='give detailed information about the process')
    return parser.parse_args()

def main():
    """
    Generate a report based on the analysis of the given alignment. The
    analysis will be done using the settings specified by the user through
    the command line.
    """

    # Deal with command line's args
    args = read_arguments()

    # Set the new text file's filename which is going to store the report
    input_filename, input_extension = splitext(args.input_filename)
    output_filename = []
    if args.output_directory:
        output_filename += args.output_directory + '/'

    output_filename += (basename(input_filename) +
                        '_' + args.frequencies_method +
                        '_' + args.conservation_method +
                        '_' + args.condition + '_' + 
                        str(args.threshold) +
                        # '{:1.2f}'.format(args.threshold) +
                        '_' + args.report_type + '.txt')
    if args.verbose:
        print('Starting analysis of the given alignment...', end='')
        sys.stdout.flush()

    if args.conservation_method == 'entropy':
        ci_method = 'shannon entropy'
    else:
        ci_method = args.conservation_method

    ci = Conservation_Index(args.sequence_type)
    freqs, cis = ci.analyze(args.input_filename, args.frequencies_method,
                            ci_method)

    if args.verbose:
        print('done')

    if args.verbose:
        print('Generating the report...', end='')
        sys.stdout.flush()

    report = Report(args.sequence_type, freqs, cis)

    if args.report_type == 'basic':
        with open(''.join(output_filename), 'w') as of:
            of.write(report.generate_basic(args.condition, args.threshold))
    else:
        with open(''.join(output_filename) , 'w') as of:
            of.write(report.generate_detailed(args.condition, args.threshold))

    if args.verbose:
        print('done')

if __name__ == '__main__':
    main()

