#! /usr/bin/python3.3
"""Generate a report based on the analysis of the given alignment."""
# Support for python2.7
from __future__ import print_function

import argparse
from errno import ENOENT, EIO
from math import ceil
from os.path import basename, splitext
import sys
from threading import Barrier
from multiprocessing import cpu_count

from Bio import AlignIO

from _Conservation_Index import CI_Thread, Report


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
                        choices=['entropy'],
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

def analyze(seqs_type, seqs_filename, freq_method, ci_method):
    """
    Analyze the set of sequences using every available cpu. Each column
    will be process by one of them.
    """
    # Initialize the stats dictionary which is going to store the CI
    # distribution for each column
    
    num_cpus = cpu_count()
    align = AlignIO.read(seqs_filename, 'fasta')
    num_rows = len(align)
    num_columns = align.get_alignment_length()
    freqs = {'-': [0.0] * num_columns, 
             'a': [0.0] * num_columns, 'g': [0.0] * num_columns,
             'c': [0.0] * num_columns, 't': [0.0] * num_columns}
    cis = [0.0 for k in range(0, num_columns)]
    align_weights = [0.0 for j in range(0, num_rows)]

    start_row = 0
    start_column = 0
    rows_section_size = ceil(num_rows / num_cpus)
    columns_section_size = ceil(num_columns / num_cpus)
    threads = []
    barrier = Barrier(num_cpus)
    for cpu in range(0, num_cpus):
        column_section = (start_column, start_column + columns_section_size)
        row_section = (start_row, start_row + rows_section_size)
        threads.append(CI_Thread(barrier, seqs_type, align, column_section,
                                 freq_method, freqs, ci_method, cis,
                                 row_section, align_weights))

        threads[cpu].start()

        start_row += rows_section_size
        start_column += columns_section_size

    # Wait for all the threads to finish their execution
    for thread in threads:
        thread.join()

    return freqs, cis

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
    output_filename = (args.output_directory + basename(input_filename) +
                       '_' + args.frequencies_method +
                       '_' + args.conservation_method +
                       '_' + args.condition + '_' + str(args.threshold) +
                       '_' + args.report_type + '.txt')
    if args.verbose:
        print('Starting analysis of the given alignment...', end='')
        sys.stdout.flush()

    if args.conservation_method == 'entropy':
        ci_method = 'shannon entropy'

    freqs, cis = analyze(args.sequence_type, args.input_filename,
                         args.frequencies_method, ci_method)

    if args.verbose:
        print('done')

    if args.verbose:
        print('Generating the report...', end='')
        sys.stdout.flush()

    report = Report(args.sequence_type, freqs, cis, 0)

    if args.report_type == 'basic':
        with open(output_filename, 'w') as of:
            of.write(report.generate_basic(args.condition, args.threshold))
    else:
        with open(output_filename , 'w') as of:
            of.write(report.generate_detailed(args.condition, args.threshold))

    if args.verbose:
        print('done')

if __name__ == '__main__':
    main()

