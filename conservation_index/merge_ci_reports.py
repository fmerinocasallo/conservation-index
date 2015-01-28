#! /usr/bin/python3.3
"""Merge nucleotides and amino acids Conservation Index (CI) reports."""
import argparse
from collections import deque
from errno import ENOENT, EINVAL, EIO
from ntpath import basename
from os import makedirs
from os.path import exists, splitext
import sys

__author__ = 'Francisco Merino'
__copyright__ = 'Copyright 2013, Final Project'
__credits__ = ['Francisco Merino', 'Jorge Alvarez', 'Elvira Mayordomo']
__license__ = 'GPL'
__version__ = '3.0'
__maintainer__ = 'Francisco Merino'
__email__ = 'fmerino@unizar.es'
__status__ = 'Development'

AMINO_ACIDS_PER_BLOCK = 4
NUCLEOTIDES_PER_BLOCK = 12
BLOCKS_PER_LINE = 3

def read_arguments():
    """
    Read command line's arguments
    """
    parser = argparse.ArgumentParser(description=('Merge nucleotides and '
                                                  'amino acids Conservation '
                                                  'Index (CI) reports'))
    parser.add_argument('-nr', '--nucleotide_report', required=True,
                        help=('filename of the file storing the nucleotides '
                              'CI report'))
    parser.add_argument('-aar', '--amino_acid_report', required=True,
                        help=('filename of the file storing the amino acids '
                              'CI report'))
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='give detailed information about the process')
    return parser.parse_args()

def extract_cs(filename):
    """
    Retrieve nucleotides or amino acids Consensus Sequence and the columns which
    meet the given condition from the given Conservation Index (CI) report

    Arguments:
        - filename: CI report's filename, required (str).

    Returns:
        - first_column: Absolute position of the first column of the sequence
        - consensus: Consensus Sequence (str)
        - special_columns: Columns which meet the given condition (str)
        - info: CI info (str)
    """
    if not isinstance(filename, str):
        raise TypeError('"filename" argument should be a string')

    consensus = ''
    special_columns = ''
    info = ''
    with open(filename, 'r') as report:
        num_line = 0
        copy = False
        while True:
            line1 = report.readline()
            line2 = report.readline()
            line3 = report.readline()

            if not line1:
                break # EOF

            # Retrieve the absolute position for the first column of the seq
            if num_line == 0:
                first_column = int(line1[:line1.find(':')].replace(' ', ''))

            if not (line1.find(':') == -1) and (line1.find('>') == -1):
                consensus += line1[line1.find(':') + 2:].replace('\n', '')
                special_columns += line2.replace(' ', '').replace('\n', '')

            # We would append at the end of the merged report the stats
            # about the special columns, so we need to know when we have to
            # start copying the content of the CI's report
            if not (line2.find('There') == -1):
                copy = True

            if copy:
                info += line1 + line2 + line3

            if not line2:
                num_line += 1
            elif not line3:
                num_line += 2
            else:
                num_line += 3

        return first_column, consensus, special_columns, info

def get_nucleotides_lines(first_column, line, consensus, special_columns):
    """
    Generate the next dna lines to be writen into the merged report

    Arguments:
        - first_column: absolute position of the first column of the
                        nucleotides sequence
        - line: line to be generated, required (int)
        - consensus: nucleotides Consensus Sequence, required (str)
        - special_columns: nucleotides' columns which meet the given condition,
                           required (str)

    Return:
        - num_column_line (str)
        - consensus_line (str)
        - special_columns_line (str)
    """
    if not isinstance(first_column, int):
        raise TypeError('"first_column" argument should be an int')
    if not isinstance(line, int):
        raise TypeError('"line" argument should be an int')
    if not isinstance(consensus, str):
        raise TypeError('"consensus" argument should be a string')
    if not isinstance(special_columns, str):
        raise TypeError('"special_columns" argument should be a string')

    start = NUCLEOTIDES_PER_BLOCK * BLOCKS_PER_LINE * line
    num_lines = int(-(-len(consensus) //
                        (NUCLEOTIDES_PER_BLOCK * BLOCKS_PER_LINE))) # ceil div

    # Last line
    if (line + 1) == num_lines:
        columns = len(consensus) - start
        blocks = int(-(-columns // NUCLEOTIDES_PER_BLOCK))

        num_column_line = ''
        consensus_line = ''
        special_columns_line = ''
        for block in range(0, blocks):
            # Last block
            if (block + 1) == blocks:
                num_column_line += '{:5d}:\n'.format(first_column + start +
                                                     NUCLEOTIDES_PER_BLOCK *
                                                       block)

                consensus_line += '      '
                special_columns_line += '      '
                start_column = start + NUCLEOTIDES_PER_BLOCK * block
                for i in range(start_column, len(consensus), 3):
                    consensus_line += ' {:s}'.format(consensus[i:i + 3])
                    special_columns_line += \
                                    ' {:s}'.format(special_columns[i:i + 3])

                consensus_line += '\n'
                special_columns_line += '\n'
            else:
                num_column_line += ('{:5d}:                '
                                    ''.format(first_column + start +
                                              NUCLEOTIDES_PER_BLOCK * block))

                i = start + NUCLEOTIDES_PER_BLOCK * block
                j = i + 3
                consensus_line += ('       {:s} {:s} {:s} {:s}'
                                   ''.format(consensus[i:j],
                                             consensus[i + 3:j + 3],
                                             consensus[i + 6:j + 6],
                                             consensus[i + 9:j + 9]))

                special_columns_line += ('       {:s} {:s} {:s} {:s}'
                                         ''.format(special_columns[i:j],
                                                   special_columns[i + 3:j + 3],
                                                   special_columns[i + 6:j + 6],
                                                   special_columns[i + 9:j + 9]))



    else:
        num_column_line = ('{:5d}:                '
                           '{:5d}:                '
                           '{:5d}:'
                           '\n'.format(first_column + start,
                                       first_column + start +
                                         NUCLEOTIDES_PER_BLOCK,
                                       first_column + start +
                                         2 * NUCLEOTIDES_PER_BLOCK))

        consensus_line = ('       {:s} {:s} {:s} {:s}'
                          '       {:s} {:s} {:s} {:s}'
                          '       {:s} {:s} {:s} {:s}'
                          '\n'.format(consensus[start:start + 3],
                                      consensus[start + 3:start + 6],
                                      consensus[start + 6:start + 9],
                                      consensus[start + 9:start + 12],
                                      consensus[start + 12:start + 15],
                                      consensus[start + 15:start + 18],
                                      consensus[start + 18:start + 21],
                                      consensus[start + 21:start + 24],
                                      consensus[start + 24:start + 27],
                                      consensus[start + 27:start + 30],
                                      consensus[start + 30:start + 33],
                                      consensus[start + 33:start + 36]))

        special_columns_line = ('       {:s} {:s} {:s} {:s}'
                                '       {:s} {:s} {:s} {:s}'
                                '       {:s} {:s} {:s} {:s}'
                                '\n'.format(special_columns[start:start + 3],
                                            special_columns[start + 3:start + 6],
                                            special_columns[start + 6:start + 9],
                                            special_columns[start + 9:start + 12],
                                            special_columns[start + 12:start + 15],
                                            special_columns[start + 15:start + 18],
                                            special_columns[start + 18:start + 21],
                                            special_columns[start + 21:start + 24],
                                            special_columns[start + 24:start + 27],
                                            special_columns[start + 27:start + 30],
                                            special_columns[start + 30:start + 33],
                                            special_columns[start + 33:start + 36]))

    return num_column_line, consensus_line, special_columns_line

def get_peptide_lines(line, consensus, special_columns):
    """
    Generate the next amino acids lines to be writen into the merged report

    Arguments:
        - line: line to be generated, required (int)
        - consensus: amino acids Consensus Sequence, required (str)
        - special_columns: amino acid's columns which meet the given condition,
                           required (str)

    Return:
        - consensus_line (str)
        - special_columns_line (str)
    """
    if not isinstance(line, int):
        raise TypeError('"line" argument should be an int')
    if not isinstance(peptide_mfs, str):
        raise TypeError('"consensus" argument should be a string')
    if not isinstance(special_columns, str):
        raise TypeError('"special_columns" argument should be a string')

    start = AMINO_ACIDS_PER_BLOCK * BLOCKS_PER_LINE * line
    num_lines = int(-(-len(peptide_mfs) //
                        (AMINO_ACIDS_PER_BLOCK * BLOCKS_PER_LINE))) # ceil div

    # Last line
    if (line + 1) == num_lines:
        columns = len(consensus) - start
        blocks = int(-(-columns // AMINO_ACIDS_PER_BLOCK))

        consensus_line = ''
        special_columns_line = ''
        for block in range(0, blocks):
            # Last block
            if (block + 1) == blocks:
                consensus_line += '     '
                special_columns_line += '     '

                start_column = start + AMINO_ACIDS_PER_BLOCK * block
                for i in range(start_column, len(consensus)):
                    consensus_line += '   {:s}'.format(consensus[i])
                    special_columns_line += '   {:s}'.format(special_columns[i])

                consensus_line += '\n'
                special_columns_line += '\n'
            else:
                i = start + AMINO_ACIDS_PER_BLOCK * block
                consensus_line += ('        {:s}   {:s}   {:s}   {:s} '
                                   ''.format(consensus[i],
                                             consensus[i + 1],
                                             consensus[i + 2],
                                             consensus[i + 3]))

                special_columns_line += ('        {:s}   {:s}   {:s}   {:s} '
                                         ''.format(special_columns[i],
                                                   special_columns[i + 1],
                                                   special_columns[i + 2],
                                                   special_columns[i + 3]))

    else:
        consensus_line = ('        {:s}   {:s}   {:s}   {:s} '
                          '        {:s}   {:s}   {:s}   {:s} '
                          '        {:s}   {:s}   {:s}   {:s}'
                          '\n'.format(consensus[start],
                                      consensus[start + 1],
                                      consensus[start + 2],
                                      consensus[start + 3],
                                      consensus[start + 4],
                                      consensus[start + 5],
                                      consensus[start + 6],
                                      consensus[start + 7],
                                      consensus[start + 8],
                                      consensus[start + 9],
                                      consensus[start + 10],
                                      consensus[start + 11]))

        special_columns_line = ('        {:s}   {:s}   {:s}   {:s} '
                                '        {:s}   {:s}   {:s}   {:s} '
                                '        {:s}   {:s}   {:s}   {:s}'
                                '\n'.format(special_columns[start],
                                            special_columns[start + 1],
                                            special_columns[start + 2],
                                            special_columns[start + 3],
                                            special_columns[start + 4],
                                            special_columns[start + 5],
                                            special_columns[start + 6],
                                            special_columns[start + 7],
                                            special_columns[start + 8],
                                            special_columns[start + 9],
                                            special_columns[start + 10],
                                            special_columns[start + 11]))

    return consensus_line, special_columns_line

def main():
    """
    Merge nucleotides and amino acids Conservation Index (CI) reports into one
    single report, for clarity purposes.
    """
    # Deal with command line's args
    args = read_arguments()

    # Extract the Consensus Sequence from the given nucleotides CI report
    if args.verbose:
        print('Extracting CS from the nucleotides CI report...', end='')
        sys.stdout.flush()

    first_column, \
    nucleotide_consensus, \
    nucleotide_special_columns, \
    nucleotide_info = extract_cs(args.nucleotide_report)

    # Append as many 'A' as needed in order to get a nucleotides' sequence with
    # a length multiple of 3
    if (len(nucleotide_consensus) % 3) != 0:
        nucleotide_consensus += 'a' * (3 - len(nucleotide_consensus) % 3)

    if args.verbose:
        print('done')

    # Extract the Consensus Sequence from the given amino acids CI report
    if args.verbose:
        print('Extracting CS from the amino acids CI report...', end='')
        sys.stdout.flush()

    first_column, \
    amino_acid_consensus, \
    amino_acid_special_columns, \
    amino_acid_info = extract_cs(args.amino_acid_report)

    if args.verbose:
        print('done')

    # Write the merged report
    input_filename, input_extension = splitext(args.nucleotide_report)
    #!FIXME Change the way we organize the info
    # output_filename = (paths.PATH_REPORTS_HMTDNA_MERGED_CI +
    #                    basename(input_filename) + '_merged.txt')
    output_filename = (basename(input_filename) + '_merged.txt')

    # if not exists(paths.PATH_REPORTS_HMTDNA_MERGED_CI):
    #     makedirs(paths.PATH_REPORTS_HMTDNA_MERGED_CI)

    if args.verbose:
        print('Writing merge CI report...', end='')
        sys.stdout.flush()

    with open(output_filename, 'w') as merged_report:
        line = 0
        num_lines = int(-(-len(nucleotide_consensus) //
                            (NUCLEOTIDES_PER_BLOCK * BLOCKS_PER_LINE)))

        while line < num_lines:
            num_column_line, \
            nucleotide_consensus_line, \
            nucleotide_special_columns_line = \
                get_nucleotide_lines(first_column, line, nucleotide_consensus,
                                     nucleotide_special_columns)

            merged_report.write(num_column_line)
            merged_report.write(nucleotide_consensus_line)
            merged_report.write(nucleotide_special_columns_line)

            merged_report.write('\n')

            amino_acid_consensus_line, \
            amino_acid_special_columns_line = \
                get_amino_acid_lines(line, amino_acid_consensus,
                                     amino_acid_special_columns)

            merged_report.write(amino_acid_consensus_line)
            merged_report.write(amino_acid_special_columns_line)

            merged_report.write('\n\n')

            line += 1

        merged_report.writelines(nucleotide_info)
        merged_report.write('\n')

        merged_report.writelines(amino_acid_info)
        merged_report.write('\n')

    if args.verbose:
        print('done')


if __name__ == '__main__':
    main()
