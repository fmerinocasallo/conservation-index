#! /usr/bin/python3.4
"""Create a new alignment randomly selecting a given number of sequences."""
# Support for python2.7
from __future__ import print_function

import argparse
from errno import ENOENT, EIO
from os.path import basename, splitext
from random import choice

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


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
    parser = argparse.ArgumentParser(description=('Create a new alignment '
                                                  'randomly selecting a given '
                                                  'number of sequences.')
                                    )
    parser.add_argument('-if', '--input_filename', required=True,
                        help=('filename of the file storing the original '
                              'alignment (FASTA format)'))

    return parser.parse_args()

def get_section(align, length):
    """
    Create a new alignment based on fragments of the one given.

    Arguments:
        - align     - alignment from which we are going to extract the
                      fragment,
                      required (MultipleSeqAlignment)
        - length    - length of the fragment to be extracted from the
                      given alignment,
                      required (int)
    """
    if not isinstance(align, MultipleSeqAlignment):
        raise TypeError('"align" argument should be a MultipleSeqAlignment')
    if not isinstance(length, int):
        raise TypeError('"length" argument should be an integer')

    # We want to pick a fragment from the HVS and other genes
    return align[:, -int(length / 2):] + align[:, :int(length / 2)]

def get_random_seqs(filename, num_seqs):
    """
    Create a new alignment randomly selecting a given number of sequences
    of an existing one.

    Arguments:
        - filename              - filename of the file which stores the
                                  original alignment,
                                  required (str)
        - num_seqs              - number of sequences to be randomly selected
                                  from the original alignment,
                                  required (int)
    """
    if not isinstance(filename, str):
        raise TypeError('"filename" argument should be a string')
    if not isinstance(num_seqs, int):
        raise TypeError('"num_seqs" argument should be an integer')

    align = AlignIO.read(filename, 'fasta')
    # Beware that choice could return the same sequence more than once
    seqs = [choice(align) for _ in range(0, num_seqs)]

    return MultipleSeqAlignment(seqs)

def main():
    seqs_sets = [100, 1000, 10000]
    size_sections = [100, 1000, 10000]

    # Deal with command line's args
    args = read_arguments()

    input_filename, input_extension = splitext(args.input_filename)
    num_dataset = 1
    for num_seqs in seqs_sets:
        align = get_random_seqs(args.input_filename, num_seqs)
        for length in size_sections:
            output_filename = ('dataset_' + str(num_dataset) + '_[' +
                               '{:d} seqs'.format(num_seqs) + ' & ' +
                               '{:d} residues'.format(length) + ']' +
                               input_extension)

            new_align = get_section(align, length)
            AlignIO.write(new_align, output_filename, 'fasta')

            num_dataset += 1


if __name__ == '__main__':
    main()

