#-------------------------------------------------------------------------------
# File :  _Protein_Translator.py
# Description :  Definition and implementation of the classes
#                'Protein_Translator' and 'PT_Thread'.
#
# Author :  F. Merino-Casallo  ( fmerino@unizar.es )
# Last version :  v1.1 ( 03/Sep/2014 )
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  03/Sep/2014
#   VERSION :  v1.1
#   AUTHOR(s) :  F. Merino-Casallo
#   MODIFICATIONS :  Speed up and simplification of the translation process
#
#   DATE :  21/Aug/2014
#   VERSION :  v1.0
#   AUTHOR(s) :  F. Merino-Casallo
#   MODIFICATIONS :  Completed 1st development iteration
#
#   DATE :  10/Aug/2014
#   VERSION :  v0.1
#   AUTHOR(s) :  F. Merino-Casallo
#
#-------------------------------------------------------------------------------
"""Translate a nucleotides sequence into an amino acids sequence."""
from itertools import product
from threading import Thread
from _thread import LockType

from Bio import Align
from Bio.Alphabet import _get_base_alphabet, ProteinAlphabet
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein, IUPACAmbiguousDNA
from Bio.SeqRecord import SeqRecord
import Bio.Seq

#FIXME Include this info in hmtDNAData
proteins = ['ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6',
            'CO1', 'CO2', 'CO3', 'ATP6', 'ATP8', 'CYB']

class Protein_Translator (object):
    """

    """

    def __init__(self):
        """
        Creates a Protein_Translator object
        """

        self._translation_table = {'AAT': 'N', 'AAC': 'N', 'AAY': 'N',
                                   'AAA': 'K', 'AAG': 'K', 'AAR': 'K',
                                   'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
                                   'ACG': 'T', 'ACN': 'T',
                                   'AGT': 'S', 'AGC': 'S', 'AGY': 'S',
                                   'AGA': 'R', 'AGG': 'R', 'AGR': 'R',
                                   'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
                                   'ATH': 'I',
                                   'ATG': 'M',
                                   'CAT': 'H', 'CAC': 'H', 'CAY': 'H',
                                   'CAA': 'Q', 'CAG': 'Q', 'CAR': 'Q',
                                   'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
                                   'CCG': 'P', 'CCN': 'P',
                                   'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
                                   'CGG': 'R', 'CGN': 'R',
                                   'CTT': 'L', 'CTC': 'L', 'CTA': 'L',
                                   'CTG': 'L', 'CTN': 'L',
                                   'GAT': 'D', 'GAC': 'D', 'GAY': 'D',
                                   'GAA': 'E', 'GAG': 'E', 'GAR': 'E',
                                   'GCT': 'A', 'GCC': 'A', 'GCA': 'A',
                                   'GCG': 'A', 'GCN': 'A',
                                   'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
                                   'GGG': 'G', 'GGN': 'G',
                                   'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
                                   'GTG': 'V', 'GTN': 'V',
                                   'TAT': 'Y', 'TAC': 'Y', 'TAY': 'Y',
                                   'TAA': '-', 'TAG': '-', 'TAR': '-',
                                   'TCT': 'S', 'TCC': 'S', 'TCA': 'S',
                                   'TCG': 'S', 'TCN': 'S',
                                   'TGT': 'C', 'TGC': 'C', 'TGY': 'C',
                                   'TGA': '-',
                                   'TGG': 'W',
                                   'TTT': 'F', 'TTC': 'F', 'TTY': 'F',
                                   'TTA': 'L', 'TTG': 'L', 'TTR': 'L',
                                   'RAY': 'B',
                                   'SAR': 'Z'}

        # We complete the translate table with those codons which do not
        # represent any known amino acid and those which do contain at least
        # one gap
        nucleotides = IUPACAmbiguousDNA.letters +'-'
        codons = frozenset([''.join(trio)
                            for trio in product(nucleotides, repeat=3)])
        included = frozenset(self._translation_table.keys())
        not_included = codons.difference(included)
        for key in not_included:
            if key.find('-') != -1:
                self._translation_table[key] = '='
            else:
                self._translation_table[key] = 'X'

    def translate_sequences(self, nucleotides_seqs, info=''):
        """
        Translate a set of nucleotides sequences into amino acids sequences
        which is then returned.

        Arguments:
            - nucleotides_seqs  - set of nucleotides sequences to be
                                  translated, required (MultipleSeqAlignment)
            - info              - additional info about the set of nucleotides
                                - sequences to be translated, optional (str)
        """
        if not isinstance(nucleotides_seqs, Align.MultipleSeqAlignment):
            raise TypeError(('"nucleotides_seqs" argument should be a '
                             'MultipleSeqAlignment'))
        if not isinstance(info, str):
            raise TypeError('"info" argument should be a string')

        amino_acids_seqs = []
        for nucleotides_seq in nucleotides_seqs:
            if isinstance(_get_base_alphabet(nucleotides_seq.seq.alphabet),
                                             ProteinAlphabet):
                raise ValueError("Proteins cannot be translated!")

            seq = str(nucleotides_seq.seq)
            # It is required that the sequence of nucleotides has a length
            # multiple of three. If it is not, add trailing A before
            # translation
            if len(seq) % 3 != 0:
                print(("Partial codon, the nucleotides sequence's length is "
                       'not a multiple of three. We add trailing A before '
                       'translation'))
                seq += 'A' * (3 - (len(seq) % 3))

            amino_acids_seq = []
            end = len(seq) - 2
            # Translate the sequence of nucleotides into amino acids
            for i in range(0, end, 3):
                try:
                    amino_acids_seq.append(self._translation_table[seq[i:i+3]])
                except KeyError:
                    amino_acids_seq.append('X')

            amino_acids_seq_name = ('{:s}, {:s}, amino acids'
                                    ''.format(nucleotides_seq.name, info))
            amino_acids_seq_desc = ('{:s}, {:s}, amino acids'
                                    ''.format(nucleotides_seq.description,
                                              info))
            amino_acids_seq = SeqRecord(Bio.Seq.Seq(''.join(amino_acids_seq), 
                                                    ExtendedIUPACProtein()),
                                        id=nucleotides_seq.id,
                                        name=amino_acids_seq_name,
                                        description=amino_acids_seq_desc
                                       )

            amino_acids_seqs.append(amino_acids_seq)

        return amino_acids_seqs

class PT_Thread(Thread):
    """
    """
    def __init__(self, nucleotides_seqs, amino_acids_seqs, lock, info=''):
        """
        Crates a PT_Thread.

        Arguments:
            - nucleotides_seqs  - set of nucleotides sequences to be
                                  translated, required (MultipleSeqAlignment)
            - amino_acids_seqs  - set of amino acids sequences,
                                  required (list)
            - lock              - lock to control the access to the shared
                                  variable, amino_acids_seqs, required (Lock)
            - info              - additional info about the set of nucleotides
                                - sequences to be translated, optional (str)
        """
        if not isinstance(nucleotides_seqs, Align.MultipleSeqAlignment):
            raise TypeError(('"nucleotides_seqs" argument should be a '
                             'MultipleSeqAlignment'))
        if not isinstance(amino_acids_seqs, list):
            raise TypeError('"amino_acids_seqs" argument should be a list')
        elif amino_acids_seqs:
            raise ValueError(('"amino_acids_seqs" argument should be an '
                              'empty list'))
        if not isinstance(lock, LockType):
            raise TypeError('"lock" argument should be a Lock')
        if not isinstance(info, str):
            raise TypeError('"info" argument should be a string')

        # We assign to each thread the set of sequences to be translated
        # by itself
        Thread.__init__(self)
        self._pt = Protein_Translator()
        self._nucleotides_seqs = nucleotides_seqs
        self._amino_acids_seqs = amino_acids_seqs
        self._lock = lock
        self._info = info

    def run(self):
        """
        Function to be called when each thread starts its execution.
        """
        # Each thread have to translate a given set of nucleotides sequences
        amino_acids_seqs = self._pt.translate_sequences(self._nucleotides_seqs,
                                                        self._info)
        self._lock.acquire()
        try:
            self._amino_acids_seqs.extend(amino_acids_seqs)
        finally:
            self._lock.release()