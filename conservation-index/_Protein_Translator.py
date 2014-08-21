#-------------------------------------------------------------------------------
# File :  _Protein_Translator.py
# Description :  Definition and implementation of the classes
#                'Protein_Translator' and 'PT_Thread'.
#
# Author :  F. Merino-Casallo  ( fmerino@unizar.es )
# Last version :  v1.0 ( 21/Aug/2014 )
#-------------------------------------------------------------------------------
# Historical report :
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
"""Translate a hmtDNA sequence into a peptide (amino acids sequence)."""
from threading import Thread
from _thread import LockType

from Bio import Align
from Bio.Alphabet import SingleLetterAlphabet
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
                                   'SAR': 'Z',
                                   'NNN': 'X'}

    def translate_columns(self, seqs):
        """
        Translate a set of mtDNA sequences into peptides (amino acids
        sequences) which is then returned.

        Arguments:
            - seqs              - set of sequences to be analyzed,
                                  required (MultipleSeqAlignment)
        """
        if not isinstance(seqs, Align.MultipleSeqAlignment):
            raise TypeError('"seqs" argument should be a MultipleSeqAlignment')

        peptides = []
        for mtDNA in seqs:
            # It is required that the mtDNA sequence has a length multiple 
            # of three. If it is not, add trailing A before translation
            #FIXME According to Bio.Seq.translate() it should be trailing N
            if len(mtDNA) % 3 != 0:
            #     print(('Last incomplete codon: {:s}'
            #            ''.format(mtDNA.seq[-(len(mtDNA) % 3):])))
                mtDNA += 'A' * (3 - len(mtDNA) % 3)
            #     print('Last incomplete and regenerated codon: {:s}'
            #           ''.format(mtDNA.seq[-3:]))

            # incomplete_codons = {}
            if mtDNA.seq.find('-') == -1:
                # There is no gap in the sequence
                peptide = Bio.Seq.translate(mtDNA.seq,
                                            stop_symbol='-')
                peptide.alphabet = SingleLetterAlphabet()
                #FIXME Change SeqRecord's attributes
                peptide = SeqRecord(peptide,
                                    id=mtDNA.id,
                                    name=mtDNA.name,
                                    description=mtDNA.description)
            else:
                peptide = ''
                sequence_size = len(mtDNA)
                end = sequence_size - 2
                # Translate the hmtDNA sequence into amino acids and
                # store the position of the codons containing a gap
                for i in range(0, end, 3):
                    codon = str(mtDNA.seq[i:i+3])
                    if codon.find('-') == -1:
                        if codon in self._translation_table:
                            peptide += self._translation_table[codon]
                        else:
                            peptide += 'X'
                            # print('Unknown codon: {:s}\n'.format(codon))
                    else:
                        peptide += '='
                        # incomplete_codons[start + i + 1] = codon


                #FIXME Change SeqRecord's attributes
                peptide = SeqRecord(Bio.Seq.Seq(peptide),
                                    id=mtDNA.id,
                                    name=mtDNA.name,
                                    description=mtDNA.description)

            peptides.append(peptide)

        return peptides

class PT_Thread(Thread):
    """
    """
    def __init__(self, mtDNA_seqs, peptides, lock):
        """
        Crates a PT_Thread.

        Arguments:
            - mtDNA_seqs    - set of mtDNA sequences to be translated,
                              required (MultipleSeqAlignment)
            - peptides      - set of amino acids sequences,
                              required (list)
            - lock          - lock to control the access to the shared
                              variable, peptides, required (Lock)
        """
        if not isinstance(mtDNA_seqs, Align.MultipleSeqAlignment):
            raise TypeError(('"mtDNA_seqs" argument should be a '
                             'MultipleSeqAlignment'))
        if not isinstance(peptides, list):
            raise TypeError('"peptides" argument should be a list')
        elif peptides:
            raise ValueError('"peptides" argument should be an empty list')
        if not isinstance(lock, LockType):
            raise TypeError('"lock" argument should be a Lock')

        # We assign to each thread the set of sequences to be translated
        # by itself
        Thread.__init__(self)
        self._pt = Protein_Translator()
        self._mtDNA_seqs = mtDNA_seqs
        self._peptides = peptides
        self._lock = lock

    def run(self):
        """
        Function to be called when each thread starts its execution.
        """
        # Each thread have to translate a given set of mtDNA sequences
        peptides = self._pt.translate_columns(self._mtDNA_seqs)
        self._lock.acquire()
        try:
            self._peptides.extend(peptides)
        finally:
            self._lock.release()