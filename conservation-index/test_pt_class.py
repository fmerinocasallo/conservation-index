#! /usr/bin/python3.3

from math import ceil
from multiprocessing import cpu_count
from operator import itemgetter
from threading import Lock

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

from _Protein_Translator import PT_Thread

def translate(seqs_filename):
    """
    Translate the set of sequences using every available cpu. Each subset of sequences
    will be processed by one of them.
    """
    mtDNA_seqs = AlignIO.read(seqs_filename, 'fasta')
    num_cpus = cpu_count()
    # Because we want the code to be compatible with python2.7, we need
    # to ensure (len_seqs / num_cpus) return a float value.
    num_seqs = ceil(len(mtDNA_seqs) / float(num_cpus))
    start_seq = 0
    lock = Lock()
    threads = []
    peptides = []
    for cpu in range(0, num_cpus):
        end_seq = start_seq + num_seqs
        threads.append(PT_Thread(mtDNA_seqs[start_seq:end_seq],
                                 peptides, lock))
        threads[cpu].start()

        start_seq += num_seqs

    # Wait for all the threads to finish their execution
    for thread in threads:
        thread.join()

    return MultipleSeqAlignment(peptides)

seqs_filename = 'hmtDNA_rCRS_12S.fasta'
peptides = translate(seqs_filename)

SeqIO.write(peptides, 'hmtDNA_rCRS_12S_peptides.fasta', 'fasta')
