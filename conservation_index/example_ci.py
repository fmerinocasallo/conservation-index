#! /usr/bin/python3.3

from math import ceil
from multiprocessing import cpu_count

from Bio import AlignIO

from PhyloDAG.Data import hmtDNAData
from _Conservation_Index import CI_Thread, Report

def analyze(seqs_type, seqs_filename):
    """
    Analyze the set of sequences using every available cpu. Each column
    will be process by one of them.
    """
    # Initialize the stats dictionary which is going to store the CI
    # distribution for each column
    
    seqs = AlignIO.read(seqs_filename, 'fasta')
    len_seqs = seqs.get_alignment_length()
    stats = {'-': [0.0] * len_seqs, 'A': [0.0] * len_seqs,
            'G': [0.0] * len_seqs, 'C': [0.0] * len_seqs,
            'T': [0.0] * len_seqs}
    num_cpus = cpu_count()
    # Because we want the code to be compatible with python2.7, we need
    # to ensure (len_seqs / num_cpus) return a float value
    num_columns = ceil(len_seqs / float(num_cpus))
    start_column = 0
    threads = []
    for cpu in range(0, num_cpus):
        section = (start_column, start_column + num_columns)
        threads.append(CI_Thread(seqs_type, seqs, section, stats))
        threads[cpu].start()

        start_column += num_columns

    # Wait for all the threads to finish their execution
    for thread in threads:
        thread.join()

    return stats

seqs_type = 'mtdna'
seqs_filename = 'hmtDNA_rCRS_12S.fasta'
stats = analyze(seqs_type, seqs_filename)

report = Report(seqs_type, stats, hmtDNAData.genes['12S'][0])

#print(report.generate_detailed('greater', 0.99))
print(report.generate_basic('less', 0.75))
