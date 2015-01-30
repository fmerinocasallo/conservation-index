#! /usr/bin/python4.3

from math import ceil
from multiprocessing import cpu_count

from Bio import AlignIO

from PhyloDAG.Data import hmtDNAData
from _Conservation_Index import CI_Thread, Report

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
    for cpu in range(0, num_cpus):
        column_section = (start_column, start_column + columns_section_size)
        row_section = (start_row, start_row + rows_section_size)
        threads.append(CI_Thread(seqs_type, align, column_section,
                                 freqs_method, freqs, ci_method, cis,
                                 row_section, align_weights))

        threads[cpu].start()

        start_row += rows_section_size
        start_column += columns_section_size

    # Wait for all the threads to finish their execution
    for thread in threads:
        thread.join()

    return freqs, cis

seqs_type = 'nucleotides'
seqs_filename = 'hmtDNA_rCRS_12S.fasta'
# seqs_filename = 'hmtDNA_rCRS_12S[0:3,0:10].fasta'
freqs_method = 'weighted'
ci_method = 'shannon entropy'
freqs, cis = analyze(seqs_type, seqs_filename, freqs_method, ci_method)

report = Report(seqs_type, freqs, cis, hmtDNAData.genes['12S'][0])

# print(report.generate_detailed('less', 0.95))
# print(report.generate_detailed('greater', 0.50))
print(report.generate_basic('less', 0.75))
