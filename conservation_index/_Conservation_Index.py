#-------------------------------------------------------------------------------
# File :  _Conservation_Index.py
# Description :  Definition and implementation of the classes
#                'Conservation_Index', 'PT_Thread' and 'Report'.
#
# Author :  F. Merino-Casallo  ( fmerino@unizar.es )
# Last version :  v4.2 ( 11/Mar/2015 )
#-------------------------------------------------------------------------------
# Historical report :
# 
#   DATE :  11/Mar/2015
#   VERSION :  v4.2
#   AUTHOR(s) :  F. Merino-Casallo
#   MODIFICATIONS :  Changed the way conservation is calculated. Now we
#                    are using independent processes instead of threads. This
#                    let us avoid the GIL which was a complete pain in the ass.
#
#   DATE :  02/Mar/2015
#   VERSION :  v4.1
#   AUTHOR(s) :  F. Merino-Casallo
#   MODIFICATIONS :  Changed the way frequencies' are calculated. Now we
#                    are using independent processes instead of threads. This
#                    let us avoid the GIL which was a complete pain in the ass.
#
#   DATE :  26/Feb/2015
#   VERSION :  v4.0
#   AUTHOR(s) :  F. Merino-Casallo
#   MODIFICATIONS :  Changed the way sequences' weights are calculated. Now we
#                    are using independent processes instead of threads. This
#                    let us avoid the GIL which was a complete pain in the ass.
#
#   DATE :  19/Feb/2015
#   VERSION :  v3.1
#   AUTHOR(s) :  F. Merino-Casallo
#   MODIFICATIONS :  Added new method to estimate conservation index, now you
#                    can use a variance-based method.
#
#   DATE :  30/Jan/2015
#   VERSION :  v3.0
#   AUTHOR(s) :  F. Merino-Casallo
#   MODIFICATIONS :  Added new method to estimate conservation index, now you
#                    can use an entropy-based method. Furthermore, the Report
#                    class has been updated using a new design.
#
#   DATE :  12/Nov/2014
#   VERSION :  v2.1
#   AUTHOR(s) :  F. Merino-Casallo
#   MODIFICATIONS :  Added performance improvements and bug fixes on the
#                    Conservation_Index class
#
#   DATE :  22/Oct/2014
#   VERSION :  v2.0
#   AUTHOR(s) :  F. Merino-Casallo
#   MODIFICATIONS :  Added new method to estimate elements' frequencies, now
#                    you can use weights associated to sequences.
#
#   DATE :  03/Sep/2014
#   VERSION :  v1.1
#   AUTHOR(s) :  F. Merino-Casallo
#   MODIFICATIONS :  Minor performance improvements on the Conservation_Index
#                    and Report classes.
#
#   DATE :  21/Aug/2014
#   VERSION :  v1.0
#   AUTHOR(s) :  F. Merino-Casallo
#   MODIFICATIONS :  Added  Report class
#
#   DATE :  11/Jul/2014
#   VERSION :  v0.1
#   AUTHOR(s) :  F. Merino-Casallo
#
#-------------------------------------------------------------------------------
"""Estimate the CI for each column of the given set of sequences."""
from collections import Counter, defaultdict
from functools import reduce
from math import ceil, log, sqrt
from multiprocessing import cpu_count, Pool
from operator import add, itemgetter

from Bio import Align, AlignIO, SeqIO

def _check_section(filename, section_type, section):
    """
    Check if the given section is well formated. Return its boundaries.

    Arguments:
        - filename              - filename of the file which stores the
                                  alignment to be analyzed,
                                  required (str)
        - section_type          - type of section to be checked,
                                  required (str)
        - section               - section of each sequence to be analyzed,
                                  required (tuple)
    """
    if not isinstance(filename, str):
        raise TypeError('"filename" argument should be a string')
    if not isinstance(section_type, str):
        raise TypeError('"section_type" argument should be a string')
    if not isinstance(section, tuple):
        raise TypeError('"section" argument should be a tuple')
    try:
        if section[0] < 0:
            raise ValueError(('"section" argument has an invalid value. '
                              "It should start with a positive integer'"))
        if section[1] < 0:
            raise ValueError(('"section" argument has an invalid value. '
                              "It should finish with a positive integer'"))
        if section[0] > section[1]:
            print(('"section" argument include an inverse range '
                   '(x, y with x > y). I am sure you know what you are doing '
                   'so I will keep it that way'))

        if section_type == 'rows':
            num_rows = sum(1 for _ in SeqIO.parse(filename, 'fasta'))
            
            return section[0], num_rows if num_rows < section[1] \
                                        else section[1]
        elif section_type == 'columns':
            num_columns = len(next(SeqIO.parse(filename, 'fasta')).seq)

            return section[0], num_columns if num_columns < section[1] \
                                          else section[1]
        else:
            raise ValueError(('"section_type" argument has an invalid value. '
                              "It should be 'rows' or 'columns'"))
    except TypeError:
       raise TypeError('"section" argument should be a tuple of integers')

class Conservation_Index(object):
    """

    """
    def __init__(self, seqs_type):
        """
        Creates a Conservation_Index.

        Arguments:
            - seqs_type         - type of sequences,
                                  required (str)

        The type of sequences should be one of those in PhyloDag.Data
        """
        if not isinstance(seqs_type, str):
            raise TypeError('"seqs_type" argument should be a string')

        self._seqs_type = seqs_type.lower()

        # _weights is a dictionary which stores the weights associated with
        # each element of the selected alphabet. Because we still does not
        # know how many sequences we are dealing with, we set all weights to
        # zero.
        if self._seqs_type == 'dna':
            self._residues = ['-', 'a', 'c', 'g', 't']
            self._elem_weights = {'-': [['-'], 0],
                                  'a': [['a'], 0],
                                  'c': [['c'], 0],
                                  'g': [['g'], 0],
                                  't': [['t'], 0],
                                  'r': [['a', 'g'], 0],
                                  'y': [['c', 't'], 0],
                                  'w': [['a', 't'], 0],
                                  's': [['c', 'g'], 0],
                                  'm': [['a', 'c'], 0],
                                  'k': [['g', 't'], 0],
                                  'h': [['a', 'c', 't'], 0],
                                  'b': [['c', 'g', 't'], 0],
                                  'v': [['a', 'c', 'g'], 0],
                                  'd': [['a', 'g', 't'], 0],
                                  'n': [['a', 'c', 'g', 't'], 0]}
        elif self._seqs_type == 'protein':
            self._residues = ['=', 'g', 'a', 'v', 'l', 'i', 'p', 'f', 'y', 'c',
                              'm', 'h', 'k', 'r', 'w', 's', 't', 'd', 'e', 'n',
                              'q', '-', 'x']
            self._elem_weights = {'=': [['='], 0],
                                  'g': [['g'], 0],
                                  'a': [['a'], 0],
                                  'v': [['v'], 0],
                                  'l': [['l'], 0],
                                  'i': [['i'], 0],
                                  'p': [['p'], 0],
                                  'f': [['f'], 0],
                                  'y': [['y'], 0],
                                  'c': [['c'], 0],
                                  'm': [['m'], 0],
                                  'h': [['h'], 0],
                                  'k': [['k'], 0],
                                  'r': [['r'], 0],
                                  'w': [['w'], 0],
                                  's': [['s'], 0],
                                  't': [['t'], 0],
                                  'd': [['d'], 0],
                                  'e': [['e'], 0],
                                  'n': [['n'], 0],
                                  'q': [['q'], 0],
                                  'b': [['d', 'n'], 0],
                                  'z': [['e', 'q'], 0],
                                  '-': [['-'], 0],
                                  'x': [['x'], 0]}
        else:
            raise ValueError(('"seqs_type" argument is not a valid sequences '
                              "type. It should be 'nucleotides' or 'amino acids'"))

    def set_elem_weights(self, filename, freq_method, weights=None):
        """
        Set the correct weight associated with each element of the selected
        alpabet.

        Arguments:
            - filename          - filename of the file which stores the
                                  alignment to be analyzed,
                                  required (str)
            - freq_method       - frequencies' calculation method,
                                  required (str)
            - weights           - sequences' weights,
                                  optional (list of float)
        """
        def unweighted_frequencies(self, num_seqs):
            """
            Set the correct element's weight for the unweighted frequencies'
            method.

            Arguments:
                - num_seqs      - number of sequences to be analyzed,
                                  required (int)
            """
            for elem, weights in self._elem_weights.items():
                num_elems = len(weights[0])
                self._elem_weights[elem][1] = 1 / (num_elems * num_seqs)

        def weighted_frequencies(self, weights):
            """
            Set the correct element's weight for the weighted frequencies'
            method.

            Arguments:
                - weights       - sequences' weights,
                                  required (list of float)
            """
            weights_sum = sum(seq_weight for seq_weight in weights)

            for elem, weights in self._elem_weights.items():
                num_elems = len(weights[0])
                self._elem_weights[elem][1] = 1 / (num_elems * weights_sum)

        if weights is None:
            weights = []

        if not isinstance(filename, str):
            raise TypeError('"filename" argument should be a string.')
        if not isinstance(freq_method, str):
            raise TypeError('"freq_method" argument should be a string')
        if not isinstance(weights, list):
            raise TypeError('"weights" argument should be a list')

        num_seqs = sum(1 for _ in SeqIO.parse(filename, 'fasta'))
        if freq_method == 'unweighted':
            unweighted_frequencies(self, num_seqs)
        elif freq_method == 'weighted':
            weighted_frequencies(self, weights)
        else:
            raise ValueError(('"freq_method" argument has an invalid value. '
                              "It should be 'unweighted' or 'weighted'"))

    def _calculate_section_weights(self, filename, section):
        """
        Weigh the given section (usually more than a single column) of a set
        of aligned sequences. Our aim is to compensate for over-
        representation among multiple aligned sequences.

        Arguments:
            - filename      - filename of the file which stores the
                                alignment to be analyzed,
                                required (str)
            - section       - section of each sequence to be analyzed,
                                required (tuple)
        """
        if not isinstance(filename, str):
            raise TypeError('"filename" argument should be a string.')
        if not isinstance(section, tuple):
            raise TypeError('"section" argument should be a tuple')

        start_column, end_column = _check_section(filename, 'columns', section)
        num_columns = end_column - start_column
        stats = [defaultdict(int) for k in range(0, num_columns)]
        # calculate stats (# diff. elems, reps each elem) for each column
        weights = []
        for record in SeqIO.parse(filename, 'fasta'):
            for num_column, residue in \
                    enumerate(record.seq[start_column:end_column].lower()):
                stats[num_column][residue] += 1
            weights.append(0.0)

        # calculate the weights associated to each column of each sequence
        for num_seq, record in enumerate(SeqIO.parse(filename, 'fasta')):
            for num_column, residue in \
                    enumerate(record.seq[start_column:end_column].lower()):
                weights[num_seq] += 1 / (len(stats[num_column]) *
                                         stats[num_column][residue])

        return weights

    def calculate_weights(self, filename):
        """
        Weigh a set of aligned sequences. Our aim is to compensate for over-
        representation among multiple aligned sequences.

        Arguments:
            - filename          - filename of the file which stores the
                                  alignment to be analyzed,
                                  required (str)
        """
        if not isinstance(filename, str):
            raise TypeError('"filename" argument should be a string.')

        num_cores = cpu_count()
        num_columns = len(next(SeqIO.parse(filename, 'fasta')).seq)
        section_size = ceil(num_columns / num_cores)

        num_section = 0
        with Pool(processes=num_cores) as pool:
            results = []
            for start_column in range(0, num_columns, section_size):
                section = (start_column, start_column + section_size)
                results.append(pool.apply_async(self._calculate_section_weights,
                                                args=(filename, section)))

            # retrieve the results of each process
            sections_weights = [result.get() for result in results]

        weights = zip(*[section_weights for section_weights in \
                                            sections_weights])

        return [sum(seq_weights) for seq_weights in weights]

    def _calculate_section_frequencies(self, num_section, filename, section,
                                       freq_method, weights=None,
                                       calc_overall_freqs=False):
        """
        Calculate frequencies of the given section (usually more than a single
        column) of a set of sequences.

        Arguments:
            - num_section   - number of section to be analyzed. It is needed to
                                be able to sort the returned results by each of
                                the processes,
                                required (int)
            - filename      - filename of the file which stores the alignment
                                to be analyzed,
                                required (str)
            - section       - section of each sequence to be analyzed,
                                required (tuple)
            - freq_method   - frequencies' estimation method,
                                required (str)
            - weights       - sequences' weights,
                                optional (list of float)
            - calc_overall_freqs
                            - if True, we'll also calculate the overall
                                frequencies for each residue,
                                optional (bool)

        Return a tuple. If we do not want to calculate the overalls frequencies
        for each residue (calc_overall_freqs should be False in this case), the
        returned tuple will not contain them. Otherwise, the tuple will include
        them.
        """
        def unweighted_frequencies(self, filename, section,
                                   calc_overall_freqs=False):
            """
            Estimate frequencies using the unweighted frequencies' method.

            Arguments:
                - filename  - filename of the file which stores the
                                alignment to be analyzed,
                                required (str)
                - section   - section of each sequence to be analyzed,
                                required (tuple)
                - calc_overall_freqs
                            - if True, we'll also calculate the overall
                                frequencies for each residue,
                                optional (bool)

            Return a tuple. If we do not want to calculate the overalls
            frequencies for each residue (calc_overall_freqs should be False in
            this case), the returned tuple will not contain them. Otherwise,
            the tuple will include them.
            """
            def calculate_without_overall_freqs(self, filename, section):
                """
                Estimate only frequencies for each column.

                Arguments:
                    - filename
                            - filename of the file which stores the alignment
                                to be analyzed,
                                required (str)
                    - section
                            - section of each sequence to be analyzed,
                                required (tuple)

                Return a tuple containing just the frequencies.
                """
                start_column, end_column = _check_section(filename, 'columns',
                                                          section)
                section_length = end_column - start_column
                freqs = {}
                for residue in self._residues:
                    freqs[residue] = [0.0] * section_length 

                for record in SeqIO.parse(filename, 'fasta'):
                    for j, residue in \
                    enumerate(record.seq[start_column:end_column].lower()):
                        for elem in self._elem_weights[residue][0]:
                            try:
                                freqs[elem][j] += self._elem_weights[residue][1]
                            except TypeError:
                                raise TypeError(('"freqs" argument should be a'
                                                 "dict of lists of 'int'"))

                return (freqs, )

            def calculate_with_overall_freqs(self, filename, section):
                """
                Estimate not only frequencies for each column, but also the
                overall frequencies for each residue.

                Arguments:
                    - filename
                            - filename of the file which stores the alignment
                                to be analyzed,
                                required (str)
                    - section
                            - section of each sequence to be analyzed,
                                required (tuple)

                Return a tuple containing both frequencies and the overalls
                frequencies for each residue.
                """
                start_column, end_column = _check_section(filename, 'columns',
                                                          section)
                section_length = end_column - start_column
                freqs = {}
                overall_freqs = defaultdict(float)
                for residue in self._residues:
                    freqs[residue] = [0.0] * section_length 

                total_columns = len(next(SeqIO.parse(filename, 'fasta')).seq)
                for record in SeqIO.parse(filename, 'fasta'):
                    for j, residue in \
                    enumerate(record.seq[start_column:end_column].lower()):
                        for elem in self._elem_weights[residue][0]:
                            try:
                                freqs[elem][j] += self._elem_weights[residue][1]
                                overall_freqs[elem] += \
                                            (self._elem_weights[residue][1] /
                                             total_columns)
                            except TypeError:
                                raise TypeError(('"freqs" argument should be a'
                                                 "dict of lists of 'int'"))

                return (freqs, overall_freqs)

            if calc_overall_freqs:
                return calculate_with_overall_freqs(self, filename, section)
            else:
                return calculate_without_overall_freqs(self, filename, section)

        def weighted_frequencies(self, filename, section, weights,
                                 calc_overall_freqs=False):
            """
            Estimate frequencies using the weighted frequencies' method.

            Arguments:
                - filename  - filename of the file which stores the alignment
                                to be analyzed,
                                required (str)
                - section   - section of each sequence to be analyzed,
                                required (tuple)
                - weights   - sequences' weights,
                                optional (list of float)
                - calc_overall_freqs
                            - if True, we'll also calculate the overall
                                frequencies for each residue,
                                optional (bool)

            Return a tuple. If we do not want to calculate the overalls
            frequencies for each residue (calc_overall_freqs should be False
            in this case), the returned tuple will not contain them. Otherwise,
            the tuple will include them.
            """
            def calculate_without_overall_freqs(self, filename, section,
                                                weights):
                """
                Estimate only frequencies for each column.

                Arguments:
                    - filename
                            - filename of the file which stores the alignment
                                to be analyzed,
                                required (str)
                    - section
                            - section of each sequence to be analyzed,
                                required (tuple)
                    - align_weights
                            - sequences' weights,
                                required (list of float)

                Return a tuple containing just the frequencies.
                """
                start_column, end_column = _check_section(filename, 'columns',
                                                          section)
                section_length = end_column - start_column
                freqs = {}
                for residue in self._residues:
                    freqs[residue] = [0.0] * section_length 

                for i, record in enumerate(SeqIO.parse(filename, 'fasta')):
                    for j, residue in \
                    enumerate(record.seq[start_column:end_column].lower()):
                        for elem in self._elem_weights[residue][0]:
                            try:
                                freqs[elem][j] += (weights[i] *
                                                   self._elem_weights[residue][1])
                            except TypeError:
                                raise TypeError(('"freqs" argument should be a'
                                                 "dict of lists of 'int'"))
                            except IndexError:
                                raise IndexError(('"align_weights" argument '
                                                  'should have the same size '
                                                  "as the alignment's length"))

                return (freqs, )

            def calculate_with_overall_freqs(self, filename, section,
                                             weights):
                """
                Estimate not only frequencies for each column, but also the
                overall frequencies for each residue.

                Arguments:
                    - filename
                            - filename of the file which stores the alignment
                                to be analyzed,
                                required (str)
                    - section
                            - section of each sequence to be analyzed,
                                required (tuple)
                    - align_weights
                            - sequences' weights,
                                required (list of float)

                Return a tuple containing both frequencies and the overalls
                frequencies for each residue.
                """
                align = AlignIO.read(filename,'fasta')
                start_column, end_column = _check_section(filename, 'columns',
                                                          section)
                section_length = end_column - start_column
                freqs = {}
                overall_freqs = defaultdict(float)
                for residue in self._residues:
                    freqs[residue] = [0.0] * section_length 

                total_columns = len(next(SeqIO.parse(filename, 'fasta')).seq)
                for i, record in enumerate(SeqIO.parse(filename, 'fasta')):
                    for j, residue in \
                    enumerate(record.seq[start_column:end_column].lower()):
                        for elem in self._elem_weights[residue][0]:
                            try:
                                freqs[elem][j] += self._elem_weights[residue][1]
                                overall_freqs[elem] += \
                                            (self._elem_weights[residue][1] /
                                             total_columns)
                            except TypeError:
                                raise TypeError(('"freqs" argument should be a'
                                                 "dict of lists of 'int'"))
                            except IndexError:
                                raise IndexError(('"align_weights" argument '
                                                  'should have the same size '
                                                  "as the alignment's length"))

                return (freqs, overall_freqs)

            if calc_overall_freqs:
                return calculate_with_overall_freqs(self, filename, section,
                                                    weights)
            else:
                return calculate_without_overall_freqs(self, filename, section,
                                                       weights)

        if weights is None:
            weights = []

        if not isinstance(num_section, int):
            raise TypeError('"num_section" argument should be a int')
        if not isinstance(filename, str):
            raise TypeError('"filename" argument should be a string')
        if not isinstance(section, tuple):
            raise TypeError('"section" argument should be a tuple')
        if not isinstance(freq_method, str):
            raise TypeError('"freq_method" argument should be a string')
        if not isinstance(weights, list):
            raise TypeError('"weights" argument should be a list')
        if not isinstance(calc_overall_freqs, bool):
            raise TypeError('"calc_overall_freqs" argument should be a bool')

        if freq_method == 'unweighted':
            return ((num_section, ) +
                    unweighted_frequencies(self, filename, section,
                                           calc_overall_freqs))
        elif freq_method == 'weighted':
            return ((num_section, ) +
                    weighted_frequencies(self, filename, section,
                                         weights, calc_overall_freqs))
        else:
            raise ValueError(('"freq_method" argument has an invalid value. '
                                "It should be 'unweighted' or 'weighted'"))


    def calculate_frequencies(self, filename, freq_method, align_weights=None,
                              calc_overall_freqs=False):
        """
        Calculate frequencies of a set of aligned sequences.

        Arguments:
            - filename          - filename of the file which stores the
                                  alignment to be analyzed,
                                  required (str)
            - freq_method       - frequencies' estimation method,
                                  required (str)
            - align_weights     - sequences' weights,
                                  optional (list of float)
            - calc_overall_freqs
                                - if True, we'll also calculate the overall
                                  frequencies for each residue,
                                  optional (bool)

        Return a tuple. If we do not want to calculate the overalls
        frequencies for each residue (calc_overall_freqs should be False
        in this case), the returned tuple will not contain them. Otherwise,
        the tuple will include them.
        """
        if align_weights is None:
            align_weights = []

        if not isinstance(filename, str):
            raise TypeError('"filename" argument should be a string')
        if not isinstance(freq_method, str):
            raise TypeError('"freq_method" argument should be a string')
        if not isinstance(align_weights, list):
            raise TypeError('"align_weights" argument should be a list')
        if not isinstance(calc_overall_freqs, bool):
            raise TypeError('"calc_overall_freqs" argument should be a bool')

        num_columns = len(next(SeqIO.parse(filename, 'fasta')).seq)
        num_cores = cpu_count()
        section_size = ceil(num_columns / num_cores)

        num_section = 0
        with Pool(processes=num_cores) as pool:
            results = []
            num_section = 0
            for start_column in range(0, num_columns, section_size):
                end_column = start_column + section_size
                if num_columns < end_column:
                    end_column = num_columns
                section = (start_column, end_column)

                results.append(pool.apply_async(
                                        self._calculate_section_frequencies,
                                        args=(num_section, filename,
                                              section, freq_method,
                                              align_weights,
                                              calc_overall_freqs)))

                num_section += 1

            # retrieve the results of each process
            sections_freqs = [result.get() for result in results]
            # because we want to access the results in order, we have to sort
            # them previously by num_section.
            sections_freqs.sort()

        # recreate the original data structure which stores the residues'
        # frequencies for each column
        columns_freqs = {}
        for residue in self._residues:
            columns_freqs[residue] = []

        for section_freqs in sections_freqs:
            for residue, freqs in section_freqs[1].items():
                columns_freqs[residue].extend(freqs)

        if calc_overall_freqs:
            overall_freqs = dict(reduce(add, map(Counter, [section_freqs[2] \
                                        for section_freqs in sections_freqs])))

            return (columns_freqs, overall_freqs)
        else:
            return (columns_freqs, )

    def _calculate_section_conservation(self, num_section, filename, section,
                                        freqs, ci_method,
                                        overall_freqs=None):
        """
        Calculate conservation index of the given section (usually more than
        a single column) of a set of sequences.

        Arguments:
            - num_section   - number of section to be analyzed. It is
                                needed to be able to sort the returned
                                results by each of the processes,
                                required (int)
            - filename      - filename of the file which stores the
                                alignment to be analyzed,
                                required (str)
            - section       - section of each sequence to be analyzed,
                                required (tuple)
            - freqs         - estimated residues' frequencies for each
                                column, required (dict of lists of int)
            - ci_method     - conservation index's calculation method,
                                required (str)
            - overall_freqs - estimated overall frequencies for each
                                residue, optional (dict of floats)

        Return a tuple containing just a list with the conservation of each
        column.
        """
        def shannon_entropy(filename, section, freqs):
            """
            Estimate conservation index using the Shannon Entropy method.

            Arguments:
                - filename  - filename of the file which stores the
                                alignment to be analyzed,
                                required (str)
                - section   - section of each sequence to be
                                analyzed,
                                required (tuple)
                - freqs     - estimated residues' frequencies for each
                                column,
                                required (dict of lists of int)

            Return a tuple containing just a list with the conservation of
            each column.
            """
            num_rows = sum(1 for _ in SeqIO.parse(filename, 'fasta'))
            lambda_t = log(min(num_rows, len(freqs)))

            start_column, end_column = _check_section(filename, 'columns',
                                                      section)

            section_length = end_column - start_column
            cis = [0.0] * section_length

            for i, freqs_column in \
                                enumerate(freqs[start_column:end_column]):
                for residue, freq in freqs_column:
                    if freq != 0.0:
                        cis[i] += freq * log(freq)

                # We want scale the entropy to range [0, 1]
                cis[i] = 1 - (-1 * cis[i] / lambda_t)

            return (cis, )

        def variance(filename, section, freqs, overall_freqs):
            """
            Estimate conservation index using a variance-based measure.

            Arguments:
                - filename  - filename of the file which stores the
                                alignment to be analyzed,
                                required (str)
                - section   - section of each sequence to be
                                analyzed,
                                required (tuple)
                - freqs     - estimated residues' frequencies for each
                                column,
                                required (dict of lists of int)
                - overall_freqs
                            - estimated overall frequencies for each
                                residue,
                                required (dict of floats)

            Return a tuple containing just a list with the conservation of
            each column.
            """
            num_columns = len(next(SeqIO.parse(filename, 'fasta')).seq)
            max_value = sqrt(2 - 3 / num_columns + 1 / (num_columns ** 2))

            start_column, end_column = _check_section(filename, 'columns',
                                                      section)

            section_length = end_column - start_column
            cis = [0.0] * section_length 

            for i, freqs_column in \
                                enumerate(freqs[start_column:end_column]):
                for residue, freq in freqs_column:
                    if freq != 0.0:
                        cis[i] += (freq - overall_freqs[residue]) ** 2

                cis[i] = sqrt(cis[i]) / max_value 

            return (cis, )

        columns_freqs = [list(zip(freqs.keys(), values)) for values in \
                                                        zip(*freqs.values())]

        if overall_freqs is None:
            overall_freqs = {}

        if ci_method == 'shannon entropy':
            return ((num_section, ) + shannon_entropy(filename, section,
                                                      columns_freqs))
        elif ci_method == 'variance':
            return ((num_section, ) + variance(filename, section, columns_freqs,
                                               overall_freqs))
        else:
            raise ValueError(('"ci_method" argument has an invalid value. '
                                "It should be 'shannon entropy' or "
                                "'variance'"))

    def calculate_conservation(self, filename, freqs, ci_method,
                               overall_freqs=None):
        """
        Calculate conservation index of a set of sequences.

        Arguments:
            - filename          - filename of the file which stores the
                                  alignment to be analyzed,
                                  required (str)
            - freqs             - estimated residues' frequencies for each
                                  column, required (dict of lists of int)
            - ci_method         - conservation index's calculation method,
                                  required (str)
            - overall_freqs     - estimated overall frequencies for each
                                  residue, optional (dict of floats)

        Return a tuple containing just a list with the conservation of each
        column.
        """
        if overall_freqs is None:
            overall_freqs = {}

        if not isinstance(filename, str):
            raise TypeError('"filename" argument should be a string')
        if not isinstance(ci_method, str):
            raise TypeError('"ci_method" argument should be a string')
        if not isinstance(freqs, dict):
            raise TypeError('"freqs" argument should be a dictionary')
        elif not all(isinstance(residue, list) for residue in freqs.values()):
            raise TypeError('"freqs" argument should be a dictionary of lists')
        if not isinstance(overall_freqs, dict):
            raise TypeError('"overall_freqs" argument should be a dictionary')

        num_columns = len(next(SeqIO.parse(filename, 'fasta')).seq)
        num_cores = cpu_count()
        section_size = ceil(num_columns / num_cores)

        with Pool(processes=num_cores) as pool:
            results = []
            num_section = 0
            for start_column in range(0, num_columns, section_size):
                end_column = start_column + section_size
                if num_columns < end_column:
                    end_column = num_columns
                section = (start_column, end_column)

                results.append(pool.apply_async(
                                        self._calculate_section_conservation,
                                        args=(num_section, filename, section,
                                              freqs, ci_method, overall_freqs)))

                num_section += 1

            # retrieve the results of each process
            sections_conservations = [result.get() for result in results]
            # because we want to access the results in order, we have to sort
            # them previously by num_section.
            sections_conservations.sort()

        columns_conservation = []
        for section_conservations in sections_conservations:
            columns_conservation.extend(section_conservations[1])

        return (columns_conservation, )

    def analyze(self, filename, freq_method, ci_method):
        """
        Analyze the set of sequences using every available cpu. Each process
        will have assigned a different set of columns (also known as section)
        to analyze on their own.

        Arguments:
            - filename          - filename of the file which stores the
                                  alignment to be analyzed,
                                  required (str)
            - freq_method       - frequencies' estimation method,
                                  required (str)
            - ci_method         - conservation index's calculation method,
                                  required (str)

        Return a tuple containing both frequencies and conservation.
        """
        if not isinstance(filename, str):
            raise TypeError('"filename" argument should be a string')
        if not isinstance(freq_method, str):
            raise TypeError('"freq_method" argument should be a string')
        if not isinstance(ci_method, str):
            raise TypeError('"ci_method" argument should be a string')

        if freq_method == 'weighted':
            weights = self.calculate_weights(filename)
        elif freq_method == 'unweighted':
            weights = None
        else:
            raise ValueError(('"freq_method" argument has an invalid value. '
                                "It should be 'unweighted' or 'weighted'"))

        self.set_elem_weights(filename, freq_method, weights)

        if ci_method == 'variance':
            freqs, overall_freqs = self.calculate_frequencies(filename,
                                                              freq_method,
                                                              weights, True)
            cis, = self.calculate_conservation(filename, freqs, ci_method,
                                               overall_freqs)
        elif ci_method == 'shannon entropy':
            freqs, = self.calculate_frequencies(filename, freq_method, weights)
            cis, = self.calculate_conservation(filename, freqs, ci_method)
        else:
            raise ValueError(('"ci_method" argument has an invalid value. '
                              "It should be 'shannon entropy' or 'variance'"))

        return (freqs, cis)


class Report:
    def __init__(self, seqs_type, freqs, cis, start_column=0, ITEMS_PER_LINE=50):
        """
        Creates a Report object.

        Arguments:
            - seqs_type         - type of sequences,
                                  required (str)
            - freqs             - estimated residues' frequencies for each
                                  column,
                                  required (dict of lists of int)
            - cis               - estimated conservation indices for each
                                  column,
                                  required (list of float)
            - start_column      - absolute position of the 1st column,
                                  optional (int)
            - ITEMS_PER_LINE    - number of columns per line,
                                  optional (int)

        The type of sequences should be one of those in PhyloDag.Data and the
        start_column should be 0 if we are dealing with a global set of
        sequences or the absolute position of the 1st column of the gene.
        """
        if not isinstance(seqs_type, str):
            raise TypeError('"seqs_type" argument should be a string')
        if seqs_type.lower() not in ['dna', 'protein']:
            raise ValueError(('"seqs_type" argument has an invalid value. '
                              "It should be 'nucleotides' or 'amino acids'"))
        if not isinstance(freqs, dict):
            raise TypeError('"freqs" argument should be a dictionary')
        elif not all(isinstance(residue, list) for residue in freqs.values()):
            raise TypeError('"freqs" argument should be a dictionary of lists')
        if not isinstance(cis, list):
            raise TypeError('"cis" argument should be a list')
        if not isinstance(start_column, int):
            raise TypeError('"start_column" argument should be an int')
        if not isinstance(ITEMS_PER_LINE, int):
            raise TypeError('"ITEMS_PER_LINE" argument should be an int')

        self._seqs_type = seqs_type.lower()
        self._freqs = freqs
        self._cis = cis
        self._start_column = start_column
        self._ITEMS_PER_LINE = ITEMS_PER_LINE

    @staticmethod
    def _check_condition(report_type, condition, threshold, stats_column,
                         wanted_columns, condition_met=None):
        """
        Check if the current column meet the given condition. If so, we
        store its frequencies distribution to include it afterwards in
        the summary of all the columns we are looking for in this
        analysis.

        Arguments:
            - report_type       - type of report,
                                  required (str)
            - condition         - condition to be met,
                                  required (str)
            - threshold         - threshold to consider a column of interest,
                                  required (float)
            - stats_column      - list containing the results of the analysis
                                  previously done,
                                  required(list)
            - wanted_columns    - summary of the columns meeting the
                                  condition given,
                                  required (list)
            - condition_met     - set of symbols indicating which of the
                                  columns of the consensus sequence meet the
                                  given condition,
                                  optional (list)
        """
        def summarize_column(stats_column, wanted_columns):
            """
            Summarize the info previously retrieve through the analysis of
            the given alignment for the current column.

            Arguments:
                - wanted_columns
                                - summary of the columns meeting the
                                  condition given,
                                  required (list)
            """
            sorted_freqs = sorted(stats_column[2], key=itemgetter(0))
            #FIXME Should I print elem.upper()?
            distribution = ("'{:s}': {:8.4f}%"
                            ''.format(elem[0], elem[1] * 100) \
                                      for elem in sorted_freqs
                                      if elem[1] * 100 > 0.0)
            #FIXME Remember to use the mod operator when printing
            # the column
            wanted_columns.append('> {:5d}: {:5.4f} ({:s})\n'
                                  ''.format(stats_column[0], stats_column[1],
                                            ', '.join(distribution)))

        if condition_met is None:
            condition_met = []

        if not isinstance(report_type, str):
            raise TypeError('"report_type" argument should be a string')
        if report_type not in ['basic', 'detailed']:
            raise ValueError(('"report_type" argument has an invalid'
                              "value. It should be 'basic' or 'detailed'"))
        if not isinstance(condition, str):
            raise TypeError('"condition" argument should be a string')
        if not isinstance(threshold, float):
            raise TypeError('"threshold" argument should be a float')
        if not isinstance(stats_column, tuple):
            raise TypeError('"stats_column" argument should be a list')
        if not isinstance(wanted_columns, list):
            raise TypeError('"wanted_columns" argument should be a list')
        if not isinstance(condition_met, list):
            raise TypeError('"condition_met" argument should be a list')

        if condition == 'greater':
            if stats_column[1] > threshold:
                if report_type == 'basic':
                    condition_met.append('+')
                summarize_column(stats_column, wanted_columns)
            else:
                if report_type == 'basic':
                    condition_met.append('-')
        elif condition == 'less':
            if stats_column[1] < threshold:
                if report_type == 'basic':
                    condition_met.append('+')
                summarize_column(stats_column, wanted_columns)
            else:
                if report_type == 'basic':
                    condition_met.append('-')
        else:
            raise ValueError(('"condition" argument has an invalid '
                              "value. It should be 'greater' or "
                              "'less'"))

    def _summarize_wanted_columns(self, summary, condition, threshold,
                                  wanted_columns):
        """
        Summarize which columns of the previously analyzed alignment meet the
        given condition. It includes not only the number of columns found, but
        also the conservation index of each one as well as their frequencies
        distribution.

        Arguments:
            - summary           - content of the report been generated,
                                  required (list)
            - condition         - condition to be met, required (str)
            - threshold         - threshold to consider a column of
                                  interest, required (float)
            - wanted_columns    - summary of the columns meeting the
                                  condition given, required (list)
        """
        if not isinstance(summary, list):
            raise TypeError('"summary" argument should be a list')
        if not isinstance(condition, str):
            raise TypeError('"condition" argument should be a string')
        if not isinstance(threshold, float):
            raise TypeError('"threshold" argument should be a float')
        if not isinstance(wanted_columns, list):
            raise TypeError('"wanted_columns" argument should be a list')

        num_wanted_columns = len(wanted_columns)
        if num_wanted_columns > 0:
            if num_wanted_columns > 1:
                if condition == 'greater':
                    summary.append('\nThere are {:d} columns with a high '
                                   'degree of conservation (>{:6.2f}%) in '
                                   'the {:s} sequences:\n\n'
                                   ''.format(num_wanted_columns, threshold * 100,
                                             self._seqs_type))
                elif condition == 'less':
                    summary.append('\nThere are {:d} columns with a high '
                                   'degree of variation (<{:6.2f}%) in the '
                                   '{:s} sequences:\n\n'
                                   ''.format(num_wanted_columns, threshold * 100,
                                             self._seqs_type))
                else:
                    raise ValueError(('"condition" argument has an '
                                      'invalid value. It should be '
                                      "'greater' or 'less'"))
            else:
                if condition == 'greater':
                    summary.append('\nThere is 1 column with a high degree '
                                   'of conservation (>{:6.2f}%) in the {:s} '
                                   'sequences:\n\n'
                                   ''.format(threshold * 100, self._seqs_type))
                elif condition == 'less':
                    summary.append('\nThere is 1 column with a high degree '
                                   'of variation (<{:6.2f}%) in the {:s} '
                                   'sequences:\n\n'
                                   ''.format(threshold * 100, self._seqs_type))
                else:
                    raise ValueError(('"condition" argument has an '
                                      'invalid value. It should be '
                                      "'greater' or 'less'"))

            summary.append(''.join(wanted_columns))
        else:
            if condition == 'greater':
                summary.append('\nThere are not any columns with a high degree '
                               'of conservation (>{:6.2f}%) in the {:s} '
                               'sequences.'
                               ''.format(threshold * 100, self._seqs_type))
            elif condition == 'less':
                summary.append('\nThere are not any columns with a high degree '
                               'of variation (<{:6.2f}%) in the {:s} sequences.'
                               ''.format(threshold * 100, self._seqs_type))
            else:
                raise ValueError(('"condition" argument has an invalid '
                                  "value. It should be 'greater' or "
                                  "'less'"))

    def generate_basic(self, condition, threshold):
        """
        Generate a basic report.

        Arguments:
            - condition         - condition to be met, required (str)
            - threshold         - threshold to consider a column of interest,
                                  required (float)

        The condition should be 'greater' if we are looking for columns with
        a high degree of conservation or 'less' if we are looking for columns
        with a high degree of variation.
        """
        if not isinstance(condition, str):
            raise TypeError('"condition" argument should be a string')
        if not isinstance(threshold, float):
            raise TypeError('"threshold" argument should be a float')

        condition = condition.lower()

        # We create a data structure containing a list of elements with each
        # of them looking like this: (num_column, ci, freqs_residues)
        num_columns = range(self._start_column,
                            self._start_column + len(self._cis))
        stats = list(zip(num_columns,
                         self._cis,
                         [list(zip(self._freqs.keys(), values))
                                  for values in zip(*(self._freqs.values()))]))

        consensus_seq = []
        condition_met = []
        wanted_columns = []
        summary = []
        for stats_column in stats:
            # If we have reach the given number of residues per line of the
            # report, we update the summary been generated accordingly
            if (stats_column[0] - self._start_column) % \
                self._ITEMS_PER_LINE == 0:
                summary += consensus_seq + condition_met
                summary.append('\n{:5d}:\t'.format(stats_column[0]))
                consensus_seq = []
                condition_met = ['\n      \t']

            consensus_seq.append(max(stats_column[2], key=itemgetter(1))[0])

            self._check_condition('basic', condition, threshold, stats_column,
                                  wanted_columns, condition_met)

        # We have to append the last section of the consensus sequence at
        # the end of the report been generated 
        summary += consensus_seq + condition_met + ['\n\n']

        self._summarize_wanted_columns(summary, condition, threshold,
                                       wanted_columns)

        # We want to get rid of the first '\n'
        summary[0] = summary[0][1:]

        return "".join(summary)

    def generate_detailed(self, condition, threshold):
        """
        Generate a detailed report.

        Arguments:
            - condition         - condition to be met, required (str)
            - threshold         - threshold to consider a column of interest,
                                  required (float)

        The condition should be 'greater' if we are looking for columns with
        a high degree of conservation or 'less' if we are looking for columns
        with a high degree of variation.
        """
        if not isinstance(condition, str):
            raise TypeError('"condition" argument should be a string')
        elif condition not in ['greater', 'less']:
            raise ValueError(('"condition" argument has an invalid value. It '
                              "should be 'greater' or 'less'"))
        if not isinstance(threshold, float):
            raise TypeError('"threshold" argument should be a float')

        condition = condition.lower()

        # We create a data structure containing a list of elements with each
        # of them looking like this: (num_column, ci, freqs_residues)
        num_columns = range(self._start_column,
                            self._start_column + len(self._cis))
        stats = list(zip(num_columns,
                         self._cis,
                         [list(zip(self._freqs.keys(), values))
                                  for values in zip(*(self._freqs.values()))]))

        wanted_columns = []
        summary = []
        for stats_column in stats:
            summary.append('\n{:5d}:\t{:5.4f}\t\t'.format(stats_column[0],
                                                          stats_column[1]))
            sorted_freqs = sorted(stats_column[2], key=itemgetter(0))
            #FIXME Should I print elem.upper()?
            summary.append('\t'.join(("'{:s}':\t{:8.4f}%"
                                      ''.format(elem[0], elem[1]*100))
                                                for elem in sorted_freqs \
                                                if elem[0] != '-'))
            self._check_condition('detailed', condition, threshold, stats_column,
                                  wanted_columns)

        summary.append('\n\n')
        self._summarize_wanted_columns(summary, condition, threshold,
                                       wanted_columns)

        # We want to get rid of the first '\n'
        summary[0] = summary[0][1:]

        return "".join(summary)

