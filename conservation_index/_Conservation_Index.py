#-------------------------------------------------------------------------------
# File :  _Conservation_Index.py
# Description :  Definition and implementation of the classes
#                'Conservation_Index', 'PT_Thread' and 'Report'.
#
# Author :  F. Merino-Casallo  ( fmerino@unizar.es )
# Last version :  v1.1 ( 03/Sep/2014 )
#-------------------------------------------------------------------------------
# Historical report :
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
from itertools import tee
from math import log
from operator import itemgetter
from threading import Thread

from Bio import Align, AlignIO

def _check_section(align, section_type, section):
    """
    Check if the given section is well formated. Return its boundaries.

    Arguments:
        - align                 - set of sequences to be analyzed,
                                  required (MultipleSeqAlignment)
        - section_type          - type of section to be checked,
                                  required (str)
        - section               - section of each sequence to be analyzed,
                                  required (tuple)
    """
    def check_max_boundary(max_align, max_section):
        """
        Check the upper boundary of the given section.

        Arguments:
        - max_align             - upper boundary of the set of sequences,
                                  required (int)
        - max_section           - upper boundary of the given section,
                                  required (int)
        """
        #FIXME Should I check this, again?
        # if not isinstance(max_align, int):
        #    raise TypeError('"max_align" argument should be an integer')
        # if not isinstance(max_section, int):
        #    raise TypeError('"max_section" argument should be an integer')
        if max_align < max_section:
            # warning(("The alignment's length is {:d}, "
            #         '"section" argument is outside its bonds. I have '
            #         'modified the section to fix this'.format(max_align)))

            return max_align
        else:
            return max_section

    if not isinstance(align, Align.MultipleSeqAlignment):
        raise TypeError('"align" argument should be a MultipleSeqAlignment')
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
            num_rows = len(align)

            return section[0], check_max_boundary(num_rows, section[1])
        elif section_type == 'columns':
            num_columns = align.get_alignment_length()

            return section[0], check_max_boundary(num_columns, section[1])
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
            - seqs_type         - type of sequences, required (str)

        The type of sequences should be one of those in PhyloDag.Data
        """
        if not isinstance(seqs_type, str):
            raise TypeError('"seqs_type" argument should be a string')

        self._seqs_type = seqs_type.lower()

        # _weights is a dictionary which stores the weights associated with
        # each element of the selected alphabet. Because we still does not
        # know how many sequences we are dealing with, we set all weights to
        # zero.
        if self._seqs_type == 'nucleotides':
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
        elif self._seqs_type == 'amino acids':
            self._elem_weights = {'=': [['='],0],
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

    def set_elem_weights(self, num_seqs, freq_method):
        """
        Set the correct weight associated with each element of the selected
        alpabet.

        Arguments:
            - num_seqs          - number of sequences to be analyzed,
                                  required (int)
            - freq_method       - frequencies' calculation method,
                                  required (str)
        """
        def unweighted_frequencies():
            """
            Set the correct element's weight for the unweighted frequencies'
            method.
            """
            for elem, weights in self._elem_weights:
                num_elems = len(weights[0])
                self._elem_weights[elem][1] = 1 / (num_elems * num_seqs)

        def weighted_frequencies():
            """
            Set the correct element's weight for the weighted frequencies'
            method.
            """
            for elem, weights in self._elem_weights.items():
                num_elems = len(weights[0])
                self._elem_weights[elem][1] = 1 / num_elems

        if not isinstance(num_seqs, int):
            raise TypeError('"num_seqs" argument should be a int')
        if not isinstance(freq_method, str):
            raise TypeError('"freq_method" argument should be a string')

        if freq_method == 'unweighted':
            unweighted_frequencies()
        elif freq_method == 'weighted':
            weighted_frequencies()
        else:
            raise ValueError(('"freq_method" argument has an invalid value. '
                              "It should be 'unweighted' or 'weighted'"))

    def weigh_align(self, align, section, align_weights):
        """
        Weigh the given section (usually more than a single row) of a set
        of sequences. Our aim is to compensate for over-representation among
        multiple aligned sequences.

        Arguments:
            - align             - set of sequences to be analyzed,
                                  required (MultipleSeqAlignment)
            - section           - section of each sequence to be analyzed,
                                  required (tuple)
            - align_weights     - sequences' weights, required (list of float)
        """
        if not isinstance(align, Align.MultipleSeqAlignment):
            raise TypeError('"align" argument should be a MultipleSeqAlignment')
        if not isinstance(section, tuple):
            raise TypeError('"section" argument should be a tuple')
        if not isinstance(align_weights, list):
            raise TypeError('"align_weights" argument should be a list')

        start_row, end_row = _check_section(align, 'rows', section)

        # freqs stores the residues (nucleotides or amino acids) which appear
        # in each column and how many times they do so
        freqs = [{} for k in range(0, align.get_alignment_length())]
        for record in align[start_row:end_row]:
            for j, residue in enumerate(record.seq.lower()):
                freqs[j][residue] = freqs[j].get(residue, 0) + 1
                # for elem in self._elem_weights[residue][0]:
                #     freqs[j][elem] = freqs[j].get(elem, 0) + \
                #                      1 / len(self._elem_weights[residue][0])

        # estimate the weight associate with each sequence
        for i, record in enumerate(align[start_row:end_row]):
            for j, residue in enumerate(record.seq.lower()):
                # reps = 0
                # for elem in self._elem_weights[residue][0]:
                #     reps += freqs[j][elem] / len(self._elem_weights[residue][0])
                try:
                    align_weights[i] += 1 / (len(freqs[j]) * freqs[j][residue])
                except TypeError:
                    raise TypeError(('"align_weights" argument should be a '
                                     'list of ints'))

    def estimate_frequencies(self, align, section, freqs,
                             freq_method, align_weights=[]):
        """
        Estimate frequencies of the given section (usually more than a single
        column) of a set of sequences.

        Arguments:
            - align             - set of sequences to be analyzed,
                                  required (MultipleSeqAlignment)
            - section           - section of each sequence to be analyzed,
                                  required (tuple)
            - freqs             - estimated residues' frequencies for each
                                  column, required (dict of lists of int)
            - freq_method       - frequencies' estimation method,
                                  required (str)
            - align_weights     - sequences' weights, optional (list of float)
        """
        def unweighted_frequencies():
            """
            Estimate frequencies using the unweighted frequencies' method.
            """
            for record in align:
                for j, residue in enumerate(record.seq[start_column: \
                                                    end_column].lower()):
                    for elem in self._elem_weights[residue][0]:
                        try:
                            freqs[elem][start_column + j] += \
                                            self._elem_weights[residue][1]
                        except TypeError:
                            raise TypeError(('"freqs" argument should be a '
                                             'dict of lists of int'))

        def weighted_frequencies():
            """
            Estimate frequencies using the weighted frequencies' method.
            """
            weights_sum = sum(seq_weight for seq_weight in align_weights)

            for i, record in enumerate(align):
                for j, residue in enumerate(record.seq[start_column: \
                                                    end_column].lower()):
                    for elem in self._elem_weights[residue][0]:
                        try:
                            freqs[elem][j + start_column] += \
                                (align_weights[i] * \
                                 self._elem_weights[residue][1]) / \
                                 weights_sum
                        except TypeError:
                            raise TypeError(('"freqs" argument should be a '
                                             'dict of lists of int'))
                        except IndexError:
                            raise IndexError(('"align_weights" argument should '
                                              'have the same size as the '
                                              "alignment's length"))

        if not isinstance(align, Align.MultipleSeqAlignment):
            raise TypeError('"align" argument should be a MultipleSeqAlignment')
        if not isinstance(section, tuple):
            raise TypeError('"section" argument should be a tuple')
        if not isinstance(freqs, dict):
            raise TypeError('"freqs" argument should be a dict')
        elif not all(isinstance(elem, list) for elem in freqs.values()):
            raise TypeError('"freqs" argument should be a dict of lists')
        if not isinstance(freq_method, str):
            raise TypeError('"freq_method" argument should be a string')
        if not isinstance(align_weights, list):
            raise TypeError('"align_weights" argument should be a list')

        start_column, end_column = _check_section(align, 'columns', section)

        if freq_method == 'unweighted':
            unweighted_frequencies()
        elif freq_method == 'weighted':
            weighted_frequencies()
        else:
            raise ValueError(('"freq_method" argument has an invalid value. '
                              "It should be 'unweighted' or 'weighted'"))

    def estimate_conservation_index(self, align, section, freqs,
                                     ci_method, cis, align_weights=[]):
        """
        Estimate conservation index of the given section (usually more than a
        single column) of a set of sequences.

        Arguments:
            - align             - set of sequences to be analyzed,
                                  required (MultipleSeqAlignment)
            - section           - section of each sequence to be analyzed,
                                  required (tuple)
            - freqs             - estimated residues' frequencies for each
                                  column, required (dict of lists of int)
            - ci_method         - conservation index's calculation method,
                                  required (str)
            - cis               - estimated residues' conservation indices
                                  for each column, required (list of float)
            - align_weights     - sequences' weights, optional (list of float)
        """
        def shannon_entropy():
            """
            Estimate conservation index using the Shannon Entropy method.
            """
            lambda_t = log(min(len(align), len(freqs.keys())))
            for freqs_residue in freqs.values():
                for j, freq_column in \
                            enumerate(freqs_residue[start_column:end_column]):
                    if freq_column != 0.0:
                        cis[j + start_column] += freq_column * log(freq_column)

            # We want scale the entropy to range [0, 1]
            for num_column in range(start_column, end_column):
                cis[num_column] = 1 - (-1 * cis[num_column] * lambda_t)

        def variance():
            """
            """
            pass

        if not isinstance(align, Align.MultipleSeqAlignment):
            raise TypeError('"align" argument should be a MultipleSeqAlignment')
        if not isinstance(section, tuple):
            raise TypeError('"section" argument should be a tuple')
        if not isinstance(freqs, dict):
            raise TypeError('"freqs" argument should be a dict')
        elif not all(isinstance(elem, list) for elem in freqs.values()):
            raise TypeError('"freqs" argument should be a dict of lists')
        if not isinstance(ci_method, str):
            raise TypeError('"ci_method" argument should be a string')
        if not isinstance(cis, list):
            raise TypeError('"cis" argument should be a list')
        if not isinstance(align_weights, list):
            raise TypeError('"align_weights" argument should be a list')

        start_column, end_column = _check_section(align, 'columns', section)

        if ci_method == 'shannon entropy':
            shannon_entropy()
        elif freq_method == 'variance':
            variance()
        else:
            raise ValueError(('"ci_method" argument has an invalid value. '
                              "It should be 'shannon entropy' or 'variance'"))


class CI_Thread(Thread):
    """

    """
    def __init__(self, seqs_type, align, columns_section,
                 freq_method, freqs, ci_method, cis,
                 rows_section=(), align_weights=[]):
        """
        Creates a CI_Thread object.

        Arguments:
            - seqs_type         - type of sequences, required (str)
            - align             - set of sequences to be analyzed,
                                  required (MultipleSeqAlignment)
            - columns_section   - section of each sequence to be analyzed
                                  for the frequencies' count, required (tuple)
            - freq_method       - frequencies' calculation method,
                                  required (str)
            - freqs             - estimated residues' frequencies for each
                                  column, required (dict of lists of int)
            - ci_method         - conservation index' calculation method,
                                  required (str)
            - cis               - estimated residues' conservation indices
                                  for each column, required (list of float)
            - rows_section      - section of each sequence to be analyzed
                                  for the align's weight, optional (tuple)
            - align_weights     - sequences' weights, optional (list of float)

        The type of sequences should be one of those in PhyloDag.Data. The
        frequencies' stats should be a dict which has as keys the elements of
        the related alphabet (IUPACUnambiguousDNA for nucleotides and
        IUPACProtein for amino acids) plus the gap. The values associated with
        each of these keys should be 0.0. The conservation index for each 
        column is going to be stored in a list of floats called 'cis' and it
        should be initialize with as many 0.0 as columns in the alignment.
        The sequences' weights should be a list of 0.0 with the same size as
        the number of sequences in the alignment
        """
        if not isinstance(seqs_type, str):
            raise TypeError('"seqs_type" argument should be a string')
        if seqs_type.lower() not in ['nucleotides', 'amino acids']:
            raise ValueError(('"seqs_type" argument has an invalid value. '
                              "It should be 'nucleotides' or 'amino acids'"))
        if not isinstance(align, Align.MultipleSeqAlignment):
            raise TypeError('"align" argument should be a MultipleSeqAlignment')
        if not isinstance(columns_section, tuple):
            raise TypeError('"columns_section" argument should be a tuple')
        if not isinstance(freq_method, str):
            raise TypeError('"freq_method" argument should be a string')
        if freq_method.lower() not in ['unweighted', 'weighted']:
            raise ValueError(('"freq_method" argument has an invalid value. '
                              "It should be 'unweighted' or 'weighted'"))
        if not isinstance(freqs, dict):
            raise TypeError('"freqs" argument should be a dict')
        elif not all(isinstance(elem, list) for elem in freqs.values()):
            raise TypeError('"freqs" argument should be a dict of lists')
        if not isinstance(ci_method, str):
            raise TypeError('"ci_method" argument should be a string')
        if ci_method.lower() not in ['shannon entropy', 'variance']:
            raise ValueError(('"ci_method" argument has an invalid value. '
                              "It should be 'shannon entropy' or 'variance'"))
        if not isinstance(cis, list):
            raise TypeError('"cis" argument should be a list')
        if not isinstance(rows_section, tuple):
            raise TypeError('"rows_section" argument should be a tuple')
        if not isinstance(align_weights, list):
            raise TypeError('"align_weights" argument should be a list')

        # We assign to each thread the section of the align to be analyzed
        # by itself and the data structure which stores the results ('freqs')
        Thread.__init__(self)
        self._ci = Conservation_Index(seqs_type)
        self._align = align
        self._columns_section = columns_section
        self._freq_method = freq_method.lower()
        self._freqs = freqs
        self._ci_method = ci_method.lower()
        self._cis = cis

        # We Initialize the weights' data structures
        self._ci.set_elem_weights(len(self._align), self._freq_method)

        if self._freq_method == 'weighted':
            self._rows_section = rows_section
            self._align_weights = align_weights
            self._ci.weigh_align(self._align, self._rows_section,
                                 self._align_weights)

    def run(self):
        """
        Function to be called when each thread starts its execution.
        """
        # Each thread have to analyze a given set of columns.
        if self._freq_method == 'weighted':
            self._ci.estimate_frequencies(self._align, self._columns_section,
                                          self._freqs, self._freq_method,
                                          self._align_weights)
            self._ci.estimate_conservation_index(self._align,
                                                 self._columns_section,
                                                 self._freqs,
                                                 self._ci_method, self._cis,
                                                 self._align_weights)
        else: # self_freq_method == 'unweighted'
            self._ci.estimate_frequencies(self._align, self._columns_section,
                                          self._freqs, self._freq_method)
            self._ci.estimate_conservation_index(self._align,
                                                 self._columns_section,
                                                 self._freqs,
                                                 self._ci_method, self._cis)

class Report:
    def __init__(self, seqs_type, freqs, cis, start_column=0, ITEMS_PER_LINE=50):
        """
        Creates a Report object.

        Arguments:
            - seqs_type         - type of sequences, required (str)
            - freqs             - estimated residues' frequencies for each
                                  column, required (dict of lists of int)
            - cis               - estimated conservation indices for each
                                  column, required (list of float)
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
        if seqs_type.lower() not in ['nucleotides', 'amino acids']:
            raise ValueError(('"seqs_type" argument has an invalid value. '
                              "It should be 'nucleotides' or 'amino acids'"))
        if not isinstance(freqs, dict):
            raise TypeError('"freqs" argument should be a dict')
        elif not all(isinstance(residue, list) for residue in freqs.values()):
            raise TypeError('"freqs" argument should be a dict of lists')
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

    def _get_data_for_column(self, condition, column, report_type):
        """
        Return the data needed to generate the chosen report for the given
        column.

        Arguments:
            - condition         - condition to be met, required (str)
            - column            - column extracted from the report,
                                  required (zip)
            - report_type       - type of report, required (str)

        The condition should be 'greater' if we are looking for columns with
        a high degree of conservation or 'less' if we are looking for columns
        with a high degree of variation. The report type should be 'basic'
        or 'detailed'.
        """
        if not isinstance(condition, str):
            raise TypeError('"condition" argument should be a str')
        if not isinstance(column, zip):
            raise TypeError('"column" argument should be a zip')
        if not isinstance(report_type, str):
            raise TypeError('"report_type" argument should be a str')

        if condition == 'greater' or condition == 'less':
            if report_type == 'basic':
                max_values, special_values = tee(column, 2)
                # Retrieve the Most Frequent Residue for this column
                mfr = max(max_values, key=itemgetter(1))

                return mfr, special_values
            elif report_type == 'detailed':
                max_values, join_values, special_values = tee(column, 3)
                mfr = max(max_values, key=itemgetter(1))

                return mfr, join_values, special_values
            else:
                raise ValueError(('"report_type" argument has an invalid value. '
                                  "It should be 'basic' or 'detailed'"))
        else:
            raise ValueError(('"condition" argument has an invalid value. '
                              "It should be 'greater' or 'less'"))

    def _generate_header_second_mod(self, condition, threshold, cis,
                                    special_columns):
        """
        Return the header for the report's 2nd module.

        Arguments:
            - condition         - condition to be met, required (str)
            - threshold         - threshold to consider a column of
                                  interest, required (float)
            - cis               - estimated conservation indices for each
                                  column, required (list of float)
            - special_columns   - columns which met the given condition
                                  required (list of tuples)

        The condition should be 'greater' if we are looking for columns with
        a high degree of conservation or 'less' if we are looking for columns
        with a high degree of variation. The special_columns argument should
        only include columns which meet the given condition.
        """
        if not isinstance(condition, str):
            raise TypeError('"condition" argument should be a str')
        if not isinstance(threshold, float):
            raise TypeError('"threshold" argument should be a float')
        if not isinstance(cis, list):
            raise TypeError('"cis" argument should be a list')
        if not isinstance(special_columns, list):
            raise TypeError('"special_columns" argument should be a list')
        elif not all(isinstance(elem, tuple) for elem in special_columns):
            raise TypeError('"freqs" argument should be a list of tuples')

        if special_columns:
            if len(special_columns) > 1:
                if condition == 'greater':
                    header = ('\nThere are {:d} columns with a high degree of '
                              'conservation (>{:6.2f}%) in the {:s} sequences:'
                              '\n\n'.format(len(special_columns),
                                            threshold * 100, self._seqs_type))
                elif condition == 'less':
                    header = ('\nThere are {:d} columns with a high degree of '
                              'variation (<{:6.2f}%) in the {:s} sequences:'
                              '\n\n'.format(len(special_columns),
                                            threshold * 100, self._seqs_type))
                else:
                    raise ValueError(('"condition" argument has an invalid '
                                      "value. It should be 'greater' or "
                                      "'less'"))
            else:
                if condition == 'greater':
                    header = ('\nThere is 1 column with a high degree of '
                              'conservation (>{:6.2f}%) in the {:s} sequences:'
                              '\n\n'.format(threshold * 100, self._seqs_type))
                elif condition == 'less':
                    header = ('\nThere is 1 column with a high degree of '
                              'variation (<{:6.2f}%) in the {:s} sequences:'
                              '\n\n'.format(threshold * 100, self._seqs_type))
                else:
                    raise ValueError(('"condition" argument has an invalid '
                                      "value. It should be 'greater' or "
                                      "'less'"))

            num = 0
            for record in special_columns:
                sorted_record = sorted(record[1], key=itemgetter(0))
                #FIXME Should I print elem.upper()?
                elem_freqs = ("'{:s}': {:8.4f}%"
                              ''.format(elem[0], elem[1] * 100)
                                        for elem in sorted_record \
                                        if elem[1] * 100 > 0.0)
                #FIXME Remember to use the mod operator when printing
                # the column
                header += ('> {:5d}: {:8.4f} ({:s})\n'
                           ''.format(record[0], cis[record[0] - self._start_column],
                                     ', '.join(elem_freqs)))
            return header
        else:
            if condition == 'greater':
                return ('\nThere are not any columns with a high degree of '
                        'conservation (>{:6.2f}%) in the {:s} sequences.'
                        ''.format(threshold * 100, self._seqs_type))
            elif condition == 'less':
                return ('\nThere are not any columns with a high degree of '
                        'variation (<{:6.2f}%) in the {:s} sequences.'
                        ''.format(threshold * 100, self._seqs_type))
            else:
                raise ValueError(('"condition" argument has an invalid value. '
                                  "It should be 'greater' or 'less'"))

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
            raise TypeError('"condition" argument should be a str')
        if not isinstance(threshold, float):
            raise TypeError('"threshold" argument should be a float')

        report = [zip(self._freqs.keys(), values)
                  for values in zip(*(self._freqs.values()))]

        # In order to traverse the freqs' data structure just once,
        # we build the different report's modules at the same time
        # 1st module:
        #   - Consensus Sequence (CS)
        # 2nd module:
        #   - Columns which meet the given condition
        first_mod_str = ''
        first_mod_str_aux = ''
        second_mod_str = ''

        # Because we want to notify the user with the CI distribution of each
        # column which meet the given condition, we store them as follows:
        # (pos, iter)
        special_columns = []

        condition = condition.lower()
        num_column = self._start_column
        for column in report:
            mfr, special_values = self._get_data_for_column(condition, column,
                                                            'basic')

            if (num_column - self._start_column) % self._ITEMS_PER_LINE == 0:
                first_mod_str += first_mod_str_aux
                first_mod_str += '\n{:5d}:\t'.format(num_column)
                first_mod_str_aux = '\n      \t'
            else:
                pass

            first_mod_str += '{:s}'.format(mfr[0])

            if condition == 'greater':
                if mfr[1] > threshold:
                    first_mod_str_aux += '+'
                    freqs = (num_column, special_values)
                    special_columns.append(freqs)
                else:
                    first_mod_str_aux += '-'
            elif condition == 'less':
                if mfr[1] < threshold:
                    first_mod_str_aux += '+'
                    freqs = (num_column, special_values)
                    special_columns.append(freqs)
                else:
                    first_mod_str_aux += '-'
            else:
                raise ValueError(('"condition" argument has an invalid value. '
                                  "It should be 'greater' or 'less'"))

            num_column += 1

        first_mod_str += first_mod_str_aux
        second_mod_str = self._generate_header_second_mod(condition,
                                                          threshold,
                                                          self._cis,
                                                          special_columns)

        # We do not want to print a new line at the beginning of the report
        # so we skip the 1st character
        return first_mod_str[1:] + '\n' + second_mod_str + '\n'

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
            raise TypeError('"condition" argument should be a str')
        elif condition not in ['greater', 'less']:
            raise ValueError(('"condition" argument has an invalid value. It '
                              "should be 'greater' or 'less'"))
        if not isinstance(threshold, float):
            raise TypeError('"threshold" argument should be a float')

        report = [zip(self._freqs.keys(), values)
                  for values in zip(*(self._freqs.values()))]

        # In order to traverse the freqs' data structure just once,
        # we build the differents modules of the report at the same time
        # 1st module:
        #   - Distribution of the CI for each column
        # 2nd module:
        #   - Columns which meet the given condition
        first_mod_str = ''
        second_mod_str = ''

        # Because we want to notify the user with the CI distribution of each
        # column which meet the given condition, we store them as follows:
        # (pos, iter)
        special_columns = []

        condition = condition.lower()
        num_column = self._start_column
        for column in report:
            mfr, join_values, \
            special_values = self._get_data_for_column(condition, column,
                                                       'detailed')

            sorted_join_values = sorted(join_values, key=itemgetter(0))
            first_mod_str += '{:5d}:\t\t'.format(num_column)
            #FIXME Should I print elem.upper()?
            first_mod_str += '\t'.join(("'{:s}':\t{:8.4f}%"
                                        ''.format(elem[0], elem[1]*100))
                                       for elem in sorted_join_values \
                                       if elem[0] != '-')
            first_mod_str += '\n'
            # We store which columns meet the given condition
            if (condition == 'greater' and \
                    mfr[1] > threshold) or \
                (condition == 'less' and \
                    mfr[1] < threshold):
                freqs = (num_column, special_values)
                special_columns.append(freqs)
            else:
                pass

            num_column += 1

        second_mod_str = self._generate_header_second_mod(condition,
                                                          threshold,
                                                          self._cis,
                                                          special_columns)

        return first_mod_str + second_mod_str + '\n'

