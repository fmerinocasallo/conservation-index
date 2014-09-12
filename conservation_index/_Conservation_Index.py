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
"""Calculate the CI for each column of the given set of sequences."""
from itertools import tee
from operator import itemgetter
from threading import Thread

from Bio import Align, AlignIO

class Conservation_Index(object):
    """

    """
    def __init__(self, seqs_type):
        """
        Creates a Conservation_Index.

        Arguments:
            - seqs_type         - sequences type, required (str)

        The sequence type should be one of those in PhyloDag.Data
        """
        if not isinstance(seqs_type, str):
            raise TypeError('"seqs_type" argument should be a string')

        self._seqs_type = seqs_type.lower()

        # _weights is a dictionary which stores the weights associated with
        # each element of the selected alphabet. Because we still does not
        # know how many sequences we are dealing with, we set all weights to
        # zero.
        if self._seqs_type == 'nucleotides':
            self._weights = {'-': [['-'], 0],
                             'A': [['A'], 0],
                             'C': [['C'], 0],
                             'G': [['G'], 0],
                             'T': [['T'], 0],
                             'R': [['A', 'G'], 0],
                             'Y': [['C', 'T'], 0],
                             'W': [['A', 'T'], 0],
                             'S': [['C', 'G'], 0],
                             'M': [['A', 'C'], 0],
                             'K': [['G', 'T'], 0],
                             'H': [['A', 'C', 'T'], 0],
                             'B': [['C', 'G', 'T'], 0],
                             'V': [['A', 'C', 'G'], 0],
                             'D': [['A', 'G', 'T'], 0],
                             'N': [['A', 'C', 'G', 'T'], 0]}
        elif self._seqs_type == 'amino acids':
            self._weights = {'=': [['='],0],
                             'G': [['G'], 0],
                             'A': [['A'], 0],
                             'V': [['V'], 0],
                             'L': [['L'], 0],
                             'I': [['I'], 0],
                             'P': [['P'], 0],
                             'F': [['F'], 0],
                             'Y': [['Y'], 0],
                             'C': [['C'], 0],
                             'M': [['M'], 0],
                             'H': [['H'], 0],
                             'K': [['K'], 0],
                             'R': [['R'], 0],
                             'W': [['W'], 0],
                             'S': [['S'], 0],
                             'T': [['T'], 0],
                             'D': [['D'], 0],
                             'E': [['E'], 0],
                             'N': [['N'], 0],
                             'Q': [['Q'], 0],
                             'B': [['D', 'N'], 0],
                             'Z': [['E', 'Q'], 0],
                             '-': [['-'], 0],
                             'X': [['X'], 0]}
        else:
            raise ValueError(('"seqs_type" argument is not a valid sequences '
                              'type'))

    def set_weights(self, num_seqs):
        """
        Set the correct weights associated with each element of the selected
        alpabet.

        Arguments:
            - num_seqs          - number of sequences to be analyzed,
                                  required (int)
        """
        if not isinstance(num_seqs, int):
            raise TypeError('"num_seqs" argument should be a int')

        for elem in self._weights:
            num_elem = len(self._weights[elem][0])
            self._weights[elem][1] = 1 / (num_elem * num_seqs)

    def analyze_columns(self, seqs, section, stats):
        """
        Analyze the given section (usually more than a single column) of a set
        of sequences.

        Arguments:
            - seqs              - set of sequences to be analyzed,
                                  required (MultipleSeqAlignment)
            - section           - section to be analyzed of each sequence,
                                  required (tuple)
            - stats             - conservation index stats,
                                  required (dict of lists of int)
        """
        if not isinstance(seqs, Align.MultipleSeqAlignment):
            raise TypeError('"seqs" argument should be a MultipleSeqAlignment')
        if not isinstance(section, tuple):
            raise TypeError('"section" argument should be a tuple')
        if not isinstance(stats, dict):
            raise TypeError('"stats" argument should be a dict')
        elif not all(isinstance(elem, list) for elem in stats.values()):
            raise TypeError('"stats" argument should be a dict of lists')

        start_column = section[0]
        if seqs.get_alignment_length() < section[1]:
            end_column = seqs.get_alignment_length()
        else:
            end_column = section[1]

        for seq in seqs:
            num_column = start_column
            for column in seq.seq[start_column:end_column].upper():
                for elem in self._weights[column][0]:
                    try:
                        stats[elem][num_column] += self._weights[column][1]
                    except TypeError:
                        raise TypeError(('"stats" argument should be a '
                                         'dict of lists of int'))
                num_column += 1

class CI_Thread(Thread):
    """

    """
    def __init__(self, seqs_type, seqs, section, stats):
        """
        Creates a CI_Thread object.

        Arguments:
            - seqs_type         - sequences type, required (str)
            - seqs              - set of sequences to be analyzed,
                                  required (MultipleSeqAlignment)
            - section           - section to be analyzed of each sequence,
                                  required (tuple)
            - stats             - conservation index stats,
                                  required (dict of lists of int)

        The sequence type should be one of those in PhyloDag.Data
        """
        if not isinstance(seqs_type, str):
            raise TypeError('"seqs_type" argument should be a string')
        if not isinstance(seqs, Align.MultipleSeqAlignment):
            raise TypeError('"seqs" argument should be a MultipleSeqAlignment')
        if not isinstance(section, tuple):
            raise TypeError('"section" argument should be a tuple')
        if not isinstance(stats, dict):
            raise TypeError('"stats" argument should be a dict')
        elif not all(isinstance(elem, list) for elem in stats.values()):
            raise TypeError('"stats" argument should be a dict of lists')
        # We assign to each thread the section of sequences to be analyzed
        # by itself and the data structure which stores the results ('stats')
        Thread.__init__(self)
        self._ci = Conservation_Index(seqs_type)
        self._seqs = seqs
        self._section = section
        self._stats = stats

        # We Initialize the weights data structures
        self._ci.set_weights(len(self._seqs))

    def run(self):
        """
        Function to be called when each thread starts its execution.
        """
        # Each thread have to analyze a given set of columns.
        self._ci.analyze_columns(self._seqs, self._section, self._stats)

class Report:
    def __init__(self, seqs_type, stats, start_column=0, ITEMS_PER_LINE=50):
        """
        Creates a Report object.

        Arguments:
            - seqs_type         - sequences type, required (str)
            - stats             - conservation index stats,
                                  required (dict of lists of int)
            - start_column      - absolute position of the 1st column,
                                  optional (int)
            - ITEMS_PER_LINE    - number of columns per line,
                                  optional (int)

        The sequence type should be one of those in PhyloDag.Data and the
        start_column should be 0 if we are dealing with a global set of
        sequences or the absolute number of the 1st column of the gene.
        """
        if not isinstance(seqs_type, str):
            raise TypeError('"seqs_type" argument should be a string')
        if seqs_type not in ['nucleotides', 'amino acids']:
            raise ValueError(('"seqs_type" argument has an invalid value '
                              "it should be 'nucleotides' or 'amino acids'"))
        if not isinstance(stats, dict):
            raise TypeError('"stats" argument should be a dict')
        elif not all(isinstance(elem, list) for elem in stats.values()):
            raise TypeError('"stats" argument should be a dict of lists')
        if not isinstance(start_column, int):
            raise TypeError('"start_column" argument should be an int')
        if not isinstance(ITEMS_PER_LINE, int):
            raise TypeError('"ITEMS_PER_LINE" argument should be an int')


        self._seqs_type = seqs_type
        self._stats = stats
        self._start_column = start_column
        self._ITEMS_PER_LINE = ITEMS_PER_LINE

    def _get_data_for_column(self, condition, column, report_type):
        """
        Return the data for the given column.

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

        condition = condition.lower()
        report_type = report_type.lower()
        if condition == 'greater' or condition == 'less':
            if report_type == 'basic':
                max_values, special_values = tee(column, 2)
                mfe = max(max_values, key=itemgetter(1))

                return mfe, special_values
            elif report_type == 'detailed':
                max_values, join_values, special_values = tee(column, 3)
                mfe = max(max_values, key=itemgetter(1))

                return mfe, join_values, special_values
            else:
                raise ValueError(('"report_type" argument has an invalid '
                                  "value it should be 'basic' or 'detailed'"))
        else:
            raise ValueError(('"condition" argument has an invalid value '
                              "it should be 'greater' or 'less'"))

    def _generate_header_second_mod(self, condition, threshold,
                                    special_columns):
        """
        Return the header for the report's 2nd module.

        Arguments:
            - condition         - condition to be met, required (str)
            - threshold         - threshold to consider a column of
                                  interest, required (float)
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
        if not isinstance(special_columns, list):
            raise TypeError('"special_columns" argument should be a list')
        elif not all(isinstance(elem, tuple) for elem in special_columns):
            raise TypeError('"stats" argument should be a list of tuples')

        condition = condition.lower()
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
                                      "value it should be 'greater' or "
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
                                      "value it should be 'greater' or "
                                      "'less'"))

            num = 0
            for record in special_columns:
                elem_stats = ("'{:s}': {:8.4f}%"
                              ''.format(elem[0], elem[1] * 100)
                                        for elem in record[1] \
                                        if elem[1] * 100 > 0.0)
                #FIXME Remember to use the mod operator when printing
                # the column
                header += ('> {:5d}: ({:s})\n'
                           ''.format(record[0], ', '.join(elem_stats)))
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
                raise ValueError(('"condition" argument has an invalid value '
                                  "it should be 'greater' or 'less'"))

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

        report = [zip(self._stats.keys(), values)
                  for values in zip(*(self._stats.values()))]

        # In order to traverse the stats' data structure just once,
        # we build the differents modules of the report at the same time
        # 1st module:
        #   - Most Frequent Sequence (MFS)
        # 2nd module:
        #   - Columns which meet the given condition
        first_mod_str = ''
        first_mod_str_aux = ''
        second_mod_str = ''

        # Because we want to notify the user with the CI distribution of each
        # column which meet the given condition, we store them as follows:
        # (pos, it)
        special_columns = []

        condition = condition.lower()
        num_column = self._start_column
        for column in report:
            mfe, special_values = self._get_data_for_column(condition, column,
                                                            'basic')

            if (num_column - self._start_column) % self._ITEMS_PER_LINE == 0:
                first_mod_str += first_mod_str_aux
                first_mod_str += '\n{:5d}:\t'.format(num_column)
                first_mod_str_aux = '\n      \t'
            else:
                pass

            first_mod_str += '{:s}'.format(mfe[0])

            if condition == 'greater':
                if mfe[1] > threshold:
                    first_mod_str_aux += '+'
                    stats = (num_column, special_values)
                    special_columns.append(stats)
                else:
                    first_mod_str_aux += '-'
            elif condition == 'less':
                if mfe[1] < threshold:
                    first_mod_str_aux += '+'
                    stats = (num_column, special_values)
                    special_columns.append(stats)
                else:
                    first_mod_str_aux += '-'
            else:
                raise ValueError(('"condition" argument has an invalid value '
                                  "it should be 'greater' or 'less'"))

            num_column += 1

        first_mod_str += first_mod_str_aux
        second_mod_str = self._generate_header_second_mod(condition,
                                                          threshold,
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
            raise ValueError(('"condition" argument has an invalid value it '
                              "should be 'greater' or 'less'"))
        if not isinstance(threshold, float):
            raise TypeError('"threshold" argument should be a float')

        report = [zip(self._stats.keys(), values)
                  for values in zip(*(self._stats.values()))]

        # In order to traverse the stats' data structure just once,
        # we build the differents modules of the report at the same time
        # 1st module:
        #   - Distribution of the CI for each column
        # 2nd module:
        #   - Columns which meet the given condition
        first_mod_str = ''
        second_mod_str = ''

        # Because we want to notify the user with the CI distribution of each
        # column which meet the given condition, we store them as follows:
        # (pos, it)
        special_columns = []

        condition = condition.lower()
        num_column = self._start_column
        for column in report:
            mfe, join_values, \
            special_values = self._get_data_for_column(condition, column,
                                                       'detailed')

            first_mod_str += '{:5d}:\t\t'.format(num_column)
            first_mod_str += '\t'.join(("'{:s}':\t{:8.4f}%"
                                        ''.format(elem[0], elem[1]*100))
                                       for elem in join_values \
                                       if elem[0] != '-')
            first_mod_str += '\n'
            # We store which columns meet the given condition
            if (condition == 'greater' and \
                    mfe[1] > threshold) or \
                (condition == 'less' and \
                    mfe[1] < threshold):
                stats = (num_column, special_values)
                special_columns.append(stats)
            else:
                pass

            num_column += 1

        second_mod_str = self._generate_header_second_mod(condition,
                                                          threshold,
                                                          special_columns)

        return first_mod_str + second_mod_str + '\n'
