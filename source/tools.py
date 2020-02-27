# import standard or third party modules
import os
from math import log
from tkinter import messagebox
from itertools import filterfalse

# import own modules
from source.configuration import cfg, adjust_dir_path

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2017 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'


def check_dir(name, path, files=False):
    """
    Checks if the directory exists and is a directory, and shows an error message if it is not.
    :param name: name of the directory
    :param path: directory path
    :param files: True if files in directory are to be checked as well, False otherwise
    :return: True if the directory is valid, False otherwise
    """
    path = adjust_dir_path(path)

    if not os.path.isdir(path):
        if not files:
            os.makedirs(path)
            return True
        else:
            messagebox.showerror('Error', 'The {0} directory does not exist.\n\n'
                                          'Directory path: {1}'.format(name, path))
            return False

    # if specified, check if the files in the directory are accessible for reading
    if files:
        for fn in os.listdir(path):
            if os.path.isfile(fn) and not os.access(path + fn, os.R_OK):
                messagebox.showerror('Error', 'A file in the {1} directory is not readable.\n\n'
                                              'File name: {0}\n\n'
                                              'Directory path: {2}'.format(fn, name, path))
                return False

    return True


def setup_directory(name, directory):
    """
    Sets up the specified directory if it does not exist yet and shows and error message if that fails.
    :param name: name of the directory
    :param directory: directory path
    :return: True if the directory already existed/could be created and is accessible, False otherwise
    """
    dir = adjust_dir_path(directory)

    # try to create the directory if it does not exist
    if not os.path.isdir(dir):
        try:
            os.makedirs(dir)
            return True
        except OSError:
            messagebox.showerror('Error', 'The {0} directory did not exist and could not be created.\n\n'
                                          'Directory path: {1}'.format(name, dir))
            return False

    # check if it is accessible for reading
    if not os.access(dir, os.R_OK):
        messagebox.showerror('Error', 'The {0} directory is not readable.\n\n'
                                      'Directory path: {1}'.format(name, dir))
        return False

    # check if it is accessible for writing
    if not os.access(dir, os.W_OK):
        messagebox.showerror('Error', 'The {0} directory is not writable.\n\n'
                                      'Directory path: {1}'.format(name, dir))
        return False

    return True


def clean_directory(name, directory, subdirs=False):
    """
    Deletes all files and subdirectories.
    :param name: name of the directory
    :param directory: directory path
    :return: True if all files/subdirectories could be delete, False otherwise
    """
    directory = adjust_dir_path(directory)

    # do nothing if it is not a valid directory
    if not check_dir(name, directory):
        return False

    failed = set()

    # try to delete all files and subdirectories unless they are exempt
    for fn in os.listdir(directory):
        # delete files
        if os.path.isfile(directory + fn):
            try:
                os.remove(directory + fn)
            except OSError:
                failed.add(fn)
        # delete subdirectories if specified
        elif os.path.isdir(directory + fn) and subdirs:
            if not clean_directory(name, directory + fn, subdirs):
                failed.add(fn)
            else:
                os.rmdir(directory + fn)

    if failed:
        msg = 'The following files in the {0} directory could not be deleted:\n- {2}\n\n' \
              'Directory path: {1}'
        messagebox.showerror('Error', msg.format(name, directory, '\n- '.join(sorted(failed))))
        return False

    return True


def check_tsv_fieldnames(name, path, supposed, actual):
    """
    Checks if all required field names are present in the .tsv file and shows an error message if not.
    :param name: name of the file
    :param path: file path
    :param supposed: field names that are supposed to be in the file
    :param actual: field names that are actually in the file
    :return: True if all field names are in the file, False otherwise
    """
    difference = set(supposed).difference(actual)

    if difference:
        msg = 'The following header entries are missing from the {0} file:\n- {2}\n\n' \
              'File path: {1}'
        messagebox.showerror('TSV Error', msg.format(name.lower().replace('file', '').strip(),
                                                     path, '\n- '.join(difference)))
        return False

    return True


def tsv_start(file, comments=tuple('#')):
    """
    Returns an iterator with all non-comment and non-empty lines.
    """
    return filterfalse(lambda line: line.startswith(comments) or not line.strip(), file)


def get_file_name(path):
    """
    Obtain the file name and file type from a path.
    :param path: file path
    :return: file name, file type
    """
    path = path.replace('\\', '/')

    # POSIX path
    if '/' in path:
        file = path.split('/')[-1]
    # file_name.file_type
    else:
        file = path

    parts = file.split('.')
    file_name = parts[0]
    file_type = '.' + '.'.join(parts[1:])

    return file, file_name, file_type


def combine_reg_info(a, b):
    """
    Combines the regulatory information.
    :param a: regulation description 1
    :param b: regulation description 2
    :return: updated regulatory information
    """
    # no difference
    if a == b:
        return a
    # b has more detailed information on the regulation
    if a == cfg.misc.unknown:
        return b
    # a has more detailed information on the regulation
    if b == cfg.misc.unknown:
        return a
    # regulation depends on environment and other factors => effector
    if {a, b} == {cfg.misc.activator, cfg.misc.repressor}:
        return cfg.misc.effector
    # regulation depends on environment and other factors => effector
    if cfg.misc.effector in {a, b}:
        return cfg.misc.effector


def parse_regulation(reg):
    """
    Parses the regulation input.
    :param reg: regulation information
    :return: processed regulation
    """
    reg = reg.lower()

    if ('activ' in reg) or (reg in cfg.isup.reg_act):
        return cfg.misc.activator

    if ('repre' in reg) or (reg in cfg.isup.reg_rep):
        return cfg.misc.repressor

    if reg in cfg.isup.reg_eff:
        return cfg.misc.effector

    return cfg.misc.unknown


def parse_strand(strand):
    """
    Parses strand information.
    """
    if strand in cfg.isup.str_minus:
        return cfg.misc.minus
    if strand in cfg.isup.str_plus or strand in cfg.misc.empty:
        return cfg.misc.plus
    return ''


def is_dna_sequence(seq):
    """
    Tests whether the sequence is a DNA sequence.
    :param seq: DNA sequence
    :return: True if the sequence is DNA, False otherwise
    """
    # sequence is not given
    if seq in cfg.misc.empty:
        return False

    for base in seq:
        if base not in ['A', 'T', 'C', 'G', 'N']:
            return False
    return True


def is_protein_sequence(seq):
    """
    Tests whether the sequence is a protein sequence.
    :param seq: amino acid sequence
    :return: True if the sequence consists of amino acids, False otherwise
    """
    # sequence is not given
    if seq in cfg.misc.empty:
        return False

    for amino_acid in seq:
        if amino_acid not in ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W',
                              'Y', 'V', 'X']:
            return False
    return True


def translate(dna):
    """
    Translates a DNA sequence into the corresponding amino acid sequence. The first codon is always translated to M to
    match the protein sequences in the gene files.
    :param dna: the DNA sequence to be translated
    :return: the amino acid sequence
    """
    dna_len = len(dna)

    # there needs to be at least one codon
    if dna_len < 3:
        return cfg.misc.none

    # the first amino acid is always methionine
    protein = ['M']

    # iterate over each codon and translate it
    for i in range(3, dna_len - 2, 3):
        codon = dna[i:i + 3]

        # translate the codon
        if codon in ('ATT', 'ATC', 'ATA'):
            amino = 'I'
        elif codon in ('CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'):
            amino = 'L'
        elif codon in ('GTT', 'GTC', 'GTA', 'GTG'):
            amino = 'V'
        elif codon in ('TTT', 'TTC'):
            amino = 'F'
        elif codon == 'ATG':
            amino = 'M'
        elif codon in ('TGT', 'TGC'):
            amino = 'C'
        elif codon in ('GCT', 'GCC', 'GCA', 'GCG'):
            amino = 'A'
        elif codon in ('GGT', 'GGC', 'GGA', 'GGG'):
            amino = 'G'
        elif codon in ('CCT', 'CCC', 'CCA', 'CCG'):
            amino = 'P'
        elif codon in ('ACT', 'ACC', 'ACA', 'ACG'):
            amino = 'T'
        elif codon in ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'):
            amino = 'S'
        elif codon in ('TAT', 'TAC'):
            amino = 'Y'
        elif codon == 'TGG':
            amino = 'W'
        elif codon in ('CAA', 'CAG'):
            amino = 'Q'
        elif codon in ('AAT', 'AAC'):
            amino = 'N'
        elif codon in ('CAT', 'CAC'):
            amino = 'H'
        elif codon in ('GAA', 'GAG'):
            amino = 'E'
        elif codon in ('GAT', 'GAC'):
            amino = 'D'
        elif codon in ('AAA', 'AAG'):
            amino = 'K'
        elif codon in ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'):
            amino = 'R'
        elif codon in ('TAA', 'TAG', 'TGA'):
            break
        else:
            amino = 'X'

        # add the amino acid to the sequence
        protein.append(amino)

    # return the amino acid sequence
    return ''.join(protein)


def complement_base(base):
    """
    Returns the complement DNA base.
    :param base: DNA base
    :return: complement DNA base.
    """
    if base == 'A':
        return 'T'
    if base == 'T':
        return 'A'
    if base == 'C':
        return 'G'
    if base == 'G':
        return 'C'
    raise ValueError


def complement_strand(seq):
    """
    Returns the DNA sequence of the complement strand. Each base is replaced with its complement and the sequence
    is reversed.
    :param seq: DNA sequence
    :return: DNA sequence on the complement strand
    """
    return ''.join([complement_base(b) for b in seq][::-1])


class Edge:
    """
    Class that represents an edge in the gene regulatory network (GRN).
    """
    def __init__(self, ID, source, target, reg):
        """
        Initialises the edge.
        :param ID: the ID of the edge
        :param source: the regulator
        :param target: the regulated gene
        :param reg: the regulation the edge represents
        """
        self.ID = ID
        self.source = source
        self.target = target
        self.reg = reg.upper()[0]

    def __eq__(self, other):
        if self.target != other.target:
            return False
        if self.source != other.source:
            return False
        if self.reg != other.reg:
            return False
        return True

    def __hash__(self):
        """
        Calculates the hash of an edge using the source, target and regulation.
        :return: hash
        """
        return hash((self.source, self.target, self.reg))

    def gml(self):
        """
        :return: GML representation of the edge
        """
        return '\tedge [\n' \
               '\t\tid {0}\n' \
               '\t\tsource {1}\n' \
               '\t\ttarget {2}\n' \
               '\t\tlabel \"{3}\"\n' \
               '\t]\n'.format(self.ID, self.source, self.target, self.reg)


class SubstitutionMatrix:
    def __init__(self, path):
        """
        Initialises the substitution matrix with the information given in the matrix file.
        :param path: path to the substitution matrix file
        """
        # setup the matrix
        self.matrix = {}
        self.max = {}
        self.min = {}
        self.error = None

        with open(path, 'r') as file:
            # initialise the column and row names
            row_names = []
            column_names = []
            column_line = True

            for line in file:
                # skip comments
                if line[0] == '#':
                    continue

                # process the column header
                if column_line:
                    column_line = False

                    # read the column names
                    column_names = line.rstrip().split('\t')[1:]

                    # test if there are enough columns for all amino acids
                    if len(set(column_names)) != 21 or len(column_names) != 21:
                        self.error = 'The number of amino acid columns is =/= 20.'
                        return

                    # test if the column names are amino acids
                    if not is_protein_sequence(''.join(column_names[:-1])):
                        self.error = 'The columns contain non-amino acids.'
                        return

                    # test if there is a column for gaps
                    if '*' not in column_names and '-' not in column_names and '/' not in column_names:
                        self.error = 'There is no column for gaps.'
                        return

                    # replace potential gap symbols with '*' for internal consistency
                    column_names = ['*' if x in ['*', '-', '/'] else x for x in column_names]

                    for key in column_names:
                        self.matrix[key] = {}

                # read the rows
                else:
                    elems = line.rstrip().split('\t')
                    amino = elems[0]
                    elems = elems[1:]

                    # test if the row has enough columns for all amino acids
                    if len(column_names) != len(elems):
                        self.error = 'Line {0} in the matrix is too short or too long.'.format(amino)
                        self.matrix = {}
                        return

                    # test if the row represents an amino acid
                    if not is_protein_sequence(amino) and amino not in ['*', '-', '/']:
                        self.error = '{0} is not an amino acid or gap.'.format(amino)
                        self.matrix = {}
                        return

                    row_names.append(amino)

                    # add the substitution scores to the matrix
                    for i in range(0, len(elems)):
                        self.matrix[amino][column_names[i]] = int(elems[i])

            # test if there is a row for every amino acid
            if len(set(row_names)) != 21 or len(row_names) != 21:
                self.error = 'The number of amino acid rows is =/= 20.'
                self.matrix = {}
                return

            # test if the matrix is symmetrical
            if not self._is_symmetrical():
                self.error = 'The matrix is not symmetrical.'
                self.matrix = {}
                return

            # compute the minimum and maximum score for every amino acid
            self._compute_min_max_scores()

    # SCORES
    def missense_score(self, ref, mut):
        """
        Computes the normalised substitution score of a single substitution (= missense mutation). 1 is the best
        possible score and 0 is the worst possible score.
        :param ref: reference amino acid
        :param mut: mutated amino acid
        :return: normalised substitution score
        """
        return (self._single_score(ref, mut) - self.min[ref]) / (self.max[ref] - self.min[ref])

    def score(self, ref, mut, pos=-1):
        """
        Computes the normalised alignment score of a mutated amino acid sequence with the reference amino acid
        sequence. 1 is the best possible score and 0 is the worst possible score.
        :param ref: reference amino acid sequence
        :param mut: mutated amino acid sequence
        :param pos: protein position after which the sequences differ
        :return: normalised alignment score of reference and mutate sequence
        """
        if ref in cfg.misc.empty and mut in cfg.misc.empty:
            return 1.0
        if ref in cfg.misc.empty:
            return 0.0
        if mut in cfg.misc.empty:
            return 0.0

        ref_len = len(ref)
        mut_len = len(mut)
        ref_old = ref

        # compute global alignment more efficiently by only aligning the changed part
        if ref_len != mut_len:
            if pos > -1:
                # otherwise ref and mut could have unequal length at the end
                pos = min(ref_len, mut_len) - 1
                ref_al, mut_al = self._compute_global_alignment(ref[pos:], mut[pos:])
                ref = ref[:pos] + ref_al
                mut = mut[:pos] + mut_al
            else:
                seq_format = '{:*<' + str(max(ref_len, mut_len)) + '}'
                ref = seq_format.format(ref)
                mut = seq_format.format(mut)

        # compute the raw score of the mutate protein sequence
        raw_score = sum([self._single_score(ref[i], mut[i]) for i in range(len(ref))])

        # compute the minimum and maximum score that are possible for the reference sequence
        max_score = 0
        min_score = 0

        for c in ref_old:
            max_score += self.max[c]
            min_score += self.min[c]

        # return the normalised score
        return (raw_score - min_score) / (max_score - min_score)

    # 'HELPER'
    def _is_symmetrical(self):
        """
        Tests if the substitution matrix is symmetrical
        :return: True if symmetrical, False otherwise
        """
        # get the amino acids in the matrix
        amino_acids = list(self.matrix.keys())

        # iterate over amino acid pairs and see if they are symmetrical
        for aa1 in amino_acids:
            for aa2 in amino_acids:
                if self.matrix[aa1][aa2] != self.matrix[aa2][aa1]:
                    return False
        return True

    def _compute_min_max_scores(self):
        """
        Computes the minimum and maximum score for every amino acid in the matrix.
        """
        for key, row in self.matrix.items():
            scores = row.values()
            self.min[key] = min(scores)
            self.max[key] = max(scores)

    def _single_score(self, aa1, aa2):
        """
        Computes the score for two single amino acids (or gaps).
        :param aa1: amino acid or gap 1
        :param aa2: amino acid or gap 2
        :return: score
        """
        if aa1 in cfg.misc.empty:
            aa1 = '*'
        if aa2 in cfg.misc.empty:
            aa2 = '*'

        return self.matrix[aa1][aa2]

    def _compute_global_alignment(self, ref, mut):
        """
        Computes the global alignment of two protein sequences using the Needleman-Wunsch algorithm with the
        substitution matrix for scoring.
        :param ref: reference protein sequence
        :param mut: mutated protein sequence
        :return: aligned reference and mutated sequence
        """
        # initialise the matrix with 0s
        # matrix[ref][mut]:
        #       *   m   u   t
        #   *
        #   r
        #   e
        #   f
        table = [[0 for _ in range(len(mut) + 1)] for _ in range(len(ref) + 1)]

        # initialise first row
        for i in range(1, len(mut) + 1):
            table[0][i] = table[0][i - 1] + self._single_score(mut[i - 1], '*')

        # initialise first column
        for i in range(1, len(ref) + 1):
            table[i][0] = table[i - 1][0] + self._single_score(ref[i - 1], '*')

        # fill the rest of the table
        for i in range(1, len(ref) + 1):
            for j in range(1, len(mut) + 1):
                match = table[i - 1][j - 1] + self._single_score(ref[i - 1], mut[j - 1])
                deletion = table[i - 1][j] + self._single_score(ref[i - 1], '*')
                insertion = table[i][j - 1] + self._single_score('*', mut[j - 1])
                table[i][j] = max(match, deletion, insertion)

        # compute the aligned sequences
        ref_al = ''
        mut_al = ''

        i = len(ref)
        j = len(mut)

        while i > 0 or j > 0:
            # match
            if i > 0 and j > 0 and table[i][j] == (table[i - 1][j - 1] + self._single_score(ref[i - 1], mut[j - 1])):
                ref_al = ref[i - 1] + ref_al
                mut_al = mut[j - 1] + mut_al
                i -= 1
                j -= 1
            # deletion
            elif i > 0 and table[i][j] == (table[i - 1][j] + self._single_score(ref[i - 1], '*')):
                ref_al = ref[i - 1] + ref_al
                mut_al = '*' + mut_al
                i -= 1
            # insertion
            else:
                ref_al = '*' + ref_al
                mut_al = mut[j - 1] + mut_al
                j -= 1

        return ref_al, mut_al


class Pwm:
    """
    This class represents a position weight matrix of a transcription factor.
    """

    def __init__(self, lt, seq_align):
        """
        Initialises the matrix from the binding site alignment. The DNA sequences in the alignment need to have the same
        length without gaps.
        :param lt: transcription factor locus tag
        :param seq_align: binding site sequence alignment of the transcription factor
        """
        # matrix
        self.A = []
        self.C = []
        self.G = []
        self.T = []
        # conservation
        self.I = []
        # minima and maxima
        self.min = []
        self.max = []
        self.min_score = 0
        self.max_score = 0
        # specs
        self.len = 0
        self.nbr = 0
        self.tf_lt = lt
        # flags
        self.success = True

        first_seq = True

        # process all available sequences of the transcription factor's binding site
        for seq in seq_align:
            seq = seq.upper()

            # skip non-DNA sequences
            if not is_dna_sequence(seq):
                cfg.a_log.pwm_non_dna.add(self.tf_lt)
                continue

            # use the first sequence to initiate the PWM with pseudo counts
            if first_seq:
                first_seq = False
                self.len = len(seq)
                self.A = [0.8] * self.len
                self.C = [0.8] * self.len
                self.G = [0.8] * self.len
                self.T = [0.8] * self.len
                self.I = [0] * self.len
                self.min = [0] * self.len
                self.max = [0] * self.len

            seq_len = len(seq)

            # skip sequences that have a different length
            if seq_len != self.len:
                cfg.a_log.pwm_length.add(self.tf_lt)
                continue

            # count the occurrences of a DNA base at each position
            for i in range(0, seq_len):
                if seq[i] == 'A':
                    self.A[i] += 1
                elif seq[i] == 'C':
                    self.C[i] += 1
                elif seq[i] == 'G':
                    self.G[i] += 1
                elif seq[i] == 'T':
                    self.T[i] += 1
                else:
                    cfg.a_log.pwm_non_dna.add(self.tf_lt)

            self.nbr += 1

        # there were no valid sequences in the alignment
        if self.nbr < cfg.user.as_pwm:
            cfg.a_log.pwm_not_enough.add(self.tf_lt)
            self.success = False
            self.A = []
            self.C = []
            self.G = []
            self.T = []
            self.I = []
            self.len = 0
            self.nbr = 0
            return

        # to PWM
        for i in range(0, self.len):
            # compute the frequency of the bases at each position
            self.A[i] /= self.nbr
            self.C[i] /= self.nbr
            self.G[i] /= self.nbr
            self.T[i] /= self.nbr

            # compute the conservation at each position
            self.I[i] = self.A[i] * log(4 * self.A[i]) + \
                        self.C[i] * log(4 * self.C[i]) + \
                        self.G[i] * log(4 * self.G[i]) + \
                        self.T[i] * log(4 * self.T[i])

            # compute the minimally and maximally possible scores for the position
            self.min[i] = self.I[i] * min(self.A[i], self.C[i], self.G[i], self.T[i])
            self.max[i] = self.I[i] * max(self.A[i], self.C[i], self.G[i], self.T[i])

        # compute the minimally and maximally possible scores for the PWM
        self.min_score = sum(self.min)
        self.max_score = sum(self.max)

    def __getitem__(self, arg):
        """
        Returns the PFM entry at a certain position.
        :param arg: (character, index)
        :return: PFM entry
        """
        c, i = arg
        if i >= self.len:
            return 0
        if c == 'A':
            return self.A[i]
        if c == 'C':
            return self.C[i]
        if c == 'G':
            return self.G[i]
        if c == 'T':
            return self.T[i]
        raise ValueError('{0} is not A, C, G or T.'.format(c))

    def alignment_score(self, seq):
        """
        Computes the normalised score of a transcription factor binding site using the position weight matrix of the
        transcription factor.
        :param seq: transcription factor binding site sequence
        :return: normalised score
        """
        try:
            raw_score = sum([self.I[i] * self[seq[i], i] for i in range(0, self.len)])
        except IndexError:
            cfg.a_log.pwm_score_length.add(self.tf_lt)
            raw_score = 0

        return (raw_score - self.min_score) / (self.max_score - self.min_score)
