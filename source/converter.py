# import standard or third party modules
import re
from time import time
from os import walk, path
from csv import DictReader
from tkinter import messagebox
from collections import defaultdict

# import own modules
from source.configuration import adjust_file_path, adjust_dir_path, cfg
from source.tools import tsv_start, is_dna_sequence, get_file_name, check_tsv_fieldnames, translate, parse_regulation, \
    parse_strand

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2017 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'


class MutationVcfMerger:
    def __init__(self, in_dir, out_path):
        self.in_dir = adjust_dir_path(path.abspath(in_dir))
        self.out_path = adjust_file_path(path.abspath(out_path))
        self.mutations = defaultdict(set)
        # error logging
        self.empty = False
        self.no_files = False
        self.io_errors = set()
        self.field_errors = set()
        self.parser_errors = defaultdict(set)
        self.validate_errors = set()
        self.write_error = False
        # merge the .vcf files iin the input directory
        self.merge()

        if self.no_files:
            messagebox.showerror('Merger Error',
                                 'There are no .vcf files in the input directory.\n\n'
                                 'Check if you entered the correct directory path and if the directory contains '
                                 'valid .vcf files.\n\n'
                                 'Input directory: {0}'.format(self.in_dir))
            return

        if self.empty:
            messagebox.showerror('Merger Error',
                                 'Merging .vcf files failed.\n\n'
                                 'Check if there are mutation VCF files in the '
                                 'directory or its sub-directories, and if these files '
                                 'contain non-empty, tab-separated columns called \'POS\', '
                                 '\'REF\' and \'ALT\'. The order of these columns does not matter '
                                 'and additional columns are not a problem.\n\n'
                                 'Input directory: {0}'.format(self.in_dir))
            return

        # write the result into a .tsv file
        self.write_error = not self.tsv()

        if self.write_error:
            messagebox.showerror('Merger Error',
                                 'Writing the results to the .tsv file failed.\n\n'
                                 'Check if the file is opened in another program and thus not'
                                 'accessible for writing.\n\n'
                                 'Result file: {0}'.format(self.out_path))
            return

    def read(self, file_path):
        """
        Processes a single .vcf file and reads all mutations into the mutation dictionary.
        :param file_path: .vcf file path
        """
        try:
            vcf = open(file_path, 'r')

            reader = DictReader(tsv_start(vcf, '##'), delimiter='\t')

            # empty file
            if not reader.fieldnames:
                self.field_errors.add(file_path)
                return

            # check if all required columns ar ein the .vcf file
            if len({'POS', 'REF', 'ALT'}.intersection(reader.fieldnames)) != 3:
                self.field_errors.add(file_path)
                return

            # process all lines in the file
            for line in reader:
                # check if the position is an integer
                try:
                    pos = int(line['POS'])
                except ValueError:
                    self.parser_errors[file_path].add(line['POS'])
                    continue

                # check if the reference is a DNA base or sequence
                if not is_dna_sequence(line['REF']):
                    self.parser_errors[file_path].add(line['REF'])
                    continue

                # check if the alternative is DNA base or sequence
                if not is_dna_sequence(line['ALT']):
                    self.parser_errors[file_path].add(line['ALT'])
                    continue

                # add the mutation
                self.mutations[line['POS']].add((line['REF'], line['ALT']))

            vcf.close()

        except IOError:
            self.io_errors.add(file_path)

    def validate(self):
        """
        Makes sure that for each mutation position, all references start with the same DNA bases.
        """
        # make sure that for each mutation position, all references start with the same DNA bases
        for key, value in self.mutations.items():
            refs = [val[0] for val in value]
            longest_ref = max(refs, key=len)

            for ref in refs:
                if not longest_ref.startswith(ref):
                    self.validate_errors.add(key)
                    break

        # remove faulty mutation positions
        for pos in self.validate_errors:
            del(self.mutations[pos])

    def merge(self):
        """
        Merges the mutations in all .vcf files in the input directory.
        """
        # get all .vcf files in the input directory
        vcf_files = set()

        for root, dirs, files in walk(self.in_dir):
            for name in files:
                _, _, suffix = get_file_name(name)
                if suffix == '.vcf':
                    vcf_files.add(adjust_file_path(path.join(root, name)))

        # no .vcf file
        if not vcf_files:
            self.no_files = True
            return

        # read each .vcf file
        for vcf in vcf_files:
            self.read(vcf)

        # make sure the mutations are valid
        self.validate()

        if len(self.io_errors.union(self.field_errors)) == len(vcf_files) or not self.mutations:
            self.empty = True
            return

    def tsv(self):
        """
        Writes the result into a .tsv file.
        :return: True if successful, False otherwise
        """
        try:
            file = open(self.out_path, 'w')
            fmt = '{0}\n'

            # header
            file.write(fmt.format('\t'.join([cfg.mut.pos, cfg.mut.ref, cfg.mut.alt])))

            # write all mutations ordered by their position
            for key, value in sorted(self.mutations.items(), key=lambda x: int(x[0])):
                for ref, alt in value:
                    file.write(fmt.format('\t'.join([key, ref, alt])))

            file.close()

            return True

        except IOError:
            return False


class Parser:
    """
    Some basic functionality for parsers.
    """
    def __init__(self, in_path, name='', type='Converter'):
        self.name = name
        self.type = type
        self.in_path = adjust_file_path(path.abspath(in_path))
        # error flags
        self.io_error = False
        self.empty = False

    def empty_msg(self):
        """
        Shows an error message when there were no valid entries in the input file.
        """
        self.empty = True
        messagebox.showerror('Converter Error',
                             'There are no valid entries in the {0} file.\n\n'
                             'Check if the file is empty or the format is incorrect.\n\n'
                             'Input file: {1}'.format(self.name, self.in_path))

    def io_msg(self):
        """
        Shows an error message when an IOError occurred while reading the input file.
        """
        self.io_error = True
        messagebox.showerror('{0} Error'.format(self.type),
                             'Reading the {0} file failed.\n\n'
                             'Check if the file is opened in another program and thus not'
                             'accessible for writing.\n\n'
                             'Input file: {1}'.format(self.name, self.in_path))

    def validate(self):
        """
        Checks for errors in the parsed input.
        """
        pass

    def process_lines(self, file):
        """
        Processes the file content.
        """
        pass

    def read(self):
        """
        Reads and processes the input file.
        """
        try:
            file = open(self.in_path, 'r')
            self.process_lines(file)
            self.validate()
            file.close()
        except IOError:
            self.io_msg()


class UniProtEntry:
    """
    Processes and stores a single UniProt entry.
    """
    def __init__(self, gn, ft):
        # input data
        self.gn = gn
        self.ft = ft
        # parsed data
        self.names = set()
        self.domains = set()
        # error flag
        self.success = True
        # parse entry
        self.parse_gn()
        self.parse_ft()
        self.validate()

    def parse_gn(self):
        """
        Parses the GN line, which contains gene name, locus tags and/or ORF names.
        """
        # join all GN lines into a single line
        gn = ' '.join(self.gn).replace('\n', '')
        # a set of protein features can be associated with several genes. the names/locus tags
        # belonging to different genes are separated by a 'GN   and' line.
        pairs = gn.split('GN   and')

        for pair in pairs:
            pair = re.sub(r'[{].*?[}]', '', pair.replace('GN', '').strip().strip(';'))
            words = pair.split('; ')

            names = set()
            lts = set()

            for word in words:
                word = word.strip().strip(';').strip()
                # parse gene name
                if word.startswith('Name='):
                    names.add(word[5:])
                # parse alternative gene names
                elif word.startswith('Synonyms='):
                    names = names.union(set(word[9:].split(', ')))
                # parse locus tags
                elif word.startswith('OrderedLocusNames='):
                    lts = lts.union(set(word[18:].split(', ')))
                # parse open reading frame names
                elif word.startswith('ORFNames='):
                    lts = lts.union(set(word[9:].split(', ')))

            # create a tuple for each pair of locus tag and name
            for lt in lts:
                lt = lt.strip()
                for name in names:
                    self.names.add((lt, name.strip()))
                if not names:
                    self.names.add((lt, ''))

            # if there is no locus tag given, only add the gene name(s)
            if not lts:
                for name in names:
                    self.names.add(('', name.strip()))

    def process_domain(self, t, loc, ds):
        """
        Processes the information of a single protein domain.
        :param t: domain type
        :param loc: location
        :param ds: list of description lines
        """
        # don't add the domain if location or type are not given
        if not t or not loc:
            return

        # extract start and end position from the location
        loc = loc.split('..')
        if len(loc) == 2:
            s, e = loc
        # sometimes there is only the start position
        elif len(loc) == 1:
            s = loc[0]
            e = s
        else:
            return

        # don't add the domain if start or end are not numbers
        if not s.isdigit() or not e.isdigit():
            return

        # a single description lines can contain several descriptions as follows: /<qualifier>="<value>"
        # this creates an adjusted list where each element is a single description
        ds_adjusted = []
        for d in ds:
            if not d:
                continue
            if '="' in d:
                ds_adjusted.append(d.split('="')[1].strip('"'))
            elif d:
                ds_adjusted[-1] += ' ' + d.strip('"')

        # final list of domain descriptions
        desc = []

        for d in ds_adjusted:
            # ignore empty descriptions
            d = d.strip().replace('|', '-')
            if d:
                desc.append(d)

        # add the protein domain
        self.domains.add((s, e, t.lower(), cfg.misc.in_sep.join(desc)))

    def parse_ft(self):
        """
        Parses the FT lines, which contain protein domain information.
        """
        # do nothing if no lines are given
        if not self.ft:
            return

        domain_type = ''
        location = ''
        desc = []

        for line in self.ft:
            # parse the columns
            t_type = line[5:21].strip()

            # new feature
            if t_type:
                # add the previous feature to the domain list
                self.process_domain(domain_type, location, desc)
                # set the type, start and end of the new feature
                domain_type = t_type
                location = line[21:].strip('?<>').strip()
                desc = []
            # previous feature is continued
            else:
                desc.append(line[21:].strip())

        # add the last feature to the domain list
        self.process_domain(domain_type, location, desc)

    def validate(self):
        """
        Tests if the entry is valid.
        """
        self.names.discard(('', ''))

        # an entry needs to contain a least one name
        if not self.names:
            self.success = False
        # an entry needs to contain at least one protein domain
        if not self.domains:
            self.success = False

    def get(self):
        """
        :return: set of tuples with all name - domain combinations
        """
        ret = set()

        for name in self.names:
            if name == ('', ''):
                continue
            for dom in self.domains:
                if dom == ('', '', '', ''):
                    continue
                ret.add(name + dom)

        return ret


class UniProtConverter(Parser):
    def __init__(self, in_path, out_path):
        Parser.__init__(self, in_path, 'UniProt input')
        self.in_path = adjust_file_path(path.abspath(in_path))
        self.out_path = adjust_file_path(path.abspath(out_path))
        # domain entries linked to (locus tag, gene name) keys
        self.entries = set()
        # read the input file and output it as a .tsv file
        self.read()
        # stop if errors occurred while reading the file
        if self.io_error or self.empty:
            return
        # write the result file
        self.tsv()

    def add_entry(self, gn, ft):
        """
        Adds an entry to the entry dictionary.
        :param gn: list of gene name lines
        :param ft: list of feature lines
        """
        entry = UniProtEntry(gn, ft)
        # add the protein domains to all associated genes if the entry is valid
        if entry.success:
            self.entries = self.entries.union(entry.get())

    def validate(self):
        """
        Checks if the entries are valid.
        """
        if not self.entries:
            self.empty_msg()

    def process_lines(self, file):
        """
        Processes the input file.
        """
        # list of gene name lines, list of feature lines
        gn = []
        ft = []

        for line in file:
            # end of an entry: add the entry and reset the gene names and features
            if line.startswith('//'):
                self.add_entry(gn, ft)

                gn = []
                ft = []
            # gene name line
            elif line.startswith('GN'):
                gn.append(line)
            # feature line
            elif line.startswith('FT'):
                ft.append(line)

        # add the last entry
        if gn or ft:
            self.add_entry(gn, ft)

    def tsv(self):
        """
        Write the protein domains to the output .tsv file.
        """
        try:
            file = open(self.out_path, 'w')
            # write the file header
            file.write('\t'.join([cfg.prot.lt, cfg.prot.name, cfg.prot.start,
                                  cfg.prot.end, cfg.prot.type, cfg.prot.desc]) + '\n')

            # sort the entries by locus tag and gene name and write them into the file
            for d in sorted(self.entries):
                file.write('\t'.join([d[0], d[1], d[2], d[3], d[4], d[5]]) + '\n')

            file.close()
        except IOError:
            messagebox.showerror('Converter Error',
                                 'Writing the results to the .tsv file failed.\n\n'
                                 'Check if the file is opened in another program and thus not'
                                 'accessible for writing.\n\n'
                                 'Result file: {0}'.format(self.out_path))


class PatricConverter(Parser):
    """
    Parses, processes and converts a PATRIC antibiotic resistance file to .tsv
    """
    def __init__(self, in_path, out_path):
        Parser.__init__(self, in_path, 'PATRIC input')
        self.in_path = adjust_file_path(path.abspath(in_path))
        self.out_path = adjust_file_path(path.abspath(out_path))
        # domain entries linked to (locus tag, gene name) keys
        self.entries = defaultdict(set)
        # error flags
        self.field_error = False
        # read the input file and output it as a .tsv file
        self.read()
        # check if an error occurred while reading the input file
        if self.io_error or self.field_error or self.empty:
            return
        # write the result file
        self.tsv()

    def validate(self):
        """
        Checks if the entries are valid.
        """
        if not self.entries:
            self.empty_msg()

    def process_lines(self, file):
        """
        Processes the input file.
        """
        reader = DictReader(tsv_start(file), delimiter=',')

        # check if the file contains the necessary columns
        if not check_tsv_fieldnames('PATRIC antibiotic resistance', self.in_path,
                                    ['Property', 'RefSeq Locus Tag', 'Gene'], reader.fieldnames):
            self.field_error = True
            return

        for line in reader:
            # skip non antibiotic resistance genes
            if line['Property'].strip().strip('"').lower() != 'antibiotic resistance':
                continue
            # add the locus tag and gene name
            self.entries[line['RefSeq Locus Tag'].strip().strip('"')].add(line['Gene'].strip().strip('"').strip())

    def tsv(self):
        """
        Write the antibiotic resistance genes to the output .tsv file.
        """
        try:
            file = open(self.out_path, 'w')

            # write the header
            file.write('\t'.join([cfg.int.lt, cfg.int.name, cfg.int.desc]) + '\n')

            # sort the genes by locus tag
            for key, val in sorted(self.entries.items(), key=lambda x: x[0]):
                # skip empty locus tag
                if not key:
                    continue
                # remove empty gene names
                val.discard('')
                # if a set of gene names becomes empty, add the empty gene name
                if not val:
                    val = {''}
                # write the line
                for name in val:
                    file.write('\t'.join([key, name, 'antibiotic resistance']) + '\n')

            file.close()
        except IOError:
            messagebox.showerror('Converter Error',
                                 'Writing the result file failed.\n\n'
                                 'Check if the file is opened in another program and thus not'
                                 'accessible for writing.\n\n'
                                 'Result file: {0}'.format(self.out_path))


class Gene:
    """
    Stores a gene for the RegulonDB converter.
    """
    def __init__(self, name, desc, dna, start, end, strand):
        self.name = re.split('\W+', name.strip())[0]
        self.desc = desc

        # make sure the start and end positions are integers and => 0
        if start.isdigit():
            start = max(int(start), 0)
        else:
            start = 0
        if end.isdigit():
            end = max(int(end), 0)
        else:
            end = 0

        self.g_start = min(start, end)
        self.g_end = max(start, end)

        # process the strand
        self.strand = parse_strand(strand)

        # check the DNA and translate it
        if is_dna_sequence(dna):
            self.dna = dna
            self.protein = translate(dna)
        else:
            self.dna = ''
            self.protein = ''

        # promoter id, start, and start
        self.prom_id = ''
        self.prom_start = 0
        self.prom_end = 0

        # operon id and genes
        self.operon_id = ''
        self.operon = []


class Error:
    """
    Stores and output errors messages for the RegulonDB converter.
    """
    def __init__(self):
        # gene errors
        self.gene_incomplete = set()
        self.gene_duplicates = set()
        self.gene_wrong_length = set()
        self.gene_multiple_operons = defaultdict(set)
        self.gene_multiple_promoter = defaultdict(set)

        # operon errors
        self.operon_unknown_genes = dict()
        self.operon_duplicates = set()
        self.operon_singletons = set()
        self.operon_incomplete = set()

        # transcription unit errors
        self.tu_duplicates = set()
        self.tu_unknown_genes = dict()
        self.tu_incomplete = set()

        # regulation errors
        self.reg_incomplete = set()
        self.reg_missing_tus = set()
        self.reg_missing_tfs = set()
        self.reg_missing_operons = set()
        self.reg_missing_genes = set()
        self.reg_missing_reg = set()

        # transcription factor MSA errors
        self.msa_incomplete = set()
        self.msa_wrong_length = set()
        self.msa_duplicates = set()
        self.msa_unknown_genes = set()
        self.msa_not_dna = set()

        # transcription factor binding site errors
        self.tfbs_incomplete = set()
        self.tfbs_missing_tus = set()
        self.tfbs_not_dna = set()
        self.tfbs_missing_tfs = set()
        self.tfbs_wrong_pos = set()

    def lst(self, desc, lst):
        """
        Help function for formatting lists.
        """
        return '{0}:\n-  {1}'.format(desc, '\n-  '.join(sorted(lst)))
    
    def sections_condensed(self):
        """
        Returns a list of statistics and errors grouped in sections.
        """
        sections = []

        # genes
        header = '[RegulonDB Genes]'
        msg = []

        if self.gene_incomplete:
            msg.append('{} gene entries were incomplete.'.format(len(self.gene_incomplete)))

        if self.gene_duplicates:
            msg.append('{} gene names occurred more than once.'.format(len(self.gene_duplicates)))

        if self.gene_wrong_length:
            msg.append('The start and end of {} gene entries did not match the '
                       'gene length.'.format(len(self.gene_wrong_length)))

        if self.gene_multiple_operons:
            msg.append('{} genes had more than one operon.'
                       ''.format(sum([len(val) for val in self.gene_multiple_operons.values()])))

        if self.gene_multiple_promoter:
            msg.append('{} genes had more than one promoter.'
                       ''.format(sum([len(val) for val in self.gene_multiple_promoter.values()])))

        if msg:
            sections.append('{0}\n{1}'.format(header, '\n'.join(msg)))

        # operons
        header = '[RegulonDB Operons]'
        msg = []

        if self.operon_incomplete:
            msg.append('{} operon entries were incomplete.'.format(len(self.operon_incomplete)))

        if self.operon_duplicates:
            msg.append('{}  operon names occurred more than once.'.format(len(self.operon_duplicates)))

        if self.operon_singletons:
            msg.append('{} operons were singletons.'.format(len(self.operon_singletons)))

        if self.operon_unknown_genes:
            msg.append('{} operons contained unknown genes.'
                       ''.format(sum([len(val) for val in self.operon_unknown_genes.values()])))

        if msg:
            sections.append('{0}\n{1}'.format(header, '\n'.join(msg)))

        # TUs
        header = '[RegulonDB Transcription Units (TUs)]'
        msg = []

        if self.tu_incomplete:
            msg.append('{} TU entries were incomplete.'.format(len(self.tu_incomplete)))

        if self.tu_duplicates:
            msg.append('{} TU names occurred more than once.'.format(len(self.tu_duplicates)))

        if self.tu_unknown_genes:
            msg.append('{} TUs were not associated with genes.'.format(len(self.tu_unknown_genes)))

        if msg:
            sections.append('{0}\n{1}'.format(header, '\n'.join(msg)))

        # regulation
        header = '[RegulonDB Regulation]'
        msg = []

        if self.reg_incomplete:
            msg.append('{} regulation entries were incomplete.'.format(len(self.reg_incomplete)))

        if self.reg_missing_tfs:
            msg.append('{} TF names did not occur.'.format(len(self.reg_missing_tfs)))

        if self.reg_missing_genes:
            msg.append('{} gene names did not occur.'.format(len(self.reg_missing_genes)))

        if self.reg_missing_operons:
            msg.append('{} operon names did not occur.'.format(len(self.reg_missing_operons)))

        if self.reg_missing_tus:
            msg.append('{} TU names did not occur.'.format(len(self.reg_missing_tus)))

        if self.reg_missing_reg:
            msg.append('For {} TFs regulation was missing.'.format(len(self.reg_missing_reg)))

        if msg:
            sections.append('{0}\n{1}'.format(header, '\n'.join(msg)))

        # MSA
        header = '[RegulonDB MSA]'
        msg = []

        if self.msa_incomplete:
            msg.append('{} MSA entries were incomplete.'.format(len(self.msa_incomplete)))

        if self.msa_unknown_genes:
            msg.append('{} TFs were unknown.'.format(len(self.msa_unknown_genes)))

        if self.msa_duplicates:
            msg.append('{} TFs occurred more than once.'.format(len(self.msa_duplicates)))

        if self.msa_wrong_length:
            msg.append('{} MSA entries had inconsistent length.'.format(len(self.msa_wrong_length)))

        if self.msa_not_dna:
            msg.append('{} MSA entries contained non-DNA sequences.'.format(len(self.msa_not_dna)))

        if msg:
            sections.append('{0}\n{1}'.format(header, '\n'.join(msg)))

        # TFBS
        header = '[RegulonDB TFBS]'
        msg = []

        if self.tfbs_incomplete:
            msg.append('{}  TFBS entries were incomplete.'.format(len(self.tfbs_incomplete)))

        if self.tfbs_not_dna:
            msg.append('{} TFBS entries had a non-DNA sequence.'.format(len(self.tfbs_not_dna)))

        if self.tfbs_wrong_pos:
            msg.append('{} TFBS entries had a wrong position.'.format(len(self.tfbs_wrong_pos)))

        if self.tfbs_missing_tfs:
            msg.append('{} TFs were unknown.'.format(len(self.tfbs_missing_tfs)))

        if self.tfbs_missing_tus:
            msg.append('{} TU names were unknown.'.format(len(self.tfbs_missing_tus)))

        if msg:
            sections.append('{0}\n{1}'.format(header, '\n'.join(msg)))

        return sections

    def sections_verbose(self):
        """
        Returns a list of statistics and errors grouped in sections.
        """
        sections = []

        # genes
        header = '[RegulonDB Genes]'
        msg = []

        if self.gene_incomplete:
            msg.append(self.lst('The following gene entries were incomplete', self.gene_incomplete))

        if self.gene_duplicates:
            msg.append(self.lst('The following gene names occurred more than once', self.gene_duplicates))

        if self.gene_wrong_length:
            msg.append(self.lst('The start and end of the following gene entries did not match the gene length',
                                self.gene_wrong_length))

        if self.gene_multiple_operons:
            lst = ['Gene name: {0}, operons: {1}'.format(key, val)
                   for key, val in sorted(self.gene_multiple_operons.items())]
            msg.append('The following genes had more than one operon:\n{}'.format('\n'.join(lst)))

        if self.gene_multiple_promoter:
            lst = ['Gene name: {0}, promoter: {1}'.format(key, val)
                   for key, val in sorted(self.gene_multiple_promoter.items())]
            msg.append('The following genes had more than one promoter:\n{}'.format('\n'.join(lst)))

        if msg:
            sections.append('{0}\n{1}'.format(header, '\n\n'.join(msg)))

        # operons
        header = '[RegulonDB Operons]'
        msg = []

        if self.operon_incomplete:
            msg.append(self.lst('The following operon entries were incomplete', self.operon_incomplete))

        if self.operon_duplicates:
            msg.append(self.lst('The following operon names occurred more than once', self.operon_duplicates))

        if self.operon_singletons:
            msg.append(self.lst('The following operons were singletons', self.operon_singletons))

        if self.operon_unknown_genes:
            lst = ['Operon ID: {0}, unknown genes: {1}'.format(key, val)
                   for key, val in sorted(self.operon_unknown_genes.items())]
            msg.append('The following operons contained unknown genes:\n{}'.format('\n'.join(lst)))

        if msg:
            sections.append('{0}\n{1}'.format(header, '\n\n'.join(msg)))

        # TUs
        header = '[RegulonDB Transcription Units (TUs)]'
        msg = []

        if self.tu_incomplete:
            msg.append(self.lst('The following TU entries were incomplete', self.tu_incomplete))

        if self.tu_duplicates:
            msg.append(self.lst('The following TU names occurred more than once', self.tu_duplicates))

        if self.tu_unknown_genes:
            msg.append(self.lst('The following TUs were not associated with genes', self.tu_unknown_genes))

        if msg:
            sections.append('{0}\n{1}'.format(header, '\n\n'.join(msg)))

        # regulation
        header = '[RegulonDB Regulation]'
        msg = []

        if self.reg_incomplete:
            msg.append(self.lst('The following regulation entries were incomplete', self.reg_incomplete))

        if self.reg_missing_tfs:
            msg.append(self.lst('The following TF names did not occur', self.reg_missing_tfs))

        if self.reg_missing_genes:
            msg.append(self.lst('The following gene names did not occur', self.reg_missing_genes))

        if self.reg_missing_operons:
            msg.append(self.lst('The following operon names did not occur', self.reg_missing_operons))

        if self.reg_missing_tus:
            msg.append(self.lst('The following TU names did not occur', self.reg_missing_tus))

        if self.reg_missing_reg:
            msg.append(self.lst('For the following TFs regulation was missing', self.reg_missing_reg))

        if msg:
            sections.append('{0}\n{1}'.format(header, '\n\n'.join(msg)))

        # MSA
        header = '[RegulonDB MSA]'
        msg = []

        if self.msa_incomplete:
            msg.append(self.lst('The following MSA entries were incomplete', self.msa_incomplete))

        if self.msa_unknown_genes:
            msg.append(self.lst('The following TFs were unknown', self.msa_unknown_genes))

        if self.msa_duplicates:
            msg.append(self.lst('The following TFs occurred more than once', self.msa_duplicates))

        if self.msa_wrong_length:
            msg.append(self.lst('The following MSA entries had inconsistent length', self.msa_wrong_length))

        if self.msa_not_dna:
            msg.append(self.lst('The following MSA entries contained non-DNA sequences', self.msa_not_dna))

        if msg:
            sections.append('{0}\n{1}'.format(header, '\n\n'.join(msg)))

        # TFBS
        header = '[RegulonDB TFBS]'
        msg = []

        if self.tfbs_incomplete:
            msg.append(self.lst('The following TFBS entries were incomplete', self.tfbs_incomplete))

        if self.tfbs_not_dna:
            msg.append(self.lst('The following TFBS entries had a non-DNA sequence', self.tfbs_not_dna))

        if self.tfbs_wrong_pos:
            lst = [str(x) for x in self.tfbs_wrong_pos]
            msg.append(self.lst('The following TFBS entries had a wrong position', lst))

        if self.tfbs_missing_tfs:
            msg.append(self.lst('The following TFs were unknown', self.tfbs_missing_tfs))

        if self.tfbs_missing_tus:
            msg.append(self.lst('The following TU names were unknown', self.tfbs_missing_tus))

        if msg:
            sections.append('{0}\n{1}'.format(header, '\n\n'.join(msg)))

        return sections


class Index:
    """
    Stores and maintains the parsed RegulonDB input.
    """
    def __init__(self):
        # error collection
        self.errors = Error()
        # data
        self.genes = dict()
        self.operons = dict()
        self.tus = dict()
        self.reg = defaultdict(set)
        self.tfbs = defaultdict(set)
        self.msa = dict()

    @staticmethod
    def tf_name(tf):
        """
        Adjust transcription factor names.
        """
        # NAME -> name
        if tf.isupper():
            return tf.lower()
        # NamE -> namE
        else:
            return tf[:1].lower() + tf[1:]

    def add_gene(self, gene):
        """
        Checks if genome information is valid and if so, adds it to the data set.
        """
        # non-existent / invalid gene name
        if not gene.name or gene.name == 'Phantom Gene':
            return False

        # invalid gene start and end
        if gene.g_start <= 0 or gene.g_end <= 0:
            gene.g_start = 0
            gene.g_end = 0
            self.errors.gene_incomplete.add(gene.name)

        # invalid DNA, strand or protein
        if not gene.dna or not gene.protein or not gene.strand:
            self.errors.gene_incomplete.add(gene.name)

        # gene start and end does not match DNA length
        if len(gene.dna) != (abs(gene.g_start - gene.g_end) + 1):
            self.errors.gene_wrong_length.add(gene.name)

        # duplicate
        if gene.name in self.genes:
            self.errors.gene_duplicates.add(gene.name)
            return False

        self.genes[gene.name] = gene
        return True

    def contains_unknown_genes(self, lst):
        """
        Tests if a list contains gene names that are not in the data set and returns a list of all unknown gene names.
        """
        unknown = set()
        for l in lst:
            if l not in self.genes:
                unknown.add(l)
        return unknown

    def add_operon_to_gene(self, id, operon):
        """
        Add the given operon to all genes in the operon.
        """
        # first check if the genes are already associated with another operon
        add = True
        for gn in operon:
            gene = self.genes[gn]
            # check if the gene is already associated with another operon
            if gene.operon:
                self.errors.gene_multiple_operons[gn].add(id)
                self.errors.gene_multiple_operons[gn].add(gene.operon_id)
                add = False
        # do not add the operon to genes if any of the genes in the operon are already associated with another operon
        if not add:
            return
        # add the operon to all genes in the operon
        for gn in operon:
            gene = self.genes[gn]
            gene.operon = operon
            gene.operon_id = id

    def add_operon(self, id, genes):
        """
        Checks if the given operon information is valid and if so, adds the operon to the data set.
        :param id: operon ID
        :param genes: list of genes in the operon
        :return: True if the operon was successfully added to the data set, False otherwise
        """
        # the operon needs an ID
        if not id:
            return False

        err = False
        # the operon needs to contain at least one gene
        if not genes:
            self.errors.operon_incomplete.add(id)
            err = True

        # duplicate: the operon with the same ID is already in the data sets
        if id in self.operons:
            self.errors.operon_duplicates.add(id)
            err = True

        # the operon can only contain genes that are in the data set
        unknown = self.contains_unknown_genes(genes)
        if unknown:
            self.errors.operon_unknown_genes[id] = unknown
            err = True

        # do not add the operon if errors occurred
        if err:
            return False

        self.operons[id] = genes

        # check if the operon consists of only one gene
        if len(genes) > 2:
            self.add_operon_to_gene(id, genes)
        else:
            self.errors.operon_singletons.add(id)

        return True

    def add_tu(self, tu_id, genes):
        """
        Checks if the given transcription unit information is valid and if so, adds the TU to the data set.
        :param tu_id: transcription unit TU
        :param genes: genes in the TU
        :return: True if the TU was successfully added to the data set, False otherwise
        """
        # the TU needs to have an ID
        if not tu_id:
            return False

        err = False
        # the TU needs to contain at least one gene
        if not genes:
            self.errors.tu_incomplete.add(tu_id)
            err = True

        # all genes in the TU needs to be in the data set
        unknown = self.contains_unknown_genes(genes)
        if unknown:
            self.errors.tu_unknown_genes[tu_id] = unknown
            err = True

        # duplicate: a different TU with that ID is already in the data set
        if tu_id in self.tus:
            if genes != self.tus[tu_id]:
                self.errors.tu_duplicates.add(tu_id)
                err = True

        # do not add the TU since errors occurred
        if err:
            return False

        # add the TU
        self.tus[tu_id] = genes

        return True

    def add_tu_reg(self, tf, tu_id, reg):
        """
        Add TF-TU regulation information
        :param tf: transcription factor name
        :param tu_id: transcription unit ID (regulated)
        :param reg: regulation description
        :return: True if successful, False otherwise
        """
        # the transcription factor and transcription unit ID need to be given
        if not tu_id or not tf:
            return False

        err = False
        # make sure the regulatory description is given
        if not reg:
            self.errors.reg_incomplete.add((tf, tu_id))
            err = True

        # the transcription unit does not occur in the dataa set
        if tu_id not in self.tus:
            self.errors.reg_missing_tus.add(tu_id)
            err = True

        # adjust the transcription factor name
        tf = self.tf_name(tf)

        # the transcription factor is not in the data set
        if tf not in self.genes:
            self.errors.reg_missing_tfs.add(tf)
            err = True

        # do not add the regulation to the data set since errors occurred
        if err:
            return False

        # add the regulation for all genes in the transcription unit
        for gn in self.tus[tu_id]:
            self.reg[tf].add((gn, parse_regulation(reg)))

        return True

    def add_operon_reg(self, tf, op_id, genes, reg):
        """
        Add TF-operon regulation information
        :param tf: transcription factor name
        :param op_id: operon ID (regulated)
        :param genes: regulated genes in the operon
        :param reg: regulation description
        :return: True if successful, False otherwise
        """
        # the transcription factor and operon ID need to be given
        if not op_id or not tf:
            return False

        err = False
        # make sure the regulatory description and genes are given
        if not reg or not genes:
            self.errors.reg_incomplete.add((tf, op_id))
            err = True

        # the operon is not in the data set
        if op_id not in self.operons:
            self.errors.reg_missing_operons.add(op_id)
            err = True

        # adjust the transcription factor name
        tf = self.tf_name(tf)

        # the transcription factor is not in the data set
        if tf not in self.genes:
            self.errors.reg_missing_tfs.add(tf)
            err = True

        # do not add the regulation since errors occurred
        if err:
            return False

        # add regulation for all genes in the operon that are also given in the regulation file
        for gn in self.operons[op_id]:
            if gn in genes:
                self.reg[tf].add((gn, parse_regulation(reg)))

        return True

    def add_gene_reg(self, tf, gn, reg):
        """
        Add TF-TF or TF-gene regulation information.
        :param tf: transcription factor name
        :param gn: gene name (regulated)
        :param reg: regulation description
        :return: True if successful, False otherwise
        """
        # the gene name and transcription factor name needs to be given
        if not gn or not tf:
            return False

        err = False
        # the regulation description needs to be given
        if not reg:
            self.errors.reg_incomplete.add((tf, gn))
            err = True

        # the regulated gene is not in the data set
        if gn not in self.genes:
            self.errors.reg_missing_genes.add(gn)
            err = True

        # adjust the transcription factor name
        tf = self.tf_name(tf)

        # the transcription factor is not in the data set
        if tf not in self.genes:
            self.errors.reg_missing_tfs.add(tf)
            err = True

        # do not add the regulation since errors occurred
        if err:
            return False

        # add the regulation
        self.reg[tf].add((gn, parse_regulation(reg)))

        return True

    def add_msa(self, tf, align):
        """
        Check if a binding site alignment for a transcription factor is valid, and if so, add it to the data set.
        :param tf: transcription factor name
        :param align: aligned binding site sequences
        :return: True if successful, False otherwise
        """
        # the transcription factor needs to be given
        if not tf:
            return False

        # adjust the transcription factor name
        tf = self.tf_name(tf)

        err = False
        # the alignment needs to be given
        if not align:
            self.errors.msa_incomplete.add(tf)
            err = True

        # all sequences in the alignment needs to have the same length
        if not all(len(x) == len(align[0]) for x in align):
            self.errors.msa_wrong_length.add(tf)
            err = True

        # all sequences in the alignment need to be DNA sequences
        if not all(is_dna_sequence(x) for x in align):
            self.errors.msa_not_dna.add(tf)
            err = True

        # the transcription factor is not in the data set
        if tf not in self.genes:
            self.errors.msa_unknown_genes.add(tf)
            err = True

        # duplicate: there is already an alignment for the transcription factor in the data set
        if tf in self.msa:
            self.errors.msa_duplicates.add(tf)
            err = True

        # do not add the alignment since errors occurred
        if err:
            return False

        # add the alignment
        self.msa[tf] = align

        return True

    def get_tfbs_genes(self, tu_id, prom_id):
        """
        Return a set of all gene names associated with a transcription factor binding site.
        :param tu_id: transcription unit ID
        :param prom_id: promoter ID
        :return: set of gene names
        """
        genes = set()

        # add all genes from the TU to the set if the TU is in the data set
        if tu_id not in self.tus:
            self.errors.tfbs_missing_tus.add(tu_id)
        else:
            genes.update(self.tus[tu_id])

        return genes

    def add_tfbs(self, tf, tu_id, prom_id, start, end, strand, dna):
        """
        Check if transcription factor binding site information is valid, and if so, add it to the data set.
        :param tf: transcription factor name
        :param tu_id: transcription unit ID (regulated)
        :param prom_id: promoter ID (regulated)
        :param start: genome start position on the forward strand
        :param end: genome end position on the forward strand
        :param strand: DNA strand o respective strand
        :param dna: binding site DNA sequence
        :return: True of successful, False otherwise
        """
        # the transcription factor needs to be given
        if not tf:
            return False

        # adjust the transcription factor name
        tf = self.tf_name(tf)

        err = False
        # check if the entry is complete
        if (not tu_id and not prom_id) or not start or not end or not strand or not dna:
            self.errors.tfbs_incomplete.add(tf)
            err = True
        # the sequence needs to be a DNA sequence
        if not is_dna_sequence(dna):
            self.errors.tfbs_not_dna.add(tf)
            err = True
        # the start and end positions need to be > 0 (= 0 is invalid/unknown)
        if start <= 0 or end <= 0:
            self.errors.tfbs_incomplete.add(tf)
            err = True

        # the transcription factor is not in the data set
        if tf not in self.genes:
            self.errors.tfbs_missing_tfs.add(tf)
            err = True

        # wrong start and end position
        if not len(dna) == abs(start - end) + 1 and start and end and dna:
            self.errors.tfbs_wrong_pos.add((tf, start, end, len(dna)))

        # get a list of all genes associated with the transcription unit and promoter ID
        genes = self.get_tfbs_genes(tu_id, prom_id)

        # only add the binding site if there are genes associated with it
        if not genes or err:
            return False

        # add the binding site to all associated genes
        for gn in genes:
            self.tfbs[tf].add((gn, start, end, dna))

        return True


class RegulonDBConverter:
    def __init__(self, out_dir, gene, operons, tus, operon_reg, gene_reg, tf_reg, tu_reg, tfbs, msa, skip,
                 status):
        # setup
        self.out_dir = adjust_dir_path(path.abspath(out_dir))
        self.index = Index()
        self.skip = skip

        # status
        status.name('RegulonDB Converter')

        # parse basic input files
        start = time()
        status.running(0, 'parsing basic input files')
        self.read_file('RegulonDB gene', gene, self.parse_gene)
        self.read_file('RegulonDB operon', operons, self.parse_operon)
        self.read_file('RegulonDB TU', tus, self.parse_tus)
        status.time(0, 'parsing basic input files', time() - start)

        # parse regulation files
        start2 = time()
        status.running(1, 'parsing regulation input files')
        self.read_file('RegulonDB TF-gene', gene_reg, self.parse_gene_reg)
        self.read_file('RegulonDB TF-TF', tf_reg, self.parse_gene_reg)
        self.read_file('RegulonDB TF-operon', operon_reg, self.parse_operon_reg)
        self.read_file('RegulonDB TF-TU', tu_reg, self.parse_tu_reg)
        status.time(1, 'parsing regulation input files', time() - start2)

        # parse TFBS files
        start2 = time()
        status.running(2, 'parsing TFBS input files')
        self.read_file('RegulonDB TFBS', tfbs, self.parse_tfbs)
        self.read_file('RegulonDB TF MSA', msa, self.parse_msa)
        status.time(2, 'parsing TFBS input files', time() - start2)

        # write files
        start2 = time()
        status.running(3, 'writing converted files')
        self.write_file('gene result', 'genes.tsv', self.gene_tsv)
        self.write_file('regulation result', 'regulation.tsv', self.reg_tsv)
        self.write_file('TF MSA result', 'TF_MSA.fasta', self.msa_text)
        self.write_file('TFBS result', 'TFBS.tsv', self.tfbs_tsv)
        status.time(3, 'writing converted files', time() - start2)

        # write log
        start2 = time()
        status.running(4, 'writing log')
        self.write_file('RegulonDB converter log', 'log.txt', self.log)
        status.time(4, 'wwriting log', time() - start2)

        # total running time
        status.total_time(time() - start)

    @staticmethod
    def in_msg(name, in_path):
        """
        Show a file reading error message.
        """
        messagebox.showerror('Converter Error',
                             'Reading the {0} file failed.\n\n'
                             'Check if the file is opened in another program and thus not'
                             'accessible for writing.\n\n'
                             'Input file: {1}'.format(name, in_path))

    @staticmethod
    def empty_msg(name, in_path):
        """
        Show a no valid entries error message.
        """
        messagebox.showerror('Converter Error',
                             'There are no valid entries in the {0} file.\n\n'
                             'Check if the file is empty or the format is incorrect.\n\n'
                             'Input file: {1}'.format(name, in_path))

    @staticmethod
    def out_msg(name, out_path):
        """
        Show a file writing error message.
        """
        messagebox.showerror('RegulonDB Converter Error',
                             'Writing the {0} file failed.\n\n'
                             'Check if the file is opened in another program and thus not'
                             'accessible for writing.\n\n'
                             'Result file: {1}'.format(name, out_path))

    def read_file(self, name, in_path, funct):
        """
        Read the specified file and process it with the provided function.
        :return: 
        """
        in_path = adjust_file_path(path.abspath(in_path))
        try:
            file = open(in_path, 'r')
            if not funct(file):
                self.empty_msg(name, in_path)
            file.close()
        except IOError:
            self.in_msg(name, in_path)

    def write_file(self, name, out_path, funct):
        """
        Write the specified file using the provided function.
        """
        out_path = self.out_dir + out_path
        try:
            file = open(out_path, 'w')
            funct(file)
            file.close()
        except IOError:
            self.out_msg(name, out_path)

    def parse_gene(self, file):
        """
        Parse the gene input file.
        """
        # remove comment lines
        lines = tsv_start(file)
        success = False

        for line in lines:
            # obtain the columns
            words = line.split('\t')

            # skip incomplete lines
            if len(words) < 10:
                continue

            gene = Gene(name=words[1].strip(), start=words[2].strip(), end=words[3].strip(),
                        strand=words[4].strip(), desc=words[6].strip(), dna=words[9].strip())

            # try to add the gene
            if self.index.add_gene(gene):
                success = True

        return success

    def parse_operon(self, file):
        """
        Parse the operon input file.
        """
        # remove comment lines
        lines = tsv_start(file)
        success = False

        for line in lines:
            # obtain the columns
            cols = line.split('\t')

            # skip incomplete lines
            if len(cols) < 8:
                continue

            # skip insufficient evidence
            if self.skip and cols[7].strip() in ['null', 'Weak', '']:
                continue

            # try to add the operon
            if self.index.add_operon(id=cols[0].strip(), genes=cols[5].strip().split(',')):
                success = True

        return success

    def parse_tus(self, file):
        """
        Parse the TU input file.
        """
        # remove comment lines
        lines = tsv_start(file)
        success = False

        for line in lines:
            # obtain the columns
            cols = line.split('\t')

            # skip incomplete lines
            if len(cols) < 7:
                continue

            # skip insufficient evidence
            if self.skip and cols[6].strip() in ['null', 'Weak', '']:
                continue

            # try to add the TU
            if self.index.add_tu(tu_id=cols[1].strip(), genes=cols[3].strip().split(',')):
                success = True

        return success

    def parse_gene_reg(self, file):
        """
        Process the gene or TF regulation input file.
        """
        # remove comment lines
        lines = tsv_start(file)
        success = False

        for line in lines:
            # obtain the columns
            words = line.split('\t')

            # skip incomplete lines
            if len(words) < 5:
                continue

            # skip entries with insufficient evidence
            if self.skip and words[4].strip() in ['null', 'Weak', '']:
                continue

            # try to add the regulatory information
            if self.index.add_gene_reg(tf=words[0].strip(), gn=words[1].strip(), reg=words[2].strip()):
                success = True

        return success

    def parse_tu_reg(self, file):
        """
        Parse the TF-TU input file.
        """
        # remove comment lines
        lines = tsv_start(file)
        success = False

        for line in lines:
            # obtain the columns
            words = line.split('\t')

            # skip incomplete lines
            if len(words) < 5:
                continue

            # skip entries with insufficient evidence
            if self.skip and words[4].strip() in ['null', 'Weak', '']:
                continue

            tu = re.sub(r'[[].*?[]]', '', words[1].strip())

            # try to add the TU regulation
            if self.index.add_tu_reg(tf=words[0].strip(), tu_id=tu, reg=words[2].strip()):
                success = True

        return success

    def parse_operon_reg(self, file):
        """
        Parse the TF-operon file.
        """
        # remove comment lines
        lines = tsv_start(file)
        success = False

        for line in lines:
            # obtain the columns
            words = line.split('\t')

            # skip incomplete lines
            if len(words) < 5:
                continue

            # skip entries with insufficient evidence
            if self.skip and words[4].strip() in ['null', 'Weak', '']:
                continue

            # parse the columns
            operon = re.split('\W+', words[1].strip())

            name = operon[0].strip()

            if len(operon) > 1:
                genes = re.split('\W+', operon[1].strip())
            else:
                genes = []

            # try to add the operon regulation
            if self.index.add_operon_reg(tf=words[0].strip(), op_id=name, genes=genes, reg=words[2].strip()):
                success = True

        return success

    def parse_tfbs(self, file):
        """
        Parse the TFBS input file.
        """
        # remove comment lines
        lines = tsv_start(file)
        success = False

        for line in lines:
            # get the columns
            cols = line.split('\t')
            # skip incomplete entries
            if len(cols) < 14:
                continue

            # skip insufficient evidence
            if self.skip and cols[13].strip() in ['null', 'Weak', '']:
                continue

            # check if the start and end positions are integers
            try:
                start = int(cols[3].strip())
                end = int(cols[4].strip())
            except ValueError:
                continue

            # process the strand
            if cols[5].strip().lower() in cfg.isup.str_minus:
                strand = cfg.misc.minus
            else:
                strand = cfg.misc.plus

            # process the DNA
            dna = ''.join([c for c in cols[11].strip() if c.isupper() and is_dna_sequence(c)])

            # try to add the TFBS
            if self.index.add_tfbs(tf=cols[1].strip(), tu_id=cols[7].strip(), prom_id=cols[9].strip(), start=start,
                                   end=end, strand=strand, dna=dna):
                success = True

        return success

    def parse_msa(self, file):
        """
        Parse the TF MSA input file.
        """
        # remove comment lines
        lines = tsv_start(file)
        success = False

        name = ''
        alignment = []

        for line in lines:
            # skip empty lines
            if not line or line == '\n':
                continue

            # skip the matrix
            if line.startswith('a') or line.startswith('c') or line.startswith('g') or line.startswith('t'):
                continue

            # skip unnecessary information
            if line.startswith('PSSM') or line.startswith('Total') or line.startswith('Transcription Factor ID'):
                continue

            # parse the gene name and add the previous entry
            if line.startswith('Transcription Factor Name'):
                if self.index.add_msa(name, alignment):
                    success = True
                # reset the alignment and gene name
                name = line.split(':')[1].strip()
                alignment = []
                continue
            # add a sequence to the alignment
            if is_dna_sequence(line.strip()):
                alignment.append(line.strip().upper())

        # add the last entry in the file
        if self.index.add_msa(name, alignment):
            success = True

        return success

    def gene_tsv(self, file):
        """
        Write the lines of the gene information result .tsv file.
        """
        # write the header line
        header = [cfg.gene.lt, cfg.gene.name, cfg.gene.desc, cfg.gene.strand, cfg.gene.start, cfg.gene.end,
                  cfg.gene.dna, cfg.gene.prot, cfg.gene.operon]
        file.write('{0}\n'.format('\t'.join(header)))

        # compile the line for each gene in the data set
        for gene in self.index.genes.values():
            line = [gene.name, gene.name, gene.desc, gene.strand, str(gene.g_start), str(gene.g_end),
                    gene.dna, gene.protein]

            if gene.operon:
                line += [cfg.misc.operon_fwd.join(gene.operon)]
            else:
                line += [cfg.misc.none]

            file.write('{0}\n'.format('\t'.join(line)))

    def reg_tsv(self, file):
        """
        Write the lines of the regulation information result .tsv file.
        """
        # write the header line
        header = [cfg.reg.lt1, cfg.reg.name1, cfg.reg.lt2, cfg.reg.name2, cfg.reg.desc]
        file.write('{0}\n'.format('\t'.join(header)))

        # compile the line for each regulation in the data set
        for key, reg in sorted(self.index.reg.items()):
            for r in sorted(reg):
                line = [key, key, r[0], r[0], r[1]]
                file.write('{0}\n'.format('\t'.join(line)))

    def msa_text(self, file):
        """
        Write the lines of the transcription factor MSA result .fasta file.
        """
        for key, align in sorted(self.index.msa.items()):
            # write the section header
            file.write('> {0}({0})\n'.format(key))
            # write the alignment
            for seq in align:
                file.write(seq + '\n')

    def tfbs_tsv(self, file):
        """
        Write the lines of the transcription factor binding site result .tsv file.
        """
        # write the header line
        header = [cfg.tf.tf_lt, cfg.tf.tf_name, cfg.tf.tfbs_lt, cfg.tf.tfbs_name, cfg.tf.astart, cfg.tf.aend,
                  cfg.tf.start, cfg.tf.seq]
        file.write('{0}\n'.format('\t'.join(header)))

        # compile the line for each binding site in the data set
        for key, lst in sorted(self.index.tfbs.items()):
            for tfbs in sorted(lst):
                line = [key, key, tfbs[0], tfbs[0], str(tfbs[1]), str(tfbs[2]), '', tfbs[3]]
                file.write('{0}\n'.format('\t'.join(line)))

    def log(self, file):
        """
        Write a log containing some statistics and error messages.
        """
        # statistics
        stats = ['[STATISTICS]',
                 'Number of genes: {}'.format(len(self.index.genes)),
                 'Number of operons (all): {}'.format(len(self.index.operons)),
                 'Number of operons (without singletons): '
                 '{}'.format(len(self.index.operons) - len(self.index.errors.operon_singletons)),
                 'Number of TFBS: {}'.format(sum([len(x) for x in self.index.tfbs.values()])),
                 'Number of regulations: {}'.format(sum([len(x) for x in self.index.reg.values()]))]
        stats = ['\n'.join(stats)]

        sections = self.index.errors.sections_condensed()

        file.write('\n\n{:-<100}\n\n'.format('').join(stats + sections))
