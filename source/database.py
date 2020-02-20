# import standard or third party modules
import os
import re
import sqlite3
from collections import defaultdict
from copy import deepcopy
from csv import DictReader
from itertools import islice
from random import randint, choice
from tkinter import messagebox
from typing import List, Set

from source.configuration import cfg, gui
from source.mutation import Mutation
from source.tools import complement_strand, translate, SubstitutionMatrix, Pwm, Edge, tsv_start, check_tsv_fieldnames, \
    is_dna_sequence

# import own modules
from source.gene import Gene

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2017 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'


class Database:
    """
    Class that represents the database containing information on all genes and mutations. It is the basis of all further
    analyses.
    """
    def __init__(self):
        self.genes = {}
        self.lt_by_gn = defaultdict(set)

        self.mutations = {}
        self.id_by_pos = defaultdict(set)

        # operons
        self.operons = set()
        # transcription factor binding site
        self.tfbs = defaultdict(lambda: {'TF': set(), 'TFBS': set(), 'mutations': set(), 'DNA': set(), 'strand': set()})
        # promoter regions
        self.prom = defaultdict(set)

        # GRN
        # IDs
        self.gene_nbr = 0
        self.mut_nbr = 0
        # genes of interest regulation and operon locus tags
        self.interest_lts = set()   # type: Set
        self.reg_lts = set()        # type: Set
        self.op_lts = set()         # type: Set
        # edges
        self.next_edge_id = 1           # type: int
        self.interest_edges = set()     # type: Set
        self.all_edges = set()          # type: Set

        # database
        self.gene_header = []   # type: List[str]
        self.mut_header = []    # type: List[str]

        # genes of interest
        self.categories = set()

        # error flags
        self.build = True
        self.coding = True

    """
    ACCESS
    """

    def by_gn(self, name):
        """
        Takes the name of a gene and tries to find it. It is possible that the gene occurs at multiple loci, in which
        case all these loci are returned.
        :param name: gene name
        :return: list with all genes with the name
        """
        # empty input field
        if name in cfg.misc.empty:
            return []
        # there is no gene with the name in the database
        if name not in self.lt_by_gn:
            cfg.a_log.db_find_gene.add(name)
            return []

        return [self.genes[lt] for lt in sorted(self.lt_by_gn[name])]

    def by_lt(self, lt):
        """
        Takes the locus tag of a gene to find it in the database. There can only be only one gene with any given
        locus tag.
        :param lt: locus tag of the gene
        :return: the gene or None if it is not in the database
        """
        # empty input field
        if lt in cfg.misc.empty:
            return None
        # there is not gene with the locus tag in the database
        if lt not in self.genes:
            cfg.a_log.db_find_gene.add(lt)
            return None

        return self.genes[lt]

    def find_gene(self, lt, name):
        """
        First tries to find the gene via the locus tag. If the locus tag is not given, tries to find the gene(s) via
        the gene name.
        :param lt: locus tag of a gene
        :param name: name of gene(s)
        :return: list of gene(s) with the locus tag or gene name
        """
        # if the locus tag is not given, try to find genes by gene name
        if lt in cfg.misc.empty or lt not in self.genes:
            return self.by_gn(name)
        # if the locus tag is given, try to find the gene with the locus tag
        gene = self.by_lt(lt)

        if gene:
            return [gene]
        return []

    def by_id(self, i):
        """
        Finds a mutation in the database by its ID.
        :param i: mutation ID
        :return: mutation if in database, None otherwise
        """
        if i in self.mutations:
            return self.mutations[i]

        cfg.a_log.db_find_mut.add(i)
        return None

    def by_pos(self, p):
        """
        Finds mutation(s) by their position.
        :param p: mutation position
        :return: list of mutation(s) if in database, empty list otherwise
        """
        if p in self.id_by_pos:
            return [self.mutations[i] for i in sorted(self.id_by_pos[p])]
        return []

    """
    BUILD DATABASE
    """

    def build_database(self):
        """
        Sets up the database by reading and processing the input files.
        """
        # process gene information
        self._read_gene_tsv()
        self._check_and_set_operons()
        if cfg.run.a_adop:
            self._adjust_promoters()
        self._build_promoter_set()
        # process mutation information
        self._read_mutation_tsv()

        if cfg.run.a_int:
            self._read_genes_of_interest_tsv()
        if cfg.run.a_reg:
            self._read_regulation_tsv()
        if cfg.run.a_prot:
            self._read_protein_domain_tsv()
        if cfg.run.a_tfbs:
            self._read_tfbs_tsv()
            self._read_tfbs_msa_fasta()
            if cfg.run.a_adop:
                self._adjust_tfbs()

        self._intersect()

        return self.build

    def _read_gene_tsv(self):
        """
        Reads and processes the .tsv file with gene information.
        """
        # open the .tsv file for reading and load it into a dictionary reader
        with open(cfg.user.ai_gene, 'r', encoding='latin-1') as file:
            reader = DictReader(tsv_start(file), delimiter='\t')
            # check if all required field names are in the .tsv file
            field_names = {cfg.gene.lt, cfg.gene.strand, cfg.gene.start, cfg.gene.end, cfg.gene.dna,
                           cfg.gene.prot}

            if not check_tsv_fieldnames(gui.desc[gui.key.ai_gene].label, cfg.user.ai_gene,
                                        field_names, reader.fieldnames):
                self.build = False
                return

            # each line represents information belonging to a single gene
            for line in reader:
                lt = line[cfg.gene.lt]

                # no locus tag
                if lt in cfg.misc.empty:
                    continue

                # locus tags need to be unique
                if lt in self.genes:
                    cfg.a_log.db_duplicate_genes.add(lt)
                    continue

                # initialise the gene with the required information and stop if it is not valid
                gene = Gene()
                if not gene.set_locus_tag(lt):
                    continue
                if not gene.set_strand(line[cfg.gene.strand]):
                    continue
                if not gene.set_gene_start_and_end(line[cfg.gene.start], line[cfg.gene.end]):
                    continue
                if not gene.set_dna(line[cfg.gene.dna]):
                    continue
                if not gene.set_coding_length():
                    continue
                if not gene.set_protein(line[cfg.gene.prot]):
                    continue

                # set the promoter region based on the available information
                if cfg.gene.prom_astart in reader.fieldnames and cfg.gene.prom_end in reader.fieldnames:
                    p_astart = line[cfg.gene.prom_astart]
                    p_end = line[cfg.gene.prom_end]
                else:
                    p_astart = cfg.misc.none
                    p_end = cfg.misc.none

                if cfg.gene.prom_rstart in reader.fieldnames:
                    p_rstart = line[cfg.gene.prom_rstart]
                else:
                    p_rstart = cfg.misc.none
                gene.set_promoter(p_astart, p_end, p_rstart)

                # set the gene name if it is in the file
                if cfg.gene.name in reader.fieldnames:
                    name = line[cfg.gene.name]
                    gene.set_gene_name(name)
                    if name not in cfg.misc.empty:
                        self.lt_by_gn[name].add(lt)

                # add operon information if it is in the file
                if cfg.gene.operon in reader.fieldnames:
                    self._add_operon(line[cfg.gene.operon])

                # set the description if it is in the file
                if cfg.gene.desc in reader.fieldnames:
                    gene.set_description(line[cfg.gene.desc])

                # set the accession numbers if they are in the file
                if cfg.gene.gene_id in reader.fieldnames:
                    gene.set_gene_id(line[cfg.gene.gene_id])
                if cfg.gene.genbank in reader.fieldnames:
                    gene.set_genbank(line[cfg.gene.genbank])
                if cfg.gene.pfam in reader.fieldnames:
                    gene.set_pfam(line[cfg.gene.pfam])
                if cfg.gene.uni in reader.fieldnames:
                    gene.set_uniprot(line[cfg.gene.uni])
                if cfg.gene.pubmed in reader.fieldnames:
                    gene.set_pubmed(line[cfg.gene.pubmed])
                if cfg.gene.chromosome in reader.fieldnames:
                    gene.set_chromosome(line[cfg.gene.chromosome])

                # set the gene regulatory network ID and add the gene to the database
                self.gene_nbr += 1
                gene.grn_id = self.gene_nbr
                self.genes[gene.locus_tag] = gene

    def _add_operon(self, field):
        """
        Parse an operon and add it to the operon set if it is valid.
        :param field: field with operon information
        """
        # input field is empty
        if field in cfg.misc.empty:
            return

        # only operons with one direction are supported
        if (cfg.misc.operon_fwd in field) and (cfg.misc.operon_bwd in field):
            cfg.a_log.db_operon_inconsistent.add(field)
            return

        # parse the operon field and add it to the operon set
        if cfg.misc.operon_fwd in field:
            genes = field.split(cfg.misc.operon_fwd)
        elif cfg.misc.operon_bwd in field:
            genes = field.split(cfg.misc.operon_bwd)[::-1]
        # operon format is not supported
        else:
            cfg.a_log.db_operon_unsupported.add(field)
            return
        # make sure that no gene occurs twice
        if len(genes) != len(set(genes)):
            cfg.a_log.db_operon_with_duplicate_genes.add(field)
            return

        self.operons.add(tuple(genes))

    def _check_and_set_operons(self):
        """
        Removes all operons that contain genes that are not in the database or that occur in more than
        one operon. Then set the operon of all genes in the valid operons.
        """
        # remove all operons that contain genes that are not in the database
        new = set(filter(lambda x: False not in [lt in self.genes for lt in x], self.operons))
        cfg.a_log.db_operon_not_in_database.set = self.operons - new
        self.operons = new

        # make sure all genes are at most in one operon
        genes = set()
        remove = set()

        for operon in self.operons:
            strand = self.genes[operon[0]].strand

            for lt in operon:
                # if a gene is already in another operon, mark this operon for removal
                if lt in genes:
                    remove.add(operon)
                    cfg.a_log.db_genes_in_several_operons.add(lt)
                # if the genes are on different strands, mark this operon for removal
                if self.genes[lt].strand != strand:
                    remove.add(operon)
                    cfg.a_log.db_operon_inconsistent_strand.add(operon)
                # add the locus tag to the gene set
                genes.add(lt)

        # remove the operons with duplicate genes
        self.operons -= remove

        # set the operon for the associated genes
        for operon in self.operons:
            for lt in operon:
                self.genes[lt].operon = operon

    def _adjust_promoters(self):
        """
        Adjusts the promoter regions of operons by setting it to the promoter region of the first gene in the operon.
        """
        for operon in self.operons:
            first = self.genes[operon[0]]   # type: Gene

            for lt in operon[1:]:
                gene = self.genes[lt]       # type: Gene

                gene.prom_start_abs = first.prom_start_abs
                gene.prom_end_abs = first.prom_end_abs
                gene.prom_length = first.prom_length

                if gene.strand == cfg.misc.plus:
                    gene.prom_start_rel = first.prom_start_abs - gene.gene_start
                else:
                    gene.prom_start_rel = gene.gene_end - first.prom_end_abs

    def _build_promoter_set(self):
        """
        Builds the set of promoters.
        """
        for gene in self.genes.values():
            if gene.prom_start_abs > 0 and gene.prom_end_abs > 0:
                self.prom[(gene.chromosome, gene.prom_start_abs, gene.prom_end_abs)].add(gene.locus_tag)

    def _read_mutation_tsv(self):
        """
        Reads and processes the .tsv file with mutation information.
        """
        # open the .tsv file for reading and load it into a dictionary reader
        with open(cfg.user.ai_mut, 'r', encoding='latin-1') as file:
            reader = DictReader(tsv_start(file), delimiter='\t')

            # check if all required field names are in the .tsv file
            field_names = {cfg.mut.pos, cfg.mut.ref, cfg.mut.alt}

            if not check_tsv_fieldnames(gui.desc[gui.key.ai_mut].label, cfg.user.ai_mut,
                                        field_names, reader.fieldnames):
                self.build = False
                return

            mutations = set()
            # read the mutation lines and eliminate duplicates
            for line in reader:
                chrom = cfg.misc.none_out if cfg.mut.chrom not in reader.fieldnames else line[cfg.mut.chrom]
                # in case there are several alternatives given, create a mutation object for each one
                for alt in re.split('\W+', line[cfg.mut.alt]):
                    mutations.add((chrom, line[cfg.mut.pos], line[cfg.mut.ref], alt))

            # add the mutations to the database
            for chrom, pos, ref, alt in sorted(mutations):
                mutation = Mutation()
                # make sure the position is valid
                if not mutation.set_pos(pos):
                    continue

                # make sure the reference and alternative base is valid
                if not mutation.set_ref_and_alt(ref, alt):
                    continue

                mutation.set_chromosome(chrom)

                self.mut_nbr += 1
                mutation.ID = self.mut_nbr
                self.mutations[self.mut_nbr] = mutation
                self.id_by_pos[mutation.pos].add(mutation.ID)

    def _read_genes_of_interest_tsv(self):
        """
        Reads and processes the .tsv file with the genes of interest.
        """
        # open the .tsv file for reading and load it into a dictionary reader
        with open(cfg.user.ai_int, 'r', encoding='latin-1') as file:
            reader = DictReader(tsv_start(file), delimiter='\t')

            # check if all required field names are in the .tsv file
            field_names = {cfg.int.lt, cfg.int.name, cfg.int.desc}

            if not check_tsv_fieldnames(gui.desc[gui.key.ai_int].label, cfg.user.ai_int,
                                        field_names, reader.fieldnames):
                self.build = False
                return

            for line in reader:
                if not line[cfg.int.desc]:
                    continue

                temp = ' '.join(line[cfg.int.desc].split()).strip()
                if temp in cfg.misc.empty:
                    interest = cfg.user.as_int
                elif re.match(r'[\w -]*$', temp):
                    interest = temp
                else:
                    continue

                # find all genes associated with the locus tag and/or gene name
                for gene in self.find_gene(line[cfg.int.lt], line[cfg.int.name]):
                    gene.set_interest(interest)
                    cfg.cat.add(interest)

                    # add genes of interest to the sub-network of interest
                    self.interest_lts.add(gene.locus_tag)

            # too many sub-categories
            if len(cfg.cat) > 10:
                self.build = False

    def _read_regulation_tsv(self):
        """
        Reads and processes the .tsv file with regulatory information.
        """
        # open the .tsv file for reading and load it into a dictionary reader
        with open(cfg.user.ai_reg, 'r', encoding='latin-1') as file:
            reader = DictReader(tsv_start(file), delimiter='\t')

            # check if all required field names are in the .tsv file
            field_names = {cfg.reg.lt1, cfg.reg.lt2, cfg.reg.name1, cfg.reg.name2, cfg.reg.desc}

            if not check_tsv_fieldnames(gui.desc[gui.key.ai_reg].label, cfg.user.ai_reg,
                                        field_names, reader.fieldnames):
                self.build = False
                return

            for line in reader:
                # find the regulator(s) and regulated gene(s) in the database
                regulators = self.find_gene(line[cfg.reg.lt1], line[cfg.reg.name1])
                regulated_genes = self.find_gene(line[cfg.reg.lt2], line[cfg.reg.name2])

                # add the regulatory information to the regulator(s) and regulated gene(s)
                for regulator in regulators:
                    for regulated in regulated_genes:
                        regulator.add_regulated_gene(regulated.locus_tag, regulated.gene_name, line[cfg.reg.desc])
                        regulated.add_regulator(regulator.locus_tag, regulator.gene_name, line[cfg.reg.desc])

    def _read_protein_domain_tsv(self):
        """
        Reads and processes the .tsv file with transcription factor binding site information.
        """
        # open the .tsv file for reading and load it into a dictionary reader
        with open(cfg.user.ai_prot, 'r', encoding='latin-1') as file:
            reader = DictReader(tsv_start(file), delimiter='\t')

            # check if all required field names are in the .tsv file
            field_names = {cfg.prot.lt, cfg.prot.name, cfg.prot.start, cfg.prot.end, cfg.prot.desc, cfg.prot.type}

            if not check_tsv_fieldnames(gui.desc[gui.key.ai_prot].label, cfg.user.ai_prot,
                                        field_names, reader.fieldnames):
                self.build = False
                return

            for line in reader:
                # find the gene(s) in the database
                genes = self.find_gene(line[cfg.prot.lt], line[cfg.prot.name])

                # add the protein domain information to the gene(s)
                for gene in genes:
                    gene.add_protein_domain(line[cfg.prot.start], line[cfg.prot.end],
                                            line[cfg.prot.type], line[cfg.prot.desc])

    def _read_tfbs_tsv(self):
        """
        Reads and processes the .tsv file with transcription factor binding site information.
        """
        # open the .tsv file for reading and load it into a dictionary reader
        with open(cfg.user.ai_tfbs, 'r', encoding='latin-1') as file:
            reader = DictReader(tsv_start(file), delimiter='\t')

            # check if all required field names are in the .tsv file
            field_names = {cfg.tf.tf_lt, cfg.tf.tf_name, cfg.tf.tfbs_lt, cfg.tf.tfbs_name, cfg.tf.seq, cfg.tf.start,
                           cfg.tf.astart, cfg.tf.aend}

            if not check_tsv_fieldnames(gui.desc[gui.key.ai_tfbs].label, cfg.user.ai_tfbs,
                                        field_names, reader.fieldnames):
                self.build = False
                return

            for line in reader:
                # find the regulator(s) and regulated gene(s)
                regulators = self.find_gene(line[cfg.tf.tf_lt], line[cfg.tf.tf_name])
                regulated_genes = self.find_gene(line[cfg.tf.tfbs_lt], line[cfg.tf.tfbs_name])

                # the regulator(s) and/or regulated gene(s) are not in the database
                if not regulators or not regulated_genes:
                    continue

                # the sequence is not a DNA sequence
                dna = line[cfg.tf.seq].upper()
                if dna in cfg.misc.empty:
                    cfg.a_log.tfbs_not_dna.add((tuple(regulators), tuple(regulated_genes)))
                    continue

                if not is_dna_sequence(dna):
                    cfg.a_log.tfbs_not_dna.add((tuple(regulators), tuple(regulated_genes)))
                    continue

                # check if the absolute position is given
                astart, aend = self._tfbs_from_abs(line[cfg.tf.astart], line[cfg.tf.aend], dna)

                # check if the relative position is valid if the absolute one is invalid
                if not astart and not aend:
                    if len(regulated_genes) > 1:
                        cfg.a_log.tfbs_no_unique_pos.add((tuple(regulators), tuple(regulated_genes)))
                        continue
                    astart, aend = self._tfbs_from_rel(line[cfg.tf.start], dna, regulated_genes[0])

                # no valid position
                if not astart and not aend:
                    cfg.a_log.tfbs_no_position.add((tuple(regulators), tuple(regulated_genes)))
                    continue

                self.tfbs[(astart, aend)]['TF'].update([x.locus_tag for x in regulators])
                self.tfbs[(astart, aend)]['TFBS'].update([x.locus_tag for x in regulated_genes])
                self.tfbs[(astart, aend)]['DNA'].add(dna)
                self.tfbs[(astart, aend)]['strand'].update([x.strand for x in regulated_genes])

        delete = set()
        # remove TFBS with more than one DNA sequence
        for pos in self.tfbs.keys():
            if len(self.tfbs[pos]['DNA']) != 1:
                cfg.a_log.tfbs_multiple_dna.add((tuple(self.tfbs[pos]['TF']), tuple(self.tfbs[pos]['TFBS'])))
                delete.add(pos)
            if len(self.tfbs[pos]['strand']) != 1:
                cfg.a_log.tfbs_multiple_strand.add((tuple(self.tfbs[pos]['TF']), tuple(self.tfbs[pos]['TFBS'])))
                delete.add(pos)
            else:
                self.tfbs[pos]['strand'] = list(self.tfbs[pos]['strand'])[0]

        for pos in delete:
            del self.tfbs[pos]

        # add the TFBS to the genes
        for pos, tfbs in self.tfbs.items():
            for r1 in tfbs['TF']:
                for r2 in tfbs['TFBS']:
                    regulated = self.genes[r2]
                    regulator = self.genes[r1]
                    regulated.add_regulator(regulator.locus_tag, regulator.gene_name, cfg.misc.unknown)
                    regulator.add_regulated_gene(regulated.locus_tag, regulated.gene_name, cfg.misc.unknown)
                    regulated.add_tfbs(regulator.locus_tag, pos[0], pos[1], list(tfbs['DNA'])[0])

    @staticmethod
    def _tfbs_from_rel(rstart, seq, gene):
        """
        Try to obtain the absolute start and end position of the transcription factor binding site from the
        relative start.
        :param rstart: relative start
        :param seq: DNA sequence of the binding site
        :param gene: regulated gene
        :return: absolute start, end
        """
        # relative start not given
        if rstart in cfg.misc.empty:
            return 0, 0

        # must be an integer
        try:
            rstart = int(rstart)
        except ValueError:
            return 0, 0

        # the coding region is not defined
        if not gene.gene_start or not gene.gene_end:
            return 0, 0

        if gene.strand == cfg.misc.minus:
            aend = gene.gene_end - rstart
            astart = aend - len(seq) + 1
        else:
            astart = gene.gene_start + rstart
            aend = astart + len(seq) - 1

        # return valid absolute position
        if aend > 0 and astart > 0:
            return min(astart, aend), max(astart, aend)

        return 0, 0

    @staticmethod
    def _tfbs_from_abs(start, end, seq):
        """
        Check if the given absolute start and end are valid.
        :param start: absolute TFBS start
        :param end: absolute TFBS end
        :param seq: TFBS DNA sequence
        :return: processed absolute start, end
        """
        # not given
        if start in cfg.misc.empty or end in cfg.misc.empty:
            return 0, 0

        # must be integers
        try:
            start = int(start)
            end = int(end)
        except ValueError:
            return 0, 0

        # outside gneome
        if start <= 0 or end <= 0:
            return 0, 0

        # does not match DNA sequence length
        if (abs(start - end) + 1) != len(seq):
            return 0, 0

        return min(start, end), max(start, end)

    def _adjust_tfbs(self):
        """
        Adjusts the TFBSs of operons.
        """
        for operon in self.operons:
            # build the set of all TFBSs in the operon
            # (start, end, DNA, TF locus tag)
            tfbs = set()

            for lt in operon:
                gene = self.genes[lt]
                tfbs.update(gene.tfbs.keys())

            for t in tfbs:
                # add the operon locus tags to the TFBS
                self.tfbs[t[:2]]['TFBS'].update(operon)
                tf = self.genes[t[3]]

                # add the TFBS to all regulated genes, and the regulatory information to the regulator and the
                # regulated genes
                for lt in operon:
                    gene = self.genes[lt]
                    gene.add_tfbs(tf.locus_tag, *t[:-1])
                    tf.add_regulated_gene(gene.locus_tag, gene.gene_name, cfg.misc.unknown)
                    gene.add_regulator(tf.locus_tag, tf.gene_name, cfg.misc.unknown)

    def _read_tfbs_msa_fasta(self):
        """
        Reads and processes the .fasta file with the transcription factor sequence alignments.
        """
        with open(cfg.user.ai_msa, 'r', encoding='latin-1') as file:
            lt = ''
            name = ''
            alignment = []

            for line in file:
                # skip empty lines
                if line in cfg.misc.empty or not line.strip():
                    continue

                # new transcription factor
                if line[0] == '>':
                    # compute the position weight matrix for the previous transcription factor
                    self._add_alignment(lt, name, alignment)

                    # set the locus tag of the next transcription factor and reset the alignment
                    words = re.split('\W+', line[1:].strip())

                    lt = words[0]

                    if len(words) > 1:
                        name = words[1]
                    else:
                        name = ''

                    alignment = []

                # add sequence to the alignment
                else:
                    alignment.append(line.strip())

            # compute the position weight matrix of the last transcription factor
            self._add_alignment(lt, name, alignment)

    def _add_alignment(self, lt, name, alignment):
        """
        Adds a TF MSA to the corresponding transcription factor(s) if the alignment is valid.
        :param lt: transcription factor locus tag
        :param name: transcription factor gene name
        :param alignment: list of transcription factor sequences
        """
        # the locus tag or the gene name must be given
        if lt in cfg.misc.empty and name in cfg.misc.empty:
            return

        genes = self.find_gene(lt, name)

        # no gene found in the database
        if not genes:
            return

        # build the PWM
        pwm = Pwm('{0} ({1})'.format(lt, name), alignment)

        # not enough valid sequence in the alignment
        if not pwm.success:
            return

        # add the PWM to all associated transcription factors
        for gene in genes:
            gene.tf_pwm = deepcopy(pwm)
            gene.tf_pwm.tf_lt = gene.locus_tag

    def _intersect(self):
        """
        Computes which mutations occur in which gene regions and cross-links the gene and mutation information.
        """
        # sort the mutations by chromosome and position
        mutation_dict = defaultdict(set)
        for mutation in self.mutations.values():
            mutation_dict[mutation.chromosome].add(mutation)

        mutation_dict = {chrom: sorted(mutations, key=lambda x: x.pos) for chrom, mutations in mutation_dict.items()}
        # build a set of gene regions
        region_dict = self._get_regions()

        for chromosome in set(mutation_dict.keys()).intersection(region_dict.keys()):
            mutations = mutation_dict[chromosome]
            regions = region_dict[chromosome]
            # index of the next gene region
            g_index = 0
            # number of genes in the database
            g_max = len(regions)

            for mut in mutations:
                ID = mut.ID
                pos = mut.pos

                for start, end, lt, interest, loc in islice(regions, g_index, g_max):
                    # mutation pos is smaller than the start of the gene region, therefore, all regions for this
                    # mutation have been found due to sorting
                    if pos < start:
                        break

                    # mutation is larger than the end of the gene region, therefore, it is also larger than the end of
                    # the previous gene regions and the next mutation is also larger than the end of this and previous
                    # regions => can be skipped
                    if pos > end:
                        g_index += 1
                        continue

                    # intersection found, cross-link the information
                    self.mutations[ID].add_gene_region(start, end, lt, interest, loc)
                    gene = self.genes[lt]
                    gene.add_mutation_location(ID, loc)

                    # add the corresponding protein domains
                    if cfg.run.a_prot and loc == cfg.misc.coding:
                        domains = gene.get_protein_domains_for_pos(pos)
                        if domains:
                            self.mutations[ID].add_protein_domain(gene.gene_start, gene.gene_end, lt, domains)

                    if loc == cfg.misc.tfbs:
                        self.tfbs[(start, end)]['mutations'].add(ID)

    def _get_regions(self):
        """
        Compute all gene regions in the database.
        :return: sorted list of all gene regions
        """
        regions = defaultdict(set)
        # promoter regions
        for key, val in self.prom.items():
            for lt in val:
                interest = self.genes[lt].interest
                regions[key[0]].add((key[1], key[2], lt, interest, cfg.misc.prom))
        # TFBS regions
        for key, val in self.tfbs.items():
            for lt in val['TFBS']:
                interest = self.genes[lt].interest
                regions[self.genes[lt].chromosome].add((key[0], key[1], lt, interest, cfg.misc.tfbs))
        # coding regions
        for key, val in self.genes.items():
            if val.gene_start <= 0 or val.gene_end <= 0:
                continue
            interest = val.interest
            regions[val.chromosome].add((val.gene_start, val.gene_end, key, interest, cfg.misc.coding))

        return {chromosome: sorted(region, key=lambda x: (x[0], x[1], x[2])) for chromosome, region in regions.items()}

    """
    WRITE DATABASE
    """

    def write_database(self):
        """
        Sets up the SQLite database and fills it with the information in the database.
        """
        try:
            connect = sqlite3.connect(cfg.op.db)
            connect.isolation_level = None
            cursor = connect.cursor()

            self._setup_database(cursor)
            self._fill_gene_table(cursor)
            self._fill_mutation_table(cursor)

            connect.close()
            return True
        except Exception as e:
            print(e)
            return False

    def _setup_database(self, cursor):
        """
        Sets up the tables in the SQLite database.
        :param cursor: database cursor
        """
        # delete the tables if they already exist to reset them
        cursor.execute('DROP TABLE IF EXISTS {0}'.format(cfg.gdb.table))
        cursor.execute('DROP TABLE IF EXISTS {0}'.format(cfg.mdb.table))

        # short hand for column types
        i = 'INTEGER'
        t = 'TEXT'

        # assemble the column names for the genes and the mutation table
        gene = [(cfg.gdb.name, t),
                (cfg.gdb.int, t),
                (cfg.gdb.desc, t),
                (cfg.gdb.chrom, t),
                (cfg.gdb.strand, t),
                (cfg.gdb.gstart, i),
                (cfg.gdb.gend, i),
                (cfg.gdb.glen, i),
                (cfg.gdb.dna, t),
                (cfg.gdb.prot_len, i),
                (cfg.gdb.prot, t),
                (cfg.gdb.id, t),
                (cfg.gdb.genbank, t),
                (cfg.gdb.pfam, t),
                (cfg.gdb.uniprot, t),
                (cfg.gdb.pubmed, t),
                (cfg.gdb.prom_rstart, i),
                (cfg.gdb.prom_astart, i),
                (cfg.gdb.prom_end, i),
                (cfg.gdb.prom_len, i),
                (cfg.gdb.op, t),
                (cfg.gdb.reg1, t),
                (cfg.gdb.reg2, t),
                (cfg.gdb.tf, t),
                (cfg.gdb.tfbs_rstart, t),
                (cfg.gdb.tfbs_astart, t),
                (cfg.gdb.tfbs_end, t),
                (cfg.gdb.tfbs_seq, t),
                (cfg.gdb.mut_id, t),
                (cfg.gdb.md_syn, i),
                (cfg.gdb.md_ns, i),
                (cfg.gdb.md_mis, i),
                (cfg.gdb.md_non, i),
                (cfg.gdb.md_rt, i),
                (cfg.gdb.md_rf, i),
                (cfg.gdb.md_prom, i),
                (cfg.gdb.md_tfbs, i),
                (cfg.gdb.pd_start, t),
                (cfg.gdb.pd_end, t),
                (cfg.gdb.pd_type, t),
                (cfg.gdb.pd_desc, t)]

        mut = [(cfg.mdb.chrom, t),
               (cfg.mdb.pos, i),
               (cfg.mdb.ref, t),
               (cfg.mdb.alt, t),
               (cfg.mdb.type, t),
               (cfg.mdb.reg, t),
               (cfg.mdb.lt, t),
               (cfg.mdb.int, t),
               (cfg.mdb.tfbs, t),
               (cfg.mdb.aa, t),
               (cfg.mdb.sub, t),
               (cfg.mdb.eff, t),
               (cfg.mdb.pd, t)]

        # create the gene table
        cursor.execute(self._setup_cmd(cfg.gdb.table, cfg.gdb.lt, t, gene))
        # create the mutation table
        cursor.execute(self._setup_cmd(cfg.mdb.table, cfg.mdb.id, i, mut))

        self.gene_header = [cfg.gdb.lt] + [x[0] for x in gene]
        self.mut_header = [cfg.mdb.id] + [x[0] for x in mut]

    @staticmethod
    def _setup_cmd(table, pk, pkt, cols):
        """
        Builds the SQLite command for setting up a table with primary key and columns.
        :param table: table name
        :param pk: primary key name
        :param pkt: primary key type
        :param cols: a list of (column name, column type) tuples
        :return: SQLite table set up command
        """
        cmd = 'CREATE TABLE {0} ({1} {2} PRIMARY KEY, {3})'
        cols_str = ', '.join([' '.join(temp) for temp in cols])

        return cmd.format(table, pk, pkt, cols_str)

    def _fill_gene_table(self, cursor):
        """
        Fills the gene table with the gene information. Parts that have not been analysed are skipped.
        :param cursor: database cursor
        """
        # input command
        values_fmt = self._get_db_value_format(len(self.gene_header))

        try:
            # open the tsv file and write the header
            file = open(cfg.op.db_gene, 'w')
            file.write('{0}\n'.format('\t'.join(self.gene_header)))

            # start the committing cycle
            cursor.execute('BEGIN')

            # iterate over all genes
            for key, gene in sorted(self.genes.items()):
                # assemble the columns
                values = [gene.locus_tag,
                          gene.gene_name,
                          gene.interest,
                          gene.description,
                          gene.chromosome,
                          gene.strand,
                          gene.gene_start,
                          gene.gene_end,
                          gene.gene_length,
                          gene.dna,
                          gene.prot_length,
                          gene.prot,
                          gene.gene_id,
                          gene.genbank,
                          gene.pfam,
                          gene.uniprot,
                          gene.get_pubmed(),
                          gene.prom_start_rel,
                          gene.prom_start_abs,
                          gene.prom_end_abs,
                          gene.prom_length,
                          gene.get_operon(),
                          gene.get_regulated_genes(),
                          gene.get_regulators(),
                          gene.get_tfs(),
                          gene.get_tfbs_start_rel(),
                          gene.get_tfbs_start_abs(),
                          gene.get_tfbs_end(),
                          gene.get_tfbs_seq(),
                          gene.get_mutation_ids(),
                          gene.syn_md,
                          gene.non_syn_md,
                          gene.missense_md,
                          gene.nonsense_md,
                          gene.readthrough_md,
                          gene.frame_shift_md,
                          gene.prom_md,
                          gene.tfbs_md,
                          gene.get_prot_dom_start(),
                          gene.get_prot_dom_end(),
                          gene.get_prot_dom_type(),
                          gene.get_prot_dom_description()]

                file.write('{}\n'.format('\t'.join([str(x) for x in values])))

                insert_cmd = 'INSERT INTO {0} VALUES {1}'.format(cfg.gdb.table, values_fmt)

                # insert their data into the database
                cursor.execute(insert_cmd, tuple(values))

            # commit the changes
            cursor.execute('COMMIT')

            file.close()

        except IOError:
            self.write_db_error = True

    def _fill_mutation_table(self, cursor):
        """
        Fills the mutation table with mutation information.
        :param cursor: database cursor
        """
        # set up the insertion argument format
        values_fmt = self._get_db_value_format(len(self.mut_header))

        try:
            # open the tsv file and write the header
            file = open(cfg.op.db_mut, 'w')
            file.write('{0}\n'.format('\t'.join(self.mut_header)))

            # start the committing cycle
            cursor.execute('BEGIN')

            # iterate over the mutations in the dictionary
            for key, mutation in sorted(self.mutations.items()):
                # set up the insertion argument values
                # general information
                values = [mutation.ID,
                          mutation.chromosome,
                          mutation.pos,
                          mutation.ref,
                          mutation.alt,
                          mutation.type,
                          mutation.get_regions(),
                          mutation.get_locus_tags(),
                          mutation.get_interest(),
                          mutation.get_tfbs_scores(),
                          mutation.get_amino_change(),
                          mutation.get_substitution_score(),
                          mutation.get_effect(),
                          mutation.get_protein_domains()]

                file.write('{}\n'.format('\t'.join([str(x) for x in values])))

                insert_cmd = 'INSERT INTO {0} VALUES {1}'.format(cfg.mdb.table, values_fmt)

                # insert their data into the database
                cursor.execute(insert_cmd, tuple(values))

            # commit the changes
            cursor.execute('COMMIT')

            file.close()

        except IOError:
            self.write_db_error = True

    @staticmethod
    def _get_db_value_format(nbr_values):
        """
        :param nbr_values: number of values
        :return: '(?,?,...,?,?)'
        """
        return '({0})'.format(','.join(['?'] * nbr_values))

    def write_gene_list(self):
        """
        Writes a comprehensive list in .tsv format of all genes of interest or their regulation.
        For each gene, the list includes the locus tag, gene name, category information, regulated genes, regulators
        and all non-synonymous, transcription factor binding site and promoter mutations in that gene.
        """
        try:
            file = open(cfg.op.gm_int, 'w')

            # set up the header
            header = [cfg.gm.lt, cfg.gm.name, cfg.gm.int, cfg.gm.chrom, cfg.gm.reg1, cfg.gm.reg2, cfg.gm.mut, cfg.gm.id,
                      cfg.gm.gpos, cfg.gm.ppos, cfg.gm.type, cfg.gm.loc, cfg.gm.tfbs, cfg.gm.eff, cfg.gm.aa,
                      cfg.gm.sub, cfg.gm.desc]
            empty_mut_cols = [''] * 11

            file.write('{0}\n'.format('\t'.join(header)))

            # list of locus tags in the sub-network of interest
            sub_network = sorted(self.interest_lts) + sorted(self.reg_lts)

            for lt in sub_network:
                gene = self.by_lt(lt)
                if not gene:
                    continue

                # get list of mutations in the gene
                mutations = sorted(list(gene.non_syn_mut) + list(gene.tfbs_mut) + list(gene.prom_mut))

                # get regulators and regulated genes in the sub-network of interest
                regulated = [lt for lt in gene.regulated.keys() if lt in sub_network]
                regulators = [lt for lt in gene.regulators.keys() if lt in sub_network]

                len_regulated = len(regulated)
                len_regulators = len(regulators)
                len_mutations = len(mutations)

                max_index = max(len_regulated, len_regulators, len_mutations, 1)

                # write the lines for the gene
                for i in range(0, max_index):
                    # the first line contains general information on the gene
                    if i == 0:
                        line = [gene.locus_tag, gene.gene_name]

                        if gene.interest == cfg.misc.none:
                            line.append('regulator')
                        else:
                            line.append(gene.interest)

                        line.append(gene.chromosome)
                    else:
                        line = ['', '', '', '']

                    # add the next regulated gene to the line
                    if i < len_regulated:
                        line.append(gene.get_regulated_genes_by_locus_tag(regulated[i]))
                    else:
                        line.append('')

                    # add the next regulator to the line
                    if i < len_regulators:
                        line.append(gene.get_regulator_by_locus_tag(regulators[i]))
                    else:
                        line.append('')

                    # all mutations processed
                    if i >= len_mutations:
                        line += empty_mut_cols
                        file.write('{0}\n'.format('\t'.join(line)))
                        continue

                    # mutation is not in the database
                    if not mutations[i] in self.mutations.keys():
                        line += empty_mut_cols
                        file.write('{0}\n'.format('\t'.join(line)))
                        continue

                    # add the next mutation information to the line
                    mut = self.mutations[mutations[i]]

                    # no gene information for the mutation
                    if lt not in mut.region_ids_by_lt.keys():
                        line += empty_mut_cols
                        file.write('{0}\n'.format('\t'.join(line)))
                        continue

                    # add mutation information
                    line += ['{0} --> {1}'.format(mut.ref, mut.alt),
                             str(mut.ID),
                             str(mut.pos),
                             str(gene.get_protein_position(mut.pos)),
                             mut.type,
                             mut.get_regions(lt).replace(cfg.misc.s_sep, ', '),
                             mut.get_tfbs_scores(lt).replace(cfg.misc.s_sep, ', '),
                             mut.get_effect(lt).replace(cfg.misc.s_sep, ', '),
                             mut.get_amino_change(lt).replace(cfg.misc.s_sep, ', '),
                             mut.get_substitution_score(lt).replace(cfg.misc.s_sep, ', '),
                             mut.get_protein_domains(lt).replace(cfg.misc.s_sep, ', ')]

                    file.write('{0}\n'.format('\t'.join(line)))
            file.close()
            return True

        except IOError:
            return False

    def write_mutation_list(self):
        """
        Writes a list of all known mutations from the protein domain data set that are in genes of interest. The
        output is in .tsv format.
        """
        try:
            file = open(cfg.op.mut, 'w')
            # write header
            file.write('{0}\n'.format('\t'.join([cfg.gm.lt, cfg.gm.name, cfg.gm.int, cfg.gm.chrom, cfg.gm.ppos,
                                                 cfg.gm.desc])))

            # iterate over all genes in the regulatory sub-network of interest
            for lt in sorted(self.interest_lts) + sorted(self.reg_lts):
                gene = self.by_lt(lt)
                if not gene:
                    continue

                first = True

                for domain in gene.mut_dom:
                    # if it is the first line for the gene, add the locus tag and gene name to the row
                    if first:
                        line = [gene.locus_tag, gene.gene_name]
                        first = False
                    else:
                        line = ['', '']

                    # if no category of interest is not given, the gene is a direct or indirect regulator
                    if gene.interest != cfg.misc.none:
                        line.append(gene.interest)
                    else:
                        line.append('regulation')

                    line.append(gene.chromosome)

                    # calculate the length of the domain
                    dif = domain[0] - domain[1]

                    # if the length is 0, the mutation is a SNP
                    if dif == 0:
                        line.append(str(domain[0]))
                    else:
                        line.append('{0} - {1}'.format(domain[0], domain[1]))

                    line.append(domain[3])
                    file.write('{0}\n'.format('\t'.join(line)))
            file.close()
            return True

        except IOError:
            return False

    def write_grn(self):
        """
        Writes the GRN and the sub-network of interest into two .gml files.
        """
        with open(cfg.op.grn_full, 'w', encoding='latin-1') as file, open(cfg.op.grn_int, 'w') as int_file:
            # write the graph specifications for both networks
            file.write('graph [\n')
            file.write('\tdirected 1\n')
            file.write('\tgraphics [\n\t]\n')
            file.write('\tLabelGraphics [\n\t]\n')

            int_file.write('graph [\n')
            int_file.write('\tdirected 1\n')
            int_file.write('\tgraphics [\n\t]\n')
            int_file.write('\tLabelGraphics [\n\t]\n')

            # set to make sure that each gene name can only be used once as a node label
            gene_names = set()

            # iterate over all genes in the gene dictionary
            for lt, gene in self.genes.items():
                # if the gene has no gene name the node name is the locus tag
                if gene.gene_name in cfg.misc.empty:
                    label = gene.locus_tag

                # if the gene name has a gene name, use the gene name as label
                else:
                    # labels have to be unique, therefore add a number to the gene name if it belongs to multiple genes
                    counter = 1
                    label = gene.gene_name
                    while label in gene_names:
                        label = '{0}_{1}'.format(gene.gene_name, counter)
                        counter += 1
                    gene_names.add(label)

                # add the number of non-synonymous mutations and number of mutations in the promoter and transcription
                # factor binding sites to the label
                label = '{0} ({1},{2},{3})'.format(label, gene.non_syn_nbr, gene.prom_nbr, gene.tfbs_nbr)

                node = '\tnode [\n' \
                       '\t\tid {0}\n' \
                       '\t\tlabel \"{1}\"\n' \
                       '\t\tcategory_of_interest \"{2}\"\n' \
                       '\t]\n'.format(gene.grn_id, label, gene.interest)

                # write the node into the corresponding GRN file(s)
                if lt in self.interest_lts or lt in self.reg_lts or lt in self.op_lts:
                    int_file.write(node)
                file.write(node)

            # write the edges into the GRN file
            for edge in self.all_edges:
                file.write(edge.gml())

            # write all edges of the sub-network of interest into the corresponding .gml file
            for edge in self.interest_edges:
                int_file.write(edge.gml())

            # write the closing tag of both files
            file.write(']')
            int_file.write(']')

        if not cfg.run.a_int:
            os.remove(cfg.op.grn_int)

    """
    ANALYSIS
    """

    def analyse_coding_region(self):
        """
        Computes the effect on the amino acid sequence, as well as the amino acid substitution score, of mutations in
        coding regions.
        :return: dictionary with the substitution scores of the mutations
        """
        # initialise substitution matrix
        matrix = SubstitutionMatrix(cfg.user.ai_sm)

        if matrix.error:
            self.coding = False
            messagebox.showerror('Substitution Matrix Error',
                                 matrix.error)
            return

        # initialise score dictionary
        scores = {cfg.user.as_int: {cfg.misc.missense: [],
                                    cfg.misc.nonsense: [],
                                    cfg.misc.readthrough: [],
                                    cfg.misc.frame_shift: [],
                                    cfg.misc.syn: []},
                  cfg.user.as_not_int: {cfg.misc.missense: [],
                                        cfg.misc.nonsense: [],
                                        cfg.misc.readthrough: [],
                                        cfg.misc.frame_shift: [],
                                        cfg.misc.syn: []}}

        # iterate over all coding region mutations in all genes
        for lt, gene in self.genes.items():
            # skip genes that have have no DNA or protein sequence
            if gene.gene_length <= 0:
                continue

            for mut_id in gene.coding_mut:
                # get the mutation from the database
                if mut_id in self.mutations.keys():
                    mut = self.mutations[mut_id]
                else:
                    cfg.a_log.db_find_mut.add(mut_id)
                    continue

                # mutation reference and alternative bases are given on the plus strand, therefore they have to be
                # adjusted if the gene is on the minus strand
                if gene.strand == cfg.misc.minus:
                    ref = complement_strand(mut.ref)
                    alt = complement_strand(mut.alt)
                else:
                    ref = mut.ref
                    alt = mut.alt

                # compute the position of the mutation in the gene
                if gene.strand == cfg.misc.minus:
                    pos = abs(mut.pos - gene.gene_end)
                else:
                    pos = abs(mut.pos - gene.gene_start)
                # compute the effect on the amino acid sequence and the amino acid substitution score
                ref_amino, alt_amino, effect, score = self._coding_effect(gene.dna, gene.prot, ref, alt, pos, matrix)

                # add the information to the mutation and the gene
                mut.add_coding(gene.gene_start, gene.gene_end, lt,
                               ref_amino, alt_amino, effect, score)
                gene.add_mutation_effect(mut_id, effect)

                # add the score to the score dictionary
                if gene.interest in cfg.cat:
                    scores[cfg.user.as_int][effect].append(score)
                else:
                    scores[cfg.user.as_not_int][effect].append(score)

        return scores

    @staticmethod
    def _coding_effect(dna_ref, prot_ref, ref, alt, pos, matrix):
        """
        Computes the effect a mutation has on the amino acid sequence of the protein. Supported effects:
        - synonymous: no change in amino acid sequence, the protein is the same
        - nonsense: an amino acid is replaced by a stop codon, the protein is shorter
        - missense: an amino acid is replaced by a different amino acid, the protein is different in one amino acid
        - readthrough: a stop codon is replaced by an amino acid, the protein is longer
        - reading-box shift: the reading box changes due to an indel, the protein is completely different
        :param dna_ref: reference DNA sequence
        :param prot_ref: reference protein amino acid sequence
        :param ref: reference DNA base(s)
        :param alt: alternative DNA base(s)
        :param pos: gene position of the mutation in the gene
        :param matrix: amino acid substitution matrix
        :return: reference amino acid, alternative amino acid, effect on amino acid sequence, substitution score
        """
        if prot_ref in cfg.misc.empty:
            return cfg.misc.none, cfg.misc.none, cfg.misc.syn, 1.0

        # compute the protein position of the mutation
        pos_prot = pos // 3
        prot_ref_len = len(prot_ref)

        # mutation outside the protein => synonymous
        if pos_prot > prot_ref_len:
            return cfg.misc.none, cfg.misc.none, cfg.misc.syn, 1.0

        # insert the mutation into the DNA sequence and compute the corresponding, mutated protein
        dna_mut = dna_ref[:pos] + alt + dna_ref[pos + len(ref):]
        prot_mut = translate(dna_mut)
        prot_mut_len = len(prot_mut)

        # change not divisible by 3 => reading box shift
        if (abs(len(ref) - len(alt)) % 3) != 0:
            ref_amino = cfg.misc.unknown
            alt_amino = cfg.misc.unknown
            effect = cfg.misc.frame_shift
        else:
            # mutated protein longer => stop coding replaced => readthrough
            if prot_mut_len > prot_ref_len:
                ref_amino = cfg.misc.stop
                alt_amino = prot_mut[prot_ref_len]
                effect = cfg.misc.readthrough
            # mutated protein shorter => new stop codon => nonsense
            elif prot_mut_len < prot_ref_len:
                ref_amino = prot_ref[prot_mut_len]
                alt_amino = cfg.misc.stop
                effect = cfg.misc.nonsense
            else:
                # proteins are the same => synonymous
                if prot_mut == prot_ref:
                    # synonymous mutation in the stop codon
                    if prot_ref_len == pos_prot:
                        ref_amino = cfg.misc.stop
                        alt_amino = cfg.misc.stop
                    # synonymous mutation in an amino acid
                    else:
                        ref_amino = prot_ref[pos_prot]
                        alt_amino = prot_mut[pos_prot]
                    return ref_amino, alt_amino, cfg.misc.syn, 1.0
                # same length but different => missense
                else:
                    ref_amino = prot_ref[pos_prot]
                    alt_amino = prot_mut[pos_prot]
                    effect = cfg.misc.missense

        # compute the mutation scores
        if effect == cfg.misc.missense:
            score = matrix.missense_score(ref_amino, alt_amino)
        # additionally perform global pairwise sequence alignment
        elif effect == cfg.misc.frame_shift:
            score = matrix.score(prot_ref, prot_mut, pos_prot)
        # readthrough needs to be the other way around to avoid negative scores
        elif effect == cfg.misc.readthrough:
            score = matrix.score(prot_mut, prot_ref)
        # no alignment needed
        else:
            score = matrix.score(prot_ref, prot_mut)

        return ref_amino, alt_amino, effect, score

    def analyse_tfbs(self):
        """
        Uses position weight matrices of transcription factors to compute the scores of mutations in transcription
        factor binding sites and compares them to the scores of randomly generated mutations.
        :return: score dictionary
        """
        scores = {cfg.user.as_int: {'actual': [], 'random': []},
                  cfg.user.as_not_int: {'actual': [], 'random': []}}

        for tfbs, values in self.tfbs.items():
            start, end = tfbs
            tf_lts = values['TF']
            tfbs_lts = values['TFBS']
            mutations = values['mutations']
            seqs = values['DNA']
            strand = values['strand']

            # check categories of interest of all associated genes
            inter = tfbs_lts.intersection(self.interest_lts)
            if len(inter) == len(tfbs_lts):
                interest = True
                non = False
            elif len(inter) > 0:
                interest = True
                non = True
            else:
                non = True
                interest = False

            for mut_id in mutations:
                # find the mutations
                mut = self.by_id(mut_id)
                if not mut:
                    continue

                # skip indels
                if mut.type not in [cfg.misc.transition, cfg.misc.transversion]:
                    continue

                # check again that the mutation is in the binding site
                if not (start <= mut.pos <= end):
                    cfg.a_log.db_tfbs_mutation_wrong.add(mut.ID)
                    continue

                for dna in seqs:
                    for tf_lt in tf_lts:
                        tf = self.by_lt(tf_lt)  # type: Gene

                        # not transcription factor
                        if not tf:
                            continue

                        # no PWM for the transcription factor
                        if not tf.tf_pwm:
                            cfg.a_log.db_tfbs_no_tf_pwm.add(tf_lt)
                            continue

                        # compute the score of the mutation and of the random mutations
                        actual, random = self._tfbs_score(start, end, strand, dna, mut, tf.tf_pwm)

                        # add the score to the mutation
                        mut.add_tfbs(start, end, actual)

                        if interest:
                            scores[cfg.user.as_int]['random'] += random
                            scores[cfg.user.as_int]['actual'].append(actual)
                        if non:
                            scores[cfg.user.as_not_int]['random'] += random
                            scores[cfg.user.as_not_int]['actual'].append(actual)

        return scores

    @staticmethod
    def _tfbs_score(start, end, strand, dna, mut, pwm):
        """
        Computes the transcription factor binding site score of a mutation, generates a number of random mutations in
        the binding site and scores them as well using the position weight matrix of the transcription factor.
        :param start: absolute start transcription factor binding site
        :param dna: DNA sequence of the transcription factor binding site
        :param mut: mutation
        :param pwm: position weight matrix of the transcription factor
        :return: score of the actual mutation, scores of random mutations
        """
        if strand == cfg.misc.minus:
            pos = abs(mut.pos - end)
            alt = complement_strand(mut.alt)
        else:
            pos = abs(mut.pos - start)
            alt = mut.alt
        scores_rnd = []

        # score the reference binding site sequence
        score_ref = pwm.alignment_score(dna)

        # score the mutated binding site
        seq_mut = dna[:pos] + alt + dna[pos + len(mut.ref):]
        score_mut = score_ref - pwm.alignment_score(seq_mut)

        # generate and score random mutations in the binding site
        for _ in range(0, cfg.user.as_rnd):
            # get a random position in the binding site
            pos = randint(0, len(dna) - 1)

            # preserve the type of the actual mutation
            if mut.type == cfg.misc.transition:
                if dna[pos] == 'A':
                    alt = 'G'
                elif dna[pos] == 'G':
                    alt = 'A'
                elif dna[pos] == 'T':
                    alt = 'C'
                elif dna[pos] == 'C':
                    alt = 'T'
            elif mut.type == cfg.misc.transversion:
                if dna[pos] == 'A':
                    alt = choice(['C', 'T'])
                elif dna[pos] == 'C':
                    alt = choice(['A', 'G'])
                elif dna[pos] == 'G':
                    alt = choice(['C', 'T'])
                elif dna[pos] == 'T':
                    alt = choice(['A', 'G'])

            # compute and score the randomly mutated binding site
            seq_rnd = dna[:pos] + alt + dna[pos + 1:]
            score_rnd = pwm.alignment_score(seq_rnd)
            scores_rnd.append(score_ref - score_rnd)

        return score_mut, scores_rnd

    def build_grn(self):
        """
        Builds the entire gene regulatory network, including the regulatory sub-network of interest, and outputs
        them in .gml format.
        """
        # keep going over the genes and their regulation as long as changes are still occurring
        change = True

        while change:
            change = False

            for locus_tag, gene in self.genes.items():
                # store locus tags of genes of interest
                if gene.interest in self.categories:
                    self.interest_lts.add(locus_tag)

                # add the regulated genes to the network
                for lt, reg in gene.regulated.items():
                    regulated = self.by_lt(lt)
                    if not regulated:
                        cfg.a_log.db_find_gene.add(lt)
                        continue

                    # add the edge to the regulated gene to the network
                    edge = Edge(0, gene.grn_id, regulated.grn_id, reg['regulation'])
                    self.all_edges.add(edge)

                    if locus_tag in self.interest_lts:
                        # if the regulator and the regulated gene are part of the sub-network, add the regulatory
                        # edge to the sub-network
                        if lt in self.interest_lts or lt in self.reg_lts:
                            self.interest_edges.add(edge)
                    elif locus_tag in self.reg_lts:
                        # if the regulator is a gene that regulates genes in the regulatory sub-network of interest
                        # and the regulated gene is a gene of interest or also regulating genes
                        # in the regulatory sub-network, add the regulatory edge to the sub-network
                        if lt in self.interest_lts or lt in self.reg_lts:
                            self.interest_edges.add(edge)
                    # if the regulator is neither, but the regulated gene is a gene of interest or a gene
                    # regulating genes in the regulatory sub-network of interest, then add the regulator
                    # and the regulatory edge to the sub-network, and also add the regulator to the set of genes
                    # regulating genes in the sub-network. also mark that the set of genes in the sub-network
                    # changed
                    elif lt in self.interest_lts or lt in self.reg_lts:
                        self.interest_edges.add(edge)
                        self.reg_lts.add(locus_tag)
                        change = True

        # adds operon information to the gene regulatory network
        self._add_operons_to_grn()

        # set edge IDs
        for edge in self.all_edges:
            edge.ID = self.next_edge_id
            self.next_edge_id += 1

        self.next_edge_id = 1
        for edge in self.interest_edges:
            edge.ID = self.next_edge_id
            self.next_edge_id += 1

        # write the
        self.write_grn()

    def _add_operons_to_grn(self):
        """
        Adds edges for operons to the gene regulatory network.
        """
        for operon in self.operons:
            # determine if the operon contains a gene of interest
            if self.interest_lts.intersection(operon):
                interest = True
            else:
                interest = False

            # add the operon edge to the grn
            for i in range(0, len(operon) - 1):
                source = self.genes[operon[i]].grn_id
                target = self.genes[operon[i + 1]].grn_id
                edge = Edge(0, source, target, cfg.misc.operon)

                if interest:
                    self.interest_edges.add(edge)
                    self.op_lts.add(operon[i])
                    self.op_lts.add(operon[i + 1])
                self.all_edges.add(edge)
