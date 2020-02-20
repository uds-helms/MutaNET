# import standard or third party modules
from collections import defaultdict

from source.tools import parse_strand, is_dna_sequence, is_protein_sequence, parse_regulation, combine_reg_info

# import own modules
from source.configuration import cfg

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2017 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'


class Gene:
    def __init__(self):
        """
        Initialises an "empty" gene.
        """
        # basic gene information
        self.locus_tag = cfg.misc.none_out          # type: str
        self.gene_name = cfg.misc.none_out          # type: str
        self.interest = cfg.misc.none_out           # type: str
        self.description = cfg.misc.none_out        # type: str

        # coding region
        self.chromosome = cfg.misc.none_out     # type: str
        self.strand = cfg.misc.plus             # type: str
        self.gene_start = 0                     # type: int
        self.gene_end = 0                       # type: int
        self.gene_length = 0                    # type: int
        self.dna = cfg.misc.none_out            # type: str
        self.prot = cfg.misc.none_out           # type: str
        self.prot_length = 0                    # type: int

        # promoter region
        self.prom_start_rel = 0         # type: int
        self.prom_start_abs = 0         # type: int
        self.prom_end_abs = 0           # type: int
        self.prom_length = 0            # type: int

        # protein domains
        self.prot_dom = set()       # type: set
        self.mut_dom = set()        # type: set

        # regulation
        self.operon = tuple()                   # type: tuple
        self.regulated = defaultdict(dict)      # type: defaultdict(dict)
        self.regulators = defaultdict(dict)     # type: defaultdict(dict)

        # transcription factor binding sites
        self.tfbs = dict()                  # type: dict

        # transcription factor position weight matrix
        self.tf_pwm = None                  # type: PWM

        # accession numbers
        self.gene_id = cfg.misc.none_out    # type: str
        self.uniprot = cfg.misc.none_out    # type: str
        self.genbank = cfg.misc.none_out    # type: str
        self.pfam = cfg.misc.none_out       # type: str
        self.pubmed = set()                 # type: set

        # mutation IDs
        self.mutations = set()          # type: set
        self.coding_mut = set()         # type: set
        self.syn_mut = set()            # type: set
        self.non_syn_mut = set()        # type: set
        self.missense_mut = set()       # type: set
        self.nonsense_mut = set()       # type: set
        self.readthrough_mut = set()    # type: set
        self.frame_shift_mut = set()    # type: set
        self.prom_mut = set()           # type: set
        self.tfbs_mut = set()           # type: set

        # mutation densities
        self.mutation_density = 0       # type: int
        self.coding_md = 0              # type: int
        self.syn_md = 0                 # type: int
        self.non_syn_md = 0             # type: int
        self.missense_md = 0            # type: int
        self.nonsense_md = 0            # type: int
        self.readthrough_md = 0         # type: int
        self.frame_shift_md = 0         # type: int
        self.prom_md = 0                # type: int
        self.tfbs_md = 0                # type: int

        # number of mutations
        self.mutation_nbr = 0       # type: int
        self.mutation_nbr = 0       # type: int
        self.coding_nbr = 0         # type: int
        self.syn_nbr = 0            # type: int
        self.non_syn_nbr = 0        # type: int
        self.missense_nbr = 0       # type: int
        self.nonsense_nbr = 0       # type: int
        self.readthrough_nbr = 0    # type: int
        self.frame_shift_nbr = 0    # type: int
        self.prom_nbr = 0           # type: int
        self.tfbs_nbr = 0           # type: int

        self.grn_id = 0

    """ SET REQUIRED GENE INFORMATION """
    def set_locus_tag(self, lt):
        """
        Sets the locus tag of the gene.
        :param lt: locus tag
        """
        if lt not in cfg.misc.empty:
            self.locus_tag = lt
            return True
        return False

    def set_strand(self, strand):
        """
        Sets the gene strand.
        :param strand: gene strand
        """
        strand = parse_strand(strand.lower())

        # the strand format is not supported
        if not strand:
            cfg.a_log.gene_strand.add(self.locus_tag)
            return False

        self.strand = strand

        return True

    def set_gene_start_and_end(self, start, end):
        """
        Sets the gene start and end (coding region).
        :param start: gene start
        :param end: gene end
        """
        # one of the input fields is empty
        if (start in cfg.misc.empty) or (end in cfg.misc.empty):
            return False

        # get the integer representation of the parameters
        try:
            start = int(start)
            end = int(end)
        except ValueError:
            cfg.a_log.gene_se_not_int.add(self.locus_tag)
            return False

        # genome positions cannot be negative
        if (start < 0) or (end < 0):
            cfg.a_log.gene_se_less_zero.add(self.locus_tag)
            return False

        # if one of them is 0, set the other to 0 as well
        if (start == 0) or (end == 0):
            cfg.a_log.gene_se_zero.add(self.locus_tag)
            self.gene_start = 0
            self.gene_end = 0
            return True

        # start <= end position
        self.gene_start = min(start, end)
        self.gene_end = max(start, end)

        return True

    def set_dna(self, dna):
        """
        Sets the DNA sequence of the gene.
        :param dna: DNA sequence
        """
        # input field is empty
        if dna in cfg.misc.empty:
            cfg.a_log.gene_dna_empty.add(self.locus_tag)
            return True

        dna = dna.upper()

        # input field is not a DNA sequence
        if not is_dna_sequence(dna):
            cfg.a_log.gene_not_dna.add(self.locus_tag)
            return False

        self.dna = dna
        return True

    def set_coding_length(self):
        """
        Computes and sets the length of the coding region, either using the DNA sequence or the gene start and end.
        """
        if self.dna == cfg.misc.none_out:
            dna_length = 0
        else:
            dna_length = len(self.dna)

        # coding region position not defined
        if (self.gene_start == 0) or (self.gene_end == 0):
            self.gene_length = dna_length
            return True

        # + 1 because the gene end is included
        start_end_length = abs(self.gene_start - self.gene_end) + 1

        # make sure the DNA length matches the coding region length
        if dna_length != start_end_length:
            cfg.a_log.gene_coding_length.add(self.locus_tag)
            return False

        self.gene_length = start_end_length

        return True

    def set_protein(self, prot):
        """
        Sets the amino acid sequence of the protein.
        :param prot: amino acid sequence
        """
        # input field is empty
        if prot in cfg.misc.empty:
            cfg.a_log.gene_protein_empty.add(self.locus_tag)
            return True

        prot = prot.upper()

        # input field is not a prot sequence
        if not is_protein_sequence(prot):
            cfg.a_log.gene_not_protein.add(self.locus_tag)
            return False

        self.prot = prot
        self.prot_length = len(prot)

        return True

    def set_promoter(self, astart, aend, rstart):
        """
        Sets the promoter. If the absolute start and end are given use those to determine the promoter fields. Otherwise
        try the relative start, and if that also fails use the default relative start.
        :param astart: start position relative to the gene start
        :param aend: 1-based end position in the genome
        :param rstart: 1-based start position in the genome
        """
        if not self.gene_start or not self.gene_end:
            return

        # first try to use the absolute start and end
        success = self._set_promoter_from_abs(astart, aend)
        # if that was not successful, try to use the relative start
        if not success:
            success = self._set_promoter_from_rel(rstart)
        # if that was not successful either, try to use the default relative start
        if not success:
            success = self._set_promoter_from_rel(cfg.user.as_prom)
        if not success:
            cfg.a_log.gene_promoter.add(self.locus_tag)

    def _set_promoter_from_rel(self, start):
        """
        Sets the promoter start, end and length based on the start position relative to the gene start.
        :param start: position relative to the gene start
        :return: True if this succeeds, False otherwise
        """
        # input field is empty
        if start in cfg.misc.empty:
            return False

        # convert start to integer
        try:
            start = int(start)
        except ValueError:
            cfg.a_log.gene_prom_not_int.add(self.locus_tag)
            return False

        # must be before the gene
        if start > 0:
            cfg.a_log.gene_prom_pos.add(self.locus_tag)
            return False

        # no promoter region
        if start == 0:
            return True

        self.prom_start_rel = start

        # account for the gene being on the minus strand
        if self.strand == cfg.misc.minus:
            self.prom_end_abs = self.gene_end - self.prom_start_rel
            self.prom_start_abs = self.gene_end + 1
        # account for the gene being on the plus strand
        else:
            self.prom_start_abs = max(1, self.gene_start + self.prom_start_rel)
            self.prom_end_abs = max(1, self.gene_start - 1)

        # set the length of the promoter, + 1 because the end is included
        if self.prom_start_abs == self.prom_end_abs:
            self.prom_length = 0
        else:
            self.prom_length = abs(self.prom_start_abs - self.prom_end_abs) + 1

        return True

    def _set_promoter_from_abs(self, start, end):
        """
        Sets the promoter start, end and length based on the start and end position in the genome.
        :param start: 1-based start position of the promoter in the genome
        :param end: 1-based end position of the promoter in the genome
        :return: True if this succeeds, False otherwise
        """
        # at least one of the input fields is empty
        if (start in cfg.misc.empty) or (end in cfg.misc.empty):
            return False

        # convert start and end to integer
        try:
            start = int(start)
            end = int(end)
        except ValueError:
            cfg.a_log.gene_prom_not_int.add(self.locus_tag)
            return False

        # genome positions cannot be negative or 0, as the positions are 1-based
        if (start <= 0) or (end <= 0):
            cfg.a_log.gene_prom_less_zero.add(self.locus_tag)
            return False

        self.prom_start_abs = min(start, end)
        self.prom_end_abs = max(start, end)

        # account for the gene being on the minus strand
        if self.strand == cfg.misc.minus:
            self.prom_start_rel = self.gene_end - self.prom_end_abs
        # account for the gene being on the plus strand
        else:
            self.prom_start_rel = self.prom_start_abs - self.gene_start

        # compute the promoter length
        self.prom_length = abs(self.prom_start_abs - self.prom_end_abs) + 1
        return True

    """ SET OPTIONAL GENE INFORMATION """
    def set_gene_name(self, name):
        """
        Sets the name of the gene.
        :param name: gene name
        """
        if name not in cfg.misc.empty:
            self.gene_name = name
            return

    def set_chromosome(self, chromosome: str):
        """
        Sets the chromosome.
        :param chromosome: chromosome
        """
        if chromosome not in cfg.misc.empty:
            self.chromosome = chromosome

    def set_description(self, desc):
        """
        Sets the gene description.
        :param desc: description
        """
        if desc not in cfg.misc.empty:
            self.description = desc

    def set_gene_id(self, acc):
        """
        Sets the gene ID of the gene.
        :param acc: gene ID
        """
        if acc not in cfg.misc.empty:
            self.gene_id = acc

    def set_pfam(self, acc):
        """
        Sets the PFAM accession of the gene.
        :param acc: PFAM accession
        """
        if acc not in cfg.misc.empty:
            self.pfam = acc

    def set_uniprot(self, acc):
        """
        Sets the UniProt accession of the gene.
        :param acc: UniProt accession
        """
        if acc not in cfg.misc.empty:
            self.uniprot = acc

    def set_genbank(self, acc):
        """
        Sets the Genbank accession of the gene.
        :param acc: Genbank accession
        """
        if acc not in cfg.misc.empty:
            self.genbank = acc

    def set_pubmed(self, acc):
        """
        Sets the pubmed IDs of the gene.
        :param acc: pubmed IDs (ID1; ID2; ...)
        """
        # input field is empty
        if acc in cfg.misc.empty:
            return

        self.pubmed = set(acc.split(cfg.misc.in_sep))

    """ SET ANTIBIOTIC RESISTANCE INFORMATION """
    def set_interest(self, interest):
        """
        Sets the of interest information of the gene.
        :param interest: of interest
        """
        # input field is empty
        if interest in cfg.misc.empty:
            return

        self.interest = interest

    """ ADD REGULATORY INFORMATION """
    def add_regulated_gene(self, lt, gn, reg):
        """
        Adds a regulated gene to the gene.
        :param lt: locus tag of the regulated gene
        :param gn: gene name of the regulated gene
        :param reg: regulation information
        """
        # parse the regulation information
        reg_processed = parse_regulation(reg)

        # unsupported regulation information
        if (reg_processed == cfg.misc.unknown) and (reg not in cfg.misc.empty):
            cfg.a_log.gene_regulation_unsupported.add((self.locus_tag, lt))

        # add the gene name
        self.regulated[lt]['gene name'] = gn

        # initialise the regulation set for the locus tag
        if 'regulation' not in self.regulated[lt].keys():
            self.regulated[lt]['regulation'] = cfg.misc.unknown

        # remove redundant information
        self.regulated[lt]['regulation'] = combine_reg_info(reg_processed, self.regulated[lt]['regulation'])

    def add_regulator(self, lt, gn, reg):
        """
        Adds a regulator to the gene.
        :param lt: locus tag of the regulator
        :param gn: gene name of the regulator
        :param reg: regulation information
        """
        # parse the regulation information
        reg_processed = parse_regulation(reg)

        # unsupported regulation information
        if (reg_processed == cfg.misc.unknown) and (reg not in cfg.misc.empty):
            cfg.a_log.gene_regulation_unsupported.add((lt, self.locus_tag))

        # add the gene name
        self.regulators[lt]['gene name'] = gn

        # initialise the regulation set for the locus tag
        if 'regulation' not in self.regulators[lt].keys():
            self.regulators[lt]['regulation'] = cfg.misc.unknown

        # combine regulatory information
        self.regulators[lt]['regulation'] = combine_reg_info(reg_processed, self.regulators[lt]['regulation'])

    """ ADD PROTEIN DOMAIN INFORMATION """
    def add_protein_domain(self, start, end, kind, desc):
        """
        Adds the specified protein domain to the list of protein domains.
        :param kind: type of the protein domain (i.e. helix)
        :param start: 1-based start position of the domain in the amino acid sequence
        :param end: 1-based end position of the domain in the amino acid sequence (included)
        :param desc: further description of the domain
        """
        # the input fields are empty
        if (start in cfg.misc.empty) or (end in cfg.misc.empty) or (kind in cfg.misc.empty and desc in cfg.misc.empty):
            return

        # convert start and end to integers
        try:
            start = int(start)
            end = int(end)
        except ValueError:
            cfg.a_log.gene_prot_int.add(self.locus_tag)
            return

        # start or end are not in the protein
        if (start > self.prot_length) or (start <= 0) or (end > self.prot_length) or (end <= 0):
            return

        if desc in cfg.misc.empty:
            desc = cfg.misc.none_out
        if kind in cfg.misc.empty:
            kind = cfg.misc.none_out

        domain = (min(start, end), max(start, end), kind, desc)
        self.prot_dom.add(domain)

        if kind.lower() in ['mut', 'mutagen', 'mutation', 'snp']:
            self.mut_dom.add(domain)

    """ ADD MUTATION INFORMATION """
    def add_mutation_location(self, id, loc):
        """
        Adds a mutation (by its ID) to the gene depending on its location.
        :param id: mutation ID
        :param loc: gene region of the mutation
        """
        # convert the ID to integer
        try:
            ID = int(id)
        except ValueError:
            cfg.a_log.gene_mut_loc_int.add(self.locus_tag)
            return

        # sort by location
        if loc == cfg.misc.coding:
            self.mutations.add(ID)
            self.coding_mut.add(ID)
        elif loc == cfg.misc.prom:
            self.mutations.add(ID)
            self.prom_mut.add(ID)
        elif loc == cfg.misc.tfbs:
            self.mutations.add(ID)
            self.tfbs_mut.add(ID)
        else:
            cfg.a_log.gene_mut_loc_unsupported.add(self.locus_tag)

    def add_mutation_effect(self, id, effect):
        """
        Adds a mutation (by its ID) to the gene depending on its effect on the coding region.
        :param id: mutation ID
        :param effect: effect on the coding region
        """
        # convert the ID to integer
        try:
            ID = int(id)
        except ValueError:
            cfg.a_log.gene_mut_eff_int.add(self.locus_tag)
            return

        self.mutations.add(ID)
        self.coding_mut.add(ID)

        # sort by effect: synonymous
        if effect == cfg.misc.syn:
            self.syn_mut.add(ID)
        # non-synonymous
        elif effect in [cfg.misc.non_syn, cfg.misc.missense, cfg.misc.nonsense, cfg.misc.readthrough,
                        cfg.misc.frame_shift]:
            self.non_syn_mut.add(ID)
            if effect == cfg.misc.missense:
                self.missense_mut.add(ID)
            elif effect == cfg.misc.nonsense:
                self.nonsense_mut.add(ID)
            elif effect == cfg.misc.readthrough:
                self.readthrough_mut.add(ID)
            elif effect == cfg.misc.frame_shift:
                self.frame_shift_mut.add(ID)
            else:
                cfg.a_log.gene_mut_eff_unsupported.add(self.locus_tag)
        else:
            cfg.a_log.gene_mut_eff_unsupported.add(self.locus_tag)

    """ TRANSCRIPTION FACTOR BINDING SITE ANALYSIS """
    def add_tfbs(self, reg, start, end, dna):
        """
        Add a transcription factor binding site to the gene.
        :param reg: locus tag of the transcription factor
        :param start: absolute start of the binding site
        :param end: absolute end of the binding site
        :param dna: DNA sequence of the binding site
        """
        # TFBS already associated with the gene
        if (start, end, dna, reg) in self.tfbs:
            return

        # compute the relative start
        if self.strand == cfg.misc.minus:
            rstart = self.gene_end - end
        else:
            rstart = start - self.gene_start

        self.tfbs[(start, end, dna, reg)] = rstart

    """ COMPUTE STATISTICS """
    def calculate_mutation_densities(self):
        """
        Calculates the mutation densities (# mutations / kbp) for various categories of mutations.
        """
        # calculate the coding and promoter region factors
        fct_coding = self.gene_length / 1000
        fct_prom = self.prom_length / 1000

        # calculate transcription factor binding site factor
        fct_tfbs = 0
        for tfbs in self.tfbs.keys():
            fct_tfbs += len(tfbs[2])
        fct_tfbs /= 1000

        # compute number of mutations
        self.mutation_nbr = len(self.mutations)
        self.coding_nbr = len(self.coding_mut)
        self.syn_nbr = len(self.syn_mut)
        self.non_syn_nbr = len(self.non_syn_mut)
        self.missense_nbr = len(self.missense_mut)
        self.nonsense_nbr = len(self.nonsense_mut)
        self.readthrough_nbr = len(self.readthrough_mut)
        self.frame_shift_nbr = len(self.frame_shift_mut)
        self.prom_nbr = len(self.prom_mut)
        self.tfbs_nbr = len(self.tfbs_mut)

        # calculate the densities
        if fct_coding > 0:
            self.mutation_density = round(self.mutation_nbr / fct_coding, 3)
            self.coding_md = round(self.coding_nbr / fct_coding, 3)
            self.syn_md = round(self.syn_nbr / fct_coding, 3)
            self.non_syn_md = round(self.non_syn_nbr / fct_coding, 3)
            self.missense_md = round(self.missense_nbr / fct_coding, 3)
            self.nonsense_md = round(self.nonsense_nbr / fct_coding, 3)
            self.readthrough_md = round(self.readthrough_nbr / fct_coding, 3)
            self.frame_shift_md = round(self.frame_shift_nbr / fct_coding, 3)

        if fct_prom > 0:
            self.prom_md = round(self.prom_nbr / fct_prom, 3)

        if fct_tfbs > 0:
            self.tfbs_md = round(self.tfbs_nbr / fct_tfbs, 3)

    """ GETTERS """
    def get_protein_domains_for_pos(self, pos):
        """
        Finds and returns all protein domains for the given genome position.
        :param pos: 1-based position in the genome
        :return: all protein domains for the given position
        """
        prot_pos = self.get_protein_position(pos)

        return [domain for domain in self.prot_dom if domain[0] <= prot_pos <= domain[1]]

    def get_protein_position(self, pos):
        """
        Takes a 1-based position in the genome and returns a 1-based position in the amino sequence (if that position
        is in the amino acid sequence to begin with).
        :param pos: 1-based position in the genome
        :return: the protein position, 0 (= not in protein) otherwise
        """
        try:
            pos = int(pos)
        except ValueError:
            cfg.a_log.gene_prot_pos_int.add(self.locus_tag)
            return 0

        # return 0 if the position is negative or equal to zero
        if pos <= 0:
            cfg.a_log.gene_prot_pos_zero.add(self.locus_tag)
            return 0

        if pos < min(self.gene_end, self.gene_start) or pos > max(self.gene_end, self.gene_start):
            return 0

        # compute the 1-based protein position
        if self.strand == cfg.misc.minus:
            pos = (abs(pos - self.gene_end) // 3) + 1
        else:
            pos = (abs(pos - self.gene_start) // 3) + 1

        # return 0 if the position is not in the protein
        if pos > self.prot_length:
            return 0

        return pos

    def get_intervals(self):
        """
        Get all coding region, promoter region and transcription factor binding site intervals.
        :return: intervals
        """
        intervals = set()

        if self.gene_start > 0 and self.gene_end > 0:
            intervals.add((self.gene_start, self.gene_end, self.locus_tag, self.interest, cfg.misc.coding))

        if self.prom_start_abs > 0 and self.prom_end_abs > 0:
            intervals.add((self.prom_start_abs, self.prom_end_abs, self.locus_tag, self.interest, cfg.misc.prom))

        for tfbs in self.tfbs.keys():
            if tfbs[0] > 0 and tfbs[1] > 0:
                intervals.add((tfbs[0], tfbs[1], self.locus_tag, self.interest, cfg.misc.tfbs))

        return intervals

    """ TO STRING """
    def get_locus_tag_and_gene_name(self):
        """
        Returns the string representation of the locus tag and gene name.
        :return: locus_tag (gene_name)
        """
        if self.gene_name == cfg.misc.none_out:
            return self.locus_tag
        return '{0} ({1})'.format(self.locus_tag, self.gene_name)

    def get_operon(self):
        """
        Returns the string representation of the operon.
        :return: gene1 > gene2 > gene3 > ...
        """
        if self.operon:
            return cfg.misc.operon_fwd.join(self.operon)
        return cfg.misc.none_out

    def get_regulator_by_locus_tag(self, lt):
        """
        Returns the string representation of the regulation information of a given regulator.
        :param lt: locus tag of the regulator
        :return: 'locus_tag (gene_name) [regulation]'
        """
        # no regulatory information for the locus tag
        if lt not in self.regulators.keys():
            return cfg.misc.none_out

        # shorthands
        name = self.regulators[lt]['gene name']
        reg = self.regulators[lt]['regulation'][0].upper()

        # format output
        if name != cfg.misc.none_out:
            return '{0} ({1}) [{2}]'.format(lt, name, reg)
        return '{0} [{1}]'.format(lt, reg)

    def get_regulators(self):
        """
        Returns string representation of all regulators and the corresponding regulatory information.
        :return: 'locus_tag (gene_name) [regulation]' | 'locus_tag (gene_name) [regulation]'
        """
        if not self.regulators.keys():
            return cfg.misc.none_out
        return cfg.misc.p_sep.join([self.get_regulator_by_locus_tag(lt) for lt in sorted(self.regulators.keys())])

    def get_regulated_genes_by_locus_tag(self, lt):
        """
        Returns the string representation of the regulation information of a given regulated gene.
        :param lt: locus tag of the regulated gene
        :return: 'locus_tag (gene_name) [regulation]'
        """
        # no regulatory information for the locus tag
        if lt not in self.regulated.keys():
            return cfg.misc.none_out

        # shorthand
        name = self.regulated[lt]['gene name']
        if name != cfg.misc.none_out:
            return '{0} ({1}) [{2}]'.format(lt, name, self.regulated[lt]['regulation'][0].upper())
        return '{0} [{1}]'.format(lt, self.regulated[lt]['regulation'][0].upper())

    def get_regulated_genes(self):
        """
        Returns string representation of all regulated genes and the corresponding regulatory information.
        :return: 'locus_tag (gene_name) [regulation]' | 'locus_tag (gene_name) [regulation]'
        """
        if not self.regulated.keys():
            return cfg.misc.none_out
        return cfg.misc.p_sep.join([self.get_regulated_genes_by_locus_tag(lt) for lt in sorted(self.regulated.keys())])

    def get_mutation_ids(self):
        """
        Returns the string representation of the mutation IDS.
        :return: 'ID1 | ID2 | ...'
        """
        if not self.mutations:
            return cfg.misc.none_out
        return cfg.misc.p_sep.join([str(x) for x in sorted(self.mutations)])

    def get_prot_dom_start(self):
        """
        Returns the string representation of all protein domain start positions.
        :return: 'start1 | start2 | ...'
        """
        if not self.prot_dom:
            return cfg.misc.none_out
        return cfg.misc.p_sep.join([str(domain[0]) for domain in sorted(self.prot_dom, key=lambda key: key[2])])

    def get_prot_dom_end(self):
        """
        Returns the string representation of all protein domain end positions.
        :return: 'end1 | end2 | ...'
        """
        if not self.prot_dom:
            return cfg.misc.none_out
        return cfg.misc.p_sep.join([str(domain[1]) for domain in sorted(self.prot_dom, key=lambda key: key[2])])

    def get_prot_dom_type(self):
        """
        Returns the string representation of all protein domain types.
        :return: 'type1 | type2 | ...'
        """
        if not self.prot_dom:
            return cfg.misc.none_out
        return cfg.misc.p_sep.join([domain[2] for domain in sorted(self.prot_dom, key=lambda key: key[2])])

    def get_prot_dom_description(self):
        """
        Returns the string representation of all protein domain description.
        :return: 'description1 | description2 | ...'
        """
        if not self.prot_dom:
            return cfg.misc.none_out
        return cfg.misc.p_sep.join([domain[3] for domain in sorted(self.prot_dom, key=lambda key: key[2])])

    def get_tfbs_seq(self):
        """
        Returns the string representation of the transcription factor binding site sequences sorted by the
        absolute start position of the binding site.
        :return: 'sequence1 | sequence2 | ...'
        """
        if not self.tfbs:
            return cfg.misc.none_out

        return cfg.misc.p_sep.join([t[2] for t in sorted(self.tfbs.keys())])

    def get_tfbs_end(self):
        """
        Returns the string representation of the transcription factor binding site end positions sorted by the
        absolute start position of the binding site.
        :return: 'end1 | end2 | ...'
        """
        if not self.tfbs:
            return cfg.misc.none_out
        return cfg.misc.p_sep.join([str(t[1]) for t in sorted(self.tfbs.keys())])

    def get_tfbs_start_rel(self):
        """
        Returns the string representation of the transcription factor binding site relative start position sorted by the
        absolute start position of the binding site.
        :return: 'start_rel1 | start_rel2 | ...'
        """
        if not self.tfbs:
            return cfg.misc.none_out
        return cfg.misc.p_sep.join([str(self.tfbs[t]) for t in sorted(self.tfbs.keys())])

    def get_tfbs_start_abs(self):
        """
        Returns the string representation of the transcription factor binding site absolute start positions sorted by
        the start position of the binding site.
        :return: 'start_abs1 | start_abs2 | ...'
        """
        if not self.tfbs:
            return cfg.misc.none_out
        return cfg.misc.p_sep.join([str(t[0]) for t in sorted(self.tfbs.keys())])

    def get_tfs(self):
        """
        Returns the string representation of the transcription factor locus tags.
        :return: 'tf1 | tf2 | ...'
        """
        if not self.tfbs:
            return cfg.misc.none_out
        return cfg.misc.p_sep.join([str(t[3]) for t in sorted(self.tfbs.keys())])

    def get_pubmed(self):
        """
        Returns the string representation of the pubmed IDs.
        :return: 'ID1 | ID2 | ...'
        """
        if not self.pubmed:
            return cfg.misc.none_out
        return cfg.misc.p_sep.join(sorted(self.pubmed))
