# standard or third party modules
from collections import defaultdict

# import own modules
from source.configuration import cfg
from source.tools import is_dna_sequence

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2017 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'


class GeneRegion:
    def __init__(self, id, pos, lt, region, interest):
        self.id = id                        # type: int
        self.pos = pos                      # type: int

        self.locus_tags = {lt}              # type: set
        
        self.protein_domains = set()        # type: set

        # transcription factor binding site score
        self.tfbs_score = cfg.misc.none     # type: str

        # effect of the mutation on the coding region of the gene
        self.coding_ref = cfg.misc.none     # type: str
        self.coding_alt = cfg.misc.none     # type: str
        self.coding_effect = cfg.misc.none  # type: str
        self.coding_score = cfg.misc.none   # type: str

        # flags
        self.prom = region == cfg.misc.prom         # type: bool
        self.tfbs = region == cfg.misc.tfbs         # type: bool
        self.coding = region == cfg.misc.coding     # type: bool
        self.int = interest in cfg.cat              # type: bool
        self.non = not self.int                     # type: bool

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    """ TO STRING """
    def get_region(self):
        """
        :return: string representation of the region information
        """
        if self.coding:
            return cfg.misc.coding
        if self.tfbs:
            return cfg.misc.tfbs
        if self.prom:
            return cfg.misc.prom
        return cfg.misc.none_out

    def get_tfbs_score(self):
        """
        :return: string representation of the TFBS score
        """
        if not self.tfbs or self.tfbs_score == cfg.misc.none:
            return cfg.misc.none_out
        return '{:.5f}'.format(round(self.tfbs_score, 5))
    
    def get_substitution(self):
        """
        :return: string representation of the amino acid substitution
        """
        if not self.coding or self.coding_alt == cfg.misc.none or self.coding_ref == cfg.misc.none:
            return cfg.misc.none_out
        return '{0} --> {1}'.format(self.coding_ref, self.coding_alt)
    
    def get_effect(self):
        """
        :return: string representation of the amino acid substitution effect
        """
        if not self.coding or self.coding_effect == cfg.misc.none:
            return cfg.misc.none_out
        return self.coding_effect
    
    def get_substitution_score(self):
        """
        :return: string representation of the amino acid substitution score
        """
        if not self.coding or self.coding_score == cfg.misc.none:
            return cfg.misc.none_out
        return '{:.5f}'.format(round(self.coding_score, 5))
    
    def get_protein_domains(self):
        """
        :return: string representation of the protein domains
        """
        if not self.protein_domains or not self.coding:
            return cfg.misc.none_out
        return cfg.misc.s_sep.join(['{0} [{1}]'.format(dom[2], dom[3]) for dom in sorted(self.protein_domains)])
        
    def get_interest(self):
        """
        :return: string representation of the of interest information
        """
        if self.int and self.non:
            return cfg.misc.both
        if self.int:
            return cfg.misc.int
        if self.non:
            return cfg.misc.non_int
        
    def get_locus_tags(self):
        """
        :return: string representation of the locus tags
        """
        lts = sorted(self.locus_tags)
        if not lts:
            return cfg.misc.none_out
        if self.coding or len(lts) == 1:
            return lts[0]
        return '({})'.format(cfg.misc.s_sep.join(lts))
        

class Mutation:
    """
    Class that represents a mutation and stores information about the type, reference and alternative DNA base, and
    basic information about the gene(s) in which it occurs.
    """
    def __init__(self):
        self.ID = 0                                 # type: int
        self.type = cfg.misc.none_out               # type: str
        self.pos = 0                                # type: int
        self.chromosome = cfg.misc.none_out         # type: str
        self.ref = cfg.misc.none_out                # type: str
        self.alt = cfg.misc.none_out                # type: str

        self.regions = dict()
        self.next_region_id = 1
        self.region_ids_by_pos = {cfg.misc.prom: dict(), cfg.misc.coding: dict(), cfg.misc.tfbs: dict()}
        self.region_ids_by_lt = defaultdict(set)

    def __eq__(self, other):
        return other.ID == self.ID

    def __hash__(self):
        return self.ID

    """ SET MUTATION INFORMATION """
    def set_pos(self, pos):
        """
        Sets the genome position of the mutation.
        :param pos: 1-based genome position
        """
        # convert pos to int
        try:
            pos = int(pos)
        except ValueError:
            cfg.a_log.mut_pos_int.add(self.ID)
            return False

        # the position is 1-based
        if pos <= 0:
            cfg.a_log.mut_pos_zero.add(self.ID)
            return False

        self.pos = pos
        return True

    def set_chromosome(self, chromosome):
        """
        Sets the chromosome of the mutation
        :param chromosome: chromosome
        """
        if chromosome not in cfg.misc.empty:
            self.chromosome = chromosome

    def set_ref_and_alt(self, ref, alt):
        """
        Sets the reference and alternative DNA base(s) of a mutation, as well as the mutation type.
        :param ref: reference base(s)
        :param alt: alternative base(s)
        """
        ref = ref.upper()
        alt = alt.upper()

        if not is_dna_sequence(ref) or not is_dna_sequence(alt):
            cfg.a_log.mut_type_dna.add(self.ID)
            return False

        # parse the mutation type
        if {ref, alt} == {'A', 'C'}:
            t = cfg.misc.transversion
        elif {ref, alt} == {'G', 'T'}:
            t = cfg.misc.transversion
        elif {ref, alt} == {'A', 'G'}:
            t = cfg.misc.transition
        elif {ref, alt} == {'C', 'T'}:
            t = cfg.misc.transition
        elif {ref, alt} == {'A', 'T'}:
            t = cfg.misc.transversion
        elif {ref, alt} == {'C', 'G'}:
            t = cfg.misc.transversion
        else:
            t = cfg.misc.indel

        self.ref = ref
        self.alt = alt
        self.type = t

        return True

    """ ADD ADDITIONAL INFORMATION """
    def add_gene_region(self, start, end, lt, interest, region):
        """
        Adds information on a gene in which the mutation occurs.
        :param start: absolute start of the gene region
        :param end: absolute end of the gene region
        :param lt: locus tag of the gene
        :param interest: gene of interest information
        :param region: gene region
        """
        # check if the region is supported
        if region not in [cfg.misc.tfbs, cfg.misc.coding, cfg.misc.prom]:
            cfg.a_log.mut_region_unsupported.add((self.ID, region))
            return

        # add the new gene region
        if (start, end) not in self.region_ids_by_pos[region]:
            self.regions[self.next_region_id] = GeneRegion(self.next_region_id, (start, end), lt, region, interest)
            self.region_ids_by_pos[region][(start, end)] = self.next_region_id
            self.region_ids_by_lt[lt].add(self.next_region_id)
            self.next_region_id += 1
        # make sure that a coding region is only associated with one locus tag
        elif region == cfg.misc.coding:
            if lt not in self.regions[self.region_ids_by_pos[region][(start, end)]].locus_tags:
                cfg.a_log.mut_coding_region.add((self.pos, lt))
                return
        # add the new information to existing promoter or TFBS information
        else:
            r = self.regions[self.region_ids_by_pos[region][(start, end)]]
            r.locus_tags.add(lt)
            r.int = r.int or (interest in cfg.cat)
            r.non = r.non or (interest not in cfg.cat)
            self.region_ids_by_lt[lt].add(r.id)

    def add_coding(self, start, end, lt, ref, alt, effect, score):
        """
        Adds information on the coding region of the gene and the mutation's effect.
        :param start: absolute coding region start
        :param end: absolute coding region end
        :param lt: locus tag of the gene
        :param ref: reference amino acid
        :param alt: alternative amino acid
        :param effect: effect on the amino acid sequence
        :param score: substitution score of the mutation
        """
        if (start, end) not in self.region_ids_by_pos[cfg.misc.coding]:
            cfg.a_log.mut_coding_not.add((self.pos, lt))
            return

        self.regions[self.region_ids_by_pos[cfg.misc.coding][(start, end)]].coding_ref = ref
        self.regions[self.region_ids_by_pos[cfg.misc.coding][(start, end)]].coding_alt = alt
        self.regions[self.region_ids_by_pos[cfg.misc.coding][(start, end)]].coding_effect = effect
        self.regions[self.region_ids_by_pos[cfg.misc.coding][(start, end)]].coding_score = score

    def add_tfbs(self, start, end, score):
        """
        Adds transcription factor binding site score.
        :param start: absolute TFBS start
        :param end: absolute TFBS end
        :param lt: locus tag of the gene with the binding site
        :param score: score of the binding site
        :param start: start of the binding site for sorting
        """
        if (start, end) not in self.region_ids_by_pos[cfg.misc.tfbs]:
            cfg.a_log.mut_tfbs_not.add(self.pos)
            return

        self.regions[self.region_ids_by_pos[cfg.misc.tfbs][(start, end)]].tfbs_score = score

    def add_protein_domain(self, start, end, lt, domains):
        """
        Adds protein domains.
        :param start: absolute coding region start
        :param end: absolute coding region end
        :param lt: locus tag of the gene with the domains
        :param domains: list with the domains
        """
        if (start, end) not in self.region_ids_by_pos[cfg.misc.coding]:
            cfg.a_log.mut_coding_not.add((self.pos, lt))
            return

        self.regions[self.region_ids_by_pos[cfg.misc.coding][(start, end)]].protein_domains = set(domains)

    """ TO STRING """
    def get_locus_tags(self):
        """
        Returns the string representation of the locus tag(s) of the gene(s) in which the mutation occurs.
        """
        temp = sorted(self.regions.values(), key=lambda x: x.pos)

        if not temp:
            return cfg.misc.none_out
        return cfg.misc.p_sep.join([r.get_locus_tags() for r in temp])

    def get_interest(self):
        """
        Returns the string representation of the of interest information of the gene(s) in which the mutation occurs.
        """
        temp = sorted(self.regions.values(), key=lambda x: x.pos)

        if not temp:
            return cfg.misc.none_out
        return cfg.misc.p_sep.join([r.get_interest() for r in temp])

    def get_tfbs_scores(self, lt=''):
        """
        Returns the string representation of th scores of the transcription factor binding site(s) in which
        the mutation occurs.
        """
        if not lt:
            temp = sorted(self.regions.values(), key=lambda x: x.pos)
        elif lt in self.region_ids_by_lt:
            temp = sorted([self.regions[i] for i in self.region_ids_by_lt[lt] if self.regions[i].tfbs],
                          key=lambda x: x.pos)
        else:
            temp = []

        if not temp:
            return cfg.misc.none_out
        if lt:
            return cfg.misc.s_sep.join([r.get_tfbs_score() for r in temp])
        return cfg.misc.p_sep.join([r.get_tfbs_score() for r in temp])

    def get_regions(self, lt=''):
        """
        Returns the string representation of the gene region(s) in which the mutation occurs.
        """
        if not lt:
            temp = sorted(self.regions.values(), key=lambda x: x.pos)
            
            if not temp:
                return cfg.misc.intergenic
            return cfg.misc.p_sep.join([r.get_region() for r in temp])
        elif lt in self.region_ids_by_lt:
            temp = set([self.regions[i].get_region() for i in self.region_ids_by_lt[lt]])
            return cfg.misc.s_sep.join(sorted(temp))
        return cfg.misc.none_out

    def get_amino_change(self, lt=''):
        """
        Returns the string representation of the amino acid substitution(s) caused by the mutation.
        """
        if not lt:
            temp = sorted(self.regions.values(), key=lambda x: x.pos)
        elif lt in self.region_ids_by_lt:
            temp = sorted([self.regions[i] for i in self.region_ids_by_lt[lt] if self.regions[i].coding],
                          key=lambda x: x.pos)
        else:
            temp = []

        if not temp:
            return cfg.misc.none_out
        if lt:
            return cfg.misc.s_sep.join([r.get_substitution() for r in temp])

        return cfg.misc.p_sep.join([r.get_substitution() for r in temp])

    def get_substitution_score(self, lt=''):
        """
        Returns the string representation of the amino acid substitution score(s) caused by the mutation.
        """
        if not lt:
            temp = sorted(self.regions.values(), key=lambda x: x.pos)
        elif lt in self.region_ids_by_lt:
            temp = sorted([self.regions[i] for i in self.region_ids_by_lt[lt] if self.regions[i].coding],
                          key=lambda x: x.pos)
        else:
            temp = []

        if not temp:
            return cfg.misc.none_out
        if lt:
            return cfg.misc.s_sep.join([r.get_substitution_score() for r in temp])
        return cfg.misc.p_sep.join([r.get_substitution_score() for r in temp])

    def get_effect(self, lt=''):
        """
        Returns the string representation of the effect of the mutation on the amino acid sequence(s).
        """
        if not lt:
            temp = sorted(self.regions.values(), key=lambda x: x.pos)
        elif lt in self.region_ids_by_lt:
            temp = sorted([self.regions[i] for i in self.region_ids_by_lt[lt] if self.regions[i].coding],
                          key=lambda x: x.pos)
        else:
            temp = []

        if not temp:
            return cfg.misc.none_out
        if lt:
            return cfg.misc.s_sep.join([r.get_effect() for r in temp])
        return cfg.misc.p_sep.join([r.get_effect() for r in temp])

    def get_protein_domains(self, lt=''):
        """
        Returns the string representation of the protein domains in which the mutation occurs.
        """
        if not lt:
            temp = sorted(self.regions.values(), key=lambda x: x.pos)
        elif lt in self.region_ids_by_lt:
            temp = sorted([self.regions[i] for i in self.region_ids_by_lt[lt] if self.regions[i].coding],
                          key=lambda x: x.pos)
        else:
            temp = []

        if not temp:
            return cfg.misc.none_out
        if lt:
            return cfg.misc.s_sep.join([r.get_protein_domains() for r in temp])
        return cfg.misc.p_sep.join([r.get_protein_domains() for r in temp])
