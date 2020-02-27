# import standard or third party modules
import os
import yaml
import traceback
from shutil import which
from sys import platform
from tkinter import messagebox

# import own modules
from source.log import ALog, NGSLog

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2017 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'


""" HELPER FUNCTIONS """


def find_exe(name, path):
    """
    Tries to find the specified executable by first checking if the given path points to a valid executable. If that
    fails, try to see if the executable name on the path is in the system environment. If that also fails,
    try to see if the name is in the system environment. If all that fails, return the given input path.
    :param name: name of the executable, i.e. varscan
    :param path: given path of the executable
    :return: actual path to the executable
    """
    # test if the given path points to a valid executable and if so, return it
    given_exe_path = which(path)

    if not given_exe_path:
        # test if executable at the end of the given path is in the system environment and if so, return it
        ending_exe_path = which(adjust_file_path(path).split('/')[-1])
        if not ending_exe_path:
            # test if the given name of the executable is in the system environment and if so, return it
            name_exe_path = which(name)
            if not name_exe_path:
                # couldn't find the executable, thus return the path
                return adjust_file_path(path)

            return adjust_file_path(name_exe_path)

        return adjust_file_path(ending_exe_path)

    return adjust_file_path(given_exe_path)


def adjust_dir_path(path):
    """
    Replaces backlashes with forward slashes and adds a forward slash at the end if missing.
    :param path: directory path
    :return: adjusted directory path
    """
    if not path:
        return path

    path = path.replace('\\', '/')

    if path[-1] != '/':
        path += '/'

    return path


def adjust_file_path(path):
    """
    Replaces backlashes with forward slashes.
    :param path: file path
    :return:  adjusted file path
    """
    return path.replace('\\', '/')


""" MISCELLANEOUS CONSTANTS """


class Misc:
    """
    Mutation type, mutation effect, regulation, resistance, separator, strand constants.
    """
    # mutation type
    transition = 'transition'                           # type: str
    transversion = 'transversion'                       # type: str
    indel = 'indel'                                     # type: str
    # mutation effect
    missense = 'missense'                               # type: str
    nonsense = 'nonsense'                               # type: str
    readthrough = 'readthrough'                         # type: str
    frame_shift = 'frame shift'                         # type: str
    non_syn = 'non-synonymous'                          # type: str
    syn = 'synonymous'                                  # type: str
    prom = 'promoter'                                   # type: str
    tfbs = 'TFBS'                                       # type: str
    coding = 'coding'                                   # type: str
    intergenic = 'intergenic'                           # type: str
    stop = 'stop'                                       # type: str
    # regulation
    effector = 'effector'                               # type: str
    activator = 'activator'                             # type: str
    repressor = 'repressor'                             # type: str
    unknown = '?'                                       # type: str
    operon = 'operon'                                   # type: str
    # strand
    plus = '+'                                          # type: str
    minus = '-'                                         # type: str
    # empty
    none = ''                                           # type: str
    none_out = '-'                                      # type: str
    empty = ['', '-', '?', '/', '*', none, none_out]    # type: list
    # separators
    p_sep = '  |  '                                     # type: str
    s_sep = ' / '                                       # type: str
    in_sep = '; '                                       # type: str
    operon_fwd = ' > '                                  # type: str
    operon_bwd = ' < '                                  # type: str
    # genes of interest
    both = 'both'                                       # type: str
    int = 'yes'                                         # type: str
    non_int = 'no'                                      # type: str


""" GUI CONFIGURATION """


class GuiEntity:
    """
    Represents an entity in the GUI: name + label + description
    """
    def __init__(self, label='', desc=''):
        self.label = label      # type: str
        self.desc = desc        # type: str


class GuiKeys:
    """
    Constant keys for the GUI.
    """
    """ ANALYSIS - MAIN """
    # section and subsection titles
    ah_title = 'Mutation Analysis'                          # type: str
    ah_int = 'Resistance Analysis'                          # type: str
    ah_coding = 'Coding Region Analysis'                    # type: str
    ah_prot = 'Protein Domain Analysis'                     # type: str
    ah_tfbs = 'Transcription Factor Binding Site Analysis'  # type: str
    ah_reg = 'Regulation Analysis'                          # type: str
    # paths
    ai_gene = 'Genes file'                          # type: str
    ai_mut = 'Mutations file'                       # type: str
    ai_int = 'Interest file'                        # type: str
    ai_sm = 'Substitution matrix file'              # type: str
    ai_prot = 'Protein domains file'                # type: str
    ai_reg = 'Regulation file'                      # type: str
    ai_tfbs = 'TFBS file'                           # type: str
    ai_msa = 'TF MSA file'                          # type: str
    ao_dir = 'Output directory (analysis)'          # type: str
    # buttons
    ar_run = 'Run (analysis)'

    """ ANALYSIS - SETTINGS """
    as_title = 'Advanced Analysis Settings'             # type: str
    as_prom = 'Relative default promoter start'         # type: str
    as_pwm = 'min. number of sequences'                 # type: str
    as_rnd = 'Number of random mutations'               # type: str
    ar_adop = 'Adjust operon promoters'                 # type: str
    ar_open = 'Open output directory (analysis)'        # type: str
    ar_del = 'Delete old output (analysis)'             # type: str
    as_int = 'interesting name'                         # type: str
    as_int_short = 'interesting name (short)'           # type: str
    as_not_int = 'not-interesting name'                 # type: str
    as_not_int_short = 'not-interesting name (short)'   # type: str

    """ NGS - MAIN """
    # section and subsection headings
    nh_title = 'NGS Pipeline'           # type: str
    # paths
    ni_reads = 'Reads directory'        # type: str
    ni_ref = 'Reference genome file'    # type: str
    no_dir = 'Output directory (NGS)'   # type: str
    # buttons
    nr_run = 'Run (NGS)'

    """ NGS - SETTINGS """
    ns_title = 'Advanced NGS Pipeline Settings'     # type: str
    ne_var = 'VarScan executable'                   # type: str
    ne_bwa = 'BWA executable'                       # type: str
    ne_sam = 'SAMTools executable'                  # type: str
    ns_qual = 'SAMTools mapping quality'            # type: str
    ns_pval = 'VarScan SNP calling p-value'         # type: str
    nr_open = 'Open output directory (NGS)'         # type: str
    nr_clean = 'Clean output directory'             # type: str
    nr_del = 'Delete old output (NGS)'              # type: str

    """ CONVERTER """
    # VCF merger
    pm_title = 'Mutation VCF Merger'                    # type: str
    pm_odir = 'Output directory (vcf merge)'            # type: str
    pm_in = '.vcf input directory'                      # type: str
    pm_run = 'Run (vcf merger)'                         # type: str
    pm_open = 'Open output directory (vcf merger)'      # type: str
    # UniProt converter
    pu_title = 'UniProt Converter'                          # type: str
    pu_odir = 'Output directory (uniprot converter)'        # type: str
    pu_in = 'UniProt input file'                            # type: str
    pu_run = 'Run (uniprot converter)'                      # type: str
    pu_open = 'Open output directory (uniprot converter)'   # type: str
    # PATRIC converter
    pp_title = 'PATRIC Converter'                           # type: str
    pp_odir = 'Output directory (patric converter)'         # type: str
    pp_in = 'PATRIC input file'                             # type: str
    pp_run = 'Run (patric converter)'                       # type: str
    pp_open = 'Open output directory (patric converter)'    # type: str
    # RegulonDB converter
    pr_title = 'RegulonDB Converter'                            # type: str
    pr_odir = 'Output directory (regulondb converter)'          # type: str
    pr_msa = 'regulondb msa input file'                         # type: str
    pr_gene = 'regulondb gene input file'                       # type: str
    pr_gene_reg = 'regulondb gene regulation file'              # type: str
    pr_tf_reg = 'regulondb tf regulation file'                  # type: str
    pr_operon = 'regulondb operon file'                         # type: str
    pr_tfbs = 'regulondb tfbs_result file'                      # type: str
    pr_run = 'Run (regulondb converter)'                        # type: str
    pr_open = 'Open output directory (regulondb converter)'     # type: str
    pr_skip = 'Skip low quality'                                # type: str
    pr_prom = 'regulondb promoters'                             # type: str
    pr_tu = 'regulondb tus'                                     # type: str
    pr_op = 'regulondb operons'                                 # type: str
    pr_tu_reg = 'regulondb tu regulation'                       # type: str

    """ MENU BAR """
    # settings
    ms_title = 'Settings'
    ms_save = 'Save current settings'                   # type: str
    ms_def = 'Restore default settings'                 # type: str
    ms_save_def = 'Restore and save default settings'   # type: str
    ms_ana = 'Advanced analysis settings'               # type: str
    ms_ngs = 'Advanced NGS settings'                    # type: str
    # converter
    mc_title = 'File Converter'                 # type: str
    mc_vmerge = 'Merge mutation .vcf files'     # type: str
    mc_uconvert = 'Convert UniProt database'    # type: str
    mc_patric = 'Convert PATRIC'                # type: str
    mc_regulon = 'Convert RegulonDB'            # type: str
    # help
    mh_title = 'Help'                       # type: str
    mh_man = 'Open user manual'             # type: str
    mh_inst = 'Open installation guide'     # type: str
    mh_install = 'Install NGS'              # type: str

    """ STATUS """
    s_title = 'Status'      # type: str

    """ OTHERS """
    g_save = 'Save'         # type: str
    g_cancel = 'Cancel'     # type: str
    g_desc = 'Description'  # type: str


""" ANALYSIS INPUT FORMATS """


class GeneTSV:
    """
    Constants for the gene.tsv format.
    """
    lt = 'locus tag'                            # type: str
    name = 'gene name'                          # type: str
    dna = 'DNA'                                 # type: str
    genbank = 'Genbank'                         # type: str
    pfam = 'PFAM'                               # type: str
    uni = 'UniProt'                             # type: str
    desc = 'description'                        # type: str
    gene_id = 'Gene ID'                         # type: str
    end = 'gene end'                            # type: str
    start = 'gene start'                        # type: str
    operon = 'operon'                           # type: str
    prot = 'protein sequence'                   # type: str
    prom_end = 'promoter end absolute'          # type: str
    prom_astart = 'promoter start absolute'     # type: str
    prom_rstart = 'promoter start relative'     # type: str
    pubmed = 'Pubmed'                           # type: str
    strand = 'strand'                           # type: str
    chromosome = 'chromosome'                   # type: str


class MutationTSV:
    """
    Constants for the mutation.tsv format.
    """
    pos = 'position'        # type: str
    ref = 'reference'       # type: str
    alt = 'alternative'     # type: str
    chrom = 'chromosome'    # type: str


class RegulationTSV:
    """
    Constants for the regulation.tsv format.
    """
    name1 = 'regulator (gene name)'         # type: str
    lt1 = 'regulator (locus tag)'           # type: str
    name2 = 'regulated gene (gene name)'    # type: str
    lt2 = 'regulated gene (locus tag)'      # type: str
    desc = 'regulation'                     # type: str


class InterestTSV:
    """
    Constants for the resistance.tsv format.
    """
    name = 'gene name'          # type: str
    lt = 'locus tag'            # type: str
    desc = 'sub-category'       # type: str


class ProteinDomainTSV:
    """
    Constants for the protein domain.tsv format.
    """
    start = 'start'         # type: str
    end = 'end'             # type: str
    name = 'gene name'      # type: str
    lt = 'locus tag'        # type: str
    type = 'type'           # type: str
    desc = 'description'    # type: str


class TfbsTSV:
    """
    Constants for the TFBS.tsv format.
    """
    tf_name = 'TF gene name'        # type: str
    tf_lt = 'TF locus tag'          # type: str
    tfbs_name = 'TFBS gene name'    # type: str
    tfbs_lt = 'TFBS locus tag'      # type: str
    seq = 'sequence'                # type: str
    start = 'relative start'        # type: str
    astart = 'absolute start'       # type: str
    aend = 'absolute end'           # type: str


class SupportedInput:
    # support resistance information
    res_gen = ['ar', 'antibiotic resistance',
               'antibiotic-resistance', 'antibiotic resistant',
               'antibiotic-resistant', 'antibiotic',
               'antimicrobial', 'resistance', 'resistant', 'res']       # type: list
    res_reg = ['regulator', 'reg', 'efflux regulator',
               'mdr regulator', 'mdr efflux regulator',
               'mdre regulator']                                        # type: list
    res_pump = ['transporter', 'pump', 'mdr efflux transporter',
                'mdr efflux pump', 'mdre transporter', 'mdre pump',
                't', 'p', 'mdr pump', 'mdr transporter',
                'efflux transporter', 'efflux pump']                    # type: list
    # supported strand information
    str_plus = ['+', 'p', 'plus', 'forward']                            # type: list
    str_minus = ['-', 'm', 'minus', 'backward', 'reverse']              # type: list
    # supported regulation information
    reg_act = ['a', 'up-regulation', 'upregulation',
               'up regulation', 'up-regulator',
               'up regulator', 'upregulation', 'up',
               'u', Misc.activator.lower(), '+', 'plus',
               'activator', 'activation']                               # type: list
    reg_rep = ['r', 'down-regulation', 'downregulation',
               'down regulation', 'down-regulator',
               'down regulator', 'downregulator',
               'down', 'd', Misc.repressor.lower(), '-',
               'minus', 'repressor', 'repression']                      # type: list
    reg_eff = ['p', 'tf', 'as_prom', 'promoter',
               'promotor', 'transcription factor',
               'transcription', 'e', 'effector',
               'transcription-factor', 't-factor',
               Misc.effector.lower(), '+-', '-+',
               'plusminus', 'minusplus', 'plus-minus',
               'minus-plus']                                            # type: list


""" ANALYSIS OUTPUT FORMATS """


class MutationDB:
    """
    Constants for the mutation table in the database.
    """
    table = 'mutations'             # type: str
    id = 'id'                       # type: str
    ref = 'ref'                     # type: str
    alt = 'alt'                     # type: str
    pos = 'pos'                     # type: str
    chrom = 'chromosome'            # type: str
    type = 'type'                   # type: str
    int = 'of_interest'             # type: str
    lt = 'locus_tags'               # type: str
    reg = 'region'                  # type: str
    tfbs = 'TFBS_scoredif'          # type: str
    pd = 'protein_domains'          # type: str
    sub = 'substitution_score'      # type: str
    aa = 'amino_acid'               # type: str
    eff = 'effect_on_prot'          # type: str


class GeneDB:
    """
    Constants for the gene table in the database.
    """
    table = 'genes'                                     # type: str
    lt = 'locus_tag'                                    # type: str
    name = 'gene_name'                                  # type: str
    chrom = 'chromosome'                                # type: str
    int = 'of_interest'                                 # type: str
    desc = 'description'                                # type: str
    gstart = 'gene_start'                               # type: str
    gend = 'gene_end'                                   # type: str
    strand = 'strand'                                   # type: str
    glen = 'gene_len'                                   # type: str
    dna = 'DNA'                                         # type: str
    prot_len = 'prot_len'                               # type: str
    prot = 'prot'                                       # type: str
    op = 'operon'                                       # type: str
    reg1 = 'regulated_genes'                            # type: str
    reg2 = 'regulators'                                 # type: str
    pfam = 'PFAM'                                       # type: str
    genbank = 'Genbank'                                 # type: str
    uniprot = 'UniProt'                                 # type: str
    pubmed = 'pubmed'                                   # type: str
    id = 'gene_ID'                                      # type: str
    prom_astart = 'promoter_start_abs'                  # type: str
    prom_rstart = 'promoter_start_rel'                  # type: str
    prom_end = 'promoter_end_abs'                       # type: str
    prom_len = 'promoter_len'                           # type: str
    tf = 'TFs'                                          # type: str
    tfbs_astart = 'TFBS_start_abs'                      # type: str
    tfbs_rstart = 'TFBS_start_rel'                      # type: str
    tfbs_end = 'TFBS_end_abs'                           # type: str
    tfbs_seq = 'TFBS_sequence'                          # type: str
    mut_id = 'mutation_ids'                             # type: str
    md_rf = 'frame_shift_mutations_per_kbp'             # type: str
    md_rt = 'readthrough_mutations_per_kbp'             # type: str
    md_mis = 'missense_mutations_per_kbp'               # type: str
    md_non = 'nonsense_mutations_per_kbp'               # type: str
    md_ns = 'non_synonymous_mutations_per_kbp'          # type: str
    md_syn = 'synonymous_mutations_per_kbp'             # type: str
    md_prom = 'promoter_mutations_per_kbp'              # type: str
    md_tfbs = 'tfbs_mutations_per_kbp'                  # type: str
    pd_type = 'prot_dom_tyoe'                           # type: str
    pd_desc = 'prot_dom_desc'                           # type: str
    pd_start = 'prot_dom_start'                         # type: str
    pd_end = 'prot_dom_end'                             # type: str


class GeneMutationList:
    """
    Constants for the gene and mutation list.tsv format.
    """
    name = 'gene name'              # type: str
    lt = 'locus tag'                # type: str
    chrom = 'chromosome'            # type: str
    int = 'of interest'             # type: str
    reg1 = 'regulated genes'        # type: str
    reg2 = 'regulators'             # type: str
    dom = 'protein domain(s)'       # type: str
    id = 'ID'                       # type: str
    mut = 'mutation'                # type: str
    type = 'type'                   # type: str
    gpos = 'genome position'        # type: str
    ppos = 'protein position'       # type: str
    aa = 'amino acid'               # type: str
    eff = 'protein effect'          # type: str
    sub = 'substitution score'      # type: str
    tfbs = 'TFBS score'             # type: str
    loc = 'location(s)'             # type: str
    desc = 'description'            # type: str


""" CONFIGURATION COLLECTIONS """


class OutputPaths:
    """
    Constant output paths.
    """
    def __init__(self, a_out_dir, n_out_dir):
        a_out_dir = adjust_dir_path(a_out_dir)
        n_out_dir = adjust_dir_path(n_out_dir)
        # sub-directories
        plots = a_out_dir + 'individual_plots/'
        grn = a_out_dir + 'GRN_files/'
        tsv = a_out_dir + 'tsv_files/'
        # setup the sub-directories, the main ones are set up during the main routines
        if not os.path.exists(plots):
            os.makedirs(plots)
        if not os.path.exists(grn):
            os.makedirs(grn)
        if not os.path.exists(tsv):
            os.makedirs(tsv)
        # directory name for NGS .vcf files
        self.vcf = 'vcf_files/'                                             # type: str
        # statistics
        self.stats = plots + 'statistics.txt'                               # type: str
        # transcription factor binding site plot 
        self.tfbs = plots + 'TFBS_score_plot'                               # type: str
        # gene regulatory network
        self.grn_full = grn + 'GRN_complete.gml'                            # type: str
        self.grn_int = grn + 'GRN_of_interest.gml'                          # type: str
        # coding region plot
        self.coding = plots + 'coding_region_score_plot'                    # type: str
        # gene and mutation lists
        self.gm_int = tsv + 'genes_of_interest_table.tsv'                   # type: str
        self.mut = tsv + 'list_known_mutations.tsv'                         # type: str
        # mutation density and gene stats
        self.md = plots + 'mutation_density_plot'                           # type: str
        # mutation distribution and mutation stats
        self.mut_type = plots + 'mutation_type_distribution_plot'           # type: str
        self.mut_dist = plots + 'mutation_distribution_plot'                # type: str
        # database
        self.db = a_out_dir + 'genes_mutations_database.db'                 # type: str
        self.db_gene = tsv + 'database_gene_table.tsv'                      # type: str
        self.db_mut = tsv + 'database_mutation_table.tsv'                   # type: str
        # result pdf
        self.result = a_out_dir + 'results.pdf'
        # analysis log
        self.a_log = a_out_dir + 'log.txt'                                  # type: str
        # NGS log
        self.n_log = n_out_dir + 'log.txt'                                  # type: str
        # NGS mutation file
        self.n_mut = n_out_dir + 'mutations.tsv'                            # type: str


class RunConfig:
    """
    Flags for running different parts of the program.
    """
    def __init__(self, cfg_dict):
        self.cfg_dict = cfg_dict.copy()                                         # type: dict

        # analysis
        self.a_coding = bool(cfg_dict['run coding region analysis'])            # type: bool
        self.a_tfbs = bool(cfg_dict['run TFBS analysis'])                       # type: bool
        self.a_int = bool(cfg_dict['run genes of interest analysis'])           # type: bool
        self.a_reg = bool(cfg_dict['run regulation analysis'])                  # type: bool
        self.a_prot = bool(cfg_dict['run protein domain analysis'])             # type: bool
        self.a_open = bool(cfg_dict['open analysis output directory'])          # type: bool
        self.a_adop = bool(cfg_dict['run adjust operon promoters'])             # type: bool
        self.a_del = bool(cfg_dict['delete old analysis output'])               # type: bool

        # NGS
        self.n_open = bool(cfg_dict['open NGS output directory'])               # type: bool
        self.n_clean = bool(cfg_dict['clean up NGS output directory'])          # type: bool

        # parser
        self.p_open = bool(cfg_dict['open parser output directory'])            # type: bool

    def update_cfg_dict(self):
        """
        Updates and returns the configuration dictionary.
        :return: updated configuration dictionary
        """
        # analysis
        self.cfg_dict['run coding region analysis'] = bool(self.a_coding)
        self.cfg_dict['run TFBS analysis'] = bool(self.a_tfbs)
        self.cfg_dict['run genes of interest analysis'] = bool(self.a_int)
        self.cfg_dict['run regulation analysis'] = bool(self.a_reg)
        self.cfg_dict['run protein domain analysis'] = bool(self.a_prot)
        self.cfg_dict['open analysis output directory'] = bool(self.a_open)
        self.cfg_dict['run adjust operon promoters'] = bool(self.a_adop)

        # NGS
        self.cfg_dict['open NGS output directory'] = bool(self.n_open)
        self.cfg_dict['clean up NGS output directory'] = bool(self.n_clean)

        # parser
        self.cfg_dict['open parser output directory'] = bool(self.p_open)

        return self.cfg_dict


class UserConfig:
    """
    Analysis, NGS and converter input paths, settings and output directories.
    """
    def __init__(self, cfg_dict):
        self.cfg_dict = cfg_dict.copy()                                                     # type: dict
        """ ANALYSIS """
        # input paths
        self.ai_gene = adjust_file_path(cfg_dict['Analysis']['gene path'])                  # type: str
        self.ai_mut = adjust_file_path(cfg_dict['Analysis']['mutation path'])               # type: str
        self.ai_tfbs = adjust_file_path(cfg_dict['Analysis']['TFBS path'])                  # type: str
        self.ai_prot = adjust_file_path(cfg_dict['Analysis']['protein domain path'])        # type: str
        self.ai_reg = adjust_file_path(cfg_dict['Analysis']['regulation path'])             # type: str
        self.ai_sm = adjust_file_path(cfg_dict['Analysis']['substitution matrix path'])     # type: str
        self.ai_int = adjust_file_path(cfg_dict['Analysis']['interest path'])               # type: str
        self.ai_msa = adjust_file_path(cfg_dict['Analysis']['TF MSA path'])                 # type: str
        # output directory
        self.ao_dir = adjust_dir_path(cfg_dict['Analysis']['output directory'])             # type: str
        # settings
        self.as_prom = int(cfg_dict['Analysis']['relative promoter start'])                 # type: int
        self.as_rnd = int(cfg_dict['Analysis']['TFBS random mutations'])                    # type: int
        self.as_pwm = int(cfg_dict['Analysis']['PWM number of sequence threshold'])         # type: int
        self.as_int = cfg_dict['Analysis']['genes of interest (name long)']                 # type: str
        self.as_int_short = cfg_dict['Analysis']['genes of interest (name short)']          # type: str
        self.as_not_int = cfg_dict['Analysis']['other genes (name long)']                   # type: str
        self.as_not_int_short = cfg_dict['Analysis']['other genes (name short)']            # type: str

        """ NGS """
        # input paths
        self.ni_reads = adjust_dir_path(cfg_dict['NGS']['reads directory'])                 # type: str
        self.ni_ref = adjust_file_path(cfg_dict['NGS']['reference genome path'])            # type: str
        # executables
        self.ne_samtools = find_exe('samtools', cfg_dict['NGS']['SAMTools executable'])     # type: str
        self.ne_varscan = find_exe('varscan', cfg_dict['NGS']['VarScan executable'])        # type: str
        self.ne_bwa = find_exe('bwa', cfg_dict['NGS']['BWA executable'])                    # type: str
        # output directory
        self.no_dir = adjust_dir_path(cfg_dict['NGS']['output directory'])                  # type: str
        # settings
        self.ns_qcut = int(cfg_dict['NGS']['quality control cut-off'])                      # type: int
        self.ns_pv = float(cfg_dict['NGS']['p-value'])                                      # type: float

        """ PARSER """
        # UniProt
        self.pu_out = adjust_file_path(cfg_dict['UniProt']['output file'])                  # type: str
        self.pu_in = adjust_file_path(cfg_dict['UniProt']['input file'])                    # type: str
        # PATRIC
        self.pp_out = adjust_file_path(cfg_dict['PATRIC']['output file'])                   # type: str
        self.pp_in = adjust_file_path(cfg_dict['PATRIC']['input file'])                     # type: str
        # VCF Merger
        self.pm_in = adjust_dir_path(cfg_dict['VCF Merger']['input directory'])             # type: str
        self.pm_out = adjust_file_path(cfg_dict['VCF Merger']['output file'])               # type: str
        # RegulonDB
        self.pr_out = adjust_dir_path(cfg_dict['RegulonDB']['output directory'])            # type: str
        self.pr_gene = adjust_file_path(cfg_dict['RegulonDB']['gene input'])                # type: str
        self.pr_operon = adjust_file_path(cfg_dict['RegulonDB']['operon input'])            # type: str
        self.pr_prom = adjust_file_path(cfg_dict['RegulonDB']['promoter input'])            # type: str
        self.pr_tfbs = adjust_file_path(cfg_dict['RegulonDB']['TFBS input'])                # type: str
        self.pr_tu = adjust_file_path(cfg_dict['RegulonDB']['TU input'])                    # type: str
        self.pr_msa = adjust_file_path(cfg_dict['RegulonDB']['TF MSA input'])               # type: str
        self.pr_tg = adjust_file_path(cfg_dict['RegulonDB']['TF-gene input'])               # type: str
        self.pr_tftf = adjust_file_path(cfg_dict['RegulonDB']['TF-TF input'])               # type: str
        self.pr_tftu = adjust_file_path(cfg_dict['RegulonDB']['TF-TU input'])               # type: str
        self.pr_to = adjust_file_path(cfg_dict['RegulonDB']['TF-operon input'])             # type: str

    def update_cfg_dict(self):
        """
        Updates and returns the configuration dictionary.
        :return: updated configuration dictionary
        """
        """ ANALYSIS """
        # input paths
        self.cfg_dict['Analysis']['gene path'] = adjust_file_path(self.ai_gene)
        self.cfg_dict['Analysis']['mutation path'] = adjust_file_path(self.ai_mut)
        self.cfg_dict['Analysis']['TFBS path'] = adjust_file_path(self.ai_tfbs)
        self.cfg_dict['Analysis']['protein domain path'] = adjust_file_path(self.ai_prot)
        self.cfg_dict['Analysis']['regulation path'] = adjust_file_path(self.ai_reg)
        self.cfg_dict['Analysis']['substitution matrix path'] = adjust_file_path(self.ai_sm)
        self.cfg_dict['Analysis']['interest path'] = adjust_file_path(self.ai_int)
        self.cfg_dict['Analysis']['TF MSA path'] = adjust_file_path(self.ai_msa)
        # output directory
        self.cfg_dict['Analysis']['output directory'] = adjust_dir_path(self.ao_dir)
        # settings
        self.cfg_dict['Analysis']['relative promoter start'] = int(self.as_prom)
        self.cfg_dict['Analysis']['TFBS random mutations'] = int(self.as_rnd)
        self.cfg_dict['Analysis']['PWM number of sequence threshold'] = int(self.as_pwm)
        self.cfg_dict['Analysis']['genes of interest (name long)'] = self.as_int
        self.cfg_dict['Analysis']['genes of interest (name short)'] = self.as_int_short
        self.cfg_dict['Analysis']['other genes (name long)'] = self.as_not_int
        self.cfg_dict['Analysis']['other genes (name short)'] = self.as_not_int_short

        """ NGS """
        # input paths
        self.cfg_dict['NGS']['reads directory'] = adjust_dir_path(self.ni_reads)
        self.cfg_dict['NGS']['reference genome path'] = adjust_file_path(self.ni_ref)
        # executables
        self.cfg_dict['NGS']['SAMTools executable'] = adjust_file_path(self.ne_samtools)
        self.cfg_dict['NGS']['VarScan executable'] = adjust_file_path(self.ne_varscan)
        self.cfg_dict['NGS']['BWA executable'] = adjust_file_path(self.ne_bwa)
        # output directory
        self.cfg_dict['NGS']['output directory'] = adjust_dir_path(self.no_dir)
        # settings
        self.cfg_dict['NGS']['quality control cut-off'] = int(self.ns_qcut)
        self.cfg_dict['NGS']['p-value'] = float(self.ns_pv)

        return self.cfg_dict.copy()

    def update_regulon(self):
        """
        Updates and returns the RegulonDB configuration dictionary.
        :return: updated configuration dictionary
        """
        self.cfg_dict['RegulonDB']['output directory'] = adjust_dir_path(self.pr_out)
        self.cfg_dict['RegulonDB']['gene input'] = adjust_file_path(self.pr_gene)
        self.cfg_dict['RegulonDB']['operon input'] = adjust_file_path(self.pr_operon)
        self.cfg_dict['RegulonDB']['promoter input'] = adjust_file_path(self.pr_prom)
        self.cfg_dict['RegulonDB']['TFBS input'] = adjust_file_path(self.pr_tfbs)
        self.cfg_dict['RegulonDB']['TU input'] = adjust_file_path(self.pr_tu)
        self.cfg_dict['RegulonDB']['TF MSA input'] = adjust_file_path(self.pr_msa)
        self.cfg_dict['RegulonDB']['TF-gene input'] = adjust_file_path(self.pr_tg)
        self.cfg_dict['RegulonDB']['TF-TF input'] = adjust_file_path(self.pr_tftf)
        self.cfg_dict['RegulonDB']['TF-TU input'] = adjust_file_path(self.pr_tftu)
        self.cfg_dict['RegulonDB']['TF-operon input'] = adjust_file_path(self.pr_to)

        return self.cfg_dict['RegulonDB'].copy()

    def update_uniprot(self):
        """
        Updates and returns the UniProt configuration dictionary.
        :return: updated configuration dictionary
        """
        self.cfg_dict['UniProt']['output file'] = adjust_file_path(self.pu_out)
        self.cfg_dict['UniProt']['input file'] = adjust_file_path(self.pu_in)

        return self.cfg_dict['UniProt'].copy()

    def update_patric(self):
        """
        Updates and returns the PATRIC configuration dictionary.
        :return: updated configuration dictionary
        """
        self.cfg_dict['PATRIC']['output file'] = adjust_file_path(self.pp_out)
        self.cfg_dict['PATRIC']['input file'] = adjust_file_path(self.pp_in)

        return self.cfg_dict['PATRIC'].copy()

    def update_merger(self):
        """
        Updates and returns the VCF merger configuration dictionary.
        :return: updated configuration dictionary
        """
        self.cfg_dict['VCF Merger']['input directory'] = adjust_dir_path(self.pm_in)
        self.cfg_dict['VCF Merger']['output file'] = adjust_file_path(self.pm_out)

        return self.cfg_dict['VCF Merger'].copy()


class GuiConfig:
    def __init__(self, cfg_dict):
        self.cfg_dict = cfg_dict    # type: dict

        """ COLOURS """
        # foreground
        self.st_colour = cfg_dict['colours']['section title']                   # type: str
        self.sst_colour = cfg_dict['colours']['sub-section title']              # type: str
        self.dt_colour = cfg_dict['colours']['description title']               # type: str
        self.dd_colour = cfg_dict['colours']['description text']                # type: str

        # background

        """ FONTS """
        self.st_font = tuple(cfg_dict['fonts']['section title'])
        self.sst_font = tuple(cfg_dict['fonts']['sub-section title'])
        self.sl_font = tuple(cfg_dict['fonts']['status line'])
        self.el_font = tuple(cfg_dict['fonts']['entry field label'])
        self.ef_font = tuple(cfg_dict['fonts']['entry field text'])
        self.dt_font = tuple(cfg_dict['fonts']['description title'])
        self.dd_font = tuple(cfg_dict['fonts']['description text'])
        self.cb_font = tuple(cfg_dict['fonts']['check button'])
        self.fb_font = tuple(cfg_dict['fonts']['file chooser button'])
        self.rb_font = tuple(cfg_dict['fonts']['run button'])

        """ STYLES """
        self.theme = cfg_dict['styles']['theme']                                # type: str
        self.st_style = 'SectionTitle.TLabel'                                   # type: str
        self.sst_style = 'SubSectionTitle.TLabel'                               # type: str
        self.sl_style = 'StatusLine.TLabel'                                     # type: str
        self.el_style = 'EntryLabel.TLabel'                                     # type: str
        self.dt_style = 'DescriptionTitle.TLabel'                               # type: str
        self.cb_style = 'Custom.TCheckbutton'                                   # type: str
        self.fb_style = 'FileChooser.TButton'                                   # type: str
        self.rb_style = 'RunButton.TButton'

        """ WIDTHS """
        # entry fields
        self.fe_width = cfg_dict['widths']['file/directory entry field']        # type: int
        self.ne_width = cfg_dict['widths']['number entry field']                # type: int
        # buttons
        self.fb_width = cfg_dict['widths']['file/directory chooser button']     # type: int
        self.rb_width = cfg_dict['widths']['run button']                        # type: int
        self.rsb_width = cfg_dict['widths']['right settings button']            # type: int
        self.lsb_width = cfg_dict['widths']['left settings button']             # type: int
        # borders
        self.b_width = cfg_dict['widths']['status box']                         # type: int
        self.sb_width = cfg_dict['widths']['status box border']                 # type: int
        self.sp_width = cfg_dict['widths']['separator line width']              # type: int
        # labels/messages
        self.dd_width = cfg_dict['widths']['description text']                  # type: int
        # frames
        self.df_width = cfg_dict['widths']['description frame']                 # type: int

        """ HEIGHT """
        # frames
        self.df_height = cfg_dict['heights']['description frame']               # type: int

        """ X-MARGINS """
        # labels
        self.st_padx = tuple(cfg_dict['x-margins']['section title'])            # type: (int, int)
        self.sst_padx = tuple(cfg_dict['x-margins']['subsection title'])        # type: (int, int)
        # window borders
        self.wl_padx = cfg_dict['x-margins']['window left']                     # type: (int, int)
        self.wr_padx = cfg_dict['x-margins']['window right']                    # type: (int, int)
        self.col_padx = cfg_dict['x-margins']['column']                         # type: (int, int)
        # separators
        self.sp_padx = tuple(cfg_dict['x-margins']['separator'])                # type: (int, int)
        self.tss_padx = tuple(cfg_dict['x-margins']['text separator'])          # type: (int, int)
        # buttons
        self.rb_padx = tuple(cfg_dict['x-margins']['run button'])               # type: (int, int)
        self.lsb_padx = tuple(cfg_dict['x-margins']['left settings button'])    # type: (int, int)
        self.rsb_padx = tuple(cfg_dict['x-margins']['right settings button'])   # type: (int, int)
        self.fb_padx = tuple(cfg_dict['x-margins']['file chooser button'])      # type: (int, int)
        # check buttons
        self.tcb_padx = tuple(cfg_dict['x-margins']['text check button'])       # type: (int, int)
        self.ecb_padx = tuple(cfg_dict['x-margins']['empty check button'])      # type: (int, int)
        # mouse-over description
        self.dt_padx = tuple(cfg_dict['x-margins']['description title'])        # type: (int, int)
        self.dd_padx = tuple(cfg_dict['x-margins']['description text'])         # type: (int, int)
        # entry fields
        self.el_padx = tuple(cfg_dict['x-margins']['entry field label'])        # type: (int, int)
        self.ef_padx = tuple(cfg_dict['x-margins']['entry field text'])         # type: (int, int)
        # status box
        self.sf_padx = tuple(cfg_dict['x-margins']['status box frame'])         # type: (int, int)

        """ Y-MARGINS """
        # labels
        self.st_pady = tuple(cfg_dict['y-margins']['section title'])            # type: (int, int)
        self.sst_pady = tuple(cfg_dict['y-margins']['subsection title'])        # type: (int, int)
        # window borders
        self.wb_pady = cfg_dict['y-margins']['window bottom']                   # type: (int, int)
        # separators
        self.sp_pady = tuple(cfg_dict['y-margins']['separator'])                # type: (int, int)
        self.tss_pady = tuple(cfg_dict['y-margins']['text separator'])          # type: (int, int)
        # buttons
        self.rb_pady = tuple(cfg_dict['y-margins']['run button'])               # type: (int, int)
        self.lsb_pady = tuple(cfg_dict['y-margins']['left settings button'])    # type: (int, int)
        self.rsb_pady = tuple(cfg_dict['y-margins']['right settings button'])   # type: (int, int)
        self.fb_pady = tuple(cfg_dict['y-margins']['file chooser button'])      # type: (int, int)
        # check buttons
        self.tcb_pady = tuple(cfg_dict['y-margins']['text check button'])       # type: (int, int)
        self.ecb_pady = tuple(cfg_dict['y-margins']['empty check button'])      # type: (int, int)
        # mouse-over description
        self.dt_pady = tuple(cfg_dict['y-margins']['description title'])        # type: (int, int)
        self.dd_pady = tuple(cfg_dict['y-margins']['description text'])         # type: (int, int)
        # entry fields
        self.el_pady = tuple(cfg_dict['y-margins']['entry field label'])        # type: (int, int)
        self.ef_pady = tuple(cfg_dict['y-margins']['entry field text'])         # type: (int, int)
        # status box
        self.sf_pady = tuple(cfg_dict['y-margins']['status box frame'])         # type: (int, int)

        """ CONSTANTS """
        self.key = GuiKeys                  # type: GuiKeys
        self.desc = self.descriptions()     # type: dict
        # application title
        self.title = 'MutaNET'              # type: str

    def descriptions(self):
        """
        Sets up the labels and descriptions of different GUI elements.
        :return: dictionary with labels and descriptions
        """
        d = dict()

        """ MENU BAR """
        # settings
        d[self.key.ms_title] = GuiEntity(label='Settings',
                                         desc='TODO')
        d[self.key.ms_save] = GuiEntity(label='Save current settings',
                                        desc='TODO')
        d[self.key.ms_def] = GuiEntity(label='Restore default settings',
                                       desc='TODO')
        d[self.key.ms_save_def] = GuiEntity(label='Restore and save default settings',
                                            desc='TODO')
        d[self.key.ms_ana] = GuiEntity(label='Advanced analysis settings',
                                       desc='TODO')
        d[self.key.ms_ngs] = GuiEntity(label='Advanced NGS settings',
                                       desc='TODO')
        # converter
        d[self.key.mc_title] = GuiEntity(label='Tools',
                                         desc='TODO')
        d[self.key.mc_vmerge] = GuiEntity(label='Mutation VCF merger',
                                          desc='TODO')
        d[self.key.mc_uconvert] = GuiEntity(label='UniProt protein domain converter',
                                            desc='TODO')
        d[self.key.mc_patric] = GuiEntity(label='PATRIC antibiotic resistance converter',
                                          desc='TODO')
        d[self.key.mc_regulon] = GuiEntity(label='RegulonDB converter',
                                           desc='')
        # help
        d[self.key.mh_title] = GuiEntity(label='Help',
                                         desc='TODO')
        d[self.key.mh_man] = GuiEntity(label='Open user manual',
                                       desc='TODO')
        d[self.key.mh_inst] = GuiEntity(label='Open installation guide',
                                        desc='TODO')
        d[self.key.mh_install] = GuiEntity(label='Install BWA, SAMTools and VarScan',
                                           desc='TODO')

        """ ANALYSIS - MAIN """
        # section and subsection titles
        d[self.key.ah_title] = GuiEntity(label='Mutation Analysis',
                                         desc='Takes a set of genes and mutations, maps the mutations onto different '
                                              'gene regions, computes various statistics and performs optional '
                                              'analysis steps.')
        d[self.key.ah_int] = GuiEntity(label='Genes of Interest Analysis',
                                       desc='Analyse the impact of mutations in genes of interest. The '
                                            'analysis becomes more in-depth the more other analysis steps are '
                                            'enabled and the more information is made available.')
        d[self.key.ah_coding] = GuiEntity(label='Coding Region Analysis',
                                          desc='Analyse the impact of mutations in coding regions by using an '
                                               'amino acid substitution matrix. If antibiotic resistance analysis is '
                                               'enable, additional statistics related to antibiotic resistance are '
                                               'computed.')
        d[self.key.ah_prot] = GuiEntity(label='Protein Domain Analysis',
                                        desc='Map mutations to protein domains.')
        d[self.key.ah_tfbs] = GuiEntity(label='Transcription Factor Binding Site Analysis',
                                        desc='Analyse the impact of mutations in transcription factor binding sites '
                                             'by using position weight matrices to compare them to randomly '
                                             'generated mutations. If antibiotic resistance analysis is enabled, '
                                             'additional statistics related to antibiotic resistance are computed.')
        d[self.key.ah_reg] = GuiEntity(label='Regulation Analysis',
                                       desc='Build a gene regulatory network of the organism. If antibiotic '
                                            'resistance analysis is enabled, the sub-network of antibiotic resistance '
                                            'genes and their regulators is computed. Additionally, a table with '
                                            'all non-synonymous mutations and their scores is created for the '
                                            'sub-network.')
        # paths
        d[self.key.ai_gene] = GuiEntity(label='Genes file',
                                        desc='Required file with columns for gene information. See the user manual for '
                                             'more information on the format, such as required and optional columns.'
                                             '\n\nFormat: tab-separated (.tsv)')
        d[self.key.ai_mut] = GuiEntity(label='Mutations file',
                                       desc='Required file columns for mutation information. See the user manual for '
                                             'more information on the format, such as the required columns.'
                                             '\n\nFormat: tab-separated (.tsv)')
        d[self.key.ai_int] = GuiEntity(label='Genes of interest file',
                                       desc='File with locus tag and/or gene name of genes of interest. It is possible '
                                            'to specify up to 15 sub-categories in this file. See the user manual for '
                                            'more information on the format.'
                                            '\n\nFormat: tab-separated (.tsv)')
        d[self.key.ai_sm] = GuiEntity(label='Substitution matrix file',
                                      desc='File with a column and row for gaps and for each amino acid with scores '
                                           'of amino acid substitutions. See the user manual for more information '
                                           'on the format.'
                                           '\n\nThe example input (/example_data/input/analysis/) contains the PAM10, '
                                           'PAM30 and PAM50 matrices.'
                                           '\n\nFormat: tab-separated (.tsv)')
        d[self.key.ai_prot] = GuiEntity(label='Protein domains file',
                                        desc='File with columns for protein domain information. See the user manual '
                                             'for information on the format, such as the required columns.\n\n'
                                             'Format: tab-separated (.tsv)')
        d[self.key.ai_reg] = GuiEntity(label='Regulation file',
                                       desc='File with information on regulators, regulated genes and the type of '
                                            'regulation. See the user manual for more information on the format '
                                            'and supported regulation descriptions.\n\n'
                                            'Format: tab-separated (.tsv)')
        d[self.key.ai_tfbs] = GuiEntity(label='TFBS file',
                                        desc='File with information on transcription factor binding sites (TFBSs). '
                                             'See the user manual for more information on the format.\n\n'
                                             'Format: tab-separated (.tsv)')
        d[self.key.ai_msa] = GuiEntity(label='TF motif MSA file',
                                       desc='File with multiple sequence alignments (MSA) of the binding motifs of '
                                            'transcription factors (TFs). See the user manual for more information '
                                            'on the format.\n\n'
                                            'Format: .fasta')
        d[self.key.ao_dir] = GuiEntity(label='Result directory',
                                       desc='Directory that is going to contain a SQLite database and .tsv file with '
                                            'the mapped gene and mutation information, as well as scores and '
                                            'statistics if the corresponding analysis steps are enabled.\n\n'
                                            'Additionally, it is going to contain plot, gene regulatory network '
                                            'files and summary tables of the different analysis steps, as well as '
                                            'a log file with information on run time errors that were not grave '
                                            'enough to stop the analysis.')
        # buttons
        d[self.key.ar_run] = GuiEntity(label='Run',
                                       desc='Takes a set of genes and mutations, maps the mutations onto different '
                                            'gene regions, computes various statistics and performs the enabled '
                                            'analysis steps. The window might stop responding while the '
                                            'analysis is running in the background; this has no impact on the '
                                            'analysis.\n\n'
                                            'Warning: If the result directory already exists, this will delete '
                                            'all its content!')

        """ ANALYSIS - SETTINGS """
        d[self.key.as_title] = GuiEntity(label='Advanced Analysis Settings',
                                         desc='Options to further customise various analysis steps.')
        d[self.key.as_prom] = GuiEntity(label='Relative default promoter start',
                                        desc='If no promoter region information is given for a gene in the genes .tsv '
                                             'file, this value is used to set the promoter region of that gene from '
                                             'the given value to the transcription start.\n\n'
                                             'Example: Let the genome position of the transcription start be 1800 and '
                                             'the this value be -100, then, the promoter region is going to be from '
                                             '17900 to 1800.')
        d[self.key.as_rnd] = GuiEntity(label='Number of random mutations',
                                       desc='The number of mutations that are randomly generated for the analysis of '
                                            'each transcription factor binding site to investigate the impact of the '
                                            'observed mutation in that binding site. This value should not be higher '
                                            'than the length of these binding sites.')
        d[self.key.ar_adop] = GuiEntity(label='Adjust operon promoters and TFBSs',
                                        desc='If enabled, sets the promoter regions of all genes belong to an operon '
                                             'to the promoter region of the first gene in the operon, and merges the '
                                             'TFBS information of all genes in the operon.')
        d[self.key.ar_open] = GuiEntity(label='Open result directory',
                                        desc='If enabled, opens the result directory in the operating '
                                             'system\'s file explorer after completing all analysis steps.')
        d[self.key.ar_del] = GuiEntity(label='Delete old output',
                                       desc='TODO')
        d[self.key.as_pwm] = GuiEntity(label='Minimum number of sequences',
                                       desc='The minimum number of valid DNA sequences with equal length that must be '
                                            'given for a single transcription factor. If the number of sequences is '
                                            'below this number, the transcription factor is not considered in the '
                                            'transcription factor binding site analysis.')
        d[self.key.as_int] = GuiEntity(label='Genes of interest (name long)',
                                       desc='The name that is going to be used in the result plots and tables for '
                                            'the genes of interest. For example, when studying antibiotic resistance, '
                                            'this might be "antibiotic resistance".')
        d[self.key.as_int_short] = GuiEntity(label='Genes of interest (name short)',
                                             desc='The shortened name that is going to be used in some result plots '
                                                  'for the genes of interest. For example, when studying antibiotic '
                                                  'resistance, this might be "AR".')
        d[self.key.as_not_int] = GuiEntity(label='Other genes (name long)',
                                           desc='The name that is going to be used in the result plots and tables for '
                                                'the genes that are not of special interest. For example, when '
                                                'studying antibiotic resistance, this might be "non-antibiotic '
                                                'resistance".')
        d[self.key.as_not_int_short] = GuiEntity(label='Other genes (name short)',
                                                 desc='The shortened name is going to be used in some result plots '
                                                      'for genes that are not of special interest. For example, when '
                                                      'studying antibiotic resistance, this might be "non-AR".')

        """ NGS - MAIN """
        # section and subsection headings
        d[self.key.nh_title] = GuiEntity(label='NGS Pipeline',
                                         desc='Takes a reference genome and NGS reads to align the reads to the '
                                              'genome, perform quality control and variant calling using BWA, '
                                              'SAMTools and VarScan.\n\n'
                                              'Installation instructions for BWA, SAMTools and VarScan on Linux, '
                                              'Mac OS X and Windows can be found under '
                                              '\'{0}\' --> \'{1}\'.'.format(d[self.key.mh_title].label,
                                                                            d[self.key.mh_inst].label))
        # paths
        d[self.key.ni_reads] = GuiEntity(label='Reads directory',
                                         desc='Directory with (paired) NGS read files. See the user manual for '
                                              'further information on the format.\n\n'
                                              'Format: name_1.fastq, name_2.fastq')
        d[self.key.ni_ref] = GuiEntity(label='Reference genome file',
                                       desc='File with the reference genome the NGS reads are aligned to.\n\n'
                                            'Format: .fasta')
        d[self.key.no_dir] = GuiEntity(label='Result directory',
                                       desc='Directory that is going to contain variant calling format (.vcf) files '
                                            'with information on mutations, as well as a log with errors that might '
                                            'have occurred during run time without being grave enough to stop the '
                                            'pipeline. If specified in the advanced settings, it '
                                            'is also going to contain files of the intermediate steps.')
        # buttons
        d[self.key.nr_run] = GuiEntity(label='Run',
                                       desc='Takes a reference genome and NGS reads to align the reads to the '
                                            'genome, perform quality control and variant calling using BWA, '
                                            'SAMTools and VarScan. The window might stop responding while the '
                                            'pipeline is running in the background,; this has no impact on the '
                                            'pipeline.\n\n'
                                            'Warning: If the result directory already exists, this will delete '
                                            'all its content!')

        """ NGS - SETTINGS """
        d[self.key.ns_title] = GuiEntity(label='Advanced NGS Pipeline Settings',
                                         desc='Options to further customise the NGS pipeline.')
        d[self.key.ne_var] = GuiEntity(label='VarScan executable',
                                       desc='Path to the VarScan executable, which is required for variant '
                                            'calling.\n\n'
                                            'Installation instructions for VarScan on Linux, '
                                            'Mac OS X and Windows can be found under '
                                            '\'{0}\' --> \'{1}\'.'.format(d[self.key.mh_title].label,
                                                                          d[self.key.mh_inst].label))
        d[self.key.ne_bwa] = GuiEntity(label='BWA executable',
                                       desc='Path to the Burrows-Wheeler-Aligner executable, which is required for '
                                            'aligning the reads to the reference genome.\n\n'
                                            'Installation instructions for BWA on Linux, '
                                            'Mac OS X and Windows can be found under '
                                            '\'{0}\' --> \'{1}\'.'.format(d[self.key.mh_title].label,
                                                                          d[self.key.mh_inst].label))
        d[self.key.ne_sam] = GuiEntity(label='SAMTools executable',
                                       desc='Path to the SAMTools executable, which is required for indexing and '
                                            'quality control.\n\n'
                                            'Installation instructions for SAMTools on Linux, '
                                            'Mac OS X and Windows can be found under '
                                            '\'{0}\' --> \'{1}\'.'.format(d[self.key.mh_title].label,
                                                                          d[self.key.mh_inst].label))
        d[self.key.ns_qual] = GuiEntity(label='SAMTools mapping quality',
                                        desc='Reads with mapping quality below this threshold are ignored.')
        d[self.key.ns_pval] = GuiEntity(label='VarScan SNP calling p-value',
                                        desc='P-value threshold for calling variants.')
        d[self.key.nr_open] = GuiEntity(label='Open result directory',
                                        desc='If enabled, opens the result directory in the operating '
                                             'system\'s file explorer after completing all analysis steps.')
        d[self.key.nr_clean] = GuiEntity(label='Clean intermediate results',
                                         desc='If enabled, deletes the intermediate results of the pipeline to save '
                                              'disc space.')
        d[self.key.nr_del] = GuiEntity(label='Delete old output',
                                       desc='TODO')

        """ CONVERTER """
        # VCF merger
        d[self.key.pm_title] = GuiEntity(label='Mutation VCF Merger',
                                         desc='Merges variant call format (.vcf) files with mutation information into '
                                              'a single tab-separated (.tsv) file. That file can be used in the '
                                              'mutation analysis of this tool and contains information on the '
                                              'genome position, reference and alternative DNA base(s) of the '
                                              'mutations.')
        d[self.key.pm_odir] = GuiEntity(label='Result file',
                                        desc='Path to the result .tsv file. It can be used in the mutaiton analysis '
                                             'and contains information on genome position, reference and alternative '
                                             'DNA base(s) of mutations.\n\n'
                                             'If the file does not exist, it is '
                                             'automatically created if the containing directory exists. It can also '
                                             'be created via the file dialog.\n\n'
                                             'File type: tab-separated (.tsv)')
        d[self.key.pm_in] = GuiEntity(label='VCF input directory',
                                      desc='Directory that contains .vcf files with mutation information, including '
                                           'in its sub-directories.\n\n'
                                           'These .vcf files need to contain non-empty, tab-separated columns called '
                                           '\'POS\', \'REF\' and \'AlT\'. The order of these columns does not matter '
                                           'and additional columns are not a problem.\n\n'
                                           'File type: variant call format (.vcf)')
        d[self.key.pm_run] = GuiEntity(label='Run',
                                       desc='Merges the .vcf files in the specified input directory and its sub-'
                                            'directories into a single .tsv file that can be used in the '
                                            'mutation analysis.\n\n'
                                            'Warning: If the result file already exists, this will override its '
                                            'content!')
        d[self.key.pm_open] = GuiEntity(label='Open result directory',
                                        desc='If enabled, opens the directory containing the result file after '
                                             'successfully merging the .vcf files.')

        # UniProt converter
        d[self.key.pu_title] = GuiEntity(label='UniProt Protein Domain Converter',
                                         desc='Extracts protein domains from a UniProt database text (.txt) file '
                                              'and returns a tab-separated (.tsv) file with these protein domains. '
                                              'That file can be used in the mutation analysis and '
                                              'contains information on locus tag, gene name, start, end, type and '
                                              'description of protein domains.')
        d[self.key.pu_odir] = GuiEntity(label='Result file',
                                        desc='Path to the result .tsv file. It can be used for protein domain '
                                             'analysis in the mutation analysis and contains '
                                             'gene name and domain information.\n\n'
                                             'If the file does not exist, it is '
                                             'automatically created if the containing directory exists. It can also '
                                             'be created via the file dialog.\n\n'
                                             'File type: tab-separated (.tsv)')
        d[self.key.pu_in] = GuiEntity(label='UniProt input file',
                                      desc='A text file that contains protein entries from the UniProt database.\n\n'
                                           'File type: text (.txt)')
        d[self.key.pu_run] = GuiEntity(label='Run',
                                       desc='Converts the UniProt file to a protein domain .tsv file that can be '
                                            'used in the mutation analysis.\n\n'
                                            'Warning: If the result file already exists, this will override its '
                                            'content!')
        d[self.key.pu_open] = GuiEntity(label='Open result directory',
                                        desc='If enabled, opens the directory containing the result file after '
                                             'successfully converting the UniProt file.')
        # PATRIC converter
        d[self.key.pp_title] = GuiEntity(label='PATRIC Antibiotic Resistance Converter',
                                         desc='Converts a comma-separated (.csv) antibiotic resistance file from '
                                              'PATRIC to a tab-separated (.tsv) antibiotic resistance file that '
                                              'can be used in the mutation analysis.')
        d[self.key.pp_odir] = GuiEntity(label='Result file',
                                        desc='Path to the result .tsv file. It can be used for antibiotic resistance '
                                             'analysis in the mutation analysis and contains '
                                             'locus tags and gene names of antibiotic resistance genes.\n\n'
                                             'If the file does not exist, it is '
                                             'automatically created if the containing directory exists. It can also '
                                             'be created via the file dialog.\n\n'
                                             'File type: tab-separated (.tsv)')
        d[self.key.pp_in] = GuiEntity(label='PATRIC input file',
                                      desc='An antibiotic resistance .csv file from PATRIC that contains '
                                           'locus tags, gene names and antibiotic resistance information.\n\n'
                                           'File type: comma-separated (.csv)')
        d[self.key.pp_run] = GuiEntity(label='Run',
                                       desc='Converts the PATRIC file to an antibiotic resistance .tsv file that can '
                                            'be used for antibiotic resistance analysis in the mutation analysis.\n\n'
                                            'Warning: If the result file already exists, this will override its '
                                            'content!')
        d[self.key.pp_open] = GuiEntity(label='Open result directory',
                                        desc='If enabled, opens the directory containing the result file after '
                                             'successfully converting the PATRIC file.')
        # RegulonDB converter
        d[self.key.pr_title] = GuiEntity(label='RegulonDB Converter',
                                         desc='Converts several files from RegulonDB, a database for E. coli data, '
                                              'to files suitable for the mutation analysis.')
        d[self.key.pr_odir] = GuiEntity(label='Result Directory',
                                        desc='Path to the result directory that contains the converted files such as '
                                             'a gene information, regulation, transcription factor binding site and '
                                             'transcription factor sequence files.')
        d[self.key.pr_gene] = GuiEntity(label='Gene file',
                                        desc='A .tsv file with gene information such as name, start, end, strand, '
                                             'description and DNA sequence from RegulonDB.\n\n'
                                             'File type: tab-separated (.tsv)')
        d[self.key.pr_gene_reg] = GuiEntity(label='TF-gene regulation file',
                                            desc='A .tsv file with transcription factors - gene regulation information '
                                                 'from RegulonDB.\n\n'
                                                 'File type: tab-separated (.tsv)')
        d[self.key.pr_tf_reg] = GuiEntity(label='TF-TF regulation file',
                                          desc='A .tsv file with transcription factor - transcription factor '
                                               'regulation information from RegulonDB.\n\n'
                                               'File type: tab-separated (.tsv)')
        d[self.key.pr_operon] = GuiEntity(label='TF-operon regulation file',
                                          desc='A .tsv file with transcription factor - operon '
                                               'regulation information from RegulonDB. This adds the operons to the '
                                               'result gene information file if the corresponding input file is '
                                               'provided as well.\n\n'
                                               'File type: tab-separated (.tsv)')
        d[self.key.pr_msa] = GuiEntity(label='TF MSA file',
                                       desc='A text file with transcription factor sequences from RegulonDB.\n\n'
                                            'File type: text (.txt)')
        d[self.key.pr_run] = GuiEntity(label='Run',
                                       desc='Converts the RegulonDB files to files suitable for the mutation analysis.'
                                            '\n\n'
                                            'Warning: If the result file already exists, this will override its '
                                            'content!')
        d[self.key.pr_open] = GuiEntity(label='Open result directory',
                                        desc='If enabled, opens the directory containing the result file after '
                                             'successfully converting the RegulonDB files.')
        d[self.key.pr_skip] = GuiEntity(label='Only consider entries with strong evidence',
                                        desc='If enabled, ignores all operons, binding sites and regulatory '
                                             'information with a weak evidence confidence level.')
        d[self.key.pr_tfbs] = GuiEntity(label='TFBS information file',
                                        desc='A .tsv file with information on transcription factor binding sites '
                                             'from RegulonDB.\n\n'
                                             'File type: tab-separated (.tsv)')
        d[self.key.pr_prom] = GuiEntity(label='Promoter file',
                                        desc='A .tsv file with information on promoters from RegulonDB.\n\n'
                                             'File type: tab-separated (.tsv)')
        d[self.key.pr_tu] = GuiEntity(label='TU file',
                                      desc='A .tsv file with information on transcription units from RegulonDB.\n\n'
                                           'File type: tab-separated (.tsv)')
        d[self.key.pr_op] = GuiEntity(label='Operon file',
                                      desc='A .tsv file with information on operons from RegulonDB.\n\n'
                                           'File type: tab-separated (.tsv)')
        d[self.key.pr_tu_reg] = GuiEntity(label='TF-TU regulation file',
                                          desc='A .tsv file with transcription factor - transcription unit regulation '
                                               'informaiton from RegulonDB.\n\n'
                                               'File type: tab-separated (.tsv)')

        """ STATUS """
        d[self.key.s_title] = GuiEntity(label='Status',
                                        desc='Gives information on the intermediate steps of the currently running '
                                             'process, or of the last completed process.\n\n'
                                             'Skipped: The step was skipped due to failure in an earlier '
                                             'step or because it was not enabled by the user.\n\n'
                                             'Failed: The step failed. See error messages or the log.txt file in the '
                                             'specified output directory, or the user manual for more information.\n\n'
                                             'Time: The step succeeded.')

        """ OTHERS """
        d[self.key.g_save] = GuiEntity(label='Save',
                                       desc='Saves the settings to the user configuration and closes the settings '
                                            'window. Before saving, the validity of the new settings is checked.')
        d[self.key.g_cancel] = GuiEntity(label='Cancel',
                                         desc='Discards the new settings and closes the settings window.')
        d[self.key.g_desc] = GuiEntity(label='Description',
                                       desc='Place the mouse on a text, button or input field '
                                            'to obtain a brief description. Open the user manual under '
                                            '\'Help\' for more detailed explanations.')

        return d


class Config:
    def __init__(self, cfg_path):
        """
        Initialises the application configuration depending on the operating system.
        """
        """ OPERATING SYSTEM """
        # OS names
        self.windows = 'Windows'    # type: str
        self.linux = 'Linux'        # type: str
        self.mac = 'Mac OS X'       # type: str
        self.cygwin = 'Cygwin'      # type: str
        # True if the OS is supported
        self.os_supported = True    # type: bool
        # determine OS
        self.os = self._get_os()    # type: str

        """ LOAD CONFIGURATION FROM FILE """
        # path to the configuration file
        self.cfg_path = cfg_path                                                            # type: str
        # load the configuration file
        self.cfg = yaml.load(open(self.cfg_path, 'r'), Loader=yaml.Loader)                  # type: dict
        # get the user, default and GUI configuration depending on the operating system
        self.user = UserConfig(self.cfg[self.os]['User'])   # type: UserConfig
        self.gui = GuiConfig(self.cfg[self.os]['GUI'])      # type: GuiConfig
        self.run = RunConfig(self.cfg['Run']['User'])       # type: RunConfig

        """ LOGS """
        self.a_log = ALog()                 # type: ALog
        self.n_log = None                   # type: NGSLog

        """ CONSTANTS """
        # user manual
        self.manual = 'user_manual.pdf'          # type: str
        self.installation = 'installation_guide.pdf'
        # output paths
        self.op = OutputPaths(self.user.ao_dir, self.user.no_dir)   # type: OutputPaths
        # file formats
        self.gene = GeneTSV                 # type: GeneTSV
        self.mut = MutationTSV              # type: MutationTSV
        self.reg = RegulationTSV            # type: RegulationTSV
        self.int = InterestTSV              # type: InterestTSV
        self.gm = GeneMutationList          # type: GeneMutationList
        self.prot = ProteinDomainTSV        # type: ProteinDomainTSV
        self.tf = TfbsTSV                   # type: TfbsTSV
        # database formats
        self.gdb = GeneDB                   # type: GeneDB
        self.mdb = MutationDB               # type: MutationDB
        # input parsing
        self.isup = SupportedInput          # type: SupportedInput
        # miscellaneous
        self.misc = Misc                    # type: Misc
        # genes of interest categories
        self.cat = set()                    # type: set

    def _get_os(self):
        """
        Determines the operating system running the program.
        :return: operating system
        """
        if platform.startswith('linux'):
            return self.linux
        elif platform == 'darwin':
            return self.mac
        elif platform == 'cygwin':
            return self.cygwin
        elif platform == 'win32':
            return self.windows
        else:
            self.os_supported = False
            return ''

    def update_output(self, a_out_dir, n_out_dir):
        """
        Updates the analysis output paths.
        :param a_out_dir: new analysis output directory
        :param n_out_dir: new NGS output directory
        """
        self.op = OutputPaths(a_out_dir, n_out_dir)

    def restore_default(self):
        """
        Restores the default settings.
        """
        # restore the default configuration
        self.cfg[self.os]['User'] = self.cfg[self.os]['Default']
        self.cfg['Run']['User'] = self.cfg['Run']['Default']
        self.user = UserConfig(self.cfg[self.os]['User'])
        self.run = RunConfig(self.cfg['Run']['User'])
        # adjust output paths
        self.update_output(self.user.ao_dir, self.user.no_dir)

    def write_config(self):
        """
        Writes the configuration dictionary into the configuration.yaml file.
        """
        # write the configuration to disk
        with open(self.cfg_path, 'w') as cfg_file:
            # prevent aliases in the YAML file
            dumper = yaml.SafeDumper
            dumper.ignore_aliases = lambda dump, data: True
            cfg_file.write(yaml.dump(self.cfg, default_flow_style=False, Dumper=dumper))

    def save_config(self):
        """
        Saves the current settings to the configuration file.
        """
        # update the user configuration dict
        self.cfg[self.os]['User'] = self.user.update_cfg_dict()
        self.cfg['Run']['User'] = self.run.update_cfg_dict()
        # write the configuration to disk
        self.write_config()

    def save_regulon(self):
        """
        Saves the current RegulonDB settings to the configuration file.
        """
        self.cfg[self.os]['User']['RegulonDB'] = self.user.update_regulon()
        self.write_config()

    def save_patric(self):
        """
        Saves the current PATRIC settings to the configuration file.
        """
        self.cfg[self.os]['User']['PATRIC'] = self.user.update_patric()
        self.write_config()

    def save_uniprot(self):
        """
        Saves the current UniProt settings to the configuration file.
        """
        self.cfg[self.os]['User']['UniProt'] = self.user.update_uniprot()
        self.write_config()

    def save_merger(self):
        """
        Saves the current VCF Merger settings to the configuration file.
        """
        self.cfg[self.os]['User']['VCF Merger'] = self.user.update_merger()
        self.write_config()


""" INITIALISATION """
try:
    # load the configuration
    cfg = Config('config.yaml')
    gui = cfg.gui
except Exception:
    print(os.getcwd())
    messagebox.showerror('Configuration Loading Error',
                         'An unhandled exception occurred while loading the user configuration.\n\n'
                         '{0}'.format(traceback.format_exc()))