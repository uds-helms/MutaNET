# import standard or third party modules
from collections import defaultdict

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2017 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'


class LogSet:
    def __init__(self, text):
        self.set = set()
        self.text = text

    def __str__(self):
        if not self.set:
            return '0'
        return str(sorted(self.set)).strip('[]')

    def add(self, elem):
        self.set.add(elem)

    def print(self, verbose):
        if verbose:
            return '{0}: {1}'.format(self.text, str(self))
        else:
            return '{0}: {1}'.format(self.text, len(self.set))

    def __bool__(self):
        return self.set


class ALog:
    """
    Groups, stores and outputs potential problems that occurred during the runtime.
    """
    def __init__(self):
        # DATABASE - LOOKUP
        self.db_find_gene = LogSet('Genes were not found in database')
        self.db_duplicate_genes = LogSet('Genes were entered more than once')
        self.db_find_mut = LogSet('Mutations not in the database')

        # DATABASE - OPERONS
        self.db_operon_inconsistent = LogSet('Inconsistent operon direction')
        self.db_operon_unsupported = LogSet('Operon format was not supported')
        self.db_operon_with_duplicate_genes = LogSet('Operons contained duplicate genes')
        self.db_operon_not_in_database = LogSet('Operons contained genes that were not in the database')
        self.db_operon_inconsistent_strand = LogSet('Operon genes were not just on one strand')
        self.db_genes_in_several_operons = LogSet('Genes occurred in more than one operon')

        # DATABASE - READ TFBS
        self.tfbs_not_dna = LogSet('TFBS sequence was not a DNA sequence')
        self.tfbs_no_position = LogSet('TFBS did not have a valid start and end position')
        self.tfbs_multiple_dna = LogSet('TFBS had more than one DNA sequence')
        self.tfbs_no_unique_pos = LogSet('TFBS position could not be determined')
        self.tfbs_multiple_strand = LogSet('TFBS strand could not be determined')

        # DATABASE - TFBS ANALYSIS
        self.db_tfbs_no_tf_pwm = LogSet('TF corresponding to TFBS did not have a PWM')
        self.db_tfbs_mutation_wrong = LogSet('Mutations falsely associated with TFBS')
        self.db_tfbs_length = LogSet('TFBS length did not match PWM length')

        # PWM
        self.pwm_not_enough = LogSet('TF MSA did not contain enough valid sequences')
        self.pwm_non_dna = LogSet('TF MSA contained non-DNA sequences')
        self.pwm_length = LogSet('TF MSA contained sequences with different lengths')
        self.pwm_score_length = LogSet('TFBS sequence length does not match PWM length')

        # GENE
        self.gene_strand = LogSet('Unsupported strand format')
        # gene start and end
        self.gene_se_not_int = LogSet('Gene start or end was not an integer')
        self.gene_se_less_zero = LogSet('Gene start or end was < 0')
        self.gene_se_zero = LogSet('Gene start or end was 0')
        # DNA
        self.gene_dna_empty = LogSet('DNA sequence was not given')
        self.gene_not_dna = LogSet('DNA sequence contained non-DNA base(s)')
        self.gene_coding_length = LogSet('Coding region length did not match DNA length')
        # protein
        self.gene_protein_empty = LogSet('Protein sequence was not given')
        self.gene_not_protein = LogSet('Protein sequence contained non-amino acids')
        # promoter
        self.gene_prom_not_int = LogSet('Promoter information was not an integer')
        self.gene_prom_pos = LogSet('Relative promoter start was > 0')
        self.gene_prom_less_zero = LogSet('Absolute promoter information was < 0')
        # add mutations
        self.gene_mut_loc_int = LogSet('Add mutation by location - mutation ID was not an integer')
        self.gene_mut_loc_unsupported = LogSet('Add mutation by location - location was not supported')
        self.gene_mut_eff_int = LogSet('Add mutation by effect - mutation ID was not an integer')
        self.gene_mut_eff_unsupported = LogSet('Add mutation by effect - effect was not supported')
        # resistance
        self.gene_res = LogSet('Resistance information was not supported')
        # protein domains
        self.gene_prot_int = LogSet('Protein domain start or end was not an integer')
        self.gene_prot_pos_int = LogSet('Compute protein position - genome position is not an integer')
        self.gene_prot_pos_zero = LogSet('Compute protein position - genome position is less or equal 0')
        # regulation
        self.gene_regulation_unsupported = LogSet('Regulation description was not supported')

        # MUTATION
        # position
        self.mut_pos_int = LogSet('Position was not an integer')
        self.mut_pos_zero = LogSet('Position was smaller or equal 0')
        # type
        self.mut_type_dna = LogSet('Reference or alternative base was not DNA')
        # gene regions
        self.mut_region_unsupported = LogSet('Gene region was not supported')
        self.mut_coding_region = LogSet('Coding region was associated with more than one locus tag')
        self.mut_coding_not = LogSet('Coding region was not associated with the mutation')
        self.mut_tfbs_not = LogSet('TFBS region was not associated with the mutation')

    def write(self, path):
        """
        Writes the analysis log into the file.
        :param path: path to the output file
        """
        try:
            file = open(path, 'w')
            v = False
            sections = [['[DATABASE - LOOKUP]',
                         self.db_find_gene.print(v),
                         self.db_duplicate_genes.print(v),
                         self.db_find_mut.print(v)],
                        ['[DATABASE - OPERON PROCESSING]',
                         self.db_operon_unsupported.print(v),
                         self.db_operon_inconsistent.print(v),
                         self.db_operon_with_duplicate_genes.print(v),
                         self.db_operon_not_in_database.print(v),
                         self.db_operon_inconsistent_strand.print(v),
                         self.db_genes_in_several_operons.print(v)],
                        ['[DATABASE - READ TFBS]',
                         self.tfbs_not_dna.print(v),
                         self.tfbs_no_position.print(v),
                         self.tfbs_multiple_dna.print(v),
                         self.tfbs_no_unique_pos.print(v),
                         self.tfbs_multiple_strand.print(v)],
                        ['[DATABASE - TFBS ANALYSIS]',
                         self.db_tfbs_no_tf_pwm.print(v),
                         self.db_tfbs_mutation_wrong.print(v),
                         self.db_tfbs_length.print(v)],
                        ['[PWM]',
                         self.pwm_not_enough.print(v),
                         self.pwm_non_dna.print(v),
                         self.pwm_length.print(v),
                         self.pwm_score_length.print(v)],
                        ['[GENES]',
                         self.gene_strand.print(v),
                         '',
                         self.gene_se_not_int.print(v),
                         self.gene_se_less_zero.print(v),
                         self.gene_se_zero.print(v),
                         '',
                         self.gene_dna_empty.print(v),
                         self.gene_not_dna.print(v),
                         self.gene_coding_length.print(v),
                         '',
                         self.gene_protein_empty.print(v),
                         self.gene_not_protein.print(v),
                         '',
                         self.gene_prom_not_int.print(v),
                         self.gene_prom_pos.print(v),
                         self.gene_prom_less_zero.print(v),
                         '',
                         self.gene_mut_loc_int.print(v),
                         self.gene_mut_loc_unsupported.print(v),
                         self.gene_mut_eff_int.print(v),
                         self.gene_mut_eff_unsupported.print(v),
                         '',
                         self.gene_res.print(v),
                         '',
                         self.gene_prot_int.print(v),
                         '',
                         self.gene_regulation_unsupported.print(v),
                         '',
                         self.gene_prot_pos_int.print(v),
                         self.gene_prot_pos_zero.print(v)],
                        ['[MUTATIONS]',
                         self.mut_pos_int.print(v),
                         self.mut_pos_zero.print(v),
                         '',
                         self.mut_type_dna.print(v),
                         '',
                         self.mut_region_unsupported.print(v),
                         self.mut_coding_region.print(v),
                         self.mut_coding_not.print(v),
                         self.mut_tfbs_not.print(v)]]

            sections = ['\n'.join(sec) for sec in sections]

            file.write('\n\n{0:-<100}\n\n'.format('').join(sections))
            file.close()

            return True

        except IOError:
            return False


class NGSLog:
    def __init__(self, ref_path, reads_dir, qual, pval):
        self.ref_path = ref_path
        self.reads_dir = reads_dir
        self.qual = qual
        self.pval = pval
        # pipeline errors
        self.strains = set()
        # merger errors
        self.validate_errors = set()
        self.field_errors = set()
        self.parser_errors = defaultdict(set)
        self.io_errors = set()

    def _reads_to_string(self):
        """
        :return: string representation of the read errors
        """
        if not self.strains:
            return ''

        msg = 'The following strains did not have the required two .fastq files:'
        return '[READ FILES]\n{0}\n- {1}\n\n'.format(msg, '\n- '.join(sorted(self.strains)))

    def _merger_to_string(self):
        """
        :return: string representation of .vcf merger errors
        """
        if not self.validate_errors and not self.parser_errors and not self.field_errors and not self.io_errors:
            return ''

        msg = '[VCF FILE MERGER]'

        if self.io_errors:
            msg += 'The following .vcf files could not be opened or read:\n- {0}\n\n'
            msg = msg.format('\n- '.join(sorted(self.io_errors)))

        if self.field_errors:
            msg += 'The following .vcf files did not have the required POS, REF and ALT columns:\n- {0}\n\n'
            msg = msg.format('\n- '.join(sorted(self.field_errors)))

        if self.parser_errors:
            msg += 'The there are errors at the following mutation positions in the .vcf files:\n{0}\n\n'

            ret = ['{0} - {1}'.format(f, pos)
                   for f in sorted(self.parser_errors.keys()) for pos in sorted(self.parser_errors[f])]
            msg = msg.format('\n'.join(ret))

        if self.validate_errors:
            msg += 'The following positions have inconsistent references:\n- {0}\n\n'
            msg = msg.format('\n- '.join(sorted(self.validate_errors)))

        return msg + '\n\n'

    def write(self, path):
        """
        Writes the NGS log into the file.
        :param path: path to the output file
        """
        try:
            file = open(path, 'w')
            file.write('Reference genome: {0}\n'.format(self.ref_path))
            file.write('Reads directory:  {0}\n'.format(self.reads_dir))
            file.write('Quality cut-off:  {0}\n'.format(self.qual))
            file.write('P-value cut-off:  {0}\n\n'.format(self.pval))
            file.write(self._reads_to_string())
            file.write(self._merger_to_string())
            file.close()
            return True

        except IOError:
            return False
