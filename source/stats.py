# import standard or third party modules
import numpy as np
from scipy.stats import ranksums
from source.plots import BarPlot, CodingPlot, TfbsPlot, MdPlot, Result

# import own modules
from source.configuration import cfg

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2017 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'


class Stats:
    """
    Class that combines and stores statistical information on the database, including distributions and
    p-values.
    """
    def __init__(self):
        self.categories = sorted(cfg.cat)
        self.total = 'total'
        self.all_interest = 'total ' + cfg.user.as_int
        self.all_categories = self.categories + [self.all_interest, cfg.user.as_not_int, self.total, ]

        """ STATS """
        # actual and randomly generated TFBS scores
        keys = ['random', 'actual']
        subcategories = [('nbr', 0), ('scores', [])]
        self.tfbs = {cfg.user.as_int: self._build_dict(keys, subcategories),
                     cfg.user.as_not_int: self._build_dict(keys, subcategories)}
        self.tfbs_p = {k: -1 for k in ['actual', cfg.user.as_int, cfg.user.as_not_int]}

        # coding region substitution scores
        keys = [cfg.misc.missense, cfg.misc.nonsense, cfg.misc.readthrough, cfg.misc.syn, cfg.misc.frame_shift]
        subcategories = [('nbr', 0), ('scores', [])]
        self.coding = {cfg.user.as_int: self._build_dict(keys, subcategories),
                       cfg.user.as_not_int: self._build_dict(keys, subcategories)}
        self.coding_p = {k: -1 for k in keys}

        # mutation densities
        keys += [cfg.misc.tfbs, cfg.misc.prom]
        self.md = {cfg.user.as_int: {key: [] for key in keys},
                   cfg.user.as_not_int: {key: [] for key in keys}}
        self.md_p = {k: -1 for k in keys}

        # mutation distribution
        keys += [cfg.misc.transition, cfg.misc.transversion, cfg.misc.indel]
        self.dist = {cfg.user.as_int: {key: [] for key in keys},
                     cfg.user.as_not_int: {key: [] for key in keys}}
        self.dist_p = {k: -1 for k in keys}

        # general statistics
        keys = self.categories + [self.all_interest, cfg.user.as_not_int, self.total]
        subcategories = [('names', set()), ('lts', set()), ('genes nbr', 0), ('tfbs nbr', 0), ('tf nbr', 0),
                         ('coding ids', set()), ('prom ids', set()), ('tfbs ids', set()), ('total ids', set()),
                         ('coding nbr', 0), ('prom nbr', 0), ('tfbs mut nbr', 0), ('total nbr', 0)]
        self.general = self._build_dict(keys, subcategories)

        """ FLAGS """
        self.run_dist = True
        self.run_tfbs = True
        self.run_type_dist = True
        self.run_md = True
        self.run_coding = True

        """ FORMATTING """
        self._length = 120
        self._line_format = '{:-<' + str(self._length) + '}\n'
        self._header_line = self._line_format.format('')
        self._title_format = '{:^' + str(self._length) + '}\n'

    @staticmethod
    def _get_new_empty(a):
        if isinstance(a, set):
            return set()
        if isinstance(a, list):
            return []
        if isinstance(a, dict):
            return {}
        return int(str(a))

    def _build_dict(self, keys, subcategories):
        return {key: {k: self._get_new_empty(v) for k, v in subcategories} for key in keys}

    def compute(self, db):
        """
        Computes number of statistics such as the number of different genes and mutations, the mutation density,
        mutation distributions and p-values.
        :param db: database with gene and mutation information
        """
        j = 0
        for gene in db.genes.values():
            gene.calculate_mutation_densities()

            # genes of interest
            if gene.interest in self.categories:
                interest = [gene.interest, self.all_interest, self.total]
                short = cfg.user.as_int
            else:
                interest = [cfg.user.as_not_int, self.total]
                short = cfg.user.as_not_int

            for i in interest:
                # gene locus tag and name
                self.general[i]['names'].add(gene.get_locus_tag_and_gene_name())
                self.general[i]['lts'].add(gene.locus_tag)

                # count genes
                self.general[i]['genes nbr'] += 1

                # count how many genes are transcription factors
                if gene.regulated:
                    self.general[i]['tf nbr'] += 1

                # number of mutations
                self.general[i]['coding ids'].update(gene.coding_mut)
                self.general[i]['prom ids'].update(gene.prom_mut)
                self.general[i]['tfbs ids'].update(gene.tfbs_mut)
                self.general[i]['total ids'].update(gene.mutations)

            # mutation densities
            self.md[short][cfg.misc.missense].append(gene.missense_md)
            self.md[short][cfg.misc.nonsense].append(gene.nonsense_md)
            self.md[short][cfg.misc.readthrough].append(gene.readthrough_md)
            self.md[short][cfg.misc.frame_shift].append(gene.frame_shift_md)
            self.md[short][cfg.misc.syn].append(gene.syn_md)
            self.md[short][cfg.misc.tfbs].append(gene.tfbs_md)
            self.md[short][cfg.misc.prom].append(gene.prom_md)

        if cfg.run.a_int:
            # p-values antibiotic vs non-antibiotic
            self.md_p[cfg.misc.missense] = ranksums(self.md[cfg.user.as_not_int][cfg.misc.missense], 
                                                    self.md[cfg.user.as_int][cfg.misc.missense])[1]
            self.md_p[cfg.misc.nonsense] = ranksums(self.md[cfg.user.as_not_int][cfg.misc.nonsense],
                                                    self.md[cfg.user.as_int][cfg.misc.nonsense])[1]
            self.md_p[cfg.misc.readthrough] = ranksums(self.md[cfg.user.as_not_int][cfg.misc.readthrough],
                                                       self.md[cfg.user.as_int][cfg.misc.readthrough])[1]
            self.md_p[cfg.misc.frame_shift] = ranksums(self.md[cfg.user.as_not_int][cfg.misc.frame_shift],
                                                       self.md[cfg.user.as_int][cfg.misc.frame_shift])[1]
            self.md_p[cfg.misc.syn] = ranksums(self.md[cfg.user.as_not_int][cfg.misc.syn],
                                               self.md[cfg.user.as_int][cfg.misc.syn])[1]
            self.md_p[cfg.misc.prom] = ranksums(self.md[cfg.user.as_not_int][cfg.misc.prom],
                                                self.md[cfg.user.as_int][cfg.misc.prom])[1]
            self.md_p[cfg.misc.tfbs] = ranksums(self.md[cfg.user.as_not_int][cfg.misc.tfbs],
                                                self.md[cfg.user.as_int][cfg.misc.tfbs])[1]
            
        # compute the number of transcription factor binding sites
        for tfbs in db.tfbs.values():
            self.general[self.total]['tfbs nbr'] += 1

            interest = False

            for c in self.categories:
                if self.general[c]['lts'].intersection(tfbs['TFBS']):
                    self.general[c]['tfbs nbr'] += 1
                    interest = True
            if interest:
                self.general[self.all_interest]['tfbs nbr'] += 1
            if self.general[cfg.user.as_not_int]['lts'].intersection(tfbs['TFBS']):
                self.general[cfg.user.as_not_int]['tfbs nbr'] += 1

        # compute the number of mutations via the mutation IDs stored in the genes
        for c in self.all_categories:
            self.general[c]['coding nbr'] = len(self.general[c]['coding ids'])
            self.general[c]['prom nbr'] = len(self.general[c]['prom ids'])
            self.general[c]['tfbs mut nbr'] = len(self.general[c]['tfbs ids'])
            self.general[c]['total nbr'] = len(self.general[c]['total ids'])

        # compute mutation distributions
        for mut in db.mutations.values():
            # skip mutations in intergenic regions
            if mut.get_regions() == cfg.misc.intergenic:
                continue
            
            # dictionary to store distribution information associated with (non-) antibiotic resistance
            subcategories = [(cfg.misc.missense, 0), (cfg.misc.nonsense, 0), (cfg.misc.readthrough, 0),
                             (cfg.misc.frame_shift, 0), (cfg.misc.syn, 0), (cfg.misc.tfbs, 0), (cfg.misc.prom, 0),
                             (cfg.misc.transition, 0), (cfg.misc.transversion, 0), (cfg.misc.indel, 0)]
            temp = self._build_dict([cfg.user.as_int, cfg.user.as_not_int], subcategories)

            # whether the mutation occurs in antibiotic or non-antibiotic resistance genes (or both)
            interest = False
            non = False

            # iterate over all gene information associated with the mutation
            for r in mut.regions.values():
                # determine if the gene information belongs to gene of interest or not
                keys = set()
                if r.int:
                    keys.add(cfg.user.as_int)
                    interest = True
                if r.non:
                    keys.add(cfg.user.as_not_int)
                    non = True

                # sore the information in the respective dictionary entry
                for key in keys:
                    if r.coding:
                        temp[key][r.coding_effect] = 1
                    elif r.tfbs:
                        temp[key][cfg.misc.tfbs] = 1
                    elif r.prom:
                        temp[key][cfg.misc.prom] = 1

            keys = set()
            if interest:
                temp[cfg.user.as_int][mut.type] = 1
                keys.add(cfg.user.as_int)
            if non:
                temp[cfg.user.as_not_int][mut.type] = 1
                keys.add(cfg.user.as_not_int)
            
            # add the information to the corresponding distribution(s)       
            for k in keys:
                self.dist[k][cfg.misc.missense].append(temp[k][cfg.misc.missense])
                self.dist[k][cfg.misc.nonsense].append(temp[k][cfg.misc.nonsense])
                self.dist[k][cfg.misc.readthrough].append(temp[k][cfg.misc.readthrough])
                self.dist[k][cfg.misc.frame_shift].append(temp[k][cfg.misc.frame_shift])
                self.dist[k][cfg.misc.syn].append(temp[k][cfg.misc.syn])
                self.dist[k][cfg.misc.prom].append(temp[k][cfg.misc.prom])
                self.dist[k][cfg.misc.tfbs].append(temp[k][cfg.misc.tfbs])
                self.dist[k][cfg.misc.transition].append(temp[k][cfg.misc.transition])
                self.dist[k][cfg.misc.transversion].append(temp[k][cfg.misc.transversion])
                self.dist[k][cfg.misc.indel].append(temp[k][cfg.misc.indel])
        if cfg.run.a_int:
            # compute mutation distribution p-values for of interest versus not of interest
            self.dist_p[cfg.misc.missense] = ranksums(self.dist[cfg.user.as_not_int][cfg.misc.missense],
                                                      self.dist[cfg.user.as_int][cfg.misc.missense])[1]
            self.dist_p[cfg.misc.nonsense] = ranksums(self.dist[cfg.user.as_not_int][cfg.misc.nonsense],
                                                      self.dist[cfg.user.as_int][cfg.misc.nonsense])[1]
            self.dist_p[cfg.misc.readthrough] = ranksums(self.dist[cfg.user.as_not_int][cfg.misc.readthrough],
                                                         self.dist[cfg.user.as_int][cfg.misc.readthrough])[1]
            self.dist_p[cfg.misc.frame_shift] = ranksums(self.dist[cfg.user.as_not_int][cfg.misc.frame_shift],
                                                         self.dist[cfg.user.as_int][cfg.misc.frame_shift])[1]
            self.dist_p[cfg.misc.syn] = ranksums(self.dist[cfg.user.as_not_int][cfg.misc.syn],
                                                 self.dist[cfg.user.as_int][cfg.misc.syn])[1]
            self.dist_p[cfg.misc.prom] = ranksums(self.dist[cfg.user.as_not_int][cfg.misc.prom],
                                                  self.dist[cfg.user.as_int][cfg.misc.prom])[1]
            self.dist_p[cfg.misc.tfbs] = ranksums(self.dist[cfg.user.as_not_int][cfg.misc.tfbs],
                                                  self.dist[cfg.user.as_int][cfg.misc.tfbs])[1]
            self.dist_p[cfg.misc.transition] = ranksums(self.dist[cfg.user.as_not_int][cfg.misc.transition],
                                                        self.dist[cfg.user.as_int][cfg.misc.transition])[1]
            self.dist_p[cfg.misc.transversion] = ranksums(self.dist[cfg.user.as_not_int][cfg.misc.transversion],
                                                          self.dist[cfg.user.as_int][cfg.misc.transversion])[1]
            self.dist_p[cfg.misc.indel] = ranksums(self.dist[cfg.user.as_not_int][cfg.misc.indel],
                                                   self.dist[cfg.user.as_int][cfg.misc.indel])[1]
        
    def set_coding_scores(self, scores):
        """
        Sets the substitution scores of mutations in coding regions.
        :param scores: dictionary with the scores
        """
        for key, val in scores.items():
            for k, v in val.items():
                self.coding[key][k]['scores'] = v
                self.coding[key][k]['nbr'] = len(v)

        if cfg.run.a_int:
            # compute p-values using the Wilcoxon rank sums test
            for k in [cfg.misc.missense, cfg.misc.nonsense, cfg.misc.readthrough, cfg.misc.frame_shift, cfg.misc.syn]:
                self.coding_p[k] = ranksums(self.coding[cfg.user.as_not_int][k]['scores'],
                                            self.coding[cfg.user.as_int][k]['scores'])[1]

    def set_tfbs_scores(self, scores):
        """
        Stores the TFBS analysis scores and computes p-values with the Wilcoxon ranksum test.
        :param scores: scores from the TFBS analysis
        :return:
        """
        for key, val in scores.items():
            for k, v in val.items():
                self.tfbs[key][k]['scores'] = v
                self.tfbs[key][k]['nbr'] = len(v)

        if cfg.run.a_int:
            # compute p-values using the Wilcoxon rank sums test
            self.tfbs_p['actual'] = ranksums(self.tfbs[cfg.user.as_not_int]['actual']['scores'],
                                             self.tfbs[cfg.user.as_int]['actual']['scores'])[1]
            self.tfbs_p[cfg.user.as_not_int] = ranksums(self.tfbs[cfg.user.as_not_int]['actual']['scores'],
                                                        self.tfbs[cfg.user.as_not_int]['random']['scores'])[1]
            self.tfbs_p[cfg.user.as_int] = ranksums(self.tfbs[cfg.user.as_int]['actual']['scores'],
                                                    self.tfbs[cfg.user.as_int]['random']['scores'])[1]

    def plot(self):
        """
        Generates the plot for coding region substitution scores, transcription factor binding site scores,
        mutation density and mutation distribution.
        """
        if not cfg.run.a_int:
            return

        self.plot_coding()
        self.plot_mutation_density()
        self.plot_mutation_distribution_type()
        self.plot_mutation_distribution_all()

        if cfg.run.a_tfbs:
            self.plot_tfbs()

    def write(self):
        """
        Write the results of the statistical analysis into the specified file.
        :return: True if no error occurred, False otherwise
        """
        try:
            result = Result('Times')

            self.gene_table(result)
            self.mutation_table(result)

            if cfg.run.a_int:
                self.gene_name_list(result)
                self.mutation_distribution_type(result)
                self.mutation_distribution_all(result)
                self.mutation_density(result)
                self.coding_region(result)

                if cfg.run.a_tfbs:
                    self.tfbs_result(result)

            result.output(cfg.op.result, 'F')

            return True
        
        except IOError:
            return False

    def plot_coding(self):
        """
        Generates the plot for coding region substitution scores.
        """
        if not cfg.run.a_coding:
            return
        
        i = cfg.user.as_int
        n = cfg.user.as_not_int
        keys = [cfg.misc.missense, cfg.misc.nonsense, cfg.misc.readthrough, cfg.misc.frame_shift]
        
        interest = [self.coding[i][k]['scores'] for k in keys]
        non = [self.coding[n][k]['scores'] for k in keys]
        pvals = [self.coding_p[k] for k in keys]

        labels = ['missense\n[{:,} and {:,}]'.format(self.coding[n][keys[0]]['nbr'], self.coding[i][keys[0]]['nbr']),
                  'nonsense\n[{:,} and {:,}]'.format(self.coding[n][keys[1]]['nbr'], self.coding[i][keys[1]]['nbr']),
                  'readthrough\n[{:,} and {:,}]'.format(self.coding[n][keys[2]]['nbr'], self.coding[i][keys[2]]['nbr']),
                  'frame shift\n[{:,} and {:,}]'.format(self.coding[n][keys[3]]['nbr'], self.coding[i][keys[3]]['nbr'])]

        if np.isnan(pvals).any():
            self.run_coding = False
            return

        CodingPlot(labels, interest, non, pvals, cfg.op.coding)

    def plot_tfbs(self):
        """
        Generates the plot for transcription factor binding site scores using position weight matrices.
        """
        i = cfg.user.as_int
        n = cfg.user.as_not_int
        
        non = [self.tfbs[n]['random']['scores'], self.tfbs[i]['random']['scores']]
        interest = [self.tfbs[n]['actual']['scores'], self.tfbs[i]['actual']['scores']]

        labels = ['{0}\n[{1:,} and {2:,}]'.format(cfg.user.as_not_int, self.tfbs[n]['random']['nbr'], 
                                                  self.tfbs[n]['actual']['nbr']),
                  '{0}\n[{1:,} and {2:,}]'.format(cfg.user.as_int, self.tfbs[i]['random']['nbr'], 
                                                  self.tfbs[i]['actual']['nbr'])]
        pvals = [self.tfbs_p[n], self.tfbs_p[i]]

        if np.isnan(pvals).any():
            self.run_tfbs = False
            return
        
        TfbsPlot(labels, interest, non, pvals, cfg.op.tfbs)

    def plot_mutation_density(self):
        """
        Generates the plot of the mutation density in antibiotic and non-antibiotic resistance genes.
        """
        if not cfg.run.a_coding:
            return
        
        i = cfg.user.as_int
        n = cfg.user.as_not_int

        labels = [cfg.misc.syn, cfg.misc.missense, cfg.misc.nonsense, cfg.misc.readthrough,
                  cfg.misc.frame_shift, cfg.misc.prom]
        
        interest = [self.md[i][k] for k in labels]
        non = [self.md[n][k] for k in labels]
        pvals = [self.md_p[k] for k in labels]

        if cfg.run.a_tfbs:
            interest += [self.md[i][cfg.misc.tfbs]]
            non += [self.md[n][cfg.misc.tfbs]]
            labels += [cfg.misc.tfbs]
            pvals += [self.md_p[cfg.misc.tfbs]]

        if np.isnan(pvals).any():
            self.run_md = False
            return

        MdPlot(labels, interest, non, pvals, cfg.op.md)

    def plot_mutation_distribution_type(self):
        """
        Generates the plot of the mutation distribution for type.
        """
        i = cfg.user.as_int
        n = cfg.user.as_not_int
        
        labels = [cfg.misc.transition, cfg.misc.transversion, cfg.misc.indel]
        
        interest = [self.dist[i][k] for k in labels]
        non = [self.dist[n][k] for k in labels]
        pvals = [self.dist_p[k] for k in labels]

        if np.isnan(pvals).any():
            self.run_type_dist = False
            return

        BarPlot(labels, interest, non, pvals, cfg.op.mut_type, self.general[self.all_interest]['total nbr'],
                self.general[cfg.user.as_not_int]['total nbr'])

    def plot_mutation_distribution_all(self):
        """
        Generates the plot of the mutation distribution for type.
        """
        # do nothing without coding region analysis
        if not cfg.run.a_coding:
            return

        i = cfg.user.as_int
        n = cfg.user.as_not_int
        
        keys = [cfg.misc.syn, cfg.misc.missense, cfg.misc.nonsense, cfg.misc.readthrough, cfg.misc.frame_shift, 
                cfg.misc.prom]

        interest = [self.dist[i][k] for k in keys]
        non = [self.dist[n][k] for k in keys]
        pvals = [self.dist_p[k] for k in keys]

        labels = ['coding\nsyn', 'coding\nmissense', 'coding\nnonsense', 'coding\nreadthrough',
                  'coding\nframe-shift', 'promoter']

        # add TFBS analysis information
        if cfg.run.a_tfbs:
            interest += [self.dist[i][cfg.misc.tfbs]]
            non += [self.dist[n][cfg.misc.tfbs]]
            labels += [cfg.misc.tfbs]
            pvals += [self.dist_p[cfg.misc.tfbs]]

        if np.isnan(pvals).any():
            self.run_dist = False
            return

        BarPlot(labels, interest, non, pvals, cfg.op.mut_dist, self.general[self.all_interest]['total nbr'],
                self.general[cfg.user.as_not_int]['total nbr'])

    @staticmethod
    def table_field(nbr, total=-1):
        """
        Computes the percentage and returns a string with the number and the percentage
        :param nbr: number 
        :param total: total number
        :return: a string with the number and the percentage
        """
        if total == -1:
            return '{0:,}'.format(nbr)

        field = '{0:,} ({1:.1f}%)'

        # avoid division by 0
        if total == 0:
            per = 0
        else:
            per = (nbr / total) * 100
            
        return field.format(nbr, per)
    
    def gene_table(self, result):
        """
        Write the gene number table.
        :param result: statistics PDF
        """
        title = 'Gene Table'

        data = [['gene type', 'genes', 'TFBSs', 'TFs']]

        for k in self.all_categories[:-1]:
            data += [[k,
                      self.table_field(self.general[k]['genes nbr'], self.general[self.total]['genes nbr']),
                      self.table_field(self.general[k]['tfbs nbr'], self.general[self.total]['tfbs nbr']),
                      self.table_field(self.general[k]['tf nbr'], self.general[self.total]['tf nbr'])]]

        data += [[self.total,
                  self.table_field(self.general[self.total]['genes nbr']),
                  self.table_field(self.general[self.total]['tfbs nbr']),
                  self.table_field(self.general[self.total]['tf nbr'])]]

        note = 'Number of genes in the database, number of transcription factors and number of transcription factor ' \
               'binding sites for different types of genes.'  \
               '\'Total {0}\' is the sum of {1}. \'Total\' is the sum of \'total {0}\' and \'{2}\'.\n\n' \
               'Transcription factor binding sites can belong to an entire operon instead of a single gene, and an ' \
               'operon can contain both \'{0}\' and \'{2}\'. Thus, the total ' \
               'number of transcription factor binding sites can be lower than the sum of the TFBS column.' \
               ''.format(cfg.user.as_int, ', '.join(self.categories), cfg.user.as_not_int)

        result.table_page(title, note, data)

    def mutation_table(self, result):
        """
        Write the mutation number table.
        :param result: statistics PDF
        """

        note = 'Number of mutations in different types of genes (rows) and gene regions (columns). ' \
               '\'Total {0}\' is the sum of {1}. \'Total\' is the sum of \'total {0}\' and \'{2}\'. TFBS stands for ' \
               'transcription factor binding site.\n\n' \
               'Gene regions can overlap, including with regions of ' \
               'other genes. As a result, some mutations are located in both \'{0}\' and in \'{2}\'. In such a case, ' \
               'only mutation information associated with \'{0}\' is counted for the \'{0}\' rows, same for mutation ' \
               'information associated with \'{2}\' genes. The total numbers can be a bit smaller than the sum for ' \
               'that reason.'.format(cfg.user.as_int, ', '.join(self.categories), cfg.user.as_not_int)

        data = [['mutations in', 'in coding', 'in promoter', 'in TFBSs', 'total']]
        
        for k in self.all_categories[:-1]:
            data += [[k,
                      self.table_field(self.general[k]['coding nbr'], self.general[self.total]['coding nbr']),
                      self.table_field(self.general[k]['prom nbr'], self.general[self.total]['prom nbr']),
                      self.table_field(self.general[k]['tfbs mut nbr'], self.general[self.total]['tfbs mut nbr']),
                      self.table_field(self.general[k]['total nbr'], self.general[self.total]['total nbr'])]]

        data += [[self.total,
                  self.table_field(self.general[self.total]['coding nbr']),
                  self.table_field(self.general[self.total]['prom nbr']),
                  self.table_field(self.general[self.total]['tfbs mut nbr']),
                  self.table_field(self.general[self.total]['total nbr'])]]

        result.table_page('Mutation Table', note, data)

    def gene_name_list(self, result):
        """
        Write the list of gene names.
        :param result: statistics PDF
        """
        lists = []

        for k in self.categories:
            if self.general[k]['names']:
                lists.append(['\n{0}'.format(k), '\n-  '.join(sorted(self.general[k]['names']))])

        if lists:
            result.lists_page('Genes of Interest', lists)

    def mutation_density(self, result):
        """
        Write mutation density p-values into the statistics file.
        :param result: statistics PDF
        """
        if not cfg.run.a_coding or not self.run_dist:
            return

        title = 'Mutation Density'

        note = ['For each gene, the mutation density (#mutations / kbp) for different gene regions was computed. '
                'Mutations in coding regions were divided into four groups depending on their effect on the amino '
                'acid sequence of the protein. \'Missense\' means a single amino acid substitution, \'nonsense\' '
                'the replacement of an amino acid by a translation stop, \'readthrough\' a translation stop replaced '
                'by an amino acid and \'frame-shift\' a reading frame shift caused by an insertion or deletion.',

                'The plot compares the mutation density between \'{0}\' and \'{1}\'. At the top are the mean/median '
                'values. The blue circles represent outliers, the red squares the mean, and the middle line the '
                'median. \'*\' indicates a p-values < 0.05 and \'**\' a p-value < 0.01.'
                ''.format(cfg.user.as_int, cfg.user.as_not_int),

                'The Wilcoxon rank-sum test was used to statistically compare the mutation density distributions '
                'between \'{0}\' and \'{1}\'. Note that the test requires a sample size > 20 to be reliable.'
                ''.format(cfg.user.as_int, cfg.user.as_not_int)]

        pvals = [['', 'p-value'],
                 ['coding: synonymous', '{0:.8f}'.format(round(self.md_p[cfg.misc.syn], 8))],
                 ['coding: missense', '{0:.8f}'.format(round(self.md_p[cfg.misc.missense], 8))],
                 ['coding: nonsense', '{0:.8f}'.format(round(self.md_p[cfg.misc.nonsense], 8))],
                 ['coding: readthrough', '{0:.8f}'.format(round(self.md_p[cfg.misc.readthrough], 8))],
                 ['coding: frame-shift', '{0:.8f}'.format(round(self.md_p[cfg.misc.frame_shift], 8))],
                 ['promoter', '{0:.8f}'.format(round(self.md_p[cfg.misc.prom], 8))]]

        if cfg.run.a_tfbs:
            pvals += [['TFBS', '{0:.8f}'.format(round(self.md_p[cfg.misc.tfbs], 8))]]

        result.plot_page(title, cfg.op.md + '.png', note, pvals)
        
    def mutation_distribution_all(self, result):
        """
        Write mutation distribution p-values into the statistics file.
        :param result: statistics PDF
        """
        if not cfg.run.a_coding or not self.run_dist:
            return

        title = 'Mutation Distribution - Gene Regions'

        note = ['Mutations were divided based on their location in different gene regions. Mutations in coding regions '
                'were further divided into four groups depending on their effect on the amino acid sequence of the '
                'protein. \'Missense\' means a single amino acid substitution, \'nonsense\' '
                'the replacement of an amino acid by a translation stop, \'readthrough\' a translation stop replaced '
                'by an amino acid and \'frame-shift\' a reading frame shift caused by an insertion or deletion. '
                'Then, the number of mutations (in percent) in all \'{0}\' and in \'{1}\' was compared.'
                ''.format(cfg.user.as_int, cfg.user.as_not_int),

                'The plot shows the number of mutations (in percentage) in \'{0}\' and \'{1}\'. Brackets denote the '
                'number of mutations in the respective category of genes. The percentages can add up to more than '
                '100% due to an overlap of gene regions. '
                '\'*\' indicates a p-values < 0.05 and \'**\' a p-value < 0.01.'
                ''.format(cfg.user.as_int, cfg.user.as_not_int),

                'The Wilcoxon rank-sum test was used to statistically compare the \'{0}\' and \'{1}\' mutation '
                'distributions. Note that the test requires a sample size > 20 to be reliable.'
                ''.format(cfg.user.as_int, cfg.user.as_not_int)]

        pvals = [['', 'p-value'],
                 ['coding: synonymous', '{0:.8f}'.format(round(self.dist_p[cfg.misc.syn], 8))],
                 ['coding: missense', '{0:.8f}'.format(round(self.dist_p[cfg.misc.missense], 8))],
                 ['coding: nonsense', '{0:.8f}'.format(round(self.dist_p[cfg.misc.nonsense], 8))],
                 ['coding: readthrough', '{0:.8f}'.format(round(self.dist_p[cfg.misc.readthrough], 8))],
                 ['coding: frame-shift', '{0:.8f}'.format(round(self.dist_p[cfg.misc.frame_shift], 8))],
                 ['promoter', '{0:.8f}'.format(round(self.dist_p[cfg.misc.prom], 8))]]

        if cfg.run.a_tfbs:
            pvals += [['TFBS', '{0:.8f}'.format(round(self.dist_p[cfg.misc.tfbs], 8))]]

        result.plot_page(title, cfg.op.mut_dist + '.png', note, pvals)

    def mutation_distribution_type(self, result):
        """
        Write mutation distribution p-values into the statistics file.
        :param result: statistics PDF
        """
        if not self.run_type_dist:
            return

        title = 'Mutation Distribution - Type'

        note = ['Mutations were divided based on their type. Then, the number of mutations (in percent) in '
                '\'{0}\' and \'{1}\' was compared.'
                ''.format(cfg.user.as_int, cfg.user.as_not_int),

                'The plot show the number of mutations (in percentage) in \'{0}\' and \'{1}\'. Brackets denote the '
                'number of mutations in the respective category of genes. \'*\' indicates a p-values < 0.05 and '
                '\'**\' a p-value < 0.01.'
                ''.format(cfg.user.as_int, cfg.user.as_not_int),

                'The Wilcoxon rank-sum test was used to statistically compare the \'{0}\' and \'{1}\' mutation '
                'distributions. Note that the test requires a sample size > 20 to be reliable.'
                ''.format(cfg.user.as_int, cfg.user.as_not_int)]

        pvals = [['', 'p-value'],
                 ['transition', '{0:.8f}'.format(round(self.dist_p[cfg.misc.transition], 8))],
                 ['transversion', '{0:.8f}'.format(round(self.dist_p[cfg.misc.transversion], 8))],
                 ['indel', '{0:.8f}'.format(round(self.dist_p[cfg.misc.indel], 8))]]

        result.plot_page(title, cfg.op.mut_type + '.png', note, pvals, 0.5)

    def coding_region(self, result):
        """
        Write coding region substitution score p-values into the statistics file.
        :param result: statistics PDF
        """
        if not cfg.run.a_coding or not self.run_coding:
            return

        title = 'Coding Region Analysis'

        note = ['Mutations in coding regions were divided into four groups depending on their effect on the amino '
                'acid sequence of the protein. \'Missense\' means a single amino acid substitution, \'nonsense\' '
                'the replacement of an amino acid by a translation stop, \'readthrough\' a translation stop replaced '
                'by an amino acid and \'frame-shift\' a reading frame shift caused by an insertion or deletion.',

                'The mutations were scored with an amino acid substitution matrix to estimate their impact on the '
                'amino acid sequence and thus potentially protein function. The raw score was then normalised with the '
                'minimum and maximum score possible for the reference amino acid sequence, resulting in a score '
                'between 0 and 1. A normalised score of 0 means the mutated amino acid sequence differs as much as '
                'possible from the reference one, whereas a score of 1 means the mutated and reference sequences are '
                'identical (synonymous). Finally, the scores of mutations in \'{0}\' and \'{1}\' were '
                'compared.'.format(cfg.user.as_int, cfg.user.as_not_int),

                'The plot compares the scores of mutations in \'{0}\' and \'{1}\'. Brackets denote the number of '
                'mutations that were scored for each category. At the top are the mean/median values. The blue '
                'circles represent outliers, the red squares the mean, and the middle line the median. \'*\' '
                'indicates a p-values < 0.05 and \'**\' a p-value < 0.01.'.format(cfg.user.as_int, cfg.user.as_not_int),

                'The Wilcoxon rank-sum test was used to compute the statistical difference between mutations in '
                '\'{0}\' and \'{1}\'. Note that the test requires a sample size > 20 to be reliable.'
                ''.format(cfg.user.as_int, cfg.user.as_not_int)]

        pvals = [['', 'p-value'],
                 ['missense', '{0:.8f}'.format(round(self.coding_p[cfg.misc.missense], 8))],
                 ['nonsense', '{0:.8f}'.format(round(self.coding_p[cfg.misc.nonsense], 8))],
                 ['readthrough', '{0:.8f}'.format(round(self.coding_p[cfg.misc.readthrough], 8))],
                 ['frame-shift', '{0:.8f}'.format(round(self.coding_p[cfg.misc.frame_shift], 8))]]

        result.plot_page(title, cfg.op.coding + '.png', note, pvals)

    def tfbs_result(self, result):
        """
        Write transcription factor binding site score p-values into the statistics file.
        :param result: statistics PDF
        """
        if not self.run_tfbs:
            return

        title = 'Transcription Factor Binding Site Analysis'

        note = ['Mutations in transcription factor binding sites were scored with the position weight matrix of the '
                'corresponding transcription factor to ascertain if the mutation likely increases or decreases the '
                'ability of the transcription factor to bind to the mutated binding site. '
                'Additionally, for each mutation, {0} mutations were randomly generated in the binding site, such '
                'that the type of the mutation was conserved. I.e., if the observed mutation is a transition, then '
                'the randomly generated mutations are also transitions. These random mutations were also scored with '
                'the position weight matrix.'.format(cfg.user.as_rnd),

                'The scores of observed and random mutations were then normalised with the minimum and maximum score '
                'possible for the transcription factor binding site, resulting in a score between 0 and 1. A '
                'normalised score of 0 means that the sequence differs as much as possible from the transcription '
                'factor motif, whereas a normalised score of 1 means that the sequence is as similar as possible to '
                'the motif.',

                'The final score of observed and random mutations was computed as the difference from the score '
                'of the reference binding site sequence (score = score_ref - score_mut). Scores > 0 suggest that a '
                'mutation decreases the ability of the transcription factor to bind to the  mutated binding site, '
                'whereas scores < 0 suggest an increase. The more the score deviates from 0, he higher the difference '
                'between the reference and mutated binding site sequence and thus the higher he potential impact of '
                'the mutation.',

                'The plot shows the difference in score of reference and mutated binding site sequences for observed '
                'and randomly generated mutations, for both \'{0}\' and \'{1}\'. Brackets denote the number of '
                'mutations that were scored for each category. At the top are the mean/median values. The blue '
                'circles represent outliers, the red squares the mean, and the middle line the median. \'*\' '
                'indicates a p-values < 0.05 and \'**\' a p-value < 0.01.'.format(cfg.user.as_int, cfg.user.as_not_int),

                'The Wilcoxon rank-sum test was used to compute the statistical difference between observed and '
                'randomly generated mutations in \'{0}\' and \'{1}\'. Furthermore, the scores of observed mutations '
                'in \'{0}\' and \'{1}\' were statistically compared. Note that the test requires a sample size > 20 '
                'to be reliable.'.format(cfg.user.as_int, cfg.user.as_not_int)]

        pvals = [['', 'p-value'],
                 ['observed vs random (in \'{0}\')', '{1:.8f}'.format(cfg.user.as_int,
                                                                      round(self.tfbs_p[cfg.user.as_int], 8))],
                 ['observed vs random (in \'{0}\')', '{1:.8f}'.format(cfg.user.as_not_int,
                                                                      round(self.tfbs_p[cfg.user.as_not_int], 8))],
                 ['observed vs observed', '{0:.8f}'.format(round(self.tfbs_p['actual'], 8))]]

        result.plot_page(title, cfg.op.tfbs + '.png', note, pvals, 0.75)
