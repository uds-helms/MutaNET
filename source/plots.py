# import standard or third party modules
import fpdf
import math
import os.path as op
from datetime import datetime
from itertools import chain

# set up matplotlib
import matplotlib
import matplotlib as mpl
import numpy as np

matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# import own modules
from source.configuration import cfg

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2017 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'


class Plot:
    def __init__(self, height, weight, labels, int_data, non_data, pvals, file):
        # initialise the plot
        self.fig = plt.figure(figsize=(weight, height))
        self.plot = plt.subplot(111)
        # data
        self.labels = labels
        self.int = int_data  # of interest
        self.non = non_data  # not of interest
        self.pvals = pvals   # p-values
        # file path
        self.file = file
        # colours
        self.col_non = '#A0A0A0'  # not of interest
        self.col_int = '#1F497D'  # of interest
        # number of columns
        self.nbr = len(self.int)
        # minimum and maximum value in the combined data
        self.min, self.max = self.min_max()
        # set the general plot style
        self.base_style()

    def base_style(self):
        """
        Sets the basic base_style of a plot.
        """
        # set overall line base_style
        mpl.rcParams['lines.color'] = '#A0A0A0'
        mpl.rcParams['lines.linestyle'] = '-'
        mpl.rcParams['lines.linewidth'] = 1
        mpl.rcParams['lines.markersize'] = 4
        mpl.rcParams['lines.markeredgewidth'] = 0

        # set overall text and font base_style
        mpl.rcParams['text.color'] = '#464646'
        mpl.rcParams['font.family'] = 'serif'
        mpl.rcParams['font.size'] = 7

        # set overall axes base_style
        mpl.rcParams['axes.edgecolor'] = '#646464'
        mpl.rcParams['axes.labelcolor'] = '#464646'
        mpl.rcParams['axes.labelsize'] = 7
        mpl.rcParams['axes.titlesize'] = 7

        # set overall xtick and ytick base_style
        mpl.rcParams['xtick.color'] = '#464646'
        mpl.rcParams['ytick.labelsize'] = 7
        mpl.rcParams['ytick.color'] = '#464646'
        mpl.rcParams['ytick.labelsize'] = 7

        # set overall legend base_style
        mpl.rcParams['legend.frameon'] = False
        mpl.rcParams['legend.fontsize'] = 7

        # remove self.fig box lines
        self.plot.spines['top'].set_visible(False)
        self.plot.spines['bottom'].set_visible(False)
        self.plot.spines['right'].set_visible(False)
        self.plot.spines['left'].set_visible(False)

        # set the grid
        self.plot.yaxis.grid(True, linestyle='-', which='major', color='#A0A0A0', alpha=0.5, zorder=0)

        # remove ticks and set tick size
        plt.tick_params(axis='both', which='both', labelsize=7,
                        bottom='off', top='off', left='off', right='off',
                        labelleft='on', labelbottom='on')

        # set y-axis label size
        self.plot.yaxis.label.set_size(8)

    def save(self):
        """
        Saves the plot as .pdf and as .png.
        """
        self.fig.savefig(self.file + '.pdf', dpi=900, format='pdf', bbox_inches='tight', pad_inches=0.01)
        self.fig.savefig(self.file + '.png', format='png', bbox_inches='tight', pad_inches=0.01, dpi=300)

    def min_max(self):
        """
        Computes the minimum and maximum value in the combined data.
        """
        if not self.int and not self.non:
            return 0, 0

        try:
            combined = list(chain.from_iterable(self.int + self.non))
        except TypeError:
            combined = self.int + self.non

        if not combined:
            return 0, 0

        return min(combined), max(combined)


class BarPlot(Plot):
    def __init__(self, labels, int_data, non_data, pvals, file, tint, tnon):
        # set the plot dimensions
        height = 2.5
        width = 0.5 + len(int_data) * 0.75

        # adjust the data
        for i in range(len(int_data)):
            int_data[i] = np.mean(int_data[i]) * 100
            non_data[i] = np.mean(non_data[i]) * 100

        # initialise the plot
        Plot.__init__(self, height, width, labels, int_data, non_data, pvals, file)

        self.w = 0.3  # bar width
        self.wa = self.w - 0.05  # bar width adjusted
        self.i = np.arange(self.nbr) + 0.1  # bar indices

        # plot the bars
        self.lb = self.plot.bar(self.i, self.non, self.wa, color=self.col_non, zorder=3, edgecolor='none')
        self.rb = self.plot.bar(self.i + self.w, self.int, self.wa, color=self.col_int, zorder=3, edgecolor='none')

        # set the x and y axis
        self.plot.set_ylabel('number of mutations [%]')
        self.ylim()
        self.adjust_labels()
        self.plot.set_xticks(self.i + self.w)
        self.plot.set_xticklabels(self.labels, ha='center')

        # legend
        self.plot.legend((self.lb[0], self.rb[0]), ('{0} [{1:,}]'.format(cfg.user.as_not_int_short, tnon),
                                                    '{0} [{1:,}]'.format(cfg.user.as_int_short, tint)))

        # label the bars with the percentage
        self.label_bars()

        # save the plot
        self.save()

    def ylim(self):
        """
        Sets the y-limit.
        """
        if self.max > 85:
            self.plot.set_ylim(0, 110)
        else:
            self.plot.set_ylim(0, 100)

    def adjust_labels(self):
        """
        Adds a p-value marker to each label.
        """
        for i in range(self.nbr):
            if self.pvals[i] < 0.001:
                self.labels[i] += '\n**'
            elif self.pvals[i] < 0.05:
                self.labels[i] += '\n*'

    def label_bars(self):
        """
        Places a label with the percentage above each bar.
        """
        for bar in self.lb + self.rb:
            value = '{:<.1f}'.format(bar.get_height())
            y = bar.get_height() + 1
            x = bar.get_x() + bar.get_width() / 2

            self.plot.text(x, y, value, ha='center', va='bottom')


class BoxPlot(Plot):
    def __init__(self, labels, int_data, non_data, pvals, file):
        # set the plot dimensions
        height = 2.5
        weight = 1.5 + len(int_data)

        # initialise the plot
        Plot.__init__(self, height, weight, labels, int_data, non_data, pvals, file)

        # merge the data sets into one
        self.data = []
        for i in range(self.nbr):
            self.data.append(self.non[i])
            self.data.append(self.int[i])

        mean_props = dict(marker='s', markerfacecolor='#C0504D', markeredgecolor='none', zorder=5)
        flier_props = dict(marker='o', markerfacecolor='#4F81BD', markersize=6, linestyle='none', alpha=0.5)

        # plot the data
        self.bp = self.plot.boxplot(self.data, patch_artist=True, showmeans=True, meanprops=mean_props,
                                    flierprops=flier_props, showfliers=True)

        # default values
        self.ymax = 0
        self.dif = 0

    def legend(self, l, r, anchor):
        """
        Plot the legend.
        """
        patch1 = mpatches.Patch(color=self.col_non, label=l)
        patch2 = mpatches.Patch(color=self.col_int, label=r)

        self.plot.legend(handles=[patch1, patch2], loc='upper center', bbox_to_anchor=(0.5, anchor), ncol=5)

    def mean_median(self, dec):
        # write the mean and median above the box plot
        # create a list of positions
        pos = np.arange(self.nbr * 2) + 1

        # creates a list of medians and means, and then labels
        medians = [np.median(s) for s in self.data]
        means = [np.mean(s) for s in self.data]
        labels = [str(np.round(means[i], dec)) + '/' + str(np.round(medians[i], dec)) for i in range(0, len(medians))]
        # write the labels
        for tick, label in zip(range(self.nbr * 2), self.plot.get_xticklabels()):
            self.plot.text(pos[tick], self.ymax - (self.ymax * 0.03), labels[tick], ha='center')

    def separators(self):
        """
        Place vertical separators between pairs of box plot.
        """
        for i in range(1, self.nbr):
            self.plot.axvline(x=(i * 2 + 0.5))

    def pval_inditators(self):
        """
        Places p-value indicators at the top of the plot.+
        """
        for i in range(0, self.nbr):
            if self.pvals[i] < 0.001:
                self.plot.text(i * 2 + 1.5, self.ymax - (self.dif * 0.05), '**', ha='center')
            elif self.pvals[i] < 0.05:
                self.plot.text(i * 2 + 1.5, self.ymax - (self.dif * 0.1), '*', ha='center')

    def shrink(self, box_val):
        # shrink current axis's height by 10% on the bottom
        box = self.plot.get_position()
        self.plot.set_position([box.x0, box.y0 + box.height * (1 - box_val), box.width, box.height * box_val])

    def boxes(self):
        """
        Set up the boxes.
        """
        # set up the boxes
        counter = 1
        for box in self.bp['boxes']:
            if counter % 2:
                box.set(color=self.col_non, zorder=3)
            else:
                box.set(color=self.col_int, zorder=3)
            box.set(facecolor='white')
            counter += 1
        # set up the whiskers
        counter = 1
        for whisker in self.bp['whiskers']:
            if counter > 4:
                counter = 1
            if counter <= 2:
                whisker.set(color=self.col_non, linestyle='-', zorder=3)
            else:
                whisker.set(color=self.col_int, linestyle='-', zorder=3)
            counter += 1
        # set up the caps
        counter = 1
        for cap in self.bp['caps']:
            if counter > 4:
                counter = 1
            if counter <= 2:
                cap.set(color=self.col_non, zorder=3)
            else:
                cap.set(color=self.col_int, zorder=3)
            counter += 1
        # set up the medians
        counter = 1
        for median in self.bp['medians']:
            if counter % 2:
                median.set(color=self.col_non, linewidth=1.5, zorder=4)
            else:
                median.set(color=self.col_int, linewidth=1.5, zorder=4)
            counter += 1


class TfbsPlot(BoxPlot):
    def __init__(self, labels, res_data, non_data, pvals, file):
        BoxPlot.__init__(self, labels, res_data, non_data, pvals, file)
        # set the y-limit and y-label
        self.ylim()
        self.plot.set_ylabel('ScoreDif\n= ScoreRef - ScoreMut')
        # vertical separator lines
        self.separators()
        # compute mean and median and place them at the top
        self.mean_median(2)
        # place the p-value indicators
        self.pval_inditators()
        # set up the boxes
        self.boxes()
        # shrink the box at the bottom of the blot
        self.shrink(0.6)
        # plot the legend
        self.legend('random', 'observed', -0.22)
        # x-labels
        self.plot.set_xticks([1.5 + i * 2 for i in range(self.nbr)])
        self.plot.set_xticklabels(self.labels, ha='center')
        # save the plot as pdf and png
        self.save()

    def ylim(self):
        """
        Computes the y-limit based on the data.
        """
        dif = abs(self.max - self.min)
        self.ymax = 0.1 * math.ceil(float(self.max + dif * 0.2) / 0.1)
        ymin = 0.1 * math.floor(float(self.min) / 0.1)

        self.plot.set_ylim(ymin, self.ymax)
        self.dif = abs(self.ymax - ymin)


class CodingPlot(BoxPlot):
    def __init__(self, labels, res_data, non_data, pvals, file):
        BoxPlot.__init__(self, labels, res_data, non_data, pvals, file)
        # set the y-limit and y-label
        self.ylim()
        self.plot.set_ylabel('normalised amino\nacid substitution score')
        # place vertical separator lines
        self.separators()
        # compute mean and median and place them at the top
        self.mean_median(2)
        # place the p-value indicators
        self.pval_inditators()
        # set up the boxes
        self.boxes()
        # shrink the box at the bottom of the blot
        self.shrink(0.75)
        # plot the legend
        self.legend(cfg.user.as_not_int, cfg.user.as_int, -0.2)
        # x-labels
        self.plot.set_xticks([1.5 + i * 2 for i in range(self.nbr)])
        self.plot.set_xticklabels(self.labels, ha='center')
        # save the plot as pdf and png
        self.save()

    def ylim(self):
        """
        Computes the y-limit based on the data.
        """
        self.ymax = 1.1
        ymin = 0.0

        self.plot.set_ylim(ymin, self.ymax)
        self.dif = abs(self.ymax - ymin)


class MdPlot(BoxPlot):
    def __init__(self, labels, res_data, non_data, pvals, file):
        BoxPlot.__init__(self, labels, res_data, non_data, pvals, file)
        # set the y-limit and y-label
        self.ylim()
        self.plot.set_ylabel('#mutations / kbp')
        # place vertical separator lines
        self.separators()
        # compute mean and median and place them at the top
        self.mean_median(1)
        # place the p-value indicators
        self.pval_inditators()
        # set up boxes
        self.boxes()
        # shrink the box at the bottom of the blot
        self.shrink(0.9)
        # plot the legend
        self.legend(cfg.user.as_not_int, cfg.user.as_int, -0.125)
        # x-labels
        self.plot.set_xticks([1.5 + i * 2 for i in range(self.nbr)])
        self.plot.set_xticklabels(self.labels, ha='center')
        # save the plot as pdf and png
        self.save()

    def ylim(self):
        """
        Computes the y-limit based on the data.
        """
        dif = abs(self.max - self.min)
        self.ymax = 20 * math.ceil(float(self.max + dif * 0.2) / 20)
        ymin = 20 * math.floor(float(self.min) / 20)

        self.plot.set_ylim(ymin, self.ymax)
        self.dif = abs(self.ymax - ymin)


class Result(fpdf.FPDF):
    """
    Result PDF file.
    """
    def __init__(self, font):
        fpdf.FPDF.__init__(self, format='A4', unit='cm')
        # margins
        self.l_margin = 2.5
        self.r_margin = 2.5
        self.t_margin = 2.5
        self.b_margin = 1.5
        self.epw = self.w - 2 * self.l_margin
        # page numbers
        self.pn = 0
        # fonts and colours
        self.blue = (31, 73, 125)
        self.black = (70, 70, 70)
        self.red = (192, 80, 77)
        self.green = (118, 146, 61)
        self.f = font
        self.title_f = (self.f, 'B', 14)
        self.text_f = (self.f, '', 10)
        # title page
        self.title_page()

    def footer(self):
        """
        Set up the footer for each page.
        """
        self.pn += 1
        # show page number on the bottom on each page apart from the title page
        if self.pn > 1:
            self.set_y(-1.5)
            self.set_font(self.f, '', 8)
            self.set_text_color(*self.black)
            self.cell(self.epw, 0, '- {} -'.format(self.pn), align='C')

    def centre_start(self, w):
        """
        Give the start position of an element that is supposed to be centred.
        :param w: width of the element
        :return: start position
        """
        return (self.epw - w) / 2 + self.l_margin

    def titlefy(self, name):
        """
        Properly capitalises titles.
        """
        words = name.split()

        for i in range(len(words)):
            if not words[i].isupper():
                words[i] = words[i].title()

        return ' '.join(words)

    def title_page(self):
        """
        Sets up the title page.
        """
        self.add_page()
        self.ln(5)
        # pre title
        self.set_font(self.f, '', 20)
        self.set_text_color(*self.black)
        self.cell(self.epw, 0.0, self.titlefy('Mutation Analysis'), align='C')
        self.ln(1.5)
        # main title
        self.set_font(self.f, 'B', 28)
        self.set_text_color(*self.blue)
        self.cell(self.epw, 0.0, self.titlefy('Statistics Results'), align='C')
        self.ln(1.5)
        # date and time of the analysis
        self.set_font(self.f, '', 12)
        self.set_text_color(*self.black)
        date = datetime.today().strftime('Generated on %B %d, %Y - %H:%M')
        self.cell(self.epw, 0, date, align='C', ln=1)
        self.ln(3)
        # show files that have been used in the analysis
        self.path_line(cfg.gui.desc[cfg.gui.key.ao_dir].label, cfg.user.ao_dir)
        self.path_line(cfg.gui.desc[cfg.gui.key.ai_gene].label, cfg.user.ai_gene)
        self.path_line(cfg.gui.desc[cfg.gui.key.ai_mut].label, cfg.user.ai_mut)
        # show optional files that have been used in the analysis
        if cfg.run.a_int:
            self.path_line(cfg.gui.desc[cfg.gui.key.ai_int].label, cfg.user.ai_int)
        if cfg.run.a_coding:
            self.path_line(cfg.gui.desc[cfg.gui.key.ai_sm].label, cfg.user.ai_sm)
        if cfg.run.a_prot:
            self.path_line(cfg.gui.desc[cfg.gui.key.ai_prot].label, cfg.user.ai_prot)
        if cfg.run.a_reg:
            self.path_line(cfg.gui.desc[cfg.gui.key.ai_reg].label, cfg.user.ai_reg)
        if cfg.run.a_tfbs:
            self.path_line(cfg.gui.desc[cfg.gui.key.ai_tfbs].label, cfg.user.ai_tfbs)
            self.path_line(cfg.gui.desc[cfg.gui.key.ai_msa].label, cfg.user.ai_msa)

    def path_line(self, name, path):
        """
        Sets up a sub title with the file name and below a line with the file path.
        :param name: file name
        :param path: file path
        """
        self.set_text_color(*self.green)
        self.set_font(self.f, 'B', 12)
        self.cell(self.epw, 0, self.titlefy(name), align='C')
        self.ln(self.font_size)
        self.set_text_color(*self.black)
        self.set_font(self.f, '', 8)
        self.cell(self.epw, 0, op.abspath(path), align='C')
        self.ln(self.font_size * 3)

    def page_title(self, title):
        """
        Sets up a chapter/page title.
        """
        self.add_page()
        self.set_text_color(*self.blue)
        self.set_font(self.f, 'B', 18)
        self.cell(self.epw, 0.0, self.titlefy(title), align='C')
        self.ln(1)
        self.set_font(*self.text_f)
        self.set_text_color(*self.black)

    def sh(self, text, size):
        """
        Sets up a sub-title.
        :param text: title text
        :param size: title size
        """
        self.set_text_color(*self.red)
        self.set_font(self.f, 'B', size)
        self.cell(self.epw, 0, self.titlefy(text))
        self.set_text_color(*self.black)
        self.ln(self.font_size * 0.75)

    def description(self, title, desc):
        """
        Sets up a description section with title followed by text.
        """
        self.sh(self.titlefy(title), 12)
        self.set_font(self.f, '', 10)
        self.multi_cell(self.epw, self.font_size, desc)

    def table_page(self, title, desc, data):
        """
        Sets up a page with a title, table and description.
        """
        self.page_title(title)

        tso = 10
        tsi = 8

        padx_2 = self.epw / 6
        padx_1 = padx_2 * 1.5

        cols = len(data[0])
        rows = len(data)
        start = self.centre_start(padx_2 * (cols - 1) + padx_1)
        h = self.font_size * 2

        self.set_text_color(*self.black)
        self.set_font(self.f, 'B', tso)

        self.set_draw_color(*self.green)

        for r in range(rows):
            # first row
            if r == 0:
                bd = 'TB'
                self.set_line_width(0.075)
                self.set_font(self.f, 'B', tso)
            # last row
            elif r == rows - 1:
                bd = 'B'
                self.set_line_width(0.075)
                self.set_font(self.f, '', tsi)
            elif r == (rows - 4) or r == (rows - 2):
                bd = 'B'
                self.set_line_width(0.03)
                self.set_font(self.f, '', tsi)
            else:
                bd = ''
                self.set_font(self.f, '', tsi)

            for c in range(cols):
                if c == 0:
                    self.set_x(start)
                    self.cell(padx_1, h, str(data[r][c]), border=bd)
                else:
                    self.cell(padx_2, h, str(data[r][c]), border=bd, align='C')
            self.ln(h)

        self.ln(self.font_size * 5)
        self.description('Description', desc)

    def lists_page(self, title, lists):
        """
        Sets up a page with title and (several) list(s).
        """
        self.page_title(title)

        for t, l in lists:
            self.sh(t, 12)
            self.set_font(self.f, '', 10)
            self.multi_cell(self.epw, self.font_size, '-  ' + l)
            self.ln(self.font_size * 2)

    def plot_page(self, title, image, desc, pvals, scale=1):
        """
        Sets up a page with title, method description, plot, plot description, p-value description and p-value table.
        :param title: page title
        :param image: image path
        :param desc: list with method, plot and p-value description
        :param pvals: p-values
        :param scale: plot scale
        """
        self.page_title(title)
        self.description('Method', desc[0])
        self.ln(self.font_size * 2)
        self.sh('Plot', 12)
        w = self.epw * scale
        x = self.centre_start(w)
        self.image(image, x=x, w=w)
        self.ln(self.font_size)
        self.set_text_color(*self.black)
        self.set_font(self.f, '', 8)
        self.ln(self.font_size * 2)
        self.description('Description', desc[1])
        self.ln(self.font_size * 2)
        self.description('P-Values', desc[2])
        self.ln(self.font_size)
        tso = 10
        tsi = 8

        padx_2 = self.epw / 6.5

        max_len = max([len(x[0]) for x in pvals])

        if max_len >= 20:
            padx_1 = padx_2 * 1.85
        elif max_len >= 15:
            padx_1 = padx_2 * 1.25
        else:
            padx_1 = padx_2

        rows = len(pvals)
        start = self.centre_start(padx_2 + padx_1)
        h = self.font_size * 2

        self.set_text_color(*self.black)
        self.set_font(self.f, 'B', tso)
        self.set_draw_color(*self.green)

        for r in range(rows):
            if r == 0:
                bd = 'TB'
                self.set_line_width(0.075)
                self.set_font(self.f, 'B', tso)
            elif r == rows - 1:
                bd = 'B'
                self.set_line_width(0.075)
                self.set_font(self.f, '', tsi)
            else:
                bd = ''
                self.set_font(self.f, '', tsi)

            self.set_x(start)
            self.cell(padx_1, h, str(pvals[r][0]), border=bd)
            self.cell(padx_2, h, str(pvals[r][1]), border=bd, align='C')
            self.ln(h)