# import standard or third party modules
import os
import shutil
from time import time
from shlex import quote
import subprocess as sp
from collections import defaultdict

# import own modules
from source.log import NGSLog
from source.converter import MutationVcfMerger
from source.configuration import cfg, adjust_dir_path, adjust_file_path
from source.tools import get_file_name, setup_directory, clean_directory

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2017 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'


class NGSPipeline:
    def __init__(self):
        # PATHS - absolute
        self.ref_path = adjust_file_path(os.path.abspath(cfg.user.ni_ref))
        self.reads_dir = adjust_dir_path(os.path.abspath(cfg.user.ni_reads))
        self.out_dir = adjust_dir_path(os.path.abspath(cfg.user.no_dir))

        # ERROR FLAGS
        self.ref_not_fasta = False
        self.reads_no_fastq = False
        self.reads_no_strains = False
        self.strain_false_nbr = set()
        self.inconsistent_suffixes = set()
        self.fail = False
        
        # INITIALISE LOG
        cfg.n_log = NGSLog(self.ref_path, self.reads_dir, cfg.user.ns_qcut, cfg.user.ns_pv)

        # PARSE INPUT - REF GENOME
        self.ref_file_name, self.ref_name, ref_suffix = get_file_name(self.ref_path)

        if ref_suffix != '.fasta':
            self.ref_not_fasta = True

        # PARSE INPUT - READS
        # initialise fields
        self.strains = dict()
        self.suffix = ''
        # obtain list of files in the reads directory
        read_file_names = os.listdir(self.reads_dir)
        # to see how many different
        suffixes = set()
        strains = defaultdict(set)

        # process the file names in the reads directory
        for path in read_file_names:
            # skip directories
            if os.path.isdir(self.reads_dir + path):
                continue

            complete_name, name, suffix = get_file_name(path)
            # skip files that are not .fastq
            if '.fastq' not in suffix:
                continue
            # name = strain_delim
            parts = name.split('_')
            strains['_'.join(parts[:-1])].add('_' + parts[-1])
            suffixes.add(suffix)

        # check if there is at least one .fastq in the reads folder and if the file type is consistent
        if len(suffixes) == 0:                      # there are no files in .fastq format in the reads directory
            self.reads_empty = True
        elif len(suffixes) > 1:
            self.inconsistent_suffixes = suffixes   # inconsistent suffixes
            self.suffix = ''
        else:
            self.suffix = list(suffixes)[0]

        # check if for each strain two .fastq files exist in the reads directory
        for key, value in strains.items():
            if len(value) != 2:
                cfg.n_log.strains.add(key)
            else:
                self.strains[key] = value
        
        # there are no strains with two .fastq files
        if not self.strains:
            self.reads_no_strains = True

    def _command_line(self, cmd_format, cmd_args):
        """
        Builds and executes the specified command on the command line.
        :param cmd_format: format of the command
        :param cmd_args: arguments
        :return: True if the call succeeded, False otherwise
        """

        # make the arguments safe for the shell and account for varscan potentially being a java executable
        if cmd_args[0][-4:].lower() == '.jar':
            cmd_args = list(map(quote, cmd_args))
            cmd_args[0] = 'java -jar {}'.format(cmd_args[0])
        else:
            cmd_args = list(map(quote, cmd_args))

        # adjust the command for different operating systems
        if cfg.os == cfg.windows:
            cmd_format = 'sh -c \"' + cmd_format + '\"'  # cygwin needs to be installed

        try:
            # build the command and execute it in the shell
            if cfg.os == cfg.windows:
                sp.check_call(args=cmd_format.format(*cmd_args), stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT,
                              shell=True)
            else:
                sp.check_call(args=cmd_format.format(*cmd_args), stdout=sp.PIPE, shell=True)
            return True
        except Exception as e:
            print(e)
            self.fail = True
            return False

    def run_pipeline(self, status):
        """
        Executes the NGS pipeline that indexes the reads, maps them onto the reference genome and calls SNPs and indels.
        :type status: Status
        """
        """ SETUP OUTPUT DIRECTORY """
        # create output directory for each strain
        for strain in self.strains.keys():
            strain_dir = self.out_dir + strain + '/'
            vcf_dir = strain_dir + cfg.op.vcf
            # this also creates out_path/strain/ if it does not exist
            if not setup_directory(strain, vcf_dir):
                return
            # copy the reference genome to the output directory
            ref_copy_path = strain_dir + self.ref_file_name
            if not os.path.exists(ref_copy_path):
                shutil.copy2(self.ref_path, strain_dir)
               
        """ PIPELINE """
        # number of strains for convenience
        strains = sorted(self.strains.keys())
        nbr_strains = len(strains)
        some_success = False

        # iterate over all strains given in the reads directory
        for i in range(nbr_strains):
            # update the GUI
            status.of(1, 'processing strain', i + 1, nbr_strains)
            status.clear(2)

            # set up a couple of file/directory paths for convenience
            strain = strains[i]
            strain_dir = self.out_dir + strain + '/'
            vcf_dir = strain_dir + cfg.op.vcf
            ref_path = strain_dir + self.ref_file_name
            ref_name = strain_dir + self.ref_name
            # list of read files belonging to the strain
            reads = [strain + nbr for nbr in self.strains[strain]]

            # build the BWA index
            if not self.fail:
                start2 = time()
                status.running(2, 'BWA (mapping)', True)
                cmd_format = '{0} index -a is -p {1} {2}'
                cmd_args = [cfg.user.ne_bwa, ref_name, ref_path]
                worked = self._command_line(cmd_format, cmd_args)

                # map the reads to the reference genome, resulting in out_directory/strain/read.sai files
                if worked:
                    for read in reads:
                        cmd_format = '{0} aln {1} {2} > {3}'
                        cmd_args = [cfg.user.ne_bwa, ref_name, self.reads_dir + read + self.suffix,
                                    strain_dir + read + '.sai']
                        worked = self._command_line(cmd_format, cmd_args) and worked
                        status.app.update()
                # update the GUI
                if worked:
                    status.time(2, 'BWA (mapping)', time() - start2, True)
                else:
                    status.fail(2, 'BWA (mapping)', True)
            else:
                status.skip(2, 'BWA (mapping)', True)

            if not self.fail:
                # join the out_directory/strain/read.sai files into .sam file
                start2 = time()
                status.running(3, 'BWA (to .sam)', True)
                join1 = []
                join2 = []
    
                # the full command is: bwa sampe ref_genome [list of .sai files] [list of .fastq files] > output.sam
                # it is constructed while constructing the lists of .sai and .fastq files
                nbr_reads = len(reads)
                cmd_format = '{0} sampe {1}'
    
                for j in range(nbr_reads):
                    read = reads[j]
                    join1.append(strain_dir + read + '.sai')
                    join2.append(self.reads_dir + read + self.suffix)
                    cmd_format += ' {' + str(j * 2 + 2) + '}'
                    cmd_format += ' {' + str(j * 2 + 3) + '}'
    
                # finish building the command and execute it
                next_index = 2 + nbr_reads * 2
                cmd_format += ' > {' + str(next_index) + '}'
                cmd_args = [cfg.user.ne_bwa, ref_name, *join1, *join2, strain_dir + strain + '.sam']
                worked = self._command_line(cmd_format, cmd_args)
            
                if worked:
                    status.time(3, 'BWA (to .sam)', time() - start2, True)
                else:
                    status.fail(3, 'BWA (to .sam)', True)
            else:
                status.skip(3, 'BWA (to .sam)', True)

            # convert to the .sam files to .bam files, sort and index
            # faidx: index and extract .fasta
            if not self.fail:
                start2 = time()
                status.running(4, 'SAMTools (to .bam + indexing)', True)
                
                cmd_format = '{0} view -bSh {1} |{2} sort -o {3}'
                cmd_args = [cfg.user.ne_samtools, strain_dir + strain + '.sam',
                            cfg.user.ne_samtools, strain_dir + strain + '_s.bam']
                worked = self._command_line(cmd_format, cmd_args)
                status.app.update()
            
                if worked:
                    cmd_format = '{0} faidx {1}'
                    cmd_args = [cfg.user.ne_samtools, ref_path]
                    worked = self._command_line(cmd_format, cmd_args)
                    status.app.update()
            
                if worked:
                    cmd_format = '{0} index {1}'
                    cmd_args = [cfg.user.ne_samtools, strain_dir + strain + '_s.bam']
                    self._command_line(cmd_format, cmd_args)

                # update the GUI
                if worked:
                    status.time(4, 'SAMTools (to .bam + indexing)', time() - start2, True)
                else:
                    status.fail(4, 'SAMTools (to .bam + indexing)', True)
            else:
                status.skip(4, 'SAMTools (to .bam + indexing)', True)
            
            """ SAMTools - QUALITY CONTROL """
            # filter reads with mapping quality < 30 and duplicated reads
            # rmdup:    remove PCR dupl_lts
            # faidx:    index and extract .fasta
            if not self.fail:
                start2 = time()
                status.running(5, 'SAMTools (quality control)', True)
    
                cmd_format = '{0} view -q' + str(cfg.user.ns_qcut) + ' -b {1} |{2} rmdup - {3}'
                cmd_args = [cfg.user.ne_samtools, strain_dir + strain + '_s.bam',
                            cfg.user.ne_samtools, strain_dir + strain + '_sfr.bam']
                worked = self._command_line(cmd_format, cmd_args)
                status.app.update()
            
                if worked:
                    cmd_format = '{0} index {1}'
                    cmd_args = [cfg.user.ne_samtools, strain_dir + strain + '_sfr.bam']
                    worked = self._command_line(cmd_format, cmd_args)

                # update the GUI
                if worked:
                    status.time(5, 'SAMTools (quality control)', time() - start2, True)
                else:
                    status.fail(5, 'SAMTools (quality control)', True)
            else:
                status.skip(5, 'SAMTools (quality control)', True)

            # generate Samtools mpileup file
            if not self.fail:
                start2 = time()
                status.running(6, 'SAMTools (generating mpileup)', True)
    
                cmd_format = '{0} mpileup -f {1} {2} > {3}'
                cmd_args = [cfg.user.ne_samtools, ref_path, strain_dir + strain + '_sfr.bam',
                            strain_dir + strain + '_sfr.mpileup']
                worked = self._command_line(cmd_format, cmd_args)
                
                if worked:
                    status.time(6, 'SAMTools (generating mpileup)', time() - start2, True)
                else:
                    status.fail(6, 'SAMTools (generating mpileup)', True)
            else:
                status.skip(6, 'SAMTools (generating mpileup)', True)

            """ VARIANT CALLING with VARSCAN """
            # SNP calling
            if not self.fail:
                start2 = time()
                status.running(7, 'VarScan (SNP + indel calling)', True)
                # compute the SNP .vcf file from the mpileup
                cmd_format = '{0} mpileup2snp {1} --p-value ' + str(cfg.user.ns_pv) + ' --output-vcf 1 > {2}'

                cmd_args = [cfg.user.ne_varscan, strain_dir + strain + '_sfr.mpileup', 
                            vcf_dir + strain + '_mpileup2snp.vcf']
                worked = self._command_line(cmd_format, cmd_args)
                status.app.update()

                # indel calling
                # compute the indel .vcf file from the mpileup
                if worked:
                    cmd_format = '{0} mpileup2indel {1} --p-value ' + str(cfg.user.ns_pv) + ' --output-vcf 1 > {2}'
                    cmd_args = [cfg.user.ne_varscan, strain_dir + strain + '_sfr.mpileup',
                                vcf_dir + strain + '_mpileup2indel.vcf']
                    worked = self._command_line(cmd_format, cmd_args)
                    status.app.update()

                # update the GUI
                if worked:
                    status.time(7, 'VarScan (SNP + indel calling)', time() - start2, True)
                else:
                    status.fail(7, 'VarScan (SNP + indel calling)', True)
            else:
                status.skip(7, 'VarScan (SNP + indel calling)', True)

            # if specified, remove unnecessary intermediate file
            if cfg.run.n_clean:
                start2 = time()
                status.running(8, 'cleaning up intermediate results', True)

                worked = clean_directory(strain, strain_dir, False)

                # update the GUI
                if worked:
                    status.time(8, 'cleaning up intermediate results', time() - start2, True)
                else:
                    status.fail(8, 'cleaning up intermediate results', True)
            else:
                status.skip(8, 'cleaning up intermediate results', True)

            some_success = some_success or not self.fail
            self.fail = False

        # merge the output files into a single mutation .tsv file that can be used for analysis
        if some_success:
            start2 = time()
            status.running(9, 'generating mutation .tsv')
            merger = MutationVcfMerger(cfg.user.no_dir, cfg.op.n_mut)
    
            if merger.empty or merger.write_error:
                status.fail(9, 'generating mutation .tsv')
            else:
                cfg.n_log.parser_errors = merger.parser_errors
                cfg.n_log.validate_errors = merger.validate_errors
                cfg.n_log.io_errors = merger.io_errors
                cfg.n_log.field_errors = merger.field_errors
                status.time(9, 'generating mutation .tsv', time() - start2)
        else:
            status.skip(9, 'generating mutation .tsv')

        # write the log
        start2 = time()
        status.running(10, 'writing log')
        worked = cfg.n_log.write(cfg.op.n_log)
        if worked:
            status.time(10, 'writing log', time() - start2)
        else:
            status.fail(10, 'writing log')

        # set fail status
        self.fail = not some_success
