#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(SCRIPT_PATH + '/../src/')

import re;

import subprocess
import multiprocessing
import getpass
from time import gmtime, strftime
from dataspec import *

ASSEMBLER_TYPE = 'hybrid';   # hybrid or nonhybrid

try:
    import basicdefines
    MODULE_BASICDEFINES = True;
    ASSEMBLERS_PATH_ROOT_ABS = basicdefines.ASSEMBLERS_PATH_ROOT_ABS;
    TOOLS_ROOT = basicdefines.TOOLS_ROOT;
    CGMEMTIME_PATH = basicdefines.CGMEMTIME_PATH;
    CGMEMTIME_FILE = basicdefines.CGMEMTIME_FILE;
    TOOLS_ROOT_ABS = basicdefines.TOOLS_ROOT_ABS;
except:
    MODULE_BASICDEFINES = False;
    ASSEMBLERS_PATH_ROOT_ABS = os.path.join(SCRIPT_PATH, 'assemblers/');
    TOOLS_ROOT = '%s' % (SCRIPT_PATH);
    CGMEMTIME_PATH = os.path.join(SCRIPT_PATH, 'tools/cgmemtime/');
    CGMEMTIME_FILE = CGMEMTIME_PATH + '/cgmemtime';
    TOOLS_ROOT_ABS = '%s/tools/' % (SCRIPT_PATH);

ASSEMBLER_URL = 'ftp://ftp.broadinstitute.org/pub/crd/ALLPATHS/Release-LG/latest_source_code/allpathslg-52488.tar.gz'
ASSEMBLER_PATH = os.path.join(ASSEMBLERS_PATH_ROOT_ABS, 'allpathslg-52488')
ASSEMBLER_BIN = os.path.join(ASSEMBLER_PATH, 'bin/RunAllPathsLG')
ASSEMBLER_NAME = 'ALLPATHS-LG'
# ASSEMBLER_RESULTS = 'contig-100.fa'
ASSEMBLY_UNPOLISHED = 'benchmark-unpolished_assembly.fasta'
ASSEMBLY_POLISHED = 'benchmark-polished_assembly.fasta'
CREATE_OUTPUT_FOLDER = True

PICARDTOOLS_URL = 'https://github.com/broadinstitute/picard/releases/download/1.140/picard-tools-1.140.zip';
PICARDTOOLS_PATH = os.path.join(ASSEMBLER_PATH, 'picard-tools-1.140');

# DRY_RUN = True;
DRY_RUN = False;

### Logs messages to STDERR and an output log file if provided (opened elsewhere).
def log(message, fp_log):
    timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());

    sys.stderr.write('[%s wrapper %s] %s\n' % (ASSEMBLER_NAME, timestamp, message))
    sys.stderr.flush();
    if (fp_log != None):
        fp_log.write('[%s wrapper %s] %s\n' % (ASSEMBLER_NAME, timestamp, message))
        fp_log.flush();

import traceback;
def execute_command(command, fp_log, dry_run=True):
    if (dry_run == True):
        log('Executing (dryrun): "%s".' % (command), fp_log);
        log('\n', fp_log);
        return 0;
    if (dry_run == False):
        log('Executing: "%s".' % (command), fp_log);
        rc = subprocess.call(command, shell=True);
        if (rc != 0):
            log('ERROR: subprocess call returned error code: %d.' % (rc), fp_log);
            log('Traceback:', fp_log);
            traceback.print_stack(fp_log);
            exit(1);
        return rc;

# def execute_command(command, fp_log, dry_run=True):
#     sys.stderr.write('Executing command: "%s"\n' % command);
#     if (dry_run == True):
#         log('Executing (dryrun): "%s".' % (command), fp_log);
#     if (dry_run == False):
#         log('Executing: "%s".' % (command), fp_log);
#         p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
#     [output, err] = p.communicate()
#     rc = p.returncode
#     sys.stderr.write('\n');
#     if (rc != 0):
#         log('ERROR: subprocess call returned error code: %d.\n' % (rc), fp_log);
#         traceback.print_stack(fp_log);
#         exit(1);
#     return [rc, output, err];

def measure_command(measure_file):
    if (MODULE_BASICDEFINES == True and os.path.exists(CGMEMTIME_FILE)):
        return basicdefines.measure_command(measure_file);
    else:
        return '/usr/bin/time --format "Command line: %%C\\nReal time: %%e s\\nCPU time: -1.0 s\\nUser time: %%U s\\nSystem time: %%S s\\nMaximum RSS: %%M kB\\nExit status: %%x" --quiet -o %s ' % measure_file;

def peek(fp, num_chars):
    data = fp.read(num_chars);
    if len(data) == 0:
        return '';
    fp.seek(num_chars * -1, 1);
    return data;



# Returns a single read from the given FASTA/FASTQ file.
# Parameter header contains only the header of the read.
# Parameter lines contains all lines of the read, which include:
# - header
# - seq
# - '+' if FASTQ
# - quals if FASTQ
# Parameter lines is an array of strings, each for one component.
# Please note that multiline FASTA/FASTQ entries (e.g. sequence line)
# will be truncated into one single line.
def get_single_read(fp):
    lines = [];
    
    line = fp.readline();
    header = line.rstrip();
    header_leading_char = '';
    if (len(header) > 0):
        sequence_separator = header[0];
        header_leading_char = header[0];
        header = header[1:];            # Strip the '>' or '@' sign from the beginning.
    else:
        return ['', []];
    
    next_char = peek(fp, 1);
    
    line_string = '';
    lines.append(header_leading_char + header);
    
    num_lines = 1;
    #while len(next_char) > 0 and next_char != sequence_separator or (next_char == '@' and num_lines < 4):
    while (len(next_char) > 0 and (next_char != sequence_separator or (next_char == '@' and num_lines < 4))):
        line = fp.readline();
        if (line.rstrip() == '+' or line.rstrip() == ('+' + header)):
        #if (line.rstrip()[0] == '+'):
            lines.append(line_string);
            lines.append(line.rstrip());
            line_string = '';
        else:
            line_string += line.rstrip();
        next_char = peek(fp, 1);
        num_lines += 1;
        
    lines.append(line_string);
    
    return [header, lines];

##############################################################
##############################################################
##############################################################
def parse_memtime(memtime_path):
    cmdline = '';
    realtime = 0;
    cputime = 0;
    usertime = 0;
    systemtime = 0;
    maxrss = 0;
    rsscache = 0;
    time_unit = '';
    mem_unit = '';

    try:
        fp = open(memtime_path, 'r');
        lines = [line.strip() for line in fp.readlines() if (len(line.strip()) > 0)];
        fp.close();
    except Exception, e:
        sys.stderr.write('Could not find memory and time statistics in file "%s".\n' % (memtime_path));
        return [cmdline, realtime, cputime, usertime, systemtime, maxrss, time_unit, mem_unit];

    for line in lines:
        if (line.startswith('Command line:')):
            cmdline = line.split(':')[1].strip();
        elif (line.startswith('Real time:')):
            split_line = line.split(':')[1].strip().split(' ');
            realtime = float(split_line[0].strip());
            time_unit = split_line[1].strip();
        elif (line.startswith('CPU time:')):
            split_line = line.split(':')[1].strip().split(' ');
            cputime = float(split_line[0].strip());
            time_unit = split_line[1].strip();
        elif (line.startswith('User time:')):
            split_line = line.split(':')[1].strip().split(' ');
            usertime = float(split_line[0].strip());
            time_unit = split_line[1].strip();
        elif (line.startswith('System time:')):
            split_line = line.split(':')[1].strip().split(' ');
            systemtime = float(split_line[0].strip());
            time_unit = split_line[1].strip();
        elif (line.startswith('Maximum RSS:')):
            split_line = line.split(':')[1].strip().split(' ');
            maxrss = float(split_line[0].strip());
            mem_unit = split_line[1].strip();
        # elif (line.startswith('')):
        #   split_line = line.split(':')[1].strip().split(' ');
        #   rsscache = float(split_line[0].strip());
        #   mem_unit = split_line[1].strip();

    return [cmdline, realtime, cputime, usertime, systemtime, maxrss, time_unit, mem_unit];

def parse_memtime_files_and_accumulate(memtime_files, final_memtime_file):
    final_command_line = '';
    final_real_time = 0.0;
    final_cpu_time = 0.0;
    final_user_time = 0.0;
    final_system_time = 0.0;
    final_time_unit = '';
    final_max_rss = 0;
    final_mem_unit = '';

    i = 0;
    for memtime_file in memtime_files:
        i += 1;
        sys.stderr.write('Parsing memtime file "%s"...\n' % (memtime_file));

        [cmdline, realtime, cputime, usertime, systemtime, maxrss, time_unit, mem_unit] = parse_memtime(memtime_file);
        if (i == 1):
            final_command_line = cmdline;
            final_real_time = realtime;
            final_cpu_time = cputime;
            final_user_time = usertime;
            final_system_time = systemtime;
            final_max_rss = maxrss;
            final_time_unit = time_unit;
            final_mem_unit = mem_unit;
        else:
            if (time_unit == final_time_unit and mem_unit == final_mem_unit):
                final_command_line += '; ' + cmdline;
                final_real_time += realtime;
                final_cpu_time += cputime;
                final_user_time += usertime;
                final_system_time += systemtime;
                final_max_rss = max(final_max_rss, maxrss);
            else:
                sys.stderr.write('Memory or time units not the same in all files! Instead of handling this, we decided to be lazy and just give up.\n');
                break;

    try:
        fp = open(final_memtime_file, 'w');
    except Exception, e:
        sys.stderr.write('ERROR: Could not open file "%s" for writing!\n' % (final_memtime_file));
        return;

    if (final_cpu_time <= 0.0):
        final_cpu_time = final_user_time + final_system_time;

    fp.write('Command line: %s\n' % (final_command_line));
    fp.write('Real time: %f %s\n' % (final_real_time, final_time_unit));
    fp.write('CPU time: %f %s\n' % (final_cpu_time, final_time_unit));
    fp.write('User time: %f %s\n' % (final_user_time, final_time_unit));
    fp.write('System time: %f %s\n' % (final_system_time, final_time_unit));
    fp.write('Maximum RSS: %f %s\n' % (final_max_rss, final_mem_unit));

    fp.close();

##############################################################
##############################################################
##############################################################

# Function 'run' should provide a standard interface for running a mapper. Given input parameters, it should run the
# alignment process, and convert any custom output results to the SAM format. Function should return a string with the
# path to the output file.
#    reads_file            Path to a FASTA/FASTQ file containing reads.
#    reference_file        Path to a reference genome FASTA file.
#    machine_name        A symbolic name to specify a set of parameters for a specific sequencing platform.
#    output_path            Folder to which the output will be placed to. Filename will be automatically generated according to the name of the mapper being run.
#    output_suffix        A custom suffix that can be added to the output filename.
# def run(reads_files, reference_file, machine_name, output_path, output_suffix=''):
def run(datasets, output_path, approx_genome_len=0, move_exiting_out_path=True):
    ##################################################################################
    ### Simple variable definitions.
    ##################################################################################        
    ### Reference for running a local job instead of using an SGE cluster:
    ### http://seqanswers.com/forums/showthread.php?t=50937
    num_threads = multiprocessing.cpu_count() / 2;
    output_path = os.path.abspath(output_path);

    ##################################################################################
    ### Backup old assembly results, and create the new output folder.
    ##################################################################################
    if (move_exiting_out_path == True):
        if (os.path.exists(output_path)):
            timestamp = strftime("%Y_%m_%d-%H_%M_%S", gmtime());
            os.rename(output_path, '%s.bak_%s' % (output_path, timestamp));
    if (not os.path.exists(output_path)):
        log('Creating a directory on path "%s".' % (output_path), None);
        os.makedirs(output_path);

    ##################################################################################
    ### Prepare a log file.
    ##################################################################################
    log_file = '%s/wrapper_log.txt' % (output_path);
    try:
        fp_log = open(log_file, 'a');
    except Exception, e:
        log('ERROR: Could not open file "%s" for writing! Using only STDERR for logging.' % (log_file), None);
        fp_log = None;

    ##################################################################################
    ### Check if permissions are given for Cgmemtime.
    ##################################################################################
    if (MODULE_BASICDEFINES == True):
        execute_command('%s date' % (measure_command('%s/test.memtime' % output_path)), fp_log, dry_run=DRY_RUN);
        if (not os.path.exists('%s/test.memtime' % (output_path))):
            command = 'sudo %s/cgmemtime/cgmemtime --setup -g %s --perm 775' % (basicdefines.TOOLS_ROOT, getpass.getuser());
            sys.stderr.write('[] %s\n' % (command));
            execute_command(command, fp_log, dry_run=DRY_RUN);

    ##################################################################################
    ### Preparing the input datasets.
    ##################################################################################
    log('Preparing the input datasets.\n', fp_log);

    ### 'type' parameter is only informative (according to the ALLPATHS-LG documentation).
    in_libs = 'library_name, project_name, organism_name, type, paired, frag_size, frag_stddev, insert_size, insert_stddev, read_orientation, genomic_start, genomic_end\n';
    in_groups = 'group_name, library_name, file_name\n';
    num_libs = 0;
    for dataset in datasets:
        num_libs += 1;
        # library_name = '%s-%d' % (os.path.splitext(os.path.basename(dataset.reads_path))[0], num_libs);
        library_name = 'lib-%d' % (num_libs);
        group_name = 'group-%d' % (num_libs);

        in_groups += '%s, %s, %s\n' % (group_name, library_name, dataset.reads_path);

        if (dataset.type == 'single'):
            in_libs += '%s, dataset, dataset, fragment, 0, , , , , inward, 0, 0\n' % (library_name);

        elif (dataset.type == 'paired'):
            # in_groups += '%s-a, %s, %s\n' % (group_name, library_name, dataset.reads_path_a);
            # in_groups += '%s-b, %s, %s\n' % (group_name, library_name, dataset.reads_path_b);
            in_libs += '%s, dataset, dataset, fragment, 1, %d, %d, , , inward, 0, 0\n' % (library_name, dataset.frag_len, dataset.frag_stddev);

        elif (dataset.type == 'mate'):
            in_libs += '%s, dataset, dataset, jumping, 1, , , %d, %d, outward, 0, 0\n' % (library_name, dataset.frag_len, dataset.frag_stddev);

        elif (dataset.type == 'pacbio'):
            in_libs += '%s, dataset, dataset, pacbio, 0, , , , , inward, 0, 0\n' % (library_name);

        elif (dataset.type == 'nanopore'):
            in_libs += '%s, dataset, dataset, nanopore, 0, , , , , inward, 0, 0\n' % (library_name);

    log(in_groups, fp_log);
    log(in_libs, fp_log);

    in_groups_csv_path = '%s/in_groups.csv' % (output_path);
    in_libs_csv_path = '%s/in_libs.csv' % (output_path);

    fp_groups = open(in_groups_csv_path, 'w');
    fp_groups.write(in_groups);
    fp_groups.close();

    fp_libs = open(in_libs_csv_path, 'w');
    fp_libs.write(in_libs);
    fp_libs.close();

    ##################################################################################
    ### Start the important work.
    ##################################################################################
    log('Running assembly using %s.' % (ASSEMBLER_NAME), fp_log);

    num_memtimes = 0;

    data_dir = '%s/data/data_test' % (output_path);
    if (not os.path.exists(data_dir)):
        log('Creating a directory on path "%s".' % (data_dir), None);
        os.makedirs(data_dir);

    ### ALLPATHS-LG crashes when it cannot parse a double value from a string, e.g. "101.0" would cause a crash if the locale is not set to US. Set correct decimal separator!! export LC_NUMERIC='en_US.utf8'
    command = 'export LC_NUMERIC=\'en_US.utf8\'; export PATH=$PATH:%s/bin; %s %s/bin/PrepareAllPathsInputs.pl DATA_DIR=%s PLOIDY=1 IN_GROUPS_CSV=%s IN_LIBS_CSV=%s PICARD_TOOLS_DIR=%s OVERWRITE=True | tee %s' % \
                (ASSEMBLER_PATH, measure_command('%s/%s-%d.memtime' % (output_path, ASSEMBLER_NAME, num_memtimes)), ASSEMBLER_PATH, data_dir, in_groups_csv_path, in_libs_csv_path, PICARDTOOLS_PATH, 'prepare.out');
    execute_command(command, fp_log, dry_run=DRY_RUN);

    num_memtimes += 1;
    ### ALLPATHS-LG crashes when it cannot parse a double value from a string, e.g. "101.0" would cause a crash if the locale is not set to US. Set correct decimal separator!! export LC_NUMERIC='en_US.utf8'
    command = 'export LC_NUMERIC=\'en_US.utf8\'; export PATH=$PATH:%s/bin; %s %s PRE=%s DATA_SUBDIR=data_test REFERENCE_NAME=data THREADS=%d RUN=run TARGETS=standard OVERWRITE=True | tee -a %s/assemble_test.out' % \
                (ASSEMBLER_PATH, measure_command('%s/%s-%d.memtime' % (output_path, ASSEMBLER_NAME, num_memtimes)), ASSEMBLER_BIN, output_path, num_threads, output_path);
    execute_command(command, fp_log, dry_run=DRY_RUN);

    ##################################################################################
    ### Accumulate the memtime stats.
    ##################################################################################
    # memtime_file = '%s/%s.memtime' % (output_path, ASSEMBLER_NAME);
    memtime_file = '%s/total.memtime' % (output_path);
    all_memtimes = ['%s-%d.memtime' % (ASSEMBLER_NAME, value) for value in xrange(1, num_memtimes)];
    parse_memtime_files_and_accumulate(all_memtimes, memtime_file);

    command = 'cp %s/data/data_test/run/ASSEMBLIES/test/final.assembly.fasta %s/%s' % (output_path, output_path, ASSEMBLY_UNPOLISHED);
    execute_command(command, fp_log, dry_run=DRY_RUN);
    
    if (fp_log != None):
        fp_log.close();



# A placeholder for a function that runs quast on assembly results
# This is the same for all wrappers at the moment and is done in the main program
def run_quast():
    pass

# A function that gets the results of cgmemtime for a wrapper
# Since some aligners are run multiple times and create multiple measurements
# this is specific for each wrapper
# Returns relevant real time, user time and maximum memory reserved
def get_memtime():
    pass


# This is a standard interface for setting up an assembler. It should assume that the assembler
# is not present localy, but needs to be retrieved, unpacked, compiled and set-up, without requireing
# root privileges.
def download_and_install():
    if os.path.exists(ASSEMBLER_BIN):
        log('Bin found at %s. Skipping installation.' % (ASSEMBLER_BIN), None);

    else:
        log('Started installation of %s.' % (ASSEMBLER_NAME), None);
        if (not os.path.exists(ASSEMBLERS_PATH_ROOT_ABS)):
            log('Creating a directory on path "%s".' % (ASSEMBLERS_PATH_ROOT_ABS), None);
            os.makedirs(ASSEMBLERS_PATH_ROOT_ABS);

        command = 'cd %s; wget %s; tar -xzvf %s' % (ASSEMBLERS_PATH_ROOT_ABS, ASSEMBLER_URL, os.path.basename(ASSEMBLER_URL));
        execute_command(command, None, dry_run=DRY_RUN);

        command = 'cd %s; ./configure --prefix=%s' % (ASSEMBLER_PATH, ASSEMBLER_PATH);
        execute_command(command, None, dry_run=DRY_RUN);

        command = 'cd %s; make -j; make install' % (ASSEMBLER_PATH);
        execute_command(command, None, dry_run=DRY_RUN);

        command = 'cd %s; wget %s; unzip %s' % (ASSEMBLER_PATH, PICARDTOOLS_URL, os.path.basename(PICARDTOOLS_URL));
        execute_command(command, None, dry_run=DRY_RUN);

    # if os.path.exists(CGMEMTIME_PATH + '/' + CGMEMTIME_BIN):
    #     sys.stderr.write('Cgmemtime already installed. Skipping...\n')
    # else:
    #     command = 'mkdir -p %s; cd %s; git clone https://github.com/isovic/cgmemtime.git' % (TOOLS_ROOT_ABS, TOOLS_ROOT_ABS)
    #     execute_command(command, None, dry_run=DRY_RUN);

def verbose_usage_and_exit():
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s mode [<output_path> approx_genome_len dataset1 [dataset2 ...]]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\t- mode - either "run" or "install". If "install" other parameters can be omitted.\n')
    sys.stderr.write('\t- dataset - specification of a dataset in the form: reads_type,<reads_path>[<reads_path_b,frag_len,frag_stddev] .\n');
    sys.stderr.write('\t            Reads_type can be nanopore/pacbio/single/paired/mate. If reads_type != "paired" or "mate", last three parameters can be omitted".\n');
    sys.stderr.write('\t            If reads_type == "paired" or "mate", other end of the pair needs to be in another file provided by reads_path_b.\n');
    sys.stderr.write('\t- approx_genome_len - approximate length of the genome to be assembled. If unknown, use "-" or "0" instead.\n');
    sys.stderr.write('\n');

    sys.stderr.write('Example:\n');
    sys.stderr.write('\t%s run results/%s - paired,datasets/frag_reads.Solexa-25396.A.fastq,datasets/frag_reads.Solexa-25396.B.fastq,180,10 mate,datasets/jump_reads.Solexa-42866.A.fastq,datasets/jump_reads.Solexa-42866.B.fastq,3000,500 mate,datasets/jump_reads.Solexa-44956.A.fastq,datasets/jump_reads.Solexa-44956.B.fastq,3000,500 nanopore,datasets/reads.fastq\n' % (os.path.basename(sys.argv[0]), ASSEMBLER_NAME));
    sys.stderr.write('\n');

    exit(0)

if __name__ == "__main__":
    if (len(sys.argv) < 2):
        verbose_usage_and_exit()

    if (sys.argv[1] == 'install'):
        download_and_install()

    elif (sys.argv[1] == 'run'):
        if (len(sys.argv) < 5):
            verbose_usage_and_exit()

        output_path = sys.argv[2]
        approx_genome_len = 0 if (sys.argv[3] == '-') else int(sys.argv[3]);
        datasets = [];
        for arg in sys.argv[4:]:
            dataset = Dataset(arg);
            datasets.append(dataset);
        run(datasets, output_path, approx_genome_len);

    else:
        verbose_usage_and_exit()
