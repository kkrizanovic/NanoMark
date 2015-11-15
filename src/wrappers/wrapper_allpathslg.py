#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(SCRIPT_PATH + '/../')

import re;

import subprocess
import multiprocessing
import getpass
from time import gmtime, strftime

try:
    import basicdefines
    MODULE_BASICDEFINES = True;
    ASSEMBLERS_PATH_ROOT_ABS = basicdefines.ASSEMBLERS_PATH_ROOT_ABS;
    TOOLS_ROOT = basicdefines.TOOLS_ROOT;
except:
    MODULE_BASICDEFINES = False;
    ASSEMBLERS_PATH_ROOT_ABS = os.path.join(SCRIPT_PATH, 'assemblers/');
    TOOLS_ROOT = '%s' % (SCRIPT_PATH);

ASSEMBLER_URL = 'http://www.cellprofiler.org/ftp/pub/crd/nightly/allpathslg/allpathslg-52488.tar.gz'
ASSEMBLER_PATH = os.path.join(ASSEMBLERS_PATH_ROOT_ABS, 'allpathslg-52488')
ASSEMBLER_BIN = os.path.join(ASSEMBLER_PATH, 'src/allpaths-lg')
ASSEMBLER_NAME = 'ALLPATHS-LG'
# ASSEMBLER_RESULTS = 'contig-100.fa'
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

def execute_command(command, fp_log, dry_run=True):
    if (dry_run == True):
        log('Executing (dryrun): "%s".' % (command), fp_log);
    if (dry_run == False):
        log('Executing: "%s".' % (command), fp_log);
        subprocess.call(command, shell=True);
    log('\n', fp_log);

def measure_command(measure_file):
    if (MODULE_BASICDEFINES == True):
        return basicdefines.measure_command(measure_file);
    else:
        sys.stderr.write('ERROR: Cgmemtime tool not found! Exiting.\n');
        exit(1);
        # return '/usr/bin/time --format "Command line: %%C\\nReal time: %%e s\\nCPU time: -1.0 s\\nUser time: %%U s\\nSystem time: %%S s\\nMaximum RSS: %%M kB\\nExit status: %%x" --quiet -o %s ' % measure_file;

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

# Function 'run' should provide a standard interface for running a mapper. Given input parameters, it should run the
# alignment process, and convert any custom output results to the SAM format. Function should return a string with the
# path to the output file.
#    reads_file            Path to a FASTA/FASTQ file containing reads.
#    reference_file        Path to a reference genome FASTA file.
#    machine_name        A symbolic name to specify a set of parameters for a specific sequencing platform.
#    output_path            Folder to which the output will be placed to. Filename will be automatically generated according to the name of the mapper being run.
#    output_suffix        A custom suffix that can be added to the output filename.
def run(reads_files, reference_file, machine_name, output_path, output_suffix=''):
    ### Reference for running a local job instead of using an SGE cluster:
    ### http://seqanswers.com/forums/showthread.php?t=50937
    num_threads = multiprocessing.cpu_count() / 2;

    reference_file = os.path.abspath(reference_file);
    output_path = os.path.abspath(output_path);

    ##################################################################################
    ### Backup old assembly results, and create the new output folder.
    ##################################################################################
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
        execute_command('%s/cgmemtime/cgmemtime -o %s/test.memtime ls -lhrt' % (basicdefines.TOOLS_ROOT, output_path), fp_log, dry_run=DRY_RUN);
        if (not os.path.exists('%s/test.memtime' % (output_path))):
            command = 'sudo %s/cgmemtime/cgmemtime --setup -g %s --perm 775' % (basicdefines.TOOLS_ROOT, getpass.getuser());
            sys.stderr.write('[] %s\n' % (command));
            execute_command(command, fp_log, dry_run=DRY_RUN);

    ##################################################################################
    ### Start the important work.
    ##################################################################################
    log('Running assembly using %s.' % (ASSEMBLER_NAME), fp_log);

    ### Create a config file for the assembly.
    if (machine_name == 'pacbio'):
        pass;
    elif (machine_name == 'nanopore'):
# ## Prepare data for allpaths, dataset1
# # Prepare nanopore reads file
# # Edit in_groups.csv and in_libs.csv files, they should contain one entery for a one fastq file
# # create folder structure
# cd; cd benchmark/results/dataset1_20x2d/allpaths
# mkdir benchmark
# cd benchmark; mkdir data
# # Run PrepareAllPathsInputs.pl
# cd; cd benchmark/results/datasets/dataset1_20x2d/allpaths

# PrepareAllPathsInputs.pl DATA_DIR=$PWD/benchmark/data_test PLOIDY=1 IN_GROUPS_CSV=in_groups.csv IN_LIBS_CSV=in_libs.csv GENOME_SIZE=4650000 PICARD_TOOLS_DIR=/home/kresimir/build/picard-tools-1.140 OVERWRITE=True | tee prepare.out

# # Prepare data with specific coverage
# PrepareAllPathsInputs.pl DATA_DIR=$PWD/benchmark/data_20_30 PLOIDY=1 FRAG_COVERAGE=20 JUMP_COVERAGE=30 IN_GROUPS_CSV=in_groups.csv IN_LIBS_CSV=in_libs.csv GENOME_SIZE=4650000 PICARD_TOOLS_DIR=/home/kresimir/build/picard-tools-1.140 OVERWRITE=True | tee prepare_20_30.out

# # Copy created files to data folder, copy jump and frag libraries to data folder
# # Rename created files to "long_reads_orig" .fastb and .qualb

# # Run allpaths on dataset1
# RunAllPathsLG PRE=$PWD REFERENCE_NAME=benchmark DATA_SUBDIR=data_test RUN=run TARGETS=standard  OVERWRITE=True | tee -a assemble_test.out

# ./dataset4_2_40x2d/Allpaths-LG/in_groups.csv
# ./dataset4_2_40x2d/Allpaths-LG/in_libs.csv
#                                                      file_name, library_name, group_name
# /home/kresimir/benchmark/datasets/dataset4_2_40x2d/reads.fastq,     dataset4,   dataset4
         pass;

    memtime_file = '%s/%s.memtime' % (output_path, ASSEMBLER_NAME);
    run_commands = [];
    command = '; '.join(run_commands);
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

def verbose_usage_and_exit():
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s mode [<reads_file1>,<reads_file2>,...,<reads_fileN> <machine_name> <output_path>]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\t- mode - either "run" or "install". Is "install" other parameters can be ommitted.\n')
    # sys.stderr.write('\n');
    # sys.stderr.write('\tMultiple reads file can')

    exit(0)

if __name__ == "__main__":
    if (len(sys.argv) < 2 or len(sys.argv) > 7):
        verbose_usage_and_exit()

    if (sys.argv[1] == 'install'):
        download_and_install()
        exit(0)

    elif (sys.argv[1] == 'run'):
        if (len(sys.argv) < 5 or len(sys.argv) > 6):
            verbose_usage_and_exit()

        reads_files = sys.argv[2].split(',')         ### Enable specifying multiple FASTQ files for input.
        machine_name = sys.argv[3]
        output_path = sys.argv[4]
        reference_file = '';
        output_suffix = ''
        
        run(reads_files, reference_file, machine_name, output_path, output_suffix)

    else:
        verbose_usage_and_exit()
