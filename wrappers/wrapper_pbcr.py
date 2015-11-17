#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(SCRIPT_PATH + '/../src/')

import subprocess
import multiprocessing

import basicdefines
from dataspec import *

ASSEMBLER_URL = 'http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.3/wgs-8.3rc2-Linux_amd64.tar.bz2'
ASSEMBLER_PATH = os.path.join(basicdefines.ASSEMBLERS_PATH_ROOT_ABS, 'wgs-8.3rc2')
ZIP_FILE = 'wgs-8.3rc2-Linux_amd64.tar.bz2'
ZIP_PATH = os.path.join(basicdefines.ASSEMBLERS_PATH_ROOT_ABS, ZIP_FILE)
ASSEMBLER_BIN = os.path.join(ASSEMBLER_PATH,'Linux-amd64/bin/runCA')
ASSEMBLER_ECBIN = os.path.join(ASSEMBLER_PATH,'Linux-amd64/bin/PBcR')
ASSEMBLER_NAME = 'PBcR'
ASSEMBLER_RESULTS = 'out/9-terminator/asm.ctg.fasta'
CREATE_OUTPUT_FOLDER = True

PACBIO_SPEC = os.path.join(ASSEMBLER_PATH, 'wgs_pacbio.spec')
OXFORD_SPEC = os.path.join(ASSEMBLER_PATH, 'wgs_oxford.spec')

DO_ERROR_CORRECTION = True

# Run parameters
P_LENGTH = 500
P_PARTITIONS = 200

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
        # sys.stderr.write('Warning: dry_run == True\n');
#        sys.stderr.write('[%s wrapper] Executing (dryrun): "%s"\n' % (ASSEMBLER_NAME, command));
        log('Executing (dryrun): "%s".' % (command), fp_log);
    if (dry_run == False):
        # sys.stderr.write('[%s wrapper] Executing: "%s"\n' % (ASSEMBLER_NAME, command));
        log('Executing: "%s".' % (command), fp_log);
        subprocess.call(command, shell=True);
    # sys.stderr.write('\n');
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
# def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):
def run(datasets, output_path):
    ##################################################################################
    ### Sanity check for input datasets.
    ##################################################################################    
    machine_name = None
    reads_file = None;
    for dataset in datasets:
        if (machine_name != None and dataset.type != machine_name):
            sys.stderr.write('ERROR: %s is not a hybrid assembler, but datasets from disparate technologies are specified! Exiting.\n');
            exit(1);
        machine_name = dataset.type;
        reads_file = dataset.reads_path;
    if (machine_name == None):
        sys.stderr.write('ERROR: Input datasets not specified correctly! Exiting.\n');
        exit(1);

    if (len(datasets) > 1):
        sys.stderr.write('ERROR: More than one input dataset specified. Only one is expected. Exiting.\n');
        exit(1);

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
    ### Start the important work.
    ##################################################################################
    used_bin = ''
    if DO_ERROR_CORRECTION:
        used_bin = ASSEMBLER_ECBIN
    else:
        used_bin = ASSEMBLER_BIN
    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '.memtime')

    spec_file = ''
    if machine_name == 'pacbio':
        spec_file = PACBIO_SPEC
        command = 'cd %s; %s %s -pbCNS -length 500 -partitions 200 -l %s -s %s -fastq %s ' % (output_path, measure_command(memtime_path), used_bin, ASSEMBLER_NAME, spec_file, reads_file)
        execute_command(command, fp_log, dry_run=DRY_RUN);
    elif machine_name == 'nanopore':
        spec_file = OXFORD_SPEC
        command = 'cd %s; %s %s -length 500 -partitions 200 -l %s -s %s -fastq %s ' % (output_path, easure_command(memtime_path), used_bin, ASSEMBLER_NAME, spec_file, reads_file)
        execute_command(command, fp_log, dry_run=DRY_RUN);
    elif machine_name == 'illumina':
        log('\nMachine name "%s" not implemented for %s.\n' % (machine_name, ASSEMBLER_NAME));
        log('Skipping ....\n', fp_log)
        return;
    else:
        log('\nMachine name "%s" not implemented for %s.\n' % (machine_name, ASSEMBLER_NAME));
        log('Skipping ....\n', fp_log)
        return;

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
        sys.stderr.write('[%s wrapper] Bin found at %s. Skipping installation ...\n' % (ASSEMBLER_NAME, ASSEMBLER_BIN))
    else:
        sys.stderr.write('[%s wrapper] Started installation of %s.\n' % (ASSEMBLER_NAME, ASSEMBLER_NAME))

        if not os.path.exists(ZIP_PATH):
            sys.stderr.write('[%s wrapper] Downloading tar.bz2...\n' % (ASSEMBLER_NAME))
            command = 'cd %s; wget %s' % (basicdefines.ASSEMBLERS_PATH_ROOT_ABS, ASSEMBLER_URL)
            sys.stderr.write('[%s wrapper] %s\n' % (ASSEMBLER_NAME, command))
            subprocess.call(command, shell='True')

        # Decompress
        command = 'cd %s; tar -xvjf %s' % (basicdefines.ASSEMBLERS_PATH_ROOT_ABS, ZIP_FILE)
        sys.stderr.write('[%s wrapper] %s\n' % (ASSEMBLER_NAME, command))
        subprocess.call(command, shell='True')

        # make
        # WGS/Celera comes as precompiled binaries
        # Doesn't require make


# def verbose_usage_and_exit():
#     sys.stderr.write('Usage:\n')
#     sys.stderr.write('\t%s mode [<reads_file> <reference_file> <machine_name> <output_path> [<output_suffix>]]\n' % sys.argv[0])
#     sys.stderr.write('\n')
#     sys.stderr.write('\t- mode - either "run" or "install". Is "install" other parameters can be ommitted.\n')

#     exit(0)

# if __name__ == "__main__":
#     if (len(sys.argv) < 2 or len(sys.argv) > 7):
#         verbose_usage_and_exit()

#     if (sys.argv[1] == 'install'):
#         download_and_install()
#         exit(0)

#     elif (sys.argv[1] == 'run'):
#         if (len(sys.argv) < 6):
#             verbose_usage_and_exit()

#         reads_file = sys.argv[2]
#         reference_file = sys.argv[3]
#         machine_name = sys.argv[4]
#         output_path = sys.argv[5]
#         output_suffix = ''

#         if (len(sys.argv) == 7):
#             output_suffix = sys.argv[6]
#         run(reads_file, reference_file, machine_name, output_path, output_suffix)

#     else:
#         verbose_usage_and_exit()

def verbose_usage_and_exit():
    sys.stderr.write('Usage:\n')
    # sys.stderr.write('\t%s mode [<reads_file1>,<reads_file2>,...,<reads_fileN> <machine_name> <output_path>]\n' % sys.argv[0])
    sys.stderr.write('\t%s mode [<output_path> dataset1 [dataset2 ...]]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\t- mode - either "run" or "install". If "install" other parameters can be omitted.\n')
    sys.stderr.write('\t- dataset - specification of a dataset in the form: reads_type,<reads_path>[<reads_path_b,frag_len,frag_stddev] .\n');
    sys.stderr.write('\t            Reads_type can be nanopore/pacbio/single/paired/mate. If reads_type == "paired" or "mate", last three parameters can be omitted".\n');
    sys.stderr.write('\t            If reads_type == "paired" or "mate", other end of the pair needs to be in another file provided by reads_path_b.\n');
    sys.stderr.write('\n');

    sys.stderr.write('Example:\n');
    sys.stderr.write('\t%s run results/%s nanopore,datasets/reads.fastq\n' % (os.path.basename(sys.argv[0]), ASSEMBLER_NAME));
    sys.stderr.write('\n');

    exit(0)

if __name__ == "__main__":
    if (len(sys.argv) < 2 or len(sys.argv) > 7):
        verbose_usage_and_exit()

    if (sys.argv[1] == 'install'):
        download_and_install()
        exit(0)

    elif (sys.argv[1] == 'run'):
        if (len(sys.argv) < 4):
            verbose_usage_and_exit()

        output_path = sys.argv[2]
        datasets = [];
        for arg in sys.argv[3:]:
            dataset = Dataset(arg);
            datasets.append(dataset);
        run(datasets, output_path);

    else:
        verbose_usage_and_exit()