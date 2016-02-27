#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(SCRIPT_PATH + '/../src/')

import subprocess
import multiprocessing
from time import gmtime, strftime

from dataspec import *

ASSEMBLER_TYPE = 'nonhybrid';   # hybrid or nonhybrid

try:
    import basicdefines
    MODULE_BASICDEFINES = True;
    ASSEMBLERS_PATH_ROOT_ABS = basicdefines.ASSEMBLERS_PATH_ROOT_ABS;
    TOOLS_ROOT = basicdefines.TOOLS_ROOT;
except:
    MODULE_BASICDEFINES = False;
    ASSEMBLERS_PATH_ROOT_ABS = os.path.join(SCRIPT_PATH, 'assemblers/');
    TOOLS_ROOT = '%s' % (SCRIPT_PATH);

ASSEMBLER_URL = 'https://github.com/lh3/miniasm.git'
ASSEMBLER_PATH = os.path.join(ASSEMBLERS_PATH_ROOT_ABS, 'miniasm')
ASSEMBLER_BIN = os.path.join(ASSEMBLER_PATH,'miniasm/miniasm')
ASSEMBLER_NAME = 'Miniasm'
# ASSEMBLER_RESULTS = 'out/9-terminator/asm.ctg.fasta'
ASSEMBLY_UNPOLISHED = 'benchmark-unpolished_assembly.fasta'
ASSEMBLY_POLISHED = 'benchmark-polished_assembly.fasta'
CREATE_OUTPUT_FOLDER = True

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
def run(datasets, output_path, approx_genome_len=0, move_exiting_out_path=True):
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

    # if (len(datasets) > 1):
    #     sys.stderr.write('ERROR: More than one input dataset specified. Only one is expected. Exiting.\n');
    #     exit(1);

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
    ### Preparing the input datasets.
    ##################################################################################
    reads_file = '%s/all_reads.fastq' % (output_path);
    for dataset in datasets:
        if (dataset.reads_path.endswith('fasta') or dataset.reads_path.endswith('fa')):
            log('ERROR: Assembler %s expects only FASTQ files for input. Trouble loading "%s". Exiting.\n' % (ASSEMBLER_NAME, dataset.reads_path));
        else:
            log('Aggregating FASTQ file "%s" to "%s".\n' % (dataset.reads_path, reads_file), fp_log);
            command = 'cat %s >> %s' % (dataset.reads_path, reads_file);
            execute_command(command, fp_log, dry_run=DRY_RUN);

    ##################################################################################
    ### Start the important work.
    ##################################################################################
    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '.memtime')

    spec_file = ''
    if machine_name == 'pacbio':
        overlaps_file = '%s/overlaps.paf.gz' % (output_path);
        command = '%s/minimap/minimap -Sw5 -L100 -m0 -t%d %s %s | gzip -1 > %s' % (ASSEMBLER_PATH, num_threads, reads_file, reads_file, overlaps_file);
        execute_command(command, fp_log, dry_run=DRY_RUN);

        gfa_file = '%s/miniasm.gfa' % (output_path);
        command = '%s -f %s %s > %s' % (ASSEMBLER_BIN, reads_file, overlaps_file, gfa_file);
        execute_command(command, fp_log, dry_run=DRY_RUN);

        command = 'awk \'$1 ~/S/ {print ">"$2"\\n"$3}\' %s > %s/%s' % (gfa_file, output_path, ASSEMBLY_UNPOLISHED);
        execute_command(command, fp_log, dry_run=DRY_RUN);

    elif machine_name == 'nanopore':
        overlaps_file = '%s/overlaps.paf.gz' % (output_path);
        command = '%s/minimap/minimap -Sw5 -L100 -m0 -t%d %s %s | gzip -1 > %s' % (ASSEMBLER_PATH, num_threads, reads_file, reads_file, overlaps_file);
        execute_command(command, fp_log, dry_run=DRY_RUN);

        gfa_file = '%s/miniasm.gfa' % (output_path);
        command = '%s -f %s %s > %s' % (ASSEMBLER_BIN, reads_file, overlaps_file, gfa_file);
        execute_command(command, fp_log, dry_run=DRY_RUN);

        command = 'awk \'$1 ~/S/ {print ">"$2"\\n"$3}\' %s > %s/%s' % (gfa_file, output_path, ASSEMBLY_UNPOLISHED);
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

        if (not os.path.exists(ASSEMBLER_PATH)):
            log('Creating a directory on path "%s".' % (ASSEMBLER_PATH), None);
            os.makedirs(ASSEMBLER_PATH);

        command = 'cd %s; git clone https://github.com/lh3/minimap && (cd minimap && make)' % (ASSEMBLER_PATH)
        execute_command(command, None, dry_run=DRY_RUN);

        command = 'cd %s; git clone %s && (cd miniasm && make)' % (ASSEMBLER_PATH, ASSEMBLER_URL)
        execute_command(command, None, dry_run=DRY_RUN);

def verbose_usage_and_exit():
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s mode [<output_path> approx_genome_len dataset1 [dataset2 ...]]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\t- mode - either "run" or "install". If "install" other parameters can be omitted.\n')
    sys.stderr.write('\t- dataset - specification of a dataset in the form: reads_type,<reads_path>[<reads_path_b,frag_len,frag_stddev] .\n');
    sys.stderr.write('\t            Reads_type can be nanopore/pacbio/single/paired/mate. If reads_type != "paired" or "mate", last three parameters can be omitted".\n');
    sys.stderr.write('\t            If reads_type == "paired" or "mate", other end of the pair needs to be in another file provided by reads_path_b.\n');
    sys.stderr.write('\t- approx_genome_len - approximate length of the genome to be assembled. This parameter is neccesary for PBcR.\n');
    sys.stderr.write('\n');

    sys.stderr.write('Example:\n');
    sys.stderr.write('\t%s run results/%s - nanopore,datasets/reads.fastq\n' % (os.path.basename(sys.argv[0]), ASSEMBLER_NAME));
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
