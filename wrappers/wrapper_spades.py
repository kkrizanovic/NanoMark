#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(SCRIPT_PATH + '/../src/')

import subprocess
import multiprocessing

from time import gmtime, strftime

from dataspec import *

ASSEMBLER_TYPE = 'hybrid';   # hybrid or nonhybrid

try:
    import basicdefines
    MODULE_BASICDEFINES = True;
    ASSEMBLERS_PATH_ROOT_ABS = basicdefines.ASSEMBLERS_PATH_ROOT_ABS;
    TOOLS_ROOT = basicdefines.TOOLS_ROOT;
except:
    MODULE_BASICDEFINES = False;
    ASSEMBLERS_PATH_ROOT_ABS = os.path.join(SCRIPT_PATH, 'assemblers/');
    TOOLS_ROOT = '%s' % (SCRIPT_PATH);

ASSEMBLER_URL = 'http://spades.bioinf.spbau.ru/release3.6.1/SPAdes-3.6.1-Linux.tar.gz'
ASSEMBLER_PATH = os.path.join(ASSEMBLERS_PATH_ROOT_ABS, 'SPAdes-3.6.1-Linux')
ZIP_FILE = os.path.join(ASSEMBLERS_PATH_ROOT_ABS, os.path.basename(ASSEMBLER_URL));
ASSEMBLER_BIN = os.path.join(ASSEMBLER_PATH,'bin/spades.py')
ASSEMBLER_NAME = 'SPADES'
# ASSEMBLER_RESULTS = 'contigs.fasta'
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
    sys.stderr.write('Executing command: "%s"\n' % command);
    if (dry_run == True):
        log('Executing (dryrun): "%s".' % (command), fp_log);
    if (dry_run == False):
        log('Executing: "%s".' % (command), fp_log);
        p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
    [output, err] = p.communicate()
    rc = p.returncode
    sys.stderr.write('\n');
    if (rc != 0):
        log('ERROR: subprocess call returned error code: %d.\n' % (rc), fp_log);
        traceback.print_stack(fp_log);
        exit(1);
    return [rc, output, err];

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
            final_max_rss += maxrss;
            final_time_unit = time_unit;
            final_mem_unit = mem_unit;
        else:
            if (time_unit == final_time_unit and mem_unit == final_mem_unit):
                final_command_line += '; ' + cmdline;
                final_real_time += realtime;
                final_cpu_time += cputime;
                final_user_time += usertime;
                final_system_time += systemtime;
                final_max_rss += maxrss;
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
# def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):
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
    ### Start the important work.
    ##################################################################################
    log('Running assembly using %s.' % (ASSEMBLER_NAME), fp_log);

    yaml_lines = [];
    yaml_lines.append('[');
    current_dataset = 0;
    for dataset in datasets:
        current_dataset += 1;
        if (dataset.type == 'single'):
            pass;

        elif (dataset.type == 'paired'):
            yaml_lines.append(' {');
            yaml_lines.append('  orientation: "fr",');
            yaml_lines.append('  type: "paired-end",');
            yaml_lines.append('  right reads: [');
            yaml_lines.append('   "%s"' % (dataset.reads_path_b));
            yaml_lines.append('  ],');
            yaml_lines.append('  left reads: [');
            yaml_lines.append('   "%s"' % (dataset.reads_path_a));
            yaml_lines.append('  ]');
            if (current_dataset == len(datasets)):
                yaml_lines.append(' }');
            else:
                yaml_lines.append(' },');

        elif (dataset.type == 'mate'):
            yaml_lines.append(' {');
            yaml_lines.append('  orientation: "rf",');
            yaml_lines.append('  type: "mate-pairs",');
            yaml_lines.append('  right reads: [');
            yaml_lines.append('   "%s"' % (dataset.reads_path_b));
            yaml_lines.append('  ],');
            yaml_lines.append('  left reads: [');
            yaml_lines.append('   "%s"' % (dataset.reads_path_a));
            yaml_lines.append('  ]');
            if (current_dataset == len(datasets)):
                yaml_lines.append(' }');
            else:
                yaml_lines.append(' },');

        elif (dataset.type == 'pacbio'):
            yaml_lines.append(' {');
            yaml_lines.append('  type: "pacbio",');
            yaml_lines.append('  single reads: [');
            yaml_lines.append('   "%s"' % (dataset.reads_path));
            yaml_lines.append('  ]');
            if (current_dataset == len(datasets)):
                yaml_lines.append(' }');
            else:
                yaml_lines.append(' },');

        elif (dataset.type == 'nanopore'):
            yaml_lines.append(' {');
            yaml_lines.append('  type: "nanopore",');
            yaml_lines.append('  single reads: [');
            yaml_lines.append('   "%s"' % (dataset.reads_path));
            yaml_lines.append('  ]');
            if (current_dataset == len(datasets)):
                yaml_lines.append(' }');
            else:
                yaml_lines.append(' },');

    yaml_lines.append(']');

    yaml = '\n'.join(yaml_lines);

    yaml_file_path = '%s/datasets.yaml' % (output_path);
    try:
        fp_yaml = open(yaml_file_path, 'w');
        fp_yaml.write(yaml + '\n');
        fp_yaml.close();
    except IOError, e:
        log('ERROR: Could not generate spec file in path: "%s"! Exiting.\n' % (yaml_file_path), fp_log);
        log(str(e), fp_log);

    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '.memtime')
    command = 'cd %s; %s %s -o %s/assembly_results --dataset %s' % (output_path, measure_command(memtime_path), ASSEMBLER_BIN, output_path, yaml_file_path)
    execute_command(command, fp_log, dry_run=DRY_RUN);

    command = 'cp %s/assembly_results/scaffolds.fasta %s/%s' % (output_path, output_path, ASSEMBLY_UNPOLISHED);
    execute_command(command, fp_log, dry_run=DRY_RUN);
    


    # if machine_name == 'pacbio':

    # elif machine_name == 'nanopore':
    # elif machine_name == 'illumina':
    #     log('\nMachine name "%s" not implemented for %s.\n' % (machine_name, ASSEMBLER_NAME));
    #     log('Skipping ....\n', fp_log)
    #     return;
    # else:
    #     log('\nMachine name "%s" not implemented for %s.\n' % (machine_name, ASSEMBLER_NAME));
    #     log('Skipping ....\n', fp_log)
    #     return;

    # if reads_file[-3:] == '.fq' or reads_file[-6:] == '.fastq':
    #     command = '%s %s -t %d --s1 %s -o %s' % (measure_command(memtime_path), ASSEMBLER_BIN, num_threads, reads_file, output_path)
    #     subprocess.call(command, shell='True')
    # # If reads file is fasta (spades cannot perform error correction)
    # elif reads_file[-3:] == '.fa' or reads_file[-6:] == '.fasta':
    #     command = '%s %s --only-assembler -t %d --s1 %s -o %s' % (measure_command(memtime_path), ASSEMBLER_BIN, num_threads, reads_file, output_path)
    #     subprocess.call(command, shell='True')
    # else:
    #     sys.stderr.write('\n[%s wrapper] Unsupported file format (%s)!\n' % (ASSEMBLER_NAME, reads_file))



    if (fp_log != None):
        fp_log.close();


# A placeholder for a function that runs quast on assembly results
# This is the same for all wrappers at the moment and is done in the main program
def run_quast():
    pass

# A function that gets the results of cgmemtime for a wrapper
# Since some aligners are run multiple times and create multiple measurements
# this is specific for each wrapper
def get_memtime():
    pass


# This is a standard interface for setting up an assembler. It should assume that the assembler
# is not present localy, but needs to be retrieved, unpacked, compiled and set-up, without requireing
# root privileges.
def download_and_install():
    if os.path.exists(ASSEMBLER_BIN):
        sys.stderr.write('[%s wrapper] Bin found at %s. Skipping installation ... \n' % (ASSEMBLER_NAME, ASSEMBLER_BIN))
    else:
        # sys.stderr.write('[%s wrapper] Couldn\'t find SPADES installation.\n')
        # sys.stderr.write('[%s wrapper] Install SPADES manually to %s.\n' % (basicdefines.ASSEMBLERS_PATH_ROOT_ABS, basicdefines.ASSEMBLERS_PATH_ROOT_ABS))
        # raw_input('Press key to continue...\n')
        if not os.path.exists(ZIP_FILE):
            sys.stderr.write('[%s wrapper] Downloading archive...\n' % (ASSEMBLER_NAME))
            command = 'cd %s; wget %s' % (ASSEMBLERS_PATH_ROOT_ABS, ASSEMBLER_URL)
            sys.stderr.write('[%s wrapper] %s\n' % (ASSEMBLER_NAME, command))
            subprocess.call(command, shell='True')

        # Decompress
        command = 'cd %s; tar -xvzf %s' % (ASSEMBLERS_PATH_ROOT_ABS, ZIP_FILE)
        sys.stderr.write('[%s wrapper] %s\n' % (ASSEMBLER_NAME, command))
        subprocess.call(command, shell='True')

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
