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

ASSEMBLER_URL = 'https://github.com/lh3/miniasm.git'
ASSEMBLER_PATH = os.path.join(ASSEMBLERS_PATH_ROOT_ABS, 'miniasmraconoverlap')
ASSEMBLER_BIN = os.path.join(ASSEMBLER_PATH,'miniasm/miniasm')
ASSEMBLER_NAME = 'MiniasmRaconOverlap'
RACON_BIN = os.path.join(ASSEMBLER_PATH,'racon/bin/racon')
# PAF2MHAP_BIN = os.path.join(ASSEMBLER_PATH,'racon/scripts/paf2mhap.pl')
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

def fastq2fasta(input_fastq_path, output_fasta_path):
    try:
        fp_in = open(input_fastq_path, 'r');
    except:
        sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
        exit(0);

    try:
        fp_out = open(output_fasta_path, 'w');
    except:
        sys.stderr.write('ERROR: Could not open file "%s" for writing! Exiting.\n' % output_fasta_path);
        exit(0);

    i = 0;
    while True:
        i += 1;
        [header, read] = get_single_read(fp_in);
        if (len(read) == 0):
            break;
        fp_out.write('>%s\n%s\n' % (read[0][1:], read[1]));

    sys.stderr.write('\n');
    fp_in.close();
    fp_out.close();



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
    reads_file_fasta = '%s/all_reads.fasta' % (output_path);
    for dataset in datasets:
        if (dataset.reads_path.endswith('fasta') or dataset.reads_path.endswith('fa')):
            # log('ERROR: Assembler %s expects only FASTQ files for input. Trouble loading "%s". Exiting.\n' % (ASSEMBLER_NAME, dataset.reads_path), fp_log);
            converted_reads_path = '%s/%s.fastq' % (output_path, os.path.splitext(os.path.basename(dataset.reads_path))[0]);
            log('Converting file "%s" to FASTQ format and aggregating to "%s".\n' % (dataset.reads_path, reads_file), fp_log);
            command = 'java -jar %s/convertFastaAndQualToFastq.jar %s >> %s' % (ASSEMBLER_PATH, dataset.reads_path, reads_file);
            execute_command(command, fp_log, dry_run=DRY_RUN);
        else:
            log('Aggregating FASTQ file "%s" to "%s".\n' % (dataset.reads_path, reads_file), fp_log);
            command = 'cat %s >> %s' % (dataset.reads_path, reads_file);
            execute_command(command, fp_log, dry_run=DRY_RUN);
    fastq2fasta(reads_file, reads_file_fasta);

    ##################################################################################
    ### Start the important work.
    ##################################################################################
    memtime_path = os.path.join(output_path, 'total.memtime')
    memtime_files_prefix =  '%s/%s' % (output_path, ASSEMBLER_NAME);
    current_memtime_id = 0;

    spec_file = ''
    if machine_name == 'pacbio' or machine_name == 'nanopore':
        overlaps_file = '%s/overlaps.paf.gz' % (output_path);
        current_memtime_id += 1;
        command = '%s %s/minimap/minimap -Sw5 -L100 -m0 -t%d %s %s | gzip -1 > %s' % (measure_command('%s-%s.memtime' % (memtime_files_prefix, current_memtime_id)), ASSEMBLER_PATH, num_threads, reads_file, reads_file, overlaps_file);
        execute_command(command, fp_log, dry_run=DRY_RUN);

        assembly_raw_gfa = '%s/miniasm.gfa' % (output_path);
        current_memtime_id += 1;
        command = '%s %s -f %s %s > %s' % (measure_command('%s-%s.memtime' % (memtime_files_prefix, current_memtime_id)), ASSEMBLER_BIN, reads_file, overlaps_file, assembly_raw_gfa);
        execute_command(command, fp_log, dry_run=DRY_RUN);

        assembly_iter0_fasta = '%s/consensus-iter0.fasta' % (output_path);
        assembly_iter1_fasta = '%s/consensus-iter1.fasta' % (output_path);
        assembly_iter2_fasta = '%s/consensus-iter2.fasta' % (output_path);
        mapping_iter1_paf = '%s/mapping-iter1.paf' % (output_path);
        mapping_iter1_mhap = '%s/mapping-iter1.mhap' % (output_path);
        mapping_iter2_paf = '%s/mapping-iter2.paf' % (output_path);
        mapping_iter2_mhap = '%s/mapping-iter2.mhap' % (output_path);

        command = "awk '$1 ~/S/ {print \">\"$2\"\\n\"$3}' %s > %s" % (assembly_raw_gfa, assembly_iter0_fasta);
        execute_command(command, fp_log, dry_run=DRY_RUN);

        # Run the first iteration of Racon.
        current_memtime_id += 1;
        command = '%s %s/minimap/minimap %s %s > %s' % (measure_command('%s-%s.memtime' % (memtime_files_prefix, current_memtime_id)), ASSEMBLER_PATH, assembly_iter0_fasta, reads_file, mapping_iter1_paf);
        execute_command(command, fp_log, dry_run=DRY_RUN);
        # command = '%s %s %s %s > %s' % (PAF2MHAP_BIN, assembly_iter0_fasta, reads_file_fasta, mapping_iter1_paf, mapping_iter1_mhap);
        # execute_command(command, fp_log, dry_run=DRY_RUN);
        current_memtime_id += 1;
        command = '%s %s -M 5 -X -4 -G -8 -E -6 --bq 10 -t %d %s %s %s %s' % (measure_command('%s-%s.memtime' % (memtime_files_prefix, current_memtime_id)), RACON_BIN, num_threads, reads_file, mapping_iter1_paf, assembly_iter0_fasta, assembly_iter1_fasta);
        execute_command(command, fp_log, dry_run=DRY_RUN);

        # Run the first iteration of Racon.
        current_memtime_id += 1;
        command = '%s %s/minimap/minimap %s %s > %s' % (measure_command('%s-%s.memtime' % (memtime_files_prefix, current_memtime_id)), ASSEMBLER_PATH, assembly_iter1_fasta, reads_file, mapping_iter2_paf);
        execute_command(command, fp_log, dry_run=DRY_RUN);
        # command = '%s %s %s %s > %s' % (PAF2MHAP_BIN, assembly_iter1_fasta, reads_file_fasta, mapping_iter2_paf, mapping_iter2_mhap);
        # execute_command(command, fp_log, dry_run=DRY_RUN);
        current_memtime_id += 1;
        command = '%s %s -M 5 -X -4 -G -8 -E -6 --bq 10 -t %d %s %s %s %s' % (measure_command('%s-%s.memtime' % (memtime_files_prefix, current_memtime_id)), RACON_BIN, num_threads, reads_file, mapping_iter2_paf, assembly_iter1_fasta, assembly_iter2_fasta);
        execute_command(command, fp_log, dry_run=DRY_RUN);

        command = 'cp %s %s/%s' % (assembly_iter2_fasta, output_path, ASSEMBLY_UNPOLISHED);
        execute_command(command, fp_log, dry_run=DRY_RUN);

        # command = 'awk \'$1 ~/S/ {print ">"$2"\\n"$3}\' %s > %s/%s' % (gfa_file, output_path, ASSEMBLY_UNPOLISHED);
        # execute_command(command, fp_log, dry_run=DRY_RUN);

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

    all_memtimes = ['%s-%s.memtime' % (memtime_files_prefix, value) for value in xrange(1, current_memtime_id+1)];
    parse_memtime_files_and_accumulate(all_memtimes, memtime_path);



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

        command = 'cd %s; git clone https://github.com/lh3/minimap && cd minimap && git checkout 1cd6ae3bc7c7a6f9e7c03c0b7a93a12647bba244 && make' % (ASSEMBLER_PATH)
        execute_command(command, None, dry_run=DRY_RUN);

        command = 'cd %s; git clone %s && cd miniasm && git checkout 17d5bd12290e0e8a48a5df5afaeaef4d171aa133 && make' % (ASSEMBLER_PATH, ASSEMBLER_URL)
        execute_command(command, None, dry_run=DRY_RUN);
        
        command = 'cd %s; git clone https://github.com/isovic/racon.git && (cd racon && git checkout overlap && make modules && make tools && make)' % (ASSEMBLER_PATH)
        execute_command(command, None, dry_run=DRY_RUN);

        command = 'cd %s; wget http://www.cbcb.umd.edu/software/PBcR/data/convertFastaAndQualToFastq.jar' % (ASSEMBLER_PATH);
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
