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


ASSEMBLER_URL = 'git://github.com/PacificBiosciences/FALCON-integrate.git'
ASSEMBLER_PATH = os.path.join(ASSEMBLERS_PATH_ROOT_ABS, 'FALCON-integrate')
ASSEMBLER_BIN = os.path.join(ASSEMBLER_PATH, 'fc_env/bin/fc_run.py')
ASSEMBLER_NAME = 'FALCON'
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
#     log('Captured STDOUT:\n%s\n\n' % (output), fp_log);
#     log('Captured STDERR:\n%s\n\n' % (err), fp_log);
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



### DALIGNER requires PacBio format of headers in reads.
### Parameter read_count_offset specifies the ID of the first read in the file.
### For every read, this ID gets updated.
### This ID is given as the first added parameter to the PacBio header.
### The new header will be formatted as: old_header/read_id/0_readlen .
### IMPORTANT: Falcon requires FASTA format of input files, so this function converts
### from FASTQ to FASTA in case reads_file was FASTQ.
def convert_reads_to_pacbio_format(reads_file, out_reads_file, fp_log, read_count_offset=0):
    try:
        fp_in = open(reads_file, 'r');
    except:
        log('ERROR: Could not open file "%s" for reading! Exiting.\n' % (reads_file), fp_log);
        exit(0);

    try:
        fp_out = open(out_reads_file, 'w');
    except:
        log('ERROR: Could not open file "%s" for writing! Exiting.\n' % (out_reads_file), fp_log);
        exit(0);

    current_read = read_count_offset;
    header_conversion_hash = {};

    while True:
        [header, read] = get_single_read(fp_in);
        
        if (len(read) == 0):
            break;

        current_read += 1;

        if (len(read[1]) <= 20):    ### DALIGNER has a lower length limit of 10bp.
            log('Found a read shorter than 20bp. Removing from the output.\n', fp_log);
            continue;

        ### Check if the read is already formatted like PacBio.
        if (re.match('^(.*?)/([\d]+)/([\d]+)_([\d]+)$', header) != None):
            read[1] = re.sub("(.{500})", "\\1\n", read[1], 0, re.DOTALL);   ### Wrap the sequence line, because DALIGNER has a 9998bp line len limit.
            fp_out.write('>%s\n' % (read[0]));
            fp_out.write('%s\n' % (read[1]));   ### Wrap the sequence line, because DALIGNER has a 9998bp line len limit.
            continue;

        ### Take only the part of the header up to the first whitespace, and replace all existing '/' with '_' to reduce chances of
        ### incompatibility with weird tools.
        trimmed_header = re.sub('[^0-9a-zA-Z]', '_', header.split()[0]); # re.sub("[|:", "_", read[0][1:]);

        # pacbio_header = '%s/%d/0_%d RQ=0.850' % (trimmed_header, current_read, len(read[1]));
        pacbio_header = 'S1/%d/0_%d' % (current_read, len(read[1]));
        header_conversion_hash[pacbio_header] = header;
        read[0] = '%s%s' % (read[0][0], pacbio_header); ### Keep the first char of the header line.
        read[1] = re.sub("(.{500})", "\\1\n", read[1], 0, re.DOTALL);   ### Wrap the sequence line, because DALIGNER has a 9998bp line len limit.
#        if (len(read) == 4):
#            read[3] = re.sub("(.{500})", "\\1\n", read[3], 0, re.DOTALL);   ### Wrap the qual line, because DALIGNER has a 9998bp line len limit.
#        fp_out.write('\n'.join(read) + '\n');
        fp_out.write('>%s\n' % (pacbio_header));
        fp_out.write('%s\n' % (read[1]));   ### Wrap the sequence line, because DALIGNER has a 9998bp line len limit.

    log('\n', fp_log);
    fp_in.close();
    fp_out.close();

    last_read_offset = current_read;

    return [header_conversion_hash, last_read_offset];



# Function 'run' should provide a standard interface for running a mapper. Given input parameters, it should run the
# alignment process, and convert any custom output results to the SAM format. Function should return a string with the
# path to the output file.
#    reads_file            Path to a FASTA/FASTQ file containing reads.
#    reference_file        Path to a reference genome FASTA file.
#    machine_name        A symbolic name to specify a set of parameters for a specific sequencing platform.
#    output_path            Folder to which the output will be placed to. Filename will be automatically generated according to the name of the mapper being run.
#    output_suffix        A custom suffix that can be added to the output filename.
# def run(reads_files, reference_file, machine_name, output_path, output_suffix=''):
#    approx_genome_len      Some assemblers need this parameter in order to run. For others, it's not necessary.
#    rerun                  If True, the specified output path will be moved to a backup location, and a new output folder created. Otherwise, the same folder will be used if it exists.
def run(datasets, output_path, approx_genome_len=0, move_exiting_out_path=True, run_type='run'):
    ##################################################################################
    ### Sanity check for input datasets.
    ##################################################################################    
    machine_name = None
    for dataset in datasets:
        if (machine_name != None and dataset.type != machine_name):
            sys.stderr.write('ERROR: %s is not a hybrid assembler, but datasets from disparate technologies are specified! Exiting.\n');
            exit(1);
        machine_name = dataset.type;
    if (machine_name == None):
        sys.stderr.write('ERROR: Input datasets not specified correctly! Exiting.\n');
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
    ### Prepare a holder for all reads paths.
    fofn_paths = [];

    ### Convert read headers to PacBio format so that DALIGNER can process them.
    # sys.stderr.write('[%s wrapper] Converting reads files:\n' % (ASSEMBLER_NAME))
    log('Converting reads files:', fp_log);
    last_read_offset = 0;
    all_header_conversion_hash = {};
    i = 0;
    for dataset in datasets:
        i += 1;
        reads_file = dataset.reads_path;
        reads_file = os.path.abspath(reads_file);
        out_reads_file = '%s/%s.pbheader.fasta' % (output_path, os.path.splitext(os.path.basename(reads_file))[0]);
        # sys.stderr.write('[%s wrapper] \t(%d) %s -> %s\n' % (ASSEMBLER_NAME, i, reads_file, out_reads_file))
        log('\t(%d) %s -> %s' % (i, reads_file, out_reads_file), fp_log);
        fofn_paths.append(out_reads_file);
        if (DRY_RUN == False):
            [header_conversion_hash, last_read_offset] = convert_reads_to_pacbio_format(reads_file, out_reads_file, fp_log, read_count_offset=last_read_offset);
            all_header_conversion_hash.update(header_conversion_hash);

    ### Create a .fofn file containing paths to all reads files.
    fofn_file = '%s/input.fofn' % (output_path);
    try:
        fp_fofn = open(fofn_file, 'w');
    except Exception, e:
        log('ERROR: Could not open file "%s" writing!\n' % (fofn_file), fp_log);
        return;
    fp_fofn.write('\n'.join(fofn_paths) + '\n');
    fp_fofn.close();

    ##################################################################################
    ### Prepare the config file.
    ##################################################################################
    log('Running assembly using %s.' % (ASSEMBLER_NAME), fp_log);

    ### Create a config file for the assembly.
    if (machine_name == 'pacbio'):
        cfg_lines = '';
        cfg_lines += '[General]\n';
        cfg_lines += 'job_type = local\n';
        cfg_lines += '# list of files of the initial bas.h5 files\n';
        cfg_lines += 'input_fofn = %s\n' % (fofn_file);
        cfg_lines += '\n';
        cfg_lines += 'input_type = raw\n';
        if (run_type == 'erc'):
                cfg_lines += 'target = pre-assembly\n';
        else:
                cfg_lines += 'target = assembly\n';
        cfg_lines += '#input_type = preads\n';
        cfg_lines += '\n';
        cfg_lines += '# The length cutoff used for seed reads used for initial mapping\n';
        cfg_lines += 'length_cutoff = 2000\n';
        cfg_lines += '\n';
        cfg_lines += '# The length cutoff used for seed reads usef for pre-assembly\n';
        cfg_lines += 'length_cutoff_pr = 2000\n';
        cfg_lines += '\n';
        cfg_lines += '# job_type= local\n';
        cfg_lines += 'jobqueue = your_queue\n';
        cfg_lines += 'sge_option_da =\n';
        cfg_lines += 'sge_option_la =\n';
        cfg_lines += 'sge_option_pda =\n';
        cfg_lines += 'sge_option_pla =\n';
        cfg_lines += 'sge_option_fc =\n';
        cfg_lines += 'sge_option_cns =\n';
        cfg_lines += '\n';
        cfg_lines += 'pa_concurrent_jobs = 24\n';
        cfg_lines += 'ovlp_concurrent_jobs = 24\n';
        cfg_lines += '\n';
        cfg_lines += 'pa_HPCdaligner_option =  -v -dal4 -t16 -e.70 -l1000 -s1000\n';
        cfg_lines += 'ovlp_HPCdaligner_option = -v -dal4 -t32 -h60 -e.96 -l500 -s1000\n';
        cfg_lines += '\n';
        cfg_lines += 'pa_DBsplit_option = -x500 -s50\n';
        cfg_lines += 'ovlp_DBsplit_option = -x500 -s50\n';
        cfg_lines += '\n';
        cfg_lines += 'falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --local_match_count_threshold 2 --max_n_read 200 --n_core %d\n' % (num_threads);
        cfg_lines += '\n';
        cfg_lines += 'overlap_filtering_setting = --max_diff 100 --max_cov 100 --min_cov 20 --bestn 10 --n_core %d\n' % (num_threads);
        cfg_lines += '\n';
        cfg_lines += '# Running a PacBio reads configuration.\n';

    elif (machine_name == 'nanopore' or machine_name == 'nanoporenocorrect'):
        cfg_lines = '';
        cfg_lines += '[General]\n';
        cfg_lines += 'job_type = local\n';
        cfg_lines += '# list of files of the initial bas.h5 files\n';
        cfg_lines += 'input_fofn = %s\n' % (fofn_file);
        cfg_lines += '\n';
        if (machine_name == 'nanoporenocorrect'):
            cfg_lines += 'input_type = preads\n';
        else:
            cfg_lines += 'input_type = raw\n';
        if (run_type == 'erc'):
                cfg_lines += 'target = pre-assembly\n';
        else:
                cfg_lines += 'target = assembly\n';
        cfg_lines += '\n';
        cfg_lines += '# The length cutoff used for seed reads used for initial mapping\n';
        cfg_lines += 'length_cutoff = 1000\n';
        cfg_lines += '\n';
        cfg_lines += '# The length cutoff used for seed reads usef for pre-assembly\n';
        cfg_lines += 'length_cutoff_pr = 1000\n';
        cfg_lines += '\n';
        cfg_lines += '# job_type= local\n';
        cfg_lines += 'jobqueue = your_queue\n';
        cfg_lines += 'sge_option_da =\n';
        cfg_lines += 'sge_option_la =\n';
        cfg_lines += 'sge_option_pda =\n';
        cfg_lines += 'sge_option_pla =\n';
        cfg_lines += 'sge_option_fc =\n';
        cfg_lines += 'sge_option_cns =\n';
        cfg_lines += '\n';
        cfg_lines += 'pa_concurrent_jobs = 24\n';
        cfg_lines += 'ovlp_concurrent_jobs = 24\n';
        cfg_lines += '\n';

        ### DALIGNER options:
        ### -dal4   The -dal option (default 4) gives the desired number of block comparisons per call to
        ###         daligner. Some must contain dal-1 comparisons, and the first dal-2 block comparisons
        ###         even less, but the HPCdaligner "planner" does the best it can to give an average load
        ###         of dal block comparisons per command.
        ### -t parameter which suppresses the use of any k-mer that occurs more than t times in either the subject or target block
        ### -h the total number of bases covered by the k-mer hits is h (default 35)
        ### -l searching for local alignments involving at least -l base pairs (default 1000)
        ### -s The local alignments found will be output in a sparse encoding where a trace point on the alignment is recorded every -s base pairs of the a-read (default 100bp)
        ### 
        ### The options -k, -h, and -w control the initial filtration search for possible matches
        ### between reads.  Specifically, our search code looks for a pair of diagonal bands of
        ### width 2^w (default 2^6 = 64) that contain a collection of exact matching k-mers
        ### (default 14) between the two reads, such that the total number of bases covered by the
        ### k-mer hits is h (default 35). k cannot be larger than 32 in the current implementation.
        cfg_lines += 'pa_HPCdaligner_option =  -v -dal4 -t100 -e.70 -l100 -s100\n';
        cfg_lines += 'ovlp_HPCdaligner_option = -v -dal4 -t100 -h60 -e.92 -l100 -s100\n';
        cfg_lines += '\n';
        cfg_lines += 'pa_DBsplit_option = -x100 -s50\n';
        cfg_lines += 'ovlp_DBsplit_option = -x100 -s50\n';
        #        cfg_lines += 'falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --local_match_count_threshold 2 --max_n_read 200 --n_core 6\n';
        #        cfg_lines += 'overlap_filtering_setting = --max_diff 100 --max_cov 100 --min_cov 20 --bestn 10 --n_core 24\n';

        # --max_n_read put a cap on the number of reads used for error correction. In high repetitive genome, you will need to put smaller --max_n_read to make sure the consensus code does not waste time aligning repeats.
        cfg_lines += 'falcon_sense_option = --output_multi --min_idt 0.50 --local_match_count_threshold 0 --max_n_read 200 --n_core %d\n' % (num_threads);

        ### The --max_diff parameter can be used to filter out the reads where one ends has much more coverage than the other end.
        ### The --max_cov and --min_cov are used for filtering reads that have too high or too low overlaps.
        ### The --bestn parameter in overlap_filtering_setting option can be used to control the maximum overlap reported for each read.
        cfg_lines += 'overlap_filtering_setting = --max_diff 100 --max_cov 100 --min_cov 5 --bestn 20 --n_core %d\n' % (num_threads);

        cfg_lines += '\n';
        cfg_lines += '# Running an Oxford Nanopore reads configuration.\n';

    else:
        log('ERROR: Unknown machine name: "%s". Skipping assembly.\n' % (machine_name), fp_log);
        return;

    cfg_file = '%s/fc_run.cfg' % (output_path);
    try:
        fp_cfg = open(cfg_file, 'w');
    except Exception, e:
        log('ERROR: Could not open file "%s" writing!\n' % (cfg_file), fp_log);
        return;
    fp_cfg.write(cfg_lines);
    fp_cfg.close();

    #####

    ##################################################################################
    ### Start the important work.
    ##################################################################################
    # if (MODULE_BASICDEFINES == True):
    #     command = 'sudo %s/cgmemtime/cgmemtime --setup -g %s --perm 775' % (basicdefines.TOOLS_ROOT, getpass.getuser());
    #     execute_command(command, fp_log, dry_run=DRY_RUN);

    # memtime_file = '%s/%s.memtime' % (output_path, ASSEMBLER_NAME);
    memtime_file = '%s/total.memtime' % (output_path);
    FC_path = '%s/fc_env' % (ASSEMBLER_PATH);
    run_commands = [];
    run_commands.append('. %s/bin/activate' % (FC_path));
    run_commands.append('cd %s' % (output_path));
    run_commands.append('%s %s/bin/fc_run.py %s' % (measure_command(memtime_file), FC_path, cfg_file));
    if (run_type != 'erc'):
        run_commands.append('cp %s/2-asm-falcon/p_ctg.fa %s/%s' % (output_path, output_path, ASSEMBLY_UNPOLISHED));
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
        # sys.stderr.write('[%s wrapper] Bin found at %s. Skipping installation.\n' % (ASSEMBLER_NAME, ASSEMBLER_BIN))
        log('Bin found at %s. Skipping installation.' % (ASSEMBLER_BIN), None);

    else:
        # sys.stderr.write('[%s wrapper] Started installation of %s.\n' % (ASSEMBLER_NAME, ASSEMBLER_NAME))
        log('Started installation of %s.' % (ASSEMBLER_NAME), None);
	try:
            import virtualenv;
	except:
            log('Virtualenv not installed! Sudo access will be required.\n', None);
            execute_command('sudo pip install virtualenv', None, dry_run=DRY_RUN);
            
        # Falcon requires the python-dev package to be installed.
        from distutils import sysconfig as s; python_path = s.get_config_vars()['INCLUDEPY'];
        if (os.path.exists('%s/Python.h' % (python_path)) == False):
            execute_command('sudo apt-get install python-dev', None, dry_run=DRY_RUN);




        if (not os.path.exists(ASSEMBLERS_PATH_ROOT_ABS)):
            # sys.stderr.write('[%s wrapper] Creating a directory on path "%s".\n' % (ASSEMBLER_NAME, ASSEMBLERS_PATH_ROOT_ABS))
            log('Creating a directory on path "%s".' % (ASSEMBLERS_PATH_ROOT_ABS), None);
            os.makedirs(ASSEMBLERS_PATH_ROOT_ABS);

        command = 'cd %s; git clone %s; cd %s; git checkout 8bb2737fd1d70e8914ade7479b434cd1223d5009' % (ASSEMBLERS_PATH_ROOT_ABS, ASSEMBLER_URL, ASSEMBLER_PATH)
        execute_command(command, None, dry_run=DRY_RUN);

        setup_commands = [];
        setup_commands.append('cd %s' % (ASSEMBLER_PATH));
        setup_commands.append('WORK=$PWD');
        setup_commands.append('FC=$WORK/fc_env');
        setup_commands.append('virtualenv --no-site-packages  --always-copy   $FC');
        setup_commands.append('. $FC/bin/activate');
        setup_commands.append('git submodule update --init');
        setup_commands.append('cd pypeFLOW');
        setup_commands.append('python setup.py install');
        setup_commands.append('cd ../FALCON');
        setup_commands.append('python setup.py install');
        setup_commands.append('cd ../DAZZ_DB/');
        setup_commands.append('make');
        setup_commands.append('cp DBrm DBshow DBsplit DBstats fasta2DB $FC/bin/');
        setup_commands.append('cd ../DALIGNER');
        setup_commands.append('make');
        setup_commands.append('cp daligner daligner_p DB2Falcon HPCdaligner LA4Falcon LAmerge LAsort  $FC/bin');
        command = '; '.join(setup_commands);
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

    elif (sys.argv[1] == 'run' or sys.argv[1] == 'erc'):
        if (len(sys.argv) < 5):
            verbose_usage_and_exit()

        output_path = sys.argv[2]
        approx_genome_len = 0 if (sys.argv[3] == '-') else int(sys.argv[3]);
        datasets = [];
        for arg in sys.argv[4:]:
            dataset = Dataset(arg);
            datasets.append(dataset);
        run(datasets, output_path, approx_genome_len, run_type=sys.argv[1]);

    else:
        verbose_usage_and_exit()
