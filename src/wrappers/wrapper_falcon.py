#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(SCRIPT_PATH + '/../')

import subprocess
import multiprocessing
import getpass

try:
    import basicdefines
    MODULE_BASICDEFINES = True;
    ASSEMBLERS_PATH_ROOT_ABS = basicdefines.ASSEMBLERS_PATH_ROOT_ABS;
    TOOLS_ROOT = basicdefines.TOOLS_ROOT;
except:
    MODULE_BASICDEFINES = False;
    ASSEMBLERS_PATH_ROOT_ABS = os.path.join(SCRIPT_PATH, 'assemblers/');
    TOOLS_ROOT = '%s' % (SCRIPT_PATH);


ASSEMBLER_URL = 'git://github.com/PacificBiosciences/FALCON-integrate.git'
ASSEMBLER_PATH = os.path.join(ASSEMBLERS_PATH_ROOT_ABS, 'FALCON-integrate')
ASSEMBLER_BIN = os.path.join(ASSEMBLER_PATH, 'fc_env/bin/fc_run.py')
ASSEMBLER_NAME = 'FALCON'
# ASSEMBLER_RESULTS = 'contig-100.fa'
CREATE_OUTPUT_FOLDER = True

# DRY_RUN = True;
DRY_RUN = False;

def execute_command(command, dry_run=True):
    if (dry_run == True):
        # sys.stderr.write('Warning: dry_run == True\n');
        sys.stderr.write('[%s wrapper] Executing (dryrun): "%s"\n' % (ASSEMBLER_NAME, command));
    if (dry_run == False):
        sys.stderr.write('[%s wrapper] Executing: "%s"\n' % (ASSEMBLER_NAME, command));
        subprocess.call(command, shell=True);
    sys.stderr.write('\n');

def measure_command(measure_file):
    if (MODULE_BASICDEFINES == True):
        return basicdefines.measure_command(measure_file);
    else:
        sys.stderr.write('ERROR: Cgmemtime tool not found! Exiting.\n');
        exit(1);
        # return '/usr/bin/time --format "Command line: %%C\\nReal time: %%e s\\nCPU time: -1.0 s\\nUser time: %%U s\\nSystem time: %%S s\\nMaximum RSS: %%M kB\\nExit status: %%x" --quiet -o %s ' % measure_file;



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

    sys.stderr.write('[%s wrapper] Running assembly using %s.\n' % (ASSEMBLER_NAME, ASSEMBLER_NAME))

    reference_file = os.path.abspath(reference_file);
    output_path = os.path.abspath(output_path);

    if (not os.path.exists(output_path)):
        sys.stderr.write('[%s wrapper] Creating a directory on path "%s".\n' % (ASSEMBLER_NAME, output_path))
        os.makedirs(output_path);

    ### Create a .fofn file containing paths to all reads files.
    fofn_file = '%s/input.fofn' % (output_path);
    try:
        fp_fofn = open(fofn_file, 'w');
    except Exception, e:
        sys.stderr.write('ERROR: Could not open file "%s" writing!\n' % (fofn_file));
        return;
    sys.stderr.write('[%s wrapper] All reads files:\n' % (ASSEMBLER_NAME))
    i = 0;
    for reads_file in reads_files:
        i += 1;
        fp_fofn.write(os.path.abspath(reads_file) + '\n');
        sys.stderr.write('[%s wrapper] \t(%d) %s\n' % (ASSEMBLER_NAME, i, reads_file))
    fp_fofn.close();

    ### Create a config file for the assembly.
    if (machine_name == 'pacbio'):
        cfg_lines = '';
        cfg_lines += '[General]\n';
        cfg_lines += 'job_type = local\n';
        cfg_lines += '# list of files of the initial bas.h5 files\n';
        cfg_lines += 'input_fofn = %s\n' % (fofn_file);
        cfg_lines += '\n';
        cfg_lines += 'input_type = raw\n';
        cfg_lines += '#input_type = preads\n';
        cfg_lines += '\n';
        cfg_lines += '# The length cutoff used for seed reads used for initial mapping\n';
        cfg_lines += 'length_cutoff = 12000\n';
        cfg_lines += '\n';
        cfg_lines += '# The length cutoff used for seed reads usef for pre-assembly\n';
        cfg_lines += 'length_cutoff_pr = 12000\n';
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
        cfg_lines += 'falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --local_match_count_threshold 2 --max_n_read 200 --n_core 6\n';
        cfg_lines += '\n';
        cfg_lines += 'overlap_filtering_setting = --max_diff 100 --max_cov 100 --min_cov 20 --bestn 10 --n_core 24\n';
        cfg_lines += '\n';
        cfg_lines += '# Running a PacBio reads configuration.\n';
    else:
        sys.stderr.write('ERROR: Unknown machine name: "%s". Skipping assembly.\n' % (machine_name));
        return;

    cfg_file = '%s/fc_run.cfg' % (output_path);
    try:
        fp_cfg = open(cfg_file, 'w');
    except Exception, e:
        sys.stderr.write('ERROR: Could not open file "%s" writing!\n' % (cfg_file));
        return;
    fp_cfg.write(cfg_lines);
    fp_cfg.close();

    #####

    if (MODULE_BASICDEFINES == True):
        command = 'sudo %s/cgmemtime/cgmemtime --setup -g %s --perm 775' % (basicdefines.TOOLS_ROOT, getpass.getuser());
        execute_command(command, dry_run=DRY_RUN);

    memtime_file = '%s/%s.memtime' % (output_path, ASSEMBLER_NAME);
    FC_path = '%s/fc_env' % (ASSEMBLER_PATH);
    setup_commands = [];
    setup_commands.append('. %s/bin/activate' % (FC_path));
    setup_commands.append('cd %s' % (output_path));
    setup_commands.append('%s fc_run.py %s' % (measure_command(memtime_file), cfg_file));
    command = '; '.join(setup_commands);
    execute_command(command, dry_run=DRY_RUN);



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
        sys.stderr.write('[%s wrapper] Bin found at %s. Skipping installation.\n' % (ASSEMBLER_NAME, ASSEMBLER_BIN))
    else:
        sys.stderr.write('[%s wrapper] Started installation of %s.\n' % (ASSEMBLER_NAME, ASSEMBLER_NAME))
        if (not os.path.exists(ASSEMBLERS_PATH_ROOT_ABS)):
            sys.stderr.write('[%s wrapper] Creating a directory on path "%s".\n' % (ASSEMBLER_NAME, ASSEMBLERS_PATH_ROOT_ABS))
            os.makedirs(ASSEMBLERS_PATH_ROOT_ABS);

        command = 'cd %s; git clone %s' % (ASSEMBLERS_PATH_ROOT_ABS, ASSEMBLER_URL)
        execute_command(command, dry_run=DRY_RUN);

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
        execute_command(command, dry_run=DRY_RUN);

def verbose_usage_and_exit():
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s mode [<reads_file1>,<reads_file2>,...,<reads_fileN> <machine_name> <output_path> <reference_file> [<output_suffix>]]\n' % sys.argv[0])
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
        if (len(sys.argv) < 5 or len(sys.argv) > 7):
            verbose_usage_and_exit()

        reads_files = sys.argv[2].split(',')         ### Enable specifying multiple FASTQ files for input.
        machine_name = sys.argv[3]
        output_path = sys.argv[4]
        reference_file = '';
        output_suffix = ''

        if (len(sys.argv) >= 6):
            reference_file = sys.argv[5]
        if (len(sys.argv) == 7):
            output_suffix = sys.argv[6]
        run(reads_files, reference_file, machine_name, output_path, output_suffix)

    else:
        verbose_usage_and_exit()
