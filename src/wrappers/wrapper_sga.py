#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

import sys;
sys.path.append(SCRIPT_PATH + '/../src');

import subprocess
import multiprocessing

import basicdefines

ASSEMBLER_URL = 'https://github.com/jts/sga.git'
ASSEMBLER_PATH = os.path.join(basicdefines.ASSEMBLERS_PATH_ROOT_ABS, 'sga')
ASSEMBLER_BIN = '/usr/local/bin/sga'
ASSEMBLER_NAME = 'SGA'
# SGA doesn't have a set result filename, we have to define it ourselves
ASSEMBLER_RESULTS = 'assemble-contigs.fa'
CREATE_OUTPUT_FOLDER = True


BAMTOOLS_URL = 'http://github.com/pezmaster31/bamtools'
BAMTOOLS_PATH = os.path.join(basicdefines.TOOLS_ROOT_ABS, 'bamtools')
BAMTOOLS_BIN = os.path.join(BAMTOOLS_PATH, 'bin/bamtools')



# Function 'run' should provide a standard interface for running a mapper. Given input parameters, it should run the
# alignment process, and convert any custom output results to the SAM format. Function should return a string with the
# path to the output file.
#    reads_file            Path to a FASTA/FASTQ file containing reads.
#    reference_file        Path to a reference genome FASTA file.
#    machine_name          A symbolic name to specify a set of parameters for a specific sequencing platform.
#    output_path           Folder to which the output will be placed to. Filename will be automatically generated according to the name of the mapper being run.
#    output_suffix         A custom suffix that can be added to the output filename.
def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):
    # SGA has a rathar long series of steps to do to run an assembly
    # TODO: Here is one sequence of programs producing one results
    #        Parameters are many and they could all influence end result
    #        Some parameters should be inferred from the reads file or set by user

    # COMMENT: changing directory every time because it seems that it is not preserved
    #          across multiple shell commands

    num_threads = multiprocessing.cpu_count() / 2

    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '.memtime')
    # command = '%s %s --num_threads %d -r %s -o %s' % (basicdefines.measure_command(memtime_path), ASSEMBLER_BIN, num_threads, reads_file, output_path)
    # subprocess.call(command, shell='True')


    # 1. Preprocess
    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '_1.memtime')
    command = 'cd %s; %s %s preprocess -o sgaccs.fasta %s' % (output_path, basicdefines.measure_command(memtime_path), ASSEMBLER_BIN, reads_file)
    subprocess.call(command, shell='True')

    # Not doing error correction
    # # 2a. Build index for error correction
    # command = 'cd %s; %s index -a ropebwt -t 32 --no-reverse sgaccs.fasta' % (output_path, ASSEMBLER_BIN)
    # subprocess.call(command, shell='True')

    # # 2b. Perform error correction
    # command = 'cd %s; %s correct -k 21 --learn -t 32 -o reads.ec.k21.fasta sgaccs.fasta' % (output_path, ASSEMBLER_BIN)
    # subprocess.call(command, shell='True')

    # 3. Contig assembly

    #3a. Index data
    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '_2.memtime')
    command = 'cd %s; %s %s index -a ropebwt -t %d --no-reverse %s' % (output_path, basicdefines.measure_command(memtime_path), ASSEMBLER_BIN, num_threads, reads_file)
    subprocess.call(command, shell='True')

    # Not doing filtering
    # #3b. Remove exact-match duplicates and reads with low frequency kmers
    # # COMMENT: In my experience filtering could filter out too much
    # #          It might be better to skip it
    # #          Not sure how to decide when to skip and when not to
    # command = 'cd %s; %s filter -x 2 -t %d --homopolymer-check --low-complexity-check reads.ec.k21.fasta' % (output_path, ASSEMBLER_BIN, num_threads)
    # subprocess.call(command, shell='True')

    #3c. Merge simple, unbranched chains of vertices
    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '_3.memtime')
    command = 'cd %s; %s %s fm-merge -m 30 -t %d -o merged.k21.fa %s' % (output_path, basicdefines.measure_command(memtime_path), ASSEMBLER_BIN, num_threads, reads_file)
    subprocess.call(command, shell='True')

    # 3d. Build an index of the merged sequences
    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '_4.memtime')
    command = 'cd %s; %s %s index -d 1000000 -t %d merged.k21.fa' % (output_path, basicdefines.measure_command(memtime_path), ASSEMBLER_BIN, num_threads)
    subprocess.call(command, shell='True')

    # 3e. Remove any substrings that were generated from the merge process
    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '_5.memtime')
    command = 'cd %s; %s %s rmdup -t %d merged.k21.fa' % (output_path, basicdefines.measure_command(memtime_path), ASSEMBLER_BIN, num_threads)
    subprocess.call(command, shell='True')

    # 3f. Compute the structure of the string graph
    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '_6.memtime')
    command = 'cd %s; %s %s overlap -m 30 -t %d merged.k21.rmdup.fa' % (output_path, basicdefines.measure_command(memtime_path), ASSEMBLER_BIN, num_threads)
    subprocess.call(command, shell='True')

    # 3g. Perform the contig assembly without bubble popping
    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '_7.memtime')
    command = 'cd %s; %s %s assemble -m 30 -o assemble merged.k21.rmdup.asqg.gz' % (output_path, basicdefines.measure_command(memtime_path), ASSEMBLER_BIN)
    subprocess.call(command, shell='True')

    # Callculate memtime summary and write it in a file
    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '.memtime')
    real_time = 0.0
    cpu_time = 0.0
    user_time = 0.0
    system_time = 0.0
    max_rss = 0.0
    with open(memtime_path, 'w') as fmemtime:
        fmemtime.write('SGA summary .memtime file (%s)' % output_path)
        for i in xrange(1, 8):
            memtime_file = '%s_%d.memtime' % (ASSEMBLER_NAME, i)
            tmemtime_path = os.path.join(output_path, memtime_file)
            with open(tmemtime_path, 'r') as tfmemtime:
                tfmemtime.readline()            # skipping 1st line
                line = tfmemtime.readline()
                treal_time = float(line.split()[2])     # real time
                line = tfmemtime.readline()
                tcpu_time = float(line.split()[2])      # cpu time
                line = tfmemtime.readline()
                tuser_time = float(line.split()[2])     # user time
                line = tfmemtime.readline()
                tsystem_time = float(line.split()[2])   # system time
                line = tfmemtime.readline()
                tmax_rss = int(line.split()[2])         # Max RSS

                real_time += treal_time
                cpu_time += tcpu_time
                user_time += tuser_time
                system_time += tsystem_time
                if tmax_rss > max_rss:
                    max_rss = tmax_rss

        # writing to summary .memtime file
        fmemtime.write('\nReal time:  %.3f s' % real_time)
        fmemtime.write('\nCPU time:  %.3f s' % cpu_time)
        fmemtime.write('\nUser time:  %.3f s' % user_time)
        fmemtime.write('\nSystem time:  %.3f s' % system_time)
        fmemtime.write('\nMaximum RSS: %d MB' % max_rss)


# A placeholder for a function that runs quast on assembly results
# This is the same for all wrappers at the moment and is done in the main program
def run_quast():
    pass

# A function that gets the results of cgmemtime for a wrapper
# Since some aligners are run multiple times and create multiple measurements
# this is specific for each wrapper
def get_memtime():
    pass



# This is a standard interface for setting up the aligner. It should assume that the aligner
# is not present localy, but needs to be retrieved, unpacked, compiled and set-up, without requireing
# root privileges.
def download_and_install():
    if os.path.exists(ASSEMBLER_BIN):
        sys.stderr.write('[%s wrapper] Bin found at %s. Skipping installation ...\n' % (ASSEMBLER_NAME, ASSEMBLER_BIN))
    else:
        sys.stderr.write('[%s wrapper] Started installation of %s.\n' % (ASSEMBLER_NAME, ASSEMBLER_NAME))
        sys.stderr.write('[%s wrapper] Cloning git repository.\n' % (ASSEMBLER_NAME))
        command = 'cd %s; git clone %s' % (basicdefines.ASSEMBLERS_PATH_ROOT_ABS, ASSEMBLER_URL)
        sys.stderr.write('[%s wrapper] %s\n' % (ASSEMBLER_NAME, command))
        subprocess.call(command, shell='True')
        sys.stderr.write('\n')

        if os.path.exists(BAMTOOLS_BIN):
            sys.stderr.write('[%s wrapper] Bamtools already installed. Skipping ...\n' % (ASSEMBLER_NAME))
        else:
            # Setting up bamtools
            sys.stderr.write('[%s wrapper] Cloning bamtools git repository.\n' % (ASSEMBLER_NAME))
            command = 'cd %s; git clone %s' % (basicdefines.ASSEMBLERS_PATH_ROOT_ABS, BAMTOOLS_URL)
            sys.stderr.write('[%s wrapper] %s\n' % (ASSEMBLER_NAME, command))
            subprocess.call(command, shell='True')
            sys.stderr.write('\n')

            sys.stderr.write('[%s wrapper] Setting up bamtools.\n' % (ASSEMBLER_NAME))
            command = 'cd %s; mkdir build; cd build; cmake ..; make' % (BAMTOOLS_PATH)
            sys.stderr.write('[%s wrapper] %s\n' % (ASSEMBLER_NAME, command))
            subprocess.call(command, shell='True')
            sys.stderr.write('\n')

        # Making SGA
        sys.stderr.write('[%s wrapper] Running make.\n' % (ASSEMBLER_NAME))
        command = 'cd %s; ./autogen.sh' % (os.path.join(ASSEMBLER_PATH, 'src'))
        sys.stderr.write('[%s wrapper] %s\n' % (ASSEMBLER_NAME, command))
        subprocess.call(command, shell='True')
        sys.stderr.write('\n')
        command = 'cd %s; ./configure --with-bamtools=%s; make' % (ASSEMBLER_PATH, BAMTOOLS_PATH)
        sys.stderr.write('[%s wrapper] %s\n' % (ASSEMBLER_NAME, command))
        subprocess.call(command, shell='True')
        sys.stderr.write('\n')



def verbose_usage_and_exit():
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s mode [<reads_file> <machine_name> <output_path> [<output_suffix>]]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\t- mode - either "run" or "install". Is "install" other parameters can be ommitted.\n')

    exit(0)

if __name__ == "__main__":
    if (len(sys.argv) < 2 or len(sys.argv) > 7):
        verbose_usage_and_exit()

    if (sys.argv[1] == 'install'):
        download_and_install()
        exit(0)

    elif (sys.argv[1] == 'run'):
        if (len(sys.argv) < 6):
            verbose_usage_and_exit()

        reads_file = sys.argv[2]
        reference_file = sys.argv[3]
        machine_name = sys.argv[4]
        output_path = sys.argv[5]
        output_suffix = ''

        if (len(sys.argv) == 7):
            output_suffix = sys.argv[6]
        run(reads_file, reference_file, machine_name, output_path, output_suffix)

    else:
        verbose_usage_and_exit()
