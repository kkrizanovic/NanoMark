#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(SCRIPT_PATH + '/../src')

import subprocess
import multiprocessing

import basicdefines

ASSEMBLER_URL = 'http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.3/wgs-8.3rc2-Linux_amd64.tar.bz2'
ASSEMBLER_PATH = os.path.join(basicdefines.ASSEMBLERS_PATH_ROOT_ABS, 'wgs-8.3rc2')
ZIP_FILE = 'wgs-8.3rc2-Linux_amd64.tar.bz2'
ZIP_PATH = os.path.join(basicdefines.ASSEMBLERS_PATH_ROOT_ABS, ZIP_FILE)
ASSEMBLER_BIN = os.path.join(ASSEMBLER_PATH,'Linux-amd64/bin/runCA')
ASSEMBLER_ECBIN = os.path.join(ASSEMBLER_PATH,'Linux-amd64/bin/PBcR')
ASSEMBLER_NAME = 'Celera'
ASSEMBLER_RESULTS = 'out/9-terminator/asm.ctg.fasta'
CREATE_OUTPUT_FOLDER = True

PACBIO_SPEC = os.path.join(basicdefines.ASSEMBLERS_PATH_ROOT_ABS, 'wgs_pacbio.spec')
OXFORD_SPEC = os.path.join(basicdefines.ASSEMBLERS_PATH_ROOT_ABS, 'wgs_oxford.spec')

DO_ERROR_CORRECTION = True

# Run parameters
P_LENGTH = 500
P_PARTITIONS = 200


# Function 'run' should provide a standard interface for running a mapper. Given input parameters, it should run the
# alignment process, and convert any custom output results to the SAM format. Function should return a string with the
# path to the output file.
#    reads_file            Path to a FASTA/FASTQ file containing reads.
#    reference_file        Path to a reference genome FASTA file.
#    machine_name        A symbolic name to specify a set of parameters for a specific sequencing platform.
#    output_path            Folder to which the output will be placed to. Filename will be automatically generated according to the name of the mapper being run.
#    output_suffix        A custom suffix that can be added to the output filename.
def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):

    # TODO: Finish this
    num_threads = multiprocessing.cpu_count() / 2

    used_bin = ''
    if DO_ERROR_CORRECTION:
        used_bin = ASSEMBLER_ECBIN
    else:
        used_bin = ASSEMBLER_BIN
    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '.memtime')

    spec_file = ''
    if machine_name == 'pacbio':
        spec_file = PACBIO_SPEC
        command = 'cd %s; %s %s -pbCNS -length 500 -partitions 200 -l %s -s %s -fastq %s ' % (output_path, basicdefines.measure_command(memtime_path), used_bin, ASSEMBLER_NAME, spec_file, reads_file)
        subprocess.call(command, shell='True')
    elif machine_name == 'nanopore':
        spec_file = OXFORD_SPEC
        command = 'cd %s; %s %s -length 500 -partitions 200 -l %s -s %s -fastq %s ' % (output_path, basicdefines.measure_command(memtime_path), used_bin, ASSEMBLER_NAME, spec_file, reads_file)
        subprocess.call(command, shell='True')
    elif machine_name == 'illumina':
        sys.stderr.write('\n\nmachine_name \'illumina\' not implemented for assembler %s' % ASSEMBLER_NAME)
        sys.stderr.write('\nSkipping ....')
    else:
        sys.stderr.write('Invalid machine_name parameter for assembler %s' % ASSEMBLER_NAME)
        sys.stderr.write('\nSkipping ....')

    # Atm, quast is run in the main program


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


def verbose_usage_and_exit():
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s mode [<reads_file> <reference_file> <machine_name> <output_path> [<output_suffix>]]\n' % sys.argv[0])
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
