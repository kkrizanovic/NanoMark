#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(SCRIPT_PATH + '/../src')

import subprocess
import multiprocessing

import basicdefines

ASSEMBLER_URL = 'http://sourceforge.net/projects/sparseassembler/files/SparseAssembler.zip'
ASSEMBLER_PATH = os.path.join(basicdefines.ASSEMBLERS_PATH_ROOT_ABS, 'SparseAssembler')
ZIP_FILE = 'SparseAssembler.zip'
ZIP_PATH = os.path.join(basicdefines.ASSEMBLERS_PATH_ROOT_ABS, ZIP_FILE)
ASSEMBLER_BIN = os.path.join(ASSEMBLER_PATH,'Compiled_Linux/SparseAssembler')
ASSEMBLER_NAME = 'SPARSE'
ASSEMBLER_RESULTS = 'Contigs.txt'
CREATE_OUTPUT_FOLDER = True



# Function 'run' should provide a standard interface for running a mapper. Given input parameters, it should run the
# alignment process, and convert any custom output results to the SAM format. Function should return a string with the
# path to the output file.
#    reads_file            Path to a FASTA/FASTQ file containing reads.
#    reference_file        Path to a reference genome FASTA file.
#    machine_name        A symbolic name to specify a set of parameters for a specific sequencing platform.
#    output_path            Folder to which the output will be placed to. Filename will be automatically generated according to the name of the mapper being run.
#    output_suffix        A custom suffix that can be added to the output filename.
def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):
    # Sparse also runs only on fasta
    # Atm parameters are hardcoded.
    # TODO: if fastq is given convert it to fasta
    #       callculate estimated genome size (GS) from reference and/or reads files
    num_threads = multiprocessing.cpu_count() / 2

    memtime_path = os.path.join(output_path, ASSEMBLER_NAME + '.memtime')
    command = 'cd %s; %s %s -t %d k 21 GS 60000000 f %s' % (output_path, basicdefines.measure_command(memtime_path), ASSEMBLER_BIN, num_threads, reads_file)
    subprocess.call(command, shell='True')

    # Atm, quast is run in the main program


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
        sys.stderr.write('[%s wrapper] Bin found at %s. Skipping installation ...\n' % (ASSEMBLER_NAME, ASSEMBLER_BIN))
    else:
        sys.stderr.write('[%s wrapper] Started installation of %s.\n' % (ASSEMBLER_NAME, ASSEMBLER_NAME))
        if not os.path.exists(ZIP_PATH):
            sys.stderr.write('[%s wrapper] Downloading zip...\n' % (ASSEMBLER_NAME))
            command = 'cd %s; wget %s' % (basicdefines.ASSEMBLERS_PATH_ROOT_ABS, ASSEMBLER_URL)
            sys.stderr.write('[%s wrapper] %s\n' % (ASSEMBLER_NAME, command))
            subprocess.call(command, shell='True')

        # Decompress
        command = 'cd %s; unzip %s' % (basicdefines.ASSEMBLERS_PATH_ROOT_ABS, ZIP_FILE)
        sys.stderr.write('[%s wrapper] %s\n' % (ASSEMBLER_NAME, command))
        subprocess.call(command, shell='True')

        # Already compiled, no need to make
        # Adding execute permissions to bin file
        command = 'cd %s; chmod 777 %s' % (basicdefines.ASSEMBLERS_PATH_ROOT_ABS, ASSEMBLER_BIN)
        sys.stderr.write('[%s wrapper] %s\n' % (ASSEMBLER_NAME, command))
        subprocess.call(command, shell='True')

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
