#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

import sys;
sys.path.append(SCRIPT_PATH + '/../src');

import subprocess;
import multiprocessing;

import basicdefines;
import fastqparser;

ALIGNER_URL = 'http://last.cbrc.jp/last-534.zip';

# ALIGNER_PATH = SCRIPT_PATH + '/../aligners/last-534/src';
ALIGNER_PATH = basicdefines.ALIGNERS_PATH_ROOT_ABS + '/last-534/src';
BIN = 'lastal';
MAPPER_NAME = 'LAST';



# LAST's MAF to SAM conversion script doesn't include the SAM header lines which are needed for
# some downstream analyses. This function extracts the reference and formats it to SAM format.
def get_sam_header(reference_file):
	[headers, seqs, quals] = fastqparser.read_fastq(reference_file);

	line = '';

	i = 0;
	while i < len(headers):
		line += '@SQ\tSN:%s\tLN:%d\n' % (headers[i], len(seqs[i]));
		i += 1;

	return line;

# Function 'run' should provide a standard interface for running a mapper. Given input parameters, it should run the
# alignment process, and convert any custom output results to the SAM format. Function should return a string with the
# path to the output file.
#	reads_file			Path to a FASTA/FASTQ file containing reads.
#	reference_file		Path to a reference genome FASTA file.
#	machine_name		A symbolic name to specify a set of parameters for a specific sequencing platform.
#	output_path			Folder to which the output will be placed to. Filename will be automatically generated according to the name of the mapper being run.
#	output_suffix		A custom suffix that can be added to the output filename.
def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):
	parameters = '';
	num_threads = multiprocessing.cpu_count();

	if ((machine_name.lower() == 'illumina') or (machine_name.lower() == 'roche')):
		parameters = ' ';

	elif ((machine_name.lower() == 'pacbio')):
		parameters = '-q 1 -r 1 -a 1 -b 1';

	elif ((machine_name.lower() == 'nanopore')):
		parameters = '-q 1 -r 1 -a 1 -b 1';

	elif ((machine_name.lower() == 'debug')):
		parameters = ' ';

	else:			# default
		parameters = ' ';

	# Running in overlap mode, trying to minimize clipping
	parameters = parameters + ' -T 1'



	if (output_suffix != ''):
		output_filename = '%s-%s' % (MAPPER_NAME, output_suffix);
	else:
		output_filename = MAPPER_NAME;



	reads_fasta = reads_file;
	reads_basename = os.path.splitext(os.path.basename(reads_file))[0];
	maf_file = '%s/%s.maf' % (output_path, output_filename);
	sam_file = '%s/%s.sam' % (output_path, output_filename);
	memtime_file = '%s/%s.memtime' % (output_path, output_filename);
	memtime_file_index = '%s/%s-index.memtime' % (output_path, output_filename);
	memtime_file_maftosam = '%s/%s-maftosam.memtime' % (output_path, output_filename);
	reference_db_file = reference_file + '.db';

	# Check if the given input file is a FASTA or FASTQ, and convert to FASTA if necessary.
	if (reads_file[-1] == 'q'):
		sys.stderr.write('[%s wrapper] Converting FASTQ to FASTA...\n' % (MAPPER_NAME));
		reads_fasta = reads_file[0:-1] + 'a';
		fastqparser.convert_to_fasta(reads_file, reads_fasta);
		sys.stderr.write('\n');

	if not os.path.exists(reference_db_file + '.suf'):
		# Run the indexing process, and measure execution time and memory.
		sys.stderr.write('[%s wrapper] Generating index...\n' % (MAPPER_NAME));
		# command = '%s %s/lastdb %s %s' % (basicdefines.measure_command(memtime_file_index), ALIGNER_PATH, reference_db_file, reference_file);
		command = '%s/lastdb %s %s' % (ALIGNER_PATH, reference_db_file, reference_file);
		sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
		subprocess.call(command, shell=True);
		sys.stderr.write('\n\n');
	else:
		sys.stderr.write('[%s wrapper] Reference index already exists. Continuing.\n' % (MAPPER_NAME));

	# Run the alignment process, and measure execution time and memory.
	sys.stderr.write('[%s wrapper] Running %s...\n' % (MAPPER_NAME, MAPPER_NAME));
	# command = '%s %s/%s %s %s %s > %s' % (basicdefines.measure_command(memtime_file), ALIGNER_PATH, BIN, parameters, reference_db_file, reads_fasta, maf_file);
	command = '%s/%s %s %s %s > %s' % (ALIGNER_PATH, BIN, parameters, reference_db_file, reads_fasta, maf_file);
	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n\n');

	# Run the alignment process, and measure execution time and memory.
	sys.stderr.write('[%s wrapper] Converting the output MAF to SAM file...\n' % (MAPPER_NAME));
	fp = open(sam_file, 'w');
	fp.write(get_sam_header(reference_file));
	fp.close();
	# command = '%s %s/../scripts/maf-convert sam %s >> %s' % (basicdefines.measure_command(memtime_file_maftosam), ALIGNER_PATH, maf_file, sam_file);
	command = '%s/../scripts/maf-convert sam %s >> %s' % (ALIGNER_PATH, maf_file, sam_file);
	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n\n');

	sys.stderr.write('[%s wrapper] %s wrapper script finished processing.\n' % (MAPPER_NAME, MAPPER_NAME));

	return sam_file


# This is a standard interface for setting up the aligner. It should assume that the aligner
# is not present localy, but needs to be retrieved, unpacked, compiled and set-up, without requireing
# root privileges.
def download_and_install():
	sys.stderr.write('[%s wrapper] Started installation of %s.\n' % (MAPPER_NAME, MAPPER_NAME));
	sys.stderr.write('[%s wrapper] Using wget to download the aligner.\n' % (MAPPER_NAME));
	command = 'cd %s; wget %s; unzip %s' % (basicdefines.ALIGNERS_PATH_ROOT_ABS, ALIGNER_URL, os.path.basename(ALIGNER_URL));
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	sys.stderr.write('[%s wrapper] Running make.\n' % (MAPPER_NAME));
	command = 'cd %s/..; make' % (ALIGNER_PATH);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');



def verbose_usage_and_exit():
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s mode [<reads_file> <reference_file> <machine_name> <output_path> [<output_suffix>]]\n' % sys.argv[0]);
	sys.stderr.write('\n');
	sys.stderr.write('\t- mode - either "run" or "install". Is "install" other parameters can be ommitted.\n');

	exit(0);

if __name__ == "__main__":
	if (len(sys.argv) < 2 or len(sys.argv) > 7):
		verbose_usage_and_exit();

	if (sys.argv[1] == 'install'):
		download_and_install();
		exit(0);

	elif (sys.argv[1] == 'run'):
		if (len(sys.argv) < 6):
			verbose_usage_and_exit();

		reads_file = sys.argv[2];
		reference_file = sys.argv[3];
		machine_name = sys.argv[4];
		output_path = sys.argv[5];
		output_suffix = '';

		if (len(sys.argv) == 7):
			output_suffix = sys.argv[6];
		run(reads_file, reference_file, machine_name, output_path, output_suffix);

	else:
		verbose_usage_and_exit();
