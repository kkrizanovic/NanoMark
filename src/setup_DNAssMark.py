#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
import sys
sys.path.append(SCRIPT_PATH + '')
sys.path.append(SCRIPT_PATH + '/wrappers');
import subprocess

import basicdefines

QUAST_URL = 'http://sourceforge.net/projects/quast/files/quast-3.0.tar.gz'
QUAST_ZIP_FILE = 'quast-3.0.tar.gz'
QUAST_ZIP_PATH = os.path.join(basicdefines.TOOLS_ROOT_ABS, QUAST_ZIP_FILE)
QUAST_BIN_FILE = 'quast.py'
QUAST_BIN_PATH = os.path.join(basicdefines.TOOLS_ROOT_ABS, 'quast-3.0', QUAST_BIN_FILE)
CGMEMTIME_FILE = 'cgmemtime'
CGMEMTIME_PATH = os.path.join(basicdefines.TOOLS_ROOT_ABS, 'cgmemtime', CGMEMTIME_FILE)

def create_folders():
	sys.stderr.write('Generating folders...\n')

	if not os.path.exists(basicdefines.TOOLS_ROOT_ABS):
		sys.stderr.write('Creating folder "%s".\n' % basicdefines.TOOLS_ROOT_ABS)
		os.makedirs(basicdefines.TOOLS_ROOT_ABS)

	if not os.path.exists(basicdefines.INTERMEDIATE_PATH_ROOT_ABS):
		sys.stderr.write('Creating folder "%s".\n' % basicdefines.INTERMEDIATE_PATH_ROOT_ABS)
		os.makedirs(basicdefines.INTERMEDIATE_PATH_ROOT_ABS)

	sys.stderr.write('\n')


def download_assemblers():
	sys.stderr.write('Installing assembly tools.\n')

	assembler_wrappers = basicdefines.find_files(basicdefines.WRAPPERS_PATH_ROOT_ABS, 'wrapper_*.py')

	for wrapper in assembler_wrappers:
		wrapper_basename = os.path.splitext(os.path.basename(wrapper))[0]
		command = 'import %s; %s.download_and_install()' % (wrapper_basename, wrapper_basename)
		exec(command)


def setup_tools():
	if os.path.exists(CGMEMTIME_PATH):
		sys.stderr.write('Cgmemtime already installed. Skipping...\n')
	else:
		sys.stderr.write('Cloning Cgmemtime Git repo. Git needs to be installed.\n')
		command = 'cd %s; git clone https://github.com/isovic/cgmemtime.git' % (basicdefines.TOOLS_ROOT_ABS)
		subprocess.call(command, shell='True')
		sys.stderr.write('\n')

		# TODO: The following needs to be run to be able to use cgmemtime
		#       sudo ./cgmemtime --setup -g myusergroup --perm 775

	if os.path.exists(QUAST_BIN_PATH):
		sys.stderr.write('Quast already installed. Skipping ....\n')
	else:
		if not os.path.exists(QUAST_ZIP_PATH):
			sys.stderr.write('Downloading Quast 3.0 from Sourceforge. Need tar to decompress.\n')
			# Download
			command = 'cd %s; wget %s' % (basicdefines.TOOLS_ROOT_ABS, QUAST_URL)
			subprocess.call(command, shell='True')

		# Decompress
		command = 'cd %s; tar -xzf %s' % (basicdefines.TOOLS_ROOT_ABS, QUAST_ZIP_FILE)
		subprocess.call(command, shell='True')
		sys.stderr.write('\n')


def setup_all():
	create_folders()
	setup_tools()
	download_assemblers()



def verbose_usage_and_exit():
	sys.stderr.write('Usage:\n')
	sys.stderr.write('\t%s\n' % sys.argv[0])
	sys.stderr.write('\n')
	exit(0)

if __name__ == '__main__':
	setup_all()

	# if (len(sys.argv) > 2):
	# 	verbose_usage_and_exit()
