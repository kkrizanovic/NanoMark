 #! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

import fnmatch

WRAPPERS_PATH = 'wrappers'
WRAPPERS_PATH_ROOT_ABS = SCRIPT_PATH + '/../src/' + WRAPPERS_PATH

TOOLS_ROOT = 'tools'
TOOLS_ROOT_ABS = SCRIPT_PATH + '/../' + TOOLS_ROOT

ASSEMBLERS_PATH = 'assemblers'
ASSEMBLERS_PATH_ROOT_ABS = SCRIPT_PATH + '/../' + ASSEMBLERS_PATH

INTERMEDIATE_PATH = 'intermediate'
INTERMEDIATE_PATH_ROOT_ABS = SCRIPT_PATH + '/../' + INTERMEDIATE_PATH

QUAST_PATH = 'quast-3.0'
QUAST_PATH_ROOT_ABS = os.path.join(TOOLS_ROOT_ABS, QUAST_PATH)
QUAST_BIN = os.path.join(QUAST_PATH_ROOT_ABS, 'quast.py')

TECH = ['illumina', 'pacbio', 'nanopore']



# Finds all files with a given filename pattern starting from the given path, recursivelly.
def find_files(start_path, file_pattern, max_depth=-1):
	matches = []
	for root, dirnames, filenames in os.walk(start_path):
		for filename in fnmatch.filter(filenames, file_pattern):
			depth = len(root.split('/'))
			if (max_depth < 0 or (max_depth >= 0 and depth <= max_depth)):
				matches.append(os.path.join(root, filename))

	matches.sort()
	return matches

# Simply lists all subfolders, recursivelly. Parameter 'depth' can be used to specify the minimum depth of the subfolders.
def find_folders(start_path, depth=0):
	matches = []
	for x in os.walk(start_path):
		if (len(x[0].split('/')) >= depth):
			matches.append(x[0])
	return matches

def measure_command(measure_file):
        return (SCRIPT_PATH + r'/../tools/cgmemtime/cgmemtime -o ' + measure_file + ' ')


if __name__ == '__main__':
	pass
