 #! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

import sys;
import subprocess;

import fnmatch
from time import gmtime, strftime

WRAPPERS_PATH = 'wrappers'
WRAPPERS_PATH_ROOT_ABS = SCRIPT_PATH + '/../' + WRAPPERS_PATH

TOOLS_ROOT = 'tools'
TOOLS_ROOT_ABS = SCRIPT_PATH + '/../' + TOOLS_ROOT

ASSEMBLERS_PATH = 'assemblers'
ASSEMBLERS_PATH_ROOT_ABS = SCRIPT_PATH + '/../' + ASSEMBLERS_PATH

INTERMEDIATE_PATH = 'intermediate'
INTERMEDIATE_PATH_ROOT_ABS = SCRIPT_PATH + '/../' + INTERMEDIATE_PATH

RESULTS_PATH = 'results'
RESULTS_PATH_ROOT_ABS = SCRIPT_PATH + '/../' + RESULTS_PATH

QUAST_PATH = 'quast-3.1'
QUAST_PATH_ROOT_ABS = os.path.join(TOOLS_ROOT_ABS, QUAST_PATH)
QUAST_BIN = os.path.join(QUAST_PATH_ROOT_ABS, 'quast.py')

DNADIFF_PATH = 'quast-3.1/libs/MUMmer3.23-linux/';
DNADIFF_PATH_ROOT_ABS = os.path.join(TOOLS_ROOT_ABS, DNADIFF_PATH)
DNADIFF_BIN = os.path.join(DNADIFF_PATH_ROOT_ABS, 'dnadiff')

DATASETS_PATH = 'intermediate'
DATASETS_PATH_ROOT_ABS = SCRIPT_PATH + '/../' + DATASETS_PATH

CGMEMTIME_FILE = 'cgmemtime'
CGMEMTIME_PATH = os.path.join(TOOLS_ROOT_ABS, 'cgmemtime')
CGMEMTIME_PATH_ABS = SCRIPT_PATH + '/../' + os.path.join(TOOLS_ROOT_ABS, 'cgmemtime')

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
        return (SCRIPT_PATH + r'/../tools/cgmemtime/cgmemtime -o ' + measure_file)

### Logs messages to STDERR and an output log file if provided (opened elsewhere).
def log_message(message):
    timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());
    return '[%s] %s\n' % (timestamp, message);

def log(message, fp_log):
    timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());

    sys.stderr.write('[%s] %s\n' % (timestamp, message))
    sys.stderr.flush();
    if (fp_log != None):
        fp_log.write('[%s] %s\n' % (timestamp, message))
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

if __name__ == '__main__':
	pass
