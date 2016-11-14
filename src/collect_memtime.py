#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

import sys
# sys.path.append(SCRIPT_PATH + '/../src/')

import subprocess
import multiprocessing
import fnmatch
from time import gmtime, strftime

from dataspec import *

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

def verbose_usage_and_exit():
	sys.stderr.write('Tool for collecting memtime results.\n');
	sys.stderr.write('Usage:\n');
#	sys.stderr.write('\t%s <memtime_folder> <out_total_memtime_file> memtime_prefix [last_memtime_name]\n' % (sys.argv[0]));
	sys.stderr.write('\t%s mode' % (sys.argv[0]));
	sys.stderr.write('mode - either "pattern" or "number". If pattern, all memtime files will be found using a given pattern. If number, a prefix is used in combination with a range of numbers spanning all file names.\n' % (sys.argv[0]));
#	sys.stderr.write('\t%s <memtime_folder> <out_total_memtime_file> memtime_prefix [last_memtime_name]\n' % (sys.argv[0]));
#	sys.stderr.write('last_memtime_name - the list of memtime files is first sorted. If there is a file specified by this parameter, all following files will be ignored.\n')
	sys.stderr.write('\n');
	exit(1);

if __name__ == "__main__":
    if (len(sys.argv) < 2):
        verbose_usage_and_exit();

#    if (len(sys.argv) != 4):
#        verbose_usage_and_exit()

    if (sys.argv[1] == 'pattern'):
        if (len(sys.argv) != 5):
            sys.stderr.write('Tool for collecting memtime results.\n');
	    sys.stderr.write('Usage:\n');
	    sys.stderr.write('\t%s mode <memtime_folder> <out_total_memtime_file> memtime_pattern [last_memtime_name]\n' % (sys.argv[0]));
	    sys.stderr.write('mode - either "pattern" or "number". If pattern, all memtime files will be found using a given pattern. If number, a prefix is used in combination with a range of numbers spanning all file names.\n' % (sys.argv[0]));
	    sys.stderr.write('\n');
    	    exit(1);

        memtime_folder = sys.argv[2];
        out_file = sys.argv[3];
        memtime_prefix = sys.argv[4];

        files = find_files(memtime_folder, memtime_prefix);
        print files;
        parse_memtime_files_and_accumulate(files, out_file);

    elif (sys.argv[1] == 'number'):
        if (len(sys.argv) != 7):
            sys.stderr.write('Tool for collecting memtime results.\n');
	    sys.stderr.write('Usage:\n');
	    sys.stderr.write('\t%s mode <memtime_folder> <out_total_memtime_file> memtime_prefix start_num end_num' % (sys.argv[0]));
	    sys.stderr.write('start_num - number of the first memtime file.\n');
	    sys.stderr.write('end_num - number of the last memtime file (inclusive).\n');
	    sys.stderr.write('mode - either "pattern" or "number". If pattern, all memtime files will be found using a given pattern. If number, a prefix is used in combination with a range of numbers spanning all file names.\n');
	    sys.stderr.write('\n');
    	    exit(1);

        memtime_folder = sys.argv[2];
        out_file = sys.argv[3];
        memtime_prefix = sys.argv[4];
        start_num = int(sys.argv[5]);
        end_num = int(sys.argv[6]) + 1;

#        files = find_files(memtime_folder, memtime_prefix);
        files = ['%s/%s-%s.memtime' % (memtime_folder, memtime_prefix, value) for value in xrange(start_num, end_num)];
        print files;
        parse_memtime_files_and_accumulate(files, out_file);
    else:
        verbose_usage_and_exit();
