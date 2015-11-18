#! /usr/bin/python

import os
import sys
import re
import shutil
import datetime
import subprocess

import setup_nanomark
import basicdefines

# To generate random number for each run, so that multiple runs don't clash
import uuid

# To be able to import wrappers more easily
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(SCRIPT_PATH + '')
sys.path.append(SCRIPT_PATH + '/../wrappers');


# def benchmark_cf(configfile):
#     sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name))
#     exit(1)
#     pass


# A helper function that return all available wrappers from 'src/wrappers' folder
def get_allwrappers():
    wrappers = []
    wrapper_basenames = [];
    for filename in os.listdir(basicdefines.WRAPPERS_PATH_ROOT_ABS):
        if filename.startswith('wrapper_') and filename.endswith('.py'):
            wrappername = filename[8:-3]
            wrappers.append(wrappername)
            wrapper_basenames.append(os.path.splitext(os.path.basename(filename))[0]);
    return [wrappers, wrapper_basenames];


# Lists all available assemblers/wrappers
# All wrappers are stored in a separate folder in python files starting with 'wrapper_'
def list_wrappers():
    sys.stdout.write('\nAvailable wrappers:')
    [wrappers, wrapper_basenames] = get_allwrappers()
    for wrapper,wrapper_basenames in zip(wrappers, wrapper_basenames):
        wrapper_path = os.path.join(basicdefines.WRAPPERS_PATH_ROOT_ABS, wrapper + '.py')
        command = 'import %s; assembler_name = %s.ASSEMBLER_NAME; assembler_type = %s.ASSEMBLER_TYPE' % (wrapper_basenames, wrapper_basenames, wrapper_basenames)
        exec(command)
        padding = ' ' * (16 - len(assembler_name));
        sys.stdout.write('\n\t- %s%s- %s' % (assembler_name, padding, assembler_type))

    sys.stdout.write('\n\n')



def benchmark(reads_file, reference_file, technology, wrapper_list = []):

    uuid_string = str(uuid.uuid4())            # Generate a random UUID so that multiple runs don't clash
    output_folder = 'benchmark_' + uuid_string
    output_path = os.path.join(basicdefines.INTERMEDIATE_PATH_ROOT_ABS, output_folder)
    log_filename = os.path.join(output_path, 'log.txt')
    creffilename = 'REF_' + uuid_string + os.path.splitext(reference_file)[1]        # preserve extension
    crefpath = os.path.join(output_path, creffilename)
    creadsfilename = 'READS_' + uuid_string + os.path.splitext(reads_file)[1]        # preserve extension
    creadspath = os.path.join(output_path, creadsfilename)


    if (not os.path.exists(output_path)):
        sys.stderr.write('Creating folder "%s".\n' % (output_path))
        os.makedirs(output_path)

    with open(log_filename, 'w') as logfile:
        cDateTime = str(datetime.datetime.now())
        logfile.write('NanoMark benchmark UUID: %s\n' % uuid_string)
        logfile.write('Date and time: %s\n' % cDateTime)
        sys.stderr.write('\n')
        sys.stderr.write('Running benchmark for one reads file!\n')
        sys.stderr.write('UUID: %s\n' % uuid_string)
        sys.stderr.write('Reads file: %s\n' % reads_file)
        sys.stderr.write('References file: %s\n' % reference_file)

        # Copy reference file to output folder while preserving extension
        shutil.copy(reference_file, crefpath)
        logfile.write('Reference file: %s\n' % reference_file)

        # Copy reads file to output folder while preserving extension
        shutil.copy(reads_file, creadspath)
        logfile.write('Reads file: %s\n' % reads_file)

        # Create global Quast folder
        gl_quast_folder = os.path.join(output_path, 'quast')
        os.makedirs(gl_quast_folder)

        # If wrappers are not specified, get all wrappers
        if not wrapper_list:
            assembler_wrappers = basicdefines.find_files(basicdefines.WRAPPERS_PATH_ROOT_ABS, 'wrapper_*.py')
            for wrapper in assembler_wrappers:
                wrapper_basename = os.path.splitext(os.path.basename(wrapper))[0]       # Removing full path and extension
                wrapper_basename = wrapper_basename[8:]                                 # Removing 'wrapper_' from the start
                wrapper_list.append(wrapper_basename)

        # Checking if all wrappers in wrapper_list actually exist
        for wrapper_basename in wrapper_list:
            wrapper = 'wrapper_' + wrapper_basename
            wrapper_path = os.path.join(basicdefines.WRAPPERS_PATH_ROOT_ABS, wrapper + '.py')
            if not os.path.exists(wrapper_path):
                sys.stderr.write('\nWrapper %s does not exist! It will not be used for the benchmark!\n' % wrapper_basename)
                wrapper_list.remove(wrapper_basename)

        # Writting wrappers to be used to the log file
        # THis will anable an interrupted run to be continued
        logfile.write('Wrappers used: %s' % ','.join(wrapper_list))

        for wrapper_basename in wrapper_list:
            assembler_name = ''
            results_file = ''
            create_output_folder = True
            wrapper = 'wrapper_' + wrapper_basename
            wrapper_path = os.path.join(basicdefines.WRAPPERS_PATH_ROOT_ABS, wrapper + '.py')
            command = 'import %s; assembler_name = %s.ASSEMBLER_NAME; create_output_folder = %s.CREATE_OUTPUT_FOLDER; assembler_type = %s.ASSEMBLER_TYPE' % (wrapper, wrapper, wrapper, wrapper, wrapper)
            sys.stderr.write('Executing command: %s\n' % (command));
            exec(command)
            # Run assembly only using non-hybrid assemblers. Using hybrid assemblers complicates uniform specification of the datasets.
            if (assembler_type != 'nonhybrid'):
                continue;
            logfile.write('%s (%s): %s\n' % (assembler_name, wrapper, wrapper_path))

            # Create folder for each assembler's results, and for quast results
            assembler_folder = os.path.join(output_path, assembler_name)
            loc_quast_folder = os.path.join(gl_quast_folder, assembler_name)
            # # Ray creates his own output folder (doesn't work if its already created!)
            # if create_output_folder:
            #     os.makedirs(assembler_folder)
            # os.makedirs(loc_quast_folder)
            # Run each wrapper
            # def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):
            sys.stderr.write('\n\nRunning assembler %s\n' % assembler_name)
            command = 'import %s; %s.run(\'%s\', \'%s\', \'%s\', \'%s\')' % (wrapper, wrapper, creadspath, crefpath, technology, assembler_folder)
            sys.stderr.write('Executing command: %s\n' % (command));
#            exec(command)

            # # Run quast on results file?
            # # This might be a part of a wrapper implementation
            # # however it seems more logical to place come here instead of in each wrapper
            # results_path = os.path.join(assembler_folder, results_file)
            # if os.path.exists(results_path):
            #     sys.stderr.write('Running quast for assembler %s\n' % assembler_name)
            #     command = '%s %s -R %s -o %s' % (basicdefines.QUAST_BIN, results_path, crefpath, loc_quast_folder)
            #     subprocess.call(command, shell='True')
            #     logfile.write('OK\n')
            # else:
            #     sys.stderr.write('Cannot find results file %s for assembler %s\n' % (results_file, assembler_name))
            #     logfile.write('ERROR\n')

    # Collect quast and cgmemtime results
    # summarize_results(output_path)


# Continue benchmark that was interrupted
# Takes benchmark results folder as an argument
def continueBenchmark(results_folder):
    sys.stderr.write('\n\nContinuing benchmark in folder: %s\n' % results_folder)
    # First check if log file exists
    log_filename = os.path.join(results_folder, 'log.txt')
    if not os.path.exists(log_filename):
        sys.stderr.write('\n\nCannot find log file (log.txt) within benchmark results folder (%s)! Exiting ...\n' % results_folder)
        verbose_usage_and_exit()

    # Global quast folder
    gl_quast_folder = os.path.join(output_path, 'quast')
    if not os.path.exists(gl_quast_folder):
        sys.stderr.write('Creating global QUAST folder for the benchmark...\n')
        os.makedirs(gl_quast_folder)

    sys.stderr.write('Reading log file and checking wrapers....\n')
    # Load wrapper names from log file, log file is only used to see which wrappers were supposed to be used
    # in a bechmark. Which wrappers completed their run is check by looking if a corresponding results file exists.
    # This could also be done read from the logfile.
    # Check if a results file exists for each wrapper. If results file doesn't exist, run the assembler again
    with open(log_filename, 'r') as logfile:
        # Skip 4 lines
        logfile.readline();logfile.readline();logfile.readline();logfile.readline()
        # Get a line with wrapper names
        line = logfile.readline()
        sys.stderr.write(line)
        line = line[15:-1]          # Remove 'Wrappers used: ' from the start and \n from the end
        wrapper_list = line.split(',')

        for wrapper_basename in wrapper_list:
            assembler_name = ''
            results_file = ''
            create_output_folder = True
            wrapper = 'wrapper_' + wrapper_basename
            wrapper_path = os.path.join(basicdefines.WRAPPERS_PATH_ROOT_ABS, wrapper + '.py')
            command = 'import %s; assembler_name = %s.ASSEMBLER_NAME; create_output_folder = %s.CREATE_OUTPUT_FOLDER; assembler_type = %s.ASSEMBLER_TYPE' % (wrapper, wrapper, wrapper, wrapper, wrapper)
            exec(command)
            # Run assembly only using non-hybrid assemblers. Using hybrid assemblers complicates uniform specification of the datasets.
            if (assembler_type != 'nonhybrid'):
                continue;
            assembler_folder = os.path.join(output_path, assembler_name)
            results_path = os.path.join(assembler_folder, results_file)
            loc_quast_folder = os.path.join(gl_quast_folder, assembler_name)

            # Get reads and reference files, they should be inside benchmark folder
            # With names that start with READS_ and REF_

            if os.path.exists(results_path):
                # The results file exists, thi means that the assembler run completed and will not be repeated
                sys.stderr.write('Assembler %s run completed.\n')

            # # If results file doesn't exist, the assembler didn't complete its run and needs to repeat it
            # # First remove whole assembler folder if it exists, and then recreate it
            # if os.path.exists(assembler_folder):
            #     os.rmdir(assembler_folder)
            # # Create output folder if specified in the wrapper
            # if create_output_folder:
            #     os.makedirs(assembler_folder)
            # # Remove local QUAST folder if it exists, and then recreate it
            # if os.path.exists(loc_quast_folder):
            #     os.rmdir(loc_quast_folder)
            # os.makedirs(loc_quast_folder)


        # Run each wrapper
        # def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):
        sys.stderr.write('\n\nRunning assembler %s\n' % assembler_name)
        command = 'import %s; %s.run(\'%s\', \'%s\', \'machine_name\', \'%s\')' % (wrapper, wrapper, creadspath, crefpath, assembler_folder)
        sys.stderr.write('Executing command: %s\n' % (command));
    #            exec(command)

        # # Run quast on results file?
        # # This might be a part of a wrapper implementation
        # # however it seems more logical to place come here instead of in each wrapper
        # results_path = os.path.join(assembler_folder, results_file)
        # if os.path.exists(results_path):
        #     sys.stderr.write('Running quast for assembler %s\n' % assembler_name)
        #     command = '%s %s -R %s -o %s' % (basicdefines.QUAST_BIN, results_path, crefpath, loc_quast_folder)
        #     subprocess.call(command, shell='True')
        #     logfile.write('OK\n')
        # else:
        #     sys.stderr.write('Cannot find results file %s for assembler %s\n' % (results_file, assembler_name))
        #     logfile.write('ERROR\n')

        # # Collect quast and cgmemtime results
        # summarize_results(output_path)



# def continueBenchmark_cf(results_folder):
#     sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name))
#     exit(1)
#     pass


# def run_quast(results_folder):
#     sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name))
#     exit(1)
#     pass



def summarize_results(results_folder):

    # get global quast folder and log file
    gl_quast_folder = os.path.join(results_folder, 'quast')
    log_filename = os.path.join(results_folder, 'log.txt')

    # Get all wrappers
    assembler_wrappers = basicdefines.find_files(basicdefines.WRAPPERS_PATH_ROOT_ABS, 'wrapper_*.py')

    # Collect quast results
    summaryfilepath = os.path.join(results_folder, 'benchmark_summary.tsv')
    report_filename = 'report.txt'

    # Names of quast fields to write in summary tsv file
    qfields = ['# contigs', 'Largest contig', 'Total length', 'GC (%)', 'N50', 'NG50', '# misassemblies', 'Genome fraction (%)'
             , 'Duplication ratio', '# N\'s per 100 kbp', '# mismatches per 100 kbp', '# indels per 100 kbp']
    title_line = 'Name\tReal time\tCPU time\tMaximum RSS\t' + '\t'.join(qfields)
    with open(summaryfilepath, 'w') as summaryfile, open(log_filename, 'a') as logfile:
        summaryfile.write(title_line + '\n')
        logfile.write('Gathering results:\n')

        for wrapper in assembler_wrappers:
            assembler_name = ''
            wrapper_basename = os.path.splitext(os.path.basename(wrapper))[0]
            command = 'import %s; assembler_name = %s.ASSEMBLER_NAME' % (wrapper_basename, wrapper_basename)
            exec(command)
            assembler_folder = os.path.join(results_folder, assembler_name)
            memtime_path = os.path.join(assembler_folder, assembler_name + '.memtime')

            # read local .memtime and quast report.txt file
            loc_quast_folder = os.path.join(gl_quast_folder, assembler_name)
            report_path = os.path.join(loc_quast_folder, report_filename)
            reportdict = {}
            line = assembler_name
            real_time = cpu_time = user_time = system_time = '0.0'
            max_rss = '0'
            # Reading .memtime file
            if os.path.exists(memtime_path):
                with open(memtime_path, 'r') as memtime_file:
                    memtime_file.readline()                # skipping 1st line
                    tline = memtime_file.readline()
                    real_time = tline.split()[2]     # real time
                    tline = memtime_file.readline()
                    cpu_time = tline.split()[2]      # cpu time
                    tline = memtime_file.readline()
                    user_time = tline.split()[2]     # user time
                    tline = memtime_file.readline()
                    system_time = tline.split()[2]   # system time
                    tline = memtime_file.readline()
                    max_rss = tline.split()[2]       # Max RSS
                logfile.write(assembler_name + ': CGMemtime report processed!\n')
            else:
                logfile.write(assembler_name + ': CGMemtime report NOT found!\n')

            line += '\t%s s\t%s s\t%s MB' % (real_time, cpu_time, max_rss)

            # Reading quast report.txt file
            if os.path.exists(report_path):
                with open(report_path, 'r') as reportfile:
                    #skipp first three lines
                    for i in range(3):
                        reportfile.readline()
                    # Read everything else into a dictionary
                    for tline in reportfile:
                        elements = tline.split()
                        key = ' '.join(elements[:-1])
                        value = elements[-1]
                        reportdict[key] = value
                    for name in qfields:
                        line += '\t' + reportdict.get(name, 'None')
                logfile.write(assembler_name + ': Quast report processed!\n')
            else:
                logfile.write(assembler_name + ': Quast report NOT found!\n')

            # Writting a line to summary file
            summaryfile.write(line + '\n')

def verbose_usage_and_exit():
    sys.stderr.write('NanoMark - a DNA assembly benchmarking tool.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    sys.stderr.write('\t\tsetup\n')
    sys.stderr.write('\t\tbenchmark\n')
    sys.stderr.write('\t\tcontinue\n')
    # sys.stderr.write('\t\tbenchmark_cf\n')
    # sys.stderr.write('\t\tcontinue_cf\n')
    # sys.stderr.write('\t\trun_quast\n')
    sys.stderr.write('\t\tsummarize_results\n')
    sys.stderr.write('\t\tlist\n')
    sys.stderr.write('\n')
    exit(0)

def main():
    if (len(sys.argv) < 2):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if mode == 'setup':
        if (len(sys.argv) != 2):
            sys.stderr.write('Setup the folder structures and install necessary tools.\n')
            sys.stderr.write('Requires no additional parameters to run.\n')
            sys.stderr.write('\n')
            exit(1)

        setup_nanomark.setup_all()

    elif mode == 'benchmark':
        if (len(sys.argv) < 5 or len(sys.argv) > 7):
            sys.stderr.write('Runs a benchmark on given reads and reference files.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('%s %s <reads_file> <reference_file> <technology> options' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            sys.stderr.write('<technology> can be one of: \n')
            sys.stderr.write('\t%s\n' % ', '.join(basicdefines.TECH))
            sys.stderr.write('options:"\n')
            sys.stderr.write('-as <assemblers> - comma separated list of assemblers to be used"\n')
            sys.stderr.write('                   wrapper files are of the form: wrapper_[aligner name].py"\n')
            sys.stderr.write('                   available wrappers can be obtained using option: list\n')
            sys.stderr.write('                   if not specified, all available assemblers are used"\n')
            sys.stderr.write('\n')
            exit(1)

        reads_file = sys.argv[2]
        reference_file = sys.argv[3]
        technology = sys.argv[4].lower()

        if technology not in basicdefines.TECH:
            sys.stderr.write('\nInvalid Technology!')
            sys.stderr.write('\nTechnology must be one of:')
            sys.stderr.write('\t%s\n' % ', '.join(basicdefines.TECH))
            verbose_usage_and_exit()

        used_assemblers = []
        [all_assemblers, all_assemblers_basenames] = get_allwrappers()

        for i in range(5, len(sys.argv), 2):
            if sys.argv[i] == '-as':
                assemblersstring = sys.argv[i+1]
                assemblers = assemblersstring.split(',')
                for assembler in assemblers:
                    if assembler not in all_assemblers:
                        sys.stderr.write('\nUnavailable assembler: %s\n\n' % assembler)
                        exit(1)
                    else:
                        used_assemblers.append(assembler)
            else:
                sys.stderr.write('Invalid parameter!.\n\n')
                exit()

        # If no aligners are specified, use all of them
        if len(used_assemblers) == 0:
            used_assemblers = all_assemblers

        benchmark(reads_file, reference_file, technology, used_assemblers)

    elif mode == 'continue':
        if (len(sys.argv) != 3):
            sys.stderr.write('Continues a benchmark on given reads and reference files.\n')
            sys.stderr.write('Reads a corresponding log file to see what is left to be done.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <results_folder>' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        results_folder = sys.argv[2]
        continueBenchmark(results_folder)


    # elif mode == 'benchmark_cf':
    #     if (len(sys.argv) != 3):
    #         sys.stderr.write('Runs a benchmark according to a given config file.\n')
    #         sys.stderr.write('Config file contains a list of assemblers to be used, reads file and reference file.\n')
    #         sys.stderr.write('\n')
    #         sys.stderr.write('Usage:\n')
    #         sys.stderr.write('\t%s %s <config_file>' % (sys.argv[0], sys.argv[1]))
    #         sys.stderr.write('\n')
    #         exit(1)


    #     configfile = sys.argv[2]
    #     benchmark_cf(configfile)

    # elif mode == 'continue_cf':
    #     if (len(sys.argv) != 3):
    #         sys.stderr.write('Continues a benchmark according to a given config file.\n')
    #         sys.stderr.write('Config file contains a list of assemblers to be used, reads file and reference file.\n')
    #         sys.stderr.write('Reads a corresponding log file to see what is left to be done.\n')
    #         sys.stderr.write('\n')
    #         sys.stderr.write('Usage:\n')
    #         sys.stderr.write('\t%s %s <results_folder>' % (sys.argv[0], sys.argv[1]))
    #         sys.stderr.write('\n')
    #         exit(1)

    #     results_folder = sys.argv[2]
    #     continueBenchmark_cf(reads_folder)

    # elif mode == 'run_quast':
    #     if (len(sys.argv) != 3):
    #         sys.stderr.write('Runs quast on a set of assembler results.\n')
    #         sys.stderr.write('\n')
    #         sys.stderr.write('Usage:\n')
    #         sys.stderr.write('\t%s %s <results_folder>' % (sys.argv[0], sys.argv[1]))
    #         sys.stderr.write('\n')
    #         exit(1)

    #     results_folder = sys.argv[2]
    #     run_quast(results_folder)

    elif mode == 'summarize_results':
        if (len(sys.argv) != 3):
            sys.stderr.write('Summarizes results from quast and cgmem timefiles.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <results_folder>' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        results_folder = sys.argv[2]
        summarize_results(results_folder)

    elif (mode == 'list'):
        if (len(sys.argv) != 2):
            sys.stderr.write('Lists available wrappers.\n')
            sys.stderr.write('Requires no additional parameters to run.\n')
            sys.stderr.write('\n')
            exit(1)

        list_wrappers()

    else:
        print 'Invalid option!'



# Standard template for the main() function
if __name__ == '__main__':
  main()
