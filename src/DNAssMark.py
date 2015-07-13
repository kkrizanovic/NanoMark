
import os
import sys
import re
import shutil
import datetime
import subprocess

import setup_DNAssMark
import basicdefines

# To generate random number for each run, so that multiple runs don't clash
import uuid

# To be able to import wrappers more easily
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(SCRIPT_PATH + '')
sys.path.append(SCRIPT_PATH + '/wrappers');

# TODO: Continue a benchmark from a place where it was interrupted

def benchmark_cf(configfile):
    sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name))
    exit(1)
    pass


def benchmark(reads_file, reference_file):

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
        logfile.write('DNAssMark benchmark UUID: %s\n' % uuid_string)
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

        # Get all wrappers
        logfile.write('Wrappers used:\n')
        assembler_wrappers = basicdefines.find_files(basicdefines.WRAPPERS_PATH_ROOT_ABS, 'wrapper_*.py')

        for wrapper in assembler_wrappers:
            assembler_name = ''
            results_file = ''
            create_output_folder = True
            wrapper_basename = os.path.splitext(os.path.basename(wrapper))[0]
            command = 'import %s; assembler_name = %s.ASSEMBLER_NAME; results_file = %s.ASSEMBLER_RESULTS; create_output_folder = %s.CREATE_OUTPUT_FOLDER' % (wrapper_basename, wrapper_basename, wrapper_basename, wrapper_basename)
            exec(command)
            logfile.write('%s (%s): %s\n' % (assembler_name, wrapper_basename, wrapper))

            # Create folder for each assembler's results, and for quast results
            assembler_folder = os.path.join(output_path, assembler_name)
            loc_quast_folder = os.path.join(gl_quast_folder, assembler_name)
            # Ray creates his own output folder (doesn't work if its already created!)
            if create_output_folder:
                os.makedirs(assembler_folder)
            os.makedirs(loc_quast_folder)
            # Run each wrapper
            # def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):
            sys.stderr.write('\n\nRunning assembler %s\n' % assembler_name)
            command = 'import %s; %s.run(\'%s\', \'%s\', \'machine_name\', \'%s\')' % (wrapper_basename, wrapper_basename, creadspath, crefpath, assembler_folder)
            exec(command)

            # Run quast on results file?
            # This might be a part of a wrapper implementation
            # however it seems more logical to place come here instead of in each wrapper
            results_path = os.path.join(assembler_folder, results_file)
            if os.path.exists(results_path):
                sys.stderr.write('Running quast for assembler %s\n' % assembler_name)
                command = '%s %s -R %s -o %s' % (basicdefines.QUAST_BIN, results_path, crefpath, loc_quast_folder)
                subprocess.call(command, shell='True')
                logfile.write('OK\n')
            else:
                sys.stderr.write('Cannot find results file %s for assembler %s\n' % (results_file, assembler_name))
                logfile.write('ERROR\n')

    # Collect quast and cgmemtime results
    summarize_results(output_path)



def continueBenchmark(results_folder):
    sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name))
    exit(1)
    pass


def continueBenchmark_cf(results_folder):
    sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name))
    exit(1)
    pass


def run_quast(results_folder):
    sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name))
    exit(1)
    pass



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
    sys.stderr.write('DNAssMark - a DNA assembly benchmarking tool.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    sys.stderr.write('\t\tsetup\n')
    sys.stderr.write('\t\tbenchmark\n')
    sys.stderr.write('\t\tcontinue\n')
    sys.stderr.write('\t\tbenchmark_cf\n')
    sys.stderr.write('\t\tcontinue_cf\n')
    sys.stderr.write('\t\trun_quast\n')
    sys.stderr.write('\t\tsummarize_results\n')
    sys.stderr.write('\n')
    exit(0)

def main():
    if (len(sys.argv) < 2):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if mode == '--setup':
        if (len(sys.argv) != 2):
            sys.stderr.write('Setup the folder structures and install necessary tools.\n')
            sys.stderr.write('Requires no additional parameters to run.\n')
            sys.stderr.write('\n')
            exit(1)

        setup_DNAssMark.setup_all()

    elif mode == '--benchmark':
        if (len(sys.argv) != 4):
            sys.stderr.write('Runs a benchmark on given reads and reference files.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <reads_file> <reference_file>' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        reads_file = sys.argv[2]
        reference_file = sys.argv[3]
        benchmark(reads_file, reference_file)

    elif mode == '--continue':
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


    elif mode == '--benchmark_cf':
        if (len(sys.argv) != 3):
            sys.stderr.write('Runs a benchmark according to a given config file.\n')
            sys.stderr.write('Config file contains a list of assemblers to be used, reads file and reference file.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <config_file>' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)


        configfile = sys.argv[2]
        benchmark_cf(configfile)

    elif mode == '--continue_cf':
        if (len(sys.argv) != 3):
            sys.stderr.write('Continues a benchmark according to a given config file.\n')
            sys.stderr.write('Config file contains a list of assemblers to be used, reads file and reference file.\n')
            sys.stderr.write('Reads a corresponding log file to see what is left to be done.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <results_folder>' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        results_folder = sys.argv[2]
        continueBenchmark_cf(reads_folder)

    elif mode == '--run_quast':
        if (len(sys.argv) != 3):
            sys.stderr.write('Runs quast on a set of assembler results.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <results_folder>' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        results_folder = sys.argv[2]
        run_quast(results_folder)

    elif mode == '--summarize_results':
        if (len(sys.argv) != 3):
            sys.stderr.write('Summarizes results from quast and cgmem timefiles.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <results_folder>' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        results_folder = sys.argv[2]
        summarize_results(results_folder)

    else:
        print 'Invalid option!'



# Standard template for the main() function
if __name__ == '__main__':
  main()
