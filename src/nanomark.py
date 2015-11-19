#! /usr/bin/python

import os
import sys
import re
import shutil
import datetime
import subprocess

import setup_nanomark
import basicdefines
import fastqparser
import parseresults

# To generate random number for each run, so that multiple runs don't clash
import uuid
from time import gmtime, strftime

# To be able to import wrappers more easily
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(SCRIPT_PATH + '')
sys.path.append(SCRIPT_PATH + '/../wrappers');

from dataspec import *



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

    sys.stdout.write('\n')



def benchmark(reads_file, reference_file, technology, wrapper_list = []):
    reads_file = os.path.abspath(reads_file);
    reference_file = os.path.abspath(reference_file);

    # uuid_string = str(uuid.uuid4())            # Generate a random UUID so that multiple runs don't clash
    timestamp = strftime("%Y_%m_%d-%H_%M_%S", gmtime());
    uuid_string = timestamp;

    output_folder = 'benchmark_' + uuid_string
    output_path = os.path.join(basicdefines.RESULTS_PATH_ROOT_ABS, output_folder)
    log_filename = os.path.join(output_path, 'log.txt')
    # creffilename = 'REF_' + uuid_string + os.path.splitext(reference_file)[1]        # preserve extension
    # crefpath = os.path.join(output_path, creffilename)
    # creadsfilename = 'READS_' + uuid_string + os.path.splitext(reads_file)[1]        # preserve extension
    # creadspath = os.path.abspath(os.path.join(output_path, creadsfilename))

    [ret_string, num_refs, total_ref_len, average_ref_len] = fastqparser.count_seq_length(reference_file);

    if (not os.path.exists(output_path)):
        sys.stderr.write('Creating folder "%s".\n' % (output_path))
        os.makedirs(output_path)

    with open(log_filename, 'w') as logfile:
        cDateTime = str(datetime.datetime.now())
        # logfile.write('NanoMark benchmark UUID: %s\n' % uuid_string)
        logfile.write('Date and time: %s\n' % cDateTime)
        sys.stderr.write('\n')
        sys.stderr.write('Running benchmark for one reads file!\n')
        sys.stderr.write('UUID: %s\n' % uuid_string)
        sys.stderr.write('Reads file: %s\n' % reads_file)
        sys.stderr.write('Reference file: %s\n' % reference_file)

        # Copy reference file to output folder while preserving extension
        # shutil.copy(reference_file, crefpath)
        logfile.write('Reference file: %s\n' % reference_file)

        # Copy reads file to output folder while preserving extension
        # shutil.copy(reads_file, creadspath)
        logfile.write('Reads file: %s\n' % reads_file)

        logfile.write('Technology: %s\n' % technology)

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
        logfile.write('Wrappers used: %s\n' % ','.join(wrapper_list))

        for wrapper_basename in wrapper_list:
            assembler_name = ''
            results_file = ''
            create_output_folder = True
            wrapper = 'wrapper_' + wrapper_basename
            wrapper_path = os.path.join(basicdefines.WRAPPERS_PATH_ROOT_ABS, wrapper + '.py')

            sys.stderr.write('Importing wrapper: "%s".\n' % (wrapper));
            exec('import %s as current_wrapper' % (wrapper))
            assembler_name = current_wrapper.ASSEMBLER_NAME;
            assembly_unpolished = current_wrapper.ASSEMBLY_UNPOLISHED;
            assembly_polished = current_wrapper.ASSEMBLY_POLISHED;
            create_output_folder = current_wrapper.CREATE_OUTPUT_FOLDER;
            assembler_type = current_wrapper.ASSEMBLER_TYPE;

            # Run assembly only using non-hybrid assemblers. Using hybrid assemblers complicates uniform specification of the datasets.
            if (assembler_type != 'nonhybrid'):
                continue;
            # logfile.write('%s (%s): %s\n' % (assembler_name, wrapper, wrapper_path))
            logfile.write(basicdefines.log_message('%s (%s) started: %s' % (assembler_name, wrapper, wrapper_path)));
        
            # Create folder for each assembler's results, and for quast results
            assembler_folder = os.path.join(output_path, assembler_name)
            sys.stderr.write('\nRunning assembler %s\n' % assembler_name)
            dataset = Dataset('%s,%s' % (technology, reads_file))
            current_wrapper.run([dataset], assembler_folder, total_ref_len, move_exiting_out_path=True);

            logfile.write(basicdefines.log_message('%s (%s) finished.' % (assembler_name, wrapper)));


def parse_log_file(log_file_path):
    try:
        fp_log = open(log_file_path, 'r');
        log_lines = fp_log.readlines();
        fp_log.close();
    except IOError, e:
        sys.stderr.write('ERROR: Could not open log file "%s"! Exiting.\n' % (log_file_path));
        exit(1);

    # Load info about the input dataset.
    reference_file = '';
    reads_file = '';
    technology = '';
    # uuid_string = '';
    wrappers_used = '';
    for log_line in log_lines:
        split_line = log_line.strip().split(':');
        if (split_line[0] == 'Reference file'):
            reference_file = split_line[-1].strip();
        elif (split_line[0] == 'Reads file'):
            reads_file = split_line[-1].strip();
        elif (split_line[0] == 'Technology'):
            technology = split_line[-1].strip();
        # elif (split_line[0] == 'UUID'):
        #     uuid_string = split_line[-1].strip();
        elif (split_line[0] == 'Wrappers used'):
            wrappers_used = split_line[-1].strip();
    if (reference_file == '' or reads_file == '' or technology == '' or wrappers_used == ''):
        sys.stderr.write('ERROR: Could not load all parameters from log file to continue the assembly.\n');
        sys.stderr.write('Reference file: "%s"\n' % (reference_file));
        sys.stderr.write('Reads file: "%s"\n' % (reads_file));
        sys.stderr.write('Technology: "%s"\n' % (technology));
        # sys.stderr.write('UUID: "%s"\n' % (uuid_string));
        exit(1);

    return [reference_file, reads_file, technology, wrappers_used];

# Continue benchmark that was interrupted
# Takes benchmark results folder as an argument
def continue_benchmark(results_folder):
    sys.stderr.write('\nContinuing benchmark in folder: %s\n' % results_folder)
    # First check if log file exists
    log_filename = os.path.join(results_folder, 'log.txt')
    if not os.path.exists(log_filename):
        sys.stderr.write('\nCannot find log file (log.txt) within benchmark results folder (%s)! Exiting ...\n' % results_folder)
        verbose_usage_and_exit()

    sys.stderr.write('Reading log file and checking wrapers....\n')
    [reference_file, reads_file, technology, wrappers_used] = parse_log_file(log_file_path);

    reference_file = os.path.abspath(reference_file);
    reads_file = os.path.abspath(reads_file);
    output_path = results_folder;

    [ret_string, num_refs, total_ref_len, average_ref_len] = fastqparser.count_seq_length(reference_file);

    wrapper_list = wrappers_used.split(',');

    logfile = open(log_filename, 'a');
    for wrapper_basename in wrapper_list:
        wrapper = 'wrapper_' + wrapper_basename
        wrapper_path = os.path.join(basicdefines.WRAPPERS_PATH_ROOT_ABS, wrapper + '.py')

        exec('import %s as current_wrapper' % (wrapper))
        assembler_name = current_wrapper.ASSEMBLER_NAME;
        assembler_folder = os.path.join(output_path, assembler_name)
        assembly_unpolished = os.path.join(assembler_folder, current_wrapper.ASSEMBLY_UNPOLISHED);
        assembly_polished = os.path.join(assembler_folder, current_wrapper.ASSEMBLY_POLISHED);
        create_output_folder = current_wrapper.CREATE_OUTPUT_FOLDER;
        assembler_type = current_wrapper.ASSEMBLER_TYPE;
        # Run assembly only using non-hybrid assemblers. Using hybrid assemblers complicates uniform specification of the datasets.
        if (assembler_type != 'nonhybrid'):
            continue;
        if (os.path.exists(assembly_unpolished) or os.path.exists(assembly_polished)):
            # The results file exists, thi means that the assembler run completed and will not be repeated
            if (os.path.exists(assembly_unpolished)):
                sys.stderr.write('Previously constructed (unpolished) assembly found at: "%s".\n' % (assembly_unpolished));
            if (os.path.exists(assembly_polished)):
                sys.stderr.write('Previously constructed (polished) assembly found at: "%s".\n' % (assembly_polished));
            sys.stderr.write('Assembler %s run previously completed.\n' % (assembler_name))
        else:

            logfile.write(basicdefines.log_message('%s (%s) started: %s' % (assembler_name, wrapper, wrapper_path)));
            sys.stderr.write('\nRunning assembler %s\n' % assembler_name)
            dataset = Dataset('%s,%s' % (technology, reads_file))
            # The wrapper internally takes care of backing up existing output folders.
            current_wrapper.run([dataset], assembler_folder, total_ref_len, move_exiting_out_path=True);
            logfile.write(basicdefines.log_message('%s (%s) finished.' % (assembler_name, wrapper)));
    logfile.close();

def summarize_results(results_folder):
    run_type = 'calc';  # 'calc' or 'collectall'.
    quast_parameters = [
                            '# contigs (>= 0 bp)', 'Largest contig', 'Total length (>= 0 bp)', 'N50', 'Unaligned length',
                            '# mismatches per 100 kbp', '# indels per 100 kbp', 'Genome fraction (%)', 'Duplication ratio'
                        ];
    dnadiff_parameters = [
                            '[Bases]-TotalBases,query', '[Bases]-AlignedBases,ref', '[Bases]-AlignedBases,query', '[Bases]-UnalignedBases,ref',
                            '[Bases]-UnalignedBases,query', '[Alignments]-AvgIdentity-1-to-1,query', '[Feature Estimates]-Breakpoints,ref',
                            '[Feature Estimates]-Breakpoints,query', '[SNPs]-TotalSNPs,query', '[SNPs]-TotalIndels,query'
                         ];
    memtime_parameters = ['Real time', 'CPU time', 'Maximum RSS'];

    # get global quast folder and log file
    gl_quast_folder = os.path.join(results_folder, 'quast')
    gl_dnadiff_folder = os.path.join(results_folder, 'dnadiff')
    log_filename = os.path.join(results_folder, 'log.txt')
    # # Collect quast results
    summaryfilepath = os.path.join(results_folder, 'benchmark_summary.tsv')
    # report_filename = 'report.txt'

    sys.stderr.write('Reading log file and checking wrapers....\n')
    [reference_file, reads_file, technology, wrappers_used] = parse_log_file(log_filename);

    wrapper_list = wrappers_used.split(',');
    reference_file = os.path.abspath(reference_file);
    reads_file = os.path.abspath(reads_file);
    output_path = results_folder;

    # Collect paths to all assembled FASTA files.
    results_files = [];
    for wrapper_basename in wrapper_list:
        wrapper = 'wrapper_' + wrapper_basename
        wrapper_path = os.path.join(basicdefines.WRAPPERS_PATH_ROOT_ABS, wrapper + '.py')

        exec('import %s as current_wrapper' % (wrapper))
        assembler_name = current_wrapper.ASSEMBLER_NAME;
        assembler_folder = os.path.join(output_path, assembler_name)
        assembly_unpolished = os.path.join(assembler_folder, current_wrapper.ASSEMBLY_UNPOLISHED);
        assembly_polished = os.path.join(assembler_folder, current_wrapper.ASSEMBLY_POLISHED);
        create_output_folder = current_wrapper.CREATE_OUTPUT_FOLDER;
        assembler_type = current_wrapper.ASSEMBLER_TYPE;

        results_files.append([assembly_unpolished, assembler_name]);
        results_files.append([assembly_polished, assembler_name]);

    try:
        fp_out = open(summaryfilepath, 'w');
    except:
        sys.stderr.write('ERROR: Could not open file "%s" for writing! Exiting.\n' % (summaryfilepath));
        exit(1);

    header = ['Assembler', 'Assembly'] + quast_parameters + dnadiff_parameters + memtime_parameters;
    fp_out.write('%s\n' % (';'.join(header)));
    fp_out.flush();

    # Run Quast and DNAdiff on the assemblies.
    for results_file in results_files:
        for results_file in results_files:
            [assembly_path, assembler_name] = results_file;
            contig_basename = os.path.splitext(os.path.basename(assembly_path))[0];

            sys.stderr.write('Processing contig file "%s".\n' % (assembly_path));
            quast_out_folder = '%s/%s/%s' % (gl_quast_folder, assembler_name, contig_basename);
            dnadiff_out_folder = '%s/%s/%s' % (gl_dnadiff_folder, assembler_name, contig_basename);

            if (os.path.exists(assembly_path)):
                if (run_type == 'calc'):
                    sys.stderr.write('Running analyses on the assembly output.\n');
                    basicdefines.execute_command('mkdir -p %s' % (quast_out_folder), None);
                    basicdefines.execute_command('mkdir -p %s' % (dnadiff_out_folder), None);
                    basicdefines.execute_command('%s %s -R %s -o %s' % (basicdefines.QUAST_BIN, assembly_path, reference_file, quast_out_folder), None);
                    basicdefines.execute_command('%s %s %s -p %s/out' % (basicdefines.DNADIFF_BIN, reference_file, assembly_path, dnadiff_out_folder), None);
                pass;
            else:   # In this case, the contig file does not exist.
                if (run_type != 'collectall'):  # This will skip processing of non-exitent file. Otherwise, '-' will be placed for that particular assembly.
                    continue;

            results_dnadiff = parseresults.parse_quast_report('%s/report.tsv' % (quast_out_folder), quast_parameters);
            results_quast = parseresults.parse_dnadiff_report('%s/out.report' % (dnadiff_out_folder), dnadiff_parameters);
            parseresults.parse_memtime_folder_and_accumulate('%s' % (os.path.dirname(assembly_path)), '%s/total.memtime' % (os.path.dirname(assembly_path)));
            results_memtime = parseresults.parse_memtime_report('%s/total.memtime' % (os.path.dirname(assembly_path)), memtime_parameters, 'h', 'GB');

            fp_out.write('%s\n' % (';'.join([assembler_name, contig_basename] + results_dnadiff + results_quast + results_memtime)));
            fp_out.flush();

    fp_out.close();

def summarize_results2(results_folder):

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
    sys.stderr.write('\t\tsummarize\n')
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
                        sys.stderr.write('\nUnavailable assembler: %s\n' % assembler)
                        exit(1)
                    else:
                        used_assemblers.append(assembler)
            else:
                sys.stderr.write('Invalid parameter!.\n')
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
        continue_benchmark(results_folder)

    elif mode == 'summarize':
        if (len(sys.argv) != 3):
            sys.stderr.write('Summarizes results from Quast, DNAdiff and cgmem timefiles.\n')
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
