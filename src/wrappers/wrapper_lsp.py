#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(SCRIPT_PATH + '/../')

import re;

import subprocess
import multiprocessing
import getpass
from time import gmtime, strftime

try:
    import basicdefines
    MODULE_BASICDEFINES = True;
    ASSEMBLERS_PATH_ROOT_ABS = basicdefines.ASSEMBLERS_PATH_ROOT_ABS;
    TOOLS_ROOT = basicdefines.TOOLS_ROOT;
except:
    MODULE_BASICDEFINES = False;
    ASSEMBLERS_PATH_ROOT_ABS = os.path.join(SCRIPT_PATH, 'assemblers/');
    TOOLS_ROOT = '%s' % (SCRIPT_PATH);

# Loman Simpson pipeline
# https://github.com/jts/nanopore-paper-analysis/blob/master/full-pipeline.make

ASSEMBLER_URL = 'https://github.com/jts/nanopore-paper-analysis.git'
ASSEMBLER_PATH = os.path.join(basicdefines.ASSEMBLERS_PATH_ROOT_ABS, 'LSP')
# ZIP_FILE = 'wgs-8.3rc2-Linux_amd64.tar.bz2'
# ZIP_PATH = os.path.join(basicdefines.ASSEMBLERS_PATH_ROOT_ABS, ZIP_FILE)
ASSEMBLER_BIN = os.path.join(ASSEMBLER_PATH, 'full-pipeline.make')
ASSEMBLER_NAME = 'LSP'
ASSEMBLER_RESULTS = 'contig-100.fa'
CREATE_OUTPUT_FOLDER = True

# FASTQ to FASTA converter supplied by IDBA
FQ2FA_BIN = os.path.join(ASSEMBLER_PATH, 'bin/fq2fa')

# DRY_RUN = True;
DRY_RUN = False;

### Logs messages to STDERR and an output log file if provided (opened elsewhere).
def log(message, fp_log):
    timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());

    sys.stderr.write('[%s wrapper %s] %s\n' % (ASSEMBLER_NAME, timestamp, message))
    sys.stderr.flush();
    if (fp_log != None):
        fp_log.write('[%s wrapper %s] %s\n' % (ASSEMBLER_NAME, timestamp, message))
        fp_log.flush();

def execute_command(command, fp_log, dry_run=True):
    if (dry_run == True):
        # sys.stderr.write('Warning: dry_run == True\n');
#        sys.stderr.write('[%s wrapper] Executing (dryrun): "%s"\n' % (ASSEMBLER_NAME, command));
        log('Executing (dryrun): "%s".' % (command), fp_log);
    if (dry_run == False):
        # sys.stderr.write('[%s wrapper] Executing: "%s"\n' % (ASSEMBLER_NAME, command));
        log('Executing: "%s".' % (command), fp_log);
        subprocess.call(command, shell=True);
    # sys.stderr.write('\n');
    log('\n', fp_log);

def measure_command(measure_file):
    if (MODULE_BASICDEFINES == True):
        return basicdefines.measure_command(measure_file);
    else:
        sys.stderr.write('ERROR: Cgmemtime tool not found! Exiting.\n');
        exit(1);
        # return '/usr/bin/time --format "Command line: %%C\\nReal time: %%e s\\nCPU time: -1.0 s\\nUser time: %%U s\\nSystem time: %%S s\\nMaximum RSS: %%M kB\\nExit status: %%x" --quiet -o %s ' % measure_file;

def peek(fp, num_chars):
    data = fp.read(num_chars);
    if len(data) == 0:
        return '';
    fp.seek(num_chars * -1, 1);
    return data;

# Returns a single read from the given FASTA/FASTQ file.
# Parameter header contains only the header of the read.
# Parameter lines contains all lines of the read, which include:
# - header
# - seq
# - '+' if FASTQ
# - quals if FASTQ
# Parameter lines is an array of strings, each for one component.
# Please note that multiline FASTA/FASTQ entries (e.g. sequence line)
# will be truncated into one single line.
def get_single_read(fp):
    lines = [];
    
    line = fp.readline();
    header = line.rstrip();
    header_leading_char = '';
    if (len(header) > 0):
        sequence_separator = header[0];
        header_leading_char = header[0];
        header = header[1:];            # Strip the '>' or '@' sign from the beginning.
    else:
        return ['', []];
    
    next_char = peek(fp, 1);
    
    line_string = '';
    lines.append(header_leading_char + header);
    
    num_lines = 1;
    #while len(next_char) > 0 and next_char != sequence_separator or (next_char == '@' and num_lines < 4):
    while (len(next_char) > 0 and (next_char != sequence_separator or (next_char == '@' and num_lines < 4))):
        line = fp.readline();
        if (line.rstrip() == '+' or line.rstrip() == ('+' + header)):
        #if (line.rstrip()[0] == '+'):
            lines.append(line_string);
            lines.append(line.rstrip());
            line_string = '';
        else:
            line_string += line.rstrip();
        next_char = peek(fp, 1);
        num_lines += 1;
        
    lines.append(line_string);
    
    return [header, lines];

def convert_to_fasta(fastq_path, out_fasta_path):
    fp_in = None;
    fp_out = None;
    
    try:
        fp_in = open(fastq_path, 'r');
    except IOError:
        print 'ERROR: Could not open file "%s" for reading!' % fastq_path;
        return;
    
    try:
        fp_out = open(out_fasta_path, 'w');
    except IOError:
        print 'ERROR: Could not open file "%s" for writing!' % out_fasta_path;
        fp_in.close();
        return;
    
    while True:
        [header, read] = get_single_read(fp_in);
        if (len(header) == 0):
            break;
        seq = read[1];
        fp_out.write('>' + header + '\n');
        fp_out.write(seq + '\n');
    fp_in.close();
    fp_out.close();

##############################################################
##############################################################
##############################################################



# Function 'run' should provide a standard interface for running a mapper. Given input parameters, it should run the
# alignment process, and convert any custom output results to the SAM format. Function should return a string with the
# path to the output file.
#    reads_file            Path to a FASTA/FASTQ file containing reads.
#    reference_file        Path to a reference genome FASTA file.
#    machine_name        A symbolic name to specify a set of parameters for a specific sequencing platform.
#    output_path            Folder to which the output will be placed to. Filename will be automatically generated according to the name of the mapper being run.
#    output_suffix        A custom suffix that can be added to the output filename.
def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):
    ### make -f full-pipeline.make CORES=64 polished_genome.fasta
    # num_threads = multiprocessing.cpu_count() / 2;
    # num_cores = 2;
    num_cores = 16;
    num_threads = 4;

    reads_file = os.path.abspath(reads_file);
    reference_file = os.path.abspath(reference_file);
    output_path = os.path.abspath(output_path);
    reads_folder = os.path.dirname(reads_file);
    reads_basename = os.path.basename(reads_file);
    raw_reads_path = '%s/raw.reads.fasta' % (output_path);
    raw_reads_filename = os.path.basename(raw_reads_path);
    raw_reads_basename = os.path.splitext(os.path.basename(raw_reads_path))[0];

    ### Backup old assembly results, and create the new output folder.
    if (os.path.exists(output_path)):
        timestamp = strftime("%Y_%m_%d-%H_%M_%S", gmtime());
        os.rename(output_path, '%s.bak_%s' % (output_path, timestamp));
    if (not os.path.exists(output_path)):
        log('Creating a directory on path "%s".' % (output_path), None);
        os.makedirs(output_path);

    log_file = '%s/wrapper_log.txt' % (output_path);
    try:
        fp_log = open(log_file, 'a');
    except Exception, e:
        log('ERROR: Could not open file "%s" for writing! Using only STDERR for logging.' % (log_file), None);
        # sys.stderr.write(str(e) + '\n');
        fp_log = None;


    log('Preparing raw reads.', fp_log);

    if (os.path.exists(raw_reads_path) == True):
        os.rename(raw_reads_path, raw_reads_path + '.bak');
    ### Generate a raw reads file in the output folder, which will be used for assembly.
    convert_to_fasta(reads_file, raw_reads_path);



    log('Running assembly using %s.' % (ASSEMBLER_NAME), fp_log);

    if (MODULE_BASICDEFINES == True):
        command = 'sudo %s/cgmemtime/cgmemtime --setup -g %s --perm 775' % (basicdefines.TOOLS_ROOT, getpass.getuser());
        execute_command(command, fp_log, dry_run=DRY_RUN);
    memtime_file = '%s/%s.memtime' % (output_path, ASSEMBLER_NAME);

    commands = [];
    commands.append('cd %s' % (output_path));
    # The programs we will install must be on the PATH
    commands.append('export PATH=%s/DAZZ_DB:%s/DALIGNER:%s/nanocorrect:%s/poaV2:%s/wgs-8.2/Linux-amd64/bin/:%s/samtools/:%s/bwa/:%s/:$PATH' % (ASSEMBLER_PATH, ASSEMBLER_PATH, ASSEMBLER_PATH, ASSEMBLER_PATH, ASSEMBLER_PATH, ASSEMBLER_PATH, ASSEMBLER_PATH, reads_folder));
    # Parameters to control execution
    commands.append('CORES=%d' % (num_cores));
    commands.append('THREADS=%d' % (num_threads));
    commands.append('NC_PROCESS=%d' % (num_cores));
    commands.append('NP_PROCESS=%d' % (num_cores / num_threads));

    commands.append('LS_ENV=%s/ls_env' % (ASSEMBLER_PATH));
    commands.append('. $LS_ENV/bin/activate');

    # Error correction, the first step.
    commands.append('make -f nanocorrect/nanocorrect-overlap.make INPUT=%s NAME=%s' % (raw_reads_path, raw_reads_basename));
    commands.append('samtools faidx %s.pp.fasta' % (raw_reads_basename));
    commands.append('python nanocorrect/makerange.py %s | parallel -v --eta -P $NC_PROCESS \'python nanocorrect/nanocorrect.py %s {} > %s.{}.corrected.fasta\'' % (raw_reads_path, raw_reads_basename, raw_reads_basename));
    commands.append('cat %s.*.corrected.fasta | python lengthsort.py > %s.corrected.fasta' % (raw_reads_basename, raw_reads_basename));
    #rm raw.reads.*.corrected.fasta
    # Error correction, the second step.
    commands.append('make -f nanocorrect/nanocorrect-overlap.make INPUT=%s.corrected.fasta NAME=%s.corrected' % (raw_reads_basename, raw_reads_basename));
    commands.append('samtools faidx %s.corrected.pp.fasta' % (raw_reads_basename));
    commands.append('python nanocorrect/makerange.py %s.corrected.fasta | parallel -v --eta -P $NC_PROCESS \'python nanocorrect/nanocorrect.py %s.corrected {} > %s.corrected.{}.corrected.fasta\'' % (raw_reads_basename, raw_reads_basename, raw_reads_basename));
    commands.append('cat %s.corrected.*.corrected.fasta | python lengthsort.py > %s.corrected.corrected.fasta' % (raw_reads_basename));
    # commands.append('rm raw.reads.corrected.*.corrected.fasta');

    reads_error_corrected = 'raw.reads.corrected.corrected.fasta';
    reads_assembly_input = 'assembly.input.fastq';
    commands.append('java -Xmx1024M -jar ./wgs-8.2/Linux-amd64/bin/convertFastaAndQualToFastq.jar %s > %s' % (reads_error_corrected, reads_assembly_input));
    commands.append('fastqToCA -technology sanger -libraryname assembly -reads %s > assembly.frg' % (reads_assembly_input));
    commands.append('runCA -d celera-assembly -p asm -s revised_ovlErrorRate0.04.spec assembly.frg');
    commands.append('ln -s celera-assembly/9-terminator/asm.scf.fasta draft_genome.fasta');

    # preprocess the fasta file for nanopolish
    commands.append('nanopolish/consensus-preprocess.pl %s > %s.np.fasta' % (raw_reads_path, raw_reads_basename));
    # index the draft assembly for bwa
    commands.append('bwa index draft_genome.fasta');
    # index the draft assembly for faidx
    commands.append('samtools faidx draft_genome.fasta');
    # align reads to draft assembly
    commands.append('bwa mem -t $THREADS -x ont2d draft_genome.fasta %s.np.fasta | samtools view -Sb - | samtools sort -f - reads_to_draft.sorted.bam' % (raw_reads_basename));
    # index the bam file
    commands.append('samtools index reads_to_draft.sorted.bam');
    # run nanopolish
    commands.append('python nanopolish/nanopolish_makerange.py draft_genome.fasta | parallel --progress -P $NP_PROCESS nanopolish/nanopolish consensus -o nanopolish.{1}.fa -r %s.np.fasta -b reads_to_draft.sorted.bam -g draft_genome.fasta -w {1} -t $THREADS python nanopolish/nanopolish_merge.py draft_genome.fasta nanopolish.scf*.fa > polished_genome.fasta' % (raw_reads_basename));

    commands.append('cp %s/polished_genome.fasta %s/benchmark-final_assembly.fasta' % (output_path, output_path));

    command = '; '.join(command);
    execute_command(command, None, dry_run=DRY_RUN);



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
def download_and_install():
    if os.path.exists(ASSEMBLER_BIN):
        # sys.stderr.write('[%s wrapper] Bin found at %s. Skipping installation.\n' % (ASSEMBLER_NAME, ASSEMBLER_BIN))
        log('Bin found at %s. Skipping installation.' % (ASSEMBLER_BIN), None);

    else:
        # sys.stderr.write('[%s wrapper] Started installation of %s.\n' % (ASSEMBLER_NAME, ASSEMBLER_NAME))
        log('Started installation of %s.' % (ASSEMBLER_NAME), None);
        if (not os.path.exists(ASSEMBLER_PATH)):
            # sys.stderr.write('[%s wrapper] Creating a directory on path "%s".\n' % (ASSEMBLER_NAME, ASSEMBLERS_PATH_ROOT_ABS))
            log('Creating a directory on path "%s".' % (ASSEMBLER_PATH), None);
            os.makedirs(ASSEMBLER_PATH);

        # command = 'cd %s; git clone %s' % (ASSEMBLERS_PATH_ROOT_ABS, ASSEMBLER_URL)
        # execute_command(command, None, dry_run=DRY_RUN);

        setup_commands = [];
        # The programs we will install must be on the PATH');
        # Install the BioPerl dependency.
        setup_commands.append('cd %s' % (ASSEMBLER_PATH));
        setup_commands.append('sudo apt-get install bioperl');
        setup_commands.append('sudo apt-get install parallel');
        setup_commands.append('sudo pip install virtualenv');
        # Install Python dependencies.
        setup_commands.append('LS_ENV=%s/ls_env' % (ASSEMBLER_PATH));
        setup_commands.append('virtualenv --no-site-packages --always-copy $LS_ENV');
        # Activate the virtual environment.
        setup_commands.append('. $LS_ENV/bin/activate');
        # Install the required Python packages in the virtual environment.
        setup_commands.append('pip install pysam > pythonlibs.version');
        setup_commands.append('pip install cython >> pythonlibs.version');
        setup_commands.append('pip install numpy==1.8.1 >> pythonlibs.version');
        setup_commands.append('pip install h5py==2.3.0 >> pythonlibs.version');
        setup_commands.append('pip install cython >> pythonlibs.version');
        setup_commands.append('pip install poretools >> pythonlibs.version');
        setup_commands.append('pip install biopython >> pythonlibs.version');
        setup_commands.append('pip freeze >> pythonlibs.version');
        # Install samtools
        setup_commands.append('git clone --recursive https://github.com/samtools/htslib.git');
        setup_commands.append('cd htslib; make; cd ..');
        setup_commands.append('git clone --recursive https://github.com/samtools/samtools.git');
        setup_commands.append('cd samtools; make; cd ..');
        setup_commands.append('cd samtools; git log | head -1 > ../samtools.version');
        setup_commands.append('cd ..');
        # Install nanocorrect & dependencies
        setup_commands.append('git clone https://github.com/jts/nanocorrect.git');
        setup_commands.append('ln -s nanocorrect/poa-blosum80.mat');
        setup_commands.append('cd nanocorrect; git checkout 47dcd7f147c; git log | head -1 > ../nanocorrect.version');
        setup_commands.append('cd ..');
        # Install nanopolish, automatically downloading libhdf5
        setup_commands.append('git clone --recursive https://github.com/jts/nanopolish.git');
        # setup_commands.append('cd nanopolish; git checkout 6440bfbfcf4fa; make libhdf5.install nanopolish; cd ..');
        setup_commands.append('cd nanopolish; make libhdf5.install nanopolish; cd ..');
        setup_commands.append('cd nanopolish; git log | head -1 > ../nanopolish.version');
        setup_commands.append('cd ..');
        # Install bwa
        setup_commands.append('git clone https://github.com/lh3/bwa.git');
        setup_commands.append('cd bwa; make');
        setup_commands.append('cd bwa; git log | head -1 > ../bwa.version');
        setup_commands.append('cd ..');
        # Install poa
        setup_commands.append('wget http://downloads.sourceforge.net/project/poamsa/poamsa/2.0/poaV2.tar.gz');
        setup_commands.append('tar -xzf poaV2.tar.gz');
        setup_commands.append('cd poaV2; make CFLAGS=\'-O3 -g -DUSE_WEIGHTED_LINKS -DUSE_PROJECT_HEADER -I.\' poa');
        setup_commands.append('cd ..');
        setup_commands.append('ln -s poaV2/poa');
        setup_commands.append('echo http://downloads.sourceforge.net/project/poamsa/poamsa/2.0/poaV2.tar.gz > poa.version');
        # Install DALIGNER
        setup_commands.append('git clone https://github.com/thegenemyers/DALIGNER.git');
        setup_commands.append('cd DALIGNER; git checkout 549da77b91395dd; make; cd ..');
        setup_commands.append('echo "549da77b91395dd" > daligner.version');
        # Install DAZZ_DB
        setup_commands.append('git clone https://github.com/thegenemyers/DAZZ_DB');
        setup_commands.append('cd DAZZ_DB; git checkout 8cb2f29c4011a2c2; make; cd ..');
        setup_commands.append('echo "8cb2f29c4011a2c2" > dazz_db.version');
        # Install Celera Assembler
        setup_commands.append('wget http://downloads.sourceforge.net/project/wgs-assembler/wgs-assembler/wgs-8.2/wgs-8.2-Linux_amd64.tar.bz2');
        setup_commands.append('tar -xjf wgs-8.2-Linux_amd64.tar.bz2');
        setup_commands.append('wget http://www.cbcb.umd.edu/software/PBcR/data/convertFastaAndQualToFastq.jar');
        setup_commands.append('mv convertFastaAndQualToFastq.jar wgs-8.2/Linux-amd64/bin/');
        setup_commands.append('echo http://downloads.sourceforge.net/project/wgs-assembler/wgs-assembler/wgs-8.2/wgs-8.2-Linux_amd64.tar.bz2 > ca.version');
        # Install lengthsort.
        setup_commands.append('wget https://raw.githubusercontent.com/jts/nanopore-paper-analysis/master/lengthsort.py');
        setup_commands.append('echo "lengthsort" > lengthsort.version');
        # Download the spec file for Celera.
        setup_commands.append('wget --no-check-certificate https://raw.githubusercontent.com/jts/nanopore-paper-analysis/c25373d93a99e51c2fedb57d8b08b81826e7c80c/revised_ovlErrorRate0.04.spec');

        command = '; '.join(setup_commands);
        execute_command(command, None, dry_run=DRY_RUN);



def verbose_usage_and_exit():
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s mode [<reads_file1>,<reads_file2>,...,<reads_fileN> <machine_name> <output_path> <reference_file> [<output_suffix>]]\n' % sys.argv[0])
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
        if (len(sys.argv) < 5 or len(sys.argv) > 7):
            verbose_usage_and_exit()

        reads_files = sys.argv[2].split(',')         ### Enable specifying multiple FASTQ files for input.
        machine_name = sys.argv[3]
        output_path = sys.argv[4]
        reference_file = '';
        output_suffix = ''

        if (len(sys.argv) >= 6):
            reference_file = sys.argv[5]
        if (len(sys.argv) == 7):
            output_suffix = sys.argv[6]
        run(reads_files, reference_file, machine_name, output_path, output_suffix)

    else:
        verbose_usage_and_exit()
