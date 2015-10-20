#! /bin/sh

##################################################
#
# Step 0. Preamble: set up paths and variable
#
##################################################

# The programs we will install must be on the PATH
export PATH := ./DAZZ_DB:./DALIGNER:./nanocorrect:./poaV2:./wgs-8.2/Linux-amd64/bin/:./samtools/:./bwa/:$(PATH)

# Download links for programs that are not on github
CA_LINK=http://downloads.sourceforge.net/project/wgs-assembler/wgs-assembler/wgs-8.2/wgs-8.2-Linux_amd64.tar.bz2
POA_LINK=http://downloads.sourceforge.net/project/poamsa/poamsa/2.0/poaV2.tar.gz

# Parameters to control execution
CORES=16
THREADS=4
NC_PROCESS=$(CORES)
NP_PROCESS=$(shell expr $(CORES) / $(THREADS) )



##################################################
#
# Step 2. Download and install the programs
# required to run the assembly
#
##################################################

# Install Python dependencies
LS_ENV=$PWD/ls_env
virtualenv --no-site-packages --always-copy $LS_ENV
. $LS_ENV/bin/activate

pip install pysam > pythonlibs.version
pip install cython >> pythonlibs.version
pip install numpy==1.8.1 >> pythonlibs.version
pip install h5py==2.3.0 >> pythonlibs.version
pip install cython >> pythonlibs.version
pip install poretools >> pythonlibs.version
pip install biopython >> pythonlibs.version
pip freeze >> pythonlibs.version

# Install samtools
git clone --recursive https://github.com/samtools/htslib.git
cd htslib; make
git clone --recursive https://github.com/samtools/samtools.git
cd samtools; make
-cd samtools; git log | head -1 > ../samtools.version

# Install nanocorrect & dependencies
git clone https://github.com/jts/nanocorrect.git
ln -s nanocorrect/poa-blosum80.mat
-cd nanocorrect; git checkout 47dcd7f147c; git log | head -1 > ../nanocorrect.version

# Install nanopolish, automatically downloading libhdf5
git clone --recursive https://github.com/jts/nanopolish.git
cd nanopolish; git checkout 6440bfbfcf4fa; make libhdf5.install nanopolish
-cd nanopolish; git log | head -1 > ../nanopolish.version

# Install bwa
git clone https://github.com/lh3/bwa.git
cd bwa; make
-cd bwa; git log | head -1 > ../bwa.version

# Install poa
wget $(POA_LINK)
tar -xzf poaV2.tar.gz
cd poaV2; make CFLAGS='-O3 -g -DUSE_WEIGHTED_LINKS -DUSE_PROJECT_HEADER -I.' poa
ln -s poaV2/poa
echo $(POA_LINK) > poa.version

# Install DALIGNER
git clone https://github.com/thegenemyers/DALIGNER.git
cd DALIGNER; git checkout 549da77b91395dd; make
echo "549da77b91395dd" > daligner.version

# Install DAZZ_DB
git clone https://github.com/thegenemyers/DAZZ_DB
cd DAZZ_DB; git checkout 8cb2f29c4011a2c2; make
echo "8cb2f29c4011a2c2" > dazz_db.version

# Install Celera Assembler
wget $(CA_LINK)
tar -xjf wgs-8.2-Linux_amd64.tar.bz2
wget http://www.cbcb.umd.edu/software/PBcR/data/convertFastaAndQualToFastq.jar
mv convertFastaAndQualToFastq.jar wgs-8.2/Linux-amd64/bin/
echo $(CA_LINK) > ca.version

# Install lengthsort.
wget https://raw.githubusercontent.com/jts/nanopore-paper-analysis/master/lengthsort.py
echo "lengthsort" > lengthsort.version



##################################################
#
# Step 1. Download the input data from the ENA
# and unpack it into fast5 files
#
##################################################

# Download the data and prepare reads.
wget ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_20_LomanLabz_PC_Ecoli_K12_R7.3.tar -O ERX708228.tar
wget ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_32_LomanLabz_K12_His_tag.tar -O ERX708229.tar
wget ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_33_LomanLabz_PC_K12_0.4SPRI_Histag.tar -O ERX708230.tar
wget ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_39.tar -O ERX708231.tar

# Untar into directories named x.fast5 and move the files into the top directory
set_name=ERX708228
mkdir -p ${SET_NAME}.fast5
cd ${SET_NAME}.fast5; tar -xf ../${SET_NAME}.tar
find ${SET_NAME}.fast5 -name "*.fast5" -exec mv {} ${SET_NAME}.fast5 \;

set_name=ERX708229
mkdir -p ${SET_NAME}.fast5
cd ${SET_NAME}.fast5; tar -xf ../${SET_NAME}.tar
find ${SET_NAME}.fast5 -name "*.fast5" -exec mv {} ${SET_NAME}.fast5 \;

set_name=ERX708230
mkdir -p ${SET_NAME}.fast5
cd ${SET_NAME}.fast5; tar -xf ../${SET_NAME}.tar
find ${SET_NAME}.fast5 -name "*.fast5" -exec mv {} ${SET_NAME}.fast5 \;

set_name=ERX708231
mkdir -p ${SET_NAME}.fast5
cd ${SET_NAME}.fast5; tar -xf ../${SET_NAME}.tar
find ${SET_NAME}.fast5 -name "*.fast5" -exec mv {} ${SET_NAME}.fast5 \;

# Export 2D reads to fasta files using poretools
poretools fasta --type 2D ERX708228.fast5/ > raw.reads.unsorted
poretools fasta --type 2D ERX708229.fast5/ >> raw.reads.unsorted
poretools fasta --type 2D ERX708230.fast5/ >> raw.reads.unsorted
poretools fasta --type 2D ERX708231.fast5/ >> raw.reads.unsorted
python lengthsort.py < raw.reads.unsorted > raw.reads.fasta



##################################################
#
# Step 3. Correct the raw nanopore data using
# nanocorrect. This step takes a long time.
#
##################################################

# Error correction, the first step.
make -f nanocorrect/nanocorrect-overlap.make INPUT=raw.reads.fasta NAME=raw.reads
samtools faidx raw.reads.pp.fasta
python nanocorrect/makerange.py raw.reads.fasta | parallel -v --eta -P $(NC_PROCESS) 'python nanocorrect/nanocorrect.py raw.reads {} > raw.reads.{}.corrected.fasta'
cat raw.reads.*.corrected.fasta | python lengthsort.py > raw.reads.corrected.fasta
#rm raw.reads.*.corrected.fasta

# Error correction, the second step.
make -f nanocorrect/nanocorrect-overlap.make INPUT=raw.reads.corrected.fasta NAME=raw.reads.corrected
samtools faidx raw.reads.corrected.pp.fasta
python nanocorrect/makerange.py raw.reads.corrected.fasta | parallel -v --eta -P $(NC_PROCESS) 'python nanocorrect/nanocorrect.py raw.reads.corrected {} > raw.reads.corrected.{}.corrected.fasta'
cat raw.reads.corrected.*.corrected.fasta | python lengthsort.py > raw.reads.corrected.corrected.fasta
#rm raw.reads.corrected.*.corrected.fasta



##################################################
#
# Step 4. Run the celera assembler on the 
# corrected reads.
#
##################################################
# prepare the input into celera assembler. 
# we want to use the twice-corrected data here, so two prereqs
java -Xmx1024M -jar ./wgs-8.2/Linux-amd64/bin/convertFastaAndQualToFastq.jar raw.reads.corrected.corrected.fasta > assembly.input.fastq
fastqToCA -technology sanger -libraryname assembly -reads assembly.input.fastq > assembly.frg

# Download the spec file we use
wget --no-check-certificate https://raw.githubusercontent.com/jts/nanopore-paper-analysis/c25373d93a99e51c2fedb57d8b08b81826e7c80c/revised_ovlErrorRate0.04.spec

# Run the assembly
# Resulting scaffolds will be in: celera-assembly/9-terminator/asm.scf.fasta
runCA -d celera-assembly -p asm -s revised_ovlErrorRate0.04.spec assembly.frg
ln -s celera-assembly/9-terminator/asm.scf.fasta draft_genome.fasta

##################################################
#
# Step 5. Polish the draft assembly with nanopolish
#
##################################################

# preprocess the fasta file for nanopolish
nanopolish/consensus-preprocess.pl raw.reads.fasta > raw.reads.np.fasta

# index the draft assembly for bwa
bwa index draft_genome.fasta

# index the draft assembly for faidx
samtools faidx draft_genome.fasta

# align reads to draft assembly
bwa mem -t $(THREADS) -x ont2d draft_genome.fasta raw.reads.np.fasta | samtools view -Sb - | samtools sort -f - reads_to_draft.sorted.bam

# index the bam file
samtools index reads_to_draft.sorted.bam

# run nanopolish
python nanopolish/nanopolish_makerange.py draft_genome.fasta | parallel --progress -P $(NP_PROCESS) \
	nanopolish/nanopolish consensus -o nanopolish.{1}.fa -r raw.reads.np.fasta -b reads_to_draft.sorted.bam -g draft_genome.fasta -w {1} -t $(THREADS)
python nanopolish/nanopolish_merge.py draft_genome.fasta nanopolish.scf*.fa > polished_genome.fasta
