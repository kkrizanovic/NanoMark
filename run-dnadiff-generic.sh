#! /bin/sh

tool=$1
# MiniasmRacon-newedlib
dataroot=results
indata=data/
memtimeprefix=$2
# MiniasmRacon-NewEdlib

#######################

dataset=lambda_ont-30x
ddir=$dataroot/$dataset/$tool
ref=$indata/$dataset/lambda-NC_001416.fa
mkdir -p $ddir/dnadiff

iter=consensus-iter1.fasta
MUMmer3.23/dnadiff -p $ddir/dnadiff/$iter $ref $ddir/${iter}
src/collect_memtime.py number $ddir $ddir/total-iter1.memtime $memtimeprefix 1 4

iter=consensus-iter2.fasta
MUMmer3.23/dnadiff -p $ddir/dnadiff/$iter $ref $ddir/${iter}
src/collect_memtime.py number $ddir $ddir/total-iter2.memtime $memtimeprefix 1 6

#######################

dataset=ecoli_map006_ont
ddir=$dataroot/$dataset/$tool
ref=$indata/$dataset/ecoli_K12_MG1655_U00096.3.fasta
mkdir -p $ddir/dnadiff

iter=consensus-iter1.fasta
MUMmer3.23/dnadiff -p $ddir/dnadiff/$iter $ref $ddir/${iter}
src/collect_memtime.py number $ddir $ddir/total-iter1.memtime $memtimeprefix 1 4

iter=consensus-iter2.fasta
MUMmer3.23/dnadiff -p $ddir/dnadiff/$iter $ref $ddir/${iter}
src/collect_memtime.py number $ddir $ddir/total-iter2.memtime $memtimeprefix 1 6

#######################

dataset=ecoli_pacbio_p6c4
ddir=$dataroot/$dataset/$tool
ref=$indata/$dataset/ecoli_K12_MG1655_U00096.3.fasta
mkdir -p $ddir/dnadiff

iter=consensus-iter1.fasta
MUMmer3.23/dnadiff -p $ddir/dnadiff/$iter $ref $ddir/${iter}
src/collect_memtime.py number $ddir $ddir/total-iter1.memtime $memtimeprefix 1 4

iter=consensus-iter2.fasta
MUMmer3.23/dnadiff -p $ddir/dnadiff/$iter $ref $ddir/${iter}
src/collect_memtime.py number $ddir $ddir/total-iter2.memtime $memtimeprefix 1 6

#######################

#exit

dataset=scerevisiae_w303
ddir=$dataroot/$dataset/$tool
ref=$indata/$dataset/saccharomyces_cerevisiae-S288C-2016.fa
mkdir -p $ddir/dnadiff-2016

iter=consensus-iter1.fasta
MUMmer3.23/dnadiff -p $ddir/dnadiff-2016/$iter $ref $ddir/${iter}
src/collect_memtime.py number $ddir $ddir/total-iter1.memtime $memtimeprefix 1 4

iter=consensus-iter2.fasta
MUMmer3.23/dnadiff -p $ddir/dnadiff-2016/$iter $ref $ddir/${iter}
src/collect_memtime.py number $ddir $ddir/total-iter2.memtime $memtimeprefix 1 6

#######################

dataset=scerevisiae_ont_r9
ddir=$dataroot/$dataset/${tool}-qv20
ref=$indata/$dataset/saccharomyces_cerevisiae-S288C-2016.fa
mkdir -p $ddir/dnadiff-2016

iter=consensus-iter1.fasta
MUMmer3.23/dnadiff -p $ddir/dnadiff-2016/$iter $ref $ddir/${iter}
src/collect_memtime.py number $ddir $ddir/total-iter1.memtime $memtimeprefix 1 4

iter=consensus-iter2.fasta
MUMmer3.23/dnadiff -p $ddir/dnadiff-2016/$iter $ref $ddir/${iter}
src/collect_memtime.py number $ddir $ddir/total-iter2.memtime $memtimeprefix 1 6

#######################

# This exit is here to prevent running Dnadiff on C. Elegans automatically
# because it takes about 24h to complete, each.
exit

dataset=celegans
ddir=$dataroot/$dataset/$tool
ref=$indata/$dataset/caenorhabditis_elegans.fa
mkdir -p $ddir/dnadiff

iter=consensus-iter1.fasta
src/collect_memtime.py number $ddir $ddir/total-iter1.memtime $memtimeprefix 1 4
MUMmer3.23/dnadiff -p $ddir/dnadiff/$iter $ref $ddir/${iter}

iter=consensus-iter2.fasta
src/collect_memtime.py number $ddir $ddir/total-iter2.memtime $memtimeprefix 1 6
MUMmer3.23/dnadiff -p $ddir/dnadiff/$iter $ref $ddir/${iter}
