#! /bin/sh

dd=MUMmer3.23/dnadiff

#######################

dataset=scerevisiae_ont_r9
ref=data/$dataset/saccharomyces_cerevisiae-S288C-2016.fa
dddir=results/$dataset/MiniasmSparc

iter=consensus-iter1.fasta
mkdir -p $dddir/dnadiff-2016
$dd -p $dddir/dnadiff-2016/$iter $ref $dddir/$iter
src/collect_memtime.py number $dddir $dddir/total-iter1.memtime MiniasmSparc 1 33

iter=consensus-iter2.fasta
mkdir -p $dddir/dnadiff-2016
$dd -p $dddir/dnadiff-2016/$iter $ref $dddir/$iter
src/collect_memtime.py number $dddir $dddir/total-iter2.memtime MiniasmSparc 1 64

src/collect_memtime.py pattern $dddir $dddir/total-new.memtime MiniasmSparc-*

exit

#######################

dataset=scerevisiae_w303
ref=data/$dataset/saccharomyces_cerevisiae-S288C-2016.fa
dddir=results/$dataset/MiniasmSparc

iter=consensus-iter1.fasta
mkdir -p $dddir/dnadiff-2016
$dd -p $dddir/dnadiff-2016/$iter $ref $dddir/$iter
src/collect_memtime.py number $dddir $dddir/total-iter1.memtime MiniasmSparc 1 33

iter=consensus-iter2.fasta
mkdir -p $dddir/dnadiff-2016
$dd -p $dddir/dnadiff-2016/$iter $ref $dddir/$iter
src/collect_memtime.py number $dddir $dddir/total-iter2.memtime MiniasmSparc 1 64

src/collect_memtime.py pattern $dddir $dddir/total-new.memtime MiniasmSparc-*

exit

#######################

dataset=ecoli_pacbio_p6c4
ref=data/$dataset/ecoli_K12_MG1655_U00096.3.fasta
dddir=results/$dataset/MiniasmSparc

iter=consensus-iter1.fasta
mkdir -p $dddir/dnadiff
$dd -p $dddir/dnadiff/$iter $ref $dddir/$iter
src/collect_memtime.py number $dddir $dddir/total-iter1.memtime MiniasmSparc 1 4

iter=consensus-iter2.fasta
mkdir -p $dddir/dnadiff
$dd -p $dddir/dnadiff/$iter $ref $dddir/$iter
src/collect_memtime.py number $dddir $dddir/total-iter2.memtime MiniasmSparc 1 6

src/collect_memtime.py pattern $dddir $dddir/total-new.memtime MiniasmSparc-*

exit

#######################

dataset=ecoli_map006_ont
ref=data/$dataset/ecoli_K12_MG1655_U00096.3.fasta
dddir=results/$dataset/MiniasmSparc

iter=consensus-iter1.fasta
mkdir -p $dddir/dnadiff
$dd -p $dddir/dnadiff/$iter $ref $dddir/$iter
src/collect_memtime.py number $dddir $dddir/total-iter1.memtime MiniasmSparc 1 4

iter=consensus-iter2.fasta
mkdir -p $dddir/dnadiff
$dd -p $dddir/dnadiff/$iter $ref $dddir/$iter
src/collect_memtime.py number $dddir $dddir/total-iter2.memtime MiniasmSparc 1 6

src/collect_memtime.py pattern $dddir $dddir/total-new.memtime MiniasmSparc-*

exit

#######################
dataset=lambda_ont-30x
ref=data/$dataset/lambda-NC_001416.fa
dddir=results/$dataset/MiniasmSparc

iter=consensus-iter1.fasta
mkdir -p $dddir/dnadiff
$dd -p $dddir/dnadiff/$iter $ref $dddir/$iter
src/collect_memtime.py number $dddir $dddir/total-iter1.memtime MiniasmSparc 1 4

iter=consensus-iter2.fasta
mkdir -p $dddir/dnadiff
$dd -p $dddir/dnadiff/$iter $ref $dddir/$iter
src/collect_memtime.py number $dddir $dddir/total-iter2.memtime MiniasmSparc 1 6

src/collect_memtime.py pattern $dddir $dddir/total-new.memtime MiniasmSparc-*

exit
#######################
