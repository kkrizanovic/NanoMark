#! /bin/sh

tool=MiniasmRacon-grpaper-ovlwin
script=wrappers/wrapper_miniasmracon-grpaper.py
runtype=ovlwin

dataset=lambda_ont-30x
$script $runtype results/$dataset/$tool - nanopore,data/$dataset/reads.fastq

dataset=ecoli_pacbio_p6c4
$script $runtype results/$dataset/$tool - pacbio,data/$dataset/reads.fastq

dataset=ecoli_map006_ont
$script $runtype results/$dataset/$tool - nanopore,data/$dataset/reads.fastq

dataset=scerevisiae_w303
$script $runtype results/$dataset/$tool - pacbio,data/$dataset/reads.fastq

dataset=scerevisiae_ont_r9
$script ${runtype}qv20 results/$dataset/${tool}-qv20 - nanopore,data/$dataset/reads.fastq

dataset=celegans
$script $runtype results/$dataset/$tool - pacbio,data/$dataset/reads.fastq
