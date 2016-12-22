#! /bin/sh

numthreads=12
reads=data/ecoli_map006_ont/reads.fastq

### Run without overlapping windows
outdir=results-erc/MiniasmRacon-grpaper-erc/
paf=$outdir/reads_to_reads.paf
correctedreads=$outdir/reads.corrected.fasta
memtime_minimap=$outdir/minimap.memtime
memtime_racon=$outdir/racon.memtime

mkdir -p $outdir

/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${memtime_minimap} \
	assemblers/miniasmracon-grpaper/minimap/minimap -L100 -Sw5 -m0 -t $numthreads $reads $reads > $paf

/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${memtime_racon} \
 	assemblers/miniasmracon-grpaper/racon/bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 --erc -t $numthreads $reads $paf $reads $correctedreads

ref=data/ecoli_map006_ont/ecoli_K12_MG1655_U00096.3.fasta
sam=$outdir/graphmap.sam
graphmap/bin/Linux-x64/graphmap align -a anchorgotoh -r $ref -d $correctedreads -o $sam
samscripts/src/errorrates.py base $ref $sam 2>&1 | tee $outdir/eval.txt
samscripts/src/fastqfilter.py info results-erc/MiniasmRacon-grpaper-erc/reads.corrected.fasta data/ecoli_map006_ont/ecoli_K12_MG1655_U00096.3.fasta > results-erc/MiniasmRacon-grpaper-erc/coverage.txt



### Run with overlapping windows
outdir=results-erc/MiniasmRacon-grpaper-erc-ovlwin/
paf=$outdir/reads_to_reads.paf
correctedreads=$outdir/reads.corrected.fasta
memtime_minimap=$outdir/minimap.memtime
memtime_racon=$outdir/racon.memtime

mkdir -p $outdir

/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${memtime_minimap} \
	assemblers/miniasmracon-grpaper/minimap/minimap -L100 -Sw5 -m0 -t $numthreads $reads $reads > $paf

/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${memtime_racon} \
	assemblers/miniasmracon-grpaper/racon/bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 --erc -t $numthreads --ovl-margin 0.1 $reads $paf $reads $correctedreads

ref=data/ecoli_map006_ont/ecoli_K12_MG1655_U00096.3.fasta
sam=$outdir/graphmap.sam
graphmap/bin/Linux-x64/graphmap align -a anchorgotoh -r $ref -d $correctedreads -o $sam
samscripts/src/errorrates.py base $ref $sam 2>&1 | tee $outdir/eval.txt
samscripts/src/fastqfilter.py info results-erc/MiniasmRacon-grpaper-erc-ovlwin/reads.corrected.fasta data/ecoli_map006_ont/ecoli_K12_MG1655_U00096.3.fasta > results-erc/MiniasmRacon-grpaper-erc-ovlwin/coverage.txt
