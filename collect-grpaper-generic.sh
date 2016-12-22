#! /bin/sh

tool=$1
# MiniasmRacon-newedlib
dataroot=results
indata=data/

totalresultsiter1=results-$tool-iter1.txt
totalresultsiter2=results-$tool-iter2.txt

mkdir -p temp

echo "" > $totalresultsiter1
echo "" > $totalresultsiter2

#######################

dataset=lambda_ont-30x
ddir=$dataroot/$dataset/$tool

iter=consensus-iter1.fasta
grep -P "TotalSeqs|TotalBases|AvgIdentity|AlignedBases" $ddir/dnadiff/$iter.report 2>&1 | tee temp/temp_collect-dnadiff.txt
tail -n 5 $ddir/total-iter1.memtime 2>&1 | tee temp/temp_collect-memtime.txt
echo $dataset >> $totalresultsiter1
echo $ddir >> $totalresultsiter1
cat temp/temp_collect-dnadiff.txt >> $totalresultsiter1
cat temp/temp_collect-memtime.txt >> $totalresultsiter1
echo "" >> $totalresultsiter1

iter=consensus-iter2.fasta
grep -P "TotalSeqs|TotalBases|AvgIdentity|AlignedBases" $ddir/dnadiff/$iter.report 2>&1 | tee temp/temp_collect-dnadiff.txt
tail -n 5 $ddir/total-iter2.memtime 2>&1 | tee temp/temp_collect-memtime.txt
echo $dataset >> $totalresultsiter2
echo $ddir >> $totalresultsiter2
cat temp/temp_collect-dnadiff.txt >> $totalresultsiter2
cat temp/temp_collect-memtime.txt >> $totalresultsiter2
echo "" >> $totalresultsiter2

#######################

dataset=ecoli_map006_ont
ddir=$dataroot/$dataset/$tool

iter=consensus-iter1.fasta
grep -P "TotalSeqs|TotalBases|AvgIdentity|AlignedBases" $ddir/dnadiff/$iter.report 2>&1 | tee temp/temp_collect-dnadiff.txt
tail -n 5 $ddir/total-iter1.memtime 2>&1 | tee temp/temp_collect-memtime.txt
echo $dataset >> $totalresultsiter1
echo $ddir >> $totalresultsiter1
cat temp/temp_collect-dnadiff.txt >> $totalresultsiter1
cat temp/temp_collect-memtime.txt >> $totalresultsiter1
echo "" >> $totalresultsiter1

iter=consensus-iter2.fasta
grep -P "TotalSeqs|TotalBases|AvgIdentity|AlignedBases" $ddir/dnadiff/$iter.report 2>&1 | tee temp/temp_collect-dnadiff.txt
tail -n 5 $ddir/total-iter2.memtime 2>&1 | tee temp/temp_collect-memtime.txt
echo $dataset >> $totalresultsiter2
echo $ddir >> $totalresultsiter2
cat temp/temp_collect-dnadiff.txt >> $totalresultsiter2
cat temp/temp_collect-memtime.txt >> $totalresultsiter2
echo "" >> $totalresultsiter2

#######################

dataset=ecoli_pacbio_p6c4
ddir=$dataroot/$dataset/$tool

iter=consensus-iter1.fasta
grep -P "TotalSeqs|TotalBases|AvgIdentity|AlignedBases" $ddir/dnadiff/$iter.report 2>&1 | tee temp/temp_collect-dnadiff.txt
tail -n 5 $ddir/total-iter1.memtime 2>&1 | tee temp/temp_collect-memtime.txt
echo $dataset >> $totalresultsiter1
echo $ddir >> $totalresultsiter1
cat temp/temp_collect-dnadiff.txt >> $totalresultsiter1
cat temp/temp_collect-memtime.txt >> $totalresultsiter1
echo "" >> $totalresultsiter1

iter=consensus-iter2.fasta
grep -P "TotalSeqs|TotalBases|AvgIdentity|AlignedBases" $ddir/dnadiff/$iter.report 2>&1 | tee temp/temp_collect-dnadiff.txt
tail -n 5 $ddir/total-iter2.memtime 2>&1 | tee temp/temp_collect-memtime.txt
echo $dataset >> $totalresultsiter2
echo $ddir >> $totalresultsiter2
cat temp/temp_collect-dnadiff.txt >> $totalresultsiter2
cat temp/temp_collect-memtime.txt >> $totalresultsiter2
echo "" >> $totalresultsiter2

#######################

dataset=scerevisiae_w303
ddir=$dataroot/$dataset/$tool

iter=consensus-iter1.fasta
grep -P "TotalSeqs|TotalBases|AvgIdentity|AlignedBases" $ddir/dnadiff-2016/$iter.report 2>&1 | tee temp/temp_collect-dnadiff.txt
tail -n 5 $ddir/total-iter1.memtime 2>&1 | tee temp/temp_collect-memtime.txt
echo $dataset >> $totalresultsiter1
echo $ddir >> $totalresultsiter1
cat temp/temp_collect-dnadiff.txt >> $totalresultsiter1
cat temp/temp_collect-memtime.txt >> $totalresultsiter1
echo "" >> $totalresultsiter1

iter=consensus-iter2.fasta
grep -P "TotalSeqs|TotalBases|AvgIdentity|AlignedBases" $ddir/dnadiff-2016/$iter.report 2>&1 | tee temp/temp_collect-dnadiff.txt
tail -n 5 $ddir/total-iter2.memtime 2>&1 | tee temp/temp_collect-memtime.txt
echo $dataset >> $totalresultsiter2
echo $ddir >> $totalresultsiter2
cat temp/temp_collect-dnadiff.txt >> $totalresultsiter2
cat temp/temp_collect-memtime.txt >> $totalresultsiter2
echo "" >> $totalresultsiter2

#######################

dataset=scerevisiae_ont_r9
ddir=$dataroot/$dataset/${tool}-qv20

iter=consensus-iter1.fasta
grep -P "TotalSeqs|TotalBases|AvgIdentity|AlignedBases" $ddir/dnadiff-2016/$iter.report 2>&1 | tee temp/temp_collect-dnadiff.txt
tail -n 5 $ddir/total-iter1.memtime 2>&1 | tee temp/temp_collect-memtime.txt
echo $dataset >> $totalresultsiter1
echo $ddir >> $totalresultsiter1
cat temp/temp_collect-dnadiff.txt >> $totalresultsiter1
cat temp/temp_collect-memtime.txt >> $totalresultsiter1
echo "" >> $totalresultsiter1

iter=consensus-iter2.fasta
grep -P "TotalSeqs|TotalBases|AvgIdentity|AlignedBases" $ddir/dnadiff-2016/$iter.report 2>&1 | tee temp/temp_collect-dnadiff.txt
tail -n 5 $ddir/total-iter2.memtime 2>&1 | tee temp/temp_collect-memtime.txt
echo $dataset >> $totalresultsiter2
echo $ddir >> $totalresultsiter2
cat temp/temp_collect-dnadiff.txt >> $totalresultsiter2
cat temp/temp_collect-memtime.txt >> $totalresultsiter2
echo "" >> $totalresultsiter2

# exit

#######################

dataset=celegans
ddir=$dataroot/$dataset/$tool

iter=consensus-iter1.fasta
grep -P "TotalSeqs|TotalBases|AvgIdentity|AlignedBases" $ddir/dnadiff/$iter.report 2>&1 | tee temp/temp_collect-dnadiff.txt
tail -n 5 $ddir/total-iter1.memtime 2>&1 | tee temp/temp_collect-memtime.txt
echo $dataset >> $totalresultsiter1
echo $ddir >> $totalresultsiter1
cat temp/temp_collect-dnadiff.txt >> $totalresultsiter1
cat temp/temp_collect-memtime.txt >> $totalresultsiter1
echo "" >> $totalresultsiter1

iter=consensus-iter2.fasta
grep -P "TotalSeqs|TotalBases|AvgIdentity|AlignedBases" $ddir/dnadiff/$iter.report 2>&1 | tee temp/temp_collect-dnadiff.txt
tail -n 5 $ddir/total-iter2.memtime 2>&1 | tee temp/temp_collect-memtime.txt
echo $dataset >> $totalresultsiter2
echo $ddir >> $totalresultsiter2
cat temp/temp_collect-dnadiff.txt >> $totalresultsiter2
cat temp/temp_collect-memtime.txt >> $totalresultsiter2
echo "" >> $totalresultsiter2
