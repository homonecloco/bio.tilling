#!/bin/bash

inputFolder=$1
#References/exons.bed
exonsFile=$2
outputFolder=$3

#source bedtools-2.25.0
mkdir -p $outputFolder
for f in `ls $inputFolder | grep -v bai`; do
    echo $f
    filename="${f%.*}"
    if [ -f $outputFolder/$filename.tab ] ; then
        echo "Already extracted"
    else
        #This line can be paralelized. 
        coverageBed -a $exonsFile -b $inputFolder/$f -sorted -g References/chr_order.txt  > $outputFolder/$filename.bedCov 
    fi
done

folder="$outputFolder"
cd $folder

for fullfile in `ls *.bedCov`; do
    filename=$(basename "$fullfile")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo $filename > cov.${filename}.txt
    awk '{print $4}' $fullfile >> cov.${filename}.txt
done


echo "" > exonNames.txt
awk '{print $1":"$2":"$3}' $fullfile >> exonNames.txt
head exonNames.txt
#This line may be needed if your user can't have many file opens. This is because to paste all the lines we have all the bed files open at the same time. 
ulimit -S -n 2048
paste exonNames.txt cov*.txt > allMergedCoverages.tab
