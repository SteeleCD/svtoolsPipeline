#!/bin/bash
# input details
dataDir="/Volumes/grp11$/CI_Pathology_Steele/data/genetic/Sanger/Behjati-2017"
segFile=$dataDir/ncomms15936-s8-seg.csv
bedpeFile=$dataDir/ncomms15936-s6-bedpe.csv

# running details
sample="PD13494a"
chrom="chr17"

# out details
outDir="/Users/c.steele/Dropbox/PostDoc/projects/1_undiff/results/WGS/chromothripsis/svtools"
outSeg=$sample"-seg.txt"
outBedpe=$sample"-bedpe.txt"

# pipeline details
svtoolsDir="/Users/c.steele/Dropbox/PostDoc/Rlibs/svtoolsPipeline/"
pythonDir=$svtoolsDir/py

# data preprocessing
Rscript --vanilla $svtoolsDir/R/dataPreprocess.R \
$segFile $bedpeFile $sample $outDir $outSeg $outBedpe

# SV diagram
python $pythonDir/runSVdiagram.py \
	$outDir/$outSeg \
	$outDir/$outBedpe \
	$chrom \
	$outDir/$sample"-"$chrom"-SVdiagram.pdf"
	
# KC tests
python $pythonDir/runKCtests.py \
	$outDir/$outBedpe \
	$chrom \
	$outDir/$sample"-"$chrom"-"



