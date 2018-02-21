#!/usr/bin/env Rscript
# get command line arguments
args = commandArgs(trailingOnly=TRUE)
segFile=args[1]
bedpeFile=args[2]
sample=args[3]
outDir=args[4]
outSeg=args[5]
outBedpe=args[6]
if(length(args)>6)
	{
	segSampleCol=as.numeric(args[7])
	segChromCol=as.numeric(args[8])
	segStartCol=as.numeric(args[9])
	segEndCol=as.numeric(args[10])
	segCNcol=as.numeric(args[11])
	bedpeSampleCol=as.numeric(args[12])
	bedpeChromCol1=as.numeric(args[13])
	bedpePosCol1=as.numeric(args[14])
	bedpeStrandCol1=as.numeric(args[15])
	bedpeChromCol2=as.numeric(args[16])
	bedpePosCol2=as.numeric(args[17])
	bedpeStrandCol2=as.numeric(args[18])
	segHead=as.logical(as.numeric(args[19]))
	bedpeHead=as.logical(as.numeric(args[20]))
	}

# read in a file
readFile = function(file,head=TRUE)
	{
	ending = rev(strsplit(file,split="[.]")[[1]])[1]
	if(ending=="csv") return(read.csv(file,head=head))
	if(ending%in%c("txt","tsv")) return(read.table(file,sep="\t",head=head,as.is=TRUE))
	return(read.table(file,sep="\t",head=head,as.is=TRUE))	
	}

# split seg file
splitSeg = function(seg,sampleCol,chromCol,startCol,endCol,CNcol,windowSize=5000)
	{
	library(plyr)
	starts = apply(seg,MARGIN=1,FUN=function(x) 
		{
		x=gsub(" ","",unlist(x))
		print(x)
		round_any(seq(from=as.numeric(x[startCol]),to=as.numeric(x[endCol]),by=windowSize),1000)
		})
	
	outSeg = sapply(1:length(starts),FUN=function(i) cbind(paste0(seg[i,sampleCol]),
						paste0(seg[i,chromCol]),
						starts[[i]],
						starts[[i]]+windowSize,
						seg[i,CNcol]))
	outSeg = do.call(rbind,outSeg)
	return(outSeg)
}

# preprocess data for SV tools
preprocessSVtools = function(segFile,bedpeFile,sample,outDir=".",
  outSeg="testSeg.txt",outBedpe="testBedpe.txt",
	segSampleCol=1,segChromCol=2,segStartCol=3,segEndCol=4,segCNcol=5,
	bedpeChromCol1=2,bedpePosCol1=4,bedpeStrandCol1=3,
	bedpeChromCol2=5,bedpePosCol2=7,bedpeStrandCol2=6,
	bedpeReadsCol=NULL,bedpeSampleCol=1,segHead=TRUE,bedpeHead=TRUE
	)
	{
	# seg file
	seg = readFile(segFile,head=segHead)
	seg = seg[which(seg[,segSampleCol]==sample),,drop=FALSE]
	if(!any(grepl("chr",seg[,segChromCol])))
		{
		seg[,segChromCol]=paste0("chr",seg[,segChromCol])
		}
	seg = splitSeg(seg,sampleCol=segSampleCol,chromCol=segChromCol,
			startCol=segStartCol,endCol=segEndCol,CNcol=segCNcol)
	seg = cbind(seg[,2:4],".",seg[,5])
	write.table(seg,paste0(outDir,"/",outSeg),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	# bedpefile
	bedpe = readFile(bedpeFile,head=bedpeHead)
	if(!is.null(bedpeSampleCol)) bedpe = bedpe[which(bedpe[,bedpeSampleCol]==sample),,drop=FALSE]
	if(!any(grepl("chr",bedpe[,bedpeChromCol1])))
		{
		bedpe[,bedpeChromCol1]=paste0("chr",bedpe[,bedpeChromCol1])
		bedpe[,bedpeChromCol2]=paste0("chr",bedpe[,bedpeChromCol2])
		}
	if(is.null(bedpeReadsCol))
		{
		reads = rep(20,nrow(bedpe))
		} 
	bedpe = cbind(bedpe[,c(bedpeChromCol1,bedpePosCol1,bedpeStrandCol1,
				bedpeChromCol2,bedpePosCol2,bedpeStrandCol2)],
			reads,0)
	colnames(bedpe) = c("chrom1","pos1","strand1","chrom2","pos2","strand2","reads","gap")
	write.table(bedpe,paste0(outDir,"/",outBedpe),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
  }

# run preprocessing
preprocessSVtools(segFile,bedpeFile,sample,outDir,
                outSeg,outBedpe,
		segSampleCol=segSampleCol,
		segChromCol=segChromCol,
		segStartCol=segStartCol,
		segEndCol=segEndCol,
		segCNcol=segCNcol,
		bedpeSampleCol=bedpeSampleCol,
		bedpeChromCol1=bedpeChromCol1,
		bedpePosCol1=bedpePosCol1,
		bedpeStrandCol1=bedpeStrandCol1,
		bedpeChromCol2=bedpeChromCol2,
		bedpePosCol2=bedpePosCol2,
		bedpeStrandCol2,
		segHead=segHead,
		bedpeHead=bedpeHead)


#segFile=segFile
#bedpeFile=bedpeFile
#sample=sample
#outDir=outDir
#outSeg=outSeg
#outBedpe=outBedpe
#		segSampleCol=segSampleCol
#		segChromCol=segChromCol
#		segStartCol=segStartCol
#		segEndCol=segEndCol
#		segCNcol=segCNcol
#		bedpeSampleCol=bedpeSampleCol
#		bedpeChromCol1=bedpeChromCol1
#		bedpePosCol1=bedpePosCol1
#		bedpeStrandCol1=bedpeStrandCol1
#		bedpeChromCol2=bedpeChromCol2
#		bedpePosCol2=bedpePosCol2
#		bedpeStrandCol2
#		segHead=segHead
#		bedpeHead=bedpeHead
