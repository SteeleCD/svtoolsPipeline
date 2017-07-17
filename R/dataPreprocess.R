# read in a file
readFile = function(file,head=TRUE)
	{
	ending = rev(strsplit(file,split="[.]")[[1]])[1]
	if(ending=="csv") return(read.csv(file,head=head))
	if(ending%in%c("txt","tsv")) return(read.table(file,sep="\t",head=head))
	return(read.table(file,sep="\t",head=head))	
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
preprocessSVtools = function(segFile,bedpeFile,sample,outDir=".",outSeg="testSeg.txt",outBedpe="testBedpe.txt",
	segSampleCol=1,segChromCol=2,segStartCol=3,segEndCol=4,segCNcol=5,
	bedpeChromCol1=2,bedpePosCol1=4,bedpeStrandCol1=3,
	bedpeChromCol2=5,bedpePosCol2=7,bedpeStrandCol2=6,
	bedpeReadsCol=NULL,bedpeSampleCol=NULL,segHead=TRUE,bedpeHead=TRUE
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
	seg = cbind(seg[,c(segChromCol,segStartCol,segEndCol)],".",seg[,segCNcol])
	write.table(seg,paste0(outDir,"/",outSeg),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	# bedpefile
	bedpe = readFile(bedpeFile,head=bedpeHead)
	if(!is.null(bedpeSampleCol)) bedpe = bedpe[which(bedpe[,bedpeSampleCol]==sample),,drop=FALSE]
	if(!any(grepl("chr",seg[,bedpeChromCol1])))
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
