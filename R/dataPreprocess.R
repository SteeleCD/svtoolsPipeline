preprocessSVtools = function(segFile,bedpeFile,sample,outDir=".",outSeg="testSeg.txt",outBedpe="testBedpe.txt",
	segSampleCol=1,segChromCol=2,segStartCol=3,segEndCol=4,segCNcol=5,
	bedpeChromCol1,bedpePosCol1,bedpeStrandCol1,
	bedpeChromCol2,bedpePosCol2,bedpeStrandCol2,
	bedpeReadsCol=NULL,bedpeSampleCol=NULL
	)
	{
	# seg file
	seg = readFile(segFile)
	seg = seg[which(seg[,segSampleCol]==sample),,drop=FALSE]
	if(!any(grepl("chr",seg[,segChromCol])))
		{
		seg[,segChromCol]=paste0("chr",segChromCol)
		}
	outSeg = cbind(seg[,c(segChromCol,segStartCol,segEndCol)],".",seg[,segCNcol])
	write.table(outSeg,paste0(outDir,"/",outSeg),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	# bedpefile
	bedpe = readFile(bedpeFile)
	if(!is.null(bedpeSampleCol)) bedpe = seg[which(bedpe[,bedpeSampleCol]==sample),,drop=FALSE]
	if(!any(grepl("chr",seg[,bedpeChromCol1])))
		{
		bedpe[,bedpeChromCol1]=paste0("chr",bedpe[,bedpeChromCol1])
		bedpe[,bedpeChromCol2]=paste0("chr",bedpe[,bedpeChromCol2])
		}
	if(is.null(bedpeReadsCol))
		{
		reads = rep(20,nrow(bedpe))
		} 
	outBedpe = cbind(bedpe[,c(bedpeChromCol1,bedpePosCol1,bedpeStrandCol1,
				bedpeChromCol2,bedpePosCol2,bedpeStrandCol2)],
			reads,0)
	colnames(outBedpe) = c("chrom1","pos1","strand1","chrom2","pos2","strand2","reads","gap")
	write.table(outBedpe,paste0(outDir,"/",outBedpe),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
	}
