plotDeletionsInScaffold <-function(geneticMapWithDeletions, chromosome, library){
	library(ggplot2)

	toPlot  <-geneticMapWithDeletions[geneticMapWithDeletions$chr == chromosome & geneticMapWithDeletions$Library == library,]
	toPlot <- toPlot[toPlot$ScaffoldCount > 1,]
	hetPlot <- toPlot$hetDels
	homPlot <- toPlot$homDels
	formatedToPlotHet <- data.frame(chr=toPlot$chr, cM=toPlot$cM, type="hetDel", value=hetPlot)
	formatedToPlotHom <- data.frame(chr=toPlot$chr, cM=toPlot$cM, type="homDel", value=homPlot)
	mixPlot<-rbind(formatedToPlotHet, formatedToPlotHom)
	ggplot(data=mixPlot, aes(x=cM, y=value, group=type, colour=type)) +
	geom_line() + geom_point()  +
	scale_y_log10(limits=c(0.001, 1), breaks=c(0.1, 0.25, 0.5, 0.75, 1)) +  
	ylab("Percentage of scaffolds with deletion") +
	xlab("centiMorgan") + 
	ggtitle(paste("Chromosome", chromosome , "\n", library))

}

plotLibraryOnScaffold <- function(contig, library, df, mat, dels, samplesSD, avgs  ){
	library(ggplot2)
	tempExons<-rownames(subset(df,Scaffold==contig))
	tmpMat<-mat[tempExons,]
	melted<-melt(as.matrix(tmpMat))
	tmpDf<-as.data.frame(melted)
	libDf <- as.data.frame( names(samplesSD))
	names(libDf)[1] <- "Library"
	libDf$sdLib <- samplesSD
	libDf$MinSDForDel <- (1-3*samplesSD)
	libDf$sampleName<-sapply(libDf$Library,getSampleName)
	
	tmpDf2<-subset(df,Scaffold==contig)
	tmpDf2$exon <-rownames(tmpDf2)
	#print(colnames(tmpDf2))
	#print(head(tmpDf))

	innerDf<-merge(tmpDf2, tmpDf, by.x="exon", by.y="Var1")
	#print(colnames(innerDf))
	#print(colnames(libDf))
	innerDf<-merge(innerDf, libDf, by.x="Var2", by.y="Library")

	innerDf$ValueType <- mapply(getValueType,innerDf$value, innerDf$sdLib, innerDf$sdExon)

	plot_Data <- ddply(innerDf, .(Start), mutate, Q1=quantile(value, 1/4), Q3=quantile(value, 3/4), IQR=Q3-Q1, upper.limit=Q3+1.5*IQR, lower.limit=Q1-1.5*IQR)
	plot_Data$Library <- plot_Data$Var2	

	libNames <-c(library)
	lines_plot_data<-subset(plot_Data, Library %in% library )
	
	lines_plot_data$MaxForHet<-0.5+lines_plot_data$sdExon
		
	gg<- ggplot(plot_Data, aes(x=factor(Start), y=value))
	gg<- gg + geom_errorbar(aes(ymax=1+3* sdExon, ymin=1-3* sdExon), colour='gray', alpha=0.75) 
	gg<- gg + geom_boxplot(outlier.shape = NA)
	#gg<-gg + stat_smooth() 
	
	gg<- gg+geom_line(data=lines_plot_data, aes(x=factor(Start), y=value,group=Library, colour=Library))
	gg<- gg+geom_line(data=lines_plot_data, aes(x=factor(Start), y=MaxForHet,group="Max. for Het", colour="Max. for Het"))
	#gg<-gg+geom_point(aes(colour=ValueType))
	gg<- gg + geom_point(data=plot_Data[plot_Data$Library == library,], aes(x=factor(Start), y=value, col=factor(ValueType)))
	#gg<- gg + coord_cartesian(ylim=c(0, 2))
	gg<- gg + ggtitle(paste("Normalised coverage\n ", contig))
	gg
}

plotDeletionsInScaffold <-function(geneticMapWithDeletions, chromosome, library){
	library(ggplot2)

	toPlot  <-geneticMapWithDeletions[geneticMapWithDeletions$chr == chromosome & geneticMapWithDeletions$Library == library,]
	toPlot <- toPlot[toPlot$ScaffoldCount > 1,]
	hetPlot <- toPlot$hetDels
	homPlot <- toPlot$homDels
	formatedToPlotHet <- data.frame(chr=toPlot$chr, cM=toPlot$cM, type="hetDel", value=hetPlot)
	formatedToPlotHom <- data.frame(chr=toPlot$chr, cM=toPlot$cM, type="homDel", value=homPlot)
	mixPlot<-rbind(formatedToPlotHet, formatedToPlotHom)
	ggplot(data=mixPlot, aes(x=cM, y=value, group=type, colour=type)) +
	geom_line() + geom_point()  +
	scale_y_log10(limits=c(0.001, 1), breaks=c(0.1, 0.25, 0.5, 0.75, 1)) +  
	ylab("Percentage of scaffolds with deletion") +
	xlab("centiMorgan") + 
	ggtitle(paste("Chromosome", chromosome , "\n", library))
}

prepareForHeatmap<-function(tableWithAllDels, column='hetDelsPer'){
	names_hm<-paste(tableWithAllDels$chr,tableWithAllDels$cm, sep="_")
	preHMHet<-cbind(names_hm,tableWithAllDels$Library, tableWithAllDels[column]) 

	colnames(preHMHet)<-c("cM","Library","Percentage")
	preHMHet<-data.frame(preHMHet)
	preHMHet$Percentage <- as.numeric(as.character(preHMHet$Percentage))
	preHMreshaped<-reshape(preHMHet, idvar=c("cM"), timevar="Library", direction="wide")
	rownames(preHMreshaped) <- preHMreshaped[,1]
	preHMreshaped[,1] <- NULL
	hmMat<-as.matrix(preHMreshaped)
	hmMat

}

plotDeletionsHeatmap <-function(geneticMap,selectedDels, libraries, minExonCount=4, prefix="heatmap"){
	tableWithAllDels <- getDeletionsPerCM(geneticMap,selectedDels, libraries, minExonCount=minExonCount)
	
	
	hmMat<-prepareForHeatmap(tableWithAllDels, column='hetDelsPer')
	filename <- paste0(prefix,"_het.png")
	png(filename = filename, height = 20000, width = 15000)
	heatmap.2(hmMat,dendrogram='none' ,keysize=0.05, Rowv=NA, Colv=NA, col = rich.colors(32), scale="none", margins=c(5,10), trace='none')
	dev.off()
	
	hmMat<-prepareForHeatmap(tableWithAllDels, column='homDelsPer')
	filename <- paste0(prefix,"_hom.png")
	png(filename = filename, height = 20000, width = 15000)
	heatmap.2(hmMat,dendrogram='none' ,keysize=0.05, Rowv=NA, Colv=NA, col = rich.colors(32), scale="none", margins=c(5,10), trace='none')
	dev.off()
}

