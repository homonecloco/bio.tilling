
plotLibraryOnScaffold <- function(contig, library, df, mat, dels,samplesSD  ){
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
	print(head(lines_plot_data))
	gg<- ggplot(plot_Data, aes(x=factor(Start), y=value))
	gg<- gg + geom_errorbar(aes(ymax=1+3* sdExon, ymin=1-3* sdExon), colour='gray', alpha=0.75) 
	gg<- gg + geom_boxplot(outlier.shape = NA)
	#gg<-gg + stat_smooth() 
	
	gg<- gg+geom_line(data=lines_plot_data, aes(x=factor(Start), y=value,group=Library, colour=Library))
	#gg<-gg+geom_point(aes(colour=ValueType))
	#gg<- gg + coord_cartesian(ylim=c(0, 2))
	#gg<- gg + ggtitle(paste("Libraries with deletions in \n ", contig))
	gg
}