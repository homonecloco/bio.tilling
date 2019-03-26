

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

	innerDf<-merge(tmpDf2, tmpDf, by.x="exon", by.y="Var1")
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
	
	gg<- gg+geom_line(data=lines_plot_data, aes(x=factor(Start), y=value,group=Library, colour=Library))
	gg<- gg+geom_line(data=lines_plot_data, aes(x=factor(Start), y=MaxForHet,group="Max. for Het", colour="Max. for Het"))
	gg<- gg + geom_point(data=plot_Data[plot_Data$Library == library,], aes(x=factor(Start), y=value, col=factor(ValueType)))
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

prepareForHeatmap<-function(tableWithAllDels, 
                            column='hetDelsPer',
                            sortBySum=False
                           ){
    
    names_hm<-paste(tableWithAllDels$chr,tableWithAllDels$cM, sep="_")
    preHMHet<-cbind(names_hm,tableWithAllDels$Library, tableWithAllDels[column]) 

    colnames(preHMHet)<-c("cM","Library","Percentage")
    preHMHet<-data.frame(preHMHet)
    preHMHet$Percentage <- as.numeric(as.character(preHMHet$Percentage))
    preHMreshaped<-reshape(preHMHet, idvar=c("cM"), timevar="Library", direction="wide")
    rownames(preHMreshaped) <- preHMreshaped[,1]
    preHMreshaped[,1] <- NULL
    hmMat<-as.matrix(preHMreshaped)
    if(sortBySum){
        hmMat<-hmMat[,order(colSums(hmMat))]
    }
    hmMat

}

plotDeletionsHeatmap <-function(geneticMap,selectedDels, libraries, minExonCount=4, prefix="heatmap"){
	library(gplots)
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


plotPerChromosome<-function(geneticMap,selectedDels, libraries,df, minExonCount=9,
                       groupByCM=T, chr='1D',column='homDelsPer', wtGrep='WT', 
                       parseLibFun=identity,  groupByCMPrecision=0){
    library(gplots)
    dels<-getDeletionsPerCM(geneticMap,selectedDels, libraries,df, minExonCount=minExonCount,
                       groupByCM=groupByCM, chr=chr, groupByCMPrecision=groupByCMPrecision)
    hmMat<-prepareForHeatmap(dels, column, sortBySum=TRUE)  
    colnames(hmMat)<-sapply( colnames(hmMat), parseLibFun)
    libNames<-colnames(hmMat)
    print(libNames)
    wtIndeces<-grep(wtGrep, libNames, fixed=T)
    mutIndeces<-grep(wtGrep, libNames, fixed=T, invert=T)
    wtHm<-hmMat[,wtIndeces]
    mutHm<-hmMat[,mutIndeces]
    
    mutHm<-mutHm[,colSums(mutHm) > 0]

    merged<-cbind(mutHm,wtHm )
    wtCount=ncol(mutHm) + 0.5
    my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)    
    merged<-log((merged + 0.0005) *100)
    #merged<-ifelse(merged == -Inf, 0, merged)
    code<-paste0("heatmap.2(merged,dendrogram='none' ,keysize=1, Rowv=NA, Colv=NA, 
              col = colorRampPalette(c('#ffffcc','#a1dab4','#41b6c4','#2c7fb8','#253494'))(n = 299),
              scale='none', margins=c(10,10),
              add.expr = abline(v = ",wtCount ," , col = 'black',lwd=3),
              trace='none')")
    eval(parse(text=code))
}


#IWGSC_CSS_1AL_scaff_875531
plotScaffoldDeletions<- function(df, mat,  contig, samplesSD ){
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
	gg<- ggplot(plot_Data, aes(x=factor(Start), y=value))
	gg<- gg + geom_errorbar(aes(ymax=1, ymin=1-3* sdExon), colour='gray', alpha=0.75)
	gg<- gg + geom_boxplot()
	gg<- gg + geom_point(data=plot_Data[plot_Data$value > plot_Data$upper.limit | plot_Data$value < plot_Data$lower.limit,], aes(x=factor(Start), y=value, col=factor(ValueType)))
	gg<- gg + coord_cartesian(ylim=c(0, 1))
	gg<- gg + geom_text(aes(label=ifelse(ValueType=="4. Deletion 4 sigma exon" , sampleName, ifelse(ValueType=="3. Deletion 3 sigma exon", sampleName, '')))  ,hjust=0,just=0, size=3, angle = 75)
	gg<- gg + ggtitle(paste("Normalized coverage for contig \n ", contig))
	gg
}




plotScaffoldDeletionsInLibrary<- function(df, mat,  contig, samplesSD ){
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
	innerDf<-merge(tmpDf2, tmpDf, by.x="exon", by.y="Var1")
	innerDf<-merge(innerDf, libDf, by.x="Var2", by.y="Library")
	innerDf$ValueType <- mapply(getValueType,innerDf$value, innerDf$sdLib, innerDf$sdExon)


	pre_plot_Data <- ddply(innerDf, .(Start), mutate, Q1=quantile(value, 1/4), Q3=quantile(value, 3/4), IQR=Q3-Q1, upper.limit=Q3+1.5*IQR, lower.limit=Q1-1.5*IQR)

	libraries<-subset(pre_plot_Data, grepl("Del", ValueType))

	libNames<-unique(libraries$sampleName)
	plot_data<-subset(pre_plot_Data, sampleName %in% libNames )
	plot_data$Library <- plot_data$Var2
	gg<- ggplot(plot_data, aes(x=Start, y=value))
	gg<- gg + geom_errorbar(aes(ymax=1, ymin=1-3* sdExon), colour='gray', alpha=0.75) 
	gg<- gg+geom_line(aes(linetype=Library))
	gg<-gg+geom_point(aes(colour=ValueType))
	gg<- gg +coord_cartesian(ylim=c(0, 1))
	gg<- gg + ggtitle(paste("Libraries with deletions in \n ", contig))
	gg
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plotLibraryOnChromosomeInterval <- function(contig, library, df, mat, samplesSD, from, to  ){
	library(ggplot2)
	tempExons<-rownames(subset(df,Scaffold==contig & Start > from & Ends < to))
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
	

	innerDf<-merge(tmpDf2, tmpDf, by.x="exon", by.y="Var1")
	innerDf<-merge(innerDf, libDf, by.x="Var2", by.y="Library")
	innerDf$ValueType <- mapply(getValueType,innerDf$value, innerDf$sdLib, innerDf$sdExon)

	plot_Data <- ddply(innerDf, .(Start), mutate, Q1=quantile(value, 1/4), Q3=quantile(value, 3/4), IQR=Q3-Q1, upper.limit=Q3+1.5*IQR, lower.limit=Q1-1.5*IQR)
	plot_Data$Library <- plot_Data$Var2	
	libNames <-c(library)
	lines_plot_data<-subset(plot_Data, Library %in% library )
	lines_plot_data$MaxForHet<-0.5+lines_plot_data$sdExon
	gg <- ggplot(plot_Data, aes(x=factor(Start), y=value))
	gg <- gg + geom_errorbar(aes(ymax=1+3* sdExon, ymin=1-3* sdExon), colour='gray', alpha=0.75) 
	gg <- gg + geom_boxplot(outlier.shape = NA)
	gg <- gg + geom_line(data=lines_plot_data, aes(x=factor(Start), y=value,group=Library, colour=Library))
	gg <- gg + geom_line(data=lines_plot_data, aes(x=factor(Start), y=MaxForHet,group="Max. for Het", colour="Max. for Het"))
	gg <- gg + geom_point(data=plot_Data[plot_Data$Library == library,], aes(x=factor(Start), y=value, col=factor(ValueType)))
	gg <- gg + ggtitle(paste("Normalised coverage\n", contig,"\n",from, "-", to))
    gg <- gg + theme_bw()
	gg
}
