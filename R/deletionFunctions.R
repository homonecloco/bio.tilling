
getFieldFromExonName <- function(str, field=1,  sp=":") {
	tmp<-unlist(strsplit(str, sp))
	tmp[field]
}

getFieldFromExonNameNumeric <- function(str, field=1,  sp=":") {
	tmp<-unlist(strsplit(str, sp))
	as.numeric(tmp[field])
}

plotHistWithDensity <- function(X, ... ){
	hist(X, prob=TRUE, col="grey", ...)# prob=TRUE for probabilities not counts
	lines(density(X), col="blue", lwd=2) # add a density estimate with defaults
	lines(density(X, adjust=2), lty="dotted", col="darkgreen", lwd=2) 
	n<-length(X)
	m<-mean(X, na.rm=TRUE)
	stdev<-sd(X, na.rm=TRUE)
	sk<-skewness(X,type=2, na.rm=TRUE)
	kt<-kurtosis(X,type=2, na.rm=TRUE)
	legend("topright", top,legend=c("n:",n,"Mean:",m,"St. dev.:", stdev, "Skewness:",sk,"Kurtosis:",kt))
}

countDeletions <- function(x, sdLibs, minSigmaExon=3, minSigmaLib=3, maxValue=0.75){
	numT<-as.numeric(x)
	stdev<-sd(numT, na.rm=TRUE)
	minForDel <- 1 - ( minSigmaExon * stdev)
	minForDelLib <- 1 - (minSigmaLib * sdLibs )
	ret <- length(numT[numT < minForDel & numT < minForDelLib & numT < maxValue ] )
	ret
}

getDeletions <- function(x, sdLibs, minSigmaExon=3, minSigmaLib=3){
	numT<-as.numeric(x)
	stdev<-sd(numT, na.rm=TRUE)
	minForDel <- 1 - ( minSigmaExon * stdev)
	minForDelLib <- 1 - (minSigmaLib * sdLibs )
	ret <- numT[numT < minForDel & numT < minForDelLib] 
	ret
}

plotWithLibFilter <-function(x, sdLibs , minSigmaExon=3, minSigmaLib=3, ...){
	numT<-as.numeric(x)
	stdev<-sd(numT, na.rm=TRUE)
	minForDel <- 1 - ( minSigmaExon * stdev)
	minForDelLib <- 1 - (minSigmaLib * sdLibs )
	numT<-numT[minForDelLib > 0] 
	plotHistWithNormal(numT, ...)
}

plotHistWithNormal <- function(x,lengendPos="topright", ...){
	myhist <- hist(x,plot=F,...)
	multiplier <- myhist$counts / myhist$density
	mydensity <- density(x)
	mydensity$y <- mydensity$y * multiplier[1]

	plot(myhist,col="grey",...)
	lines(mydensity)

	myx <- seq(min(x), max(x), length.out= 100)
	mymean <- mean(x)
	mysd <- sd(x)

	normal <- dnorm(x = myx, mean = mymean, sd = mysd)
	lines(myx, normal * multiplier[1], col = "blue", lwd = 2)

	sd_x <- seq(mymean - 3 * mysd, mymean + 3 * mysd, by = mysd)
	sd_y <- dnorm(x = sd_x, mean = mymean, sd = mysd) * multiplier[1]

	segments(x0 = sd_x, y0= 0, x1 = sd_x, y1 = sd_y, col = "firebrick4", lwd = 2)

	n<-length(x)
	m<-mean(x, na.rm=TRUE)
	stdev<-sd(x, na.rm=TRUE)
	sk<-skewness(x,type=2, na.rm=TRUE)
	kt<-kurtosis(x,type=2, na.rm=TRUE)
	legend(lengendPos, top,legend=c("n:",n,"Mean:",m,"St. dev.:", stdev, "Skewness:",sk,"Kurtosis:",kt))
}



plotExonDetail <- function(x, sdLibs, title){
	numT<-as.numeric(x)
	stdev<-sd(numT, na.rm=TRUE)
	minForDel3s <- 1 - ( 3 * stdev)
	minForDel4s <- 1 - ( 4 * stdev)
	minForDelLib <- 1 - (3 * sdLibs )
	
	numT <- ifelse(numT > 1 , 1, numT)

	#Plot the minimum deletion for library
	minsLibTrimmed <- ifelse(numT < minForDelLib, numT, minForDelLib)
	minsLibTrimmed <- ifelse(minsLibTrimmed < 0, 0, minsLibTrimmed)
	vals <- numT - minsLibTrimmed

	#Plot the 4 sigma
	minForDel4sAdj <- minForDel4s - minsLibTrimmed
	minForDel4sTrimmed <- ifelse(vals < minForDel4sAdj, vals, minForDel4sAdj)
	minForDel4sTrimmed <- ifelse(minForDel4sTrimmed < 0, 0, minForDel4sTrimmed)
	vals <- vals - minForDel4sTrimmed
	vals <- ifelse(vals < 0, 0, vals)
	
	minForDel4sAdj <- ifelse(minForDel4sAdj < 0, 0, minForDel4sAdj)
	#Plot the 3 sigma
	minForDel3sAdj <- minForDel3s - minsLibTrimmed - minForDel4sAdj
	minForDel3sTrimmed <- ifelse(vals < minForDel3sAdj, vals, minForDel3sAdj)
	minForDel3sTrimmed <- ifelse(minForDel3sTrimmed < 0, 0, minForDel3sTrimmed)
	vals <- vals - minForDel3sTrimmed
	vals <- ifelse(vals < 0, 0, vals)
	
	value<- c(minsLibTrimmed, minForDel4sTrimmed, minForDel3sTrimmed, vals )
	lib<-c(names(x), names(x), names(x), names(x))
	value 
	type<-c(rep("1. 3 Sigma Library", length(x)), rep("2. 4 sigma exon", length(x)), rep("3. 3 sigma exon", length(x)),rep("4. Within 3 sigma Library", length(x)))
	tmpdf<-data.frame(lib,value,type)
	ggplot(tmpdf, aes(x=lib, y=value)) + 
	geom_bar(stat="identity", aes(fill = type)) +
	ggtitle(title)
}



plotDetailsOfContig <- function(contig){

	tempExons<-rownames(subset(df,Scaffold==contig))
	pdf(file=paste("Dist_", contig, ".pdf", sep =""), onefile=TRUE)
	for(exon in tempExons){
	   tmp<-mat[exon,]
	   plotHistWithNormal(tmp, breaks=100,main=exon, xlim=c(0,2.5), ylim=c(0,40))
	}
	dev.off()

	pdf(file=paste("NormalizedCoverage_", contig, ".pdf", sep =""), onefile=TRUE)
	for(exon in tempExons){
   		tmp<-mat[exon,]
   		print(plotExonDetail(tmp,samplesSD, paste(exon, "Normalized coverage")))
	}
	dev.off()
}

getSampleName<-function(library){
	splitted<-strsplit(as.character(library), "_")
	second<-splitted[[1]][2]
	first<-splitted[[1]][1]
	libShortName<-ifelse(substr(second,1,3)=="Cad", second,first)
	libShortName
}

getValueType <- function(value, sampleSD, exonSD){
	minForDel3s <- 1 - ( 3 * exonSD)
	minForDel4s <- 1 - ( 4 * exonSD)
	minForDelLib <- 1 - (3 * sampleSD)
	
	ret<-"0. Not set up"
	if(value < minForDelLib){
		if(value > minForDel3s){
			ret <- "2. Covered (3 sigma exon)"
		}
		if(value <= minForDel3s){
			ret <- "3. Deletion 3 sigma exon"
		}
		if(value <= minForDel4s){
			ret <- "4. Deletion 4 sigma exon"
		}
	}else{
		ret<-"1. Covered (3 sigma library)"
	}
	ret
}

longestAdjacentDeletionInLibrary <- function(library, tempDf, libValues, minForDelLib, maxForDel=0){
	tempExons <- rownames(tempDf)
	longestStretch<-0
	currentStretch<-0
	tempDf$minForDel <- 1 - ( 3 * tempDf$StdDev)
	prevDeletion<-F
	for (i in tempExons) {
		val<-libValues[i]
		minForDel <- tempDf[i,"minForDel"]
		if( minForDel < 0 ){
			next
		}
		if(val < minForDel && val < minForDelLib && val < maxForDel ){
			currentStretch <- currentStretch + 1
		}else{
			currentStretch <- 0
		}
		if(currentStretch > longestStretch){
			longestStretch<-currentStretch 
		}
			
	}
	longestStretch
}


longestWithGoodFlanking <- function(library, tempDf, libValues, minForDelLib, maxForDel=0.25, minForCovered=0.50){
	tempExons <- rownames(tempDf)
	minLibValue<- min(libValues)
	minForDel <- 1 - ( 3 * tempDf$StdDev)

	if(min(minLibValue) > maxForDel){
		return(0)
	}
	if(minForDel < 0){
		return(0)
	}
	longestStretch<-0
	currentStretch<-0
	
	prevDeletion<-F
	prevCovered<-T
	prevVal1<-1
	prevVal2<-1
	nextVal1<-1
	nextVal2<-1
	
	i<-1
	deletionsInScaffold <- data.frame(scaffold= character(0), start= integer(0), length = integer(0))
	while(i < length(tempExons)  ) {
	 	exonName<-tempExons[i]
	  	currentDeletion <- F

	  	if(i - 1 > 0){
	  		prevVal1 = libValues[i - 1]
	  	}
	  	if(i - 2 > 0){
	  		prevVal2 = libValues[i - 2]
	  	}
	  	val <- libValues[i]

	  	if(prevVal2 > minForCovered && prevVal2 > minForCovered){
	  		prevCovered <- T
	  	}
	  	if(val < maxForDel && val < minForDel){
	  		currentDeletion <- T
	  	}

	  	if(currentDeletion && prevCovered){
	  		j<-i
	  		delLength<-1
	  		while(j < length(tempExons) - 2 && currentDeletion){
	  			nextVal1<-libValues[j+1]
	  			nextVal2<-libValues[j+2]

	  		}
	  		i<-j
	  	}

	  	print(exonName)
	  	print(libValues[i])  	
		i<-i+1		
	}
	longestStretch
}



longestClearDeletions<- function(contig, df, mat, samplesSD, maxForDel=0.25, minForCovered=0.50){
	#contig<-"IWGSC_CSS_2BL_scaff_7939376"
	#library<-"LIB10929_Cadenza0106"	
	tempDf <- subset(df,Scaffold==contig)
	if(max(tempDf$Deletions3Lib3Exon) == 0){
		return(0)
	}
	longestStretch <- 0 
	tempExons <- rownames(tempDf)
	tmpMat<-mat[tempExons,]

  
	for(library in colnames(tmpMat) ){	
		libValues<-tmpMat[,library]
		minForDelLib <- 1 - (3 * samplesSD[library])
		currentStretch <- longestAdjacentDeletionInLibrary(library, tempDf, libValues,minForDelLib)
		if(currentStretch > longestStretch){
			longestStretch <- currentStretch
		}
	}
	longestStretch
} 

longestAdjacentDeletions<- function(contig, df, mat, samplesSD){
	tempDf <- subset(df,Scaffold==contig)
	if(max(tempDf$Deletions3Lib3Exon) == 0){
		return(0)
	}
	longestStretch <- 0 
	tempExons <- rownames(tempDf)
	tmpMat<-mat[tempExons,]
	#library<-"LIB11002_Cadenza0432"
	for(library in colnames(tmpMat) ){	
		libValues<-tmpMat[,library]
		minForDelLib <- 1 - (3 * samplesSD[library])
		currentStretch <- longestAdjacentDeletionInLibrary(library, tempDf, libValues,minForDelLib)
		if(currentStretch > longestStretch){
			longestStretch <- currentStretch
		}
	}
	longestStretch
} 


getMinimumDeletionValueInScaffold<-function(contig, df, mat, samplesSD ){
	tempExons<-rownames(subset(df,Scaffold==contig))

	
	tmpMat<-mat[tempExons,]
	melted<-melt(tmpMat)
	tmpDf<-as.data.frame(melted)
	if(length(tempExons) < 2){
		tmpDf$Var1 <- contig
		tmpDf$Var2 <- rownames(tmpDf)
	}
	libDf <- as.data.frame( names(samplesSD))
	names(libDf)[1] <- "Library"
	libDf$sdLib <- samplesSD
	libDf$MinSDForDel 	
	libDf$sampleName<-sapply(libDf$Library,getSampleName)
	
	tmpDf2<-subset(df,Scaffold==contig)
	tmpDf2$exon <-rownames(tmpDf2)
	innerDf<-merge(tmpDf2, tmpDf, by.x="exon", by.y="Var1")
	innerDf<-merge(innerDf, libDf, by.x="Var2", by.y="Library")
	innerDf$ValueType <- mapply(getValueType,innerDf$value, innerDf$sdLib, innerDf$StdDev)
	deletions<-subset(innerDf, grepl ("Del", ValueType))
	ret<-ifelse(nrow(deletions) > 0, deletions$value, 1)
	ret
}


#IWGSC_CSS_1AL_scaff_875531
plotScaffoldDeletions<- function(df, mat,  contig, samplesSD ){
	tempExons<-rownames(subset(df,Scaffold==contig))
	tmpMat<-mat[tempExons,]
	melted<-melt(tmpMat)
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
	innerDf$ValueType <- mapply(getValueType,innerDf$value, innerDf$sdLib, innerDf$StdDev)

	plot_Data <- ddply(innerDf, .(Start), mutate, Q1=quantile(value, 1/4), Q3=quantile(value, 3/4), IQR=Q3-Q1, upper.limit=Q3+1.5*IQR, lower.limit=Q1-1.5*IQR)
	gg<- ggplot(plot_Data, aes(x=factor(Start), y=value))
	gg<- gg + geom_errorbar(aes(ymax=1, ymin=1-3* StdDev), colour='gray', alpha=0.75)
	gg<- gg + geom_boxplot(aes (fill= factor(Deletions3Lib3Exon)) )
	gg<- gg + geom_point(data=plot_Data[plot_Data$value > plot_Data$upper.limit | plot_Data$value < plot_Data$lower.limit,], aes(x=factor(Start), y=value, col=factor(ValueType)))
	gg<- gg + coord_cartesian(ylim=c(0, 1))
	gg<- gg + geom_text(aes(label=ifelse(ValueType=="4. Deletion 4 sigma exon" , sampleName, ifelse(ValueType=="3. Deletion 3 sigma exon", sampleName, '')))  ,hjust=0,just=0, size=3, angle = 75)
	gg<- gg + ggtitle(paste("Normalized coverage for contig \n ", contig))
	gg
}


plotScaffoldDeletionsInLibrary<- function(df, mat,  contig, samplesSD ){
	tempExons<-rownames(subset(df,Scaffold==contig))
    tmpMat<-mat[tempExons,]
	melted<-melt(tmpMat)
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
	innerDf$ValueType <- mapply(getValueType,innerDf$value, innerDf$sdLib, innerDf$StdDev)


	pre_plot_Data <- ddply(innerDf, .(Start), mutate, Q1=quantile(value, 1/4), Q3=quantile(value, 3/4), IQR=Q3-Q1, upper.limit=Q3+1.5*IQR, lower.limit=Q1-1.5*IQR)

	libraries<-subset(pre_plot_Data, grepl("Del", ValueType))

	libNames<-unique(libraries$sampleName)
	plot_data<-subset(pre_plot_Data, sampleName %in% libNames )
	plot_data$Library <- plot_data$Var2
	gg<- ggplot(plot_data, aes(x=Start, y=value))
	gg<- gg + geom_errorbar(aes(ymax=1, ymin=1-3* StdDev), colour='gray', alpha=0.75) 
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








