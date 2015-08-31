#' Coutns how many missing values are in the the list. 
#'
#' @param x: The list of values. It is transformed to a numneric vector. 
#'
#' @return The count of missing values (x == 0.0)

countMissing <-function(x){
	numT<-as.numeric(x)
	ret <- length(numT[numT==0.0])
	ret
}


#' Remove the libraries with bad standard deviation
#'
#' 
#' @param mat:  The matrix with the normalized coverage. The row names have the format "Scaffold:start:end". 
#'				Each column represents a library 
#' @param maxSD: The maximum standard deviation in the library to be considered as a good deletion. (default: 0.3)
#'
#'

filterSamplesPerSD<-function(covs, maxSD=0.3) {
	df<-getExonsDF(covs)
	mat<-as.matrix(normalizeCovs(covs, df))
	missing<-apply(mat,1,countMissing)
	libSD<-apply(mat[missing<1,],2,sd)
	fineSD<-libSD < maxSD
	covs[,fineSD]
}


#' Filters all the exons that have low quality. 
filterLowQualityExons <-function(mat, maxSD=0.3){
	exonSD<-apply(mat, 1, sd)
 	fineSD<-exonSD < maxSD
 	mat[fineSD,]
}



#' Reads the table with the coverages. The table has headers containing the libary name. 
#' the rows must be named with the format "Scaffold:start:end". 
#' 
#' @param filename: The filename with the table. 
readCoverageTable<-function(filename){
	counts<-read.table(filename, header=TRUE, sep="\t")
	counts$remove <-NULL
	counts<-counts[,colSums(counts,na.rm=T)>0]
	counts
}



#' Prepares the dataframe with the information per exon. This parses the rownames to get the scaffold
#' Start and end. It also calculates the exon length
#' @param mat: The matrix with the data
#' @return A data frame with the following rows: Exon, Scaffold, Start, End, ExonL (exon length)
getExonsDF<-function(mat){
	names<-rownames(mat)
	exons<-strsplit(as.character(names), ":")
	scaffolds<-sapply(exons, "[[", 1)
	starts<-as.integer(sapply(exons, "[[", 2))
    ends<-as.integer(sapply(exons, "[[", 3))
	exonLengths<-ends-starts
	sdExons<-apply(mat,1,sd)
    df<-data.frame(Exon=names, Scaffold=scaffolds, Start=starts, Ends=ends, ExonL=exonLengths, sdExon=sdExons)
    df
}

getExonsDetails<-function(names){
	exons<-strsplit(as.character(names), ":")
	scaffolds<-sapply(exons, "[[", 1)
	starts<-as.integer(sapply(exons, "[[", 2))
    df<-data.frame(Exon=names, Scaffold=scaffolds, Start=starts)
    df
}



#' Function that normalizes the coverages, in an RPKM-like values, and normalized to have a mean of 1, at the exon level.
normalizeCovs<-function(counts, exonsDF) {
	exonLengths<-exonsDF$ExonL
	totalReadsPerSample<-apply(counts,2,sum)
	multiplier<-outer(exonLengths, totalReadsPerSample )
	multiplier<-1/multiplier
	multiplier<-1000000000*multiplier
	#"RPKM-like covs"
	multiplier<-counts*multiplier
	rowsToRemove<-rowSums(multiplier)>0
	
	#Filtering the rows with AVG=0 
	multiplier<-multiplier[rowsToRemove,]
	
	#Normalize by multiplying by the mean of each row. 	
	mean.mult <- apply(multiplier, 1, mean)
	covs<-sweep(multiplier,MARGIN=1,mean.mult,'/')
	covs
}


homDeletionFilter<-function(localMat, exonsDF, maxValueForDeletion=0.1, minSigmaExon=3, minSigmaLib=3 ) {
	sdExons<-exonsDF$sdExon
	missing<-apply(localMat,1,countMissing)
	sdLibs<-apply(localMat[missing<1,],2,sd)

	sdExons <- 1 - ( minSigmaExon * sdExons)#This is a vector
	sdLibs <- 1 - (minSigmaLib * sdLibs )

	sdLibs<- ifelse(sdLibs  < maxValueForDeletion, sdLibs , maxValueForDeletion)
	sdExons<-ifelse(sdExons < maxValueForDeletion, sdExons, maxValueForDeletion)
	
	sdLibsMat  <- matrix(sdLibs ,nrow=length(sdExons),ncol=length(sdLibs), byrow=FALSE)
	sdExonsMat <- matrix(sdExons,nrow=length(sdExons),ncol=length(sdLibs), byrow=TRUE)
	ret <- pmin(sdExonsMat,sdLibsMat)
	ret
}

getHomExoDeletions<-function(localMat, exonsDF, maxValueForDeletion=0.1, ...){
	library(reshape2)
	library(plyr)
	filterMat<-homDeletionFilter(localMat, exonsDF, maxValueForDeletion=maxValueForDeletion,...)
	dels <- localMat < filterMat
	meltedDels<-melt(dels)
	meltedDels<-rename(meltedDels, c("Var1"="Exon", "Var2"="Library", "value"="HomDel"))
	meltedDels
}

getHetExoDeletions<-function(localMat, exonsDF, minSigmaExonHet=1) {
	library(reshape2)
	library(plyr)
	sdExons<-exonsDF$sdExon
	sdExons <- 0.5 + (minSigmaExonHet * sdExons)#This is a vector
	sdExonsMat <- matrix(sdExons,nrow=nrow(localMat),ncol=ncol(localMat), byrow=TRUE)
	dels <- localMat < sdExonsMat
	meltedDels<-melt(dels)
	#meltedDels<- meltedDels[meltedDels$value,]
	meltedDels<-rename(meltedDels, c("Var1"="Exon", "Var2"="Library", "value"=paste("HetDel",minSigmaExonHet,sep="")))
	meltedDels
}

getAllDeletedExons<-function(localMat, exonsDF){
	library(reshape2)
	library(plyr)
	meltedMat<-melt(as.matrix(localMat),varnames =c("Exon", "Library"), value.name="NormCov")
	#meltedMat<-rename(meltedMat, c("Var1"="Exon", "Var2"="Library", "value"="NormCov"))
	
	dat<-getHomExoDeletions(localMat, exonsDF)
	meltedMat$HomDel <- dat$HomDel
	
	dat<-getHetExoDeletions(localMat, exonsDF, 1)
	meltedMat$HetDel1 <- dat$HetDel1
	meltedMat$Scaffold<-exonsDF$Scaffold
	meltedMat$Start<-exonsDF$Start
	meltedMat[meltedMat$HetDel1 | meltedMat$HomDel,]
}



getAllScaffoldsWithDeletions<-function(exonsDF, delsDF,  minValidExonLength=9){
	library(sqldf)
	allScaffcounts <- table(exonsDF$Scaffold)
	scaffoldsDF <- as.data.frame(rownames(allScaffcounts))
	names(scaffoldsDF)[1]<-"Scaffold"
	scaffoldsDF$validExons <- allScaffcounts
	tempMat<-delsDF
	tempMat<-delsDF[delsDF$HetDel1 | delsDF$HomDel, ]
	tempMat<-delsDF
	countsHom<-sqldf("SELECT Scaffold, Library, COUNT(*) as HomDels FROM tempMat WHERE HomDel GROUP BY Scaffold, Library ")
	countsHet<-sqldf("SELECT Scaffold, Library, COUNT(*) as HetDels FROM tempMat WHERE HetDel1 GROUP BY Scaffold, Library ")
	mergedScaffDF <- merge(scaffoldsDF, countsHet, all.x=T)
	mergedScaffDF <- merge(mergedScaffDF, countsHom, all.x=T)
	mergedScaffDF$HomDels[is.na(mergedScaffDF$HomDels)] <- 0
	mergedScaffDF$HetDels[is.na(mergedScaffDF$HetDels)] <- 0
	mergedScaffDF<-na.omit(mergedScaffDF)
	mergedScaffDF<-mergedScaffDF[mergedScaffDF$validExons > minValidExonLength,]
	mergedScaffDF
}

getScaffoldAverages<-function(scaffold, library, localMat, exonsDF, minSigmaExonHet=1, maxValueForDeletion=0.1, minSigmaExon=3) {
	exonsToProcessDF<-exonsDF[exonsDF$Scaffold==scaffold,]
	exonsToProcess<-rownames(exonsToProcessDF)
	scaffMat<-localMat[exonsToProcess,library]
	allAvg<-mean(scaffMat)
	allNoHom<-mean(scaffMat[scaffMat > maxValueForDeletion])
	allNoHom<-ifelse(is.na(allNoHom), 1, allNoHom)
	maxForHetDel<-	0.5 + exonsToProcessDF$sdExon * minSigmaExonHet
	maxForDel <- 1 - 3*exonsToProcessDF$sdExon
	allNoDels <- mean(scaffMat[scaffMat > maxForDel])
	allNoHet<-mean(scaffMat[scaffMat > maxForHetDel])
	allNoHet<-ifelse(is.na(allNoHet), 1, allNoHet)
	allNoDels<-ifelse(is.na(allNoDels), 1, allNoDels)
	arg0<-list(AllAvg=allAvg, AllNoHomAvg=allNoHom, AllNoHetAvg=allNoHet,AllNoDelsAvg=allNoDels)
	df<-data.frame(t(unlist(arg0)))
	df
}


getAllScaffoldAverages<-function(delsMat, localMat, exonsDF){

	delsMat$AllAvg<-1.0
	delsMat$AllNoHomAvg<-1.0 
	delsMat$AllNoHetAvg <-1.0
	delsMat$AllNoDelsAvg<-1.0
	for(i in 1:nrow(delsMat)) {
	    row <- delsMat[i,]
	    print(row)
	    tmpDF <- getScaffoldAverages(row$Scaffold, row$Library, localMat, exonsDF)
	    print(tmpDF)
	    delsMat[i,"AllAvg"]      <-tmpDF$AllAvg 		
		delsMat[i,"AllNoHomAvg"] <-tmpDF$AllNoHomAvg  
		delsMat[i,"AllNoHetAvg"] <-tmpDF$AllNoHetAvg 
		delsMat[i,"AllNoDelsAvg"]<-tmpDF$AllNoDelsAvg
	}
}


getWholeExonDeletions<-function(scaffoldsDF, exonsDF, localMat, sdLibs, deletions, maxValueForDeletion=0.1, minSigmaExon=3, minSigmaLib=3){
	library(reshape2)
	library(plyr)

	libraries<-colnames(localMat)
	deletions<-sapply(libraries,getWholeExonDeletionsLibrary, scaffoldsDF, exonsDF, localMat, sdLibs, maxValueForDeletion=maxValueForDeletion, minSigmaExon=minSigmaExon, minSigmaLib=minSigmaLib)
	rownames(deletions)<-scaffoldsDF$Scaffold
	deletionsMelted<-melt(deletions)
	
	deletionsMelted<-deletionsMelted[deletionsMelted$value >0,]
	deletionsMelted$Category <- "Whole Exon"
	deletionsMelted<-rename(deletionsMelted, c("Var1"="Scaffold", "Var2"="Library"))
	delsJoined <- merge(x=deletionsMelted, y=scaffoldsDF, by = "Scaffold", all.x = TRUE )
	delsJoined
}




categorizeDeletions<-function(deletionsMelted, exonsDF){

	all_scaffold_counts <- sort(table(exonsDF$Scaffold))
	scaffoldsDF <- as.data.frame(rownames(all_scaffold_counts))
	names(scaffoldsDF)[1]<-"Scaffold"
	missing<-apply(localMat,1,countMissing)
	scaffoldsDF$validExons <- all_scaffold_counts
	scaffoldsDF<-subset(scaffoldsDF, validExons > 0)
	sdLibs<-apply(localMat[missing<1,],2,sd)
	deletions<-getWholeExonDeletions(scaffoldsDF, exonsDF, localMat, sdLibs, maxValueForDeletion=0.1, minSigmaExon=3, minSigmaLib=3)
	deletions
}

