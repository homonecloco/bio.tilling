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

getWholeExonDeletionsLibraryInner <- function(scaffold, library, exonsDF, mat,libSD, maxValueForDeletion=0.1, minSigmaExon=3, minSigmaLib=3 ){
	tempDf <- subset(exonsDF,Scaffold == scaffold)
	tempExons <- rownames(tempDf)
	tmpMat<-as.vector(t(mat[tempExons,library]))
	numT<-as.numeric(tmpMat)
	minForDel <- 1 - ( minSigmaExon * tempDf$sdExon)#This is a vector
	minForDelLib <- 1 - (minSigmaLib * libSD )
	ret <- numT[numT < minForDel & numT < minForDelLib & numT < maxValueForDeletion] 
	length(ret)
}

getWholeExonDeletionsLibrary <- function(library, scaffoldDF, exonsDF, mat, sdLibs, maxValueForDeletion=0.1, minSigmaExon=3, minSigmaLib=3){
	libSD<-sdLibs[library]
	delCounts<-sapply(scaffoldDF$Scaffold,getWholeExonDeletionsLibraryInner, library, exonsDF, mat, libSD,  maxValueForDeletion=maxValueForDeletion, minSigmaExon=minSigmaExon, minSigmaLib=minSigmaLib )
	delCounts	
}

getWholeExonDeletions<-function(scaffoldsDF, exonsDF, localMat, sdLibs, maxValueForDeletion=0.1, minSigmaExon=3, minSigmaLib=3){
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


getWholeScaffoldHetDeletionsLibraryInner <- function(scaffold, library, exonsDF, mat,libSD, maxValueForDeletion=0.1, minSigmaExon=3, minSigmaLib=3 ){
	tempDf <- subset(exonsDF,Scaffold == scaffold)
	tempExons <- rownames(tempDf)
	tmpMat<-as.vector(t(mat[tempExons,library]))
	numT<-as.numeric(tmpMat)
	minForDel <- 0.5 + ( minSigmaExon * tempDf$sdExon)#This is a vector
	minForDelLib <- 0.5 + (minSigmaLib * libSD )
	ret <- numT[numT < minForDel & numT < minForDelLib & numT < maxValueForDeletion] 
	length(ret)
}

getWholeScaffoldHetDeletionsLibrary <- function(library, scaffoldDF, exonsDF, mat, sdLibs, maxValueForDeletion=0.1, minSigmaExon=3, minSigmaLib=3){
	libSD<-sdLibs[library]
	delCounts<-sapply(scaffoldDF$Scaffold,getWholeScaffoldHetDeletionsLibraryInner, library, exonsDF, mat, libSD,  maxValueForDeletion=maxValueForDeletion, minSigmaExon=minSigmaExon, minSigmaLib=minSigmaLib )
	delCounts	
}

getWholeScaffoldHetDeletions<-function(scaffoldsDF, exonsDF, localMat, sdLibs, maxValueForDeletion=0.1, minSigmaExon=3, minSigmaLib=3){
	library(reshape2)
	library(plyr)

	libraries<-colnames(localMat)
	deletions<-sapply(libraries,getWholeScaffoldHetDeletionsLibrary, scaffoldsDF, exonsDF, localMat, sdLibs, maxValueForDeletion=maxValueForDeletion, minSigmaExon=minSigmaExon, minSigmaLib=minSigmaLib)
	rownames(deletions)<-scaffoldsDF$Scaffold
	deletionsMelted<-melt(deletions)
	
	deletionsMelted<-deletionsMelted[deletionsMelted$value >0,]
	deletionsMelted$Category <- "Whole Het scaffold"
	deletionsMelted<-rename(deletionsMelted, c("Var1"="Scaffold", "Var2"="Library"))
	delsJoined <- merge(x=deletionsMelted, y=scaffoldsDF, by = "Scaffold", all.x = TRUE )
	delsJoined
}


categorizeDeletions<-function(localMat, exonsDF){
	library(parallel)
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

