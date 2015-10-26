

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
