
getFieldFromExonName <- function(str, field=1,  sp=":") {
	tmp<-unlist(strsplit(str, sp))
	tmp[field]
}

getFieldFromExonNameNumeric <- function(str, field=1,  sp=":") {
	tmp<-unlist(strsplit(str, sp))
	as.numeric(tmp[field])
}

plotHistWithDensity <- function(X, ... ){
	library(e1071)
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
	library(e1071)
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

getDeletionsInChromosome<-function(exons_df, dels,
                                   chr="IWGSC_3BSEQ_3B_traes3bPseudomoleculeV1", 
                                   line="Cadenza0154_1158_LIB10947_LDI9034",
                                   max_gap = 2,
                                   window = ""
                                  ){
    exons<-exons_df[exons_df$Scaffold == chr,  c("Exon", "Start", "Ends")]
    dels<-dels[dels$Library== line & dels$HomDel == TRUE & dels$Scaffold == chr , c("Exon", "Library", "HomDel" ),]
    dels_df <- sqldf("SELECT DISTINCT exons.*, dels.Library, dels.HomDel FROM exons LEFT JOIN dels ON exons.Exon == dels.Exon ORDER BY Start")
    dels_df$Library <- line
    dels_df$HomDel <- ifelse(is.na(dels_df$HomDel), FALSE, dels_df$HomDel)
    tempExons <- rownames(dels_df)
    current_strech <- list(chr=chr, start=0, ends=0, library=line, length=0, index_start=0, index_end=0, gap_exons=0)
    found_deletions<-data.frame(row.names = c("chr","start","ends","library","length","index_start", "index_end", "gap_exons"))
    current_gap <- 0
    for (i in tempExons) {
        val<-dels_df[i,]
        if( val$HomDel == FALSE ){
            if(current_strech$length > 0){
                current_gap <- current_gap + 1
                if(current_gap > max_gap){
                    current_strech <- data.frame(current_strech)
                    found_deletions<-rbind(found_deletions, current_strech)   
                    current_strech <- list(chr=chr, start=0, ends=0, library=line, length=0, index_start=0, index_end=0, gap_exons=0)
                    current_gap <- 0
                }
            }     
            next
        }        
        if(current_strech$length == 0){
            current_strech$start <- val$Start
            current_strech$index_start <- i
        }
        current_strech$gap_exons <- current_strech$gap_exons + current_gap
        current_strech$ends = val$Ends
        current_strech$length <- current_strech$length + 1
        current_strech$index_end <- i
        current_gap <- 0 
    }
    if(current_strech$length > 0){
        current_strech <- data.frame(current_strech)
        found_deletions<-rbind(found_deletions, current_strech)  
    }
    if(length(found_deletions) > 0){
        found_deletions$window <- window
    }
    found_deletions
}


validate_deletions <- function(top, bottom, min_overlap=0.5){
    hits     <- findOverlaps(bottom,top)
    overlaps <- pintersect(bottom[queryHits(hits)], top[subjectHits(hits)])
    hits_df  <- data.frame(hits)
    hits_df$topPercentOverlap <- width(overlaps) / width( top[subjectHits(hits)])
    top_overlap <- sqldf("SELECT subjectHits as top, sum(topPercentOverlap) as overlap FROM hits_df GROUP BY subjectHits")
    top$validated <- F
    top[top_overlap$top ]$validated <- ifelse(top_overlap$overlap > min_overlap,
                                              TRUE, 
                                              top[top_overlap$top ]$validated )
    top
}

merge_deletions<-function(top_in, bottom_in, min_overlap=0.5, min_validate_cov=1){
    bottom_window<-unique(bottom_in$window)
    top   <-reduce(top_in)
    top$library   <- top_in$library
    top$window    <- top_in$window
    top$validated <- top_in$validated

    pre_validated<-top[top$validated==TRUE, ]
    
    bottom<-reduce(bottom_in)
    bottom$library   <- bottom_in$library
    bottom$window    <- bottom_in$window
    bottom$validated <- bottom_in$validated

    merged<-unlist(as(list(top,bottom), "GRangesList"))
    merged<-reduce(merged)

    #We are getting which windows are validated on either side. Also, we merge to the validated set the deletions previously validated
    top_validated   <-validate_deletions(top, bottom, min_overlap=min_overlap)
    bottom_validated<-validate_deletions(bottom, top, min_overlap=min_overlap)
    validated<-unlist(as(list(top_validated,bottom_validated,pre_validated), "GRangesList"))
    validated<-validated[validated$validated==TRUE, ]    
    
    #We validate now from the large merge to the validated windows, but now the overlap needs to be more than 1
    merged<-validate_deletions(merged,validated, min_overlap=min_validate_cov)
    merged$library <- unique(top$library)
    
    #This section gets the topmost window where there is at least some evidence of the deletion. 
    hits     <- findOverlaps(top,merged)
    hits_df  <- data.frame(hits)
    colnames(hits_df)<-c("top","merged")
    hits_df$window   <-  top[hits_df$top,]$window
    merged$window <- bottom_window
    merged[hits_df$merged,]$window <-  top[hits_df$top,]$window    
    merged
}

#The folder should contain different subfolders with the different window sizes. 
#The deletions are run in the order of the array. 
get_all_deletions<-function(folder=".",
                            deletions=c("200k","100k","075k","050k","025k","010k"),
                           chr="chr5D",
                           line="J1.33_k80d50Mc5", 
                           max_gap=5, 
                           df_h  = hash(),
                           dels_h = hash()){
    top <- NA
    first<-TRUE
    raw_dels<-NA
    for(i in deletions){
        if(!all(has.key(i, dels_h))){
            path<-paste0(folder,"/",i,"/df.csv.gz")
            df_h[i] <- read.csv(gzfile(path))
            path<-paste0(folder,"/",i,"/dels.csv.gz")
            tmp<-read.csv(gzfile(path))
            dels_h[i] <- tmp[tmp$HomDel, ]
        }
        df   <- df_h[[i]]
        dels <- dels_h[[i]]
        chr_dels <- getDeletionsInChromosome(df, dels,
                                             chr=chr,
                                             line=line,
                                             max_gap = max_gap,
                                             window=i)
        
        if(nrow(chr_dels) == 0 || ncol(chr_dels) == 0){
            next
        }
        range_chr_dels <- makeGRangesFromDataFrame(chr_dels, end.field="ends", keep.extra.columns = T)
        range_chr_dels$validated <- F
        if(first){
            range_chr_dels$length      <- NULL
            range_chr_dels$index_start <- NULL
            range_chr_dels$index_end   <- NULL
            range_chr_dels$gap_exons   <- NULL
            top      <-range_chr_dels
            raw_dels <- chr_dels
            first    <- FALSE
            next
        }
        top<-merge_deletions(top, range_chr_dels)
        raw_dels <- rbind(raw_dels,chr_dels)
    } 
    #top
    hash("raw_deletions" = raw_dels, "merged_deletions" = top)
}







