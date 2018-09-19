#!/usr/bin/env Rscript

options(echo=TRUE) 
#devtools::install_github("homonecloco/bio.tilling")
#require('bio.tilling')
library(sqldf)



thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

print(thisFile())
path<- dirname(thisFile())
print(path)

deletions_file<-paste0(path,"/../R/deletionCategory.R")
plot_file<-paste0(path,"/../R/plotFunctions.R")

source(deletions_file)
source(plot_file)

args <- commandArgs(trailingOnly = TRUE)

folder<-args[1]
parser<-ifelse(args[2], args[2], "identity")
parser<-args[2]

setwd(folder)

f<-file('dels.csv')
dels<-sqldf('select * from f', dbname = tempfile(), file.format = list(sep=',',header = T, row.names = F))

f<-file('mat.csv')
mat<-sqldf('select * from f', dbname = tempfile(), file.format = list(sep=',',header = T, row.names = F))


deslWithAVGs<-read.csv('deslWithAVGs.csv')
df<-read.csv('df.csv')
geneticMap<-read.csv('geneticMap.csv')
libSD<-getLibSD(mat)
#libSD<-read.csv('libSD.csv')
scaffsWithDels<-read.csv('scaffsWithDels.csv')
selectedDels<-read.csv('selectedDels.csv')

chromosomes<-unique(geneticMap$chr)

libraries<-data.frame(Library=colnames(mat))

cadenzaParseName <- function(x) {
    vect<-strsplit(x,"_")
    val <- vect[[1]][2] 
    #print(val)
    #print(vect)
    return (val)

}

kronosParseNames<- function(x){
    t2<-strsplit(x,"[.]")[[1]][2]
    val<-strsplit(t2,"_")[[1]][1]
    return(val)
}

for (i in 1:length(chromosomes)) {
    chr <- chromosomes[i]
    
    plotFile=paste0("chr",chr,"_min5.pdf")
    pdf( plotFile, height = 10, width = 8)
    code<-paste0("plotPerChromosome(geneticMap,selectedDels, libraries,df,column='homDelsPer', 
        groupByCM=T,groupByCMPrecision=1, chr=chr, minExonCount=4, parseLibFun=",parser,")")  
    print(code)
    eval(parse(text=code))

    dev.off()

}

homSelDel<-getHomSelectedDeletions(selectedDels,geneticMap)
write.csv(homSelDel, file="SelectedHomDeletions.csv")

