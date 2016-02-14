
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

setwd(folder)

f<-file('dels.csv')
dels<-sqldf('select * from f', dbname = tempfile(), file.format = list(sep=',',header = T, row.names = F))

f<-file('mat.csv')
mat<-sqldf('select * from f', dbname = tempfile(), file.format = list(sep=',',header = T, row.names = F))


deslWithAVGs<-read.csv('deslWithAVGs.csv')
df<-read.csv('df.csv')
geneticMap<-read.csv('geneticMap.csv')
libSD<-read.csv('libSD.csv')
scaffsWithDels<-read.csv('scaffsWithDels.csv')
selectedDels<-read.csv('selectedDels.csv')

chromosomes<-unique(geneticMap$chr)

libraries<-data.frame(Library=colnames(mat))

for (i in 1:length(chromosomes)) {
    chr <- chromosomes[i]
    
    plotFile=paste0("chr",chr,"_min5_grouped.pdf")
    pdf(plotFile, height = 10, width = 8)
    plotPerChromosome(geneticMap,selectedDels, libraries,df,column='homDelsPer', groupByCM=T, chr=chr, minExonCount=4)  
    dev.off() 
    
    plotFile=paste0("chr",chr,"_min5.pdf")
    pdf( plotFile, height = 10, width = 8)
    plotPerChromosome(geneticMap,selectedDels, libraries,df,column='homDelsPer', groupByCM=F, chr=chr, minExonCount=4)  
    dev.off()
    
    
    
}

homSelDel<-getHomSelectedDeletions(selectedDels,geneticMap)
write.csv(homSelDel, file="SelectedHomDeletions.csv")

