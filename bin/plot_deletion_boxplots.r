library(bio.tilling)
library(sqldf)



options(echo=TRUE) 
#devtools::install_github("homonecloco/bio.tilling")
require('bio.tilling')
args <- commandArgs(trailingOnly = TRUE)

#/Volumes/ramirezrVMs/deletions/deletionTablesV1
#Contains dels.csv, mat.csv, df.csv, libSD.csv and geneticMap.csv
#All produced by detect_deletions.r
inputFolder<-args[1]

outputFolder<-args[2]
#CSV with two named columns: Scaffold and Library
fileWithDeletions<-args[3]


setwd(inputFolder)
f<-file('dels.csv')
dels<-sqldf('select * from f', dbname = tempfile(), file.format = list(sep=',',header = T, row.names = F))
f<-file('mat.csv')

mat<-sqldf('select * from f', dbname = tempfile(), file.format = list(sep=',',header = T, row.names = F))
mat2 <- mat[,-1]
newRowNames<-lapply( mat[,1], function(x) gsub("\"", "", x))
rownames(mat2) <-newRowNames
head(mat2)
mat<-mat2

deslWithAVGs<-read.csv('deslWithAVGs.csv')
df<-read.csv('df.csv')
rownames(df)<-df$X

df<-df[,!(names(df) %in% c("X"))]
geneticMap<-read.csv('geneticMap.csv')

libSD<-read.csv('libSD.csv')
names(libSD)<-libSD$X
vect<-libSD$x
names(vect)<-libSD$X
head(vect)
libSD<-vect

scaffsWithDels<-read.csv('scaffsWithDels.csv')
selectedDels<-read.csv('selectedDels.csv')
