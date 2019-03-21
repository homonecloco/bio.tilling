#!/usr/bin/env Rscript
library("optparse")
options(echo=TRUE) 
require('bio.tilling')

option_list = list(
  make_option(c("-f", "--coverage_file"), 
  	type="character", 
  	default=NULL, 
  	help="Coverage file", 
  	metavar="character"),
  make_option(c("-s", "--sampleSD"), 
  	type="double", 
  	default=0.3, 
  	help="Maximum Standard Deviation allowed in sample [default= %default]"),
  make_option(c("-w", "--windowSD"), 
  	type="double", 
  	default=0.3,
  	help="Maximum Standard Deviation allowed in each window [default= %default]"),
  make_option(c("-z", "--gzip"), 
  	action="store_true", 
  	default=FALSE,
  	help="The coverage file is gzipped [default]"),
  make_option(c("-o", "--out"), 
  	type="character", 
  	default="./deletions_out", 
  	help="output folder [default= %default]", 
  	metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt

if (is.null(opt$coverage_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (coverage_file).",
  	call.=FALSE)
}

filename<-unlist(opt$coverage_file)
output_folder<-opt$out
is_gz <- opt$gzip

covs<-readCoverageTable(filename, is_gz=is_gz)
covs<-filterSamplesPerSD(covs, maxSD=opt$sampleSD)

df<-getExonsDF(covs)
mat<-normalizeCovs(covs, df)
rm(covs)

dir.create(output_folder, recursive = TRUE)
setwd(output_folder)

mat<-filterLowQualityExons(mat, maxSD=opt$windowSD)
gc()
write.csv(mat, file='mat.csv')

df<-getExonsDF(mat)
write.csv(df, file='df.csv')

dels<-getAllDeletedExons(mat,df)
write.csv(dels, file="dels.csv")

libSD<-getLibSD(mat)
write.csv(libSD, file='libSD.csv')