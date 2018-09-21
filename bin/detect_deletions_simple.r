#!/usr/bin/env Rscript
library("optparse")

#options(echo=TRUE) 
#devtools::install_github("homonecloco/bio.tilling")
#dirname(sys.frame(0)$ofile)
#programDir <- dirname(sys.frame(1)$ofile)

#devtools::install_local("/Users/ramirezr/Documents/public_code/bio.tilling")
require('bio.tilling')

option_list = list(
  make_option(c("-f", "--coverage_file"), type="character", default=NULL, 
              help="Coverage file", metavar="character"),
  make_option(c("-z", "--gzip"), action="store_true", default=FALSE,
            help="The coverage file is gzipped [default]"),
  make_option(c("-o", "--out"), type="character", default="./deletions_out", 
              help="output folder [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt

if (is.null(opt$coverage_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (coverage_file).", call.=FALSE)
}



filename<-unlist(opt$coverage_file)
output_folder<-opt$out
is_gz <- opt$gzip

dir.create(output_folder, recursive = TRUE)
#cores<-as.integer(args[2])

covs<-readCoverageTable(filename, is_gz=is_gz)
covs<-filterSamplesPerSD(covs, maxSD=0.3)

df<-getExonsDF(covs)
mat<-normalizeCovs(covs, df)
rm(covs)

mat<-filterLowQualityExons(mat, maxSD=0.3)
df<-getExonsDF(mat)

dels<-getAllDeletedExons(mat,df)
libSD<-getLibSD(mat)

#scaffsWithDels<-getAllScaffoldsWithDeletions(df,dels)

#deslWithAVGs<-getAllScaffoldAveragesParallel(scaffsWithDels, mat,df, cores=cores)

#selectedDels<-deslWithAVGs[deslWithAVGs$Score>=1,]

write.csv(dels, file="dels.csv")
write.csv(deslWithAVGs, file='deslWithAVGs.csv')
write.csv(df, file='df.csv')
write.csv(libSD, file='libSD.csv')
write.csv(mat, file='mat.csv')

