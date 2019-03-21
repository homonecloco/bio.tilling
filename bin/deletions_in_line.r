#!/usr/bin/env Rscript
library("optparse")
library("sqldf")
options(echo=TRUE) 
#devtools::install_github("homonecloco/bio.tilling")
#dirname(sys.frame(0)$ofile)
#programDir <- dirname(sys.frame(1)$ofile)

#devtools::install_local("/Users/ramirezr/Documents/public_code/bio.tilling")
require('bio.tilling')

option_list = list(
  make_option(c("-l", "--line"),
   type="character", 
   default=NULL, 
   help="Coverage file", metavar="CHARACTER"),
  make_option(c("-m", "--max_gap"), 
    type="integer",
    default=3,
    help="Maximum gap length"
    ),
  make_option(c("-o", "--out"), type="character", default="./deletions_out", 
              help="output folder [default= %default]. This folder most contains the files 'dels.csv' and 'df.csv'", 
              metavar="DIR")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt

if (is.null(opt$line)){
  print_help(opt_parser)
  stop("Please select a line to extract.", call.=FALSE)
}



output_folder<-opt$out
setwd(output_folder)

df   <- read.csv("df.csv")
dels <- read.csv("dels.csv")

all_chr_dels <- NULL
for (chr in unique(df$Scaffold)) {
  chr_dels <- getDeletionsInChromosome(df, dels,
                                   chr=chr, 
                                   line=opt$line,
                                   max_gap = opt$max_gap
                                  ) 
  chr_dels$chr <- chr
  print(chr)
  if(ncol(chr_dels) > 6 && nrow(chr_dels) > 0){
    if(is.null(all_chr_dels)){
      all_chr_dels<-chr_dels
    }else{
      all_chr_dels<-rbind(all_chr_dels,chr_dels)
    }
  }
}


dir.create("line_deletions", recursive = TRUE)
setwd("line_deletions")

write.csv(all_chr_dels, paste0(opt$line, ".csv") )

