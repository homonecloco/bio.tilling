#!/usr/bin/env Rscript
library("optparse")
library("sqldf")
library('hash')
library(GenomicRanges)
options(echo=TRUE) 
#devtools::install_github("homonecloco/bio.tilling")
#dirname(sys.frame(0)$ofile)
#programDir <- dirname(sys.frame(1)$ofile)

#devtools::install_local("/Users/ramirezr/Documents/public_code/bio.tilling")
require('bio.tilling')

option_list = list(
  make_option(c("-l", "--lines"),
   type="character", 
   default=NULL, 
   help="Coverage file", metavar="FILE"),
  make_option(c("-w", "--window_size"), 
    type="character", 
    default="200k,100k,075k,050k,025k,010k,005k",
    help='List of window sizes inside the main DIR ","',
    metavar="CHARACTER"),
  make_option(c("-m", "--max_gap"), 
    type="integer",
    default=3,
    help="Maximum gap length"
    ),
  make_option(c("-c", "--chromosomes"),
    type="character",
    default=NULL,
    help="File with the list of chromosomes",
    metavar="FILE"
    ),
  make_option(c("-i", "--index"),
    type="integer",
    default=NULL,
    help="The index of the line to calculate (optional)."
    ),
  make_option(c("-o", "--out"), 
    type="character", 
    default="./deletions_out",
    help="output folder [default= %default]. This folder most contain all the follders in 'window_size', each containing the files 'dels.csv.gz' and 'df.csv.gz'", 
    metavar="DIR")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt

if (is.null(opt$lines)){
  print_help(opt_parser)
  stop("Please select a line to extract.", call.=FALSE)
} 
windows<-strsplit(opt$window_size, ",")[[1]]
dels_h <- hash()
df_h   <- hash()
out    <- opt$out

lines = read.csv(opt$lines)[,1]
chromosomes = read.csv(opt$chromosomes)[,1]

if(!is.null(opt$index)){
  lines<-c(as.character(lines[opt$index]))
}

merged_path = paste0(out, "/line_deletions_merged")
dels_path = paste0(out, "/line_deletions")
dir.create(merged_path, recursive = TRUE)
dir.create(dels_path, recursive = TRUE)

for(line in lines){
  all_chr_dels <- NULL
  all_raw_dels <- NULL
  for (chr in chromosomes) { 
    print(paste0(line, ":",chr))
    chr_dels_h <- get_all_deletions(folder=out,
      deletions=windows,
      chr=chr,
      line=line, 
      max_gap=opt$max_gap, 
      dels_h=dels_h,
      df_h=df_h)
    chr_dels<-chr_dels_h[["merged_deletions"]]
    raw_dels<-chr_dels_h[["raw_deletions"]]      
    if(length(raw_dels) == 0 || is.na(raw_dels)){
      next
    }
    chr_dels <- data.frame(chr_dels)
    chr_dels <- chr_dels[,c("seqnames", "start", "end", "width", "strand","library", "window", "validated")]
    all_raw_dels <- all_raw_dels[,c("chr", "start", "ends", "library", "length",  "index_start", "index_end", "gap_exons","window")]
    
    if(is.null(all_raw_dels)){
      all_raw_dels<-raw_dels
    }else{
      all_raw_dels<-rbind(all_raw_dels,raw_dels)
    }
    
    if(is.null(all_chr_dels)){
      all_chr_dels<-chr_dels
    }else{
      all_chr_dels<-rbind(all_chr_dels,chr_dels)        
    }
    gc()
  }
  
 
  write.csv(all_chr_dels, paste0(merged_path, "/", line, ".csv"), row.names=F )
  write.csv(all_raw_dels, paste0(dels_path,   "/", line, ".csv"), row.names=F )
}
