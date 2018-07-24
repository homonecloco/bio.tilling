options(echo=TRUE) 
#devtools::install_github("homonecloco/bio.tilling")
require('bio.tilling')
args <- commandArgs(trailingOnly = TRUE)

filename<-args[1]
cores<-as.integer(args[2])

covs<-readCoverageTable(filename)
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

