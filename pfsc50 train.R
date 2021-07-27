#vcfR training with pcf files
setwd("~/R training/vcfR")
library(vcfR)
library(pinfsc50)
library(tidyverse)
pkg <- "pinfsc50"

#Location for data files from pinfsc50
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = pkg)
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = pkg)
gff_file <- system.file("extdata", "pinf_sc50.gff", package = pkg)

#read files into r
vcf <- read.vcfR( vcf_file, verbose = FALSE )
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="")
str(dna)

#Create chromR object
chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)
plot(chrom)

#Masking variants with qual<1, DP <300/>700, MQ btw 59.9 and 60.1
chrom <- masker(chrom, min_QUAL = 1, min_DP = 300, max_DP = 700, 
                min_MQ = 59.9, max_MQ = 60.1)

#See table of variants
variant.table(chrom)

#Plot chrom object after masking
plot(chrom)
head(vcf@gt)

#Process chrom object after masking. 
chrom <- proc.chromR(chrom, verbose=TRUE)
plot(chrom)

#Quality control
chromoqc(chrom, dp.alpha=255)

#Zoom on 5e5 to 6e5
chromoqc(chrom, xlim=c(5e+05, 6e+05))

##NYT EKSEMPEL. Kigger på Vcf file format. 3 regioner, meta, fix og genotyp
## META
data(vcfR_example)
strwrap(vcf@meta[1:8])
queryMETA(vcf)

##Fokus on one element
queryMETA(vcf, element = 'DP')

##Fokus on one element of the element
queryMETA(vcf, element = 'FORMAT=<ID=DP')

## FIX, there are more commands to get from the fix position look under ?getFIX
head(getFIX(vcf))

## GT region
vcf@gt[1:6, 1:4]

#MANIPULATION OF vcfR files
vcf <- read.vcfR( vcf_file, verbose = FALSE )
head(is.polymorphic(vcf, na.omit = TRUE))
head(is.biallelic(vcf))
vcf2 <- extract.indels(vcf) ## Stores all non-indels (SNV) in another vcf object

#Subsetting
vcf[,1:19]

#Memory load calculator... Can be used for real data
library('memuse')
nvar <- 10^{2:8}
nMb <- howbig(nrow = 10, ncol = 10, unit = "MB")@size

for(i in nvar){
  nMb <- c(nMb, howbig(nrow = i, ncol = 10, unit = "MB")@size)
}
nvar <- (c(10, nvar))
par(mar=c(5,5,4,2))
plot(log10(nvar), log10(nMb), xaxt="n", yaxt="n", type='b', xlab="Number of variants", ylab = "")
axis(side = 1, at = log10(nvar), labels=nvar)
axis(side = 2, at = log10(nMb), labels=nMb, las=2)
title(ylab="Memory use (Mb)", line=4)
abline(h=log10(nMb), lwd=2, col="#C0C0C066")
abline(v=log10(nvar), lwd=2, col="#C0C0C066")
