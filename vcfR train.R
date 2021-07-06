setwd("~/R training/vcfR")
library(tidyverse)
library(vcfR)
library(pinfsc50)
list.files()

vcf_21 <- read.vcfR(list.files()[[1]], verbose = TRUE, nrows=20000, cols=1:11)

#Explore subsets of loaded vcf file
fix_own <- vcf_21@fix
gt_own <- vcf_21@gt

## Reading annotation file (GFF) which was downloaded from USCS
## Filter Chr21
gff_own <- read.table("hg19.knownGene.gtf", sep="\t", quote="")
chr21_ann <- subset(gff_own, V1 == "chr21")

##Reading DNA (FASTA) reference file for chr 21, from USCS browse, using ape?
dna <- ape::read.dna("chr21.fa", format = "fasta")

## Create Chrom file
chrom <- create.chromR(name='chr21', vcf=vcf_21, ann=chr21_ann, seq = dna)
chrom <- masker(chrom, min_DP = 10, max_DP = 50000)
plot(chrom)

# Plot read depth and annotation, MQ and QUAL is missing in vcf from 1000 genomes
chromoqc(chrom)
chrom <- proc.chromR(chrom, verbose = TRUE)
chromoqc(chrom, xlim=c(0.5e+07, 1.2e+07)) ## Boundaries of loaded data
plot(chrom)

## Checking data
queryMETA(vcf_21)
vcf_21
head(vcf_21)
head(is.polymorphic(vcf_21, na.omit = TRUE))
vcf2 <- extract.indels(vcf_21, return.indels = TRUE)
?extract.indels
head(vcf2)

## INto tibble
Z <- vcfR2tidy(vcf_21, format_fields = c("GT", "DP"))
head(Z)
plot(chrom)

head(getFIX(vcf_own))
head(is.polymorphic(vcf_own, na.omit = TRUE))
vcf_own[1:4,]
browseVignettes('vcfR')

## ----------------------------------- READ FILES-------------------
pkg <- "pinfsc50"

vcf <- read.vcfR(vcf_file, verbose = FALSE )
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="")

head(vcf)

chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)
chrom <- masker(chrom, min_QUAL = 1, min_DP = 300, max_DP = 700, min_MQ = 59.9,  max_MQ = 60.1)


data(vcfR_example)
strwrap(vcf@meta[1:7])
queryMETA(vcf)
queryMETA(vcf, element = 'DP')

vcf@gt[1:6, 1:4]
vcf <- read.vcfR("myVCFdata.vcf.gz")

##-------------------------------------------------------------------
pkg <- "pinfsc50"
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = pkg)
vcf <- read.vcfR(vcf_file, verbose = TRUE )
show(vcf)

head(is.polymorphic(vcf, na.omit = TRUE))
head(is.biallelic(vcf))
vcf2 <- extract.indels(vcf, return.indels=TRUE)
head(vcf2)
View(vcf2)
View(as_tibble(getFIX(vcf2)))

vcf[,1:5]
View(vcf)

## -----------------------------MEMORY-----------------------------

library('memuse')
nvar <- 10^{2:8}
nMb <- howbig(nrow=10, ncol=10, unit="MB")@size
for(i in nvar){
  nMb <- c(nMb, howbig(nrow=i, ncol=10, unit="MB")@size)
}
nvar <- c(10, nvar)
par(mar=c(5,5,4,2))
plot(log10(nvar), log10(nMb), xaxt="n", yaxt="n", type='b', xlab="Number of variants", ylab = "")
axis(side = 1, at = log10(nvar), labels=nvar)
axis(side = 2, at = log10(nMb), labels=nMb, las=2)
title(ylab="Memory use (Mb)", line=4)
abline(h=log10(nMb), lwd=2, col="#C0C0C066")
abline(v=log10(nvar), lwd=2, col="#C0C0C066")

zgrep -v "^#" myVcfFile.vcf.gz | head -n 1000000 | wc -l
?read.vcfR
