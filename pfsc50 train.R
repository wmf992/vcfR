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

## Extract information
## GT extract
data("vcfR_test")
gt <- extract.gt(vcf, IDtoRowNames = TRUE)
head(gt)

## Extracting seq depth for each variant, watch not to convert 
## genotype to numeric as this doesnt make sense
gt <- extract.gt(vcf, element = 'DP', as.numeric = TRUE)

## Extract Haplo quality
gt <- extract.gt(vcfR_test, element = 'HQ')
myHQ1 <- masplit(gt[,1:3], sort = 0)

## Tidy vcfR
## Checking field names
vcf_field_names(vcfR_test, tag = "FORMAT")
vcf_field_names(vcf, tag = "FORMAT")
Z <- vcfR2tidy(vcfR_test, format_fields = c("GT", "DP"))

#ChromR objects

#Location for data files from pinfsc50
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = pkg)
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = pkg)
gff_file <- system.file("extdata", "pinf_sc50.gff", package = pkg)

#read files into r
vcf <- read.vcfR( vcf_file, verbose = FALSE )
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="")
str(dna)

chrom <- create.chromR(name='Supercontig', vcf=vcf,
                       seq=dna, ann=gff, verbose = TRUE)
plot(chrom)
chromoqc(chrom)
chrom <- proc.chromR(chrom, verbose = TRUE)
plot(chrom)
chromoqc(chrom)
chrom <- masker(chrom, 
                min_QUAL=0, 
                min_DP=350, 
                max_DP=650, 
                min_MQ=59.5, 
                max_MQ=60.5)
chrom <- proc.chromR(chrom, verbose = TRUE)
chromoqc(chrom)

#Chaning windowsize for variants pr site
chrom <- proc.chromR(chrom, verbose = FALSE,
                     win.size = 1e3)
chromoqc(chrom)


## Extracting seq depth
dp <- extract.gt(chrom, element="DP",
                 as.numeric=TRUE)
rownames(dp) <- 1:nrow(dp)
heatmap.bp(dp[1:1000,])
is.na(dp[na.omit(dp==0)]) <- TRUE
heatmap.bp(dp[1001:1500,])

par(mar=c(8,4,4,2))
barplot(apply(dp, MARGIN = 2,
              mean, na.rm= TRUE,
              las=3))
par(mar=c(5,4,4,2))

##Reading VCF data
vcf <- read.vcfR(vcf_file, verbose = FALSE)
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)

#boxplot
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Depth", las=2)
abline(h=seq(0,1e4, by=100), col="#C0C0C088")
par(mar=c(5,4,4,2))

##Violin plot
if( require(reshape2) & require(ggplot2) ){
  dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
  dpf <- dpf[ dpf$Depth > 0,]
  p <- ggplot(dpf, aes(x=Sample, y=Depth)) + geom_violin(fill="#C0C0C0", adjust=1.0,
                                                         scale = "count", trim=TRUE)
  p <- p + theme_bw()
  p <- p + theme(axis.title.x = element_blank(), 
                 axis.text.x = element_text(angle = 60, hjust = 1, size=12))
  #  p <- p + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black")
  p <- p + scale_y_continuous(trans=scales::log2_trans(), 
                              breaks=c(1, 10, 100, 800),
                              minor_breaks=c(1:10, 2:10*10, 2:8*100))
  p <- p + theme(axis.title.y = element_text(size=12))
  p <- p + theme( panel.grid.major.y=element_line(color = "#A9A9A9", size=0.6) )
  p <- p + theme( panel.grid.minor.y=element_line(color = "#C0C0C0", size=0.2) )
  p <- p + stat_summary(fun.y=median, geom="point", shape=23, size=2)
  p
} else {
  message("The packages reshape2 and ggplot2 are required for this example but do not appear
          to be installed.  Please use install.packages(c('reshape2', 'ggplot2', 'scales')) if you would
          like to install them.")
}
head(vcf@gt[,-1])
