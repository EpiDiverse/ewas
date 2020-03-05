#!/usr/bin/env R

suppressMessages(require(ggplot2))
suppressMessages(require(Rcpp))
suppressMessages(require(digest))
suppressMessages(require(devtools))
suppressMessages(require(usethis))
#suppressMessages(require(GEM))
#loadNamespace(GEM)

##################################
args <- commandArgs(trailingOnly=T)
source(paste(args[1],"GEM_model.R",sep="/"))
environment(my.GEM_GxEmodel) <- asNamespace("GEM")
##################################

snp= paste(".", args[2], sep="/")
cov= paste(".", args[3], sep="/")
meth= paste(".", args[4], sep="/")
p= as.numeric(args[5])

#main.nf takes p_value as an integer, but GEM process converts it into string. before this, this value be converted into int again. 
gxe_txt = paste("./", args[6], ".txt", sep="")

my.GEM_GxEmodel(snp, cov, meth, p, gxe_txt, topKplot=0, noFDR=TRUE, savePlot=FALSE)

