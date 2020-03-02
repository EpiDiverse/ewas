#!/usr/bin/env R

suppressMessages(require(ggplot2))
suppressMessages(require(Rcpp))
suppressMessages(require(digest))
suppressMessages(require(devtools))
suppressMessages(require(usethis))
suppressMessages(require(GEM))

args <- commandArgs(trailingOnly=T)


snp= paste(".", args[1], sep="/")
cov= paste(".", args[2], sep="/")
meth= paste(".", args[3], sep="/")
p= as.numeric(args[4])

#main.nf takes p_value as an integer, but GEM process converts it into string. before this, this value be converted into int again. 
gmodel_txt = paste("./", args[5], ".txt", sep="")


GEM_Gmodel(snp, cov, meth, p, gmodel_txt)
