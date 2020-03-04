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
source(args[1])
environment(my.GEM_Gmodel) <- asNamespace("GEM")
##################################

snp= paste(".", args[2], sep="/")
cov= paste(".", args[3], sep="/")
meth= paste(".", args[4], sep="/")
p= as.numeric(args[5])

#main.nf takes p_value as an integer, but GEM process converts it into string. before this, this value be converted into int again. 
gmodel_txt = paste("./", args[6], ".txt", sep="")

my.GEM_Gmodel(snp, cov, meth, p, gmodel_txt, noFDR=FALSE)
