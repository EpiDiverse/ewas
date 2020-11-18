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
environment(my.GEM_GWASmodel) <- asNamespace("GEM")
##################################

snp= paste(".", args[2], sep="/")
cov= paste(".", args[3], sep="/")
env= paste(".", args[2], sep="/")
p= as.numeric(args[5])

#main.nf takes p_value as an integer, but GEM process converts it into string. before this, this value be converted into int again. 
gwas_txt = paste("./output/", args[6], "GWAS.txt", sep="")
gwas_png = paste("./output/", args[7], "GWAS.png", sep="")

my.GEM_GWASmodel(env, snp, cov, p, gwas_txt, gwas_png)
