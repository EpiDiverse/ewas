#!/usr/bin/env R

library(stats)

args <- commandArgs(trailingOnly=T)

model <- read.table(args[1],header=F)
colnames(model) <- c("snp", "cpg", "beta", "stats", "pvalue")

model$FDR <- p.adjust(p = model[-1,5], method = "BH", n = length(model[-1,5]))

write.table(model[,c("cpg","snp","beta","stats","pvalue","FDR")], args[2], sep = "\t", row.names = FALSE, quote = FALSE)