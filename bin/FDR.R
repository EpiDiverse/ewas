#!/usr/bin/env R

library(stats)

args <- commandArgs(trailingOnly=T)

raw.model <- read.table(args[1],header=F)
adj.model <- p.adjust(p = raw.model[-1,5], method = "BH", n = length(raw.model[-1,5]))

write.table(adj.model, args[2], sep = "\t", row.names = FALSE, quote = FALSE)