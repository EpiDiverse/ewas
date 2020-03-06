#!/usr/bin/env R

library(ggplot2)

args <- commandArgs(trailingOnly=T)

# args[1]   path to input data
# args[2]   number of chromosomes

df <- read.table(args[1],header=F)
colnames(df)= c("cpg","cpos","snp","spos","dist")

# set the size of scaffolds and rank
chr <- aggregate(c(df$cpos, df$spos), by = list(c(as.character(df$cpg), as.character(df$snp))), max)
chr <- chr[order(-chr$x),]
chr$x1 <- cumsum(chr$x)
chr$x0 <- chr$x1 - chr$x
chr.nrow <- nrow(chr)

# set the number of grids to draw
if(chr.nrow < as.numeric(args[2])){
    grids <- chr.nrow
} else {
    grids <- as.numeric(args[2])
}

# modify the values in df based on ranking of scaffolds
df$cpos <- df$cpos + chr$x0[match(df$cpg, chr$Group.1)]
df$spos <- df$spos + chr$x0[match(df$snp, chr$Group.1)]



pdf("scatter.pdf")
p <- ggplot(data = df) +

    geom_point(aes(cpos,spos, color=dist), show.legend = FALSE) +
    scale_x_continuous(breaks = chr$x1[1:grids],minor_breaks=NULL,labels=chr$Group.1[1:grids]) +
    theme(axis.text.y=element_text(angle=0, vjust=1),
    text=element_text(size=5),axis.text.x = element_text(angle = 90, hjust=1)) +
    xlab("SNP position / bp") +
    ylab("Meth. position / bp") 
    
dev.off()
