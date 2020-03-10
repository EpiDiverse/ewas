#!/usr/bin/env R

library(ggplot2)

args <- commandArgs(trailingOnly=T)

# args[1]   path to input data
# args[2]   number of chromosomes
# args[3]   path to output data

df <- read.table(args[1],header=F)
colnames(df)= c("cpg","cpos","snp","spos","dist")

# set the size of scaffolds and rank
chr <- aggregate(c(df$cpos, df$spos), by = list(c(as.character(df$cpg), as.character(df$snp))), max)
chr <- chr[order(-chr$x),]
chr$x1 <- cumsum(chr$x)
chr$x0 <- chr$x1 - chr$x
chr.nrow <- nrow(chr)

# set the number of grids to draw
if(chr.nrow < as.numeric(args[3])){
    grids <- chr.nrow
} else {
    grids <- as.numeric(args[3])
}

chr$mid <- chr$x1 - (chr$x/2)
chr.grid.values <- chr$mid[1:grids]
chr.grid.labels <- chr$Group.1[1:grids]

chr.min <- 0
chr.max <- chr$x1[-1]

# modify the values in df based on ranking of scaffolds
df$cpos <- df$cpos + chr$x0[match(df$cpg, chr$Group.1)]
df$spos <- df$spos + chr$x0[match(df$snp, chr$Group.1)]

p <- ggplot(df, aes(cpos,spos,color=dist)) +

    geom_point(size=0.5) +
    scale_x_continuous(breaks = chr.grid.values,labels=chr.grid.labels) +
    scale_y_continuous(breaks = chr.grid.values,labels=chr.grid.labels) +
    expand_limits(x=c(chr.min,chr.max), y=c(chr.min,chr.max)) +
    theme_minimal() +
    theme(axis.text.y=element_text(angle=90,hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_line(size=1)) +
    xlab("SNP position / bp") +
    ylab("Meth. position / bp")

ggsave(paste0(args[2], ".png"), p, width=12, height=7)
