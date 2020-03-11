#!/usr/bin/env R

library(ggplot2)
library(plotly)

args <- commandArgs(trailingOnly=T)

# args[1]   path to input data
# args[2]   path to output data
# args[3]   number of scaffolds

df <- read.table(args[1],header=F)
colnames(df)= c("cpg","cpos","snp","spos","dist")

df$cpos <- ceiling(df$cpos)
df$spos <- ceiling(df$spos)

# set the size of scaffolds and rank
chr <- aggregate(c(df$cpos, df$spos), by = list(c(as.character(df$cpg), as.character(df$snp))), max)
chr <- chr[order(-chr$x),]
chr$x1 <- cumsum(chr$x)
chr$x0 <- chr$x1 - chr$x
chr.nrow <- nrow(chr)

# set the number of gridlines to draw
grids <- ifelse(chr.nrow < as.numeric(args[3]), chr.nrow, as.numeric(args[3]))

# get the mid position of scaffolds for axis breaks
chr$mid <- chr$x1 - (chr$x/2)
chr.grid.values <- chr$mid[1:grids]
chr.grid.labels <- chr$Group.1[1:grids]

# set min and max
chr.min <- 0
chr.max <- chr$x1[-1]

# modify the values in df based on ranking of scaffolds
df$tip <- paste0("Meth (", df$cpg, ":", df$cpos-1, "-", df$cpos, ")<br>SNP (", df$snp, ":", df$spos-1, "-", df$spos, ")")
df$cpos <- df$cpos + chr$x0[match(df$cpg, chr$Group.1)]
df$spos <- df$spos + chr$x0[match(df$snp, chr$Group.1)]

p <- ggplot(df, aes(spos,cpos,color=dist,text=tip)) +
    # show all points
    geom_point(size=0.5) +

    # scale x and y axis
    scale_x_continuous(breaks = chr.grid.values, labels=chr.grid.labels) +
    scale_y_continuous(breaks = chr.grid.values, labels=chr.grid.labels) +
    expand_limits(x=c(chr.min,chr.max), y=c(chr.min,chr.max)) +

    # add plot and axis titles
    ggtitle(paste0(basename(args[2]))) +
    xlab("SNP position / bp") +
    ylab("Meth. position / bp") +

    # customise the theme
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.y=element_text(angle=90,hjust=0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
    )

# add gridlines
for(i in chr.grid.labels){
p <- p + geom_hline(yintercept=chr[chr$Group.1 == i,]$x1, color="grey", alpha=0.5)
p <- p + geom_vline(xintercept=chr[chr$Group.1 == i,]$x1, color="grey", alpha=0.5)
}

# save plot image
p.png <- paste0(args[2], ".png")
ggsave(p.png, p, width=12, height=7)

# save interactive plot
p <- ggplotly(p, tooltip="text")
p.html <- paste0(args[2], ".html")
p.zip <- paste0(args[2], ".zip")
htmlwidgets::saveWidget(p, p.html, selfcontained=F, libdir=args[2], title=basename(args[2]))
zip::zipr(p.zip, c(p.html,args[2]))