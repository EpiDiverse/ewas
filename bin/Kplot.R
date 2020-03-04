#!/usr/bin/env R

library(ggplot2)

args <- commandArgs(trailingOnly=T)

# args[1]   path to input data
# args[2]   meth.txt
# args[3]   snp.txt
# args[4]   cov.txt
# args[5]   number of Kplots

model <- read.table(args[1],header=T)
cpg <- read.table(args[2],header=T)
snp <- read.table(args[3],header=T)
cvrt <- 
nSamples <- 
nCov <- 


if(!is.null(topKplot) && is.numeric(topKplot) && topKplot > 0){
    if(topKplot > length(GxEmodel$all$eqtl$pvalue))
        topKplot <- length(GxEmodel$all$eqtl$pvalue)
    for(i in 1:topKplot){
        snpiname <- GxEmodel$all$eqtl$snps[i]
        cpginame <- GxEmodel$all$eqtl$gene[i]

        snpi<- snp$FindRow(snpiname)$row[1, ]
        cpgi<- cpg$FindRow(cpginame)$row[1, ]

        sceData <- data.frame(snp = factor(snpi, levels=c(1,2,3), labels=c("AA", "AB", "BB")),
                                cpg = as.numeric(cpgi),
                                env = as.numeric(cvrt$getSlice(1)[dim(cvrt)[1], ]))

        dp <- ggplot(sceData, aes_string(x="env", y="cpg")) + geom_point() +
            geom_smooth(method=lm, se=FALSE, colour = "red") +
            facet_grid(. ~ snp) + theme_bw() + xlab("Environmental factor") + ylab(paste0(cpginame, " profile")) +
            theme(strip.text.x = element_text(size=12, face="bold"))
        
        if(savePlot){
            ggsave(paste0(dirname(output_file_name), .Platform$file.sep,
                        snpiname, " x ", cpginame, ".png"),
                    dp, width=12, height=7)
        }else{
            print(dp)
        }
        
    }
}