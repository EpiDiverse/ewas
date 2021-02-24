#!/usr/bin/env R

library(ggplot2)

args <- commandArgs(trailingOnly=T)

# args[1]   path to input data
# args[2]   meth.txt
# args[3]   snp.txt
# args[4]   gxe.txt

model <- read.table(args[1],header=T)
cpg <- read.table(args[2],header=T)
snp <- read.table(args[3],header=T)
gxe <- read.table(args[4],header=T)

topKplot <- nrow(model)

if(!is.null(topKplot) && is.numeric(topKplot) && topKplot > 0){
    if(topKplot > length(model$pvalue))
        topKplot <- length(model$pvalue)
    for(i in 1:topKplot){
        snpiname <- as.character(model$snp[i])
        cpginame <- as.character(model$ID[i])
        #cpginame <- as.character(model$cpg[i])

        snpi <- snp[snp[,1] == snpiname,2:ncol(snp)]
        cpgi <- cpg[cpg[,1] == cpginame,2:ncol(snp)]
        covi <- dim(gxe)

        sceData <- data.frame(snp = factor(snpi, levels=c(1,2,3), labels=c("AA", "AB", "BB")),
                                cpg = as.numeric(cpgi),
                                env = as.numeric(gxe[covi[1],2:covi[2]]))

        sceData <- sceData[complete.cases(sceData),]
        dp.name <- paste0(snpiname, "_vs_", cpginame)

        if(nrow(sceData) > 0){

            dp <- ggplot(sceData, aes_string(x="env", y="cpg")) + geom_point() +
                geom_smooth(method=lm, se=FALSE, colour = "red") +
                facet_grid(. ~ snp) + theme_bw() + xlab("Environmental factor") + ylab(paste0(cpginame, " profile")) +
                theme(strip.text.x = element_text(size=12, face="bold"))

            ggsave(paste0(dirname(args[1]), "/", dp.name, ".png"), dp, width=12, height=7)

        }else{
            write(paste0(dp.name, " comparison resulted in zero rows"), stderr())
        }
    }
}