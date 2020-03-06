#!/usr/bin/env R

library(ggplot2)

args <- commandArgs(trailingOnly=T)

# args[1]   path to input data

model <- read.table(args[1],header=T)

dotplot <- nrow(model)

if(!is.null(dotplot) && is.numeric(dotplot) && dotplot > 0){


    

    snpiname <- as.character(model$snp[i])
    cpginame <- as.character(model$cpg[i])

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