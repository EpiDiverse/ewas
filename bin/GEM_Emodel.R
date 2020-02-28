#!/usr/bin/env R

suppressMessages(require(ggplot2))
suppressMessages(require(Rcpp))
suppressMessages(require(digest))
suppressMessages(require(devtools))
suppressMessages(require(usethis))
suppressMessages(require(GEM))

args <- commandArgs(trailingOnly=T)


env= paste(".", args[1], sep="/")
cov= paste(".", args[2], sep="/")
meth= paste(".", args[3], sep="/")
#meth= paste(".", "methylation.txt", sep="/")
p= as.numeric(args[4])
#main.nf takes p_value as an integer, but GEM process converts it into string. before this, this value be converted into int again. 
emodel_txt = paste("./", args[4], ".txt", sep="")
e_model_qq = paste("./", args[4], ".jpg", sep="")


#cov= paste(args[1], "cov.txt", sep = .Platform$file.sep)
#meth= paste(args[1], "methylation.txt", sep = .Platform$file.sep)
#Emodel_result= paste(args[1], "Result_Emodel.txt", sep = .Platform$file.sep)
#Emodel_result_qqplot= paste(args[1], "QQplot_Emodel.png", sep = .Platform$file.sep)
GEM_Emodel(env, cov, meth, p , emodel_txt, e_model_qq, savePlot=TRUE)

#write.table(cov,"cov.txt", quote=FALSE, row.names = FALSE, sep = "\t")
#write.table(covar,"covar.txt", quote=FALSE, row.names = FALSE, sep = "\t")
