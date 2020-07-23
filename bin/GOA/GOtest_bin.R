args<-commandArgs(TRUE)

library("AnnotationDbi")
BiocManager::install("GSEABase")
library("GSEABase")
BiocManager::install("GOstats")
library("GOstats")
BiocManager::install("Rgraphviz")
library("Rgraphviz")

input_file <- as.character(args[1])
cat(paste("input file: ",args[1],"\n",sep=""))

species <- as.character(args[2])
cat(paste("species: ",args[2],"\n",sep=""))

output_OverTree_file <- as.character(args[3])
cat(paste("output OverTree file: ",args[3],"\n",sep=""))

output_OverTable_file <- as.character(args[4])
cat(paste("output OverTable file: ",args[4],"\n",sep=""))

output_UnderTree_file <- as.character(args[5])
cat(paste("output UnderTree file: ",args[5],"\n",sep=""))

output_UnderTable_file <- as.character(args[6])
cat(paste("output UnderTable file: ",args[6],"\n",sep=""))

ontology_option <- as.character(args[7])
cat(paste("ontology option: ",args[7],"\n",sep=""))

pvalueCutoff_option <- as.numeric(args[8])
cat(paste("pvalueCutoff option: ",args[8],"\n",sep=""))

parser= ArgumentParser(description='enter species name')
parser$add_argument('--species', type=string)


load(paste0("${baseDir}/db",species,"/frame.RData"))
load(paste0("${baseDir}/db",species,"/goframeData.RData"))
load(paste0("${baseDir}/db",species,"/goFrame.RData"))
load(paste0("${baseDir}/db",species,"/goAllFrame.RData"))
load(paste0("$${baseDir}/db",species,"/gsc.RData"))
load(paste0("/${baseDir}/db",species,"/universe.RData"))
source("${baseDir}/bin/GOA/GOtest.R")

#parse options
data <- read.table(input_file,header=F,sep="\t")
file_tree_name_over<-output_OverTree_file
file_name_over<-output_OverTable_file
file_tree_name_under<-output_UnderTree_file
file_name_under<-output_UnderTable_file
ontology<-ontology_option
pvalueCutoff<-pvalueCutoff_option

#make the GOtest
make_GO_tests(data,file_tree_name_over,file_name_over,file_tree_name_under,file_name_under,ontology,pvalueCutoff)

quit()
