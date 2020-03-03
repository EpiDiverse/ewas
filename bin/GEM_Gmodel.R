#!/usr/bin/env R

suppressMessages(require(ggplot2))
suppressMessages(require(Rcpp))
suppressMessages(require(digest))
suppressMessages(require(devtools))
suppressMessages(require(usethis))
suppressMessages(require(GEM))

##################################
GEM_Gmodel <-
    function(snp_file_name, covariate_file_name, methylation_file_name,
             Gmodel_pv, output_file_name, noFDR = FALSE)
    {

        errorCovariance = numeric();

        snp = SlicedData$new();
        snp$fileDelimiter = "\t";      # the TAB character
        snp$fileOmitCharacters = "NA"; # denote missing values;
        snp$fileSkipRows = 1;          # one row of column labels
        snp$fileSkipColumns = 1;       # one column of row labels
        snp$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        snp$LoadFile(snp_file_name);

        cvrt = SlicedData$new();
        cvrt$fileDelimiter = "\t";      # the TAB character
        cvrt$fileOmitCharacters = "NA"; # denote missing values;
        cvrt$fileSkipRows = 1;          # one row of column labels
        cvrt$fileSkipColumns = 1;       # one column of row labels
        if (length(covariate_file_name) > 0) {
            cvrt$LoadFile(covariate_file_name);
        }

        cpg = SlicedData$new();
        cpg$fileDelimiter = "\t";      # the TAB character
        cpg$fileOmitCharacters = "NA"; # denote missing values;
        cpg$fileSkipRows = 1;          # one row of column labels
        cpg$fileSkipColumns = 1;       # one column of row labels
        cpg$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        cpg$LoadFile(methylation_file_name);


        ## Run the analysis
        Gmodel = Matrix_eQTL_engine2(
            snps = snp,
            gene = cpg,
            cvrt = cvrt,
            output_file_name = NULL,
            pvOutputThreshold = Gmodel_pv,
            useModel = modelLINEAR,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            pvalue.hist = "qqplot",
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = noFDR,
            addInfo = "methQTL"
        );

        unlink(output_file_name);
        ## Results:
        cat('Analysis done in: ', Gmodel$time.in.sec, ' seconds', '\n');
        #R2 = Gmodel$all$eqtls$statistic ^ 2 / (Gmodel$all$eqtls$statistic ^ 2 + Gmodel$param$dfFull);
        
        result_Gmodel <- cbind(
          as.character(Gmodel$all$eqtls$gene),
          as.character(Gmodel$all$eqtls$snps),
          Gmodel$all$eqtls$beta,
          Gmodel$all$eqtls$statistic,
          Gmodel$all$eqtls$pvalue,
          Gmodel$all$eqtls$FDR
        )
        colnames(result_Gmodel) <- c("cpg", "snp", "beta", "stats", "pvalue", "FDR")
        
        write.table(
            result_Gmodel, output_file_name, sep = "\t", row.names = FALSE, quote = FALSE
        )
    }

assignInNamespace("GEM_Gmodel", GEM_Gmodel, "GEM")
##################################


args <- commandArgs(trailingOnly=T)
#source(args[6])

snp= paste(".", args[1], sep="/")
cov= paste(".", args[2], sep="/")
meth= paste(".", args[3], sep="/")
p= as.numeric(args[4])

#main.nf takes p_value as an integer, but GEM process converts it into string. before this, this value be converted into int again. 
gmodel_txt = paste("./", args[5], ".txt", sep="")


GEM_Gmodel(snp, cov, meth, p, gmodel_txt, noFDR=FALSE)
