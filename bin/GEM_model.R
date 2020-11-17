#' GEM_Emodel Analysis
#'
#' GEM_Emodel is to find the association between methylation and environmental factor genome widely.
#'
#' GEM_Emodel finds the association between methylation and environment genome-wide by performing matrix 
#' based iterative correlation and memory-efficient data analysis instead of millions of linear regressions 
#' (N = number_of_CpGs). The methylation data are the measurements for CpG probes, for example, 450,000 CpGs 
#' from Illumina Infinium HumanMethylation450 Array. The environmental factor can be a particular phenotype or environment 
#' factor from,for example, birth outcomes, maternal conditions or disease traits. The output of GEM_Emodel
#' for particular environmental factor is a list of CpGs that are potential epigenetic biomarkers.
#' GEM_Emodel runs linear regression like lm (M ~ E + covt), where M is a matrix with methylation data,
#' E is a matrix with environment factor and covt is a matrix with covariates, and all read from the
#' formatted text data file.
#'
#'
#' @param env_file_name Text file with rows representing environment factor and columns representing samples, such as the example data file "env.txt".
#' @param covariate_file_name Text file with rows representing covariate factors and columns representing samples, such as the example data file "cov.txt".
#' @param methylation_file_name Text file with rows representing methylation profiles for CpGs, and columns representing samples,  such as the example data file "methylation.txt".
#' @param Emodel_pv The pvalue cut off. Associations with significances at Emodel_pv level or below are saved to output_file_name, with corresponding estimate of effect size (slope coefficient), test statistics and p-value. Default value is 1.0.
#' @param outfile The result file with each row presenting a CpG and its association with environment, which contains CpGID, estimate of effect size (slope coefficient), test statistics, pvalue and FDR at each column.
#' @param noFDR toggle FDR calculation.
#'
#' @return save results automatically
#' @export
#'
#' @examples
#' DATADIR = system.file('extdata',package='GEM')
#' RESULTDIR = getwd()
#' env_file = paste(DATADIR, "env.txt", sep = .Platform$file.sep)
#' covariate_file = paste(DATADIR, "cov.txt", sep = .Platform$file.sep)
#' methylation_file = paste(DATADIR, "methylation.txt", sep = .Platform$file.sep)
#' Emodel_pv = 1
#' output_file = paste(RESULTDIR, "Result_Emodel.txt", sep = .Platform$file.sep)
#' qqplot_file = paste(RESULTDIR, "QQplot_Emodel.jpg", sep = .Platform$file.sep)
#' GEM_Emodel(env_file, covariate_file, methylation_file, Emodel_pv, output_file, qqplot_file)
my.GEM_Emodel <-
    function(env_file_name, covariate_file_name, methylation_file_name,
             Emodel_pv, outfile, noFDR = FALSE) {

        errorCovariance = numeric();

        env <- SlicedData$new();
        env$fileDelimiter = "\t";      # the TAB character
        env$fileOmitCharacters = "NA"; # denote missing values;
        env$fileSkipRows = 1;          # one row of column labels
        env$fileSkipColumns = 1;       # one column of row labels
        env$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        env$LoadFile(env_file_name);

        cvrt <- SlicedData$new();
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
        Emodel <- Matrix_eQTL_engine2(
            snps = env,
            gene = cpg,
            cvrt = cvrt,
            output_file_name = outfile,
            pvOutputThreshold = Emodel_pv,
            useModel = modelLINEAR,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            pvalue.hist = FALSE,
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = noFDR,
            addInfo = "CpGs"
        )

        ## Results:
        cat('Analysis done in: ', Emodel$time.in.sec, ' seconds', '\n');
        #show(Emodel$all$eqtls)
        #R2 = Emodel$all$eqtls$statistic ^ 2 / (Emodel$all$eqtls$statistic ^ 2 + Emodel$param$dfFull);
        
    }


/*
#use a kinship matrix from methylation data to conduct a mixed model analysis
my.mixed_GEM_Emodel <-
    function(env_file_name, covariate_file_name, methylation_file_name,
             Emodel_pv, outfile, noFDR = FALSE) {

        errorCovariance = numeric();

        env <- SlicedData$new();
        env$fileDelimiter = "\t";      # the TAB character
        env$fileOmitCharacters = "NA"; # denote missing values;
        env$fileSkipRows = 1;          # one row of column labels
        env$fileSkipColumns = 1;       # one column of row labels
        env$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        env$LoadFile(env_file_name);

        cvrt <- SlicedData$new();
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
        mixed_Emodel <- Matrix_eQTL_engine2(
            snps = env,
            gene = cpg,
            cvrt = cvrt,
            output_file_name = outfile,
            pvOutputThreshold = Emodel_pv,
            useModel = modelANOVA,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            pvalue.hist = FALSE,
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = noFDR,
            addInfo = "CpGs"
        )
        
        fstat = mixed_Emodel$all$eqtls$statistic;
        pvalue = mixed_Emodel$all$eqtls$pvalue;
        rez = c( Fstat = fstat, pvalue = pvalue)
        # And compare to those from ANOVA in R
        {
            #cat("\n\n Matrix eQTL: \n");
            #print(rez);
            #cat("\n R anova(lm()) output: \n")
            lmodel = lm( cpg.mat ~ cvrt.mat + factor(snps.mat), weights = 1/errorCovariance);
            lmout = anova(lmodel)[2, c("F value", "Pr(>F)")];
            #print( lmout )
}
        
# Results from Matrix eQTL and "lm" must agree
stopifnot(all.equal(lmout, rez, check.attributes = FALSE));
        
        
        ## Results:
        cat('Analysis done in: ', mixed_Emodel$time.in.sec, ' seconds', '\n');
        #show(Emodel$all$eqtls)
        #R2 = Emodel$all$eqtls$statistic ^ 2 / (Emodel$all$eqtls$statistic ^ 2 + Emodel$param$dfFull);
        
    }

*/

#' GEM_Gmodel Analysis
#'
#' GEM_Gmodel creates a methQTL genome-wide map.
#'
#' GEM_Gmodel creates a methQTL genome-wide map by performing matrix based iterative correlation and memory-efficient 
#' data analysis instead of millions of linear regressions (N = number_of_CpGs x number_of_SNPs)
#' between methylation and genotyping. Polymorphisms close to CpGs in the same chromosome (cis-) or different chromosome (trans-)
#' often form methylation quantitative trait loci (methQTLs) with CpGs. In GEM_Gmodel, MethQTLs can be discovered by correlating
#' single nucleotide polymorphism (SNP) data with CpG methylation from the same samples, by linear regression lm (M ~ G + covt),
#' where M is a matrix with methylation data, G is a matrix with genotype data and covt is a matrix with covariates,
#' and all read from the formatted text data file. The methylation data are the measurements for CpG probes, for example,
#' 450,000 CpGs from Illumina Infinium HumanMethylation450 Array. The genotype data are encoded as 1,2,3 or any three 
#' distinct values for major allele homozygote (AA),
#' heterozygote (AB) and minor allele homozygote (BB). The linear regression is adjusted by covariates read from covariate data file.
#' The output of GEM_Gmodel is a list of CpG-SNP pairs, where the SNP is the best fit to explain the particular CpG. The significant
#' association between CpG-SNP pair suggests the methylation driven by genotyping variants, which is so called methylation quantitative
#' trait loci (methQTL).
#'
#' @param snp_file_name Text file with rows representing genotype encoded as 1,2,3 or any three distinct values for major allele homozygote (AA), heterozygote (AB) and minor allele homozygote (BB) and columns representing samples, such as the example data file "snp.txt".
#' @param covariate_file_name Text file with rows representing covariate factors and columns representing samples, such as the example data file "cov.txt".
#' @param methylation_file_name Text file with rows representing methylation profiles for CpGs, and columns representing samples,  such as the example data file "methylation.txt".
#' @param Gmodel_pv The pvalue cut off. Associations with significances at Gmodel_pv level or below are saved to output_file_name, with corresponding estimate of effect size (slope coefficient), test statistics and p-value. Default value is 5.0E-08.
#' @param outfile The result file with each row presenting a CpG and its association with SNP, which contains CpGID, SNPID, estimate of effect size (slope coefficient), test statistics, pvalue and FDR at each column.
#' @param noFDR toggle FDR calculation.
#' 
#' @return save results automatically
#' @export
#'
#' @examples
#' DATADIR = system.file('extdata',package='GEM')
#' RESULTDIR = getwd()
#' snp_file = paste(DATADIR, "snp.txt", sep = .Platform$file.sep)
#' covariate_file = paste(DATADIR, "cov.txt", sep = .Platform$file.sep)
#' methylation_file = paste(DATADIR, "methylation.txt", sep = .Platform$file.sep)
#' Gmodel_pv = 1e-04
#' output_file = paste(RESULTDIR, "Result_Gmodel.txt", sep = .Platform$file.sep)
#' GEM_Gmodel(snp_file, covariate_file, methylation_file, Gmodel_pv, output_file)
my.GEM_Gmodel <-
    function(snp_file_name, covariate_file_name, methylation_file_name,
             Gmodel_pv, outfile, noFDR = FALSE)
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
            output_file_name = outfile,
            pvOutputThreshold = Gmodel_pv,
            useModel = modelLINEAR,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            pvalue.hist = FALSE,
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = noFDR,
            addInfo = "methQTL"
        );

        ## Results:
        cat('Analysis done in: ', Gmodel$time.in.sec, ' seconds', '\n');
        #R2 = Gmodel$all$eqtls$statistic ^ 2 / (Gmodel$all$eqtls$statistic ^ 2 + Gmodel$param$dfFull);

    }


#' GEM_GxEmodel
#'
#' GEM_GxEmodel tests the ability of the interaction of gene and environmental factor to predict DNA methylation level.
#'
#' GEM_GxEmodel explores how the genotype can work in interaction with environment (GxE) to influence specific DNA 
#' methylation level, by performing matrix based iterative correlation and memory-efficient data analysis 
#' among methylation, genotyping and environment. 
#' This has greatly released the computational burden for GxE study from billions of linear
#' regression (N = number_of_CpGs x number_of_SNPs x number_of_environment) and made it possible to be accomplished in an 
#' efficient way.
#' The linear regression is lm (M ~ G x E + covt), where M is a matrix with methylation data, G is a matrix with genotype data,
#' E is environment data and covt is covariate matrix. E values is combined into covariate file as the last row, and all read from the formatted text data file. The output of
#' GEM_GxEmodel is a list of CpG-SNP-Env triplets, where the environment factor segregated in genotype group fits to explain
#' the particular CpG. The significant association suggests the association between methylation and environment can be better
#' explained by segregation in genotype groups (GxE).
#'
#' @param snp_file_name Text file with rows representing genotype encoded as 1,2,3 or any three distinct values for major allele homozygote (AA), heterozygote (AB) and minor allele homozygote (BB) and columns representing samples, such as the example data file "snp.txt".
#' @param covariate_file_name Text file with rows representing covariate factors and the envirnoment value, and the environment value should be put in the last row, and columns representing samples, such as the example data file "gxe.txt".
#' @param methylation_file_name Text file with rows representing methylation profiles for CpGs, and columns representing samples,  such as the example data file "methylation.txt".
#' @param GxEmodel_pv The pvalue cut off. Associations with significances at GxEmodel_pv level or below are saved to output_file_name, with corresponding estimate of effect size (slope coefficient), test statistics and p-value. Default value is 5.0E-08.
#' @param output_file_name The result file with each row presenting a CpG and its association with SNPxEnv, which contains CpGID, SNPID, estimate of effect size (slope coefficient), test statistics, pvalue and FDR at each column.
#' @param topKplot The top number of topKplot CpG-SNP-Env triplets will be presented into charts to demonstrate how environment values segregated by SNP groups can explain methylation. 
#' @param noFDR toggle FDR calculation.
#' @param savePlot if save the plot.
#'    
#' @return save results automatically
#' @import ggplot2
#' @export
#'
#' @examples
#' DATADIR = system.file('extdata',package='GEM')
#' RESULTDIR = getwd()
#' snp_file = paste(DATADIR, "snp.txt", sep = .Platform$file.sep)
#' covariate_file = paste(DATADIR, "gxe.txt", sep = .Platform$file.sep)
#' methylation_file = paste(DATADIR, "methylation.txt", sep = .Platform$file.sep)
#' GxEmodel_pv = 1e-4
#' output_file = paste(RESULTDIR, "Result_GxEmodel.txt", sep = .Platform$file.sep)
#' GEM_GxEmodel(snp_file, covariate_file, methylation_file, GxEmodel_pv, output_file)
my.GEM_GxEmodel <-
    function(snp_file_name, covariate_file_name, methylation_file_name,
             GxEmodel_pv, outfile, topKplot = 10, noFDR = FALSE, savePlot=TRUE)
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
        GxEmodel = Matrix_eQTL_engine2(
            snps = snp,
            gene = cpg,
            cvrt = cvrt,
            output_file_name = outfile,
            pvOutputThreshold = GxEmodel_pv,
            useModel = modelLINEAR_CROSS,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            pvalue.hist = FALSE,
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = noFDR,
            addInfo = "cpg-snp pairs"
        );

        ## Results:
        cat('Analysis done in: ', GxEmodel$time.in.sec, ' seconds', '\n');
        #R2 = GxEmodel$all$eqtls$statistic ^ 2 / (GxEmodel$all$eqtls$statistic ^ 2 + GxEmodel$param$dfFull);

        #jpeg(qqplot_file_name, width = 2000, height = 2000, res = 300)
        #plot(GxEmodel, pch = 16, cex = 0.7);
        #dev.off()
}


#' GEM_GWASmodel
#' 
#' 
#' GEM_GWASmodel performs genome wide association study (GWAS).
#' 
#' GEM_GWASmodel finds the association between genetic variants and environment genome-wide by performing matrix 
#' based iterative correlation and memory-efficient data analysis instead of millions of linear regressions 
#' (N = number_of_SNPs). The environmental factor can be a particular phenotype or environment 
#' factor from,for example, birth outcomes, maternal conditions or disease traits. 
#' The genotype data are encoded as 1,2,3 or any three 
#' distinct values for major allele homozygote (AA),
#' heterozygote (AB) and minor allele homozygote (BB). The linear regression is adjusted by covariates read from covariate data file.
#' The output of GEM_GWASmodel is a list of SNPs and their association with environment. 
#' GEM_GWASmodel runs linear regression like lm (E ~ G + covt), where G is a matrix with genotype data,
#' E is a matrix with environment factor and covt is a matrix with covariates, and all read from the
#' formatted text data file.
#' 
#' 
#' @export
#' @param env_file_name Text file with rows representing environment factor and columns representing samples, such as the example data file "env.txt".
#' @param snp_file_name Text file with rows representing genotype encoded as 1,2,3 or any three distinct values for major allele homozygote (AA), heterozygote (AB) and minor allele homozygote (BB) and columns representing samples, such as the example data file "snp.txt".
#' @param covariate_file_name Text file with rows representing covariate factors,  and columns representing samples, such as the example data file "cov.txt".
#' @param GWASmodel_pv The pvalue cut off. Associations with significances at GWASmodel_pv level or below are saved to output_file_name, with corresponding estimate of effect size (slope coefficient), test statistics and p-value. Default value is 5.0E-08.
#' @param output_file_name The result file with each row presenting a SNP and its association with environment, which contains SNPID, estimate of effect size (slope coefficient), test statistics, pvalue and FDR at each column.
#' @param qqplot_file_name Output QQ plot for all pvalues.
#'
#' @return save results automatically
#' @examples
#' DATADIR = system.file('extdata',package='GEM')
#' RESULTDIR = getwd()
#' snp_file = paste(DATADIR, "snp.txt", sep = .Platform$file.sep)
#' covariate_file = paste(DATADIR, "cov.txt", sep = .Platform$file.sep)
#' env_file = paste(DATADIR, "env.txt", sep = .Platform$file.sep)
#' GWASmodel_pv = 1
#' GEM_GWASmodel(env_file, snp_file, covariate_file, GWASmodel_pv, output_file, qqplot_file)

GEM_GWASmodel <-
    function(env_file_name, snp_file_name, covariate_file_name,
             GWAS_pv, outfile, noFDR = FALSE)
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

        env = SlicedData$new();
        env$fileDelimiter = "\t";      # the TAB character
        env$fileOmitCharacters = "NA"; # denote missing values;
        env$fileSkipRows = 1;          # one row of column labels
        env$fileSkipColumns = 1;       # one column of row labels
        env$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        env$LoadFile(env_file_name);


        ## Run the analysis
        GWASmodel = Matrix_eQTL_engine2(
            snps = snp,
            gene = env,
            cvrt = cvrt,
            output_file_name = outfile,
            pvOutputThreshold = GWASmodel_pv,
            useModel = modelLINEAR,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            pvalue.hist = FALSE,
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = noFDR,
            addInfo = "SNPs"
        );

        ## Results:
        cat('Analysis done in: ', GWAS$time.in.sec, ' seconds', '\n');
