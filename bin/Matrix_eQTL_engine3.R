#!/usr/bin/env R

Matrix_eQTL_engine3 = function(
    snps,
    gene,
    cvrt = SlicedData$new(),
    output_file_name = "",
    pvOutputThreshold = 1e-5,
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = TRUE,
    output_file_name.cis = "",
    pvOutputThreshold.cis = 0,
    snpspos = NULL,
    genepos = NULL,
    cisDist = 1e6,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE,
    addInfo = "") {
    ################################# Basic variable checks #################################
    {
        # status("Performing basic checks of the input variables");
        stopifnot( "SlicedData" %in% class(gene) );
        stopifnot( any(c("SlicedData","SlicedData.fmt") %in% class(snps)) );
        stopifnot( "SlicedData" %in% class(cvrt) );

        # Check dimensions
        if( min(snps$nRows(),snps$nCols()) == 0 )
            stop("Empty genotype dataset");
        if( min(gene$nRows(),gene$nCols()) == 0 )
            stop("Empty expression dataset");
        if( snps$nCols() != gene$nCols() )
            stop("Different number of samples in the genotype and gene expression files");
        if( cvrt$nRows()>0 ) {
            if( snps$nCols() != cvrt$nCols() )
                stop("Wrong number of samples in the matrix of covariates");
        }

        stopifnot( class(pvOutputThreshold) == "numeric" );
        stopifnot( length(pvOutputThreshold) == 1 );
        stopifnot( pvOutputThreshold >= 0 );
        stopifnot( pvOutputThreshold <= 1 );

        stopifnot(  class(noFDRsaveMemory) == "logical" );
        stopifnot( length(noFDRsaveMemory) == 1 );

        #if( pvOutputThreshold > 0 ) {
            #stopifnot( !((length(output_file_name) == 0) && noFDRsaveMemory) )
            #stopifnot( length(output_file_name) <= 1 );
            #if( length(output_file_name) == 1 ) {
                #stopifnot( class(output_file_name) %in% c("character","connection") );
            #}
        #}

        stopifnot( class(pvOutputThreshold.cis) == "numeric" );
        stopifnot( length(pvOutputThreshold.cis) == 1 );
        stopifnot( pvOutputThreshold.cis >= 0 );
        stopifnot( pvOutputThreshold.cis <= 1 );
        stopifnot( !((pvOutputThreshold > 0) & (pvOutputThreshold.cis > 0) & (pvOutputThreshold > pvOutputThreshold.cis)) );
        stopifnot( (pvOutputThreshold > 0) | (pvOutputThreshold.cis > 0) );

        stopifnot( class(useModel) == class(modelLINEAR) );
        stopifnot( length(useModel) == 1 );
        stopifnot( useModel %in% c(modelLINEAR, modelANOVA, modelLINEAR_CROSS) );
        if( useModel %in%  c(modelLINEAR, modelLINEAR_CROSS) ) {
            if( snps$nCols() <= cvrt$nRows() + 1 + 1) {
                stop("The number of covariates exceeds the number of samples.\nLinear regression can not be fit.")
            }
        }
        if( useModel == modelLINEAR_CROSS ) {
            if( cvrt$nRows() == 0 ) {
                stop( "Model \"modelLINEAR_CROSS\" requires at least one covariate" );
            }
        }
        if( useModel == modelANOVA ) {
            n.anova.groups = getOption("MatrixEQTL.ANOVA.categories", 3);
            stopifnot( n.anova.groups == floor(n.anova.groups) );
            stopifnot( n.anova.groups >= 3 );
            # 			stopifnot( n.anova.groups < snps$nCols() - cvrt$nRows() - 2 );
            if( snps$nCols() <= cvrt$nRows() + n.anova.groups) {
                stop("The number of covariates exceeds the number of samples.\nLinear regression (ANOVA) can not be fit.")
            }
        }

        stopifnot(  class(verbose) == "logical" );
        stopifnot( length(verbose) == 1 );

        stopifnot(  class(min.pv.by.genesnp) == "logical" );
        stopifnot( length(min.pv.by.genesnp) == 1 );

        if( pvOutputThreshold.cis > 0 ) {
            stopifnot( !((length(output_file_name.cis) == 0) && noFDRsaveMemory) )
            stopifnot( length(output_file_name.cis) <= 1 );
            if( length(output_file_name.cis) == 1 ) {
                stopifnot( class(output_file_name.cis) %in% c("character","connection") );
            }

            # 			stopifnot( class(output_file_name.cis) == "character" );
            # 			stopifnot( length(output_file_name.cis) == 1 );
            stopifnot( class(snpspos) == "data.frame" );
            stopifnot( ncol(snpspos) == 3 );
            stopifnot( nrow(snpspos) > 0 );
            stopifnot( class(snpspos[1,3]) %in% c("integer", "numeric") )
            stopifnot( !any(is.na(snpspos[,3])) )
            stopifnot( class(genepos) == "data.frame" );
            stopifnot( ncol(genepos) == 4 );
            stopifnot( nrow(genepos) > 0 );
            stopifnot( class(genepos[1,3]) %in% c("integer", "numeric") )
            stopifnot( class(genepos[1,4]) %in% c("integer", "numeric") )
            stopifnot( !any(is.na(genepos[,3])) )
            stopifnot( !any(is.na(genepos[,4])) )
            stopifnot( nzchar(output_file_name.cis) )
        }

        if( pvOutputThreshold > 0 ) {
            stopifnot( nzchar(output_file_name) )
        }

        stopifnot( class(errorCovariance) %in% c("numeric", "matrix") );
        errorCovariance = as.matrix(errorCovariance);
        if(length(errorCovariance)>0) {
            if( nrow(errorCovariance) != ncol(errorCovariance) ) {
                stop("The covariance matrix is not square");
            }
            if( nrow(errorCovariance) != snps$nCols() ) {
                stop("The covariance matrix size does not match the number of samples");
            }
            if( !all(errorCovariance == t(errorCovariance)) ) {
                stop("The covariance matrix is not symmetric");
            }
        }
    }
    ################################# Initial setup #########################################
    {
        gene.std = .listBuilder$new();
        snps.std = .listBuilder$new();

        dont.clone.gene = getOption("MatrixEQTL.dont.preserve.gene.object", FALSE)
        if(is.null(dont.clone.gene))
            dont.clone.gene = FALSE;

        if( !dont.clone.gene )
            gene = gene$Clone();
        # snps = snps$Clone(); # snps is read only
        cvrt = cvrt$Clone();

        params = list(
            output_file_name = output_file_name,
            pvOutputThreshold = pvOutputThreshold,
            useModel = useModel,
            errorCovariance = errorCovariance ,
            verbose = verbose,
            output_file_name.cis = output_file_name.cis,
            pvOutputThreshold.cis = pvOutputThreshold.cis,
            cisDist = cisDist ,
            pvalue.hist = pvalue.hist,
            min.pv.by.genesnp = min.pv.by.genesnp);

        if( verbose ) {
            lastTime = 0;
            status <- function(text) {
                # gc();
                newTime = proc.time()[3];
                if(lastTime != 0) {
                    cat("Task finished in ", newTime-lastTime, " seconds\n");
                }
                #cat(text,"\n");
                lastTime <- newTime;
                unused = flush.console();
            }
        } else {
            status = function(text){}
        }
        start.time = proc.time()[3];
    }
    ################################# Error covariance matrix processing ####################
    {
        if( length(errorCovariance) > 0 ) {
            status("Processing the error covariance matrix");
            eig = eigen(errorCovariance, symmetric = TRUE)
            d = eig$values;
            v = eig$vectors;
            #  errorCovariance == v %*% diag(d) %*% t(v)
            #  errorCovariance^0.5 == v*sqrt(d)*v" (no single quotes anymore)
            #  errorCovariance^(-0.5) == v*diag(1./sqrt(diag(d)))*v"
            if( any(d<=0) ) {
                stop("The covariance matrix is not positive definite");
            }
            correctionMatrix = v %*% diag(1./sqrt(d)) %*% t(v);
            rm( eig, v, d, errorCovariance )
        } else {
            rm( errorCovariance );
            correctionMatrix = numeric();
        }
    }
    ################################# Matching gene and SNPs locations ######################
    if( pvOutputThreshold.cis > 0 ) {
        status("Matching data files and location files")

        # names in the input data
        gene_names = rownames(gene);
        snps_names = rownames(snps);

        # gene range, set: left<right
        if(any(genepos[,3] > genepos[,4])) {
            temp3 = genepos[,3];
            temp4 = genepos[,4];
            genepos[,3] = pmin(temp3,temp4);
            genepos[,4] = pmax(temp3,temp4);
            rm(temp3, temp4);
        }

        # match with the location data
        genematch = match( gene_names, genepos[ ,1],  nomatch = 0L);
        usedgene = matrix(FALSE, nrow(genepos), 1); # genes in "genepos" that are matching  "gene_names"
        usedgene[ genematch ] = TRUE;
        if( !any(genematch) ) {
            stop("Gene names do not match those in the gene location file.");
        }
        cat( sum(genematch>0), "of", length(gene_names), " genes matched\n");


        snpsmatch = match( snps_names, snpspos[ ,1],  nomatch = 0L);
        usedsnps = matrix(FALSE, nrow(snpspos),1);
        usedsnps[ snpsmatch ] = TRUE;
        if( !any(snpsmatch) ) {
            stop("SNP names do not match those in the SNP location file.");
        }
        cat( sum(snpsmatch>0), "of", length(snps_names), " SNPs matched\n");

        # list used chr names
        chrNames = unique(c( as.character(unique(snpspos[usedsnps,2])),
                             as.character(unique(genepos[usedgene,2])) ))
        chrNames = chrNames[ sort.list( suppressWarnings(as.integer(chrNames)),
                                        method = "radix", na.last = TRUE ) ];
        # match chr names
        genechr = match(genepos[,2],chrNames);
        snpschr = match(snpspos[,2],chrNames);

        # max length of a chromosome
        chrMax = max( snpspos[usedsnps, 3], genepos[usedgene, 4], na.rm = TRUE) + cisDist;

        # Single number location for all rows in "genepos" and "snpspos"
        genepos2 = as.matrix(genepos[ ,3:4, drop = FALSE] + (genechr-1)*chrMax);
        snpspos2 = as.matrix(snpspos[ ,3  , drop = FALSE] + (snpschr-1)*chrMax);

        # the final location arrays;
        snps_pos = matrix(0,length(snps_names),1);
        snps_pos[snpsmatch>0, ] = snpspos2[snpsmatch, , drop = FALSE];
        snps_pos[rowSums(is.na(snps_pos))>0, ] = 0;
        snps_pos[snps_pos==0] = (length(chrNames)+1) * (chrMax+cisDist);
        rm(snps_names, snpsmatch, usedsnps, snpschr, snpspos2)

        gene_pos = matrix(0,length(gene_names),2);
        gene_pos[genematch>0, ] = genepos2[genematch, , drop = FALSE];
        gene_pos[rowSums(is.na(gene_pos))>0, ] = 0;
        gene_pos[gene_pos==0] = (length(chrNames)+2) * (chrMax+cisDist);
        rm(gene_names, genematch, usedgene, genechr, genepos2)
        rm(chrNames, chrMax);

        if( is.unsorted(snps_pos) ) {
            status("Reordering SNPs\n");
            ordr = sort.list(snps_pos);
            snps$RowReorder(ordr);
            snps_pos = snps_pos[ordr, , drop = FALSE];
            rm(ordr);
        }
        if( is.unsorted(rowSums(gene_pos)) ) {
            status("Reordering genes\n");
            ordr = sort.list(rowSums(gene_pos));
            gene$RowReorder(ordr);
            gene_pos = gene_pos[ordr, , drop = FALSE];
            rm(ordr);
        }

        # Slice it back.
        geneloc = vector("list", gene$nSlices())
        gene_offset = 0;
        for(gc in 1:gene$nSlices()) {
            nr = gene$GetNRowsInSlice(gc);
            geneloc[[gc]] = gene_pos[gene_offset + (1:nr), , drop = FALSE];
            gene_offset = gene_offset + nr;
        }
        rm(gc, gene_offset, gene_pos);

        snpsloc = vector("list", snps$nSlices())
        snps_offset = 0;
        for(sc in 1:snps$nSlices()) {
            nr = snps$GetNRowsInSlice(sc);
            snpsloc[[sc]] = snps_pos[snps_offset + (1:nr), , drop = FALSE];
            snps_offset = snps_offset + nr;
        }
        rm(nr, sc, snps_offset, snps_pos);
    }
    ################################# Covariates processing #################################
    {
        status("Processing covariates");
        if( useModel == modelLINEAR_CROSS ) {
            last.covariate = as.vector(tail( cvrt$getSlice(cvrt$nSlices()), n = 1));
        }
        if( cvrt$nRows()>0 ) {
            cvrt$SetNanRowMean();
            cvrt$CombineInOneSlice();
            cvrt = rbind(matrix(1,1,snps$nCols()),cvrt$getSlice(1));
        } else {
            cvrt = matrix(1,1,snps$nCols());
        }
        # Correct for the error covariance structure
        if( length(correctionMatrix)>0 ) {
            cvrt = cvrt %*% correctionMatrix;
        }
        # Orthonormalize covariates
        # status("Orthonormalizing covariates");
        q = qr(t(cvrt));
        if( min(abs(diag(qr.R(q)))) < .Machine$double.eps * snps$nCols() ) {
            stop("Colinear or zero covariates detected");
        }
        cvrt = t( qr.Q(q) );
        rm( q );
    }
    ################################# Gene expression processing ############################
    {
        status("Processing gene expression data (imputation, residualization, etc.)");
        # Impute gene expression
        gene$SetNanRowMean();
        # Correct for the error covariance structure
        if( length(correctionMatrix)>0 ) {
            gene$RowMatrixMultiply(correctionMatrix);
        }
        # Orthogonolize expression w.r.t. covariates
        # status("Orthogonolizing expression w.r.t. covariates");
        gene_offsets = double(gene$nSlices()+1);
        for( sl in 1:gene$nSlices() ) {
            slice = gene$getSlice(sl);
            gene_offsets[sl+1] = gene_offsets[sl] + nrow(slice);
            rowsq1 = rowSums(slice^2);
            slice = slice - tcrossprod(slice,cvrt) %*% cvrt;
            rowsq2 = rowSums(slice^2);
            # kill rows colinear with the covariates
            delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
            slice[delete.rows,] = 0;
            rowsq2[delete.rows] = 1;
            div = sqrt(rowsq2); #sqrt( rowSums(slice^2) );
            # 			div[ div == 0 ] = 1;
            gene.std$set(sl, div);
            gene$setSlice(sl, slice / div);
        }
        rm(rowsq1, rowsq2, delete.rows, div);
        rm( sl, slice );
        #gene$RowRemoveZeroEps();
    }
    ################################# snps_process, testfun, pvfun, threshfun, afun  ########
    {
        # snps_process - preprocess SNPs slice
        #
        # afun --- abs for signed stats, identity for non-negative
        # threshfun --- internal stat threshold for given p-value
        # testfun --- t or F statistic from the internal one
        # pvfun --- p-value from the t or F statistic

        nSamples = snps$nCols();
        nGenes = gene$nRows();
        nSnps  = snps$nRows();
        nCov = nrow(cvrt);
        # nVarTested = length(snps_list); # set in case(useModel)
        # dfNull = nSamples - nCov;
        # d.f. of the full model
        betafun = NULL;

        if( useModel == modelLINEAR ) {
            snps_process = function(x) {
                return( list(.SetNanRowMean(x)) );
            };
            nVarTested = 1;
            dfFull = nSamples - nCov - nVarTested;
            statistic.fun = function(mat_list) {
                return( mat_list[[1]] );
            }
            afun = function(x) {return(abs(x))};
            threshfun = function(pv) {
                thr = qt(pv/2, dfFull, lower.tail = FALSE);
                thr = thr^2;
                thr = sqrt(  thr / (dfFull + thr) );
                thr[pv >= 1] = 0;
                thr[pv <= 0] = 1;
                return( thr );
            }
            testfun = function(x) { return( x * sqrt( dfFull / (1 - .my.pmin(x^2,1))));	}
            pvfun = function(x) { return( .pv.nz(pt(-abs(x),dfFull)*2)); }
            thresh.cis = threshfun(pvOutputThreshold.cis);
            thresh = threshfun(pvOutputThreshold);
            betafun = function(stat, ss, gg, select) {
                return(stat * gene.std$get(gg)[select[,1]] / snps.std$get(ss)[select[,2]]);
            }
        } else
            if( useModel == modelANOVA ) {
                snps_process = function(x).SNP_process_split_for_ANOVA(x,n.anova.groups);
                nVarTested = n.anova.groups - 1;
                dfFull = nSamples - nCov - nVarTested;
                # 			statistic.fun = function(mat_list) {
                # 				return( mat_list[[1]]^2 + mat_list[[2]]^2 );
                # 			}
                statistic.fun = function(mat_list) {
                    x = mat_list[[1]]^2;
                    for( j in 2:length(mat_list) )
                        x = x + mat_list[[j]]^2;
                    return( x );
                }
                afun = identity;
                threshfun = function(pv) {
                    thr = qf(pv, nVarTested, dfFull, lower.tail = FALSE);
                    thr = thr / (dfFull/nVarTested + thr);
                    thr[pv >= 1] = 0;
                    thr[pv <= 0] = 1;
                    return( thr );
                }
                testfun = function(x) { return( x / (1 - .my.pmin(x,1)) * (dfFull/nVarTested) ); }
                pvfun = function(x) { return( .pv.nz(pf(x, nVarTested, dfFull, lower.tail = FALSE)) ); }
                thresh.cis = threshfun(pvOutputThreshold.cis);
                thresh = threshfun(pvOutputThreshold);
            } else
                if( useModel == modelLINEAR_CROSS ) {
                    last.covariate = as.vector( last.covariate );
                    snps_process = .SNP_process_split_for_LINEAR_CROSS = function(x) {
                        out = vector("list", 2);
                        out[[1]] = .SetNanRowMean(x);
                        out[[2]] = t( t(out[[1]]) * last.covariate );
                        return( out );
                    };
                    nVarTested = 1;
                    dfFull = nSamples - nCov - nVarTested - 1;
                    statistic.fun = function(mat_list) {
                        return( mat_list[[2]] / sqrt(1 - mat_list[[1]]^2) );
                    }
                    afun = function(x) {return(abs(x))};
                    threshfun = function(pv) {
                        thr = qt(pv/2, dfFull, lower.tail = FALSE);
                        thr = thr^2;
                        thr = sqrt(  thr / (dfFull + thr) );
                        thr[pv >= 1] = 0;
                        thr[pv <= 0] = 1;
                        return( thr );
                    }
                    testfun = function(x) { return( x * sqrt( dfFull / (1 - .my.pmin(x^2,1))));	}
                    pvfun = function(x) { return( .pv.nz(pt(-abs(x),dfFull)*2 )); }
                    thresh.cis = threshfun(pvOutputThreshold.cis);
                    thresh = threshfun(pvOutputThreshold);
                    betafun = function(stat, ss, gg, select) {
                        return(stat * gene.std$get(gg)[select[,1]] / snps.std$get(ss)[select[,2]]);
                    }
                }
        params$dfFull = dfFull;
    }
    ################################# Saver class(es) creation ##############################
    {
        status("Creating output file(s)");
        if(noFDRsaveMemory) {
            if( pvOutputThreshold > 0 ) {
                saver.tra = .OutputSaver_direct$new();
            }
            if( pvOutputThreshold.cis > 0 ) {
                saver.cis = .OutputSaver_direct$new();
            }
        } else {
            if( pvOutputThreshold > 0 ) {
                saver.tra = .OutputSaver_FRD$new();
            }
            if( pvOutputThreshold.cis > 0 ) {
                saver.cis = .OutputSaver_FRD$new();
            }
        }
        if( pvOutputThreshold > 0 )
            if( pvOutputThreshold * gene$nRows() * snps$nRows() > 1000000 )
                if(!noFDRsaveMemory)
                    cat("Warning: pvOutputThreshold may be too large.\nExpected number of findings > ",
                        pvOutputThreshold * gene$nRows() * snps$nRows(),"\n");
        if( (useModel == modelLINEAR) || (useModel == modelLINEAR_CROSS) ) {
            statistic_name = "t-stat";
        } else if( useModel == modelANOVA ) {
            statistic_name = "F-test";
        }
        if(!is.null(betafun))
            statistic_name = paste("beta\t",statistic_name, sep="");
        if( pvOutputThreshold > 0 )
            saver.tra$start("temp.txt",     statistic_name, snps, gene, testfun, pvfun);
        if( pvOutputThreshold.cis > 0 )
            saver.cis$start(output_file_name.cis, statistic_name, snps, gene, testfun, pvfun);
        rm( statistic_name );
    }
    ################################# Some useful functions #################################
    {
        orthonormalize.snps = function(cursnps, ss) {
            for(p in 1:length(cursnps)) {
                if(length(correctionMatrix)>0) {
                    cursnps[[p]] = cursnps[[p]] %*% correctionMatrix;
                }
                rowsq1 = rowSums(cursnps[[p]]^2);
                cursnps[[p]] = cursnps[[p]] - tcrossprod(cursnps[[p]],cvrt) %*% cvrt;
                for(w in .seq(1L,p-1L))
                    cursnps[[p]] = cursnps[[p]] - rowSums(cursnps[[p]]*cursnps[[w]]) * cursnps[[w]];
                rowsq2 = rowSums(cursnps[[p]]^2);
                delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
                cursnps[[p]][delete.rows,] = 0;
                div = sqrt( rowsq2 );
                div[ delete.rows ] = 1;
                # 				show(c(rowsq2,rowsq1, div));
                cursnps[[p]] = cursnps[[p]]/div;
            }
            snps.std$set(ss, div);
            return(cursnps);
        }
        # 		if( pvOutputThreshold.cis > 0 ) {
        # 			is.cis.pair = function(gg,ss) {
        # 				return(!( ( snpsloc[[ss]][1, 1] - tail( geneloc[[gg]][ , 2], n = 1L) > cisDist) |
        # 					    ( geneloc[[gg]][1, 1] - tail( snpsloc[[ss]]      , n = 1L) > cisDist) ) );
        # 			}
        # 		}
        if( pvOutputThreshold.cis > 0 ) {
            # 			sn.l = sapply(snpsloc, function(x)x[1] );
            # 			sn.r = sapply(snpsloc, function(x)tail(x,1) );
            # 			ge.l = sapply(geneloc, function(x)x[1,1] );
            # 			ge.r = sapply(geneloc, function(x)x[nrow(x) , 2] );
            sn.l = sapply(snpsloc, "[", 1 );
            sn.r = sapply(snpsloc, tail, 1 );
            ge.l = sapply(geneloc, "[", 1, 1 );
            ge.r = sapply( lapply(geneloc, tail.matrix, 1 ), "[", 2);
            gg.1 = findInterval( sn.l , ge.r + cisDist +1) + 1;
            # 			cat(gg.1,"\n")
            gg.2 = findInterval( sn.r , ge.l - cisDist );
            # 			cat(gg.2,"\n")
            rm(sn.l, sn.r, ge.l, ge.r);
        }

    }
    ################################# Prepare counters and histogram bins ###################
    {
        pvbins = NULL; # bin edges for p-values
        statbins = 0;  # bin edges for the test statistic (|t| or F)
        do.hist = FALSE;
        if( length(pvalue.hist) == 1 ) {
            if(pvalue.hist == "qqplot") {
                pvbins = c(0, 10^rev(seq(0, log10(.Machine$double.xmin)-1, -0.05)));
            } else
                if( is.numeric(pvalue.hist) ) {
                    pvbins = seq(from = 0, to = 1, length.out = pvalue.hist+1);
                } else
                    if( pvalue.hist == TRUE ) {
                        pvbins = seq(from = 0, to = 1, length.out = 100+1);
                    }
        } else
            if( is.numeric(pvalue.hist) && (length(pvalue.hist) > 1) ) {
                pvbins = pvalue.hist;
            }
        if( is.null(pvbins) && (pvalue.hist != FALSE) ) {
            stop("Wrong value of pvalue.hist. Must be FALSE, TRUE, \"qqplot\", or numerical");
        }
        do.hist = !is.null(pvbins);
        if( do.hist ) {
            pvbins = sort(pvbins);
            statbins = threshfun(pvbins);
            if( pvOutputThreshold > 0) {
                hist.all = .histogrammer$new(pvbins, statbins);
            }
            if( pvOutputThreshold.cis > 0) {
                hist.cis = .histogrammer$new(pvbins, statbins);
            }
        }
        rm( pvbins, statbins);
        if(min.pv.by.genesnp) {
            if( pvOutputThreshold > 0) {
                minpv.tra = .minpvalue$new(snps,gene);
            }
            if( pvOutputThreshold.cis > 0) {
                minpv.cis = .minpvalue$new(snps,gene);
            }
        }
    }
    ################################# Main loop #############################################
    {
        beta = NULL;
        n.tests.all = 0;
        n.tests.cis = 0;
        n.eqtls.tra = 0;
        n.eqtls.cis = 0;

        status("Performing eQTL analysis");
        # ss = 1; gg = 1;
        # ss = snps$nSlices(); gg = gene$nSlices();

        snps_offset = 0;
        for(ss in 1:snps$nSlices()) {
            # 		for(ss in 1:min(2,snps$nSlices())) { #for debug
            cursnps = NULL;
            nrcs = snps$GetNRowsInSlice(ss);

            # loop only through the useful stuff
            for(gg in if(pvOutputThreshold>0){1:gene$nSlices()}else{.seq(gg.1[ss],gg.2[ss])} ) {
                gene_offset = gene_offsets[gg];
                curgene = gene$getSlice(gg);
                nrcg = nrow(curgene);
                if(nrcg == 0) next;

                rp = "";

                statistic = NULL;
                select.cis.raw = NULL;
                ## do cis analysis
                # 				if( (pvOutputThreshold.cis > 0) && ( is.cis.pair(gg, ss) ) ) {
                if( (pvOutputThreshold.cis > 0) && (gg >= gg.1[ss]) && (gg <= gg.2[ss]) ) {

                    if( is.null( statistic ) ) {
                        if( is.null( cursnps ) ) {
                            cursnps = orthonormalize.snps( snps_process( snps$getSlice(ss) ), ss );
                        }
                        mat = vector("list", length(cursnps));
                        for(d in 1:length(cursnps)) {
                            mat[[d]] = tcrossprod(curgene, cursnps[[d]]);
                        }
                        statistic = statistic.fun( mat );
                        astatistic = afun(statistic);
                        # 						rm(mat);
                    }

                    # 					sn.l = findInterval(geneloc[[gg]][ ,1] - cisDist-1  +1   , snpsloc[[ss]]);
                    # 					sn.r = findInterval(geneloc[[gg]][ ,2] + cisDist    -1   , snpsloc[[ss]]);
                    sn.l = findInterval(geneloc[[gg]][ ,1] - cisDist-1, snpsloc[[ss]]);
                    sn.r = findInterval(geneloc[[gg]][ ,2] + cisDist, snpsloc[[ss]]);
                    xx = unlist(lapply(which(sn.r>sn.l),FUN=function(x){(sn.l[x]:(sn.r[x]-1))*nrow(statistic)+x}))
                    select.cis.raw = xx[ astatistic[xx] >= thresh.cis ];
                    select.cis = arrayInd(select.cis.raw, dim(statistic))

                    n.tests.cis = n.tests.cis + length(xx);
                    n.eqtls.cis = n.eqtls.cis + length(select.cis.raw);

                    if( do.hist )
                        hist.cis$update(astatistic[xx]);

                    if( min.pv.by.genesnp ) {
                        # 					minpv.cis$updatecis(ss, gg, arrayInd(xx, dim(statistic)), astatistic[xx])
                        temp = double(length(astatistic));
                        dim(temp) = dim(astatistic);
                        temp[xx] = astatistic[xx];
                        minpv.cis$update(ss, gg, temp);
                    }

                    if(!is.null(betafun))
                        beta = betafun(mat[[length(mat)]][select.cis.raw], ss, gg, select.cis);

                    saver.cis$update( snps_offset + select.cis[ , 2],
                                      gene_offset + select.cis[ , 1],
                                      statistic[select.cis.raw],
                                      beta);

                    # 				statistic.select.cis  = statistic[ select.cis ];
                    # 				test = testfun( statistic.select.cis );
                    # 				pv = pvfun(test);
                    # 				Saver.cis$WriteBlock( cbind(snps_offset + select.cis[ , 2], gene_offset + select.cis[ , 1], test, pv) );
                    # 				counter.cis$Update(gg, ss, select.cis, pv, n.tests = length(xx), if(do.hist) afun(statistic[xx]) )
                    rp = paste(rp, ", ", formatC(n.eqtls.cis, big.mark=",", format = "f", digits = 0), " cis-eQTLs", sep = "");
                }
                ## do trans/all analysis
                if(pvOutputThreshold>0) {
                    if( is.null( statistic ) ) {
                        if( is.null( cursnps ) ) {
                            cursnps = orthonormalize.snps( snps_process( snps$getSlice(ss) ), ss );
                        }
                        mat = vector("list", length(cursnps));
                        for(d in 1:length(cursnps)) {
                            mat[[d]] = tcrossprod(curgene, cursnps[[d]]);
                        }
                        statistic = statistic.fun( mat );
                        astatistic = afun(statistic);
                        # 						rm(mat);
                    }

                    if( do.hist )
                        hist.all$update(astatistic);

                    if(!is.null(select.cis.raw))
                        astatistic[xx] = -1;
                    # 					select.tra.raw = select.tra.raw[!(select.tra.raw %in% select.cis.raw)];

                    select.tra.raw = which( astatistic >= thresh);
                    select.tra = arrayInd(select.tra.raw, dim(statistic))

                    n.eqtls.tra = n.eqtls.tra + length(select.tra.raw);
                    n.tests.all = n.tests.all + length(statistic);

                    if(!is.null(betafun))
                        beta = betafun(mat[[length(mat)]][select.tra.raw], ss, gg, select.tra);

                    saver.tra$update( snps_offset + select.tra[ , 2],
                                      gene_offset + select.tra[ , 1],
                                      statistic[select.tra.raw],
                                      beta);

                    if( min.pv.by.genesnp )
                        minpv.tra$update(ss, gg, astatistic)

                    # 				statistic.select.tra = statistic[ select.tra ];
                    # 				test = testfun( statistic.select.tra );
                    # 				pv = pvfun( test );
                    # 				Saver$WriteBlock( cbind( snps_offset + select.tra[ , 2], gene_offset + select.tra[ , 1], test, pv) );
                    # 				counter$Update(gg, ss, select.tra, pv, n.tests = nrcs*nrcg, if(do.hist) afun(statistic) )
                    rp = paste(rp, ", ", formatC(n.eqtls.tra, big.mark=",", format = "f", digits = 0), " ", addInfo, sep = "")
                }

                #gene_offset = gene_offset + nrcg;
                if( !is.null(statistic) ) {
                    per = 100*(gg/gene$nSlices() + ss-1) / snps$nSlices();
                    cat( formatC(floor(per*100)/100, format = "f", width = 5, digits = 2), "% done" , rp, "\n", sep = "");
                    flush.console();
                }
            } # gg in 1:gene$nSlices()
            snps_offset = snps_offset + nrcs;
        } # ss in 1:snps$nSlices()
    }
    ################################# Results collection ####################################
    {
        rez = list(time.in.sec = proc.time()[3] - start.time);
        rez$param = params;

        if(pvOutputThreshold.cis > 0) {
            rez.cis = list(ntests = n.tests.cis, neqtls = n.eqtls.cis);
            rez.cis = c(rez.cis, saver.cis$getResults( gene, snps, n.tests.cis) );
            if(do.hist)
                rez.cis = c(rez.cis, hist.cis$getResults() );
            if(min.pv.by.genesnp)
                rez.cis = c(rez.cis, minpv.cis$getResults(snps, gene, pvfun = function(x){pvfun(testfun(x))}) );
        }

        if(pvOutputThreshold>0) {
            rez.all = list(ntests = n.tests.all, neqtls = n.eqtls.tra + n.eqtls.cis);
            if(pvOutputThreshold.cis > 0) {
                rez.tra = list(ntests = n.tests.all - n.tests.cis, neqtls = n.eqtls.tra);
                rez.tra = c(rez.tra, saver.tra$getResults( gene, snps, n.tests.all - n.tests.cis) );
            } else {
                rez.all = c(rez.all, saver.tra$getResults( gene, snps, n.tests.all              ) );
            }
            if(do.hist) {
                rez.all = c(rez.all, hist.all$getResults() );
                if(pvOutputThreshold.cis > 0) {
                    rez.tra$hist.bins = rez.all$hist.bins;
                    rez.tra$hist.counts = rez.all$hist.counts - rez.cis$hist.counts;
                }
            }
            if(min.pv.by.genesnp) {
                if(pvOutputThreshold.cis > 0) {
                    rez.tra = c(rez.tra, minpv.tra$getResults(snps, gene, pvfun = function(x){pvfun(testfun(x))}) );
                } else {
                    rez.all = c(rez.all, minpv.tra$getResults(snps, gene, pvfun = function(x){pvfun(testfun(x))}) );
                }
            }
        }

        if(exists("rez.all")>0)
            rez$all = rez.all;
        if(exists("rez.tra")>0)
            rez$trans = rez.tra;
        if(exists("rez.cis")>0)
            rez$cis = rez.cis;

        class(rez) = c(class(rez),"MatrixEQTL");
        status("");
    }
    # 	cat("s std ",snps.std$get(1),"\n");
    # 	cat("g std ",gene.std$get(1),"\n");
    ################################# Results collection ####################################
    return(rez);
}


.listBuilder_NEW <- setRefClass(".listBuilder_NEW",
                            fields = list(
                                dataEnv = "environment",
                                n = "integer"
                            ),
                            methods = list(
                                initialize = function() {
                                    dataEnv <<- new.env(hash = TRUE);
                                    n <<- 0L;
                                    # 			cumlength <<- 0;
                                    return(.self);
                                },
                                add = function(x) {
                                    if(length(x) > 0) {
                                        n <<- n + 1L;
                                        # 				cumlength <<- cumlength + length(x);
                                        assign(paste(n), x, dataEnv );
                                    }
                                    return(.self);
                                },
                                set = function(i,x) {
                                    i = as.integer(i);
                                    if(length(x) > 0) {
                                        if(i>n)
                                            n <<- i;
                                        assign(paste(i), x, dataEnv );
                                    }
                                    return(.self);
                                },
                                get = function(i) {
                                    return(base::get(paste(i),dataEnv));
                                },
                                list = function() {
                                    if(n==0)	return(list());
                                    result = vector("list",n);
                                    for( i in 1:n) {
                                        result[[i]] = .self$get(i);
                                    }
                                    return(result);
                                },
                                unlist = function() {
                                    return(base::unlist(.self$list(), recursive=FALSE, use.names = FALSE));
                                },
                                show = function() {
                                    #cat(".listBuilder object.\nIternal object in MatrixEQTL package.\n");
                                    #cat("Number of elements:", .self$n, "\n");
                                }
                            ))


.OutputSaver_NEW <- setRefClass(".OutputSaver_NEW",
                                fields = list(
                                    sdata = ".listBuilder_NEW",
                                    gdata = ".listBuilder_NEW",
                                    cdata = ".listBuilder_NEW",
                                    bdata = ".listBuilder_NEW",
                                    fid = "list",
                                    testfun1 = "list",
                                    pvfun1 = "list"
                                ),
                                methods = list(
                                    initialize = function () {
                                        sdata <<- .listBuilder_NEW$new();
                                        gdata <<- .listBuilder_NEW$new();
                                        cdata <<- .listBuilder_NEW$new();
                                        bdata <<- .listBuilder_NEW$new();
                                        fid <<- list(0);
                                        testfun1 <<- list(0);
                                        pvfun1 <<- list(0);
                                        return(.self);
                                    },
                                    start = function(filename, statistic_name, unused1, unused2, testfun, pvfun) {
                                        testfun1 <<- list(testfun);
                                        pvfun1 <<- list(pvfun);
                                        if(length(filename) > 0) {
                                            if(class(filename) == "character") {
                                                fid <<- list(file(description = filename, open = "wt", blocking = FALSE, raw = FALSE), TRUE);
                                            } else {
                                                fid <<- list(filename, FALSE)
                                            }
                                            writeLines( paste("SNP\tgene\t",statistic_name,"\tp-value\tFDR", sep = ""), fid[[1]]);
                                        } else {
                                            fid <<- list();
                                        }
                                    },
                                    update = function(spos, gpos, sta, beta = NULL) {
                                        if(length(sta)>0) {
                                            sdata$add(spos);
                                            gdata$add(gpos);
                                            cdata$add(sta );
                                            if(!is.null(beta ))
                                                bdata$add(beta );
                                        }
                                        return(.self);
                                    },
                                    getResults = function( gene, snps, FDR_total_count) {
                                        pvalues = NULL;
                                        if(cdata$n > 0) {
                                            tests = testfun1[[1]](cdata$unlist());
                                            cdata <<- .listBuilder_NEW$new();

                                            pvalues = pvfun1[[1]](tests);
                                            ord = order(pvalues);

                                            tests = tests[ord];
                                            pvalues = pvalues[ord];

                                            snps_names  = rownames(snps)[sdata$unlist()[ord]];
                                            sdata <<- .listBuilder_NEW$new();
                                            gene_names  = rownames(gene)[gdata$unlist()[ord]];
                                            gdata <<- .listBuilder_NEW$new();

                                            beta = NULL;
                                            if(bdata$n > 0)
                                                beta = bdata$unlist()[ord];

                                        } else {
                                            cat("No significant associations were found.\n", file = if(length(fid)>0){fid[[1]]}else{""});
                                        }
                                        if(length(fid)>0)	{
                                            if(fid[[2]]) {
                                                close(fid[[1]]);
                                            }
                                            fid <<- list();
                                        }

                                        if(!is.null(pvalues)) {
                                            eqtls = list( snps = snps_names,
                                                          gene = gene_names,
                                                          statistic = tests,
                                                          pvalue = pvalues);
                                            if(!is.null(beta))
                                                eqtls$beta = beta;
                                        } else {
                                            eqtls = list( snps = character(),
                                                          gene = character(),
                                                          beta = numeric(),
                                                          statistic = numeric(),
                                                          pvalue = numeric());
                                        }
                                        return(list(eqtls = data.frame(eqtls)));
                                    }
                                )
)