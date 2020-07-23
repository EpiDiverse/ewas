make_GO_tests=function(data,file_tree_name_over,file_name_over,file_tree_name_under,file_name_under,ontology="BP",pvalueCutoff=0.01){
        gene<-as.character(data[,1])
        paramsOver <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",geneSetCollection=gsc,geneIds = gene,universeGeneIds = universe,ontology=ontology,pvalueCutoff = pvalueCutoff,conditional = FALSE,testDirection = "over")
        Over<-hyperGTest(paramsOver)
        SummaryOver=summary(Over)
        if(!is.na(SummaryOver[1,1])){
                SummaryOver[1:7,"q_value"]<-"NA"
                SummaryOver$q_value<-p.adjust(SummaryOver[,2],method="BH")
                write.table(SummaryOver, file=file_name_over,append=FALSE,sep="\t",col.names=TRUE)
                tgsOver = termGraphs(Over,use.terms=TRUE , pvalue=.01)
                for (i in names(tgsOver)){
                        g1 = tgsOver[[i]]
                        #bigg = inducedTermGraph(Over, id=nodes(g1))
                        bigg = inducedTermGraph(Over, id=nodes(g1),children = FALSE, parents = TRUE)
                        nodeDataDefaults(bigg) <- c(nodeDataDefaults(bigg),list(term=as.character(NA)))
                        theTerms = as.character(lapply(mget(nodes(bigg), GOTERM), Term))
                        nodeData(bigg, attr="term") = theTerms
                        write(rev(theTerms), file=file_tree_name_over,append=TRUE,sep=" ")
                        write("\t\t",file=file_tree_name_over,append=TRUE,sep=" ")
                        pdf(sprintf('%s_%s_tree.pdf',file_tree_name_over,i))
                        plotGOTermGraph(bigg, Over, node.shape="ellipse",add.counts=TRUE,node.colors=c(sig="green", not="yellow"),max.nchar = 100)
                        dev.off()
                }
        }
        paramsUnder <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",geneSetCollection=gsc,geneIds = gene,universeGeneIds = universe,ontology=ontology,pvalueCutoff = pvalueCutoff,conditional = FALSE,testDirection = "under")
        Under<-hyperGTest(paramsUnder)
        SummaryUnder=summary(Under)
        if(!is.na(SummaryUnder[1,1])){
                SummaryUnder[1:7,"q_value"]<-"NA"
                SummaryUnder$q_value<-p.adjust(SummaryUnder[,2],method="BH")
                write.table(SummaryUnder, file=file_name_under,append=FALSE,sep="\t",col.names=TRUE)
                tgsUnder = termGraphs(Under,use.terms=TRUE , pvalue=.01)
                for (j in names(tgsUnder)){
                        g1 = tgsUnder[[j]]
                        #bigg = inducedTermGraph(Over, id=nodes(g1))
                        biggUnder = inducedTermGraph(Under, id=nodes(g1),children = FALSE, parents = TRUE)
                        nodeDataDefaults(biggUnder) <- c(nodeDataDefaults(biggUnder),list(term=as.character(NA)))
                        theTerms = as.character(lapply(mget(nodes(biggUnder), GOTERM), Term))
                        nodeData(biggUnder, attr="term") = theTerms
                        write(rev(theTerms), file=file_tree_name_under,append=TRUE,sep=" ")
                        write("\t\t",file=file_tree_name_under,append=TRUE,sep=" ")
                        pdf(sprintf('%s_%s_tree.pdf',file_tree_name_under,j))
                        plotGOTermGraph(biggUnder, Under, node.shape="ellipse",add.counts=TRUE,node.colors=c(sig="red", not="yellow"),max.nchar = 100)
                        dev.off()
                }
        }
}
