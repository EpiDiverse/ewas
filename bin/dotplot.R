#!/usr/bin/env R

library(ggplot2)


file = read.table(file.choose(),header=F)
colnames(file)= c("cpg","cpg_pos","snp", "snp_pos","beta", "stats", "pvalue", "FDR")

test= read.table(file.choose(),header=F)
colnames(test)= c("cpg","cpg_pos","snp", "snp_pos","beta", "stats", "pvalue", "FDR")


df <- data.frame(file)
df_test <- data.frame(test)


min_test_cg=min(df_test$cpg_pos)
max_test_cg=max(df_test$cpg_pos)

min_test_snp=min(df_test$snp_pos)
max_test_snp =max(df_test$snp_pos)

min_cpg=min(df$cpg_pos)
max_cpg=max(df$cpg_pos)

min_snp=min(df$snp_pos)
max_snp =max(df$snp_pos)


if(min_cpg>min_snp) {
   min=min_snp
} else {
   min= min_cpg
} 

if(max_cpg>max_snp) {
   max=max_cpg
} else {
   max= max_cpg
} 

if(min_test_cg>min_test_snp) {
   min=min_test_snp
} else {
   min= min_test_cg
} 

if(max_test_cg>max_test_snp) {
   max=max_test_cg
} else {
   max= max_test_snp
} 



#I will get the average of start and start+1 positions to attain them to single points
#bash
# awk '{OFS="\t"}{$5=sprintf("%.6f",$6); $6=sprintf("%.6f",$6)}1' small.txt |  awk -F"-" '$1=$1' | awk -F":" '$1=$1' | awk '{print $1":"$2"-"$3"\t"($2+$3)/2"\t"$4":"$5"-"$6"\t"($5+$6)/2"\t"$7"\t"$8"\t"$9"\t"$10}' |  tail -n+2 > small2.txt
#convert e values to decimal since I will use "-" delimiter and number-evalue should be avoided
#then print the old files with extra columns with average info next to them

#factoring the snp values based on snp_pos column
#df$new_snp <- factor(df$snp, levels=unique(df$snp[order(df$snp_pos)]), ordered=TRUE)


pdf("scatter.pdf")
 p =ggplot(data = df) +
  geom_point(aes(snp_pos,cpg_pos, color=ifelse(abs((df$snp_pos)-(df$cpg_pos)) <= 2000,'blue','red'), size = 0.0001), show.legend = FALSE)  + theme(axis.text.y=element_text(angle=0, vjust=1),
           text=element_text(size=5),axis.text.x = element_text(angle = 90, hjust=1)) +
   ggtitle("Scatter plot of  SNPs vs CpGs")  + xlab("SNP") + ylab("CpG") 
 
 p+ xlim(min,max) + ylim(min,max) + theme(
       plot.title = element_text(color="black", size=12, face="bold.italic"),
       axis.title.x = element_text(color="black", size=12, face="bold"),
       axis.title.y = element_text(color="black", size=12, face="bold") ) +
    theme(axis.line = element_line(color="black", size = 0.5))  + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.grid.major = element_line(colour = "black"))
dev.off()


#col = ifelse(abs((df$snp_pos)-(df$cpg_pos)) < 2000,'blue','red')
#102718   102668     
#102718   102710  
#132138   103610  
pdf("sig_scatter.pdf")
p =ggplot(data = df_test) +
   geom_point(aes(snp_pos,cpg_pos, color=ifelse(abs((df_test$snp_pos)-(df_test$cpg_pos)) <= 2000,'blue','red')), show.legend = FALSE, size=1, shape=".")  + theme(axis.text.y=element_text(angle=0, vjust=1),
                                                                                                                                                    text=element_text(size=5),axis.text.x = element_text(angle = 90, hjust=1)) +
   ggtitle("Scatter plot of  SNPs vs CpGs")  + xlab("SNP") + ylab("CpG") 

p+ xlim(min,max) + ylim(min,max) + theme(
   plot.title = element_text(color="black", size=12, face="bold.italic"),
   axis.title.x = element_text(color="black", size=12, face="bold"),
   axis.title.y = element_text(color="black", size=12, face="bold") ) +
   theme(axis.line = element_line(color="black", size = 0.5))  + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.grid.major = element_line(colour = "black"))
dev.off()

geom_vline(data = df_test, aes(xintercept = as.numeric(time)))
