#!/usr/bin/env R

library(ggplot2)

args <- commandArgs(trailingOnly=T)

# args[1]   path to input data
# args[2]   path to output

df <- read.table(args[1], header=T)
df$log_pvalue <- -log10(df$pvalue)

title <- paste("Q-Q plot for all",nrow(df),"p-values")
qq <- ggplot(data = df, mapping = aes(sample = log_pvalue)) +
  stat_qq_line(col="red") +
  stat_qq() + labs(x = "-Log10 (p-value), theoretical", y = "-Log10 (p-value), observed") +
  theme_bw() +
  ggtitle(title) + theme(plot.title = element_text(color="black", size=14, face="bold.italic",hjust=0.5))

qq.name <- paste0(args[2], ".qqplot.png")
ggsave(qq.name, qq, width=12, height=7)

ph <- ggplot(df, aes(log_pvalue)) +
  geom_histogram(binwidth = 0.05)  +  
  ggtitle("Histogram of -Log10 (p-values)") +
  xlab("-Log10 (p-value)") + ylab("Frequency") +
  theme_bw() +
  theme(plot.title = element_text(color="black", size=14, face="bold.italic",hjust=0.5))

ph.name <- paste0(args[2], ".pval_hist.png")
ggsave(ph.name, ph, width=12, height=7)
