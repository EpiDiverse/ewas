#!/usr/bin/env R

#library(ggplot2)

args <- commandArgs(trailingOnly=T)

# args[1]   path to input data
# args[2]   path to output

df <- read.table(args[1], header=T)
df$pvalue <- as.numeric(df$pvalue)
df$log_pvalue <- -log10(df$pvalue)

n <- length(df$log_pvalue)    # number of observations
r <- 1:n                      # order of values, i.e. ranks without averaged ties
p <- (r - 1/2) / n            # assign to ranks using Blom's method
x <- qnorm(p)                 # theoretical standard normal quantiles for p values

title <- paste("Q-Q plot for all",nrow(df),"p-values")
qq.name <- paste0(args[2], ".qqplot.png")

png(qq.name, width=12, height=7)
plot(x, df$log_pvalue, xlab= "-Log10 (p-value), theoretical", ylab= "-Log10 (p-value), observed", main= title)
qqline(df$log_pvalue,col = "red")
dev.off()

hg <- ggplot(df, aes(log_pvalue)) +
  geom_histogram(binwidth = 0.05)  +  
  ggtitle("Histogram of -Log10 (p-values)") +
  xlab("-Log10 (p-value)") + ylab("Frequency") +
  theme_bw() +
  theme(plot.title = element_text(color="black", size=14, face="bold.italic",hjust=0.5))

hg.name <- paste0(args[2], ".pval_hist.png")
ggsave(hg.name, hg, width=12, height=7)
