#!/usr/bin/env R

suppressMessages(require(readr))
suppressMessages(require(ggrepel))
suppressMessages(require(ggplot2))
suppressMessages(require(dplyr))
suppressMessages(require(RColorBrewer))

args <- commandArgs(trailingOnly=T)

gg.manhattan <- function(df, threshold, hlight, col, ylims, title){

	# format df
	df.tmp <- df %>% 
		
		# Compute chromosome size
		group_by(CHR) %>% 
		summarise(chr_len=max(BP)) %>% 
		
		# Calculate cumulative position of each chromosome
		mutate(tot=cumsum(chr_len)-chr_len) %>%
		select(-chr_len) %>%
		
		# Add this info to the initial dataset
		left_join(df, ., by=c("CHR"="CHR")) %>%
		
		# Add a cumulative position of each SNP
		arrange(CHR, BP) %>%
		mutate( BPcum=BP+tot) %>%
		
		# Add highlight and annotation information
		mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
		mutate( is_annotate=ifelse(P < threshold, "yes", "no"))

	
	
	# get chromosome center positions for x-axis
	axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
	
	ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
		# Show all points
		geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
		scale_color_manual(values = rep(col, 22 )) +
		
		# custom X axis:
		scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
		scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
		
		# add plot and axis titles
		ggtitle(paste0(title)) +
		
		labs(x = "Scaffold") +
		
		# add genome-wide sig and sugg lines
		geom_hline(yintercept = -log10(sig), color="blue") +
		geom_hline(yintercept = -log10(sugg), linetype="dashed", color="red") +
		
		# Add highlighted points
		#geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
		
		# Add label using ggrepel to avoid overlapping
		geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
		
		# Custom the theme:
		theme_bw(base_size = 15) +
		theme( 
		plot.title = element_text(hjust = 0.5),
		legend.position="none",
		panel.border = element_blank(),
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
				axis.title.x=element_blank(),
				axis.text.x=element_blank(),
				axis.ticks.x=element_blank()
		)
}

sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line

mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B",
               "#FF817E", "#E9534F", "#D92B26", "#AE1612", "#870300",
               "#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B",
               "#4EAFAF", "#2C9696", "#0F8F8F", "#057272", "#005A5A",
               "#5D82BB", "#3B64A5", "#1E4F9E", "#103B7E", "#082B64",
               "#6E65C2", "#4D43AE", "#3226A6", "#211785", "#150D69",
               "#DDF67A", "#C2E04C", "#ADD024", "#88A711", "#678100",
               "#A7E672", "#85D047", "#6AC122", "#4F9B10", "#367800",
               "#A6C965", "#7B9F34", "#567714", "#375201", "#203000")

gmodel= args[1]
gg.manhattan(gmodel, threshold = 1e-6, hlight = NA, ylims=c(0,10),col=mypalette, title="Manhattan Plot")
