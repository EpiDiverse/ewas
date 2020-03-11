#!/usr/bin/env R

suppressMessages(require(readr))
suppressMessages(require(ggrepel))
suppressMessages(require(ggplot2))
suppressMessages(require(dplyr))
suppressMessages(require(RColorBrewer))
suppressMessages(require(plotly))

args <- commandArgs(trailingOnly=T)

# args[1]   path to input data
# args[2]   path to output data
# args[3]   significance threshold

gg.manhattan <- function(df, threshold, hlight, col, ylims, title){

	# format df
	df.tmp <- df %>% 
		
		# Compute chromosome size
		group_by(CHR) %>% 
		summarise(chr_len=max(BP)) %>% 
		arrange(desc(chr_len)) %>%
		
		# Calculate cumulative position of each chromosome
		mutate(tot=cumsum(chr_len)-chr_len) %>%
		select(-chr_len) %>%
		
		# Add this info to the initial dataset
		left_join(df, ., by=c("CHR"="CHR")) %>%
		
		# Add a cumulative position of each SNP
		arrange(CHR, BP) %>%
		mutate( BPcum=BP+tot) %>%
		arrange(BPcum) %>%
		
		# Add highlight and annotation information
		mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
		mutate( is_annotate=ifelse(P < threshold, "yes", "no")) %>%
		mutate( P=ifelse(P < 1e-10, 1e-10, P))
	
	# get chromosome center positions for x-axis
	axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
	
	chr.ord <- unique(df.tmp$CHR)
	chr.len <- length(chr.ord)

	p <- ggplot(df.tmp, aes(x=BPcum, y=-log10(P), text=as.factor(SNP))) +
		# Show all points
		geom_point(aes(color=factor(CHR, levels=chr.ord)), alpha=0.8, size=1) +
		scale_color_manual(values = rep(col, chr.len)) +
		
		# custom X axis:
		scale_x_continuous( label = levels(axisdf$CHR), breaks= axisdf$center ) +
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
		geom_text_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=2, force=5) +
		
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

	# save plot image
	p.png <- paste0(title, ".png")
	ggsave(p.png, p, width=12, height=7)

	# save interactive plot
	p <- ggplotly(p,tooltip=c("y","text"))
	p.html <- paste0(title, ".html")
	p.zip <- paste0(title, ".zip")
	htmlwidgets::saveWidget(p, p.html, selfcontained=F, libdir=title, title=title)
	zip::zipr(p.zip, c(p.html,title))

}

sig <- as.numeric(args[3]) # 5e-8 # significant threshold line
sugg <- 1e-3 # 1e-6 # suggestive threshold line

mypalette <- c("#FF817E", "#E9534F", "#D92B26", "#AE1612", "#870300")

emodel <- read.table(args[1],header=T)
if(nrow(emodel) > 0){
	gg.manhattan(emodel, threshold= sig, hlight= NA, ylims=c(0,12), col=mypalette, title=args[2])
}else{
	write(paste0(args[2], " resulted in zero rows"), stderr())
}

