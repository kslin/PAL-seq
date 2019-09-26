library(tidyverse)
source('/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R')
args <- commandArgs(trailingOnly = TRUE)

StdMedians <- read_tsv(paste0(args[1],'/standard_medians.txt'))
ModelParams <- read_tsv(paste0(args[1],'/model_params.txt'))


StdMediansGat <- gather(StdMedians, key = concentration, value = intensity, -nucleotide_length)

BestFit <- as_tibble(apply(ModelParams[,c(1,2)],1,function(x){x[1]*(0:500) + x[2]}),.name_repair='unique')
colnames(BestFit) = colnames(StdMedians[,-1])

BestFit$nucleotide_length <- 0:500
BestFitGat <- gather(BestFit, key = concentration, value = intensity, -nucleotide_length)

##Plotting
pFit <- ggplot(StdMediansGat, aes(x = nucleotide_length, y = intensity, label = nucleotide_length)) + 
	geom_point() + 
	geom_line(data = BestFitGat,  aes(x = nucleotide_length, y = intensity)) +
	facet_wrap(~concentration, ncol = 5) + 
	geom_text(position=position_jitter(width=10,height=10)) + 
	theme_tim_label()

ggsave(plot = pFit, file = "pFit.pdf", width = 6, height = 2)

