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
	geom_point(size = 0.5) + 
	geom_line(data = BestFitGat,  aes(x = nucleotide_length, y = intensity)) +
	facet_wrap(~concentration, ncol = 5) + 
	geom_text(nudge_x = 5, nudge_y = 0.01, size = 2) + 
	theme_tim_label()

ggsave(plot = pFit, file = paste0(args[1],'/BestFitModels.pdf'), width = 6, height = 2)

