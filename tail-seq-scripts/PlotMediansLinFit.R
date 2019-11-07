library(tidyverse)
source('/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R')
args <- commandArgs(trailingOnly = TRUE)

StdMedians <- read_tsv(paste0(args[1],'/standard_medians.txt'))
ModelParams <- read_tsv(paste0(args[1],'/model_params.txt'))
AllTagsStds <- read_tsv(paste0(args[1],'/standard_tail_lengths.txt'))
AllTags <- read_tsv(paste0(args[1],'/mapped_tail_lengths_stds_removed.txt'))


StdMediansGat <- gather(StdMedians, key = concentration, value = intensity, -nucleotide_length)

BestFit <- as_tibble(apply(ModelParams[,c(1,2)],1,function(x){x[1]*(0:350) + x[2]}),.name_repair='unique')
colnames(BestFit) = colnames(StdMedians[,-1])

BestFit$nucleotide_length <- 0:350
BestFitGat <- gather(BestFit, key = concentration, value = intensity, -nucleotide_length)

#pearson data frame
ModelParamsPlot = tibble(
	concentration = colnames(StdMedians[,-1]),
	pearson = paste0('R = ',round(ModelParams$r_value,3)))
print(ModelParamsPlot)
##Plotting
pFit <- ggplot(data = StdMediansGat) + 
	geom_line(data = BestFitGat,  aes(x = nucleotide_length + 12, y = intensity), linetype = 'dashed', color = 'grey') +
	geom_point(data = StdMediansGat, aes(x = nucleotide_length + 12, y = intensity), size = 0.5) + 
	geom_text(data = StdMediansGat, aes(x = nucleotide_length + 12, y = intensity, label = nucleotide_length + 12),nudge_x = 0, nudge_y = 0.5, size = 1) +
	geom_text(data = ModelParamsPlot, aes(x = 40, y = 1.5, label = pearson), size = 2) + 
	facet_wrap(~concentration, ncol = 1) + 
	theme_tim_label()

pStdCDF <- ggplot(data = AllTagsStds, aes(x = tail_length, color = as.factor(accession))) + 
	stat_ecdf() + 
	geom_vline(xintercept = unique(as.numeric(as.character(AllTagsStds$accession))) + 12, linetype = 'dashed',color = 'grey') +
	theme_tim_label()

AllTags[AllTags$tail_length > 300,]$tail_length = 300
pReadHist <- ggplot(data = AllTags, aes(x = tail_length)) + 
	geom_histogram(binwidth = 1) + 
	scale_x_continuous(expand = c(0,0), breaks = c(0:6)*50, labels = c(0:6)*50) +
	scale_y_continuous(expand = c(0,0)) +
	theme_tim_label()

ggsave(plot = pFit, file = paste0(args[1],'/BestFitModels.pdf'), width = 2, height = 8)
ggsave(plot = pStdCDF, file = paste0(args[1],'/AllStandardTagsCDF.pdf'), width = 5, height = 5)
ggsave(plot = pReadHist, file = paste0(args[1],'/AllTagsHist.pdf'), width = 5, height = 5)

