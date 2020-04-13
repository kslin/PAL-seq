library(tidyverse)
source('/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R')
args <- commandArgs(trailingOnly = TRUE)
# args = '/lab/solexa_bartel/teisen/Tail-seq/PALseqV3/Run20200213/data/total'
StdMedians <- read_tsv(paste0(args[1],'/standard_medians.txt'))
ModelParams <- read_tsv(paste0(args[1],'/model_params.txt'))
AllTagsStds <- read_tsv(paste0(args[1],'/standard_tail_lengths.txt'))
AllTags <- read_tsv(paste0(args[1],'/mapped_tail_lengths_stds_removed.txt'))


StdMediansGat <- gather(StdMedians, key = concentration, value = intensity, -nucleotide_length)

BestFit <- as_tibble(
	apply(ModelParams[,c(1,2)],1,function(x){
			x[1]*seq(0,1.5,length.out = 1000) + x[2]
		}
	),.name_repair='unique')
colnames(BestFit) = colnames(StdMedians[,-1])
BestFit$intensity <- seq(0,1.5,length.out = 1000)
BestFitGat <- gather(BestFit, key = concentration, value = nucleotide_length, -intensity)

#pearson data frame
ModelParamsPlot = tibble(
	concentration = colnames(StdMedians[,-1]),
	pearson = paste0('R = ',round(ModelParams$r_value,3)))
print(ModelParamsPlot)
##Plotting
pFit <- ggplot(data = filter(StdMediansGat, concentration == 'conc_1')) + 
	geom_line(data = filter(BestFitGat, concentration == 'conc_1'),  aes(y = nucleotide_length + 11, x = intensity), linetype = 'dashed', color = 'grey') +
	geom_point(data = filter(StdMediansGat, concentration == 'conc_1'), aes(y = nucleotide_length + 11, x = intensity), size = 0.5) + 
	# geom_text(data = StdMediansGat, aes(y = nucleotide_length + 11, x = intensity, label = nucleotide_length + 11),nudge_y = 0, nudge_x = 0.5, size = 1) +
	# geom_text(data = ModelParamsPlot, aes(y = 40, x = 1.5, label = pearson), size = 2) + 
	facet_wrap(~concentration, ncol = 1) + 
	scale_x_continuous(expand = c(0,0), name = 'Normalized intensity') +
	scale_y_continuous(expand = c(0,0), name = 'Calculated tail length') +
	theme_tim_label()
# break
pStdCDF <- ggplot(data = AllTagsStds, aes(x = tail_length, color = as.factor(accession))) + 
	stat_ecdf() + 
	geom_vline(xintercept = unique(as.numeric(as.character(AllTagsStds$accession))), linetype = 'dashed',color = 'grey') +
	scale_x_continuous(name = 'Tail length, nt', expand = c(0,0)) +
	scale_y_continuous(name = 'Cumulative fraction', expand = c(0,0)) +
	labs(color = 'Standard') + 
	theme_tim_label()

AllTags[AllTags$tail_length > 300,]$tail_length = 300

pReadHist <- ggplot(data = AllTags, aes(x = tail_length)) + 
	geom_histogram(binwidth = 1) + 
	scale_x_continuous(name = 'Tail length, nt', expand = c(0,0), breaks = c(0:6)*50, labels = c(0:6)*50) +
	scale_y_continuous(name = 'Read counts', expand = c(0,0)) +
	theme_tim_label()

# ggsave(plot = pFit, file = paste0(args[1],'/BestFitModels.pdf'), width = 2, height = 8)
ggsave(plot = pFit, file = paste0(args[1],'/BestFitModels.pdf'), width = 4, height = 4)

ggsave(plot = pStdCDF, file = paste0(args[1],'/AllStandardTagsCDF.pdf'), width = 5, height = 5)
ggsave(plot = pReadHist, file = paste0(args[1],'/AllTagsHist.pdf'), width = 5, height = 5)

