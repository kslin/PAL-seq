.DEFAULT_GOAL := help

#####################################################
#  Makefile for PALseq V3 pipeline
#  This pipeline is based on the PALseqV3 pipeline 
#  It is modified to use streptavidin fluorescence 
#  as a readout of tail length, and uses the V3.1 version of
#  PALseq library design (based on a design with Coffee)
#  Timothy J Eisen
#  2019 09 19
#####################################################


# Generates a help message. Borrowed from https://github.com/pydanny/cookiecutter-djangopackage.
help: ## Display this help message
	@echo "Please run \`make <inputs as environment variables> <target>\` where <target> is one of"
	@perl -nle'print $& if m{^[\.a-zA-Z_-]+:.*?## .*$$}' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m  %-25s\033[0m %s\n", $$1, $$2}'
	@echo "See README for examples.\n"
	@echo "Messages:\n"
	@echo "Did you set the read length in the config file?"
	@echo "Are you using an NN quality filtering for read two, set in the config file?\n\n"

#Parse the arguments for the datatype
###TJE 20181209 Need to change the lines concerning parsing of PALseq and Tail-seq data to give warnings and exit if requirements aren't met
strand=S
state=True
# readCommand=zcat #Note that this must be set along with the same argument in the config file.
readCommand= - #Note that this must be set along with the same argument in the config file.
clip3pR2=10 #This should probably be 5. 
clip5pR1=14 #???

parseArgs:
	mkdir ${outdir}

demultiplex: 
	#barcode to file name as output
	python /lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqV3/tail-seq-scripts/DemultiplexFastqIntensity.py $(barcodesFile) $(fastq1_combined) $(fastq2_combined) $(intensity_combined) $(outdir)/fastq

align-to-genome: ## Align rest of reads to genome and intersect with gff file
	#Now aligns read 2, not read 1.
	STAR --genomeDir $(genomeDir) --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM jM --alignIntronMax 1 --runThreadN 24 --outFilterMultimapNmax 1 --clip3pNbases $(clip3pR2) --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --readFilesCommand $(readCommand) --outSJfilterReads Unique --readFilesIn $(fastq2) $ --outFileNamePrefix $(outdir)/STAR_ > $(outdir)/stdOut_logFile.txt

intersect-gff:
	bedtools intersect -abam $(outdir)/STAR_Aligned.sortedByCoord.out.bam -b $(gff) -bed -wb -wa -${strand} > $(outdir)/read1.bed

#ARE YOU USING A BEDFILE?
signal-from-raw: ## generate normalized fluorscent signal for each read. 
	python3 /lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqV3/tail-seq-scripts/get_signal_from_raw_PALV3.py --f1 $(fastq1) --f2 $(fastq2) -i $(intensity) -s $(standard_file) -o $(outdir) --strand ${strand} -b

pal-seq: ## 
	python3 /lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqV3/tail-seq-scripts/tail_length_lin_mod.py -o $(outdir)

plot-model:
	head -1 $(outdir)/tail_lengths.txt > $(outdir)/standard_tail_lengths.txt
	head -1 $(outdir)/tail_lengths.txt > $(outdir)/mapped_tail_lengths_stds_removed.txt
	tail -n +2 $(outdir)/tail_lengths.txt | grep standard >> $(outdir)/standard_tail_lengths.txt
	tail -n +2 $(outdir)/tail_lengths.txt | grep -v standard >> $(outdir)/mapped_tail_lengths_stds_removed.txt

	Rscript /lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqV3/tail-seq-scripts/PlotMediansLinFit.R $(outdir)

summary: ## Aggregate individual tail lengths by accession and plot standards
	python3 /lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqV3/tail-seq-scripts/summarize_results.py -o $(outdir)

clean:
	rm $(outdir)/read2temp.bed
	rm $(outdir)/read2filtered.bed
	rm $(outdir)/Read2Filtered_STAR_Aligned.sortedByCoord.out.sam


# all: parseArgs align-to-genome intersect-gff filter_bam signal-from-raw signal-plot tail-seq summary clean ## Run all at once

testingall: parseArgs align-to-genome intersect-gff signal-from-raw pal-seq plot-model summary

testing: signal-from-raw pal-seq plot-model summary

# all_from_bam: intersect-gff filter_bam signal-from-raw signal-plot tail-seq summary ## Run without STAR aligning

# all_hmm: tail-seq summary 

# all-plots: signal-plot summary ## Run all plotting
