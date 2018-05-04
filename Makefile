.DEFAULT_GOAL := help

# Generates a help message. Borrowed from https://github.com/pydanny/cookiecutter-djangopackage.
help: ## Display this help message
	@echo "Please run \`make <inputs as environment variables> <target>\` where <target> is one of"
	@perl -nle'print $& if m{^[\.a-zA-Z_-]+:.*?## .*$$}' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m  %-25s\033[0m %s\n", $$1, $$2}'
	@echo "See README for examples."

align-to-genome: ## Align rest of reads to genome and intersect with gff file
	STAR --genomeDir $(genomeDir) --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM jM --alignIntronMax 1 --runThreadN 24 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSJfilterReads Unique  --readFilesCommand zcat --readFilesIn $(fastq1) $ --outFileNamePrefix $(outdir)/STAR_ > $(outdir)/stdOut_logFile.txt

intersect-gff:
	bedtools intersect -abam $(outdir)/STAR_Aligned.sortedByCoord.out.bam -b $(gff) -bed -wb -S > $(outdir)/read1.bed

signal-from-raw: ## Extract intensities for mapping reads and calculate normalized T-signal
	python tail-seq-scripts/get_signal_from_raw.py --f1 $(fastq1) --f2 $(fastq2) -i $(intensity) -s $(standard_file) -o $(outdir)

signal-plot: ## Plot T-signal density
	python tail-seq-scripts/plot_t_signals.py -o $(outdir) -b 100

tail-seq: ## Train and run HMM for calling tail lengths
	python tail-seq-scripts/tail_length_hmm.py -o $(outdir)

summary: ## Aggregate individual tail lengths by accession and plot standards
	python tail-seq-scripts/summarize_results.py -o $(outdir)

all: align-to-genome intersect-gff signal-from-raw signal-plot tail-seq summary ## Run all at once

all_from_bam: intersect-gff signal-from-raw signal-plot tail-seq summary ## Run without STAR aligning

all-plots: signal-plot summary ## Run all plotting
