.DEFAULT_GOAL := help

# Generates a help message. Borrowed from https://github.com/pydanny/cookiecutter-djangopackage.
help: ## Display this help message
	@echo "Please run \`make <inputs as environment variables> <target>\` where <target> is one of"
	@perl -nle'print $& if m{^[\.a-zA-Z_-]+:.*?## .*$$}' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m  %-25s\033[0m %s\n", $$1, $$2}'
	@echo "See README for examples.\n"
	@echo "Messages:\n"
	@echo "Did you set the read length in the config file?"
	@echo "Are you using an NN quality filtering for read two, set in the config file?\n\n"


DataType?=Tail-seq
Phasing?=False
clipNo?=0

#Parse the arguments for the datatype
###TJE 20181209 Need to change the lines concerning parsing of PALseq and Tail-seq data to give warnings and exit if requirements aren't met
ifeq (${DataType}, PALseq)
	strand=S
	state=True
else
	strand=s
	state=True
endif

#Parse the arguments for the Phasing
ifeq (${Phasing}, True)
	clipNo=12
endif

#NEED TO CHANGE TO .GZ, NOT .TAR.GZ
parseArgs:
	mkdir ${outdir}
	@echo ${strand}
	@echo ${DataType}
	@echo ${Phasing}

align-to-genome: ## Align rest of reads to genome and intersect with gff file
	STAR --genomeDir $(genomeDir) --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM jM --alignIntronMax 1 --runThreadN 24 --outFilterMultimapNmax 1 --clip5pNbases ${clipNo} --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --readFilesCommand zcat --outSJfilterReads Unique --readFilesIn $(fastq1) $ --outFileNamePrefix $(outdir)/STAR_ > $(outdir)/stdOut_logFile.txt

intersect-gff:
	bedtools intersect -abam $(outdir)/STAR_Aligned.sortedByCoord.out.bam -b $(gff) -bed -wb -wa -${strand} > $(outdir)/read1.bed

signal-from-raw: ## Extract intensities for mapping reads and calculate normalized T-signal
	python /lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqKlinMod/tail-seq-scripts/get_signal_from_raw.py --f1 $(fastq1) --f2 $(fastq2) -i $(intensity) -s $(standard_file) -o $(outdir) --strand ${strand}

signal-plot: ## Plot T-signal density
	python /lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqKlinMod/tail-seq-scripts/plot_t_signals.py -o $(outdir) -b 100

tail-seq: ## Train and run HMM for calling tail lengths
	python /lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqKlinMod/tail-seq-scripts/tail_length_hmm.py -o $(outdir) --twostate ${state}

summary: ## Aggregate individual tail lengths by accession and plot standards
	python /lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqKlinMod/tail-seq-scripts/summarize_results.py -o $(outdir)

all: parseArgs align-to-genome intersect-gff signal-from-raw signal-plot tail-seq summary ## Run all at once

all_from_bam: intersect-gff signal-from-raw signal-plot tail-seq summary ## Run without STAR aligning

all_from_bed: signal-from-raw signal-plot tail-seq summary ##Run without intersect

all_hmm: tail-seq summary 

all-plots: signal-plot summary ## Run all plotting
