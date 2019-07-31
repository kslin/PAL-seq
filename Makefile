  
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
clip5pR1?=0

#Parse the arguments for the datatype
###TJE 20181209 Need to change the lines concerning parsing of PALseq and Tail-seq data to give warnings and exit if requirements aren't met
ifeq (${DataType}, PALseq)
	strand=S
	state=True
	readCommand=tar xzfO #Note that this must be set along with the same argument in the config file.
# 	readCommand=zcat #Note that this must be set along with the same argument in the config file.
	clip5pR2=0

else
	strand=s
	state=True
endif

#Parse the arguments for the Phasing
ifeq (${Phasing}, True)
	clip5pR1=12
endif

#NEED TO CHANGE TO .GZ, NOT .TAR.GZ
parseArgs:
	mkdir ${outdir}
	@echo ${strand}
	@echo ${DataType}
	@echo ${Phasing}

align-to-genome: ## Align rest of reads to genome and intersect with gff file
	STAR --genomeDir $(genomeDir) --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM jM --alignIntronMax 1 --runThreadN 24 --outFilterMultimapNmax 1 --clip5pNbases ${clip5pR1} --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --readFilesCommand $(readCommand) --outSJfilterReads Unique --readFilesIn $(fastq1) $ --outFileNamePrefix $(outdir)/STAR_ > $(outdir)/stdOut_logFile.txt
	
	#read2, only 50 nt, clipping the first 4 if direct ligation. 
	STAR --genomeDir $(genomeDir) --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM jM --alignIntronMax 1 --runThreadN 24 --outFilterMultimapNmax 1 --clip5pNbases ${clip5pR2} --clip3pNbases 200 --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --readFilesCommand $(readCommand) --outSJfilterReads Unique --readFilesIn $(fastq2) --outFileNamePrefix $(outdir)/Read2STAR_ > $(outdir)/stdOut_logFile_read2.txt

intersect-gff:
	bedtools intersect -abam $(outdir)/STAR_Aligned.sortedByCoord.out.bam -b $(gff) -bed -wb -wa -${strand} > $(outdir)/read1.bed
	bedtools intersect -abam $(outdir)/Read2STAR_Aligned.sortedByCoord.out.bam -b $(gff) -bed -wb -wa -${strand} > $(outdir)/read2.bed

filter_bam:
	## Implements a filtering step for the Read2STAR_Aligned.sortedByCoord.out.bam file.
	cut -f 4,21 $(outdir)/read2.bed > $(outdir)/read2temp.bed
	cut -f 4,21 $(outdir)/read1.bed | grep -Ff - $(outdir)/read2temp.bed | cut -f 1 > $(outdir)/read2filtered.bed

	#get the bam file header. 
	samtools view -H $(outdir)/Read2STAR_Aligned.sortedByCoord.out.bam > $(outdir)/Read2Filtered_STAR_Aligned.sortedByCoord.out.sam

	#filter the bam file, add it to the sam file.
	samtools view -h $(outdir)/Read2STAR_Aligned.sortedByCoord.out.bam | grep -Ff $(outdir)/read2filtered.bed >> $(outdir)/Read2Filtered_STAR_Aligned.sortedByCoord.out.sam

	#convert back to bam.
	samtools view -b $(outdir)/Read2Filtered_STAR_Aligned.sortedByCoord.out.sam > $(outdir)/Read2Filtered_STAR_Aligned.sortedByCoord.out.bam

signal-from-raw: ## Extract intensities for mapping reads and calculate normalized T-signal
	python3 /lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqKlinMod/tail-seq-scripts/get_signal_from_raw.py --f1 $(fastq1) --f2 $(fastq2) -i $(intensity) -s $(standard_file) -o $(outdir) --strand ${strand}

signal-plot: ## Plot T-signal density
	python3 /lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqKlinMod/tail-seq-scripts/plot_t_signals.py -o $(outdir) -b 100

tail-seq: ## Train and run HMM for calling tail lengths
	python3 /lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqKlinMod/tail-seq-scripts/tail_length_hmm.py -o $(outdir) --twostate ${state} --futures 24

summary: ## Aggregate individual tail lengths by accession and plot standards
	python3 /lab/solexa_bartel/teisen/RNAseq/Scripts/PALseqKlinMod/tail-seq-scripts/summarize_results.py -o $(outdir)

clean:
	rm $(outdir)/read2temp.bed
	rm $(outdir)/read2filtered.bed
	rm $(outdir)/Read2Filtered_STAR_Aligned.sortedByCoord.out.sam


all: parseArgs align-to-genome intersect-gff filter_bam signal-from-raw signal-plot tail-seq summary clean ## Run all at once

all_no_parse: align-to-genome intersect-gff filter_bam signal-from-raw signal-plot tail-seq summary clean ## Run all at once

all_from_bam: intersect-gff filter_bam signal-from-raw signal-plot tail-seq summary ## Run without STAR aligning

all_from_bed: signal-from-raw signal-plot tail-seq summary ##Run without intersect

all_hmm: tail-seq summary 

all-plots: signal-plot summary ## Run all plotting
