move-files:
	mkdir -p $(outdir)
	tar -xOzf $(fastq1_tar) | gzip > $(fastq1)
	tar -xOzf $(fastq2_tar) | gzip > $(fastq2)

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

all-plots: signal-plot summary ## Run all plotting
