fastq1="/lab/solexa_public/Bartel/171221_WIGTC-HISEQB_HVJN7BCXY/QualityScore/TAGTGC-s_2_1_sequence.txt"
fastq2="/lab/solexa_public/Bartel/171221_WIGTC-HISEQB_HVJN7BCXY/QualityScore/TAGTGC-s_2_2_sequence.txt"
genomeDir="/lab/solexa_bartel/teisen/RNAseq/indices/STAR_mm10/"
gff="/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/Annotations/mm10_tpUTR.gff"
outdir="/lab/bartel4_ata/kathyl/Tail_Seq/testing_outputs"
intensity="/lab/solexa_bartel/tail-seq_data_links/171221_WIGTC-HISEQB_0779_HVJN7BCXY/TAGTGC-s_2_hits.txt.gz"

move_files:
	tar -xvzf $(fastq1).tar.gz -C $(outdir)/
	sed 's/\/1;0//' $(outdir)/$(fastq1) | sed 's/\/1;1//' | sed 's/WIGTC-HISEQ:2://' | fastx_reverse_complement -o $(outdir)/fastq1_revcomp.txt
	rm -rf $(outdir)/lab
	tar -xvzf $(fastq2).tar.gz -C $(outdir)/
	sed 's/\/2;0//' $(outdir)/$(fastq2) | sed 's/\/2;1//' | sed 's/WIGTC-HISEQ:2://' > $(outdir)/fastq2.txt
	rm -rf $(outdir)/lab

# pull out standard sequences from read1
separate_standards:
	python tail-seq-scripts/separate_standards.py -i $(outdir)/fastq1_revcomp.txt -s tail-seq-scripts/standard_sequences.txt -o $(outdir)/fastq1_revcomp_no_standards.txt -p $(outdir)/standards.txt

# align rest of reads to genome and intersect with gff file
align_to_genome:
	STAR --genomeDir $(genomeDir) --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM jM --alignIntronMax 1 --runThreadN 24 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSJfilterReads Unique --readFilesIn $(outdir)/fastq1_revcomp_no_standards.txt $ --outFileNamePrefix $(outdir)/STAR_ > $(outdir)/stdOut_logFile.txt
	bedtools intersect -abam $(outdir)/STAR_Aligned.sortedByCoord.out.bam -b $(gff) -bed -wb -s > $(outdir)/read1.bed

# extract relevant columns
# remove reads that map to multiple genes
# deduplicate reads that map to multiple exons
# define 3' end of read
organize_dedup_read1:
	awk '$$19=="+"{print $$1,$$3-1,$$3,$$19,$$21,$$4}$$19=="-"{print $$1,$$2,$$2+1,$$19,$$21,$$4}' OFS="\t" $(outdir)/read1.bed > $(outdir)/read1_3p_ends.txt
	python tail-seq-scripts/dedup.py -i $(outdir)/read1_3p_ends.txt -o $(outdir)/read1_3p_ends_dedup.txt -d $(outdir)/dropped_reads.txt

# concat standard read information with the rest of the reads
cat_standards:
	cat $(outdir)/read1_3p_ends_dedup.txt $(outdir)/standards.txt | sort -k6,6 > $(outdir)/read1_3p_ends_dedup_with_standards.txt

# call short tails manually and remove from HMM pipeline
filter_short_tails:
	python tail-seq-scripts/filter_read2.py -i $(outdir)/fastq2.txt -m $(outdir)/short_tails.txt -o $(outdir)/matched_filtered.txt

join:
	sort -k1,1 $(outdir)/short_tails.txt | join -1 6 -2 1 $(outdir)/read1_3p_ends_dedup_with_standards.txt - > $(outdir)/short_tails_merged.txt
	sort -k1,1 $(outdir)/matched_filtered.txt | join -1 6 -2 1 $(outdir)/read1_3p_ends_dedup_with_standards.txt - > $(outdir)/matched_filtered_merged.txt

# join filtered list of reads with intensity file
join_intensities:
	zcat $(intensity) | awk '$$4=="TAGTGC"{$$4=$$1":"$$2":"$$3"#"$$4; print $$0}' OFS="\t" | cut -f 4- | sort -k1,1 | join $(outdir)/matched_filtered_merged.txt - > $(outdir)/matched_filtered_intensities.txt


tim_pipeline:
	bash /lab/solexa_bartel/teisen/RNAseq/Scripts/tail_seq/TAILseq_pipeline_SWE.sh -1 /lab/solexa_public/Bartel/171221_WIGTC-HISEQB_HVJN7BCXY/QualityScore/TAGTGC-s_2_1_sequence.txt.tar.gz -s /lab/solexa_bartel/teisen/RNAseq/indices/STAR_mm10/ -b /lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/Annotations/mm10_tpUTR.gff -i /lab/solexa_bartel/tail-seq_data_links/171221_WIGTC-HISEQB_0779_HVJN7BCXY/ -d ../tim_output


# make fastq1=/lab/solexa_public/Bartel/171221_WIGTC-HISEQB_HVJN7BCXY/QualityScore/TAGTGC-s_2_1_sequence.txt fastq2=/lab/solexa_public/Bartel/171221_WIGTC-HISEQB_HVJN7BCXY/QualityScore/TAGTGC-s_2_2_sequence.txt genomeDir=/lab/solexa_bartel/teisen/RNAseq/indices/STAR_mm10/ gff=/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/Annotations/mm10_tpUTR.gff outdir=/lab/bartel4_ata/kathyl/Tail_Seq/testing_outputs move_files align_to_genome assign_genes organize_read1 dedup

# awk '{$4=$1":"$2":"$3"#"$4}1' OFS="\t" ../testing_outputs/mini_intensity.txt 

# /lab/solexa_bartel/tail-seq_data_links/171221_WIGTC-HISEQB_0779_HVJN7BCXY/TAGTGC-s_2_hits.txt.gz 


# .DEFAULT_GOAL := help

# Generates a help message. Borrowed from https://github.com/pydanny/cookiecutter-djangopackage.
# help: ## Display this help message
# 	@echo "Please run \`make <target>\` where <target> is one of"
# 	@perl -nle'print $& if m{^[\.a-zA-Z_-]+:.*?## .*$$}' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m  %-25s\033[0m %s\n", $$1, $$2}'

# signal: ## Normalize T-signal
# 	python tail-seq-scripts/get_signal.py --f1 $(dir)$(tag)_1.txt.gz --l1 50 --f2 $(dir)$(tag)_2.txt.gz --l2 250 -i $(dir)$(tag)_intensity.txt.gz -a $(dir)$(tag)_bowtie.txt -o ../data/outputs/testing2$(tag) -f 24

# signal-plot: ## Plot T-signal density
# 	python tail-seq-scripts/plot_t_signals.py -f ../data/outputs/testing2$(tag)/normalized_t_signal.txt -o ../data/outputs/testing2$(tag)/signals.pdf -b 100 -l 250

# tail-seq: ## Train and run HMM for calling tail lengths
# 	python tail-seq-scripts/tail_length_hmm.py -s ../data/outputs/testing$(tag)/normalized_t_signal.txt -l 250 -t 10000 -m 100 --tol 100 -o ../data/outputs/testing2$(tag) -f 24

# all: signal signal-plot tail-seq ## Run all at once

# make tag=TAGTGC-2 dir=/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/filtered_analysis/matched_fastq/ all