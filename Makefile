stephen:
	python stephen_script/Tail-seq_gmhmm_v6_50by250.py ../data/TAGTGC-1_1.txt.gz ../data/TAGTGC-1_bowtie.txt ../data/TAGTGC-1_2.txt.gz ../data/TAGTGC-1_intensity.txt.gz ../data/TAGTGC-1_individual_lengths_test.txt ../data/TAGTGC-1_median_lengths_test.txt

signal:
	python tail-seq-scripts/get_signal.py --f1 $(dir)$(tag)_1.txt.gz --l1 50 --f2 $(dir)$(tag)_2.txt.gz --l2 250 -i $(dir)$(tag)_intensity.txt.gz -a $(dir)$(tag)_bowtie.txt -o ../data/outputs/testing$(tag) -f 10

tail-seq:
	python tail-seq-scripts/tail_length_hmm.py -s ../data/outputs/testing$(tag)/normalized_t_signal.txt -l 250 -t 10000 -m 100 --tol 100 -o ../data/outputs/testing$(tag)

compare:
	python compare_results.py ../data/outputs/AATCCG-3/AATCCG-3_individual_lengths.txt ../data/outputs/AATCCG-3/tail_lengths.txt

sliding:
	python tail-seq-scripts/tail_length_sliding_window.py --f1 $(dir)$(tag)_1.txt.gz --l1 50 --f2 $(dir)$(tag)_2.txt.gz --l2 250 -i $(dir)$(tag)_intensity.txt.gz -a $(dir)$(tag)_bowtie.txt -o ../data/outputs/sliding_window -f 0

sliding2:
	python tail-seq-scripts/tail_length_sliding_window.py --f1 $(dir)$(tag)_1.txt.gz --l1 50 --f2 $(dir)$(tag)_2.txt.gz --l2 250 -i ../data/outputs/sliding_window/standards.txt -a $(dir)$(tag)_bowtie.txt -o ../data/outputs/sliding_window -f 0


# bsub make tag=AATCCG-3 dir=/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_timecourse/170120/matched_fastq_revised_2/ tail-seq
# bsub make tag=TGAACT-2 dir=/lab/solexa_bartel/eichhorn/Pail-seq/170411/matched_fastq/ signal

# make tag=TAGTGC-2 dir=/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_timecourse/170530/matched_fastq/ sliding