signal:
	python tail-seq-scripts/get_signal.py --f1 $(dir)$(tag)_1.txt.gz --l1 50 --f2 $(dir)$(tag)_2.txt.gz --l2 250 -i $(dir)$(tag)_intensity.txt.gz -a $(dir)$(tag)_bowtie.txt -o ../data/outputs/testing2$(tag) -f 24

tail-seq:
	python tail-seq-scripts/tail_length_hmm.py -s ../data/outputs/testing$(tag)/normalized_t_signal.txt -l 250 -t 10000 -m 100 --tol 100 -o ../data/outputs/testing2$(tag) -f 24

signal-plot:
	python tail-seq-scripts/plot_t_signals.py -f ../data/outputs/testing2$(tag)/normalized_t_signal.txt -o ../data/outputs/testing2$(tag)/signals.pdf -b 100 -l 250

all:
	python tail-seq-scripts/get_signal.py --f1 $(dir)$(tag)_1.txt.gz --l1 50 --f2 $(dir)$(tag)_2.txt.gz --l2 250 -i $(dir)$(tag)_intensity.txt.gz -a $(dir)$(tag)_bowtie.txt -o ../data/outputs/testing2$(tag) -f 24 \
	&& python tail-seq-scripts/tail_length_hmm.py -s ../data/outputs/testing$(tag)/normalized_t_signal.txt -l 250 -t 10000 -m 100 --tol 100 -o ../data/outputs/testing2$(tag) -f 24 \
	&& python tail-seq-scripts/plot_t_signals.py -f ../data/outputs/testing2$(tag)/normalized_t_signal.txt -o ../data/outputs/testing2$(tag)/signals.pdf -b 100 -l 250


# make tag=TAGTGC-2 dir=/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/filtered_analysis/matched_fastq/ all