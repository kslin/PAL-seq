stephen:
	python stephen_script/Tail-seq_gmhmm_v6_50by250.py data/TAGTGC-1_1.txt.gz data/TAGTGC-1_bowtie.txt data/TAGTGC-1_2.txt.gz data/TAGTGC-1_intensity.txt.gz data/TAGTGC-1_individual_lengths_test.txt data/TAGTGC-1_median_lengths_test.txt

tail-seq:
	python tail-seq-scripts/tail_seq_hmm.py --f1 data/TAGTGC-1_1.txt.gz --l1 50 --f2 data/TAGTGC-1_2.txt.gz --l2 250 -m data/TAGTGC-1_bowtie.txt -i data/TAGTGC-1_intensity.txt.gz -o data/outputs -s data/listofstandards.txt
