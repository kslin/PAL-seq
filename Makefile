stephen:
	python stephen_script/Tail-seq_gmhmm_v6_50by250.py ../data/TAGTGC-1_1.txt.gz ../data/TAGTGC-1_bowtie.txt ../data/TAGTGC-1_2.txt.gz ../data/TAGTGC-1_intensity.txt.gz ../data/TAGTGC-1_individual_lengths_test.txt ../data/TAGTGC-1_median_lengths_test.txt

tail-seq-pipeline:
	python tail-seq-scripts/tail_seq_hmm.py --f1 ../data/TAGTGC-1_1.txt.gz --l1 50 --f2 ../data/TAGTGC-1_2.txt.gz --l2 250 -m ../data/TAGTGC-1_bowtie.txt -i ../data/TAGTGC-1_intensity.txt.gz -o ../data/outputs -s ../data/listofstandards.txt
# run-mirna:
# 	bsub $(mirna).py

	# make mirna=let7 run-mirna

signal:
	python tail-seq-scripts/get_signal.py --f1 ../data/with_standards/TAGTGC-1_1.txt.gz --l1 50 --f2 ../data/with_standards/TAGTGC-1_2.txt.gz --l2 250 -i ../data/with_standards/TAGTGC-1_intensity.txt.gz -o ../data/output/single_sigma/signals2.csv

tail-seq-EM:
	python tail-seq-scripts/EM.py -s ../data/output/single_sigma/signals2.csv -l 250 -m 10 -tol 10 -o ../data/output/single_sigma2/EM_tails.txt -p ../data/output/single_sigma2/ -f