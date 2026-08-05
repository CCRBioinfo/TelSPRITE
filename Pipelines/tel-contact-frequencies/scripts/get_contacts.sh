#!/bin/bash

module load python

############
#Script locations
############

get_sprite_contacts=scripts/python/get_sprite_contacts.py
hicorrector=scripts/HiCorrector_1.2/bin/ic
plot_heatmap=scripts/r/plot_heatmap.R


############
#Options
############

max_cluster_size=1000000
min_cluster_size=2
iterations=100
downweighting=two_over_n
max_val=255

############
#Make contact matrices
############

output=$2_final.txt

python $get_sprite_contacts \
	--clusters $1 \
	--raw_contacts $2_raw_contacts.txt \
	--biases $2_biases.txt \
	--iced $2_iced.txt \
	--output $output \
	--assembly $3 \
	--chromosome $4 \
	--max_cluster_size $max_cluster_size \
	--min_cluster_size $min_cluster_size \
	--resolution $5 \
	--iterations $iterations \
	--downweighting two_over_n \
	--hicorrector $hicorrector
