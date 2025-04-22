#!/bin/bash
#SBATCH --nodes=1
#SBATCH -n 12
#SBATCH --time=04:00:00
#SBATCH --job-name=DADA2
#SBATCH -A microbiome
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 160 \
  --p-trunc-len-r 200 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza
