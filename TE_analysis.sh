#!/bin/bash
#SBATCH --job-name=TEtrans            # Job name
#SBATCH --partition=batch             # Partition (queue) name                  
#SBATCH --ntasks=1                    # Run a single task       
#SBATCH --cpus-per-task=4             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=TE.%j.out            # Standard output log
#SBATCH --error=TE.%j.err             # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mm82881@uga.edu   # Where to send mail 

cd $SLURM_SUBMIT_DIR
ml TEtranscripts
TEtranscripts --format BAM -t morc3ab-24hpf_rep1_sorted.bam \
morc3ab-24hpf_rep2_sorted.bam \
morc3ab-24hpf_rep3_sorted.bam \
-c wt_24hpf_rep1_sorted.bam \
wt_24hpf_rep2_sorted.bam \
wt_24hpf_rep3_sorted.bam \
--GTF refann.gtf \
--TE TE_annotations_DEclusters_new.gtf \
--mode multi --project morc3ab_24hpf --minread 1 -i 10 --padj 0.05 --sortByPos 
