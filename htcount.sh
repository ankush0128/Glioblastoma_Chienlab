#!/bin/bash

#SBATCH --job-name=htcount1
# Project:
#SBATCH --account=NN9632K
# Wall clock limit:
#SBATCH --time=18:00:00
#SBATCH --output=/work/users/ankushs/software/htseq/htseq-%j.out
#SBATCH --error=/work/users/ankushs/software/htseq/htseq-%j.err
# Number of cores:
# Max memory usage:
#SBATCH --mem-per-cpu=6000
# Number of cores:
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ankush.sharma@ibv.uio.no
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END



#If you need hugemem...
#--partition=hugemem
RUNPATH=/work/users/ankushs/software/htseq/
cd $RUNPATH

source /work/users/ankushs/software/htseq/htseq/bin/activate


#transcriptome_path=/usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/bowtie2-index/trans_th2_idx
#ensembl_genome_ref=/usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/bowtie2-index/
#read_path=work/users/ankushs/Miao_glioblastoma/
GTF_path=/usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/gtf/





htseq-count -i gene_id -f bam --additional-attr=gene_name /work/users/ankushs/Miao_glioblastoma/2_tophat/Chops_DvS003Glioma1RNA-chops1/accepted_hits_sorted.bam /usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/gtf/Homo_sapiens.GRCh38.93.gtf -o /work/users/ankushs/Miao_glioblastoma/3_cuffdiff/HTCOUNT/DvS003Glioma1RNA-chops1.txt
htseq-count -i gene_id -f bam --additional-attr=gene_name /work/users/ankushs/Miao_glioblastoma/2_tophat/Chops_DvS003Glioma1RNA-chops2/accepted_hits_sorted.bam /usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/gtf/Homo_sapiens.GRCh38.93.gtf -o /work/users/ankushs/Miao_glioblastoma/3_cuffdiff/HTCOUNT/DvS003Glioma1RNA-chops2.txt
htseq-count -i gene_id -f bam --additional-attr=gene_name /work/users/ankushs/Miao_glioblastoma/2_tophat/Chops_DvS003Glioma1RNA-chops3/accepted_hits_sorted.bam /usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/gtf/Homo_sapiens.GRCh38.93.gtf -o /work/users/ankushs/Miao_glioblastoma/3_cuffdiff/HTCOUNT/DvS003Glioma1RNA-chops3.txt
htseq-count -i gene_id -f bam --additional-attr=gene_name /work/users/ankushs/Miao_glioblastoma/2_tophat/day10_DvS003Glioma1RNA-D-day10/accepted_hits_sorted.bam /usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/gtf/Homo_sapiens.GRCh38.93.gtf -o /work/users/ankushs/Miao_glioblastoma/3_cuffdiff/HTCOUNT/DvS003Glioma1RNA-D-day10.txt
htseq-count -i gene_id -f bam --additional-attr=gene_name /work/users/ankushs/Miao_glioblastoma/2_tophat/day10_DvS003Glioma1RNA-E-day10/accepted_hits_sorted.bam /usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/gtf/Homo_sapiens.GRCh38.93.gtf -o /work/users/ankushs/Miao_glioblastoma/3_cuffdiff/HTCOUNT/DvS003Glioma1RNA-E-day10.txt
htseq-count -i gene_id -f bam --additional-attr=gene_name /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-A-day3/accepted_hits_sorted.bam /usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/gtf/Homo_sapiens.GRCh38.93.gtf -o /work/users/ankushs/Miao_glioblastoma/3_cuffdiff/HTCOUNT/DvS003Glioma1RNA-A-day3.txt
htseq-count -i gene_id -f bam --additional-attr=gene_name /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-B-day3/accepted_hits_sorted.bam /usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/gtf/Homo_sapiens.GRCh38.93.gtf -o /work/users/ankushs/Miao_glioblastoma/3_cuffdiff/HTCOUNT/DvS003Glioma1RNA-B-day3
htseq-count -i gene_id -f bam --additional-attr=gene_name /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-C-day3/accepted_hits_sorted.bam /usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/gtf/Homo_sapiens.GRCh38.93.gtf -o /work/users/ankushs/Miao_glioblastoma/3_cuffdiff/HTCOUNT/DvS003Glioma1RNA-C-day3
htseq-count -i gene_id -f bam --additional-attr=gene_name /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-D-day3/accepted_hits_sorted.bam /usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/gtf/Homo_sapiens.GRCh38.93.gtf -o /work/users/ankushs/Miao_glioblastoma/3_cuffdiff/HTCOUNT/DvS003Glioma1RNA-D-day3
htseq-count -i gene_id -f bam --additional-attr=gene_name /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-G-day3/accepted_hits_sorted.bam /usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/gtf/Homo_sapiens.GRCh38.93.gtf -o /work/users/ankushs/Miao_glioblastoma/3_cuffdiff/HTCOUNT/DvS003Glioma1RNA-G-day3
htseq-count -i gene_id -f bam --additional-attr=gene_name /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-H-day3/accepted_hits_sorted.bam /usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/gtf/Homo_sapiens.GRCh38.93.gtf -o /work/users/ankushs/Miao_glioblastoma/3_cuffdiff/HTCOUNT/DvS003Glioma1RNA-H-day3
htseq-count -i gene_id -f bam --additional-attr=gene_name /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-B-day3/accepted_hits_sorted.bam /usit/abel/u1/ankushs/nobackup/references/ensembl_build93/ensemble/gtf/Homo_sapiens.GRCh38.93.gtf -o /work/users/ankushs/Miao_glioblastoma/3_cuffdiff/HTCOUNT/DvS003Glioma1RNA-B-day3