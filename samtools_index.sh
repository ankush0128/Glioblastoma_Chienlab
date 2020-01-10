
#!/bin/bash
#SBATCH --job-name=Glioma1_tophat2
# Project:
#SBATCH --account=
# Wall clock limit:
#SBATCH --time=12:00:00
# Max memory usage:
#SBATCH --mem-per-cpu=4000
# Number of cores:
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ankush.sharma@ibv.uio.no
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

#If you need hugemem...
#--partition=hugemem


## Set up the job environment
source /cluster/bin/jobsetup
module purge
module load samtools

samtools index /work/users/ankushs/Miao_glioblastoma/2_tophat/DvS003Glioma1RNA-A-day3/accepted_hits.bam
samtools index /work/users/ankushs/Miao_glioblastoma/2_tophat/DvS003Glioma1RNA-D-day10/accepted_hits.bam
samtools index /work/users/ankushs/Miao_glioblastoma/2_tophat/DvS003Glioma1RNA-B-day3/accepted_hits.bam
samtools index /work/users/ankushs/Miao_glioblastoma/2_tophat/DvS003Glioma1RNA-D-day3/accepted_hits.bam
samtools index /work/users/ankushs/Miao_glioblastoma/2_tophat/DvS003Glioma1RNA-C-day3/accepted_hits.bam 
samtools index /work/users/ankushs/Miao_glioblastoma/2_tophat/DvS003Glioma1RNA-E-day10/accepted_hits.bam
samtools index /work/users/ankushs/Miao_glioblastoma/2_tophat/DvS003Glioma1RNA-chops1/accepted_hits.bam
samtools index /work/users/ankushs/Miao_glioblastoma/2_tophat/DvS003Glioma1RNA-G-day3/accepted_hits.bam
samtools index /work/users/ankushs/Miao_glioblastoma/2_tophat/DvS003Glioma1RNA-H-day3/accepted_hits.bam
samtools index /work/users/ankushs/Miao_glioblastoma/2_tophat/DvS003Glioma1RNA-chops2/accepted_hits.bam
samtools index /work/users/ankushs/Miao_glioblastoma/2_tophat/DvS003Glioma1RNA-chops3/accepted_hits.bam

