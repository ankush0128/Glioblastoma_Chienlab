

#!/bin/bash
#SBATCH --job-name=GLIOMA_cuffdiff
# Project:
#SBATCH --account=NN9632K
# Wall clock limit:
#SBATCH --time=12:00:00
# Max memory usage:
#SBATCH --mem-per-cpu=3000MB
# Number of cores:
#SBATCH --cpus-per-task=12
#SBATCH --mail-user=ankush.sharma@ibv.uio.no
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

#If you need hugemem...
#--partition=hugemem


## Set up the job environment
source /cluster/bin/jobsetup
module purge
module load samtools 
set -o errexit



#samtools sort /work/users/ankushs/Miao_glioblastoma/2_tophat/Chops_DvS003Glioma1RNA-chops1/accepted_hits.bam -o /work/users/ankushs/Miao_glioblastoma/2_tophat/Chops_DvS003Glioma1RNA-chops1/accepted_hits_sorted.bam 

samtools sort /work/users/ankushs/Miao_glioblastoma/2_tophat/Chops_DvS003Glioma1RNA-chops2/accepted_hits.bam -o /work/users/ankushs/Miao_glioblastoma/2_tophat/Chops_DvS003Glioma1RNA-chops2/accepted_hits_sorted.bam 

samtools sort /work/users/ankushs/Miao_glioblastoma/2_tophat/Chops_DvS003Glioma1RNA-chops3/accepted_hits.bam -o /work/users/ankushs/Miao_glioblastoma/2_tophat/Chops_DvS003Glioma1RNA-chops3/accepted_hits_sorted.bam 

samtools sort /work/users/ankushs/Miao_glioblastoma/2_tophat/day10_DvS003Glioma1RNA-D-day10/accepted_hits.bam -o /work/users/ankushs/Miao_glioblastoma/2_tophat/day10_DvS003Glioma1RNA-D-day10/accepted_hits_sorted.bam 

samtools sort /work/users/ankushs/Miao_glioblastoma/2_tophat/day10_DvS003Glioma1RNA-E-day10/accepted_hits.bam -o /work/users/ankushs/Miao_glioblastoma/2_tophat/day10_DvS003Glioma1RNA-E-day10/accepted_hits_sorted.bam 

samtools sort /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-A-day3/accepted_hits.bam -o /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-A-day3/accepted_hits_sorted.bam 

samtools sort /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-B-day3/accepted_hits.bam -o /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-B-day3/accepted_hits_sorted.bam 

samtools sort /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-C-day3/accepted_hits.bam -o /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-C-day3/accepted_hits_sorted.bam 

samtools sort /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-D-day3/accepted_hits.bam -o /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-D-day3/accepted_hits_sorted.bam 

samtools sort /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-G-day3/accepted_hits.bam -o /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-G-day3/accepted_hits_sorted.bam 

samtools sort /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-H-day3/accepted_hits.bam -o /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-H-day3/accepted_hits_sorted.bam 

samtools sort /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-B-day3/accepted_hits.bam -o /work/users/ankushs/Miao_glioblastoma/2_tophat/day3_DvS003Glioma1RNA-B-day3/accepted_hits_sorted.bam 



