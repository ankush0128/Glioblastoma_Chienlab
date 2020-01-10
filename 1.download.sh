#!/bin/bash
#SBATCH --job-name=glioblastoma_1
# Project:
#SBATCH --account=nn9374k
# Wall clock limit:
#SBATCH --time=8:00:00
# Max memory usage:
#SBATCH --mem-per-cpu=2G
# Number of cores:
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=ankush.sharma@ibv.uio.no
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

#If you need hugemem...
#--partition=hugemem


## Set up the job environment
source /cluster/bin/jobsetup
module purge
#module load openmpi.intel
#module load tophat
#module load bowtie2
#module load samtools
#module load python3



cd DvS003Glioma1RNA-A-day3
#wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-A-day3/V300009963_L3_DKRT190110OligoTB19-19_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-A-day3/V300009963_L3_DKRT190110OligoTB19-19_2.fq.gz
cd .. 
mkdir DvS003Glioma1RNA-B-day3
cd DvS003Glioma1RNA-B-day3
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-B-day3/V300009963_L3_DKRT190110OligoTB20-20_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-B-day3/V300009963_L3_DKRT190110OligoTB20-20_2.fq.gz

cd ..
mkdir DvS003Glioma1RNA-C-day3
cd DvS003Glioma1RNA-C-day3
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-C-day3/V300009963_L3_DKRT190110OligoTB21-21_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-C-day3/V300009963_L3_DKRT190110OligoTB21-21_2.fq.gz

cd ..
mkdir DvS003Glioma1RNA-chops1
cd DvS003Glioma1RNA-chops1
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-chops1/V300009963_L4_DKRT190110OligoTB27-27_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-chops1/V300009963_L4_DKRT190110OligoTB27-27_2.fq.gz

cd ..
mkdir DvS003Glioma1RNA-chops2
cd DvS003Glioma1RNA-chops1
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-chops2/V300009963_L4_DKRT190110OligoTB28-28_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-chops2/V300009963_L4_DKRT190110OligoTB28-28_2.fq.gz

cd ..
mkdir DvS003Glioma1RNA-chops3
cd DvS003Glioma1RNA-chops3
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-chops3/V300009963_L4_DKRT190110OligoTB29-29_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-chops3/V300009963_L4_DKRT190110OligoTB29-29_2.fq.gz

cd ..
mkdir DvS003Glioma1RNA-D-day10
cd DvS003Glioma1RNA-D-day10
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-D-day10/V300009963_L4_DKRT190110OligoTB25-25_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-D-day10/V300009963_L4_DKRT190110OligoTB25-25_2.fq.gz

echo completed Successfully
