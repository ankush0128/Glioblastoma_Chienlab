#!/bin/bash
#SBATCH --job-name=glioblastoma_2
# Project:
#SBATCH --account=NN9632K
# Wall clock limit:
#SBATCH --time=8:00:00
# Max memory usage:
#SBATCH --mem-per-cpu=2G
# Number of cores:
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=ankush.sharma@ibv.uio.no
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

#If you need hugemem...
#--partition=hugemem


## Set up the job environment
source /cluster/bin/jobsetup
module purge


mkdir DvS003Glioma1RNA-D-day3
cd DvS003Glioma1RNA-D-day3
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-D-day3/V300009963_L3_DKRT190110OligoTB22-22_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-D-day3/V300009963_L3_DKRT190110OligoTB22-22_2.fq.gz

cd ..
mkdir DvS003Glioma1RNA-E-day10
cd DvS003Glioma1RNA-E-day10
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-E-day10/V300009963_L4_DKRT190110OligoTB26-26_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-E-day10/V300009963_L4_DKRT190110OligoTB26-26_2.fq.gz

cd ..
mkdir DvS003Glioma1RNA-G-day3
cd 
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-G-day3/V300009963_L3_DKRT190110OligoTB23-23_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-G-day3/V300009963_L3_DKRT190110OligoTB23-23_2.fq.gz

cd ..
mkdir DvS003Glioma1RNA-H-day3
cd
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-H-day3/V300009963_L3_DKRT190110OligoTB24-24_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-H-day3/V300009963_L3_DKRT190110OligoTB24-24_2.fq.gzcd ..
mkdir DvS003Glioma1RNA-D-day3
cd DvS003Glioma1RNA-D-day3
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-D-day3/V300009963_L3_DKRT190110OligoTB22-22_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-D-day3/V300009963_L3_DKRT190110OligoTB22-22_2.fq.gz

cd ..
mkdir DvS003Glioma1RNA-E-day10
cd DvS003Glioma1RNA-E-day10
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-E-day10/V300009963_L4_DKRT190110OligoTB26-26_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-E-day10/V300009963_L4_DKRT190110OligoTB26-26_2.fq.gz

cd ..
mkdir DvS003Glioma1RNA-G-day3
cd 
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-G-day3/V300009963_L3_DKRT190110OligoTB23-23_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-G-day3/V300009963_L3_DKRT190110OligoTB23-23_2.fq.gz

cd ..
mkdir DvS003Glioma1RNA-H-day3
cd
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-H-day3/V300009963_L3_DKRT190110OligoTB24-24_1.fq.gz
wget --user=F18FTSEUHT1779_HUMhgnR --password='ekL3b4' ftp://5.57.48.133/Clean/DvS003Glioma1RNA-H-day3/V300009963_L3_DKRT190110OligoTB24-24_2.fq.gz


echo completed
