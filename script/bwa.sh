#!/bin/env bash
#SBATCH --partition=rimlsfnwi
#SBATCH --job-name=bwa
#SBATCH --output=./log/arr_%x-%A-%a.out
#SBATCH --error=./log/arr_%x-%A-%a.err
#SBATCH --time=1:00:00
#SBATCH --mem 100G
#SBATCH -c 16
#SBATCH --array=1

wd=/ceph/rimlsfnwi/data/cellbio/mhlanga/thsieh
sub=microC

inputDir=$wd/$sub/fastq
outputDir=$wd/$sub/mapped

hg38=$wd/GRCh38/GRCh38.primary_assembly.genome.fa

cd $inputDir

inputfile_list=($inputDir/*.gz)
inputfile=${inputfile_list[$SLURM_ARRAY_TASK_ID-1]}

basename_temp=${inputfile%_R1.fastq.gz}
basename=${basename_temp##*/}

# 43034_THP1_LPS_mc1_R1.fastq.gz 

if [ -d "$outputDir/$basename" ]; then
        echo "outputDir/$basename exists."
        rm -r $outputDir/$basename
fi
mkdir $outputDir/$basename

bwa mem -5SP -T0 -t16 $hg38 $inputDir/${basename}_R1.fastq.gz $inputDir/${basename}_R2.fastq.gz -o $outputDir/$basename/aligned.sam
