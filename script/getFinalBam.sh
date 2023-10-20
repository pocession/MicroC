#!/bin/env bash
#SBATCH --partition=rimlsfnwi
#SBATCH --job-name=getFinalBam
#SBATCH --output=./log/arr_%x-%A-%a.out
#SBATCH --error=./log/arr_%x-%A-%a.err
#SBATCH --time=1:00:00
#SBATCH --mem 100G
#SBATCH -c 16
#SBATCH --array=1

wd=/ceph/rimlsfnwi/data/cellbio/mhlanga/thsieh
sub=microC

inputDir=$wd/$sub/fastq
outputDir=$wd/$sub/finalBam
split=$wd/$sub/split
temp=$wd/$sub/temp

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

samtools sort -@16 -T $temp/$basename/temp.bam -o $outputDir/$basename/mapped.PT.bam $split/$basename/unsorted.bam
samtools index $outputDir/$basename/mapped.PT.bam
