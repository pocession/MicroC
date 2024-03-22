#!/bin/env bash
#SBATCH --partition=rimlsfnwi
#SBATCH --job-name=process
#SBATCH --output=./log/arr_%x-%A-%a.out
#SBATCH --error=./log/arr_%x-%A-%a.err
#SBATCH --time=1:00:00
#SBATCH --mem 100G
#SBATCH -c 16
#SBATCH --array=1

# This script is used to process microC data from fastq to bam.

# Set path
wd=/ceph/rimlsfnwi/data/cellbio/mhlanga/thsieh
sub=microC

inputDir=$wd/$sub/fastq
outputDir=$wd/$sub/output

# You should have hg38 genome installed somewhere
hg38=$wd/GRCh38/GRCh38.primary_assembly.genome.fa
hg38genome=$wd/GRCh38/GRCh38.primary_assembly.genome

# Specifiy a temp path for temporary files
temp=$wd/$sub/temp

cd $inputDir

# Get the base name of input files
inputfile_list=($inputDir/*.gz)
inputfile=${inputfile_list[$SLURM_ARRAY_TASK_ID-1]}

basename_temp=${inputfile%_R1.fastq.gz}
basename=${basename_temp##*/}

# Example input files: 43034_THP1_LPS_mc1_R1.fastq.gz 

# Check whether if there is already output dir
# If no, create a new one
if [ -d "$outputDir/$basename" ]; then
        echo "outputDir/$basename exists."
        rm -r $outputDir/$basename
fi
mkdir $outputDir/$basename

if [ -d "$temp/$basename" ]; then
        echo "temp/$basename exists."
        rm -r $temp/$basename
fi
mkdir $temp/$basename

# The real process pipeline
bwa mem -5SP -T0 -t16 $hg38 $inputDir/${basename}_R1.fastq.gz $inputDir/${basename}_R2.fastq.gz| \
pairtools parse --min-mapq 40 --walks-policy 5unique \
--max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path $hg38genome | \
pairtools sort --tmpdir=$temp --nproc 16|pairtools dedup --nproc-in 8 \
--nproc-out 8 --mark-dups --output-stats $outputDir/$basename/stats.txt|pairtools split --nproc-in 8 \
--nproc-out 8 --output-pairs $outputDir/$basename/mapped.pairs --output-sam -|samtools view -bS -@16 | \
samtools sort -@16 -o $outputDir/$basename/mapped.PT.bam;samtools index $outputDir/$basename/mapped.PT.bam
