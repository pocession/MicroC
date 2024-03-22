#!/bin/env bash
#SBATCH --partition=rimlsfnwi
#SBATCH --job-name=parse
#SBATCH --output=./log/arr_%x-%A-%a.out
#SBATCH --error=./log/arr_%x-%A-%a.err
#SBATCH --time=1:00:00
#SBATCH --mem 100G
#SBATCH -c 16
#SBATCH --array=1

wd=/ceph/rimlsfnwi/data/cellbio/mhlanga/thsieh
sub=microC

inputDir=$wd/$sub/fastq
outputDir=$wd/$sub/parsed

mapped=$wd/$sub/mapped
hg38genome=$wd/GRCh38/GRCh38.primary_assembly.genome

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

pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 \
 --chroms-path $hg38genome $mapped/$basename/aligned.sam > $outputDir/$basename/parsed.pairsam
