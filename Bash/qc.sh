#!/bin/env bash
#SBATCH --partition=rimlsfnwi
#SBATCH --job-name=qc
#SBATCH --output=./log/arr_%x-%A-%a.out
#SBATCH --error=./log/arr_%x-%A-%a.err
#SBATCH --time=1:00:00
#SBATCH --mem 10G
#SBATCH -c 16
#SBATCH --array=1

# This script is used to generate a QC report for the microC library with less reads (<2M)

wd=/ceph/rimlsfnwi/data/cellbio/mhlanga/thsieh
sub=microC

inputDir=$wd/$sub/fastq
outputDir=$wd/$sub/shallow_QC
finalBam=$wd/$sub/finalBam
deduped=$wd/$sub/deduped

cd $outputDir

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

get_qc=/ceph/rimlsfnwi/data/cellbio/mhlanga/thsieh/App/Micro-C

python3 $get_qc/get_qc.py -p $deduped/$basename/stats.txt > $outputDir/$basename/qc.txt

LD_LIBRARY_PATH=/ceph/rimlsfnwi/data/cellbio/mhlanga/thsieh/App/htslib-1.17/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 -seg_len 1000000000 -output $outputDir/$basename/out.preseq $finalBam/$basename/mapped.PT.bam
