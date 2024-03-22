#!/bin/env bash
#SBATCH --partition=rimlsfnwi
#SBATCH --job-name=fastqc
#SBATCH --output=./log/arr_%x-%A-%a.out
#SBATCH --error=./log/arr_%x-%A-%a.err
#SBATCH --time=1:00:00
#SBATCH --mem 10G
#SBATCH -c 2

wd=/ceph/rimlsfnwi/data/cellbio/mhlanga/thsieh
sub=microC

inputDir=$wd/$sub/fastq

cd $inputDir

fastqc *
