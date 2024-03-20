# MicroC
This repository is a resource and a toolbox for microC-seq analysis.

## About MicroC
Please read the [original paper](https://www.sciencedirect.com/science/article/pii/S0092867415006388) for the principle and details.
- PS: My name is as same as the first author of this paper, but I have nothing to do with this paper.

MicroC is a sequencing technology that captures the genomic interaction in high resolution. In theory, the resolution could reach to single nucleosome (~150 base pair). In practice, the resolution usually ranges within 1000 - 5000 bp.

## The experiment material
There are multiple wasy to perfrom MicroC experiment. We have tested [Dovetail MicroC kit](https://cantatabio.com/dovetail-genomics/products/micro-c-sequencing/) and obtain a good result.

## Protocols
Please refer to the [Dovetail official website](https://cantatabio.com/dovetail-genomics/products/micro-c-sequencing/) to have the original protocols. 

We have tested this kit for LPS-/beta-glucan-stimulated human monocytes (primary cells). To collect cells after stimulation, please add 10% FBS in PBS during all centrifugation steps.

## Analysis
For a step-by-step guide, please refer to [Dovetail analsyis documents](https://micro-c.readthedocs.io/en/latest/) for more information. Alternatively, you can also check [Juicebox](https://github.com/aidenlab/Juicebox) for a more generalized pipeline.

In [our bash repository](./Bash), you can find example bash scripts for the pre-processing step in high-performance computation (HPC) cluster.

## Usages of our homebrew tools in this repository
We have developed several Python and R functions for further processing microC data and the downstream analysis. To use those functions, the easiest way is to download this repoitory and unzip it. Then navigate to this folder and you are ready to go.

```
## In linux environment
unzip MicroC.zip -d -path-to-this-folder
cd -path-to-this-folder
```

Example commands for our homebrew tools:

- [subset_hic_data.py](./Python/subset_hic_data.py): This function helps to subset the interested regions from hic file and generate the position index for further downstream analysis. Note the position index is the mid-point of the genomic regions. For example, if there is a interaction between two genomic regions: 11000-16000 and 18000-23000, then the index will be 13500 - 20500.

```
## Example
python ./Python/subset_hic_data.py --inputDir /data/44114_C_lps_43617_mc8contact_map.hic --outputDir /Results/processing/ --rowChr chr2 --columnChr chr2 --rowStart 112735986 --rowEnd 113204585 --columnStart 112735986 --columnEnd 112735986 --res 5000
```

- [getInteractionPosIndex.R](./R/getInteractionPosIndex.R): This function is an independent R function for generating the interaction position index.

```
## Example:
df <- getInteractionPosIndex(
chr = "chr2",
start = 112735986,
end = 113204585,
res = 5000,
output = here::here("./Results/processing/start_position_index.txt")
)
```

- [getDIRWithNoReplicate.R](./R/getDIRWithNoReplicate.R): After running subset_hic_data.py, run this function to read interaction counts and perform differential analysis between two conditions. Please note this function is specifically for experiments with no biological replicate, so the result has no statistical significance. Please treat the log2 fold change and p value as descriptive values.

```
## Example
df <- GetDIRWithNoReplicate(
chr = "chr2",
treat = here::here("./Results/processing/44112_A_bg_43615_mc6contact_map_extracted.csv"),
ctrl = here::here("./Results/processing/44111_A_ctrl_43614_mc5contact_map_extracted.csv"),
bcv = 0.4,
start_position_index = here::here("./Results/processing/start_position_index.txt"),
output = NULL)
```

- [annotateDIR.R](./R/annotateDIR.R): This function annotates the differentially-interacted regions (DIRs) with hg38 genome information.

```
## Example
df <- annotateDIR(
input = here::here("./Results/DIR/DIR_44112_A_bg_43615_mc6contact_map_extracted.csv"),
output = NULL)
```

