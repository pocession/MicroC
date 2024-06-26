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

## Analysis pipeline
We have developed a pipeline that identifies the differentially-interacted regions (DIRs) between two samples. As microC is an expensive assay, we don't do biological repeat for each condition. Instead, we set the biological covariance (BCV) as 0.4 to get the DIRs. We removed interactions with very low abudance (the last 5%) when performing the differential analysis. For more details, please see [this file](./Doc/Il1B_DIR_identification.Rmd).

## Usage of our homebrew tools in this repository
We have developed several Python and R functions for further processing microC data and the downstream analysis. To use those functions, the easiest way is to download this repoitory and unzip it. 

- Please refer to this [Jupyter Notebook](./Python/subset_hic.ipynb) for subsetting hic data. 
- Please refer to this [Rmd file](./Doc/Il1B_DIR_identification.Rmd) for identifying differentially-interacted regions (DIRs).

```
## Use R functions
## In R studio

devtools::document()
devtools::load_all()

library(edgeR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationDbi)
```

Example commands for our homebrew tools:

- [subset_hic_data_v2.py](./Python/subset_hic_data_v2.py): This function helps to subset the interested regions from hic file and generate the position index for further downstream analysis. Please refer to this [tutorial](./Python/subset_hic_tutorial.ipynb) for more information. 

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

