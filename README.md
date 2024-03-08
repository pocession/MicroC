# MicroC
This repository is a resource and a toolbox for microC-seq analysis.

## About MicroC
Please read the [original paper](https://www.sciencedirect.com/science/article/pii/S0092867415006388) for the principle and details.

MicroC is a sequencing technology that captures the genomic interaction in high resolution. In theory, the resolution could reach to single nucleosome (~150 base pair). In practice, the resolution usually ranges within 1000 - 5000 bp.

## The experiment material
There are multiple wasy to perfrom MicroC experiment. We have tested [Dovetail MicroC kit](https://cantatabio.com/dovetail-genomics/products/micro-c-sequencing/) and obtain a good result.

## Protocols
Please refer to the [Dovetail official website](https://cantatabio.com/dovetail-genomics/products/micro-c-sequencing/) to have the original protocols. 

We have tested this kit for LPS-/beta-glucan-stimulated human monocytes (primary cells). To collect cells after stimulation, please add 10% FBS in PBS during all centrifugation steps.

## Analysis
For a step-by-step guide, please refer to [Dovetail analsyis documents](https://micro-c.readthedocs.io/en/latest/) for more information. Alternatively, you can also check [Juicebox](https://github.com/aidenlab/Juicebox) for a more generalized pipeline.

In [our bash repository](./Bash), you can find example bash scripts for the pre-processing step in high-performance computation (HPC) cluster.

## Other tools in this repository
We have developed several Python and R functions for further processing microC data and the downstream analysis.

- [subset_hic_data.py](./Python/subset_hic_data.py): This function helps subset the interested regions from hic file for further downstream analysis.

```
## Example
python subset_hic_data.py --inputDir ./data/input.hiC ./Results/output.csv --outputDir --startChr chr2 --endChr chr2 --startPosStartChr 112735986 --endPosStartChr 113204585 --startPosEndChr 112735986 --endPosEndChr 112735986 --res 5000
```

- [differential_interacted_region_analysis.Rmd](): This notebook illustrates how we use [getDIR.R]() and [annotateDIR.R]() to identify the differentially-interacted regions (DIRs) between the stimulated and the control conditions.


