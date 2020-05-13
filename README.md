# DComboNet

## Introduction

` DComboNet` is an R package for personalized anti-cancer drug combination prediction based on multi-level data integration. There are two main prediction models contained in the package. The level one model is for generalized anti-cancer drug combination effectiveness prediction and level two model is for cancer sample specific drug combination prediction. The two-level model based on a network-based method which integrates five subnetwork including drug-drug, drug-gene, drug-pathway, gene-gene and pathway-pathway association networks. Random walk with restart(RWR) algorithm is used to capture global proximity between drugs and the result is based on the rank returned via RWR algorithm. ` DComboNet` also provide clues for the potential mechanisms of drug combinations by extracting the top ranked genes/drugs between predicted drug combinations.

This tutorial provides the instruction of the main usage of this package fitting different scenario. This tutorial will lead to know the basic usage of the prediction, the description of prepared data how to extend the network construction to fit your own dataset. This tutorial will not present detail description of all functions contain in the packagem but you can easily learn those in help document with the R package. 

## Package installation

` DComboNet` has been upload in Github and can be install in R as follow:

```{r, eval=FALSE}
  install.packages("devtools")
  devtools::install_github("VeronicaFung/DComboNet")

  library(DComboNet)
```

Due to the large size of drug-induced gene expression profiles for Level 2 model, we upload only cell line data for repeating the results in our manuscripts, for more cell lines, we prepare an extra data folder in [links]<links>. 

After install the package, you should download this compressed data set, unzip it and put in a reachable path. This path should be used as input for parameter ` load_dir` for functions in ` DComboNet`. For more details about using level two model with provided cancer sample specific data or your own data table, you can check the tutorial below.

## Tutorial  

After installing package, you can use sample code in folder "test" for  testing.
The detail usages of ` DComboNet`, can be found here: [Instruction](https://veronicafung.github.io/DComboNet/DComboNet-vignette.html) 


## Environment

` DComboNet` is developed in R version 3.5.1. For running the test code in Instruction, you can check the "sessionInfo" chapter for more detailed session information.

## Citation
<!-- <div style = "width:120%; height:auto; margin: auto;"> -->

<!-- <p style="text-indent:16px;">If you use `DComboNet` in your publication(s), please cite:</p> -->
<!-- </div> -->

## Acknowledgement
