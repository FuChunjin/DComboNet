# DComboNet

## Introduction

` DComboNet` is an R package for personalized anti-cancer drug combination prediction based on multi-level data integration. There are two main prediction models contained in the package. The level one model is for generalized anti-cancer drug combination effectiveness prediction and level two model is for cancer sample specific drug combination prediction. The two-level model based on a network-based method which integrates five subnetwork including drug-drug, drug-gene, drug-pathway, gene-gene and pathway-pathway association networks. Random walk with restart(RWR) algorithm is used to capture global proximity between drugs and the result is based on the rank returned via RWR algorithm. ` DComboNet` also provide clues for the potential mechanisms of drug combinations by extracting the top ranked genes/drugs between predicted drug combinations.

This tutorial provides the instruction of the main usage of this package fitting different scenario. This tutorial will lead to know the basic usage of the prediction, the description of prepared data how to extend the network construction to fit your own dataset. This tutorial will not present detail description of all functions contain in the packagem but you can easily learn those in help document with the R package. 

## Package installation

` DComboNet` has been upload in Github and can be install as follow:



```{r, eval=FALSE}
  install.packages("devtools")
  devtools::install_github(".../DComboNet")
  library(DComboNet)
```
  
## Environment


## Quick start  

To quickly get the usage of `DComboNet`, prediction assignment for one drug is taken as an example to go through the basic functions. Here, `Sorafenib` is taken as example, note that `Sorafenib` is within our pre-prepared data set, if you want to know how to do prediction for new drugs, please check session 5 for more details.


Before start any prediction assignment, please make sure `data` folder is fully downloaded and put in a fetchable path. `data` folder contains data files for construct networks as well as testing datasets used in our publication.


```{r, eval=TRUE}
# 1. Make sure the path of `data` folder assign to variable correctly
## Here is an example path, you should switch to your own
load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
# 2. Specify the path to save result files
## Here is an example path, you should switch to your own
resultdir = "G:/lab/DCcomboNet/Rpackage/tryout_result/"
```

As we introduced in session 3, `DComboNet()` function incorporate two level models, you can choose level one model which is a generalized anti-cancer drug combination model and level two model which is a cancer-sample specific drug combination prediction model.

### Level one model

To start generalized prediction assignment, you should choose `L1` in parameter `model`. Then drug that you are interested in predict combinable drugs should be input. Here, we provide two possible solutions to input the drug you are interested in `DComboNet`: 
    
1\) Input one or more drug name(s) together with their drugbank ID in a `data.frame`. In this case, you should shield manual input by `manual_input = FALSE`.


```{r, eval=FALSE}
# runing level one model (L1)
drugcandidate = data.frame(drug = "Sorafenib", drugbankID = "DB00398")

DComboNet(load_dir = load_dir, 
          resultdir = resultdir, 
          model = "L1", # To choose level one model
          manual_input = FALSE, # To shield manually input drug name
          drugcandidate = drugcandidate, 
          drugnetWeight = TRUE, # Confirm if drug network should be weighted
          featuretype = 'integrated_score') # Select which drug-drug similarity should be use

```

Note that if neither has `drugcandidate` been inputed nor has `manual_input` been allowed, all drug nodes in constructed network will be traversed for prediction.

2\) Submit the name of drug in an interactive manner. Drug for prediction can be entered interactively after the question 'Please type in the drug you are interested: ' pop out. Note that, if `manual_input = TRUE`, even with input `drugcandidate`, the final result will be generates for the manually input drug.



```{r, eval=FALSE}

DComboNet(load_dir = load_dir, 
          resultdir = resultdir, 
          model = "L1",
          manual_input = TRUE, # Select manually input function
          drugcandidate = NULL, 
          drugnetWeight = TRUE,
          featuretype = 'integrated_score')

# Don't run here
# Please type in the drug you are interested: 
# Sorafenib

```

Three results will be generated and saved under the folder `drugrank`, `generank` and `pathwayrank` in `resultdir` path. Prediction result for drug candidates is saved in `drugrank` folder. The name of result file indicate `drugseed` and inside the file, drug candidates within network is ranked based on global similarity.

```{r eval=FALSE, results='hide'}
# Take Sorafenib as an example
model = "L1"
drugrank = read.csv(paste0(resultdir,model,"_result/drugrank/Sorafenib_rank.csv"))
DT::datatable(drugrank, options = list(pageLength = 5))
```

Genes and pathways are also ranked according to global similarity. This two files will be used for subnetwork extraction and visualization.


```{r eval = FALSE, results='hide'}
# Take Sorafenib as an example
generank = read.csv(paste0(resultdir,model,"_result/generank/Sorafenib_rank.csv"))
DT::datatable(generank, options = list(pageLength = 5))
pathwayrank = read.csv(paste0(resultdir,model,"_result/pathwayrank/Sorafenib_rank.csv"))
DT::datatable(pathwayrank, options = list(pageLength = 5))
```

### Level two model


To run cancer sample specific model requires more input parameters, especially to define cancer cell line and drug treatment related variables. Drug induced gene expression changes are from LINCS database where two datasets (`GSE70138` and `GSE92742`) are provided from GEO database. Multiple cancer cell lines together with different treatment time are also included in the two datasets. Depends on the specific purpose, you can choose cell line, treatment time and filtering criteria to generate your own gene set.

Still take `Sorafenib` as example, if aiming at predicting combinable drug for `Sorafenib` in Hepatocellular carcinoma cell line `HEPG2`, since corresponding dataset for `HEPG2` is provided in LINCS and CCLE database, you can choose `dataset`, `cellline` and `treatment_time` for calling the drug-induced differential expressed gene table and HEPG2 specific expressed gene table, then the model can be set as follow:



```{r, eval=FALSE, warning=FALSE}
# Prepare input drugseed
drugcandidate = data.frame(drug = "Sorafenib", drugbankID = "DB00398")

# The setting of other parameters can be the same as how we run level one model
DComboNet(load_dir = load_dir, 
          resultdir = resultdir, 
          model = "L2", # Choose level two model
          manual_input = FALSE, # You can use manual input function (see level one model example)
          drugcandidate = drugcandidate, 
          drugnetWeight = TRUE,
          featuretype = 'integrated_score', 
          dataset = "92742",
          cellline = "HEPG2",
          treatment_time = 6,
          # the absolute value of fold-change between drug-treated group and control group will be above 0.5
          foldchange = 0.5,
          # the p-value from the significent test (t.test) between drug-treated group and control group will be below 0.05
          pvalue = 0.05)

```

If you have other gene list related to the cancer cell line you are interested in and want to included in drug-gene assocaition network and/or gene-gene association network, you can input drugDEG and cancergene table (in data.frame format) manually.

Similar to level one model, three results will be generated and saved under the folder `drugrank`, `generank` and `pathwayrank` in `resultdir` path. Prediction result for drug candidates is saved in `drugrank` folder. The name of result file indicate `drugseed` and inside the file, drug candidates within network is ranked based on global similarity.


```{r eval = FALSE, results='hide'}
# Take Sorafenib as an example
model = "L2"
drugrank = read.csv(paste0(resultdir,model,"_result/drugrank/Sorafenib_rank.csv"))
# DT::datatable(drugrank, options = list(pageLength = 5))
```

Genes and pathways are also ranked according to global similarity. This two files will be used for subnetwork extraction and visualization.


```{r eval = FALSE, results='hide'}
# Take Sorafenib as an example
generank = read.csv(paste0(resultdir,model,"_result/generank/Sorafenib_rank.csv"))
# DT::datatable(generank, options = list(pageLength = 5))
pathwayrank = read.csv(paste0(resultdir,model,"_result/pathwayrank/Sorafenib_rank.csv"))
# DT::datatable(pathwayrank, options = list(pageLength = 5))
```

### Network visualization

Other than prediction function, `DComboNet` package provides functions to extract and visualize subnetwork between user-interested drug seed and its corresponding predicted combinable drug candidate, which may help infer the possible mechanism of drug combinations. This step only works when prediction assignment finished. 

To run network visualization function, you need to check if the prediction result has been saved to given path `resultdir`. Note that the result files including drugrank table, generank table and pathwayrank table. Then capitalized drug names that you are interested can be inputed as `drugseed` and `drugcandidate`. `drugtarget` table should also be prepared as a dataframe and take as input. Genes and pathways for subnetwork construction are according to their rank corresponding to drugseed, the default value of gene rank (controlled via parametre `generank_filter`)is 0.01 meaning that only genes ranked on top 1% will be kept in subnetwork while default value of pathway rank(controlled via parametre `pathwayrank_filter`) is 0.1 meaning that only pathway ranked on top 10% will be kept in subnetwork. 

After preparation of network visualization input, `network_visualization()` is used for visualizing network. This function was built based on `visNetwork` package. Drugseed, drugcandidate and their targeted genes are colorred coded. User can either click to select interesting nodes and drag around for better layout or select by id or group. An `.html` file will be generated and saved in under the result path. Furthermore, a `.graphml` file will also be saved. This file can be easily import as a network in Cytoscape and network style can be shown via choosing **column** type (e.g. **color**) and **Mapping type** set as "Passthrough Mapping".


```{r eval = FALSE, warning=FALSE}
library(visNetwork)

drugseed = "Sorafenib"
drugcandidate = "Vorinostat"
drugtarget = data.frame(drug = c(rep(drugseed,10),rep(drugcandidate,5)),
                        target = c("BRAF", "FGFR1", "FLT1", "FLT3", "FLT4", "KDR", "KIT", "PDGFRB", "RAF1", "RET", "HDAC1", "HDAC2", "HDAC3", "HDAC6", "HDAC8"))

model = "L2"

network_extract(drugseed = drugseed, 
                drugcandidate = drugcandidate, 
                drugtarget = drugtarget, 
                generank_filter = 0.01, 
                pathwayrank_filter = 0.1,
                model = model,
                cellline = "HEPG2",
                load_dir = load_dir, 
                resultdir = resultdir)

network_visualization(drugseed = drugseed, 
                      drugcandidate =  drugcandidate, 
                      drugtarget = drugtarget, 
                      model = model, 
                      cellline = "HEPG2",
                      load_dir = load_dir,
                      resultdir = resultdir)
```


## Citation
<!-- <div style = "width:120%; height:auto; margin: auto;"> -->

<!-- <p style="text-indent:16px;">If you use `DComboNet` in your publication(s), please cite:</p> -->
<!-- </div> -->

## Acknowledgement

