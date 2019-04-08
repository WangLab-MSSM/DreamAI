# DreamAI
- [Authors](#authors)
- [Overview](#Overview)
- [Installation](#installation)
- [Usage](#Usage)
- [Arguments](#Arguments)
- [Value](#Value)
- [Note](#Note)
- [Examples](#Examples)
- [Contributions](#contributions)

## Authors
Shrabanti Chowdhury<sup>1</sup>, Weiping Ma<sup>1</sup>, Sunkyu Kim<sup>2</sup>; Zhi Li<sup>3</sup>, Thomas Yu<sup>4</sup>, Mi Yang<sup>5,6</sup>, Francesca Petralia<sup>1</sup>, Jeremy Jacobsen<sup>7</sup>, Jingyi Jessica Li<sup>8</sup>, Xinzhou Ge<sup>8</sup>, Kexin Li<sup>9</sup>, Nathan Edwards<sup>10</sup>, Samuel Payne<sup>11</sup>, Henry Rodriguez<sup>12</sup>, Paul Boutros<sup>13</sup>, Gustavo Stolovitzky<sup>14</sup>, Jaewoo Kang<sup>2</sup>, David Fenyo<sup>3</sup>, Julio Saez-Rodriguez,<sup>6,15</sup>, Pei Wang<sup>1</sup>

<sup>1</sup>Icahn School of Medicine at Mount Sinai (USA), <sup>2</sup>Department of Computer Science and Engineering, Korea University (South Korea), <sup>3</sup>New York University (USA), <sup>4</sup>Sage Bionetworks (USA), <sup>5</sup>Heidelberg University, Faculty of Biosciences (Germany), <sup>6</sup>RWTH Aachen University (Germany), <sup>7</sup>University of Colorado (USA), <sup>8</sup>Department of Statistics, University of California (USA), <sup>9</sup>Department of Mathematics, Tsinghua University (China), <sup>10</sup>Georgetown University (USA), <sup>11</sup>Pacific Northwest National Laboratory (USA), <sup>12</sup>National Cancer Institute (USA), <sup>13</sup>Ontario Institute of Cancer Research (Canada), <sup>14</sup>IBM Research & Mount Sinai (USA), <sup>15</sup>European Molecular Biology Laboratory-European Bioinformatics Institute (UK)

## Overview

To develop powerful computational tools to extract the most information from the proteome, Clinical Proteomic Tumor Analysis Consortium (CPTAC) and DREAM organization launched The NCI-CPTAC DREAM Proteogenomics Challenge in 2016, one of the subchallange: impute missing values in proteomics data given observed proteins.

In this challenge, participants were invited to develop proper imputation algorithms for proteomics data. And with their help an optimal imputation method: DreamAI was ensembled as an outcome of this challenge.

Specifically in DreamAI, ensemble imputation matrix is obtained from averaging results of six imputation algorithms: top 3 teams in challenge (spectroFM: Team DMIS_PTG; RegImpute: Team Jeremy Jacobsen; Birnn: Team BruinGo) and 3 baseline algorithms (KNN, missForest, ADMIN). Bootstrap aggregating (bagging) is also adopted to improve unstable estimation and accuracy of machine learning algorithms.

In the output option of this function, it provides user the flexibility to select imputation matrix from the ensemble method or each individual algorithm:
 - "KNN": k nearest neighbor imputation
 - "MissForest": nonparametric Missing Value Imputation using Random Forest 
  - "ADMIN": abundance dependent missing imputation
   - "Brinn": imputation using IRNN-SCAD algorithm 
   - "SpectroFM": imputation using matrix factorization 
   -  "RegImpute": imputation using Glmnet ridge regression  
   -  "Ensemble": average of the 6 methods



## Installation

Packages required prior to installing DreamAI
```
require("cluster")
require("survival")
require("randomForest")
require("missForest")
require("glmnet")
require("foreach")
require("itertools")
require("iterators")
require("Matrix")
require("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("impute", version = "3.8")
require("impute")

```

Install DreamAI
```
install_github("WangLab-MSSM/DreamAI/Code")
```


## Usage
```
DreamAI(data, k = 10, maxiter_MF = 10, ntree = 100,
  maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
  gamma_ADMIN = NA, gamma = 50, CV = FALSE,
  fillmethod = "row_mean", maxiter_RegImpute = 10,
  conv_nrmse = 1e-06, iter_SpectroFM = 40, method = c("KNN",
  "MissForest", "ADMIN", "Brinn", "SpectroFM", "RegImpute"),
  out = c("Ensemble"))
```
## Arguments
  
| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| data	       |	           |dataset in the form of a matrix or dataframe with missing values or NA's. The function throws an error message and stops if any row or column in the dataset is missing all values
| k         | 10           |number of neighbors to be used in the imputation by KNN and ADMIN
| maxiter_MF 	       |	10	            |maximum number of iteration to be performed in the imputation by "MissForest" if the stopping criteria is not met beforehand
| ntree  		       | 100	           | number of trees to grow in each forest in "MissForest"
| maxnodes		           | NULL             | maximum number of terminal nodes for trees in the forest in "MissForest", has to equal at least the number of columns in the given data
| maxiter_ADMIN 	        | 30           | maximum number of iteration to be performed in the imputation by "ADMIN" if the stopping criteria is not met beforehand
| tol	         | 10^(-2)             | convergence threshold for "ADMIN"
| gamma_ADMIN          | NA           | parameter for ADMIN to control abundance dependent missing. Set gamma_ADMIN=0 for log ratio intensity data. For abundance data put gamma_ADMIN=NA, and it will be estimated accordingly
| gamma       | 50  | parameter of the supergradients of popular nonconvex surrogate functions, e.g. SCAD and MCP of L0-norm for Brinn
| CV   | FALSE         | a logical value indicating whether to fit the best gamma with cross validation for "Brinn". If CV=FALSE, default gamma=50 is used, while if CV=TRUE gamma is calculated using cross-validation.
| fillmethod			             | "row_mean" 	           | a string identifying the method to be used to initially filling the missing values using simple imputation for "RegImpute". That could be "row_mean" or "zeros", with "row_mean" being the default. It throws an warning if "row_median" is used.
| maxiter_RegImpute			     | 10         | maximum number of iterations to reach convergence in the imputation by "RegImpute"
| conv_nrmse			             | 1e-06     	     | convergence threshold for "RegImpute"
| iter_SpectroFM		    | 40     	     | number of iterations for "SpectroFM"
| method		      | c("KNN","MissForest", "ADMIN", "Brinn", "SpectroFM", "RegImpute", "Ensemble")     	   | a vector of imputation methods: ("KNN", "MissForest", "ADMIN", "Brinn", "SpectroFM, "RegImpute", "Ensemble"). Default is "Ensemble" if nothing is specified
| out		      | c("Ensemble")     	   | a vector of imputation methods for which the function will output the imputed matrices

	
## Value
a list of imputed datasets by different methods as specified by the user. Always returns imputed data by "Ensemble"

## Note
If all methods are specified for obtaining "Ensemble" imputed matrix, the approximate time required to output the imputed matrix for a dataset of dimension 26000 x 200 is ~50 hours.

## Examples
```
data(datapnnl)
data<-datapnnl.rm.ref[1:100,1:21]
impute<- DreamAI(data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40, method = c("KNN", "MissForest", "ADMIN", "Brinn", "SpectroFM", "RegImpute"),out="Ensemble")
impute$Ensemble
```

## Contributions

If you find small bugs, larger issues, or have suggestions, please email the maintainer at <shrabanti.chowdhury@mssm.edu> or <weiping.ma@mssm.edu>. Contributions (via pull requests or otherwise) are welcome.
