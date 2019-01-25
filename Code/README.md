# DreamAI
- [DreamAI::DreamAI](#dreamaidreamai)
   - Imputation of Missing Protein Abundances with Iterative Prediction Model
- [DreamAI::DreamAI_Bagging](#dreamaidreamai_bagging)
   - Bag Imputation of Missing Protein Abundances with Iterative Prediction Model
- [DreamAI::bag.summary](#dreamaibagsummary)
   - Wrapper function for summarizing the outputs from DreamAI_bagging

## DreamAI::DreamAI
- [Description](#description)
- [Usage](#usage)
- [Arguments](#arguments)
- [Value](#value)
- [Notes](#notes)
- [Example](#example)

### Description

The function DreamAI imputes a dataset with missing values or NA's using 7 different methods: 

 - "KNN": k nearest neighbor 
 - "MissForest": nonparametric Missing Value
   Imputation using Random Forest 
  - "ADMIN": abundance dependent missing
   imputation
   - "Brinn": imputation using IRNN-SCAD algorithm 
   - "SpectroFM": imputation using matrix factorization 
   -  "RegImpute": imputation using Glmnet ridge regression  
   -  "Ensemble": aggregation of the 6 methods
   using simple average.

### Usage
```
DreamAI(data, k = 10, maxiter_MF = 10, ntree = 100,
  maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
  gamma_ADMIN = NA, gamma = 50, CV = FALSE,
  fillmethod = "row_mean", maxiter_RegImpute = 10,
  conv_nrmse = 1e-06, iter_SpectroFM = 40, method = c("KNN",
  "MissForest", "ADMIN", "Brinn", "SpectroFM", "RegImpute"),
  out = c("KNN", "MissForest", "ADMIN", "Brinn", "SpectroFM", "RegImpute""Ensemble"))
```
### Arguments
  
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
| out		      | c("KNN", "MissForest", "ADMIN", "Brinn", "SpectroFM", "RegImpute""Ensemble")     	   | a vector of imputation methods for which the function will output the imputed matrices

	
### Value
a list of imputed datasets by different methods as specified by the user. Always returns imputed data by "Ensemble"

### Notes
If all methods are specified for obtaining "Ensemble" imputed matrix, the approximate time required to output the imputed matrix for a dataset of dimension 26000 x 200 is ~50 hours.

### Example
```
data(datapnnl)
data<-datapnnl.rm.ref[1:100,1:21]
impute<- DreamAI(data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40, method = c("KNN", "MissForest", "ADMIN", "Brinn", "SpectroFM", "RegImpute"),out="Ensemble")
impute$Ensemble
```

## DreamAI::DreamAI_Bagging
- [Description](#description-1)
- [Usage](#usage-1)
- [Arguments](#arguments-1)
- [Value](#value-1)
- [Notes](#notes-1)
- [Example](#example-1)

### Description

The function DreamAI_bagging imputes a dataset with missing values or NA's by bag imputaion with help of parallel processing. Pseudo datasets are generated having true missing (as in the original dataset) and pseudo missing and every such pseudo dataset is imputed by 7 different methods: KNN, MissForest, ADMIN, Brinn, SpectroFM, RegImpute and Ensemble (descriptions are included in the documentation of the function DreamAI).

### Usage
```
DreamAI_Bagging(data, k = 10, maxiter_MF = 10, ntree = 100,
  maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
  gamma_ADMIN = NA, gamma = 50, CV = FALSE,
  fillmethod = "row_mean", maxiter_RegImpute = 10,
  conv_nrmse = 1e-06, iter_SpectroFM = 40, method = c("KNN",
  "MissForest", "ADMIN", "Brinn", "SpectroFM", "RegImpute", "Ensemble"),
  SamplesPerBatch, n.bag, save.out = TRUE, path = NULL, ProcessNum)
```
### Arguments
  
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
| SamplesPerBatch			             |      	     | convergence threshold for "RegImpute"
| n.bag		    |      	     | number of pseudo datasets to generate and impute in the current process save.out logical indicator whether or not to save the output. When TRUE output is saved, when FALSE output is returned
| path		      | NULL	   | location to save the output file from the curent process. Path only needs to be specified when save.out=TRUE
| ProcessNum		      |     	   | process number starting from 1 when run in cluster, e.g. 1 - 10, 1 - 100 etc. Needs to be specified only if the output is saved
| out		      | c("Ensemble")     	   | a vector of imputation methods for which the function will output the imputed matrices

	
### Value
list of imputed dataset, n.bag and a matrix containing gene name, sample name, true and imputed values of every pseudo missing combined from n.bag datasets. Impute is a list of imputed datasets (average over all pseudo imputed data matrices) by different methods as specified by the user. Always returns imputed data by "Ensemble"

### Notes
This function can be run as parallel job in cluster. It generates and saves a .RData file containing the output from the current process in the location provided by the user, with the process number in the file name. If the user runs it in local computer multiple times, then changing the ProcessNumber everytime will generate and save .RData file with the given ProcessNumber.

### Example
```
data(datapnnl)
data<-datapnnl.rm.ref[1:100,1:21]
impute<- DreamAI_Bagging(data=data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40,method=c("KNN","MissForest","ADMIN","Brinn","SpectroFM","RegImpute","Ensemble"),SamplesPerBatch=3,n.bag=2,save.out=TRUE,path="C:\\Users\\chowds14\\Desktop\\test_package\\",ProcessNum=1)
impute$Ensemble
```

## DreamAI::bag.summary
- [Description](#description-2)
- [Usage](#usage-2)
- [Arguments](#arguments-2)
- [Value](#value-2)
- [Example](#example-2)

### Description

Wrapper function for summarizing the outputs from DreamAI_bagging

### Usage
```
bag.summary(method = c("KNN", "MissForest", "ADMIN", "Brinn",
  "SpectroFM", "RegImpute", "Ensemble"), nNodes = 2, path = NULL)
```
### Arguments
  
| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| method	       |c("KNN", "MissForest", "ADMIN", "Brinn",
  "SpectroFM", "RegImpute", "Ensemble")	           |a vector of imputation methods. Default is "Ensemble" if nothing is specified
| nNodes         | 2           |number of parallel processes
| path 	       |NULL	            |location where the bagging output is saved
	
### Value
list of final imputed data and confidence score for every gene using pseudo missing

### Example
```
data(datapnnl)
data<-datapnnl.rm.ref[1:100,1:21]
impute<- DreamAI_Bagging(data=data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40,method=c("KNN","MissForest","ADMIN","Brinn","SpectroFM","RegImpute","Ensemble"),SamplesPerBatch=3,n.bag=2,save.out=TRUE,path="C:\\Users\\chowds14\\Desktop\\test_package\\",ProcessNum=1)
final.out<-bag.summary(method=c("Ensemble"),nNodes=2,path="C:\\Users\\chowds14\\Desktop\\test_package\\")
final.out$score
final.out$imputed_data
```

