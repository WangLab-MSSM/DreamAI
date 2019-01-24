#' Imputation of Missing Protein Abundances with Iterative Prediction Model
#' @description The function DreamAI imputes a dataset with missing values or NA's using 7 different methods: "KNN": k nearest neighbor, "MissForest": nonparametric Missing Value Imputation using Random Forest, "ADMIN": abundance dependent missing imputation, "Brinn": imputation using IRNN-SCAD algorithm, "SpectroFM": imputation using matrix factorization, "RegImpute": imputation using Glmnet ridge regression and "Ensemble": aggregation of the 6 methods using simple average.
#' @param data dataset in the form of a matrix or dataframe with missing values or NA's. The function throws an error message and stops if any row or column in the dataset is missing all values
#' @param k number of neighbors to be used in the imputation by KNN and ADMIN (default is 10)
#' @param maxiter_MF maximum number of iteration to be performed in the imputation by "MissForest" if the stopping criteria is not met beforehand
#' @param ntree number of trees to grow in each forest in "MissForest"
#' @param maxnodes maximum number of terminal nodes for trees in the forest in "MissForest", has to equal at least the number of columns in the given data
#' @param maxiter_ADMIN maximum number of iteration to be performed in the imputation by "ADMIN" if the stopping criteria is not met beforehand
#' @param tol convergence threshold for "ADMIN"
#' @param gamma_ADMIN parameter for ADMIN to control abundance dependent missing. Default is 0 (for log ratio intensity data). For abundance data put gamma=NA, and it will be estimated accordingly 
#' @param gamma parameter of the supergradients of popular nonconvex surrogate functions, e.g. SCAD and MCP of L0-norm for Brinn
#' @param CV a logical value indicating whether to fit the best gamma with cross validation for "Brinn". If CV=FALSE, default gamma=50 is used, while if CV=TRUE gamma is calculated using cross-validation.
#' @param fillmethod a string identifying the method to be used to initially filling the missing values using simple imputation for "RegImpute". That could be "row_mean" or "zeros", with "row_mean" being the default. It throws an warning if "row_median" is used.
#' @param maxiter_RegImpute maximum number of iterations to reach convergence in the imputation by "RegImpute"
#' @param conv_nrmse convergence threshold for "RegImpute"
#' @param iter_SpectroFM number of iterations for "SpectroFM"
#' @param method a vector of imputation methods: ("KNN", "MissForest", "ADMIN", "Brinn", "SpectroFM, "RegImpute", "Ensemble"). Default is "Ensemble" if nothing is specified
#'
#' @return a list of imputed datasets by different methods as specified by the user. Always returns imputed data by "Ensemble"
#' @export
#' @examples
#' \dontrun{
#' data(datapnnl)
#' data<-datapnnl.rm.ref[1:100,1:21]
#' impute<- DreamAI(data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40, method = c("Ensemble"))
#' impute$Ensemble
#' }
DreamAI<-function(data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40,method=c("KNN","MissForest","ADMIN","Brinn","SpectroFM","RegImpute","Ensemble"))
{
  missing_rows = (which(rowSums(is.na(data))==dim(data)[2]))
  if(length(missing_rows)>0){
    stop("Some Rows Are Missing All Values !!\n")
  }
  missing_cols = (which(colSums(is.na(data))==dim(data)[1]))
  if(length(missing_cols)>0){
    stop("Some Columns Are Missing All Values !!\n")
  }
  options(warn=-1)
  ### subsetting data without rows/cols with all missing ###

  # data = subset(data, (rowSums(is.na(data))!=dim(data)[2]))
  # data = subset(data, (colSums(is.na(data))!=dim(data)[1]))

  ### method of imputation ###
  methods<-c("KNN","MissForest","ADMIN","Brinn","SpectroFM","RegImpute")
  
  ## KNN ##
  sink("NULL")
  d.impute.knn = impute.KNN(data = as.matrix(data),k = k)
  sink()
  print("Method 1/7 complete")
  if(is.null(maxnodes))
  {
    maxnodes<-ncol(as.matrix(data))
  }
  sink("NULL")

  # ## MF ##
  d.impute.MF = impute.MF(data = as.matrix(data),maxiter_MF=maxiter_MF,ntree=ntree,maxnodes=maxnodes)
  sink()
  print("Method 2/7 complete")

  # ## ADMIN ##
  sink("NULL")
  d.impute.ADMIN = impute.ADMIN(data = as.matrix(data),k = k, gamma=gamma_ADMIN,maxiter_ADMIN = maxiter_ADMIN, tol = tol)
  sink()
  print("Method 3/7 complete")

  # ## BruinGo ##
  sink("NULL")
  d.impute.BruinGo=impute.Brinn(as.matrix(data), gamma = 50, CV = CV)
  sink()
  print("Method 4/7 complete")

  ## DMIS ##
  sink("NULL")
  d.impute.DMIS=impute.SpectroFM(input_table=as.data.frame(data),iter=iter_SpectroFM, verbose=FALSE)
  sink()
  print("Method 5/7 complete")

  ## Jeremy ##
  sink("NULL")
  d.impute.Jeremy=impute.RegImpute(data=as.matrix(data),fillmethod=fillmethod,maxiter_RegImpute = maxiter_RegImpute,conv_nrmse = conv_nrmse)
  sink()
  print("Method 6/7 complete")

  ## Ensemble ##
  sink("NULL")
  ensemble<- (d.impute.knn+d.impute.MF+d.impute.ADMIN+d.impute.BruinGo+d.impute.DMIS+
                d.impute.Jeremy)/length(methods)
  sink()
  print("Method 7/7 complete")
  imputed_matrix=list("KNN"=as.matrix(d.impute.knn),"MissForest"=as.matrix(d.impute.MF),"ADMIN"=as.matrix(d.impute.ADMIN),"Brinn"=as.matrix(d.impute.BruinGo),"SpectroFM"=as.matrix(d.impute.DMIS),"RegImpute"=as.matrix(d.impute.Jeremy),"Ensemble"=as.matrix(ensemble))

  num<-which(methods %in% method)

  output<-imputed_matrix[c(num,7)]

  return(output)
}

