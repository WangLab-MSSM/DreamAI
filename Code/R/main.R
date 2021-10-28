#' Imputation of Missing Protein Abundances with Iterative Prediction Model
#' @description The function DreamAI imputes a dataset with missing values or NA's using 7 different methods: "KNN": k nearest neighbor, "MissForest": nonparametric Missing Value Imputation using Random Forest, "ADMIN": abundance dependent missing imputation, "Birnn": imputation using IRNN-SCAD algorithm, "SpectroFM": imputation using matrix factorization, "RegImpute": imputation using Glmnet ridge regression and "Ensemble": aggregation of the 6 methods using simple average.
#'
#' @param data dataset in the form of a matrix or dataframe with missing values or NA's. The function throws an error message and stops if any row or column in the dataset is missing all values
#' @param k number of neighbors to be used in the imputation by KNN and ADMIN (default is 10)
#' @param maxiter_MF maximum number of iteration to be performed in the imputation by "MissForest" if the stopping criteria is not met beforehand
#' @param ntree number of trees to grow in each forest in "MissForest"
#' @param maxnodes maximum number of terminal nodes for trees in the forest in "MissForest", has to equal at least the number of columns in the given data
#' @param maxiter_ADMIN maximum number of iteration to be performed in the imputation by "ADMIN" if the stopping criteria is not met beforehand
#' @param tol convergence threshold for "ADMIN"
#' @param gamma_ADMIN parameter for ADMIN to control abundance dependent missing. Set gamma_ADMIN=0 for log ratio intensity data. For abundance data put gamma_ADMIN=NA, and it will be estimated accordingly
#' @param gamma parameter of the supergradients of popular nonconvex surrogate functions, e.g. SCAD and MCP of L0-norm for Birnn
#' @param CV a logical value indicating whether to fit the best gamma with cross validation for "Birnn". If CV=FALSE, default gamma=50 is used, while if CV=TRUE gamma is calculated using cross-validation.
#' @param fillmethod a string identifying the method to be used to initially filling the missing values using simple imputation for "RegImpute". That could be "row_mean" or "zeros", with "row_mean" being the default. It throws an warning if "row_median" is used.
#' @param maxiter_RegImpute maximum number of iterations to reach convergence in the imputation by "RegImpute"
#' @param conv_nrmse convergence threshold for "RegImpute"
#' @param iter_SpectroFM number of iterations for "SpectroFM"
#' @param method a vector of imputation methods: ("KNN", "MissForest", "ADMIN", "Birnn", "SpectroFM, "RegImpute") based on which "Ensemble" imputed matrix will be obtained.  
#' @param out a vector of imputation methods for which the function will output the imputed matrices. Default is "Ensemble".  
#' 
#'@note If all methods are specified for obtaining "Ensemble" imputed matrix, the approximate time required to output the imputed matrix for a dataset of dimension 26000 x 200 is ~50 hours.
#' 
#' @return a list of imputed datasets by different methods as specified by the user. Always returns imputed data by "Ensemble"
#' @export
#' @examples
#' \dontrun{
#' data(datapnnl)
#' data<-datapnnl.rm.ref[1:100,1:21]
#' impute<- DreamAI(data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40, method = c("KNN", "MissForest", "ADMIN", "Birnn", "SpectroFM", "RegImpute"),out="Ensemble")
#' impute$Ensemble
#' }
DreamAI<-function(data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40,method=c("KNN","MissForest","ADMIN","Birnn","SpectroFM","RegImpute"),out=c("Ensemble"))
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
  methods<-c("KNN","MissForest","ADMIN","Birnn","SpectroFM","RegImpute")
  
  methods.match<- methods[which(methods %in% method)]
  
  if(length(methods.match)==0)
  {
    sink()
    return(print("specify method"))
  }
    
  ensemble<-matrix(0,nrow(data),ncol(data))
  method.idx<-1
  imputed_matrix=list()
  out_matrix=list()

  ## KNN ##
  
  if("KNN" %in% method)
  {
    sink("NULL")
    d.impute.knn = impute.KNN(data = as.matrix(data),k = k)
    ensemble<-ensemble+d.impute.knn
    sink()
    print(paste("Method",method.idx,"complete"))
    method.idx<-method.idx+1
    imputed_matrix<-c(imputed_matrix,list("KNN"=as.matrix(d.impute.knn)))
  }
  
  
  # ## MF ##
  if("MissForest" %in% method)
  {
    if(is.null(maxnodes))
    {
      maxnodes<-ncol(as.matrix(data))
    }
    sink("NULL")
    d.impute.MF = impute.MF(data = as.matrix(data),maxiter_MF=maxiter_MF,ntree=ntree,maxnodes=maxnodes)
    ensemble<-ensemble+d.impute.MF
    sink()
    print(paste("Method",method.idx,"complete"))
    method.idx<-method.idx+1
    imputed_matrix<-c(imputed_matrix,list("MissForest"=as.matrix(d.impute.MF)))
  }

  # ## ADMIN ##
  if("ADMIN" %in% method)
  {
    sink("NULL")
    d.impute.ADMIN = impute.ADMIN(data = as.matrix(data),k = k, gamma=gamma_ADMIN,maxiter_ADMIN = maxiter_ADMIN, tol = tol)
    ensemble<-ensemble+d.impute.ADMIN
    sink()
    print(paste("Method",method.idx,"complete"))
    method.idx<-method.idx+1
    imputed_matrix<-c(imputed_matrix,list("ADMIN"=as.matrix(d.impute.ADMIN)))
  }

  
  # ## Birnn ##
  if("Birnn" %in% method)
  {
    sink("NULL")
    d.impute.Birnn=impute.Birnn(as.matrix(data), gamma = 50, CV = CV)
    ensemble<-ensemble+d.impute.Birnn
    sink()
    print(paste("Method",method.idx,"complete"))
    method.idx<-method.idx+1
    imputed_matrix<-c(imputed_matrix,list("Birnn"=as.matrix(d.impute.Birnn)))
  }

  ## SpectroFM ##
  if("SpectroFM" %in% method)
  {
    sink("NULL")
    d.impute.SpectroFM=impute.SpectroFM(input_table=as.data.frame(data),iter=iter_SpectroFM, verbose=FALSE)
    ensemble<-ensemble+d.impute.SpectroFM
    sink()
    print(paste("Method",method.idx,"complete"))
    method.idx<-method.idx+1
    imputed_matrix<-c(imputed_matrix,list("SpectroFM"=as.matrix(d.impute.SpectroFM)))
  }


  ## RegImpute ##
  if("RegImpute" %in% method)
  {
    sink("NULL")
    d.impute.RegImpute=impute.RegImpute(data=as.matrix(data),fillmethod=fillmethod,maxiter_RegImpute = maxiter_RegImpute,conv_nrmse = conv_nrmse)
    ensemble<-ensemble+d.impute.RegImpute
    sink()
    print(paste("Method",method.idx,"complete"))
    method.idx<-method.idx+1
    imputed_matrix<-c(imputed_matrix,list("RegImpute"=as.matrix(d.impute.RegImpute)))
  }


  ## Ensemble ##
  
    ensemble<-ensemble/length(method)
    imputed_matrix<-c(imputed_matrix,list("Ensemble"=as.matrix(ensemble)))
    
    
    if("KNN" %in% out)
    {
      if("KNN" %in% method){
        out_matrix<-c(out_matrix,list("KNN"=as.matrix(d.impute.knn)))}else{
          print("KNN is not specified")
        }
    }
    
    if("MissForest" %in% out)
    {
      if("MissForest" %in% method){
        out_matrix<-c(out_matrix,list("MissForest"=as.matrix(d.impute.MF)))}else{
          print("MissForest is not specified")
        }
    }
    
    if("ADMIN" %in% out)
    {
      if("ADMIN" %in% method){
        out_matrix<-c(out_matrix,list("ADMIN"=as.matrix(d.impute.ADMIN)))}else{
          print("ADMIN is not specified")
        }
    }
    
    if("Birnn" %in% out)
    {
      if("Birnn" %in% method){
        out_matrix<-c(out_matrix,list("Birnn"=as.matrix(d.impute.Birnn)))}else{
          print("Birnn is not specified")
        }
    }  
    
    if("SpectroFM" %in% out)
    {
      if("SpectroFM" %in% method){
        out_matrix<-c(out_matrix,list("SpectroFM"=as.matrix(d.impute.SpectroFM)))}else{
          print("SpectroFM is not specified")
        }
    } 
    
    if("RegImpute" %in% out)
    {
      if("RegImpute" %in% method){
        out_matrix<-c(out_matrix,list("RegImpute"=as.matrix(d.impute.RegImpute)))}else{
          print("RegImpute is not specified")
        }
    } 
    
    out_matrix<-c(out_matrix,list("Ensemble"=as.matrix(ensemble)))
    
    sink()
    # print(paste("Method",method.idx,"complete"))
  

  # sink("NULL")
  # ensemble<- (d.impute.knn+d.impute.MF+d.impute.ADMIN+d.impute.Birnn+d.impute.SpectroFM+d.impute.RegImpute)/length(methods)
  # sink()
  # 
  # imputed_matrix=list("KNN"=as.matrix(d.impute.knn),"MissForest"=as.matrix(d.impute.MF),"ADMIN"=as.matrix(d.impute.ADMIN),"Birnn"=as.matrix(d.impute.Birnn),"SpectroFM"=as.matrix(d.impute.SpectroFM),"RegImpute"=as.matrix(d.impute.RegImpute),"Ensemble"=as.matrix(ensemble))
  # 
  #num<-which(method %in% out)

  output<-out_matrix

  return(output)
}

