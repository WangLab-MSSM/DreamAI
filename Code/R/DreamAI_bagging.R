gamma.est = function(Ym)
{
  ### Ym dim Lxn
  Y.mean = apply(Ym,1,mean,na.rm=T);
  q = apply(is.na(Ym),1,mean);
  par = lsfit(Y.mean[q>0],log(1/q)[q>0])[[1]];
  names(par) = c('gamma0','gamma');
  return(par);
}

avg.batch = function(data,SamplesPerBatch)
{
  t(apply(data,1,function(x){apply(matrix(x,SamplesPerBatch),2,mean,na.rm=T)}));
}

diff.mr = function(x,gm.pnnl,data.pnnl.b,N){
  p.temp = exp(x-gm.pnnl[2]*data.pnnl.b);
  p.temp[is.na(p.temp)] = 1;
  return(abs(mean(p.temp)-(sum(is.na(data.pnnl.b))+N)/sum(is.na(data.pnnl.b))*mean(is.na(data.pnnl.b))));
}



#' Bag Imputation of Missing Protein Abundances with Iterative Prediction Model
#' @description The function DreamAI_bagging imputes a dataset with missing values or NA's by bag imputaion with help of parallel processing. Pseudo datasets are generated having true missing (as in the original dataset) and pseudo missing and every such pseudo dataset is imputed by 7 different methods: KNN, MissForest, ADMIN, Birnn, SpectroFM, RegImpute and Ensemble (descriptions are included in the documentation of the function DreamAI). 
#' @details This function can be run as parallel job in cluster. It generates and saves a .RData file containing the output from the current process in the location provided by the user, with the process number in the file name. If the user runs it in local computer multiple times, then changing the ProcessNumber everytime will generate and save .RData file with the given ProcessNumber. 
#'
#' @param data dataset in the form of a matrix or dataframe with missing values or NA's. The function throws an error message and stops if any row or column in the dataset is missing all values
#' @param k number of neighbors to be used in the imputation by KNN and ADMIN (default is 10)
#' @param maxiter_MF maximum number of iteration to be performed in the imputation by "MissForest" if the stopping criteria is not met beforehand
#' @param ntree number of trees to grow in each forest in "MissForest"
#' @param maxnodes maximum number of terminal nodes for trees in the forest in "MissForest", has to equal at least the number of columns in the given data
#' @param maxiter_ADMIN maximum number of iteration to be performed in the imputation by "ADMIN" if the stopping criteria is not met beforehand
#' @param tol convergence threshold for "ADMIN"
#' @param gamma_bagging parameter to control abundance-dependence for pseudo missing values generation in bagging sets. Set gamma_bagging = NA, to learn the abundance-dependence from the observed data matrix. Set gamma_bagging = 0 to generate abundance-independent pseudo missing values.
#' @param gamma_ADMIN parameter for ADMIN to control abundance dependent missing. Set gamma_ADMIN=0 for log ratio intensity data. For abundance data put gamma_ADMIN=NA, and it will be estimated accordingly
#' @param gamma parameter of the supergradients of popular nonconvex surrogate functions, e.g. SCAD and MCP of L0-norm for Birnn
#' @param CV a logical value indicating whether to fit the best gamma with cross validation for "Birnn". If CV=FALSE, default gamma=50 is used, while if CV=TRUE gamma is calculated using cross-validation.
#' @param fillmethod a string identifying the method to be used to initially filling the missing values using simple imputation for "RegImpute". That could be "row_mean" or "zeros", with "row_mean" being the default. It throws an warning if "row_median" is used.
#' @param maxiter_RegImpute maximum number of iterations to reach convergence in the imputation by "RegImpute"
#' @param conv_nrmse convergence threshold for "RegImpute"
#' @param iter_SpectroFM number of iterations for "SpectroFM"
#' @param method a vector of imputation methods: ("KNN", "MissForest", "ADMIN", "Birnn", "SpectroFM, "RegImpute") based on which "Ensemble" imputed matrix will be obtained
#' @param out a vector of imputation methods for which the function will output the imputed matrices. Default is "Ensemble".  
#' @param SamplesPerBatch number of samples per batch (batch size in the original data)
#' @param n.bag number of pseudo datasets to generate and impute in the current process
#' @param path location to save the output file from the curent process.  Path only needs to be specified when save.out=TRUE
#' @param ProcessNum process number starting from 1 when run in cluster, e.g. 1 - 10, 1 - 100 etc. Needs to be specified only if the output is saved
#' @param save.out logical indicator whether or not to save the output. When TRUE output is saved, when FALSE output is only returned
#'
#' @return list of imputed dataset(average over all pseudo imputed data matrices) by different methods as specified by the user, n.bag and a matrix containing gene name, sample name, true and imputed values of every pseudo missing combined from n.bag datasets
#' @export
#' @examples
#' \dontrun{
#' data(datapnnl)
#' data<-datapnnl.rm.ref[1:100,1:21]
#' impute<- DreamAI_Bagging(data=data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_bagging=NA,gamma_ADMIN=NA,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40,method = c("KNN", "MissForest", "ADMIN", "Birnn", "SpectroFM", "RegImpute"),out=c("Ensemble"),SamplesPerBatch=3,n.bag=2,save.out=TRUE,path="C:\\Users\\chowds14\\Desktop\\test_package\\",ProcessNum=1)
#' impute$Ensemble
#' }
DreamAI_Bagging<-function(data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),
                          gamma_bagging = NA, gamma_ADMIN=NA,gamma=50,
                          CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40,
                          method = c("KNN", "MissForest", "ADMIN", "Birnn", "SpectroFM", "RegImpute"),out=c("Ensemble"),
                          seed.bags = NULL,SamplesPerBatch,n.bag,save.out=TRUE,path=NULL,ProcessNum)
{
 
  data.pnnl = data[rowMeans(is.na(data))!=1,];
  
  # data.pnnl.b <<- avg.batch(data.pnnl,SamplesPerBatch=SamplesPerBatch)
  data.pnnl.b <- avg.batch(data.pnnl,SamplesPerBatch=SamplesPerBatch)
  if(is.na(gamma_bagging)){gm.pnnl <- gamma.est(data.pnnl.b)}else{gm.pnnl = c(NA,gamma_bagging)}
  
  N <- sum(is.na(data.pnnl.b[!(apply(is.na(data.pnnl.b), 1, mean)==1),]))
  
  
  gm1.new = optimize(f = diff.mr,interval = c(-10,10),gm.pnnl=gm.pnnl,data.pnnl.b=data.pnnl.b,N=N)$minimum
  
  p.pnnl = exp(gm1.new-gm.pnnl[2]*data.pnnl.b);
  p.pnnl[is.na(p.pnnl)] = 1;
  
  # results<-NULL
  data.v<-c(as.matrix(data))
  true.miss <- which(is.na(data.v))
  bag.knn.sum<-0
  bag.MissForest.sum<-0
  bag.ADMIN.sum<-0
  bag.BruinGo.sum<-0
  bag.DMIS.sum<-0
  bag.Jeremy.sum<-0
  bag.ensemble.sum<-0
  summary.all<-list()

  if(length(seed.bags)!=n.bag)
  seed.bags = 1:n.bag+1000000
    
  for (i in 1:n.bag)
  {
    # TimeStart<-proc.time()
    
    set.seed(seed.bags[i]);
    m.missing.b = matrix(rbinom(length(p.pnnl), 1, c(p.pnnl)),dim(p.pnnl)[1]);
    m.missing = t(apply(m.missing.b,1,rep,each=SamplesPerBatch));
    # sum(m.missing.b-is.na(data.pnnl.b))
    #print(sum(m.missing!=t(is.na(apply(data.pnnl.b,1,rep,each=SamplesPerBatch)))));
    
    data.pnnl.new = data.pnnl;
    data.pnnl.new<-as.matrix(data.pnnl.new)
    data.pnnl.new[m.missing!=t(is.na(apply(data.pnnl.b,1,rep,each=SamplesPerBatch)))] <- NA;
    
    
    ResultDreamImputation<-DreamAI(data=data.pnnl.new,k=k,maxiter_MF = maxiter_MF, ntree = ntree,maxnodes = maxnodes,maxiter_ADMIN=maxiter_ADMIN,tol=tol,gamma_ADMIN=gamma_ADMIN,gamma=gamma,CV=CV,fillmethod=fillmethod,maxiter_RegImpute=maxiter_RegImpute,conv_nrmse = conv_nrmse,iter_SpectroFM=iter_SpectroFM,method=method,out=out)
    
    ########## summary ###############
    d.o<-as.matrix(data)
    d.p<-data.pnnl.new
    d.i<-ResultDreamImputation$Ensemble
    
    d.s = cbind(c(d.p),c(d.i))
    
    ind.pseudo = as.logical(is.na(c(d.p))-is.na(c(d.o)));
    
    d.s = cbind.data.frame(gene = rep(rownames(d.o),ncol(d.o))[ind.pseudo],
                           sample = rep(colnames(d.o),each = nrow(d.o))[ind.pseudo],
                           true = (c(d.o))[ind.pseudo],
                           imputed = c(d.i)[ind.pseudo]);
    
    summary.all[[i]]<-d.s
    
    methods<-c("KNN","MissForest","ADMIN","Birnn","SpectroFM","RegImpute")
    
    if("KNN" %in% out){
      bag.knn.sum <- bag.knn.sum+(c(ResultDreamImputation$KNN))[true.miss]
    }
    
    if("MissForest" %in% out){
    bag.MissForest.sum <- bag.MissForest.sum+(c(ResultDreamImputation$MissForest))[true.miss]
    }
    
    if("ADMIN" %in% out){
    bag.ADMIN.sum <- bag.ADMIN.sum+(c(ResultDreamImputation$ADMIN))[true.miss]
    }
    
    if("Birnn" %in% out){
    bag.BruinGo.sum <- bag.BruinGo.sum+(c(ResultDreamImputation$Birnn))[true.miss]
    }
    
    if("SpectroFM" %in% out){
    bag.DMIS.sum <- bag.DMIS.sum+(c(ResultDreamImputation$SpectroFM))[true.miss]
    }
    
    if("RegImpute" %in% out){
    bag.Jeremy.sum <- bag.Jeremy.sum+(c(ResultDreamImputation$RegImpute))[true.miss]
    }
    
    bag.ensemble.sum <- bag.ensemble.sum+(c(ResultDreamImputation$Ensemble))[true.miss]
        
    cat("\n\n")
    print(paste("imputation on bagging data ",i," completed",sep=""))
    cat("\n\n")
  }
  
  summary.all1 = do.call(rbind,summary.all);
  
  imputed_matrix<-list()
  
  if("KNN" %in% out){
    bag.knn.v<- bag.knn.sum/n.bag
    data.v[true.miss]<-bag.knn.v
    bag.knn<-matrix(data.v,nrow(data))
    imputed_matrix<-c(imputed_matrix,list("KNN"=bag.knn))
  }

  if("MissForest" %in% out){
  bag.MissForest.v<- bag.MissForest.sum/n.bag
  data.v[true.miss]<-bag.MissForest.v
  bag.MissForest<-matrix(data.v,nrow(data))
  imputed_matrix<-c(imputed_matrix,list("MissForest"=bag.MissForest))
  }
  
  if("ADMIN" %in% out){
  bag.ADMIN.v<- bag.ADMIN.sum/n.bag
  data.v[true.miss]<-bag.ADMIN.v
  bag.ADMIN<-matrix(data.v,nrow(data))
  imputed_matrix<-c(imputed_matrix,list("ADMIN"=bag.ADMIN))
  }
  
  if("Birnn" %in% out){
  bag.BruinGo.v<- bag.BruinGo.sum/n.bag
  data.v[true.miss]<-bag.BruinGo.v
  bag.BruinGo<-matrix(data.v,nrow(data))
  imputed_matrix<-c(imputed_matrix,list("Birnn"=bag.BruinGo))
  }
  
  if("SpectroFM" %in% out){
  bag.DMIS.v<- bag.DMIS.sum/n.bag
  data.v[true.miss]<-bag.DMIS.v
  bag.DMIS<-matrix(data.v,nrow(data))
  imputed_matrix<-c(imputed_matrix,list("SpectroFM"=bag.DMIS))
  }
  
  if("RegImpute" %in% out){
  bag.Jeremy.v<- bag.Jeremy.sum/n.bag
  data.v[true.miss]<-bag.Jeremy.v
  bag.Jeremy<-matrix(data.v,nrow(data))
  imputed_matrix<-c(imputed_matrix,list("RegImpute"=bag.Jeremy))
  }
  
  bag.ensemble.v<- bag.ensemble.sum/n.bag
  data.v[true.miss]<-bag.ensemble.v
  bag.ensemble<-matrix(data.v,nrow(data))
  imputed_matrix<-c(imputed_matrix,list("Ensemble"=bag.ensemble))
  

  bag.output<-list(impute=imputed_matrix,n.bag=n.bag,summary=summary.all1,out.method=out)
  
  # num<-which(methods %in% method)
  # bag.imputed_matrix[[1]]<-bag.imputed_matrix[[1]][c(num,7)]
  # bag.output<-bag.imputed_matrix
  
  if(save.out){
    save(bag.output,file=paste(path,"bag_imputed_",sprintf("%03d",ProcessNum),".RData",sep=""))
  }else{
    return(bag.output)
  }
}

