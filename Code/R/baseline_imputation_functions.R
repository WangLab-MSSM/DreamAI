
my.normlize = function(data)
{
  median.temp = apply(data,1,median,na.rm=T);
  sd.temp = apply(data,1,sd,na.rm=T);
  # apply(t(t(data)-median.temp),2,median,na.rm=T);
  # apply(t((t(data)-median.temp)/sd.temp),2,sd,na.rm=T);
  data.new = ((data-median.temp)/sd.temp);
  return(list(data = data.new, median = median.temp, sd=sd.temp));
}

my.normlize.rev = function(data, median.old, sd.old)
{
  data.new = (data*sd.old+median.old);
}

gamma.est = function(data)
{
  ### data dim Lxn
  Y.mean = apply(data,1,mean,na.rm=T);
  q = apply(is.na(data),1,mean);
  par = lsfit(Y.mean[q>0],log(1/q)[q>0])[[1]];
  names(par) = c('gamma0','gamma');
  return(par);
}


### 2 knn linear regression
knn.est.it2 = function(data, k, m.ind=T)
{
  L = dim(data)[1];
  ccc = cor(t(data));
  order.k = apply(ccc,1,function(x){order(x,decreasing = T)[1:k+1]});

  data.est = data;
  data.est[m.ind,] = t(sapply((1:L)[m.ind],function(l){predict(lm(data[l,]~t(data[order.k[,l],])))}));
  dimnames(data.est) = dimnames(data);
  data.est[data.est<0] = 0;
  return(data.est);
}

#########################################

#' @title Imputation of Missing Protein Abundances using K Nearest Neighbour
#' @description The function impute.KNN imputes a dataset with missing values or NA's using k nearest neighbour
#' @param data dataset in the form of a matrix or data frame with either NAs or 0s as missings
#' @param k number of neighbors to be used in the imputation (default=10)
#'
#' @return the imputed version of the dataset
#' @export
#'
#' @examples
#' \dontrun{
#' data(datapnnl)
#' data<-datapnnl.rm.ref[1:100,1:21]
#' impute.KNN(data=as.matrix(data), k=10)
#' }
impute.KNN = function(data,k)
{
  if (!requireNamespace("impute", quietly = TRUE)) {
    stop("\n package ", "impute", " is not yet installed \n", 
         "To install: \n", "BiocManager::install("impute") \n",  
         call. = FALSE)  
  }
  
  library(impute)
  
  norm.temp = my.normlize((data));

  data.new = (norm.temp[[1]]);
  median.new = norm.temp[[2]];
  sd.new = norm.temp[[3]];

  impu.temp = impute.knn(data.new,k,rowmax = 0.9,colmax = 0.9,maxp=dim(data)[1]);

  X1b.new = my.normlize.rev((impu.temp[[1]]), median.new, sd.new);
  # X1b.new[X1b.new<0] = 0;
  return((X1b.new));

}



#' @title Imputation of Missing Protein Abundances using MissForest
#' @description The function impute.MF imputes a dataset with missing values or NA's using MissForest
#' @param data dataset in the form of a matrix or data frame with either NAs or 0s as missing values
#' @param maxiter_MF maximum number of iteration to be performed if the stopping criteria is not met beforehand
#' @param ntree number of trees to grow in each forest
#' @param maxnodes maximum number of terminal nodes for trees in the forest, has to equal at least the number of columns in the given data
#' @return the imputed version of the dataset
#' @export
#'
#' @examples
#' #' #' \dontrun{
#' data(datapnnl)
#' data<-datapnnl.rm.ref[1:100,1:21]
#' impute.MF(data=as.matrix(data), maxiter_MF = 10, ntree = 100,maxnodes = NULL)
#' }
impute.MF = function(data = as.matrix(data),maxiter_MF, ntree, maxnodes)
{
  if(is.null(maxnodes))
  {
    maxnodes<-ncol(as.matrix(data))
  }
  mF.temp = missForest(t(data), maxiter = maxiter_MF, ntree = ntree, variablewise = FALSE,
                       decreasing = FALSE, verbose = FALSE,
                       mtry = ceiling(sqrt(dim(data)[1])), replace = TRUE,
                       classwt = NULL, cutoff = NULL, strata = NULL,
                       sampsize = NULL, nodesize = NULL, maxnodes = maxnodes,
                       xtrue = NA, parallelize = c('no', 'variables', 'forests'));
  X.temp = t(mF.temp[[1]]);
  return(X.temp);
}




#' @title Imputation of Missing Protein Abundances using ADMIN (Abundance  Dependent  Missing Imputation Mechanism)
#'@description The function impute.ADMIN imputes a dataset with missing values or NA's using ADMIN
#' @param data dataset in the form of a matrix or data frame with either NAs or 0s as missings
#' @param data.ini initial dataset, set to be NA as default
#' @param gamma parameter of the non-ignorable missing mechanism
#' @param k number of neighbors to be used in the imputation (default=10)
#' @param maxiter_ADMIN maximum number of iteration to be performed if the stopping criteria is not met beforehand
#' @param tol convergence threshold
#'
#' @return the imputed version of the dataset
#' @export
#'
#' @examples
#' \dontrun{
#' data(datapnnl)
#' data<-datapnnl.rm.ref[1:100,1:21]
#' impute.ADMIN(data = as.matrix(data),k = 10, maxiter_ADMIN = 30, tol = 10^(-2))
#' }
impute.ADMIN = function(data,data.ini=NA,gamma, k, maxiter_ADMIN,tol)
{
   if (!requireNamespace("impute", quietly = TRUE)) {
    stop("\n package ", "impute", " is not yet installed \n", 
         "To install: \n", "BiocManager::install("impute") \n",  
         call. = FALSE)  
  }
  
  library(impute)
  
  L = dim(data)[1];
  iter = 0;
  diff = 999;
  llh = Inf;
  if(is.na(gamma))
  {gamma = gamma.est(data)[2]}

  M.t = data.ini;
  if(sum(is.na(data.ini))>0){
    set.seed(1234);
    fit.0 = impute.knn(data,k=5,rowmax = 0.9,colmax = 0.9,maxp=dim(data)[1])[[1]];
    M.t = fit.0;
  }

  d.t = 0;

  while ((iter<maxiter_ADMIN) & (diff>tol))
  {
    gc();
    iter = iter+1;

    X.temp = (M.t-is.na(data)*gamma*d.t);

    M.t = knn.est.it2(X.temp, k);

    d.t = apply(M.t-X.temp,1,var);

    # X.temp = M.t;
    M.diff = (M.t-data);
    M.diff[is.na(M.diff)]=0;

    llh = c(llh,-norm(M.diff,type = 'F')^2-sum((M.t*gamma*d.t)[is.na(data)]));

    # diff = abs(diff(llh)[iter]/llh[iter+1]);
    diff = abs(diff(llh)[iter]);

    M.t[!is.na(data)] = data[!is.na(data)];

    print(c(iter,llh[iter+1]));
  }
  dimnames(M.t) = dimnames(data);
  # return(list(M.t,llh))
  return(M.t);
}



