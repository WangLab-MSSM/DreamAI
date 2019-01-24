#Iterative ridge regression imputation for biological data sets
# set.seed(1)

#default args
# filepath = "./SC1_CPTAC2_OV_retrospective_proteome_PNNL_for_R.csv"
# delim = "\t"
# outpath = "./imputed_out.tab"
# fillmethod = "row_mean"
# n_iter = 100


# args = commandArgs(trailingOnly=TRUE)
#
# if (length(args)==1) {
#   stop("inpath is a required argument..", call.=FALSE)
# }
#
#
# for(arg in args){
#
#   if(grepl('sep',arg)){
#     delim = strsplit(arg,"=")[[1]][2]
#   }
#   if(grepl('filepath',arg)){
#     filepath = strsplit(arg,"=")[[1]][2]
#   }
#   if(grepl('outpath',arg)){
#     outpath = strsplit(arg,"=")[[1]][2]
#   }
#   if(grepl('fillmethod',arg)){
#     fillmethod = strsplit(arg,"=")[[1]][2]
#   }
#   if(grepl('n_iter',arg)){
#     n_iter = as.numeric(strsplit(arg,"=")[[1]][2])
#   }
#   if(grepl('conv_nrmse',arg)){
#     conv_nrmse = as.double(strsplit(arg,"=")[[1]][2])
#   }
# }





################################################## BACK FILLING METHOD #############

#' Fill-in missing values using simple imputation.
#'
#' @param data A matrix of values.
#' @param fillmethod A string identifying the method to be used.
#' @return The imputed matrix of \code{data}.
backfill <- function(data,fillmethod=fillmethod){

  if ((fillmethod!="zeros") & (fillmethod!="row_mean")) {
    stop("Accepted Fill Methods Include \"zeros\",\"row_mean\"", call.=FALSE)
  }

  print(paste("Back Filling With Method: ",fillmethod))

  if(fillmethod == "zeros"){
    missing_indices = which(is.na(data), arr.ind=TRUE)
    data[missing_indices] = 0
  }
  if(fillmethod == "row_mean"){
    ind <- which(is.na(data), arr.ind=TRUE)
    data[ind] <- rowMeans(data,  na.rm = TRUE)[ind[,1]]
  }

  return(data)
}



#' Return current column indices where values are not missing.
#' @param col A number indicating current columns index.
#' @param data A matrix, The original data matrix with NaNs.
#' @param filled A matrix, The data matrix after running backfill function.
#' @return A list of indices corresponding to current working column values.
returnPresentIndices<-function(col,data,filled){
  #Indices @ current column, only where values aren't missing
  present_indices = which(!is.na(data[,col]), arr.ind=TRUE)
  return(present_indices)
}

#' Return training set for current column.
#' @param col A number indicating current columns index.
#' @param filled A matrix, The data matrix after running backfill function.
#' @param present_indices A list containing indices in current column where values are not missing.
#' @return A matrix.

returnTrainingSet<-function(col,filled,present_indices){
  #Only rows where value is present at current column.  All columns except for current.
  train = filled[,-col][present_indices,]
  return(train)
}

#' Return target set for current column.
#' @param col A number indicating current columns index.
#' @param data A matrix, The original data matrix with NaNs.
#' @param present_indices A list containing indices in current column where values are not missing.
#' @return A list of length n.
returnTargets<-function(col,data,present_indices){
  #Current column, only indices where values were present in original data
  target = data[present_indices,col]
  return(target)
}

#' Return current column indices where values are missing.
#' @param col A number indicating current columns index.
#' @param data A matrix, The original data matrix with NaNs.
#' @return A list of indices corresponding to current working column cells that contain NaN.
returnTestIndices<-function(col,data){
  #Indices @ current column, only where values are missing
  test_indices = which(is.na(data[,col]), arr.ind=TRUE)
  return(test_indices)
}

#' Return back filled data for test indices.
#' @param col A number indicating current columns index.
#' @param test_indices A list of indices corresponding to cells in current column that contain NaN.
#' @return A matrix
returnTestSet<-function(col,filled,test_indices){
  test = filled[,-col][test_indices,]
  return(test)
}


#' Imputation of Missing Protein Abundance using Glmnet Ridge Regression.
#'@description The function impute.RegImpute imputes a dataset with missing values or NA's using Glmnet Ridge Regression
#' @param data dataset in the form of a matrix or data frame with either NAs or 0s as missings.
#' @param fillmethod a string identifying the method to be used that could be "row_mean" or "zeros", with "row_mean" being the default. It throws an warning if "row_median" is used.
#' @param maxiter_RegImpute a integer identifying maximum number of iterations to reach convergence
#' @return the imputed version of the dataset
#' @export
#' @examples
#' \dontrun{
#' data(datapnnl)
#' data<-datapnnl.rm.ref[1:100,1:21]
#' impute.RegImpute(data=as.matrix(data), fillmethod = "row_mean", maxiter_RegImpute = 10,conv_nrmse = 1e-06)
#' }
impute.RegImpute <- function(data,fillmethod,maxiter_RegImpute,conv_nrmse){
  filled = backfill(data,fillmethod)
  print(paste("Starting Imputation With ",toString(maxiter_RegImpute)," Max. Iterations"))
  missing_indices = which(is.na(data), arr.ind=TRUE)


  #Begin iterative imputation.  Missing value slots in data are updated on each iteration
  for(i in 1:maxiter_RegImpute){
    print(paste("Working on Iteration: ",toString(i),"/",toString(maxiter_RegImpute)))

        
    #Iterate over columns
    for(col in 1:dim(data)[2]){

      if(sum(is.na(data[,col]))==0){next} #continue if there aren't any missing values in the current column

      #otherwise, begin imputation
      present_indices = returnPresentIndices(col,data,filled)
      train = returnTrainingSet(col,filled,present_indices)
      target = returnTargets(col,data,present_indices)
      test_indices = returnTestIndices(col,data)
      test = returnTestSet(col,filled,test_indices)

        #used http://ricardoscr.github.io/how-to-use-ridge-and-lasso-in-r.html as tutorial on ridge regression in R
        cv_fit <- cv.glmnet(train, target, alpha=0, standardize=TRUE)
        opt_lambda = cv_fit$lambda.min
        fit = cv_fit$glmnet.fit
        
        test<-as.matrix(test)
        if(dim(test)[2]==1)
        {
          predicted = predict(fit, s = opt_lambda, newx=t(test))
        }else{
          predicted = predict(fit, s = opt_lambda, newx=test)
        }        
        filled[test_indices,col] = predicted
    }

	#test for convergence at each iteration
   	if(i==1){
	impvals = filled[missing_indices]
	}
	else{
	NRMSE = sqrt(mean((impvals - filled[missing_indices])^2))
	print(paste("NRMSE = ",NRMSE))
	impvals = filled[missing_indices]
	if (NRMSE<conv_nrmse){
	return(filled)
	}
	}
  }
  return(filled)
}

