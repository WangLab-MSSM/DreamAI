#' Wrapper function for summarizing the outputs from DreamAI_bagging
#'
#' @param method a vector of imputation methods: ("KNN", "MissForest", "ADMIN", "Birnn", "SpectroFM, "RegImpute", "Ensemble"). This vector should be a subset or equal to the vector out in DreamAI_bagging.
#' @param nNodes number of parallel processes 
#' @param path location where the bagging output is saved
#'
#' @return  list of final imputed data and confidence score for every gene using pseudo missing
#' @export 
#'
#' @examples 
#' \dontrun{
#' data(datapnnl)
#' data<-datapnnl.rm.ref[1:100,1:21]
#' impute<- DreamAI_Bagging(data=data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40,method=c("KNN","MissForest","ADMIN","Birnn","SpectroFM","RegImpute"),out=c("Ensemble"),SamplesPerBatch=3,n.bag=2,save.out=TRUE,path="C:\\Users\\chowds14\\Desktop\\test_package\\",ProcessNum=1)
#' final.out<-bag.summary(method=c("Ensemble"),nNodes=2,path="C:\\Users\\chowds14\\Desktop\\test_package\\")
#' final.out$score
#' final.out$imputed_data
#' }
bag.summary<-function(method=c("KNN", "MissForest", "ADMIN", "Birnn", "SpectroFM", "RegImpute","Ensemble"),nNodes=2,path=NULL)
{
  load(paste(path,"bag_imputed_",sprintf("%03d",1),".RData",sep=""))
  out<-bag.output$out.method
  
  if(identical(method,intersect(out,method))==FALSE){
       return(print("method does not match out"))
  }
  
  summary.out<-list()
  
  KNN.out<-list()
  MF.out<-list()
  ADMIN.out<-list()
  Reg_Impute.out<-list()
  Birnn.out<-list()
  SpectroFM.out<-list()
  Ensemble.out<-list()
  
  n.bag.out<-list()

  for(i in 1:nNodes)
  {
    load(paste(path,"bag_imputed_",sprintf("%03d",i),".RData",sep=""))
    
    summary.out[[i]]<-bag.output$summary
    
    if("KNN" %in% method){
      KNN.out[[i]]<-bag.output$impute$KNN
    }else{
      sink("NULL")
      print("No output for KNN")
      sink()
    }
    
    if("MissForest" %in% method){
      MF.out[[i]]<-bag.output$impute$MissForest
    }else{
      sink("NULL")
      print("No output for MissForest")
      sink()
      }
    
    
    if("ADMIN" %in% method){
      ADMIN.out[[i]]<-bag.output$impute$ADMIN
    }else{
      sink("NULL")
      print("No output for ADMIN")
      sink()
      }

    if("RegImpute" %in% method){
      Reg_Impute.out[[i]]<-bag.output$impute$RegImpute
    }else{
      sink("NULL")
      print("No output for RegImpute")
      sink()
      }
    
    
    if("Birnn" %in% method){
      Birnn.out[[i]]<-bag.output$impute$Birnn
    }else{
      sink("NULL")
      print("No output for Birnn")
      sink()
      }
    
    
    if("SpectroFM" %in% method){
      SpectroFM.out[[i]]<-bag.output$impute$SpectroFM
    }else{
      sink("NULL")
      print("No output for SpectroFM")
      sink()
      }
    
    
    Ensemble.out[[i]]<-bag.output$impute$Ensemble
    
    n.bag.out[[i]]<-bag.output$n.bag
  }
  
  summary.out.all<-do.call(rbind,summary.out)
  summary.out.all.agg = aggregate(cbind(true,imputed)~gene ,data = summary.out.all,function(x){x});
  
  get.score = function(x)
  {
    dt = x$true;
    di = x$imputed;
    cor.temp = cor(dt,di,use = 'pairwise.complete.obs',method = 'spearman');
    nrmsd.temp = sqrt(mean((dt-di)^2))/diff(range(dt));
    return(c(cor = cor.temp,nrmsd = nrmsd.temp));
  }
  
  summary.score = t(apply(summary.out.all.agg,1,get.score));
  summary.score = as.data.frame(summary.score);
  rownames(summary.score) = summary.out.all.agg$gene
  
  if("KNN" %in% method){
  d.impute.knn<-matrix(0,nrow(KNN.out[[1]]),ncol(KNN.out[[1]]))
  }
  
  if("MissForest" %in% method){
  d.impute.MF<-matrix(0,nrow(MF.out[[1]]),ncol(MF.out[[1]]))
  }
  
  if("ADMIN" %in% method){
  d.impute.ADMIN<-matrix(0,nrow(ADMIN.out[[1]]),ncol(ADMIN.out[[1]]))
  }
  
  if("Birnn" %in% method){
  d.impute.Birnn<-matrix(0,nrow(Birnn.out[[1]]),ncol(Birnn.out[[1]]))
  }
  
  if("RegImpute" %in% method){
  d.impute.Reg_Impute<-matrix(0,nrow(Reg_Impute.out[[1]]),ncol(RegImpute.out[[1]]))
  }
  
  if("SpectroFM" %in% method){
  d.impute.SpectroFM<-matrix(0,nrow(SpectroFM.out[[1]]),ncol(SpectroFM.out[[1]]))
  }
  
  d.impute.Ensemble<-matrix(0,nrow(Ensemble.out[[1]]),ncol(Ensemble.out[[1]]))
  n.bag.tot<-0
  
  for(i in 1:nNodes)
  {
    if("KNN" %in% method){
    d.impute.knn<-d.impute.knn+n.bag.out[[i]]*KNN.out[[i]]
    }
    
    if("MissForest" %in% method){
    d.impute.MF<-d.impute.MF+n.bag.out[[i]]*MF.out[[i]]
    }
    
    if("ADMIN" %in% method){
    d.impute.ADMIN<-d.impute.ADMIN+n.bag.out[[i]]*ADMIN.out[[i]]
    }
    
    if("Birnn" %in% method){
    d.impute.Birnn<-d.impute.Birnn+n.bag.out[[i]]*Birnn.out[[i]]
    }
    
    if("RegImpute" %in% method){
    d.impute.Reg_Impute<-d.impute.Reg_Impute+n.bag.out[[i]]*Reg_Impute.out[[i]]
    }
    
    if("SpectroFM" %in% method){
    d.impute.SpectroFM<-d.impute.SpectroFM+n.bag.out[[i]]*SpectroFM.out[[i]]
    }
    
    d.impute.Ensemble<-d.impute.Ensemble+n.bag.out[[i]]*Ensemble.out[[i]]
    n.bag.tot<-n.bag.tot + n.bag.out[[i]]
  }
  
  imputed_matrix<-list()
  
  if("KNN" %in% method){
  d.impute.knn.final<-d.impute.knn/n.bag.tot
  imputed_matrix<-c(imputed_matrix,list("KNN"=as.matrix(d.impute.knn.final)))
  }
  
  if("MissForest" %in% method){
  d.impute.MF.final<-d.impute.MF/n.bag.tot
  imputed_matrix<-c(imputed_matrix,list("MissForest"=as.matrix(d.impute.MF.final)))
  }
  
  if("ADMIN" %in% method){
  d.impute.ADMIN.final<-d.impute.ADMIN/n.bag.tot
  imputed_matrix<-c(imputed_matrix,list("ADMIN"=as.matrix(d.impute.ADMIN.final)))
  }
  
  if("Birnn" %in% method){
  d.impute.Birnn.final<-d.impute.Birnn/n.bag.tot
  imputed_matrix<-c(imputed_matrix,list("Birnn"=as.matrix(d.impute.Birnn.final)))
  }
  
  if("RegImpute" %in% method){
  d.impute.RegImpute.final<-d.impute.Reg_Impute/n.bag.tot
  imputed_matrix<-c(imputed_matrix,list("RegImpute"=as.matrix(d.impute.RegImpute.final)))
}
  
  if("SpectroFM" %in% method){
  d.impute.SpectroFM.final<-d.impute.SpectroFM/n.bag.tot
  imputed_matrix<-c(imputed_matrix,list("SpectroFM"=as.matrix(d.impute.SpectroFM.final)))
  }
  
  d.impute.Ensemble.final<-d.impute.Ensemble/n.bag.tot
  imputed_matrix<-c(imputed_matrix,list("Ensemble"=as.matrix(d.impute.Ensemble.final)))
  

  # methods<-c("KNN","MissForest","ADMIN","Birnn","SpectroFM","RegImpute")
  # 
  # num<-which(methods %in% method)
  #
  sink()
  
  output<-imputed_matrix
  
  out<-list(score=summary.score,imputed_data=output)
  return(out)
}
  
