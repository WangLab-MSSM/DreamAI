#' Wrapper function for summarizing the outputs from DreamAI_bagging
#'
#' @param method a vector of imputation methods: ("KNN", "MissForest", "ADMIN", "Brinn", "SpectroFM, "RegImpute", "Ensemble"). Default is "Ensemble" if nothing is specified
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
#' impute<- DreamAI_Bagging(data=data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40,method=c("KNN","MissForest","ADMIN","Brinn","SpectroFM","RegImpute","Ensemble"),SamplesPerBatch=3,n.bag=2,save.out=TRUE,path="C:\\Users\\chowds14\\Desktop\\test_package\\",ProcessNum=1)
#' final.out<-bag.summary(method=c("Ensemble"),nNodes=2,path="C:\\Users\\chowds14\\Desktop\\test_package\\")
#' final.out$score
#' final.out$imputed_data
#' }
bag.summary<-function(method=c("KNN","MissForest","ADMIN","Brinn","SpectroFM","RegImpute","Ensemble"),nNodes=2,path=NULL)
{
  summary.out<-list()
  KNN.out<-list()
  MF.out<-list()
  ADMIN.out<-list()
  Reg_Impute.out<-list()
  Brinn.out<-list()
  SpectroFM.out<-list()
  Ensemble.out<-list()
  n.bag.out<-list()
  for(i in 1:nNodes)
  {
    load(paste(path,"bag_imputed_",sprintf("%03d",i),".RData",sep=""))
    
    summary.out[[i]]<-bag.output$summary
    KNN.out[[i]]<-bag.output$impute$KNN
    MF.out[[i]]<-bag.output$impute$MissForest
    ADMIN.out[[i]]<-bag.output$impute$ADMIN
    Reg_Impute.out[[i]]<-bag.output$impute$RegImpute
    Brinn.out[[i]]<-bag.output$impute$Brinn
    SpectroFM.out[[i]]<-bag.output$impute$SpectroFM
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
  
  d.impute.knn<-matrix(0,nrow(KNN.out[[1]]),ncol(KNN.out[[1]]))
  d.impute.MF<-matrix(0,nrow(KNN.out[[1]]),ncol(MF.out[[1]]))
  d.impute.ADMIN<-matrix(0,nrow(KNN.out[[1]]),ncol(ADMIN.out[[1]]))
  d.impute.Brinn<-matrix(0,nrow(KNN.out[[1]]),ncol(Brinn.out[[1]]))
  d.impute.Reg_Impute<-matrix(0,nrow(Reg_Impute.out[[1]]),ncol(KNN.out[[1]]))
  d.impute.SpectroFM<-matrix(0,nrow(SpectroFM.out[[1]]),ncol(KNN.out[[1]]))
  d.impute.Ensemble<-matrix(0,nrow(Ensemble.out[[1]]),ncol(KNN.out[[1]]))
  n.bag.tot<-0
  
  for(i in 1:nNodes)
  {
    d.impute.knn<-d.impute.knn+n.bag.out[[i]]*KNN.out[[i]]
    d.impute.MF<-d.impute.MF+n.bag.out[[i]]*MF.out[[i]]
    d.impute.ADMIN<-d.impute.ADMIN+n.bag.out[[i]]*ADMIN.out[[i]]
    d.impute.Brinn<-d.impute.Brinn+n.bag.out[[i]]*Brinn.out[[i]]
    d.impute.Reg_Impute<-d.impute.Reg_Impute+n.bag.out[[i]]*Reg_Impute.out[[i]]
    d.impute.SpectroFM<-d.impute.SpectroFM+n.bag.out[[i]]*SpectroFM.out[[i]]
    d.impute.Ensemble<-d.impute.Ensemble+n.bag.out[[i]]*Ensemble.out[[i]]
    n.bag.tot<-n.bag.tot + n.bag.out[[i]]
  }
  
  d.impute.knn.final<-d.impute.knn/n.bag.tot
  d.impute.MF.final<-d.impute.MF/n.bag.tot
  d.impute.ADMIN.final<-d.impute.ADMIN/n.bag.tot
  d.impute.Brinn.final<-d.impute.Brinn/n.bag.tot
  d.impute.RegImpute.final<-d.impute.Reg_Impute/n.bag.tot
  d.impute.SpectroFM.final<-d.impute.SpectroFM/n.bag.tot
  d.impute.Ensemble.final<-d.impute.Ensemble/n.bag.tot
  
  imputed_matrix<-list("KNN"=as.matrix(d.impute.knn.final),"MissForest"=as.matrix(d.impute.MF.final),"ADMIN"=as.matrix(d.impute.ADMIN.final),"Brinn"=as.matrix(d.impute.Brinn.final),"SpectroFM"=as.matrix(d.impute.SpectroFM.final),"RegImpute"=as.matrix(d.impute.RegImpute.final),"Ensemble"=as.matrix(d.impute.Ensemble.final))
  
  methods<-c("KNN","MissForest","ADMIN","Brinn","SpectroFM","RegImpute")
  
  num<-which(methods %in% method)
  
  output<-imputed_matrix[c(num,7)]
  
  out<-list(score=summary.score,imputed_data=output)
  return(out)
}
  
