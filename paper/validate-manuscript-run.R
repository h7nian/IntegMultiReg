# Inspect the scientific outputs, not merely the replication completion log.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(1L,2L))
library(IntegMultiReg)
library(survival)
validate_run <- function(root) {
 root <- normalizePath(root,mustWork=TRUE)
 result <- readRDS(file.path(root,'replication_results.rds'))
 stopifnot(!result$quick,result$sampler_method=='paper',result$cv_method=='refit')
 outcomes <- list(binary=simIMR$outcome.binary,continuous=simIMR$outcome.continuous,
  right.censored=simIMR$outcome.survival)
 rows <- list()
 for(type in names(outcomes))for(method in c('imr','bms')) {
  x<-readRDS(file.path(root,paste0('validation-',type,'-',method,'.rds')))
  stopifnot(x$fit$control$mcmc$draws==4000L,x$fit$control$mcmc$burnin==1000L,
   x$fit$control$sampler_method=='paper',x$fit$control$method==method,
   x$cv$validation=='refit',x$cv$control$k==5L,x$cv$control$rounds==10L,
   x$cv$control$sampler_method=='paper',x$cv$control$score_method=='standard')
  predictions<-x$cv$predictions;folds<-x$cv$control$folds
  stopifnot(nrow(predictions)==3000L,all(is.finite(predictions$prediction)))
  key<-function(x)paste(x$id,x$round,sep=':')
  stopifnot(!anyDuplicated(key(predictions)),!anyDuplicated(key(folds)),
   setequal(key(predictions),key(folds)),
   identical(predictions$fold,folds$fold[match(key(predictions),key(folds))]))
  y<-outcomes[[type]]
  calculate<-function(p) {
   obs<-y[match(p$id,y$id),];v<-p$prediction
   if(type=='continuous')return(mean((obs$y-v)^2))
   if(type=='binary') {
    n1<-sum(obs$y==1);n0<-sum(obs$y==0)
    return((sum(rank(v)[obs$y==1])-n1*(n1+1)/2)/(n1*n0))
   }
   survival::concordance(survival::Surv(obs$time,obs$status)~v)$concordance
  }
  computed<-x$cv$pooled
  for(round in 1:10) {
   p<-predictions[predictions$round==round,]
   stopifnot(nrow(p)==300L,setequal(p$id,y$id),setequal(p$fold,1:5))
   for(group in colnames(computed))computed[round,group]<-
    calculate(if(group=='all')p else p[p$subgroup==group,])
  }
  stopifnot(isTRUE(all.equal(computed,x$cv$pooled,tolerance=1e-12)))
  table<-result$prediction_table
  index<-which(table$outcome==type & table$method==toupper(method))
  stopifnot(length(index)==1L,isTRUE(all.equal(unname(colMeans(computed)),
   unname(as.numeric(table[index,colnames(computed)])),tolerance=1e-12)))
  rows[[length(rows)+1L]]<-data.frame(outcome=type,method=method,
   rounds=10L,folds=5L,predictions=nrow(predictions),
   max_score_difference=max(abs(computed-x$cv$pooled)))
 }
 kirc<-readRDS(file.path(root,'kirc_reduced_fit.rds'))
 stopifnot(kirc$control$mcmc$draws==4000L,kirc$control$mcmc$burnin==1000L,
  kirc$control$sampler_method=='paper',result$simulated_fit$control$mcmc$draws==8000L,
  result$simulated_fit$control$mcmc$burnin==2000L)
 comparison <- readRDS(file.path(root,'covariate-comparison/comparison.rds'))
 stopifnot(!comparison$settings$quick,comparison$settings$draws==2000,comparison$settings$burnin==500,
  comparison$settings$k==3,comparison$settings$rounds==2)
 outer <- comparison$outer_predictions
 stopifnot(nrow(outer)==300L,!anyDuplicated(outer$id),setequal(outer$id,simIMR$outcome.continuous$id),
  all(is.finite(outer$prediction)),
  isTRUE(all.equal(mean((outer$observed-outer$prediction)^2),
   comparison$nested_summary$pooled_outer_mse,tolerance=1e-12)))
 for(fold in 1:3) {
  test <- outer[outer$outer_fold==fold,]
  inner <- comparison$nested_inner_folds[comparison$nested_inner_folds$outer_fold==fold,]
  stopifnot(nrow(inner)==200L,!anyDuplicated(inner$id),
   setequal(inner$id,setdiff(outer$id,test$id)),setequal(inner$fold,1:3),
   !any(inner$id %in% test$id))
  selected <- comparison$selected_by_fold[comparison$selected_by_fold$outer_fold==fold,]
  chosen <- c('age_only','age_sex_stage')[which.min(c(selected$inner_mse_age_only,selected$inner_mse_age_sex_stage))]
  stopifnot(selected$selected==chosen,selected$n_train==nrow(inner),
   selected$n_test==nrow(test),isTRUE(all.equal(mean((test$observed-test$prediction)^2),
    selected$outer_mse,tolerance=1e-12)))
 }
 for(candidate in c('age_only','age_sex_stage')) {
  selected <- comparison$paired_summary[comparison$paired_summary$candidate==candidate,]
  v <- comparison$paired_rounds[[candidate]]
  stopifnot(isTRUE(all.equal(mean(v),selected$mean_mse,tolerance=1e-12)),
   isTRUE(all.equal(sd(v),selected$sd_across_rounds,tolerance=1e-12)))
 }
 report<-do.call(rbind,rows)
 write.csv(report,file.path(root,'manuscript-acceptance.csv'),row.names=FALSE)
 cat('PASS: six full refit configurations, independent scores, nested holdout separation and summary aggregation:',root,'\n')
 invisible(report)
}
for(root in args)validate_run(root)
if(length(args)==2L) {
 a<-normalizePath(args[1]);b<-normalizePath(args[2])
 csv<-setdiff(list.files(a,pattern='[.]csv$',recursive=TRUE),'manuscript-acceptance.csv')
 rds<-c(paste0('validation-',rep(c('binary','continuous','right.censored'),each=2),'-',c('imr','bms'),'.rds'),
  'kirc_reduced_fit.rds','replication_results.rds','covariate-comparison/comparison.rds')
 report<-lapply(c(csv,rds),function(name){
  read<-if(endsWith(name,'.csv'))function(p)read.csv(p,check.names=FALSE) else readRDS
  x<-read(file.path(a,name));y<-read(file.path(b,name))
  same<-isTRUE(all.equal(x,y,tolerance=1e-12))
  data.frame(file=name,identical=identical(x,y),equal_1e12=same)
 })
 report<-do.call(rbind,report)
 write.csv(report,file.path(dirname(b),'AB-numeric-comparison.csv'),row.names=FALSE)
 stopifnot(all(report$equal_1e12))
 cat('PASS:',nrow(report),'CSV/RDS artifacts agree across A/B\n')
}
