# Audit the executed simulation manifest; a one-replicate run is not a full study.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% 1:2)
root <- normalizePath(args[1])
full_study <- length(args)==2L && identical(args[2],'--require-complete-study')
if(length(args)==2L && !full_study) stop('Unknown option')
s <- readRDS(file.path(root,'settings.rds'))
stopifnot(s$experiment %in% c('simulation','correlated'))
if(full_study) stopifnot(!s$quick,s$draws==350000L,s$burnin==50000L,
 s$replicates==if(s$experiment=='simulation')50L else 30L)
stopifnot(identical(tools::md5sum(names(s$source_hashes)),s$source_hashes),
 identical(tools::md5sum(names(s$data_hash)),s$data_hash))
source(names(s$source_hashes)[basename(names(s$source_hashes))=='original-experiment-helpers.R'])
read_result <- function(p) {
 x<-readRDS(p)
 if(inherits(x,'imr_experiment_result_v1')) {
  x<-x$result;x$fit<-unpack_experiment_fit(x$fit)
 }
 x
}
# Independent pairwise AUC: probability that a positive exceeds a negative.
auc_pairwise <- function(score,truth) {
 a<-score[truth==1];b<-score[truth==0]
 if(!length(a)||!length(b))return(NA_real_)
 cmp<-outer(a,b,'-');mean((cmp>0)+.5*(cmp==0))
}
reports <- list()
configs <- if(s$experiment=='simulation')3L else 6L
selection <- list(configurations = seq_len(configs), replicates = seq_len(s$replicates))
if (!full_study && file.exists(file.path(root, "task.rds"))) {
 selection <- readRDS(file.path(root, "task.rds"))
 stopifnot(identical(selection$kind, s$experiment),
   all(selection$configurations %in% seq_len(configs)),
   all(selection$replicates %in% seq_len(s$replicates)))
}
for(config in selection$configurations) for(rep in selection$replicates) {
 job<-file.path(root,sprintf('configuration-%02d-replicate-%03d',config,rep))
 data<-readRDS(file.path(job,'data.rds'))$data
 for(method in c('imr-molecular','bms-molecular','l1-cph')) {
  path<-file.path(job,paste0(method,'.rds'));status<-readRDS(paste0(path,'.status.rds'))
  stopifnot(identical(status$status,'completed'))
  result<-read_result(path)
  if(method!='l1-cph') {
   fit<-result$fit;cv<-result$cv
   stopifnot(identical(cv$validation,s$reference_arguments$cv$cv_method))
   stopifnot(fit$control$mcmc$draws==s$draws,fit$control$mcmc$burnin==s$burnin,
     length(fit$posterior$selection_draws)==s$draws,
     fit$control$seed==s$seed+rep-1L,cv$control$k==s$k,cv$control$rounds==s$rounds,
     identical(cv$control$folds,experiment_result_folds(path)))
   for(name in names(s$reference_arguments$fit)) {
    actual<-if(name=='sampler_method')fit$control[[name]] else fit$control$numerical[[name]]
    stopifnot(identical(actual,s$reference_arguments$fit[[name]]))
   }
   for(name in setdiff(names(s$reference_arguments$cv),'cv_method'))
    if(!(name=='fold_rng' && method=='bms-molecular' &&
         !identical(s$reference_arguments$cv$fold_rng,'continue')))
      stopifnot(identical(cv$control[[name]],s$reference_arguments$cv[[name]]))
   selection<-result$selection
   for(p in seq_along(data$truth)) {
    groups<-fit$model$subgroup_names[fit$model$platform_subgroups[[p]]]
    truth<-data$truth[[p]][groups,,drop=FALSE]
    prob<-fit$posterior$inclusion_probabilities[[p]]
    expected<-vapply(seq_along(groups),function(i)auc_pairwise(prob[i,],truth[i,]),numeric(1))
    got<-selection$auc[match(paste(fit$model$platform_names[p],groups),paste(selection$platform,selection$subgroup))]
    stopifnot(isTRUE(all.equal(unname(got),unname(expected),tolerance=1e-12)))
   }
   pred<-cv$predictions;folds<-cv$control$folds;score_method<-cv$control$score_method
   if(method=='imr-molecular')first_folds<-folds
   result$fit<-NULL;rm(fit,cv);gc()
  } else {pred<-result$predictions;folds<-first_folds;score_method<-'standard'}
  stopifnot(all(is.finite(pred$prediction)),nrow(pred)==nrow(data$outcome)*s$rounds)
  key<-function(x)paste(x$id,x$round,sep=':')
  stopifnot(!anyDuplicated(key(pred)),setequal(key(pred),key(folds)),
    identical(pred$fold,folds$fold[match(key(pred),key(folds))]))
  for(round in seq_len(s$rounds))stopifnot(setequal(pred$id[pred$round==round],data$outcome$id))
  calculated<-fold_scores(pred,data$outcome,score_method,weighted_overall=method=='l1-cph')
  stopifnot(isTRUE(all.equal(result$scores,calculated,tolerance=1e-12)),
    isTRUE(all.equal(result$summary,summarize_scores(calculated),tolerance=1e-12)))
  reports[[length(reports)+1L]]<-data.frame(config,rep,method,draws=s$draws,
    rounds=s$rounds,k=s$k,predictions=nrow(pred),warnings=length(status$warnings))
  cat('PASS:',config,rep,method,'\n');rm(result);gc()
 }
 path<-file.path(job,'uni-cph-selection.csv')
 stopifnot(readRDS(paste0(path,'.status.rds'))$status=='completed')
 uni<-read.csv(path,colClasses=c(subgroup='character'))
 stopifnot(all(is.finite(uni$auc)),all(uni$auc>=0 & uni$auc<=1))
}
report<-do.call(rbind,reports)
write.csv(report,file.path(root,'simulation-acceptance.csv'),row.names=FALSE)
saveRDS(list(settings=s,report=report,full_study=full_study,verified_at=Sys.time(),
 scope='Executed manifest, fit settings, pairwise IMR/BMS selection AUC, folds and predictive scores. Univariate/L1 scientific stability needs separate warning review.'),
 file.path(root,'simulation-acceptance.rds'))
cat('PASS executed manifest. Complete-study requirement:',full_study,'\n')
