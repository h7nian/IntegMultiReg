native_chain <- function(groups, feature_platform, nu, draws=30000L, burnin=2000L,
                         seed=1L, model_variant='imr', initial='empty', keep_latent=TRUE,
                         interaction_prior=c(shape=2,rate=2)) {
 n_groups<-length(groups);n_platforms<-length(nu)
 spec<-lapply(groups,function(g) list(design=unname(g$X),response=as.double(g$y),
  status=if(is.null(g$status))integer()else as.integer(g$status),
  prior_scale=as.double(g$prior_scale),residual_prior=as.double(g$residual_prior),
  forced=as.integer(g$n_forced),outcome=as.integer(match(g$outcome_type,c('continuous','binary','right.censored'))-1L)))
 block<-lapply(seq_len(n_platforms),function(l) {
  features<-which(feature_platform==l)
  members<-which(vapply(groups,function(g)all(features%in%g$feature_index),TRUE))
  columns<-matrix(unlist(lapply(features,function(f)vapply(members,function(s)
    groups[[s]]$n_forced+match(f,groups[[s]]$feature_index)-1L,1L))),nrow=length(members))
  storage.mode(columns)<-'integer'
  list(groups=as.integer(members-1L),columns=columns,nu=as.double(nu[l]))
 })
 beta<-lapply(seq_along(groups),function(s){g<-groups[[s]]
  selected<-switch(initial,empty=rep(0,length(g$feature_index)),full=rep(1,length(g$feature_index)),alternating=as.integer((g$feature_index+s)%%2==0))
  c(rep(.1,g$n_forced),.1*selected)})
 theta<-lapply(block,function(b){x<-matrix(if(model_variant=='bms')0 else .25,length(b$groups),length(b$groups));diag(x)<-0;x})
 settings<-list(draws=as.integer(draws),burnin=as.integer(burnin),thin=1L,
   sharing=as.integer(model_variant=='imr'),keep_latent=as.integer(keep_latent),verbose=0L,
   theta_step=.4,swap_rate=.5,interaction_prior=as.double(interaction_prior))
 initial_state<-list(coefficients=beta,variance=vapply(groups,function(g)g$residual_prior[['rate']]/(g$residual_prior[['shape']]+1),0),interaction=theta)
 set.seed(seed)
 result<-.Call('imr_joint_sample',spec,block,settings,initial_state,PACKAGE='IntegMultiReg')
 names(result$coefficients)<-names(groups)
 selection<-matrix(NA_integer_,draws,n_groups*length(feature_platform))
 for(s in seq_along(groups)) {
  index<-(s-1L)*length(feature_platform)+groups[[s]]$feature_index
  selection[,index]<-result$coefficients[[s]][,groups[[s]]$n_forced+seq_along(groups[[s]]$feature_index),drop=FALSE]!=0
 }
 list(beta=result$coefficients,variance=result$variance,selection=selection,
      theta=do.call(cbind,result$interaction),latent=result$latent,log_density=result$log_density,
      diagnostics=result$acceptance)
}
