source('exact-reference.R')
gaussian_availability_reference <- function(groups, feature_platform, nu,
 interaction_prior=c(shape=2,rate=2), model_variant=c('imr','bms')) {
 model_variant<-match.arg(model_variant)
 G<-length(groups);P<-length(feature_platform);L<-length(nu)
 available<-matrix(FALSE,G,P)
 for(g in seq_len(G))available[g,groups[[g]]$feature_index]<-TRUE
 free<-which(as.vector(t(available)))
 stopifnot(length(free)<=12)
 states<-as.matrix(expand.grid(rep(list(0:1),length(free))))
 cache<-lapply(groups,function(dat) {
  local<-as.matrix(expand.grid(rep(list(0:1),length(dat$feature_index))))
  lapply(seq_len(nrow(local)),function(i){active<-c(seq_len(dat$n_forced),dat$n_forced+which(local[i,]==1))
   list(active=active,moments=gaussian_pmom_evidence(dat$X[,active,drop=FALSE],dat$y,dat$prior_scale[active],dat$residual_prior))})
 })
 members<-lapply(seq_len(L),function(l)which(rowSums(available[,feature_platform==l,drop=FALSE])>0))
 stopifnot(all(lengths(members)<=2))
 estimated_theta<-which(lengths(members)==2)
 theta<-matrix(0,nrow(states),length(estimated_theta))
 beta_mean<-beta_second<-lapply(groups,function(g)matrix(0,nrow(states),ncol(g$X)))
 variance_mean<-matrix(0,nrow(states),G)
 log_weight<-numeric(nrow(states))
 for(row in seq_len(nrow(states))) {
  gamma<-matrix(0,G,P);flat<-as.vector(t(gamma));flat[free]<-states[row,];gamma<-matrix(flat,G,byrow=TRUE)
  prior<-0
  for(l in seq_len(L)) {
   block<-gamma[members[[l]],feature_platform==l,drop=FALSE]
   if(model_variant=='bms'||nrow(block)==1) {
    prior<-prior+nu[l]*sum(block)-length(block)*log1p(exp(nu[l]))
   } else {
    density<-function(t) {
     terms<-cbind(0,nu[l],nu[l],2*nu[l]+2*t)
     maxv<-apply(terms,1,max);logz<-maxv+log(rowSums(exp(terms-maxv)))
     exp(nu[l]*sum(block)+2*t*sum(block[1,]*block[2,])-ncol(block)*logz+
      dgamma(t,interaction_prior[['shape']],rate=interaction_prior[['rate']],log=TRUE))
    }
    mass<-integrate(density,0,Inf,rel.tol=1e-10,subdivisions=500)$value
    prior<-prior+log(mass)
    theta[row,match(l,estimated_theta)]<-integrate(function(t)t*density(t),0,Inf,rel.tol=1e-10,subdivisions=500)$value/mass
   }
  }
  log_weight[row]<-prior
  for(g in seq_len(G)) {
   index<-1L+sum(gamma[g,groups[[g]]$feature_index]*2^(seq_along(groups[[g]]$feature_index)-1L))
   item<-cache[[g]][[index]];m<-item$moments
   log_weight[row]<-log_weight[row]+m$log_evidence
   beta_mean[[g]][row,item$active]<-m$beta_mean
   beta_second[[g]][row,item$active]<-m$beta_second
   variance_mean[row,g]<-m$variance_mean
  }
 }
 probability<-exp(log_weight-max(log_weight));probability<-probability/sum(probability)
 list(free=free,states=states,probability=probability,inclusion=drop(crossprod(probability,states)),
      beta_mean=lapply(beta_mean,function(x)drop(crossprod(probability,x))),
      beta_second=lapply(beta_second,function(x)drop(crossprod(probability,x))),
      variance_mean=drop(crossprod(probability,variance_mean)),theta_mean=drop(crossprod(probability,theta)))
}
