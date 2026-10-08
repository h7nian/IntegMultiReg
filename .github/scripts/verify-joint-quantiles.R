# Independent marginal CDFs for the small continuous joint posterior.
# Fixing one active coefficient t changes the IG prior to
# shape a+3/2 and rate b+t^2/(2*tau), with an explicit normalizer.
source('.github/reference/exact-reference.R')
root <- file.path(Sys.getenv('RUNNER_TEMP', tempdir()), 'joint-reference-evidence')
groups<-readRDS('.github/reference/fixtures/continuous-data.rds')
reference<-readRDS('.github/reference/fixtures/continuous-reference.rds')
chains<-readRDS(file.path(root, 'continuous/chains.rds'))
local_states<-as.matrix(expand.grid(marker1=0:1,marker2=0:1))
local_weights<-lapply(1:2,function(g) {
 pattern<-1+drop(reference$states[,2*(g-1)+1:2,drop=FALSE]%*%c(1,2))
 vapply(1:4,function(i)sum(reference$probability[pattern==i]),0)
})
log_evidence<-function(X,y,tau,prior) {
 if(ncol(X)) return(gaussian_pmom_evidence(X,y,tau,prior)$log_evidence)
 a<-prior[['shape']];b<-prior[['rate']];n<-length(y)
 -n/2*log(2*pi)+a*log(b)-lgamma(a)+lgamma(a+n/2)-(a+n/2)*log(b+sum(y^2)/2)
}
models<-lapply(1:2,function(g)lapply(1:4,function(i) {
 d<-groups[[g]];active<-c(1,1+which(local_states[i,]==1));X<-d$X[,active,drop=FALSE]
 tau<-d$prior_scale[active];a<-d$residual_prior[['shape']];b<-d$residual_prior[['rate']]
 A<-crossprod(X)+diag(1/tau,nrow=length(tau));covariance<-solve(A)
 mu<-drop(covariance%*%crossprod(X,d$y));shape<-a+nrow(X)/2+ncol(X)
 rate<-b+(sum(d$y^2)-sum(drop(crossprod(X,d$y))*mu))/2
 polynomial<-normal_moment_polynomial(rep(seq_along(tau),each=2L),mu,covariance)
 k<-which(polynomial!=0);power<-k-1
 logw<-log(abs(polynomial[k]))+lgamma(shape-power)-(shape-power)*log(rate)
 w<-sign(polynomial[k])*exp(logw-max(logw));w<-w/sum(w)
 list(active=active,X=X,y=d$y,tau=tau,a=a,b=b,
  log_evidence=log_evidence(X,d$y,tau,d$residual_prior),
  variance_cdf=function(v)sum(w*pgamma(rate/v,shape=shape-power,lower.tail=FALSE)))
}))
beta_density<-function(t,model,j) {
 if(t==0)return(0)
 pos<-match(j,model$active);tau<-model$tau[pos]
 a2<-model$a+1.5;b2<-model$b+t^2/(2*tau)
 shift<-model$y-model$X[,pos]*t
 log_constant<-2*log(abs(t))-.5*log(2*pi)-1.5*log(tau)+
  model$a*log(model$b)-lgamma(model$a)+lgamma(a2)-a2*log(b2)
 exp(log_constant+log_evidence(model$X[,-pos,drop=FALSE],shift,model$tau[-pos],
  c(shape=a2,rate=b2))-model$log_evidence)
}
beta_cdf<-function(q,g,j,left=FALSE) {
 sum(vapply(1:4,function(i) {
  model<-models[[g]][[i]]
  probability<-if(!j%in%model$active) as.numeric(if(left)q>0 else q>=0) else
   integrate(function(t)vapply(t,beta_density,0,model=model,j=j),-Inf,q,
     rel.tol=1e-7,abs.tol=1e-9,subdivisions=300L)$value
  local_weights[[g]][i]*probability
 },0))
}
# Check each active coefficient density integrates to one independently of its CDF.
normalizations<-list()
for(g in 1:2)for(i in 1:4)for(j in models[[g]][[i]]$active) {
 m<-models[[g]][[i]]
 value<-integrate(function(t)vapply(t,beta_density,0,model=m,j=j),-Inf,Inf,
  rel.tol=1e-7,abs.tol=1e-9,subdivisions=300L)$value
 stopifnot(abs(value-1)<1e-6)
 normalizations[[length(normalizations)+1L]]<-data.frame(group=g,model=i,coefficient=j,integral=value)
}
rows<-list();probs<-c(.025,.5,.975)
record<-function(parameter,draws,cdf) {
 qs<-quantile(draws,probs=probs,names=FALSE,type=1)
 for(k in seq_along(probs)) {
  q<-qs[k];lower<-cdf(q,TRUE);upper<-cdf(q,FALSE)
  error<-max(0,lower-probs[k],probs[k]-upper)
  rows[[length(rows)+1L]]<<-data.frame(parameter,probability=probs[k],empirical_quantile=q,
   exact_cdf_left=lower,exact_cdf_right=upper,probability_error=error,passed=error<.01)
 }
}
for(g in 1:2) {
 for(j in 1:3) record(paste0('beta:',g,':',j),unlist(lapply(chains,function(x)x$beta[[g]][,j])),
  function(q,left)beta_cdf(q,g,j,left))
 record(paste0('variance:',g),unlist(lapply(chains,function(x)x$variance[,g])),
  function(q,left)sum(vapply(1:4,function(i)local_weights[[g]][i]*models[[g]][[i]]$variance_cdf(q),0)))
}
# Theta CDF: integrate the normalized MRF prior at each complete selection state.
prior_density<-function(theta,state) {
 selection<-matrix(state,2,byrow=TRUE);nu<--1
 logs<-cbind(0,nu,nu,2*nu+2*theta);largest<-apply(logs,1,max)
 logZ<-largest+log(rowSums(exp(logs-largest)))
 exp(nu*sum(selection)+2*theta*sum(selection[1,]*selection[2,])-2*logZ+
  dgamma(theta,shape=2,rate=2,log=TRUE))
}
prior_mass<-apply(reference$states,1,function(state)integrate(prior_density,0,Inf,state=state,rel.tol=1e-9)$value)
record('theta',unlist(lapply(chains,`[[`,'theta')),function(q,left)
 sum(vapply(seq_len(nrow(reference$states)),function(i)reference$probability[i]*
  integrate(prior_density,0,q,state=reference$states[i,],rel.tol=1e-9)$value/prior_mass[i],0)))
result<-do.call(rbind,rows)
write.csv(result,file.path(root, 'credible-quantile-reference.csv'),row.names=FALSE)
write.csv(do.call(rbind,normalizations),file.path(root, 'coefficient-density-normalizations.csv'),row.names=FALSE)
print(result,row.names=FALSE)
stopifnot(all(result$passed))
cat('All',nrow(result),'quantile checks passed; maximum CDF-probability discrepancy',max(result$probability_error),'\n')
