# Standalone external-library diagnosis: never loads IntegMultiReg or knitr.
stopifnot(!"IntegMultiReg" %in% loadedNamespaces())
for (i in 1:12) {
  png(tempfile(fileext=".png"), width=672, height=403, type="cairo", res=96)
  par(mfrow=c(2,2))
  plot(1:30, sin(1:30), type="l", main="Posterior log survival", xlab="Iteration", ylab="Value", font.main=1+(i%%4))
  barplot(1:12, names.arg=paste0("clinical:marker",1:12), las=2, main="Selection probabilities")
  plot(1:10, pch=1:10, main="Model diagnostics", sub="Intercept and coefficients")
  text(1:10,1:10,labels=c("(Intercept)","age","sex","stage","G01","G02","M01","M02","theta","95%"),font=1+(i%%4))
  legend("topleft",c("Clinical","Molecular","Future outcome"),lty=1:3)
  plot(1:10, log(1:10), main=expression(theta+beta),xlab="Log time",ylab="Outcome")
  dev.off()
}
gc()
stopifnot(!"IntegMultiReg" %in% loadedNamespaces())
print(loadedNamespaces())
cat("Completed standalone graphics; no IntegMultiReg code was loaded.\n")
