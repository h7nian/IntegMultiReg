# Match vignette rendering without loading IntegMultiReg or using its data.
stopifnot(!'IntegMultiReg' %in% loadedNamespaces())
work <- tempfile('base-graphics-vignette-'); dir.create(work)
input <- file.path(work, 'control.Rmd')
writeLines(c('---', 'title: "Base graphics control"', 'output: html_document', '---',
  '```{r setup, include=FALSE}',
  'stopifnot(!"IntegMultiReg" %in% loadedNamespaces())',
  '```',
  '```{r graphics}',
  'for (i in 1:12) {',
  '  plot(1:20, sin((1:20)/i), main=paste("Font and axis control", i),',
  '       xlab="Iteration", ylab="Value", type="l")',
  '  legend("topright", c("Group A", "Group B"), lty=1:2)',
  '}',
  '```'), input)
rmarkdown::render(input, quiet=TRUE, envir=new.env(parent=globalenv()))
stopifnot(!'IntegMultiReg' %in% loadedNamespaces())
gc()
