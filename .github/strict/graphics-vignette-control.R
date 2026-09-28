# Match vignette rendering without loading IntegMultiReg or using its data.
stopifnot(!'IntegMultiReg' %in% loadedNamespaces())
probe <- Sys.getenv('IMR_GRAPHICS_PROBE', 'default')
stopifnot(probe %in% c('default', 'width'))
cat('Package-free graphics probe:', probe, '\n')
first_graphic <- if (probe == 'width') c(
  'plot.new()',
  'invisible(strwidth(c("Group A", "Group B", "Iteration", "Value")))'
) else character()
work <- tempfile('base-graphics-vignette-'); dir.create(work)
input <- file.path(work, 'control.Rmd')
writeLines(c('---', 'title: "Base graphics control"', 'output: html_document', '---',
  '```{r setup, include=FALSE}',
  'stopifnot(!"IntegMultiReg" %in% loadedNamespaces())',
  '```',
  '```{r graphics}',
  first_graphic,
  'for (i in 1:12) {',
  '  plot(1:20, sin((1:20)/i), main=paste("Font and axis control", i),',
  '       xlab="Iteration", ylab="Value", type="l")',
  '  legend("topright", c("Group A", "Group B"), lty=1:2)',
  '}',
  '```'), input)
rmarkdown::render(input, quiet=TRUE, envir=new.env(parent=globalenv()))
stopifnot(!'IntegMultiReg' %in% loadedNamespaces())
gc()
