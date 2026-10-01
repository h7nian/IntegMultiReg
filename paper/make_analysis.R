# Compatibility entry point; the maintained manuscript workflow lives alongside.
# Pass the same arguments as IntegMultiReg-replication.R (for example --quick).
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                   value = TRUE)[1])
entry <- file.path(dirname(normalizePath(script)),
                   "IntegMultiReg-replication.R")
status <- system2(file.path(R.home("bin"), "Rscript"),
  c(shQuote(entry), shQuote(commandArgs(TRUE))))
quit(status = status)
