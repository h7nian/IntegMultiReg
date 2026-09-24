# Test adapter for the frozen serial/split fixtures; no collection logic here.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(dirname(script), "collection.R"))
settings <- readRDS(file.path(args[2], "001/settings.rds"))
n <- if (args[1] == "chains") 8L else settings$replicates * if (args[1] == "simulation") 3L else 6L
collect_tasks(args[1], file.path(args[2], sprintf("%03d", seq_len(n))), args[3])
