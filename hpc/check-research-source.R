args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, args[1] %in% c("capture", "verify"))
root <- normalizePath(Sys.getenv("IMR_ROOT"), mustWork = TRUE)
files <- c(list.files(file.path(root, "paper"), "[.]R$", full.names = TRUE),
  list.files(file.path(root, "hpc"), "[.](R|sh|py)$", full.names = TRUE),
  file.path(root, c("paper/reference/original-generator.c",
                    "paper/data/kirc_table1_full.rda",
                  "Ref/biom12587-sup-0002-suppdata_code.zip")))
hashes <- tools::md5sum(sort(files))
if (args[1] == "capture") {
  stopifnot(!file.exists(args[2]), !anyNA(hashes))
  saveRDS(hashes, args[2])
} else stopifnot(identical(hashes, readRDS(args[2])))
