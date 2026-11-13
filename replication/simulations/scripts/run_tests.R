here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                      value = TRUE)[1])))
here <- dirname(here)
source(file.path(here, "R", "execution.R"))
load_kgmselect()
res <- tinytest::run_test_dir(file.path(here, "tests", "R"), verbose = 0)
print(res)
if (!all(as.logical(res))) quit(status = 1L)
