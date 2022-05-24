# Package names
packages <- c("mlbench", "mnormt", "VGAM", "miscTools", "monomvn", "SSLASSO", "glmnet", "flare", "resample")
## monomvn package - blasso function
## docs:  https://cran.r-project.org/web/packages/monomvn/monomvn.pdf

# Packages loading
suppressMessages(invisible(lapply(packages, library, character.only = TRUE)))
error_func <- function(){
    # Handles errors that arise
    traceback(2, max.lines=10)
    q()
}
options(error=error_func, keep.source=TRUE)
source("main_func.r")

#DATA OPTIONS: "bh", "bb_syndrome", "synthetic"
#  Must provide filename if "synthetic"
# dataset_name = "bb_syndrome" 
# dataset_name = "bh" 
dataset_name = "synthetic" 
dataset_fname = "../dep_data_high2.RData"

reps = 100
main(dataset_name, dataset_fname, reps=reps)
