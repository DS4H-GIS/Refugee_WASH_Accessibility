set_script_dir_as_wd <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    script_path <- sub("--file=", "", args[grep("--file=", args)])
    
    if (length(script_path) == 0) {
        stop("Cannot determine script path. Are you running this via Rscript?")
    }
    
    script_dir <- dirname(script_path)
    setwd(script_dir)
    message("Working directory set to script location: ", script_dir)
}