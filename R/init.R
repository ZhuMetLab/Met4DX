.set_package_options <- function(pkgname, run_env) {

    ## add met4dx specific options
    ## (not unlike what is done in 'affy')
    if (is.null(getOption("BioC"))) {
        BioC <- list()
        class(BioC) <- "BioCOptions"
        options("BioC"=BioC)
    }


    pkg_option <- list('run_env' = run_env)

    class(pkg_option) <- "BioCPkg"

    BioC <- getOption("BioC")
    BioC[[pkgname]] <- pkg_option
    options("BioC"=BioC)
}