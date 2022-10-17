# project: Met4DX
# File name: zzz.R.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/1/17 9:51
# Copyright (c) 2022- ZhuMSLab ALL right reserved

#' MetAnalyzer2_QE
#'
#' Processing MS data
#'
#' This package processes MS data for peak detection and MSMS spectra assignment
#' to generate feature table and corresponding MSMS spectra.
#' @docType package
#' @author Mingdu Luo(luomd@sioc.ac.cn), Yandong Yin (\email{yinyandong@@sioc.ac.cn})
#' @import SpectraTools opentimsr BiocParallel ggpmisc smoother Rcpp splus2R fastmatch
#' @importFrom Rcpp evalCpp
#' @useDynLib Met4DX
#' @name Met4DX
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\n========== WELCOME to ",
                        getPackageName(), " v", packageVersion(getPackageName()),
                        " ===========",
                        "\nDeveloped by zhulab for MS Data Processing",
                        "\nImported packages:",
                        "\n  SpectraTools v", packageVersion("SpectraTools"),
                        "\n  OpenTIMSR v", packageVersion("opentimsr"),
                        "\n=============================================\n")
}

.onLoad <- function(libname, pkgname) {
  # require(methods)
  .set_package_options(pkgname, 'production')
  # .set_package_options(pkgname, 'development')
  getOption('BioC')[[pkgname]]
}
