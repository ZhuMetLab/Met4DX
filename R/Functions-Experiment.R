# project: Met4DX
# File name: Functions-Experiment.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/1/17 13:09
# Copyright (c) 2022- ZhuMSLab ALL right reserved

.get_sample_groups <- function(paths, keep_suffix = FALSE) {
  # from xcms (phenoDataFromPaths)
  ## create factors from filesystem hierarchy
  sclass <- gsub("^\\.$", "sample", dirname(paths))
  lev <- strsplit(sclass, "/")
  levlen <- sapply(lev, length)
  if(length(lev) > 1 && !all(levlen[1] == levlen))
    stop("Directory tree must be level")
  pdata <- as.data.frame(matrix(unlist(lev), nrow=length(lev), byrow=TRUE))
  redundant <- apply(pdata, 2, function(col) length(unique(col)) == 1)
  if (!any(!redundant)) {
    redundant[length(redundant)] <- FALSE
  }
  pdata <- pdata[,!redundant,drop=FALSE]
  if (ncol(pdata) == 1) { ## if not multiple factors, behave as before
    ## Make the default group names less redundant
    scomp <- strsplit(substr(sclass, 1, min(nchar(sclass))), "")
    scomp <- matrix(c(scomp, recursive = TRUE), ncol = length(scomp))
    i <- 1
    while(all(scomp[i,1] == scomp[i,-1]) && i < nrow(scomp))
      i <- i + 1
    i <- min(i, tail(c(0, which(scomp[1:i,1] == .Platform$file.sep)), n = 1) + 1)
    if (i > 1 && i <= nrow(scomp))
      sclass <- substr(sclass, i, max(nchar(sclass)))
    pdata <- data.frame(factor(sclass))
    colnames(pdata) <- "class"
  }
  if (keep_suffix) {
    rownames(pdata) <- basename(paths)
  } else {
    rownames(pdata) <- gsub("\\.[^.]*$", "", basename(paths))
  }
  pdata
}

.sample_groups_to_class <- function(x) {
  cls <- as.character(x$class)
  factor(cls, levels = unique(cls))
}

.sample_class <- function(paths, keep_suffix = FALSE) {
  .sample_groups_to_class(.get_sample_groups(paths, keep_suffix))
}

.get_tims_files <- function(wd) {
  res <- lapply(grep("\\.d$", list.dirs(wd, recursive = TRUE), value = TRUE), function(dr) {
    if (length(list.files(dr, pattern = "^analysis.tdf$", recursive = FALSE)) == 0) {
      return(NA)
    } else {
      return(dr)
    }
  })
  return(do.call(c, res[!is.na(res)]))
}

# .get_tims_data <- function(file, colum_names) {
#   D <- opentimsr::OpenTIMS(file)
#   MS1_frame <- opentimsr::MS1(D)
#
#   precursors_table <- opentimsr::table2df(D, "Precursors")$Precursors
#   frame_id_d <- opentimsr::table2df(D, "Frames")$Frames
# }


## Experiment
#'
#' Params for MS experiment
#' @param wd \code{character} Working directory for experiment to be processed
#' @param injection_order \code{character} Injection order file name (relative to wd).
#'  Two/Three columns are required: 'file', 'order', (optional) 'sample_type'
#' @param ion_mode \code{charactor} ionization mode for MS experiment data acquisition
#' @param res_dir \code{character} directory for saving results
#' @param tmp_dir \code{character} directory for saving temporary files
#' @param ce \code{character} CE for MS experiment
#' @param lc_column \code{character} LC column for MS experiment
#' @param lc_method \code{character} LC method for MS experiment
#' @param ms1range \code{numeric(2) or NULL} Mass range setup for MS1 acquisition
#' @param ms2range \code{numeric(2) or NULL} Mass range setup for MS2 acquisition
#' @param res_define_at numeric m/z threshold for using ppm tolerance for MS1 or MS2 match,
#'  for smaller m/z values, using the tolerance of ppm.ms1 * res.defineat to mathing experimental
#'  and librarial fragments.
#' @param nSlaves \code{numeric} Number of threads to be used for multiprocessing
#' @param BPPARAM \code{BiocParallelParam or NULL} BPPARAM for multiprocessing data processing
#'
#' @export
Experiment <- function(
  wd = ".",
  experiment_type = c("DDA", "DIA"),
  injection_order = NULL,
  ion_mode = c("positive", "negative"),
  res_dir = "results",
  tmp_dir = file.path(res_dir, "tmp"),
  ce = "30",
  ms1range = NULL,
  ms2range = NULL,
  rt_range = c(0, 720),
  lc_column = c("RP", "HILIC"),
  lc_method = "RP18",
  res_define_at = 200,
  nSlaves = 4,
  BPPARAM = NULL
) {
  if (!dir.exists(wd)) {
    stop("Working directory does not exist")
  }
  experiment_type <- match.arg(experiment_type)
  ion_mode <- match.arg(ion_mode)
  lc_column <- match.arg(lc_column)
  lc_method <- match.arg(lc_method)

  options(mc.cores = nSlaves)

  if (!dir.exists(res_dir)) {
    dir.create(res_dir)
  }

  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir)
  }

  if (!missing(injection_order)) {
    if (!is.null(injection_order) & !file.exists(file.path(wd, injection_order))) {
      stop('File not exist: ', injection_order)
    }
  }

  if (missing(BPPARAM) || is.null(BPPARAM)) {
    os <- Sys.info()["sysname"]
    BPPARAM <- switch(os,
                      "Darwin" = {
                        BiocParallel::MulticoreParam(workers = nSlaves, progressbar = TRUE)
                      },
                      "Linux" = {
                        BiocParallel::SnowParam(workers = nSlaves, progressbar = TRUE)
                      },
                      "Windows" = {
                        BiocParallel::SnowParam(workers = nSlaves, progressbar = TRUE)
                      }
    )
  }

  new(paste0(experiment_type, "Experiment"),
      wd = wd,
      injection_order = injection_order,
      ion_mode = ion_mode,
      res_dir = res_dir,
      tmp_dir = tmp_dir,
      ce = ce,
      ms1range = ms1range,
      ms2range = ms2range,
      rt_range = rt_range,
      lc_column = lc_column,
      lc_method = lc_method,
      res_define_at = res_define_at,
      BPPARAM = BPPARAM
  )
}


setMethod(
  "show",
  signature = c("Experiment"),
  function(object) {
    cat("A", class(object), "Object\n")
})