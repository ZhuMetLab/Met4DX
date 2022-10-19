# project: Met4DX
# File name: Methods-TimsData.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/1/17 10:45
# Copyright (c) 2022- ZhuMSLab ALL right reserved

#' @export
setMethod(
  "show",
  signature = "TimsData",
  function(object) {
    cat("An \"TimsData\" object with", length(object@files), "samples.\n")
  }
)

#' @export
setMethod(
  "ReadSpectraData",
  signature = c("TimsData", "ReadSpectraParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    message("Reading spectra data...")
    object@.processHistory <- c(object@.processHistory, param)
    param_list <- as.list(param)
    par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)

    files <- .tmp_files(object@files, object@experiment@tmp_dir, 'spectra')
    names(files) <- object@files

    spec_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(object@files[idxs], function(data_file) {
        spec_file <- list.files(data_file, pattern = "(?i)\\.mgf$", recursive = FALSE, full.names = TRUE)
        if (length(spec_file) == 0) {
          stop("No mgf file found in " + dr)
        } else if (length(spec_file) > 1) {
          warning("More than one mgf file found in " +
                    dr +
                    "; using the first one")
        } else {
          spec_file <- spec_file[1]
        }
        c("data_file" = spec_file,
          "res_define_at" = object@experiment@res_define_at,
          param_list)
      })
      .parallel_parser(".read_spectra", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))

    names(spec_files) <- object@files
    object@tmp_data_files$spectra_files <- spec_files

    setwd(wd0)
    return(object)
  })

#' @export
setMethod(
  "BinPrecursors",
  signature = c("TimsData", "BinPrecursorParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    message("Binning precursor data...")

    param_list <- as.list(param)
    par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)

    files <- .tmp_files(object@files, object@experiment@tmp_dir, 'precursor_bins')
    names(files) <- object@files

    precursor_bin_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(object@files[idxs], function(data_file) {
        c("data_file" = unname(object@tmp_data_files$spectra_files[data_file]),
          "res_define_at" = object@experiment@res_define_at,
          "ms2_range" = object@experiment@ms2range,
          param_list)
      })
      .parallel_parser(".bin_precursors", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))

    names(precursor_bin_files) <- object@files
    object@tmp_data_files$precursor_bin_files <- precursor_bin_files

    setwd(wd0)
    return(object)
  })


#' @export
setMethod(
  "QueryTimsData",
  signature = c("TimsData", "QueryTimsDataParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    
    message("Querying tims data...")
    
    param_list <- as.list(param)
    par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)
    
    files <- .tmp_files(object@files, object@experiment@tmp_dir, 'tims_data')
    names(files) <- object@files
    
    tims_data_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(object@files[idxs], function(data_file) {
        c("data_file" = data_file,
          param_list)
      })
      .parallel_parser(".query_tims_data", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))
    names(tims_data_files) <- object@files
    object@tmp_data_files$tims_data_files <- tims_data_files
    
    setwd(wd0)
    return(object)
  })



#' @export
setMethod(
  "ExtractIMData",
  "signature" = c("TimsData", "ExtractIMDataParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    message("Extracting IM data...")

    sub_dir <- ifelse(is.null(param@interpolate_method), "im_data",
                      file.path("im_data", param@interpolate_method))

    sub_dir <- file.path(sub_dir, paste0('smooth_', param@smooth_method))


    param_list <- as.list(param)
    par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)

    files <- .tmp_files(object@files, object@experiment@tmp_dir, sub_dir)
    names(files) <- object@files

    im_data_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(object@files[idxs], function(data_file) {
        c(list("tims_data_file" = unname(object@tmp_data_files$tims_data_files[data_file]),
               "info" = readRDS(object@tmp_data_files$spectra_files[data_file])@info,
               "precursor_bin_file" = unname(object@tmp_data_files$precursor_bin_files[data_file]),
               "res_define_at" = object@experiment@res_define_at),
          param_list)
      })
      .parallel_parser(".extract_im_data", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))

    names(im_data_files) <- object@files
    object@tmp_data_files$im_data_files <- im_data_files

    files <- .tmp_files(paste0(object@files, '_peaks'), object@experiment@tmp_dir, sub_dir)
    names(files) <- object@files

    peak_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(idxs, function(smp_idx) {
        data_file <- object@files[smp_idx]
        c("spectra_file" = unname(object@tmp_data_files$spectra_files[data_file]),
          "im_data_file" = unname(object@tmp_data_files$im_data_files[data_file]),
          "smp_idx" = smp_idx,
          "tmp_dir" = unname(files[data_file]),
          param_list)
      })
      .parallel_parser(".get_peaks", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))
    object@peaks <- do.call(rbind, lapply(peak_files, function(peak_file) {
      peaks <- readRDS(peak_file)
      rownames(peaks) <- NULL
      peaks
    }))
    rownames(object@peaks) <- .gen_indexes(object@peaks)
    setwd(wd0)
    return(object)
  })

#' @export
setMethod(
  "DereplicatePeaks",
  signature = c("TimsData", "DereplicatePeaksParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    message("Dereplicating peaks...")

    sub_dir <- file.path('dereplicate_peaks', ifelse(param@match_msms, 'w_msms_match', 'wo_msms_match'))
    param_list <- as.list(param)

    par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)

    files <- .tmp_files(object@files, object@experiment@tmp_dir, sub_dir)
    names(files) <- object@files

    derep_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(idxs, function(idx) {
        peaks <- object@peaks[object@peaks$smp_idx == idx, , drop = FALSE]
        spectra_data <- readRDS(object@tmp_data_files$spectra_files[idx])
        precursor_bins <- readRDS(object@tmp_data_files$precursor_bin_files[idx])
        precursor_bins <- lapply(precursor_bins, `[`, 1)[peaks$spec_idx]
        names(precursor_bins) <- rownames(peaks)
        c(list("peaks" = peaks,
               "spectra_data" = spectra_data,
               "precursor_bins" = precursor_bins,
               "res_define_at" = object@experiment@res_define_at,
               "ms2_range" = object@experiment@ms2range),
          param_list)
      })
      .parallel_parser(".dereplicate_peaks", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))

    names(derep_files) <- object@files
    object@tmp_data_files$derepicate_files <- derep_files

    files <- .tmp_files(paste0(object@files, "_spectra"),
                        object@experiment@tmp_dir,
                        sub_dir)

    object@peaks <- do.call(rbind, lapply(seq_along(derep_files), function(idx) {
      # browser()
      precursor_bins <- readRDS(derep_files[idx])
      spectra_data <- readRDS(object@tmp_data_files$spectra_files[idx])
      peaks <- object@peaks[object@peaks$smp_idx == idx, , drop = FALSE]
      peaks <- peaks[rownames(peaks) %in% names(precursor_bins), , drop = FALSE]
      nms <- rownames(peaks)
      spectra_data <- spectra_data[sapply(precursor_bins[nms], `[`, 1)]
      spectra_data <- UpdateNames(spectra_data, nms)
      saveRDS(spectra_data, file = files[idx], version = 2)
      peaks
    }))
    object@tmp_data_files$spectra_files_dereplicaters <- files
    setwd(wd0)
    return(object)
  }
)

#' @export
setMethod(
  "AlignPeaks",
  signature = c("TimsData", "AlignPeakParamRaw", "ANY", "ANY", "ANY"),
  function(object, param, ref_index = NULL, ref_sample = NULL, num_pools = 1) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    message("Aligning and grouping peaks...")

    if (missing(ref_index)) {
      ref_index <- round(length(object@files) / 2)
    }

    if (!missing(ref_sample)) {
      ref_index <- match(ref_sample, basename(object@files))
      if (is.na(ref_index)) {
        stop("ref_sample '", ref_sample, "' was not found!")
      }
    }

    files <- .tmp_files(paste("alignment_result", 'ref', ref_index, sep = '_'),
                        object@experiment@tmp_dir, 'peak_alignment')
    param_list <- as.list(param)

    arg_list <- c(list("ref_index" = ref_index,
                       "peaks" = object@peaks,
                       "spectra_files" = object@tmp_data_files$spectra_files,
                       "precursor_bin_file" = object@tmp_data_files$precursor_bin_files,
                       "rt_range" = object@experiment@rt_range,
                       "ms1_range" = object@experiment@ms1range,
                       "ms2_range" = object@experiment@ms2range,
                       "res_define_at" = object@experiment@res_define_at),
                  param_list)
    align_file <- .analysis_parser(".align_peaks", arg_list, files)
    object@tmp_data_files$align_file <- align_file
    tmp <- readRDS(align_file)
    object@features <- tmp$features
    object@peak_groups <- tmp$peak_groups
    setwd(wd0)
    return(object)
  })

#' @export
setMethod(
  "AlignPeaks",
  signature = c("TimsData", "AlignPeakParam", "ANY", "ANY", "ANY"),
  function(object, param, ref_index = NULL, ref_sample = NULL, num_pools = NULL) {
    print('pool')
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    param_list <- as.list(param)

    if (missing(ref_index)) {
      ref_index <- round(length(object@files) / 2)
    }

    if (!missing(ref_sample)) {
      ref_index <- match(ref_sample, basename(object@files))
      if (is.na(ref_index)) {
        stop("ref_sample '", ref_sample, "' was not found!")
      }
    }

    if (missing(num_pools)) {
      num_pools <- ceiling(length(object@files) / object@experiment@BPPARAM$workers)
    }
    query_indexes <- seq_along(object@files)[-ref_index]
    peaks <- object@peaks
    if (param@strict_rt_constrains) {
      peaks <- peaks[peaks$rt >= object@experiment@rt_range[1] & peaks$rt <= object@experiment@rt_range[2], ,
                     drop = FALSE]
    }
    peaks$ccs_tol <- .percentage2value(peaks$ccs, param@ccs_tol)
    if (param@mz_tol > 1) {
      peaks$mz_tol <- sapply(peaks$mz, function(mz) .ppm2dalton(mz, param@mz_tol, object@experiment@res_define_at))
    } else {
      peaks$mz_tol <- param@mz_tol
    }
    peaks$rt_align <- peaks$rt

    spectra_all <- .get_peak_spectra(peaks, object@tmp_data_files$spectra_files,
                                     param@min_num_fragments, update_info = TRUE)
    files <- .tmp_files(paste(object@files, 'ref', ref_index, sep = '_'),
                        object@experiment@tmp_dir, 'rt_correction')
    names(files) <- object@files
    par_idx <- .gen_parallel_indexes(length(object@files), ceiling(length(query_indexes) / num_pools))
    # idx_from <- seq(from = 1, to = length(query_indexes), by = ceiling(length(query_indexes) / num_pools))
    # idx_to <- c(idx_from[-1] - 1, length(query_indexes))
    rt_correction_data <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', query_indexes[idxs], '\n')
      arg_list <- lapply(query_indexes[idxs], function(query_index) {
        c(list("peak_ref" = peaks[peaks$smp_idx == ref_index, , drop = FALSE],
               "peak_query" = peaks[peaks$smp_idx == query_index, , drop = FALSE],
               "spectra" = spectra_all[names(spectra_all) %in% rownames(peaks[peaks$smp_idx %in% c(ref_index,
                                                                                                   query_index),
                                                                          , drop = FALSE])],
               "rt_range" = object@experiment@rt_range,
               "rt_tol" = param@rt_tol_landmark,
               "res_define_at" = object@experiment@res_define_at),
          param_list)
      })
      res <- .parallel_parser(".correct_rt_landmarks", arg_list, files[query_indexes[idxs]],
                              object@experiment@BPPARAM, TRUE, TRUE)
      names(res) <- object@files[query_indexes[idxs]]
      return(res)
    }, simplify = FALSE))
    # arg_list <- lapply(query_indexes, function(query_index) {
    #   c(list("peak_ref" = peaks[peaks$smp_idx == ref_index, , drop = FALSE],
    #          "peak_query" = peaks[peaks$smp_idx == query_index, , drop = FALSE],
    #          "spectra" = spectra_all[names(spectra_all) %in% rownames(peaks[peaks$smp_idx %in% c(ref_index,
    #                                                                                              query_index),
    #                                                                     , drop = FALSE])],
    #          "rt_range" = object@experiment@rt_range,
    #          "rt_tol" = param@rt_tol_landmark,
    #          "res_define_at" = object@experiment@res_define_at),
    #     param_list)
    # })
    # names(arg_list) <- names(files[query_indexes])
    # rt_correction_data <- .parallel_parser(".correct_rt_landmarks", arg_list, files,
    #                                        object@experiment@BPPARAM, TRUE, TRUE)

    for (smp in names(rt_correction_data)) {
      peaks$rt_align[peaks$smp_idx == match(smp, names(files))] <- rt_correction_data[[smp]]$rt
    }

    files <- .tmp_files(paste('pool_data', 'ref', ref_index, num_pools, 'pools', seq(num_pools), sep = '_'),
                        object@experiment@tmp_dir, 'peak_alignment_pools')

    # idx_from <- seq(from = 1, to = length(query_indexes), by = ceiling(length(query_indexes) / num_pools))
    # idx_to <- c(idx_from[-1] - 1, length(query_indexes))
    arg_list <- apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      c(list("ref_index" = ref_index,
             "query_index" = query_indexes[idxs],
             "peaks" = peaks,
             "query_land_marks" = lapply(rt_correction_data[idxs], `[[`, 'land_marks'),
             "query_peak_cads" = lapply(rt_correction_data[idxs], `[[`, 'peak_cads'),
             "rt_tol" = param@rt_tol_match
      ),
        param_list)
    }, simplify = FALSE)

    raw_workers <- object@experiment@BPPARAM$workers
    bpparam_tmp <- object@experiment@BPPARAM
    bpparam_tmp$workers <- min(num_pools, raw_workers)

    peak_align_data <- .parallel_parser(".align_peaks_raw", arg_list, files, bpparam_tmp, TRUE, TRUE)

    files <- .tmp_files(paste('pool_alignment', 'ref', ref_index, num_pools, 'pools', sep = '_'),
                        object@experiment@tmp_dir, 'peak_alignment_final')

    arg_list <- c(list("ref_index" = ref_index,
                       "ref_names" = rownames(peaks[peaks$smp_idx == ref_index, , drop = FALSE]),
                       "peak_align_data" = peak_align_data,
                       "rt_tol" = param@rt_tol_match),
                  param_list)
    res <- .analysis_parser(".align_peaks_pool", arg_list, files, TRUE)
    object@tmp_data_files$align_file <- files
    object@features <- res$peaks
    object@peak_groups <- res$peak_groups

    object@experiment@BPPARAM$workers <- raw_workers
    setwd(wd0)
    return(object)
  })

#' @export
setMethod(
  "CorrectRT",
  signature = c("TimsData", "CorrectRTLandmarksParam", "ANY", "ANY"),
  function(object, param, ref_index = NULL, ref_sample = NULL) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    param_list <- as.list(param)

    if (missing(ref_index)) {
      ref_index <- round(length(object@files) / 2)
    }

    if (!missing(ref_sample)) {
      ref_index <- match(ref_sample, basename(object@files))
      if (is.na(ref_index)) {
        stop("ref_sample '", ref_sample, "' was not found!")
      }
    }

    message("Correcting RT with reference sample: ", object@files[[ref_index]])

    query_indexes <- seq_along(object@files)[-ref_index]
    peaks <- object@peaks
    if (param@strict_rt_constrains) {
      peaks <- peaks[peaks$rt >= object@experiment@rt_range[1] & peaks$rt <= object@experiment@rt_range[2], ,
                     drop = FALSE]
    }
    peaks$ccs_tol <- .percentage2value(peaks$ccs, param@ccs_tol)
    if (param@mz_tol > 1) {
      peaks$mz_tol <- sapply(peaks$mz, function(mz) .ppm2dalton(mz, param@mz_tol, object@experiment@res_define_at))
    } else {
      peaks$mz_tol <- param@mz_tol
    }
    peaks$rt_align <- peaks$rt
    peaks$peak_idx <- rownames(peaks)

    # spectra_all <- .get_peak_spectra(peaks, object@tmp_data_files$spectra_files_dereplicaters,
    #                                  param@min_num_fragments, update_info = TRUE)
    files <- .tmp_files(paste(object@files, 'ref', ref_index, sep = '_'),
                        object@experiment@tmp_dir, 'rt_correction/landmarks')
    names(files) <- object@files

    par_idx <- .gen_parallel_indexes(length(query_indexes), object@experiment@BPPARAM$workers)

    rt_correction_data <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', query_indexes[idxs], '\n')
      arg_list <- lapply(query_indexes[idxs], function(query_index) {
        c(list("peak_ref" = peaks[peaks$smp_idx == ref_index, , drop = FALSE],
               "peak_query" = peaks[peaks$smp_idx == query_index, , drop = FALSE],
               "spectra" = .get_peak_spectra(peaks[peaks$smp_idx %in% c(ref_index, query_index), ],
                                             object@tmp_data_files$spectra_files_dereplicaters[c(ref_index, query_index)],
                                             param@min_num_fragments, update_info = TRUE),
               "rt_range" = object@experiment@rt_range,
               "rt_tol" = param@rt_tol_landmark,
               "res_define_at" = object@experiment@res_define_at),
          param_list)
      })
      res <- .parallel_parser(".correct_rt_landmarks", arg_list, files[query_indexes[idxs]],
                              object@experiment@BPPARAM, TRUE, TRUE)
      names(res) <- object@files[query_indexes[idxs]]
      return(res)
    }, simplify = FALSE))

    for (smp in names(rt_correction_data)) {
      peaks$rt_align[peaks$smp_idx == match(smp, names(files))] <- rt_correction_data[[smp]]$rt
    }

    object@peaks$rt_align <- peaks$rt_align
    object@tmp_data_files$correct_rt_files <- files[query_indexes]

    files <- .tmp_files(paste('RTShift_ref', ref_index, sep = '_'),
                        object@experiment@res_dir)
    pdf(paste0(files, '.pdf'))
    plot(0,0,
         xlim = c(0, max(object@peaks$rt)), ylim = c(-1, 1) * max(abs(object@peaks$rt_align - object@peaks$rt)),
         type = 'n',
         xlab = 'Retention time (s)',
         ylab = 'RT shift (s)',
         main = 'RT Correction')
    abline(h=0)
    for (idx in seq(length(object@files))[-ref_index]) {
      x = tims_data@peaks[tims_data@peaks$smp_idx == idx, 'rt']
      x_diff =tims_data@peaks[tims_data@peaks$smp_idx == idx, 'rt'] - tims_data@peaks[tims_data@peaks$smp_idx == idx, 'rt_align']
      lines(sort(x), x_diff[order(x)], type = 'l', col = idx)
    }
    dev.off()
    setwd(wd0)
    return(object)
  })

#' @export
setMethod(
  "GroupPeaks",
  signature = c("TimsData", "GroupDensityParam"),
  function(object, param, dereplication_param = NULL) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    message("Grouping peaks using density method...")

    sub_dir <- file.path('group_peaks', 'density')
    files <- .tmp_files('group_data', object@experiment@tmp_dir, sub_dir)
    arglist <- c(list("peaks" = object@peaks,
                      "plot_dir" = dirname(files)),
                 as.list(param))

    res <- .analysis_parser('.group_peaks_density',
                            arglist,
                            files, TRUE)
    if (!is.null(dereplication_param)) {
      if (class(dereplication_param) != "DereplicatePeaksParam") {
        stop("'dereplication_param' must be of class DereplicatePeaksParam")
      }
      message('Dereplicating and merging overlapped features ...')
      spectra_data <- .get_peak_spectra(object@peaks, object@tmp_data_files$spectra_files_dereplicaters,
                                        update_info = TRUE)
      spectra_data <- spectra_data[sapply(res$peak_groups, `[`, 1)]
      spectra_data <- UpdateNames(spectra_data, rownames(res$features))
      arg_list <- c(list("peaks" = res$features,
                         "spectra_data" = spectra_data,
                         "precursor_bins" = res$peak_groups),
                    as.list(dereplication_param))
      files <- .tmp_files(paste0('group_data_merged_', ifelse(dereplication_param@match_msms, 'w_msms_match', 'wo_msms_match')),
                          object@experiment@tmp_dir, sub_dir)
      merged_indexes <- .analysis_parser('.dereplicate_peaks',
                                         arg_list,
                                         files, TRUE)
      col_features <- colnames(res$features)
      pk_nms <- rownames(object@peaks)
      res$peak_groups <- merged_indexes[rownames(res$features)[rownames(res$features) %in% names(merged_indexes)]]
      res$features <- as.data.frame(t(sapply(res$peak_groups, function(idx) {
        peaks <- object@peaks[fastmatch::fmatch(idx, pk_nms), , drop = FALSE]
        c(median(peaks[, "mz"]),
          range(peaks[, "mz"]),
          median(peaks[, "mobility"]),
          range(peaks[, "mobility"]),
          median(peaks[, "rt_align"]),
          range(peaks[, "rt_align"]),
          nrow(peaks),
          colMedians(peaks[, c("ccs", "target_intensity", "height", "height_fit", "area"), drop = FALSE]))
      })), stringsAsFactors = FALSE)
      colnames(res$features) <- col_features
    }

    object@.processHistory
    object@tmp_data_files$align_file <- files
    object@features <- res$features
    object@peak_groups <- res$peak_groups
    setwd(wd0)
    return(object)
  })

#' @export
setMethod(
  "MatchBetweenRuns",
  signature = c("TimsData", "MatchBetweenRunParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    message("Matching features between runs...")

    sub_dir <- ifelse(is.null(param@interpolate_method), "match_between_runs",
                      file.path("match_between_runs", param@interpolate_method))
    files <- .tmp_files(object@files, object@experiment@tmp_dir, sub_dir)
    names(files) <- object@files

    param_list <- as.list(param)
    rev_models <- readRDS(object@tmp_data_files$align_file)$rev_models
    rt_correct_files <- object@tmp_data_files$correct_rt_files

    par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)
    pks <- object@peaks[unname(do.call(c, object@peak_groups)), "smp_idx", drop = FALSE]
    pks$ft_idx <- do.call(c, mapply(rep, names(object@peak_groups), sapply(object@peak_groups, length), simplify=FALSE))
    match_between_run_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(idxs, function(idx) {
        is_fill <- !rownames(object@features) %in% pks[pks$smp_idx == idx, "ft_idx"]
        info <- object@features[is_fill, , drop = FALSE]
        if (!is.null(rev_models[[idx]])) {
          info$rt <- predict(rev_models[[idx]], info$rt_align)
        } else if (file.exists(rt_correct_file <- rt_correct_files[object@files[idx]])) {
          land_marks <- readRDS(rt_correct_file)$land_marks
          dt_rt <- data.frame('ref' = c(object@experiment@rt_range[1],
                                        object@peaks[land_marks[, 'ref'], 'rt'],
                                        object@experiment@rt_range[2]),
                              'query' = c(object@experiment@rt_range[1],
                                          object@peaks[land_marks[, 'query'], 'rt'],
                                          object@experiment@rt_range[2])
          )
          rev_model <- loess(query ~ ref, data = dt_rt[order(dt_rt$query), , drop = FALSE], span = 0.1, degree = 1L)
          info$rt <- predict(rev_model, info$rt)
        }

        c(list("tims_data_file" = unname(object@tmp_data_files$tims_data_files[idx]),
               "info" = info,
               "precursor_bin_file" = unname(object@tmp_data_files$precursor_bin_files[idx]),
               "res_define_at" = object@experiment@res_define_at),
          param_list)
      })
      .parallel_parser(".match_between_runs", arg_list, files[idxs],
                       object@experiment@BPPARAM, TRUE)
    }, simplify = FALSE))
    names(match_between_run_files) <- object@files
    object@tmp_data_files$match_between_run_files <- match_between_run_files

    files <- .tmp_files(paste0(object@files, '_peaks'), object@experiment@tmp_dir, sub_dir)
    names(files) <- object@files
    message("Finalizing matched peaks...")
    peak_files <- do.call(c, apply(par_idx, 1, function(dr) {
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      idxs <- seq(dr[1], dr[2])
      arg_list <- lapply(idxs, function(idx) {
        rt_model <- NULL
        if (file.exists(rt_correct_file <- rt_correct_files[object@files[idx]])) {
          land_marks <- readRDS(rt_correct_file)$land_marks
          dt_rt <- data.frame('ref' = c(object@experiment@rt_range[1],
                                        object@peaks[land_marks[, 'ref'], 'rt'],
                                        object@experiment@rt_range[2]),
                              'query' = c(object@experiment@rt_range[1],
                                          object@peaks[land_marks[, 'query'], 'rt'],
                                          object@experiment@rt_range[2])
          )
          rt_model <- loess(ref ~ query, data = dt_rt[order(dt_rt$ref), , drop = FALSE], span = 0.1, degree = 1L)
        }
        list("smp_index" = idx,
             "features" = object@features,
             "match_between_run_file" = unname(object@tmp_data_files$match_between_run_files[idx]),
             "rt_model" = rt_model,
             "rerun" = param@rerun)
      })
      .parallel_parser(".get_matched_peaks", arg_list, files[idxs],
                       object@experiment@BPPARAM, TRUE)
    }, simplify = FALSE))

    object@filled_peaks <- do.call(rbind, lapply(peak_files, function(peak_file) {
      peaks <- readRDS(peak_file)
      rownames(peaks) <- NULL
      peaks
    }))
    rownames(object@filled_peaks) <- .gen_indexes(object@filled_peaks, '#F')
    setwd(wd0)
    return(object)
  }
)

#' @export
setMethod(
  "FinalizeFeatures",
  signature = c("TimsData", "FinalizeFeatureParam"),
  function(object, param, par_analysis = FALSE) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    files <- .tmp_files("finalize_data", object@experiment@tmp_dir, 'finalize_data')
    message("Finalizing features...")

    param_list <- as.list(param)
    arg_list <- c(list("peaks" = object@peaks,
                       "filled_peaks" = object@filled_peaks,
                       "features" = object@features,
                       "sample_groups" = object@sample_groups,
                       "spectra_files" = object@tmp_data_files$spectra_files_dereplicaters,
                       "peak_groups" = object@peak_groups,
                       "res_define_at" = object@experiment@res_define_at,
                       "bpparam" = splus2R::ifelse1(par_analysis, object@experiment@BPPARAM, NULL)),
                  param_list)

    finalize_data_file <- .analysis_parser(".finalize_data", arg_list, files)
    object@tmp_data_files$finalize_data_file <- finalize_data_file
    tmp <- readRDS(finalize_data_file)
    object@features <- tmp$features
    object@spectra <- tmp$spectra@spectra
    object@peak_groups <- tmp$peak_groups
    write.csv(tmp$features, file.path(object@experiment@res_dir, "features.csv"))
    SpectraTools::ExportSpectra(SpectraTools::SpectraData(object@features, object@spectra),
                                file.path(object@experiment@res_dir, "spectra.msp"),
                                "msp")
    setwd(wd0)
    return(object)
  }
)

#' @export
setMethod(
  "FillPeaks",
  signature = c("TimsData", "FillPeakParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    message("Filling peaks...")
    files <- .tmp_files(object@files, object@experiment@tmp_dir, "fill_peaks")
    names(files) <- object@files

    param_list <- as.list(param)
    rev_models <- readRDS(object@tmp_data_files$align_file)$rev_models
    rt_correct_files <- object@tmp_data_files$correct_rt_files
    smp_names <- rownames(object@sample_groups)
    par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)
    fill_peak_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(idxs, function(idx) {
        is_fill <- is.na(object@features[, smp_names[idx]])
        info <- object@features[is_fill, , drop = FALSE]
        if (!is.null(rev_models[[idx]])) {
          info$rt <- predict(rev_models[[idx]], info$rt_align)
        } else if (file.exists(rt_correct_file <- rt_correct_files[object@files[idx]])) {
          land_marks <- readRDS(rt_correct_file)$land_marks
          dt_rt <- data.frame('ref' = c(object@experiment@rt_range[1],
                                        object@peaks[land_marks[, 'ref'], 'rt'],
                                        object@experiment@rt_range[2]),
                              'query' = c(object@experiment@rt_range[1],
                                          object@peaks[land_marks[, 'query'], 'rt'],
                                          object@experiment@rt_range[2])
          )
          rev_model <- loess(query ~ ref, data = dt_rt[order(dt_rt$query), , drop = FALSE], span = 0.1, degree = 1L)
          info$rt <- predict(rev_model, info$rt)
        }

        c(list("tims_data_file" = unname(object@tmp_data_files$tims_data_files[idx]),
               "info" = info,
               "res_define_at" = object@experiment@res_define_at),
          param_list)
      })
      .parallel_parser(".fill_peaks", arg_list, files[idxs],
                       object@experiment@BPPARAM, TRUE)
    }, simplify = FALSE))
    names(fill_peak_files) <- object@files
    object@tmp_data_files$fill_peak_files <- fill_peak_files

    for (idx in seq_along(object@files)) {
      is_fill <- is.na(object@features[, smp_names[idx]])
      object@features[is_fill, smp_names[idx]] <- readRDS(fill_peak_files[idx])
    }
    write.csv(object@features, file.path(object@experiment@res_dir, "features_filled.csv"))
    setwd(wd0)
    return(object)
  })

#' @export
setMethod(
  "IdentifyPeaks",
  signature = c("TimsData", "SearchParam", "MatchParam", "CombineParam"),
  function(object,
           search_param, match_param, combine_param,
           lib_file = NULL,
           rt_exp_file=NULL, rt_ref_file=NULL, 
           level3_lib_file = NULL, level3_lib_info = NULL,
           demo_mode = FALSE
           ) {
    # browser()
    wd0 <- getwd()
    object <- tims_data
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    message("Identify peaks...")

    files <- .tmp_files(object@files, object@experiment@tmp_dir, "identify_peaks")
    names(files) <- object@files

    if (!demo_mode && is.null(lib_file)) {
      stop("Either lib_file must be specified!")
    }

    pkg <- getPackageName()
    if (pkg == ".GlobalEnv") {
      pkg <- "Met4DX"
    }

    if (demo_mode) {

      lib_data <- SpectraTools:::ReadSpectraLib(system.file(package = pkg, "library",
                                                           paste0('tims', object@experiment@ion_mode, '.lib')))
    } else {
      lib_data <- SpectraTools::ParseSpectra(SpectraTools::ParseSpectraParam('msp',
                                                                             labelKeep = NULL,
                                                                             labelName = NULL,
                                                                             autoRename = TRUE,
                                                                             resDefineAt = 200,
                                                                             thrIntensityAbs = 0,
                                                                             ppmPrecursorFilter = 20),
                                             lib_file)

    }
    # browser()

    lib_data <- SpectraTools::setRT(lib_data, object@experiment@lc_column)

    exp_data <- SpectraTools::SpectraData(object@features, object@spectra)
    # browser()
    rt_ref <- rt_exp <- NULL
    if (search_param@scoreRT) {
      if (is.null(rt_ref_file)) {
        rt_ref_file <- system.file('rt_calibration',
                                   paste0(object@experiment@lc_column, '_', object@experiment@ion_mode, '.csv'),
                                   package = pkg)
      }
      rt_ref <- read.csv(rt_ref_file, stringsAsFactors = FALSE)
      rt_exp <- read.csv(rt_exp_file, stringsAsFactors = FALSE)
    }
    if (search_param@scoreCCS) {
      adduct_table <- read.csv(system.file('adducts',
                                           paste0("adducts_", object@experiment@ion_mode, ".csv"),
                                           package = pkg), stringsAsFactors = FALSE)
    } else {
      adduct_table <- NULL
    }
    lib_data@ccsInfo$`[M-H]-` <- as.numeric(lib_data@ccsInfo$`[M-H]-`)
    lib_data@ccsInfo$`[M+Na-2H]-` <- as.numeric(lib_data@ccsInfo$`[M+Na-2H]-`)
    lib_data@ccsInfo$`[M+HCOO]-` <- as.numeric(lib_data@ccsInfo$`[M+HCOO]-`)
    # browser()
    
    spec_searched <- SpectraTools::SearchSpectra(exp_data, lib_data, search_param,
                                                rtcalExp = rt_exp, rtcalRef = rt_ref,
                                                adductTable = adduct_table)
    # score_match <- BiocParallel::bplapply(spec_searched, function(specData) {
    score_match <- lapply(spec_searched, function(specData) { 
      # browser()
      # cat(names(specData))
      # specData <- spec_searched[["#989"]]
      dataExp <- specData$dataExp
      dataRef <- specData$dataRef
      dataRef <- SpectraTools::SpectraData(info = dataRef@info,
                                           spectra = dataRef@spectra)
      SpectraTools::MatchSpectra(dataExp, dataRef, match_param)
    })
    # browser()
    score_match <- score_match[!sapply(score_match, is.null)]
    # browser()
    cat("Plotting MSMS match figures ...\n")
    dirPlot <- file.path(object@experiment@tmp_dir, "MSMSMatchPlot")
    PlotMatchResult(score_match, expinfo = exp_data@info, dirPlot, addname = FALSE)
    scTable <- do.call(rbind, lapply(score_match, SpectraTools::GenOutputScore,
                                     match_param@cutoff, type = "metabolites"))
    
    pkTable <- MergeResTable(exp_data@info, scTable)
    pkTable <- add_level_class(pkTable, lib_data)
    write.csv(pkTable, file.path(object@experiment@res_dir, "result1_MSMSmatch.csv"),
              row.names = FALSE)
    
    cat("Finalizing scores ...\n")
    # score_match <- BiocParallel::bplapply(score_match, function(sc) {
    score_match <- lapply(score_match, function(sc) {
      SpectraTools::CombineScore(sc, combine_param)
    })
    score_match <- score_match[!sapply(score_match, is.null)]
    # browser()
    scTable <- do.call(rbind, lapply(score_match, SpectraTools::GenOutputScore,
                                     type = "metabolites"))
    pkTable <- MergeResTable(exp_data@info, scTable)
    pkTable <- add_level_class(pkTable, lib_data)
    write.csv(pkTable, file.path(object@experiment@res_dir, "result2_ScoreCombine.csv"),
              row.names = FALSE)
    
    ### the final table : remain the highest level for each feature #####
    cat("Refine confidence level ...\n")
    scTable <- remain_the_highest_level_res(scTable, 
                                            lib_data)
    # browser()
    pkTable <- MergeResTable(exp_data@info, scTable)
    pkTable <- add_level_class(pkTable, lib_data)
    write.csv(pkTable, file.path(object@experiment@res_dir, "result3_ScoreCombine_refined_level.csv"),
              row.names = FALSE)
    
    cat("Generate ms1 match result for MS-FIDNER ...\n")
    # browser()
    idx_not_id <- which(!row.names(exp_data@info) %in% row.names(scTable))
    
    new_exp_info <- exp_data@info[idx_not_id, ]
    new_exp_spec <- exp_data@spectra[idx_not_id]
    names(new_exp_spec) <- row.names(new_exp_info)
    new_expe_data <- SpectraTools::SpectraData(new_exp_info, 
                                               new_exp_spec)
    new_adduct_table <- adduct_table[1, ]
    
    level3_db <- gen_lib(lib_file = level3_lib_file, 
                         info_file = level3_lib_info, 
                         col_lib = c('Ion_mode', 'ExactMass', 'level'))
    
    level3_db@ccsInfo$`[M-H]-` <- as.numeric(level3_db@ccsInfo$`[M-H]-`)
    level3_db <- SpectraTools::setRT(level3_db, object@experiment@lc_column)
    
    new_spec_searched <- SpectraTools::SearchSpectra(new_expe_data, level3_db, search_param,
                                                     rtcalExp = rt_exp, rtcalRef = rt_ref,
                                                     adductTable = new_adduct_table)
    saveRDS(new_spec_searched, 
            file.path(object@experiment@res_dir, "spec_searched_for_msfinder"), 
            version = 2)
    
    setwd(wd0)
    return(object)
  })



PlotMatchResult <- function(scoreMatch, expinfo, dirPlot,
                            addname = FALSE, plotPNG = FALSE) {
  if (!dir.exists(dirPlot)) {
    dir.create(dirPlot)
  }
  
  # BiocParallel::bplapply(names(scoreMatch), function(nm) {
  lapply(names(scoreMatch), function(nm) {
    # cat(nm, '\t')
    require(ggplot2)
    pkname <- expinfo[nm, "name"]
    matchScore <- scoreMatch[[nm]]
    if (!is.null(matchScore)) {
      filePlot = file.path(dirPlot, pkname)
      SpectraTools::PlotMirror(matchScore, pkname, addname = addname,
                               plotPNG = plotPNG,
                               plotPDF = TRUE,
                               filePlot = filePlot,
                               direction = 'both')
    }
  })
  invisible()
}

MergeResTable <- function(info, res) {
  if (!is.data.frame(res)) {
    res <- SpectraTools:::.Col2Numeric(res)
  }
  tmp <- data.frame(matrix(ncol = ncol(res), nrow = nrow(info)))
  colnames(tmp) <- colnames(res)
  rownames(tmp) <- rownames(info)
  tmp[, sapply(res, is.character)] <- ""
  tmp[, sapply(res, is.numeric)] <- 0
  tmp[rownames(res), ] <- res
  info <- cbind(info, tmp)
  return(info)
}


add_level_class <- function(pkTable, 
                            lib_data){
  # parse the labid and match with the lib_data_info
  extra_info <- lapply(pkTable$labids, function(ids){
    ids <- strsplit(ids, ';')[[1]]
    idx <- match(ids, lib_data@info$labid)
    name <- paste0(lib_data@info$name[idx], collapse = ';')
    database <- paste0(lib_data@info$raw_id[idx], collapse = ';')
    level <- paste0(lib_data@info$level[idx], collapse = ';')
    formula <- paste0(lib_data@info$formula[idx], collapse = ';')
    exact_mass <- paste0(lib_data@info$mz[idx], collapse = ';')
    smiles <- paste0(lib_data@info$smiles[idx], collapse = ';')
    inchi <- paste0(lib_data@info$inchi[idx], collapse = ';')
    inchikey <- paste0(lib_data@info$inchikey[idx], collapse = ';')
    kingdom <- paste0(lib_data@info$kingdom[idx], collapse = ';')
    superclass <- paste0(lib_data@info$superclass[idx], collapse = ';')
    class <- paste0(lib_data@info$class[idx], collapse = ';')
    subclass <- paste0(lib_data@info$subclass[idx], collapse = ';')
    
    return(data.frame(name = name, 
                      database = database,
                      formula = formula, 
                      exact_mass = exact_mass, 
                      smiles = smiles, 
                      inchi = inchi, 
                      inchikey = inchikey,
                      confidence_level = level, 
                      kingdom = kingdom, 
                      superclass = superclass, 
                      class = class, 
                      subclass = subclass))
  })
  extra_info <- do.call(rbind, extra_info)
  return(cbind(pkTable, extra_info))
}

remain_the_highest_level_res <- function(scTable, 
                                         lib_data){
  res_refined <- lapply(seq(nrow(scTable)), function(t){
    temp_row <- scTable[t, ]
    id <- strsplit(temp_row$labids, ';')[[1]]
    n_id <- length(id)
    
    if(n_id <= 1){
      return(temp_row)
    }else{
      # browser()
      level <- as.numeric(lib_data@info$level[match(id,
                                                    lib_data@info$labid)])
      highest_level <- min(level)
      final_idx <- which(level == highest_level)
      
      final_res <- apply(temp_row, 2, function(v){
        remain_v <- strsplit(v, ';')[[1]][final_idx]
        return(paste0(remain_v, collapse = ';'))
      })
      temp_res <- as.data.frame(t(final_res))
      row.names(temp_res) <- row.names(temp_row)
      return(temp_res)
    }
  })
  return(do.call(rbind, res_refined))
}

