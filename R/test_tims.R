# setGeneric("ExtractTimDataIonList",
#            function(object, param,
#                     ion_of_interest_path, ...)
#              standardGeneric("ExtractTimDataIonList"))
# 
# setGeneric("DereplicatePeaksIonList",
#            function(object, param)
#              standardGeneric("DereplicatePeaksIonList"))


setGeneric("MatchBetweenRuns_IOI",
           function(object, param, ion_of_interest_path, ...)
             standardGeneric("MatchBetweenRuns_IOI"))

setGeneric("FinalizeFeatures_IOI",
           function(object, param, ...)
             standardGeneric("FinalizeFeatures_IOI"))


setClass("FinalizeFeatureIOIParam",
         slots = c(
           min_fraction = "numeric",
           msms_assign_method = "character",
           min_num_samples = "numeric",
           valid_eic_peak = "logical",
           valid_eim_peak = "logical",
           snthreshold_eic = "nullOrNumeric",
           snthreshold_eim = "nullOrNumeric",
           quant_method = "character",
           col_max = "character",
           col_quant = "character",
           mz_diff_spectral_assgined = 'numeric', 
           rt_diff_spectral_assgined = 'numeric', 
           mobility_diff_spectral_assgined = 'numeric',
           rerun = "logical"
         ),
         contains = "Met4DXParam"
)

setGeneric("FillPeaks_IOI",
           function(object, param)
             standardGeneric("FillPeaks_IOI"))


## FinalizeFeatureParam
#‘ parameters for finalize data
#' @export
FinalizeFeatureIOIParam <- function(
    min_fraction = 0.5,
    msms_assign_method = "highest",
    min_num_samples = 1,
    valid_eic_peak = FALSE,
    valid_eim_peak = TRUE,
    snthreshold_eic = NULL,
    snthreshold_eim = NULL,
    quant_method = c("ccs_median", "rt_median", "max"),
    col_max = c("area", "target_intensity", "height", "height_fit"),
    col_quant = c("area", "height", "height_fit"),
    mz_diff_spectral_assgined = 20, 
    rt_diff_spectral_assgined = 10, 
    mobility_diff_spectral_assgined = 0.02,
    rerun = FALSE
) {
  quant_method <- match.arg(quant_method)
  col_max <- match.arg(col_max)
  col_quant <- match.arg(col_quant)
  
  new("FinalizeFeatureIOIParam",
      min_fraction = min_fraction,
      msms_assign_method = msms_assign_method,
      min_num_samples = min_num_samples,
      valid_eic_peak = valid_eic_peak,
      valid_eim_peak = valid_eim_peak,
      snthreshold_eic = snthreshold_eic,
      snthreshold_eim = snthreshold_eim,
      quant_method = quant_method,
      col_max = col_max,
      col_quant = col_quant,
      mz_diff_spectral_assgined = mz_diff_spectral_assgined, 
      rt_diff_spectral_assgined = rt_diff_spectral_assgined, 
      mobility_diff_spectral_assgined = mobility_diff_spectral_assgined,
      rerun = rerun
  )
}





.ccs2mobility <- function(ccs, mz, t = 305, z = 1){
  factor <- 18509.863216340458
  # em_n2 <- 28.0134
  em_n2 <- 28.00615
  k0 <- ccs/(
    factor *
    z *
    sqrt((mz + em_n2) / (t * em_n2 * mz))
    )
  return(k0)
}



.extract_im_data_roi <- function(
    info,
    tims_data_file,
    precursor_bin_file = NULL,
    order_column = "mz",
    mz_tol = 20,
    frame_range = 30,
    frame_integration_range = 5,
    mobility_range = 0.1,
    mobility_intgration_range = 0.015,
    min_points = 4,
    min_intensity = 0,
    n_skip = 0,
    interpolate_method = NULL,
    peak_span_eim = 21,
    peak_span_eic = 7,
    snthreshold = 3,
    smooth_window_eim = 16,
    smooth_window_eic = 8,
    keep_profile = FALSE,
    res_define_at = 200,
    use_cmpd_id = TRUE,
    smooth_method = 'gaussian',
    skip_invalid_eic_peaks = TRUE,
    skip_invalid_eim_peaks = TRUE,
    filter_outlier_peaks = FALSE,
    allowed_mobility_shift = 0.02,
    allowed_rt_shift = 30,
    ...
) {
  # browser()
  # ioi <- read.csv('./ioi_tims.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  # ioi$mobility <- .ccs2mobility(ccs = ioi$ccs, mz = ioi$mz)
  # ioi$ccs <- NULL
  # info <- ioi
  remove_profile <- !keep_profile
  
  query_data <- readRDS(tims_data_file)
  num_frames <- nrow(query_data$ms1_frame_info)
  
  if (!missing(precursor_bin_file) && !is.null(precursor_bin_file)) {
    precursor_bins <- readRDS(precursor_bin_file)
    info <- info[sapply(precursor_bins, `[`, 1),]
  }
  
  if (mz_tol > 1) {
    info$mz_tol <- sapply(info$mz, function(mz) .ppm2dalton(mz, mz_tol, res_define_at))
  } else {
    info$mz_tol <- mz_tol
  }
  
  mobility_col <- ifelse("k0" %in% colnames(info), "k0", "mobility")
  intensity_col <- ifelse("intensity" %in% colnames(info), "intensity", "target_intensity")
  info <- info[order(info[[order_column]], decreasing = TRUE),]
  switch(smooth_method,
         "gaussian" = {
           smooth_param_eim <- GaussianSmoothParam(window = smooth_window_eim)
           smooth_param_eic <- GaussianSmoothParam(window = smooth_window_eic)
         },
         "loess" = {
           smooth_param_eim <- LOESSSmoothParam(window = smooth_window_eim)
           smooth_param_eic <- LOESSSmoothParam(window = smooth_window_eic)
         },
         stop("Unknown smooth method")
  )
  
  if (use_cmpd_id && ("cmpd" %in% colnames(info))) {
    info$target_frame <- sapply(info[, "cmpd"], function(cmpd) {
      query_data$precursor_info$Parent[query_data$precursor_info$Id == cmpd]
    })
  } else {
    info$target_frame <- sapply(info[, "rt"], function(rt) {
      query_data$ms1_frame_info$Id[which.min(abs(query_data$ms1_frame_info$Time - rt))]
    })
  }
  
  peak_quality_tmp <- as.list(c(rep(FALSE, 2), rep(NA, 4)))
  names(peak_quality_tmp) <- c("eim_peak", "eic_peak", "eim_sn", "eic_sn", "eim_baseline", "eic_baseline")
  
  
  
  # mz <- info$mz[1]
  # mobility <- info$mobility[1]
  # mz_tol <- info$mz_tol[1]
  # target_frame <- info$target_frame[1]
  
  res <- apply(info[, c("mz", intensity_col, mobility_col, "mz_tol", "target_frame")], 1, function(dr) {
    # browser()
    mz <- dr["mz"]
    mobility <- dr[mobility_col]
    mz_tol <- dr["mz_tol"]
    target_frame <- dr["target_frame"]
    
    peak_quality <- peak_quality_tmp
    
    target_idx <- which(query_data$ms1_frame_info$Id == target_frame)
    frame_start <- max(0, target_idx - frame_range)
    extract_frames <- query_data$ms1_frame_info$Id[frame_start:min(num_frames, target_idx + frame_range)]
    num_extract_frames <- length(extract_frames)
    
    # get eim from the target frame
    # eim <- .get_eims2(query_data$all_frames[as.character(target_frame)], query_data$all_mobility,
    #                   mz, mz_tol, mobility, mobility_range)
    # # eim <- .get_eims(D, target_frame, mz, mz_tol, mobility, mobility_range, all_columns)
    # 
    # if (sum(eim[, 1]) == 0) {
    #   return(NULL)
    # }
    # 
    # if (is.null(eim)) {
    #   return(NULL)
    # }
    # colnames(eim)[1] <- "intensity"
    # eim <- .interpolate_data(eim, c(-1, 1) * mobility_range + mobility, interpolate_method)
    # 
    # eim$intensity_smooth <- splus2R::ifelse1(smooth_method == "gaussian",
    #                                          .smooth_gaussian(data = eim$intensity, window = smooth_param_eim@window),
    #                                          .smooth_loess(data = eim$intensity, degree = smooth_param_eim@degree, window = smooth_param_eim@window)
    # )
    # 
    # # if (.get_smooth_sd(eim$intensity, eim$intensity_smooth) > 0.35) {
    # #   return(NULL)
    # # }
    # 
    # # deterim the eim apex for optimized eic-eim profile extraction
    # # ref_index <- which.min(abs(eim$k0 - mobility))
    # apex_eim <- eim$k0[.find_apex(eim$intensity, eim$intensity_smooth,
    #                               span = peak_span_eim,
    #                               min_points = min_points,
    #                               min_intensity = min_intensity,
    #                               n_skip = n_skip,
    #                               ref_index = ref_index,
    #                               find_roi = FALSE)]
    # if (length(apex_eim) == 0) {
    #   return(NULL)
    # }
    
    # .find_roi(eim$intensity,
    #           min_points = 2,
    #           min_intensity = 0,
    #           n_skip = 3,
    #           ref = 2)
    
    
    # extract eic-eim profile around the target frame and eim apex
    eims <- .get_eims2(query_data$all_frames[as.character(extract_frames)], query_data$all_mobility,
                       mz, mz_tol, mobility, mobility_range)
    check_roi <- apply(eims[,-ncol(eims)], 2 , function(fr){
      res_roi <- .find_roi(fr,
                           min_points = 8,
                           min_intensity = 0,
                           n_skip = 3,
                           ref = which.min(abs(eims$k0 - mobility)))
    })
    temp_ref <- which.min(abs(eims$k0 - mobility))
    check_roi <- check_roi[which(!sapply(check_roi, is.null))]
    # check_roi_sum <- sapply(names(check_roi), function(fr){
    #   # browser()
    #   scan_min <- check_roi[[fr]][1, 'scan_min']
    #   scan_max <- check_roi[[fr]][1, 'scan_max']
    #   sum(eims[temp_ref - ceiling(15/2):temp_ref + ceiling(15/2), as.character(fr)])
    # })
    # target_frame <- names(which.max(check_roi_sum))
    # target_frame <- query_data$ms1_frame_info$Id[which.min(abs(query_data$ms1_frame_info[names(check_roi),]$Time - query_data$ms1_frame_info[as.character(dr["target_frame"]),]$Time))]
    target_frame <- query_data$ms1_frame_info[names(check_roi),]$Id[which.min(abs(query_data$ms1_frame_info[names(check_roi),]$Time - query_data$ms1_frame_info[as.character(dr["target_frame"]),]$Time))]
    # target_frame <- query_data$ms1_frame_info$Id[which.min(abs(query_data$ms1_frame_info$Time - rt))]
    
    
    target_idx <- which(query_data$ms1_frame_info$Id == target_frame)
    frame_start <- max(0, target_idx - frame_range)
    extract_frames <- query_data$ms1_frame_info$Id[frame_start:min(num_frames, target_idx + frame_range)]
    num_extract_frames <- length(extract_frames)
    
    # get eim from the target frame
    eim <- .get_eims2(query_data$all_frames[as.character(target_frame)], query_data$all_mobility,
                      mz, mz_tol, mobility, mobility_range)
    # eim <- .get_eims(D, target_frame, mz, mz_tol, mobility, mobility_range, all_columns)
    
    if (sum(eim[, 1]) == 0) {
      return(NULL)
    }
    
    if (is.null(eim)) {
      return(NULL)
    }
    colnames(eim)[1] <- "intensity"
    eim <- .interpolate_data(eim, c(-1, 1) * mobility_range + mobility, interpolate_method)
    
    eim$intensity_smooth <- splus2R::ifelse1(smooth_method == "gaussian",
                                             .smooth_gaussian(data = eim$intensity, window = smooth_param_eim@window),
                                             .smooth_loess(data = eim$intensity, degree = smooth_param_eim@degree, window = smooth_param_eim@window)
    )
    
    # if (.get_smooth_sd(eim$intensity, eim$intensity_smooth) > 0.35) {
    #   return(NULL)
    # }
    
    # deterim the eim apex for optimized eic-eim profile extraction
    ref_index <- which.min(abs(eim$k0 - mobility))
    apex_eim <- eim$k0[.find_apex(eim$intensity, eim$intensity_smooth,
                                  span = peak_span_eim,
                                  min_points = min_points,
                                  min_intensity = min_intensity,
                                  n_skip = n_skip,
                                  ref_index = ref_index,
                                  find_roi = FALSE)]
    if (length(apex_eim) == 0) {
      return(NULL)
    }
    # re-extract the eims2
    
    
    
    eims <- .get_eims2(query_data$all_frames[as.character(extract_frames)], query_data$all_mobility,
                       mz, mz_tol, apex_eim, mobility_range)
    
    # eims <- .get_eims(D, extract_frames, mz, mz_tol, apex_eim, mobility_range, all_columns)
    # colnames(eims)[1:num_extract_frames] <- extract_frames
    eic_rt <- query_data$ms1_frame_info[as.character(extract_frames), "Time"]
    im_data <- new("IMData",
                   eic_rt = eic_rt,
                   eic_mz = attributes(eims)$mz,
                   eim_mobility = eims[, num_extract_frames + 1],
                   profile_data = eims[, 1:num_extract_frames])
    im_data <- SetData(im_data, "eic", apex_eim + c(-1, 1) * mobility_intgration_range)
    im_data <- SmoothData(im_data, "eic", smooth_param_eic)
    apex_eic <- .find_apex(im_data@eic$raw, im_data@eic$smooth,
                           span = peak_span_eic,
                           min_points = min_points,
                           min_intensity = min_intensity,
                           n_skip = n_skip,
                           ref_index = frame_range + 1,
                           force_max = TRUE)
    if (is.null(apex_eic)) {
      return(NULL)
    }
    
    if (.get_smooth_sd(im_data@eic$raw, im_data@eic$smooth) > 0.35) {
      if (skip_invalid_eic_peaks) {
        return(NULL)
      }
    } else {
      peak_quality$eic_peak <- TRUE
    }
    
    eic_quality <- .get_peak_quality(im_data@eic$raw,
                                     apex_eic,
                                     min_intensity = 0,
                                     min_points = min_points,
                                     n_skip = n_skip,
                                     snthreshold = snthreshold,
                                     skip_invalid_peaks = skip_invalid_eic_peaks)
    
    if (is.null(eic_quality) && skip_invalid_eic_peaks) {
      return(NULL)
    }
    
    peak_quality$eic_sn <- eic_quality$snr
    peak_quality$eic_baseline <- eic_quality$baseline
    
    im_data <- SetData(im_data, "eim", seq(max(1, apex_eic - 5), min(apex_eic + 5, num_extract_frames)))
    im_data <- InterpolateEIM(im_data, c(-1, 1) * mobility_range + mobility, interpolate_method)
    im_data <- SmoothData(im_data, "eim", smooth_param_eim)
    
    if (.get_smooth_sd(im_data@eim$raw, im_data@eim$smooth) > 0.35) {
      if (skip_invalid_eim_peaks) {
        return(NULL)
      }
    } else {
      peak_quality$eim_peak <- TRUE
    }
    
    ref_index <- which.min(abs(im_data@eim_mobility_interpolate - mobility))
    
    apex_eim_index <- .find_apex(im_data@eim$interpolate, im_data@eim$smooth,
                                 span = peak_span_eim,
                                 min_points = min_points,
                                 min_intensity = min_intensity,
                                 n_skip = n_skip,
                                 ref_index = ref_index,
                                 find_roi = FALSE)
    
    eim_quality <- .get_peak_quality(im_data@eim$interpolate,
                                     apex_eim_index,
                                     min_intensity = 0,
                                     min_points = min_points, 
                                     n_skip = n_skip,
                                     snthreshold = snthreshold,
                                     skip_invalid_peaks = skip_invalid_eim_peaks)
    
    if (is.null(eim_quality) && skip_invalid_eim_peaks) {
      return(NULL)
    }
    
    peak_quality$eim_sn <- eim_quality$snr
    peak_quality$eim_baseline <- eim_quality$baseline
    
    apex_eim <- im_data@eim_mobility_interpolate[apex_eim_index]
    if (length(apex_eim) == 0) {
      return(NULL)
    }
    if (filter_outlier_peaks) {
      if (abs(mobility - apex_eim) > allowed_mobility_shift ||
          abs(query_data$ms1_frame_info[as.character(target_frame), "Time"] - eic_rt[apex_eic]) > allowed_rt_shift) {
        return(NULL)
      }
    }
    im_data <- SetData(im_data, "eic", apex_eim + c(-1, 1) * mobility_intgration_range)
    im_data <- SetInfo(im_data, apex_eic, apex_eim, unname(dr[intensity_col]), frame_integration_range)
    im_data@peak_quality <- peak_quality
    if (remove_profile) {
      im_data <- DropProfile(im_data)
    }
    return(im_data)
  })
  # browser()
  names(res) <- rownames(info)
  return(res)
}




#' @export
setMethod(
  "MatchBetweenRuns_IOI",
  signature = c("TimsData", "MatchBetweenRunParam", 'ANY'),
  function(object, param, ion_of_interest_path) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    message("Matching features between runs...")

    sub_dir <- ifelse(is.null(param@interpolate_method), "match_between_runs",
                      file.path("match_between_runs", param@interpolate_method))
    files <- .tmp_files(object@files, object@experiment@tmp_dir, sub_dir)
    names(files) <- object@files

    ### read ioi and converted ccs into mobility ####
    ioi <- read.csv(ion_of_interest_path, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
    ioi$mobility <- .ccs2mobility(ccs = ioi$ccs, mz = ioi$mz)
    ioi$ccs <- NULL
    object@features <- ioi
    param_list <- as.list(param)
    # rev_models <- readRDS(object@tmp_data_files$align_file)$rev_models
    # rt_correct_files <- object@tmp_data_files$correct_rt_files

    par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)
    # pks <- object@peaks[unname(do.call(c, object@peak_groups)), "smp_idx", drop = FALSE]
    # pks$ft_idx <- do.call(c, mapply(rep, names(object@peak_groups), sapply(object@peak_groups, length), simplify=FALSE))
    match_between_run_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(idxs, function(idx) {
        # is_fill <- !rownames(object@features) %in% pks[pks$smp_idx == idx, "ft_idx"]
        # info <- object@features[is_fill, , drop = FALSE]
        # if (!is.null(rev_models[[idx]])) {
        #   info$rt <- predict(rev_models[[idx]], info$rt_align)
        # } else if (file.exists(rt_correct_file <- rt_correct_files[object@files[idx]])) {
        #   land_marks <- readRDS(rt_correct_file)$land_marks
        #   dt_rt <- data.frame('ref' = c(object@experiment@rt_range[1],
        #                                 object@peaks[land_marks[, 'ref'], 'rt'],
        #                                 object@experiment@rt_range[2]),
        #                       'query' = c(object@experiment@rt_range[1],
        #                                   object@peaks[land_marks[, 'query'], 'rt'],
        #                                   object@experiment@rt_range[2])
        #   )
        #   rev_model <- loess(query ~ ref, data = dt_rt[order(dt_rt$query), , drop = FALSE], span = 0.1, degree = 1L)
        #   info$rt <- predict(rev_model, info$rt)
        # }

        c(list("tims_data_file" = unname(object@tmp_data_files$tims_data_files[idx]),
               "info" = ioi,
               "precursor_bin_file" = NULL,
               "res_define_at" = object@experiment@res_define_at),
          param_list)
      })
      .parallel_parser(".match_between_runs_ioi", arg_list, files[idxs],
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
        # if (file.exists(rt_correct_file <- rt_correct_files[object@files[idx]])) {
        #   land_marks <- readRDS(rt_correct_file)$land_marks
        #   dt_rt <- data.frame('ref' = c(object@experiment@rt_range[1],
        #                                 object@peaks[land_marks[, 'ref'], 'rt'],
        #                                 object@experiment@rt_range[2]),
        #                       'query' = c(object@experiment@rt_range[1],
        #                                   object@peaks[land_marks[, 'query'], 'rt'],
        #                                   object@experiment@rt_range[2])
        #   )
        #   rt_model <- loess(ref ~ query, data = dt_rt[order(dt_rt$ref), , drop = FALSE], span = 0.1, degree = 1L)
        # }
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

.match_between_runs_ioi <- function(
    info,
    tims_data_file,
    precursor_bin_file = NULL,
    mz_tol = 20,
    frame_range = 30,
    frame_integration_range = 5,
    mobility_range = 0.1,
    mobility_intgration_range = 0.015,
    min_points = 4,
    min_intensity = 0,
    n_skip = 0,
    interpolate_method = NULL,
    keep_profile = FALSE,
    filter_outlier_peaks = TRUE,
    allowed_mobility_shift = 0.015,
    allowed_rt_shift = 10,
    res_define_at = 200,
    data_file = NULL,
    ...
) {
  im_data <- .extract_im_data_roi(
    info = info,
    tims_data_file = tims_data_file,
    precursor_bin_file = NULL,
    order_column = "mz",
    mz_tol = mz_tol,
    frame_range = frame_range,
    frame_integration_range = frame_integration_range,
    mobility_range = mobility_range,
    mobility_intgration_range = mobility_intgration_range,
    min_points = min_points,
    min_intensity = min_intensity,
    n_skip = n_skip,
    interpolate_method = interpolate_method,
    keep_profile = keep_profile,
    res_define_at = res_define_at,
    use_cmpd_id = FALSE,
    skip_invalid_eic_peaks= FALSE,
    skip_invalid_eim_peaks= FALSE,
    filter_outlier_peaks = filter_outlier_peaks,
    allowed_mobility_shift = allowed_mobility_shift,
    allowed_rt_shift = allowed_rt_shift,
  )
  if (!is.null(data_file)) {
    saveRDS(im_data, file = data_file, version = 2)
  } else {
    return(im_data)
  }
}

#' @export
setMethod(
  "FinalizeFeatures_IOI",
  signature = c("TimsData", "FinalizeFeatureIOIParam"),
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
                       "spectra_files" = object@tmp_data_files$spectra_files,
                       "peak_groups" = object@peak_groups,
                       "res_define_at" = object@experiment@res_define_at,
                       "bpparam" = splus2R::ifelse1(par_analysis, object@experiment@BPPARAM, NULL)),
                  param_list)

    finalize_data_file <- .analysis_parser(".finalize_data_ioi", arg_list, files)
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

# peaks = object@peaks
# features = object@features
# sample_groups = object@sample_groups
# spectra_files = object@tmp_data_files$spectra_files
# peak_groups = object@peak_groups
# valid_eic_peak = FALSE
# valid_eim_peak = FALSE
# snthreshold_eic <- NULL
# snthreshold_eim <- NULL
# min_fraction = 0.5
# min_num_samples = 1
# col_quant <- 'area'
# col_quant = "area"
# mz_diff_spectral_assgined = 20 
# rt_diff_spectral_assgined = 10
# mobility_diff_spectral_assgined = 0.02
# bpparam <- NULL

.finalize_data_ioi <- function(
    peaks,
    filled_peaks,
    features,
    sample_groups,
    spectra_files,
    peak_groups,
    min_fraction = 0.5,
    min_num_samples = 1,
    valid_eic_peak = TRUE,
    valid_eim_peak = TRUE,
    snthreshold_eic = NULL,
    snthreshold_eim = NULL,
    quant_method = "max",
    col_max = "area",
    col_quant = "area",
    mz_diff_spectral_assgined, 
    rt_diff_spectral_assgined, 
    mobility_diff_spectral_assgined,
    res_define_at = 200,
    bpparam = NULL,
    ...
) {
  # browser()
  # features <- features[, -match(c("npeaks", "area"), colnames(features)), drop = FALSE]

  sample_classes <- .sample_groups_to_class(sample_groups)

  # spectra_all <- .get_peak_spectra_ioi(peaks, spectra_files, update_info = FALSE)

  filled_groups <- sapply(unique(filled_peaks$feature_idx), function(feature_idx) {
    rownames(filled_peaks[filled_peaks$feature_idx == feature_idx, "smp_idx", drop = FALSE])
  }, simplify = FALSE)
  all_groups <- sapply(rownames(features), function(nm) {
    c(peak_groups[[nm]], filled_groups[[nm]])
  }, simplify = FALSE)

  # cols_keep <- intersect(colnames(peaks), colnames(filled_peaks))
  cols_keep <- colnames(filled_peaks)
  all_peaks <- filled_peaks[, cols_keep]
  is_keep <- rep(TRUE, nrow(all_peaks))
  if (valid_eic_peak) {
    is_keep <- is_keep & all_peaks$eic_peak
  }
  if (valid_eim_peak) {
    is_keep <- is_keep & all_peaks$eim_peak
  }
  if (!is.null(snthreshold_eic)) {
    is_keep <- is_keep & all_peaks$eic_sn > snthreshold_eic
  }
  if (!is.null(snthreshold_eim)) {
    is_keep <- is_keep & all_peaks$eim_sn > snthreshold_eim
  }
  all_peaks <- all_peaks[is_keep, , drop = FALSE]
  all_peak_names <- rownames(all_peaks)
  cat(' extracting peaks for each feature ...')
  if (is.null(bpparam)) {
    group_sample_peaks <- lapply(all_groups, function(grp) {
      kkkk <- all_peaks[na.omit(fastmatch::fmatch(grp, all_peak_names)), , drop = FALSE]
      kkkk$rt_align <- kkkk$rt
      return(kkkk)
    })
  } else {
    group_sample_peaks <- BiocParallel::bplapply(all_groups, function(grp) {
      kkkk <- all_peaks[na.omit(fastmatch::fmatch(grp, all_peak_names)), , drop = FALSE]
      kkkk$rt_align <- kkkk$rt
      return(kkkk)
    }, BPPARAM = bpparam)
  }
  feature_sample_number <- sapply(group_sample_peaks, nrow)
  group_sample_num <- do.call(rbind, lapply(group_sample_peaks, function(pks) {
    table(sample_classes[unique(pks$smp_idx)])
  }))

  features$npeaks <- feature_sample_number
  features <- cbind(features, group_sample_num)

  group_num <- table(sample_classes)
  group_fraction <- do.call(rbind, apply(group_sample_num, 1, function(dr) dr / group_num, simplify = F))

  is_keep <- group_fraction >= min_fraction & apply(group_sample_num, 1, function(dr) any(dr >= min_num_samples))
  is_keep <- apply(is_keep, 1, any)
  features <- features[is_keep, , drop = FALSE]
  rec_col <- c("mz", "rt_align", "mobility")

  tmp <- vector("numeric", 9)
  names(tmp) <- c("mz", "mzmin", "mzmax", "rt_align", "rt_alignmin", "rt_alignmax", "mobility", "mobilitymin", "mobilitymax")
  update_info <- t(sapply(rownames(features), function(nms) {
    # cat(nms, '\t')
    dt <- group_sample_peaks[[nms]][, rec_col, drop = FALSE]
    rgs <- colRanges(dt)
    res <- tmp
    res[rec_col] <- colMedians(dt)
    for (col in rec_col) {
      res[paste0(col, c("min", "max"))] <- rgs[, col]
    }
    res
  }, simplify = TRUE))
  colnames(update_info) <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "mobility", "mobilitymin", "mobilitymax")
  features[, colnames(update_info)] <- update_info

  spectra_all <- .get_peak_spectra_ioi(features = features, spectra_files = spectra_files, min_num_fragments = NULL, update_info = FALSE, 
                                       mz_diff_spectral_assgined, 
                                       rt_diff_spectral_assgined, 
                                       mobility_diff_spectral_assgined)

  # spectra_assign_nms <- sapply(rownames(features), function(nm) {
  #   nms_spec <- intersect(peak_groups[[nm]], rownames(group_sample_peaks[[nm]]))
  #   idx_max <- which.max(spectra_all[nms_spec]@info[[col_max]])
  #   nms_spec[idx_max]
  # })
  # peak_groups <- peak_groups[is_keep]
  features$ccs <- .mobility2ccs(features$mobility, features$mz)
  features$name <- .get_peak_name(features[, c("mz", "rt", "ccs")])
  features <- features[, c("name", colnames(features)[1:(ncol(features) - 1)])]
  # browser()
  smp_idxs <- seq(nrow(sample_groups))
  # all_peaks <- rbind(peaks[, cols_keep], filled_peaks[, cols_keep])
  all_peaks <- filled_peaks[, cols_keep]
  cat(' quantifying peaks...')
  if (is.null(bpparam)) {

    quant_data <- do.call(rbind,lapply(rownames(features), function(nm) {
      valid_peaks <- group_sample_peaks[[nm]]
      extra_idx <- na.omit(fastmatch::fmatch(setdiff(all_groups[[nm]], rownames(valid_peaks)), all_peak_names))
      extra_peaks <- all_peaks[extra_idx, , drop = FALSE]
      sapply(smp_idxs, function(idx) {
        res <- NA
        pks <- valid_peaks[valid_peaks$smp_idx == idx, c("mobility", "rt_align", col_quant)]
        if (nrow(pks) == 0) {
          pks <- extra_peaks[extra_peaks$smp_idx == idx, c("mobility", "rt_align", col_quant)]
        }
        if ((pkrows <- nrow(pks)) > 1) {
          res <- switch(quant_method,
                        "max" = pks[which.max(pks[, col_quant]), col_quant],
                        "rt_median" = pks[which.min(pks$rt_align - features[nm, 'rt']), col_quant],
                        "ccs_median" = pks[which.min(pks$mobility - features[nm, 'mobility']), col_quant],
                        stop("Unsupported quant method: ", quant_method))
        } else if (pkrows == 1) {
          res <- pks[, col_quant]
        }
        res
      })
    }))
  } else {
    quant_data <- do.call(rbind, BiocParallel::bplapply(rownames(features), function(nm) {
      valid_peaks <- group_sample_peaks[[nm]]
      extra_idx <- na.omit(fastmatch::fmatch(setdiff(all_groups[[nm]], rownames(valid_peaks)), all_peak_names))
      extra_peaks <- all_peaks[extra_idx, , drop = FALSE]
      sapply(smp_idxs, function(idx) {
        res <- NA
        pks <- valid_peaks[valid_peaks$smp_idx == idx, c("mobility", "rt_align", col_quant)]
        if (nrow(pks) == 0) {
          pks <- extra_peaks[extra_peaks$smp_idx == idx, c("mobility", "rt_align", col_quant)]
        }
        if ((pkrows <- nrow(pks)) > 1) {
          res <- switch(quant_method,
                        "max" = pks[which.max(pks[, col_quant]), col_quant],
                        "rt_median" = pks[which.min(pks$rt_align - features[nm, 'rt']), col_quant],
                        "ccs_median" = pks[which.min(pks$mobility - features[nm, 'mobility']), col_quant],
                        stop("Unsupported quant method: ", quant_method))
        } else if (pkrows == 1) {
          res <- pks[, col_quant]
        }
        res
      })
    }, BPPARAM = bpparam))
  }
  colnames(quant_data) <- rownames(sample_groups)
  # spec <- spectra_all[spectra_assign_nms]
  spec <- spectra_all
  # spec <- SpectraTools::UpdateNames(spec, rownames(spectra_all@info$))
  
  return(list("features" = cbind(features, quant_data), "spectra" = spec, "peak_groups" = peak_groups))
}
#### a function was needed to assign ms2 spectra to features ######
.get_peak_spectra_ioi <- function(features, spectra_files, min_num_fragments = NULL, update_info = FALSE, 
                                  mz_diff_spectral_assgined = mz_diff_spectral_assgined, 
                                  rt_diff_spectral_assgined = rt_diff_spectral_assgined, 
                                  mobility_diff_spectral_assgined = mobility_diff_spectral_assgined) {
  # browser()
  spectra_data <- lapply(seq_along(spectra_files), function(idx) {
    readRDS(spectra_files[idx])
  })

  info <- do.call(rbind, lapply(spectra_data, function(x) {x@info}))
  rownames(info) <- .gen_indexes(info)

  spec <- do.call(c, lapply(seq_along(spectra_files), function(idx) {
    spectra_data[[idx]]@spectra
  }))

  ### assign the spec into the features according to the spectral intensity

  spec_index <- sapply(seq(nrow(features)), function(temp_row){
    # browser()
    ft <- features[temp_row, ]
    
    mz_diff <- abs(ft$mz - info$mz)/ft$mz * 10 ^ 6
    rt_diff <- abs(ft$rt - info$rt)
    mobility_diff <- abs(ft$mobility - info$k0)

    idx_mz <- which(mz_diff <= mz_diff_spectral_assgined)
    idx_rt <- which(rt_diff <= rt_diff_spectral_assgined)
    idx_mobility <- which(mobility_diff <= mobility_diff_spectral_assgined)

    final_idx <- intersect(idx_mz, intersect(idx_rt, idx_mobility))

    if(length(final_idx) > 1){
      final_idx <- final_idx[which.max(info$intensity[final_idx])]
    }else if (length(final_idx) == 1){
      final_idx <- final_idx
    }else{
      final_idx <- NA
    }
    
    return(final_idx)
  }, simplify = TRUE)


  # idxxxxx <- fastmatch::fmatch(rownames(info), names(spec))
  spec <- spec[spec_index]
  names(spec) <- rownames(features)
  
  check_null <- sapply(spec, is.null)
  spec <- spec[-which(check_null)]
  new_info <- features[-which(check_null), ]
  
  # if (!is.null(min_num_fragments)) {
  #   is_keep <- sapply(spec, nrow) >= min_num_fragments
  #   info <- info[is_keep, , drop = FALSE]
  #   spec <- spec[is_keep]
  # }
  
  SpectraTools::SpectraData(info = new_info,
                            spectra = spec)
}



##### final fill gap #####
#' @export
setMethod(
  "FillPeaks_IOI",
  signature = c("TimsData", "FillPeakParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    
    message("Filling peaks...")
    files <- .tmp_files(object@files, object@experiment@tmp_dir, "fill_peaks")
    names(files) <- object@files
    
    param_list <- as.list(param)
    # rev_models <- readRDS(object@tmp_data_files$align_file)$rev_models
    # rt_correct_files <- object@tmp_data_files$correct_rt_files
    smp_names <- rownames(object@sample_groups)
    par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)
    fill_peak_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(idxs, function(idx) {
        # browser()
        is_fill <- is.na(object@features[, smp_names[idx]])
        info <- object@features[is_fill, , drop = FALSE]
        # if (!is.null(rev_models[[idx]])) {
        #   info$rt <- predict(rev_models[[idx]], info$rt_align)
        # } else if (file.exists(rt_correct_file <- rt_correct_files[object@files[idx]])) {
        #   land_marks <- readRDS(rt_correct_file)$land_marks
        #   dt_rt <- data.frame('ref' = c(object@experiment@rt_range[1],
        #                                 object@peaks[land_marks[, 'ref'], 'rt'],
        #                                 object@experiment@rt_range[2]),
        #                       'query' = c(object@experiment@rt_range[1],
        #                                   object@peaks[land_marks[, 'query'], 'rt'],
        #                                   object@experiment@rt_range[2])
        #   )
        #   rev_model <- loess(query ~ ref, data = dt_rt[order(dt_rt$query), , drop = FALSE], span = 0.1, degree = 1L)
        #   info$rt <- predict(rev_model, info$rt)
        # }
        
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



