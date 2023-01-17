################################
# all_classes
# all_generics
# ImsData
setGeneric("ImsData",
           function(experiment, ...)
             standardGeneric("ImsData"))

# QueryImsDataParam
setClass("QueryImsDataParam",
         slots = c(
           rerun = "logical"
         ),
         contains = "Met4DXParam"
)

setGeneric("MatchBetweenRuns_DTIM_IOI",
           function(object, param, calibration_table_path, ion_of_interest_path)
             standardGeneric("MatchBetweenRuns_DTIM_IOI"))


# QueryImsData
setGeneric("QueryImsData",
           function(object, param, ...)
             standardGeneric("QueryImsData"))


# QueryImsmsData
setGeneric("QueryImsmsData",
           function(object, param, ...)
             standardGeneric("QueryImsmsData"))


setGeneric("ExtractIMMSMS_DTIM",
           function(object, param, ...)
             standardGeneric("ExtractIMMSMS_DTIM"))

setClass("ExtractIMMSMSParam",
         slots = c(
           order_column = "nullOrCharacter",
           mz_tol = "numeric",
           frame_range = "numeric",
           frame_integration_range = "numeric",
           mobility_range = "numeric",
           mobility_intgration_range = "nullOrNumeric",
           min_points = "numeric",
           min_intensity = "numeric",
           n_skip = "numeric",
           interpolate_method = "nullOrCharacter",
           smooth_method = "character",
           snthreshold = "numeric",
           peak_span_eim = "numeric",
           peak_span_eic = "numeric",
           smooth_window_eim = "numeric",
           smooth_window_eic = "numeric",
           skip_invalid_eic_peaks = "logical",
           skip_invalid_eim_peaks = "logical",
           keep_profile = "logical",
           rerun = "logical"
         ),
         contains = "Met4DXParam"
)

# FinalizeFeatures_DTIM
setGeneric("FinalizeFeatures_DTIM",
           function(object, param, ...)
             standardGeneric("FinalizeFeatures_DTIM"))
# # ExtractIMDataIonList
# setGeneric("ExtractIMDataIonList",
#            function(object, param, 
#                     calibration_table_path, 
#                     ion_of_interest_path, ...)
#              standardGeneric("ExtractIMDataIonList"))


###### build a experiment for ims data ########
#' @export
setMethod(
  "ImsData",
  signature = "Experiment",
  function(experiment) {
    # browser()
    wd0 <- getwd()
    setwd(experiment@wd)
    
    # list raw data files
    raw_files <- .get_ims_files(experiment@wd)
    if (!is.null(experiment@injection_order)) {
      injection_order <- read.csv(experiment@injection_order, stringsAsFactors = FALSE)
      colnames(injection_order) <- tolower(colnames(injection_order))
      injection_order <- injection_order[order(injection_order[,'order']), ]
      sub_pattern <- ifelse(grepl('/$', experiment@wd), experiment@wd, paste0(experiment@wd, '/'))
      file_match_index <- match(injection_order[, 'file'], gsub(sub_pattern, '', raw_files))
      if (any(is.na(file_match_index))) {
        stop('Some files are not correct timsPro data. Please check file: ', experiment@injection_order)
      }
      raw_files <- raw_files[file_match_index]
    }
    
    tims_data <- new("TimsData",
                     files = raw_files,
                     experiment = experiment,
                     sample_groups = .get_sample_groups(raw_files, keep_suffix = TRUE))
    
    setwd(wd0)
    return(tims_data)
  }
)

.get_ims_files <- function(wd) {
  res <- lapply(grep("\\.mzML$", list.files(wd, recursive = TRUE), value = TRUE), function(dr) {
      return(dr)
  })
  return(do.call(c, res[!is.na(res)]))
}
##### query raw data ########
## QueryImsDataParam
#' parameters for querying IMS data
#' @export
QueryImsDataParam <- function(
  rerun = FALSE
) {
  new("QueryImsDataParam",
      rerun = rerun
  )}


#' @export
setMethod(
  "QueryImsData",
  signature = c("TimsData", "QueryImsDataParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    
    message("Querying ims data...")
    
    param_list <- as.list(param)
    # calibration_table <- read.csv(calibration_table_path, stringsAsFactors = FALSE)
    # ion_of_interest <- read.csv(ion_of_interest_path, stringsAsFactors = FALSE)
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
      .parallel_parser(".query_ims_data", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))
    names(tims_data_files) <- object@files
    object@tmp_data_files$tims_data_files <- tims_data_files
    
    setwd(wd0)
    return(object)
  })


.query_ims_data <- function(
    data_file,
    ...
) {
  raw_file <- mzR::openMSfile(data_file, backend = 'pwiz')
  spec_info <- mzR::header(raw_file)
  spec <- mzR::spectra(raw_file)
  
  new_spec <- lapply(seq(nrow(spec_info)), function(i){
    temp_spec <- as.data.frame(spec[[i]])
    temp_spec$ionMobilityDriftTime <- spec_info$ionMobilityDriftTime[i]
    temp_spec$retentionTime <- spec_info$retentionTime[i]
    return(temp_spec)
  })
  # browser()
  ms1_info <- spec_info[which(spec_info$msLevel == 1), ]
  all_unique_rt <- unique(ms1_info$retentionTime)
  
  ms1_info_combined <- lapply(seq(length(all_unique_rt)), function(i){
    temp_rt <- all_unique_rt[i]
    temp_ms1_info <- ms1_info[ms1_info$retentionTime == all_unique_rt[i], ]
    
    spec <- do.call(rbind, new_spec[temp_ms1_info$seqNum])
    final <- list()
    final$temp_ms1_info <- temp_ms1_info[1, ]
    final$spec <- spec
    return(final)
  })
  
  final_ms1_info <- lapply(ms1_info_combined, function(i){
    i$temp_ms1_info
  })
  final_ms1_info <- do.call(rbind, final_ms1_info)
  row.names(final_ms1_info) <- paste0('#', row.names(final_ms1_info))
  final_spec <- lapply(ms1_info_combined, function(i){
    i$spec
  })
  names(final_spec) <- row.names(final_ms1_info)
  
  template_data <- list()
  template_data$ms1_frame_info <- final_ms1_info[, c('seqNum', 'retentionTime')]
  row.names(template_data$ms1_frame_info) <- NULL
  
  all_frames <- final_spec
  
  all_mobility <- unique(ms1_info$ionMobilityDriftTime)
  all_mobility <- all_mobility[order(all_mobility)]
  names(all_mobility) <- seq_along(all_mobility)
  
  all_frames <- lapply(seq(length(all_frames)), function(i){
    temp_spec <- all_frames[[i]]
    temp_spec$scan <- names(all_mobility)[match(temp_spec$ionMobilityDriftTime, all_mobility)]
    temp_spec <- temp_spec[, c('scan', 'mz', 'intensity')]
    temp_spec$scan <- as.integer(temp_spec$scan)
    temp_spec$mz <- as.numeric(temp_spec$mz)
    temp_spec$intensity <- as.numeric(temp_spec$intensity)
    temp_spec <- temp_spec[which(temp_spec$intensity > 0), ]
    return(temp_spec)
  })
  
  names(all_frames) <- template_data$ms1_frame_info$seqNum
  template_data$all_frames <- all_frames
  template_data$all_mobility <- all_mobility
  template_data$precursor_info <- data.frame()
  colnames(template_data$ms1_frame_info) <- c('Id', 'Time')
  row.names(template_data$ms1_frame_info) <- template_data$ms1_frame_info$Id
  
  # browser()
  # res <- list("ms1_frame_info" = ms1_frame_info,
  #             "all_frames" = all_frames,
  #             "all_mobility" = all_mobility,
  #             "precursor_info" = precursor_info)
  
  res <- template_data
  # opentimsr::CloseTIMS(D)
  return(res)
}


###### core methods and function to extract the raw data #############

#### calculate the ccs to drift time with provide tfix and beta ######
## .CCS2DT
#' Convert CCS to drift time
.CCS2DT <- function(ccs, mz, 
                    tfix, beta, 
                    mass_n2 = 28.006148, 
                    ...){
  gama <- (mz/(mz+mass_n2))^0.5
  dt <- (ccs * gama * beta) + tfix
  return(dt)
}
## .update_ioi_with_calibration_table
#' update ioi table with calibration table
.update_ioi_with_calibration_table <- function(file_name, 
                                               calibration_table, 
                                               ioi_table){
  idx_file <- match(file_name, calibration_table$file)
  tfix <- calibration_table$tfix[idx_file]
  beta <- calibration_table$beta[idx_file]
  # row.names(ioi_table) <- paste0('#', seq(nrow(ioi_table)))
  ioi_table$mobility <- sapply(seq(nrow(ioi_table)), function(ttt){
    .CCS2DT(ccs = ioi_table$ccs[ttt], 
            mz = ioi_table$mz[ttt],
            tfix = tfix, 
            beta = beta
            )
  },simplify = TRUE)
  return(ioi_table)
}

# .get_peaks_dtim <- function(
#     spectra_file,
#     im_data_file,
#     smp_idx,
#     tfix, 
#     beta,
#     ...
# ) {
#   # browser()
#   im_data <- readRDS(im_data_file)
#   im_data <- im_data[!sapply(im_data, is.null)]
#   peaks <- do.call(rbind, lapply(im_data, GetPeakInfo))
#   # peak_list <- lapply(im_data, slot, "info")
#   #
#   # peaks <- do.call(rbind, peak_list)
#   # info <- readRDS(spectra_file)@info
#   # colnames(info)[match(c("mz", "rt"), colnames(info))] <- c("target_mz", "target_rt")
#   # peaks <- cbind(info[match(rownames(peaks), rownames(info)),], peaks)
#   peaks$ccs <- .dt2ccs(dt = peaks$mobility, mz = peaks$mz, tfix = tfix, beta = beta)
#   peaks$spec_idx <- rownames(peaks)
#   peaks$smp_idx <- smp_idx
#   return(peaks)
# }

## .dt2ccs
#' calculate drift time into ccs
.dt2ccs <- function(dt, mz, tfix, beta, 
                    mass_n2 = 28.006148, 
                    ...){
  gama <- (mz/(mz+mass_n2))^0.5
  # dt <- (ccs * gama * beta) + tfix
  ccs <- (dt - tfix)/(gama * beta)
  return(ccs)
}

###### match between run with DTIM #######
#' @export
setMethod(
  "MatchBetweenRuns_DTIM_IOI",
  signature = c("TimsData", "MatchBetweenRunParam", 'ANY', 'ANY'),
  function(object, param, calibration_table_path, ion_of_interest_path) {
    # browser()
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    message("Matching features between runs...")
    
    sub_dir <- ifelse(is.null(param@interpolate_method), "match_between_runs",
                      file.path("match_between_runs", param@interpolate_method))
    files <- .tmp_files(object@files, object@experiment@tmp_dir, sub_dir)
    names(files) <- object@files
    
    ### read ioi and converted ccs into mobility ####
    # ioi <- read.csv(ion_of_interest_path, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
    # ioi$mobility <- .ccs2mobility(ccs = ioi$ccs, mz = ioi$mz)
    # ioi$ccs <- NULL
    
    calibration_table <- read.csv(calibration_table_path, header = TRUE, stringsAsFactors = FALSE)
    object@tmp_data_files$calibration_table <- calibration_table
    ioi <- read.csv(ion_of_interest_path, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
    
    ioi_as_feature <- ioi
    ioi_as_feature$mobility <- .CCS2DT(ccs = ioi$ccs,
                                       mz = ioi$mz, 
                                       tfix = calibration_table$tfix[1], 
                                       beta = calibration_table$beta[1])
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
        # browser()
        ioi_for_this_file <- .update_ioi_with_calibration_table(file_name = object@files[idx],
                                                  calibration_table = calibration_table,
                                                  ioi_table = ioi)
        
        c(list("tims_data_file" = unname(object@tmp_data_files$tims_data_files[idx]),
               "info" = ioi_for_this_file,
               "precursor_bin_file" = NULL,
               "res_define_at" = object@experiment@res_define_at),
          param_list)
      })
      .parallel_parser(".match_between_runs_dtim_ioi", arg_list, files[idxs],
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
        
        # list("smp_index" = idx,
        #      "features" = object@features,
        #      "match_between_run_file" = unname(object@tmp_data_files$match_between_run_files[idx]),
        #      "rt_model" = rt_model,
        #      "rerun" = param@rerun)
        
        data_file <- object@files[idx]
        idx_file <- match(data_file, calibration_table$file)
        tfix <- calibration_table$tfix[idx_file]
        beta <- calibration_table$beta[idx_file]
        # c("spectra_file" = unname(object@tmp_data_files$spectra_files[data_file]),
        #   "im_data_file" = unname(object@tmp_data_files$im_data_files[data_file]),
        #   "smp_idx" = idx,
        #   "tmp_dir" = unname(files[data_file]),
        #   'tfix' = tfix,
        #   'beta' = beta,
        #   param_list)
        list("smp_index" = idx,
             "features" = object@features,
             "match_between_run_file" = unname(object@tmp_data_files$match_between_run_files[idx]),
             "rt_model" = rt_model,
             'tfix' = tfix, 
             'beta' = beta,
             "rerun" = param@rerun)
        
        
      })
      .parallel_parser(".get_matched_peaks_dtim", arg_list, files[idxs],
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

.match_between_runs_dtim_ioi <- function(
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
  im_data <- .extract_im_data_dtim_roi(
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


.extract_im_data_dtim_roi <- function(
    info,
    tims_data_file,
    precursor_bin_file = NULL,
    order_column = "target_intensity",
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
    allowed_mobility_shift = 0.5,
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
                           min_points = 5,
                           min_intensity = 0,
                           n_skip = 1,
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


.get_matched_peaks_dtim <- function(
    smp_index,
    features,
    match_between_run_file,
    rt_model = NULL,
    data_file = NULL,
    tfix, 
    beta, 
    ...
) {
  fill_gaps_data <- readRDS(match_between_run_file)
  fill_gaps_data <- fill_gaps_data[!sapply(fill_gaps_data, is.null)]
  # peak_quality <- data.frame(t(sapply(fill_gaps_data, slot, "peak_quality")), stringsAsFactors = FALSE)
  # res <- data.frame(t(sapply(fill_gaps_data, slot, "info")), stringsAsFactors = FALSE)
  res <- do.call(rbind, lapply(fill_gaps_data, GetPeakInfo))
  res$ccs <- .dt2ccs(dt = res$mobility, mz = res$mz, tfix = tfix, beta = beta)
  res$smp_idx <- smp_index
  fts_filled <- rownames(features)[rownames(features) %in% rownames(res)]
  res <- res[fts_filled, , drop = FALSE]
  res$feature_idx <- fts_filled
  res$rt_align <- splus2R::ifelse1(is.null(rt_model), features[rownames(res), "rt"], predict(rt_model, res[, "rt"]))
  if (!is.null(data_file)) {
    saveRDS(res, file = data_file, version = 2)
  } else {
    return(res)
  }
}



##### extract IMMS2 ######
#' @export
setMethod(
  "QueryImsmsData",
  signature = c("TimsData", "QueryImsDataParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    
    message("Querying imsms data...")
    
    param_list <- as.list(param)
    # calibration_table <- read.csv(calibration_table_path, stringsAsFactors = FALSE)
    # ion_of_interest <- read.csv(ion_of_interest_path, stringsAsFactors = FALSE)
    par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)
    
    files <- .tmp_files(object@files, object@experiment@tmp_dir, 'tims_data_ms2')
    names(files) <- object@files
    
    tims_data_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(object@files[idxs], function(data_file) {
        c("data_file" = data_file,
          param_list)
      })
      .parallel_parser(".query_imsms_data", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))
    names(tims_data_files) <- object@files
    object@tmp_data_files$tims_data_files_ms2 <- tims_data_files
    
    setwd(wd0)
    return(object)
  })


.query_imsms_data <- function(
    data_file,
    ...
) {
  # data_file <- './S1_MA-d3-c3_SR.mzML'
  raw_file <- mzR::openMSfile(data_file, backend = 'pwiz')
  spec_info <- mzR::header(raw_file)
  spec <- mzR::spectra(raw_file)
  
  new_spec <- lapply(seq(nrow(spec_info)), function(i){
    temp_spec <- as.data.frame(spec[[i]])
    temp_spec$ionMobilityDriftTime <- spec_info$ionMobilityDriftTime[i]
    temp_spec$retentionTime <- spec_info$retentionTime[i]
    return(temp_spec)
  })
  
  ms1_info <- spec_info[which(spec_info$msLevel == 2), ]
  all_unique_rt <- unique(ms1_info$retentionTime)
  
  ms1_info_combined <- lapply(seq(length(all_unique_rt)), function(i){
    temp_rt <- all_unique_rt[i]
    temp_ms1_info <- ms1_info[ms1_info$retentionTime == all_unique_rt[i], ]
    
    spec <- do.call(rbind, new_spec[temp_ms1_info$seqNum])
    final <- list()
    final$temp_ms1_info <- temp_ms1_info[1, ]
    final$spec <- spec
    return(final)
  })
  
  final_ms1_info <- lapply(ms1_info_combined, function(i){
    i$temp_ms1_info
  })
  final_ms1_info <- do.call(rbind, final_ms1_info)
  row.names(final_ms1_info) <- paste0('#', row.names(final_ms1_info))
  final_spec <- lapply(ms1_info_combined, function(i){
    i$spec
  })
  names(final_spec) <- row.names(final_ms1_info)
  
  template_data <- list()
  template_data$ms1_frame_info <- final_ms1_info[, c('seqNum', 'retentionTime')]
  row.names(template_data$ms1_frame_info) <- NULL
  
  all_frames <- final_spec
  
  all_mobility <- unique(ms1_info$ionMobilityDriftTime)
  all_mobility <- all_mobility[order(all_mobility)]
  names(all_mobility) <- seq_along(all_mobility)
  
  all_frames <- lapply(seq(length(all_frames)), function(i){
    temp_spec <- all_frames[[i]]
    temp_spec$scan <- names(all_mobility)[match(temp_spec$ionMobilityDriftTime, all_mobility)]
    temp_spec <- temp_spec[, c('scan', 'mz', 'intensity')]
    temp_spec$scan <- as.integer(temp_spec$scan)
    temp_spec$mz <- as.numeric(temp_spec$mz)
    temp_spec$intensity <- as.numeric(temp_spec$intensity)
    temp_spec <- temp_spec[which(temp_spec$intensity > 0), ]
    return(temp_spec)
  })
  
  names(all_frames) <- template_data$ms1_frame_info$seqNum
  
  ### add msms info 
  IsolationMz <-  mean(c(min(spec_info$lowMZ), max(spec_info$highMZ)))
  IsolationWidth <- max(spec_info$highMZ) - IsolationMz
  msms_info <- data.frame(WindowGroup = 1, 
                          Frame = template_data$ms1_frame_info$seqNum, 
                          ScanNumBegin = as.integer(names(all_mobility)[1]),
                          ScanNumEnd = as.integer(names(all_mobility)[length(all_mobility)]), 
                          IsolationMz = IsolationMz, 
                          IsolationWidth = IsolationWidth, 
                          CollisionEnergy = 20)
  
  template_data$msms_info <- msms_info
  template_data$all_frames <- all_frames
  template_data$all_mobility <- all_mobility
  # template_data$precursor_info <- data.frame()
  colnames(template_data$ms1_frame_info) <- c('Id', 'Time')
  row.names(template_data$ms1_frame_info) <- template_data$ms1_frame_info$Id
  names(template_data)[1] <- 'frame_info'
  
  # res <- list("ms1_frame_info" = ms1_frame_info,
  #             "all_frames" = all_frames,
  #             "all_mobility" = all_mobility,
  #             "precursor_info" = precursor_info)
  
  res <- template_data
  # opentimsr::CloseTIMS(D)
  return(res)
}





####### extract ms2 spectar ####
## ExtractIMMSMSParam
#' parameters for extracting MSMS spectra from IM-MSMS Data
#' @export
ExtractIMMSMSParam <- function(
    order_column = "target_intensity",
    mz_tol = 20,
    frame_range = 30,
    frame_integration_range = 5,
    mobility_range = 0.1,
    mobility_intgration_range = 0.015,
    min_intensity = 0,
    n_skip = 0,
    interpolate_method = NULL,
    smooth_method = "loess",
    snthreshold = 3,
    min_points = 2,
    peak_span_eim = 21,
    peak_span_eic = 7,
    smooth_window_eim = 16,
    smooth_window_eic = 8,
    skip_invalid_eic_peaks = FALSE,
    skip_invalid_eim_peaks = TRUE,
    keep_profile = FALSE,
    rerun = FALSE
) {
  if (!is.null(interpolate_method)) {
    interpolate_method <- match.arg(interpolate_method, c("linear", "scans"))
  }
  new("ExtractIMMSMSParam",
      order_column = order_column,
      mz_tol = mz_tol,
      frame_range = frame_range,
      frame_integration_range = frame_integration_range,
      mobility_range = mobility_range,
      mobility_intgration_range = mobility_intgration_range,
      min_points = min_points,
      min_intensity = min_intensity,
      n_skip = n_skip,
      interpolate_method = interpolate_method,
      smooth_method = smooth_method,
      snthreshold = snthreshold,
      peak_span_eim = peak_span_eim,
      peak_span_eic = peak_span_eic,
      smooth_window_eim = smooth_window_eim,
      smooth_window_eic = smooth_window_eic,
      skip_invalid_eic_peaks = skip_invalid_eic_peaks,
      skip_invalid_eim_peaks = skip_invalid_eim_peaks,
      keep_profile = keep_profile,
      rerun = rerun
  )
}




#' @export
setMethod(
  "ExtractIMMSMS_DTIM",
  signature = c("TimsData", "ExtractIMMSMSParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)

    message("Extracting IM Spectra...")

    # extract MSMS for target precursors
    par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)
    files <- .tmp_files(object@files, object@experiment@tmp_dir, 'spectra')
    names(files) <- object@files

    param_list <- as.list(param)
    spec_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(idxs, function(smp_idx) {
        # browser()
        c(list("precursor_info" = object@filled_peaks[object@filled_peaks$smp_idx == smp_idx, , drop = FALSE],
               "tims_data_file" = unname(object@tmp_data_files$tims_data_files_ms2[object@files[smp_idx]])),
          param_list)
      })
      .parallel_parser('.extract_im_msms_dtim',
                       arg_list, files[idxs],
                       object@experiment@BPPARAM,
                       save_in_analysis = TRUE)
    }, simplify = FALSE))

    names(spec_files) <- object@files
    object@tmp_data_files$spectra_files <- spec_files

    setwd(wd0)
    return(object)
  })


.extract_im_msms_dtim <- function(
    precursor_info,
    tims_data_file,
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
    smooth_method = 'loess',
    skip_invalid_eic_peaks = TRUE,
    skip_invalid_eim_peaks = TRUE,
    filter_outlier_peaks = FALSE,
    allowed_mobility_shift = 0.015,
    allowed_rt_shift = 10,
    data_file = NULL,
    ...
) {
  query_data <- readRDS(tims_data_file)
  query_window <- seq(unique(query_data$msms_info$WindowGroup))

  scans <- names(query_data$all_mobility)
  num_scan <- length(scans)
  target_scan_range <- t(sapply(precursor_info$mobility, function(x) {
    mobility_diff <- abs(query_data$all_mobility - x)
    c(range(as.numeric(scans[mobility_diff <=  mobility_range])),
      splus2R::ifelse1(is.null(mobility_intgration_range),
                       c(NA, NA),
                       range(as.numeric(scans[mobility_diff <=  mobility_intgration_range]))),
      as.numeric(scans[which.min(mobility_diff)])
    )
  }))
  colnames(target_scan_range) <- col_ranges <- c('scan_from', 'scan_to', 'scan_int_from', 'scan_int_to', 'scan_ref')
  precursor_info <- cbind(precursor_info, target_scan_range)
  res <- apply(precursor_info[, c('mz', 'rt', 'mobility', col_ranges)], 1, function(dr) {
    # browser()
    mz <- dr['mz']
    rt <- dr['rt']
    scan_from <- dr['scan_from']
    scan_to <- dr['scan_to']
    scan_ref <- dr['scan_ref']
    scan_int_from <- dr['scan_int_from']
    scan_int_to <- dr['scan_int_to']

    target_frames <- query_data$frame_info[sum(query_data$frame_info$Time < rt) + query_window, "Id"]
    target_window <- query_data$msms_info[query_data$msms_info$Frame %in% target_frames, , drop = FALSE]
    mz_diff <- abs(target_window$IsolationMz - mz)
    if (all(mz_diff > target_window$IsolationWidth)) {
      return()
    }
    scan_range <- target_window[which.min(mz_diff), c('ScanNumBegin', 'ScanNumEnd', 'Frame')]
    # # lmd 0927
    # scan_range <- data.frame(ScanNumBegin = as.integer(rep(scan_range$ScanNumBegin, 5)),
    #                          ScanNumEnd = as.integer(rep(scan_range$ScanNumEnd, 5)),
    #                          Frame = as.integer(seq(-2, 2, 1) * length(query_window) + scan_range$Frame))
    # scan_range <- scan_range[scan_range$Frame > as.integer(names(query_data$all_frames)[1]) & scan_range$Frame < as.integer(names(query_data$all_frames)[length(query_data$all_frames)]), ]
    target_spec <- query_data$all_frames[[as.character(scan_range$Frame)]]
    # lmd 0927
    # target_spec <- do.call(rbind, lapply(scan_range$Frame, function(fr){
    #   query_data$all_frames[[as.character(fr)]]
    # }))
    # spec_profile <- target_spec[target_spec$scan >= max(scan_range$ScanNumBegin, scan_from) &
    #                               target_spec$scan <= min(scan_range$ScanNumEnd, scan_to) &
    #                               target_spec$mz <= .ppm2dalton(mz, mz_tol, res_define_at) + mz,
    #                             c('mz', 'intensity', 'inv_ion_mobility', 'scan')]
    spec_profile <- target_spec[target_spec$scan >= max(scan_range$ScanNumBegin, scan_from) &
                                  target_spec$scan <= min(scan_range$ScanNumEnd, scan_to) &
                                  target_spec$mz <= .ppm2dalton(mz, mz_tol, res_define_at) + mz,
                                c('mz', 'intensity', 'scan')]
    if (nrow(spec_profile) == 0) {
      return()
    }
    rois <- .get_rois(spec_profile,
                      mz_tol = mz_tol,
                      min_points = min_points,
                      min_intensity = min_intensity,
                      n_skip = 3,
                      roi_along = 'scan',
                      res_define_at = res_define_at,
                      ref_scan = scan_ref)
    if (length(rois) == 0) {
      return(NULL)
    }
    if (is.null(mobility_intgration_range)) {
      spec <- do.call(rbind, lapply(rois, function(roi) {
        roi_spec <- roi[roi[, 'scan'] == scan_ref, c('mz', 'intensity')]
        if (nrow(roi_spec) > 1) {
          roi_spec <- roi_spec[which.min(abs(roi_spec[, 'mz'] - .get_weighted_mz(roi[, 'mz'], roi[, 'intensity']))), ]
        }
        roi_spec
      }))
    } else {
      spec <- do.call(rbind, lapply(rois, function(roi) {
        roi <- roi[roi[, 'scan'] >= scan_int_from & roi[, 'scan'] <= scan_int_to, , drop = FALSE]
        mz <- .get_weighted_mz(roi[, 'mz'], roi[, 'intensity'])
        intensity <- sum(roi[, 'intensity'])
        c(mz, intensity)
      }))
    }
    colnames(spec) <- c('mz', 'intensity')
    spec <- spec[order(spec[, 'mz']), , drop = FALSE]
    rownames(spec) <- NULL
    return(spec)
  }, simplify = FALSE)
  # browser()
  names(res) <- row.names(precursor_info)
  res <- res[!sapply(res, is.null)]
  if (!is.null(data_file)) {
    saveRDS(res, file = data_file, version = 2)
  } else {
    return(res)
  }
}





.get_rois <- function(
    data,
    mz_tol,
    min_points = 4,
    min_intensity = 0,
    n_skip = 0,
    roi_along = 'scan',
    res_define_at = 200,
    ref_scan = NULL
) {
  data <- data[order(data[, "intensity"], decreasing = TRUE), , drop = FALSE]
  tmp_data <- cbind(seq(min(data[, roi_along]), max(data[, roi_along])), 0)
  rownames(tmp_data) <- tmp_data[, 1]
  rois <- vector('list', nrow(data))
  idx <- 1
  if (is.null(ref_scan)) {
    ref_index <- sapply(data[, roi_along], function(x) {
      which.min(abs(tmp_data[, 1] - x))
    })
  } else {
    ref_index <- rep(which(tmp_data[,1] == ref_scan), nrow(data))
  }
  if (length(ref_index)) { # incase that ref_scan is out of tmp_data
    while(nrow(data)) {
      mz <- data$mz[1]
      is_mz_match <- abs(data$mz - mz) <= .ppm2dalton(mz, mz_tol, res_define_at)
      sort(data[is_mz_match, roi_along])
      roi_data <- tmp_data
      roi_data[as.character(data[is_mz_match, roi_along]), 2] <-  data[is_mz_match, 'intensity']
      roi <- .find_roi(roi_data[, 2],
                       min_points = min_points,
                       min_intensity = min_intensity,
                       n_skip = n_skip,
                       ref = ref_index[1]
      )
      if (!is.null(roi) && nrow(roi) > 0) {
        in_roi <- is_mz_match & data[, 'scan'] %in% roi_data[seq(roi[1, 1], roi[1, 2]), 1]
        roi_data <- data[in_roi, , drop = FALSE]
        rois[[idx]] <- roi_data[order(roi_data[, roi_along]), , drop = FALSE]
        data <- data[!in_roi, ]
        ref_index <- ref_index[!in_roi]
      } else {
        data <- data[-1, , drop = FALSE]
        ref_index <- ref_index[-1]
      }
      idx <- idx + 1
    }
  }
  rois <- rois[!sapply(rois, is.null)]
  return(rois)
}


##### finalize data #####
#' @export
setMethod(
  "FinalizeFeatures_DTIM",
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
                       "spectra_files" = object@tmp_data_files$spectra_files,
                       "peak_groups" = object@peak_groups,
                       "res_define_at" = object@experiment@res_define_at,
                       'tfix' = object@tmp_data_files$calibration_table$tfix[1],
                       'beta' = object@tmp_data_files$calibration_table$beta[1],
                       "bpparam" = splus2R::ifelse1(par_analysis, object@experiment@BPPARAM, NULL)),
                  param_list)
    
    finalize_data_file <- .analysis_parser(".finalize_data_dtim", arg_list, files)
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

.finalize_data_dtim <- function(
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
    res_define_at = 200,
    bpparam = NULL,
    tfix, 
    beta,
    ...
) {
  # browser()
  # features <- features[, -match(c("npeaks", "area"), colnames(features)), drop = FALSE]
  
  sample_classes <- .sample_groups_to_class(sample_groups)
  
  spectra_all <- .get_peak_spectra_dtim(filled_peaks, spectra_files, update_info = TRUE)
  
  filled_groups <- sapply(unique(filled_peaks$feature_idx), function(feature_idx) {
    rownames(filled_peaks[filled_peaks$feature_idx == feature_idx, "smp_idx", drop = FALSE])
  }, simplify = FALSE)
  all_groups <- sapply(rownames(features), function(nm) {
    c(peak_groups[[nm]], filled_groups[[nm]])
  }, simplify = FALSE)
  
  # cols_keep <- intersect(colnames(peaks), colnames(filled_peaks))
  # all_peaks <- rbind(peaks[, cols_keep], filled_peaks[, cols_keep])
  all_peaks <- filled_peaks
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
      all_peaks[na.omit(fastmatch::fmatch(grp, all_peak_names)), , drop = FALSE]
    })
  } else {
    group_sample_peaks <- BiocParallel::bplapply(all_groups, function(grp) {
      all_peaks[na.omit(fastmatch::fmatch(grp, all_peak_names)), , drop = FALSE]
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
  # browser()
  spectra_assign_nms <- sapply(rownames(features), function(nm) {
    # browser()
    # nms_spec <- intersect(peak_groups[[nm]], rownames(group_sample_peaks[[nm]]))
    nms_spec <- rownames(group_sample_peaks[[nm]])
    idx_max <- which.max(spectra_all[nms_spec]@info[[col_max]])
    if(length(idx_max) == 0){
      return(NA)
    }else{
      return(nms_spec[idx_max])
    }
  }, simplify = TRUE)
  
  peak_groups <- peak_groups[is_keep]
  features$ccs <- .dt2ccs(features$mobility, features$mz, tfix = tfix, beta = beta)
  features$name <- .get_peak_name(features[, c("mz", "rt", "ccs")])
  features <- features[, c("name", colnames(features)[1:(ncol(features) - 1)])]
  
  smp_idxs <- seq(nrow(sample_groups))
  # all_peaks <- rbind(peaks[, cols_keep], filled_peaks[, cols_keep])
  all_peaks <-  filled_peaks
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
  spec <- spectra_all[spectra_assign_nms]
  spec <- SpectraTools::UpdateNames(spec, rownames(features)[which(row.names(features) %in% spec@info$feature_idx)])
  # browser()
  return(list("features" = cbind(features, quant_data), "spectra" = spec, "peak_groups" = peak_groups))
}


.get_peak_spectra_dtim <- function(peaks, spectra_files, min_num_fragments = NULL, update_info = FALSE) {
  spectra_data <- lapply(seq_along(spectra_files), function(idx) {
    readRDS(spectra_files[idx])
  })
  if (update_info) {
    info <- peaks
  } else {
    info <- do.call(rbind, lapply(spectra_data, function(x) {
      x@info
    }))
    rownames(info) <- .gen_indexes(info)
  }
  
  spec <- do.call(c, lapply(seq_along(spectra_files), function(idx) {
    spectra_data[[idx]]
  }))
  idxxxxx <- fastmatch::fmatch(row.names(info), names(spec))
  spec <- spec[idxxxxx]
  names(spec) <- rownames(info)
  if (!is.null(min_num_fragments)) {
    is_keep <- sapply(spec, nrow) >= min_num_fragments
    info <- info[is_keep, , drop = FALSE]
    spec <- spec[is_keep]
  }
  SpectraTools::SpectraData(info = info,
                            spectra = spec)
}

##### fill gaps #########













