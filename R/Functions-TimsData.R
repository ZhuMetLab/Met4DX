# project: Met4DX
# File name: Functions-TimsData.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/1/17 12:01
# Copyright (c) 2022- ZhuMSLab ALL right reserved

.read_spectra <- function(
    data_file,
    intensity_from = c("pepmass", "ms2_intensity"),
    include_precursor = TRUE,
    mz_tol = 20,
    denoise = TRUE,
    int_threshold_abs = 30,
    int_threshold_rel = 0.01,
    ms2range = NULL,
    res_define_at = 200,
    ...
) {
  intensity_from <- match.arg(intensity_from)
  param <- SpectraTools::ParseSpectraParam(
    type = "mgf",
    denoise = FALSE,
    ms2range = ms2range,
    includePrecursor = include_precursor,
    ppmPrecursorFilter = mz_tol,
    thrIntensityAbs = int_threshold_abs,
    thrIntensityRel = int_threshold_rel,
    labelKeep = c("PEPMASS", "RTINSECONDS", "TITLE", "RAWSCANS"),
    labelName = c("precursor_info", "rt", "info", "raw_scans"),
    resDefineAt = res_define_at
  )
  spectra <- SpectraTools::ParseSpectra(param, data_file)
  mz_info_names <- c('mz', 'intensity')
  info <- t(apply(spectra@info, 1, function(dr) {
    .info <- dr["info"]
    one_over_k0 <- regmatches(.info, regexpr('(?<=1/K0=)\\d+\\.\\d+(?=,)', .info, perl = TRUE))
    cmpd <- regmatches(.info, regexpr('(?<=Cmpd\\s)\\d+(?=,)', .info, perl = TRUE))
    scans <- strsplit(dr["raw_scans"], "-|,")[[1]]
    mz_info <- head(strsplit(dr["precursor_info"], "\t| ")[[1]], 2)
    names(mz_info) <- mz_info_names[seq_along(mz_info)]
    c(dr,
      mz_info,
      "k0" = unname(one_over_k0),
      "cmpd" = unname(cmpd),
      "ms1_scans" = scans[1],
      "ms2_scans" = paste0(scans[-1], collapse = ","))
  }))
  spectra@rtInfo <- NULL
  spectra@ccsInfo <- NULL
  spectra@info <- SpectraTools:::.Col2Numeric(as.data.frame(info, stringsAsFactors = FALSE))
  if (denoise) {
    spectra <- SpectraTools::DenoiseSpectra(spectra,
                                            resetIdx = TRUE,
                                            ms2range = ms2range,
                                            includePrecursor = include_precursor,
                                            ppmPrecursorFilter = mz_tol,
                                            detectIntensityAbs = TRUE,
                                            thrDetectCounts = 3,
                                            thrIntensityAbs = int_threshold_abs,
                                            thrIntensityRel = int_threshold_rel,
                                            resDefineAt = res_define_at)
  }
  
  if (intensity_from == "ms2_intensity") {
    spectra@info[, "intensity"] <- unname(sapply(spectra@spectra, function(spec) {
      sum(sort(spec$intensity, decreasing = TRUE)[1:min(10, nrow(spec))])
    }))
  }
  return(spectra)
}

.bin_precursors <- function(
  data_file,
  mz_tol = 20,
  rt_tol = c(10, 20),
  mobility_tol = c(0.015, 0.03),
  distance_cutoff = 1,
  weight_rt = 1,
  weight_mobility = 1,
  weight_msms = 1,
  ms2_range = NULL,
  method_match = c("direct", "bootstrap"),
  method_score = c("dp", "msdial", ""),
  weight_intensity = 1,
  weight_mz = 0,
  int_threshold_abs = 0,
  int_threshold_rel = 0,
  include_precursor = TRUE,
  mz_tol_ms1 = 20,
  mz_tol_ms2 = 35,
  res_define_at = 200,
  ...
) {
  if (!file.exists(data_file)) {
    stop("File not found: %s", data_file)
  }

  spec_data <- readRDS(data_file)
  info <- spec_data@info

  info <- info[order(info$intensity, decreasing = TRUE), , drop = FALSE]
  if (mz_tol > 1) {
    info$mz_tol <- sapply(info$mz, function(mz) .ppm2dalton(mz, mz_tol, res_define_at))
  } else {
    info$mz_tol <- mz_tol
  }

  # create bin pool and init with highest precursor to store target precursors
  bin_pool <- data.frame(info[1, c("mz", "rt", "k0")])

  # create cadidates list to store all precursors
  cad_list <- vector("list")
  # init with the highest precursor
  cad_list[[1]] <- rownames(info)[1]

  # traversal all precursors and add it to corresponding bin pool and cadidates list
  for (nm in rownames(info)[-1]) {
    dr <- info[nm,]
    mobility_diff <- abs(bin_pool$k0 - dr$k0)
    cad_idx <- which(abs(bin_pool$mz - dr$mz) <= dr$mz_tol &
                       abs(bin_pool$rt - dr$rt) <= rt_tol[2] &
                       mobility_diff <= mobility_tol[2])
    if (length(cad_idx) == 0) {
      bin_pool <- rbind(bin_pool, dr[c("mz", "rt", "k0")])
      cad_list <- c(cad_list, nm)
    } else {
      best_idx <- cad_idx[which.min(mobility_diff[cad_idx])]
      cad_list[[best_idx]] <- c(cad_list[[best_idx]], nm)
    }
  }

  match_param <- SpectraTools::MatchParam(ppm = mz_tol_ms2,
                                          cutoff = 0,
                                          methodMatch = method_match,
                                          methodScore = method_score,
                                          includePrecursor = include_precursor,
                                          ppmPrecursorFilter = mz_tol_ms1,
                                          ms2range = ms2_range,
                                          intensityExpNormed = FALSE,
                                          intensityLibNormed = FALSE,
                                          tuneLibSpectra = TRUE,
                                          thrIntensityAbs = int_threshold_abs,
                                          thrIntensityRel = int_threshold_rel,
                                          weightIntensity = weight_intensity,
                                          weightMZ = weight_mz,
                                          useMS1ResDefine = TRUE,
                                          resDefineAt = res_define_at)
  # cad_len <- sort(unique(sapply(cad_list, length)))

  weights <- c(weight_rt, weight_mobility, weight_msms)

  idx_remove <- vector("numeric")

  pb <- txtProgressBar(min = 0, max = length(cad_list), initial = 0, style = 3)

  # iterate all binned cadiates and cluster the cadidates according their distances to each other
  for (idx in seq_along(cad_list)) {
    cads <- cad_list[[idx]]
    setTxtProgressBar(pb, idx)
    if (length(cads) > 1) {
      num_cads <- length(cads)
      distance_matrix <- matrix(0, ncol = num_cads, nrow = num_cads)
      colnames(distance_matrix) <- rownames(distance_matrix) <- cads

      for (idx_cad in seq(num_cads-1)) {
        idx_cad_score <- seq(idx_cad+1, num_cads)
        cads_score <- cads[idx_cad_score]
        cad <- cads[idx_cad]
        rt_score <- SpectraTools:::.TrapezoidalScore(info[cad, "rt"], info[cads_score, "rt"],
                                                     rt_tol)
        mobility_score <- SpectraTools:::.TrapezoidalScore(info[cad, "k0"], info[cads_score, "k0"],
                                                           mobility_tol)
        msms_score <- .match_spectra(spec_data[cad],
                                     spec_data[cads_score],
                                     match_param)
        distance <- .get_peak_distance(cbind(rt_score, mobility_score, msms_score), weights)
        distance_matrix[idx_cad, idx_cad_score] <- distance_matrix[idx_cad_score, idx_cad] <- distance
      }

      if (all(distance_matrix <= distance_cutoff)) {
        # do nothing if all distances are smaller than cutoff as they are all similar
        next
      } else {
        # cluster the cadidates according their distance
        idx_remove <- c(idx_remove, idx)
        new_list <- .clust_peaks(distance_matrix, distance_cutoff)
        new_list <- lapply(new_list, function(peak_names) {
          peak_names[order(info[peak_names, "intensity"], decreasing = TRUE)]
        })
        cad_list <- c(cad_list, new_list)
      }
    }
  }
  close(pb)
  if (length(idx_remove) > 0) {
    cad_list <- cad_list[-idx_remove]
  }
  names(cad_list) <- sapply(cad_list, `[`, 1)
  return(cad_list)
}

.query_tims_data <- function(
    data_file,
    opentims_thread = 2L,
    ...
) {
  opentimsr::opentims_set_threads(opentims_thread)
  
  all_columns <- .set_opentims()
  D <- opentimsr::OpenTIMS(data_file)
  
  # precursor_info <- opentimsr::table2df(D, 'Precursors')$Precursors[, c("Id", "Parent")]
  precursor_info <- NA
  if ('Precursors' %in% opentimsr::tables_names(D)) {
    precursor_info <- opentimsr::table2df(D, 'Precursors')$Precursors[, c("Id", "Parent")]
  }
  frame_id <- opentimsr::table2df(D, 'Frames')$Frames
  
  ms1_frame_info <- frame_id[frame_id$MsMsType == 0, c("Id", "Time")]
  
  all_frames <- lapply(ms1_frame_info$Id, function(frame) {
    query(D, frames = frame, columns = c("scan", "mz", "intensity"))
  })
  names(all_frames) <- ms1_frame_info$Id
  
  all_mobility <- .query_scan_mobility(D, ms1_frame_info$Id)
  res <- list("ms1_frame_info" = ms1_frame_info,
              "all_frames" = all_frames,
              "all_mobility" = all_mobility,
              "precursor_info" = precursor_info)
  opentimsr::CloseTIMS(D)
  return(res)
}




# .extract_im_data <- function(
#   info,
#   tims_data_file,
#   precursor_bin_file = NULL,
#   order_column = "target_intensity",
#   mz_tol = 20,
#   frame_range = 30,
#   frame_integration_range = 5,
#   mobility_range = 0.1,
#   mobility_intgration_range = 0.015,
#   min_points = 4,
#   min_intensity = 0,
#   n_skip = 0,
#   interpolate_method = NULL,
#   peak_span_eim = 21,
#   peak_span_eic = 7,
#   snthreshold = 3,
#   smooth_window_eim = 16,
#   smooth_window_eic = 8,
#   keep_profile = FALSE,
#   res_define_at = 200,
#   use_cmpd_id = TRUE,
#   smooth_method = 'gaussian',
#   skip_invalid_eic_peaks = TRUE,
#   skip_invalid_eim_peaks = TRUE,
#   filter_outlier_peaks = FALSE,
#   allowed_mobility_shift = 0.015,
#   allowed_rt_shift = 10,
#   ...
# ) {
#   remove_profile <- !keep_profile
# 
#   query_data <- readRDS(tims_data_file)
#   num_frames <- nrow(query_data$frame_info)
# 
#   if (!missing(precursor_bin_file) && !is.null(precursor_bin_file)) {
#     precursor_bins <- readRDS(precursor_bin_file)
#     info <- info[sapply(precursor_bins, `[`, 1),]
#   }
# 
#   if (mz_tol > 1) {
#     info$mz_tol <- sapply(info$mz, function(mz) .ppm2dalton(mz, mz_tol, res_define_at))
#   } else {
#     info$mz_tol <- mz_tol
#   }
# 
#   mobility_col <- ifelse("k0" %in% colnames(info), "k0", "mobility")
#   intensity_col <- ifelse("intensity" %in% colnames(info), "intensity", "target_intensity")
#   info <- info[order(info[[order_column]], decreasing = TRUE),]
#   switch(smooth_method,
#          "gaussian" = {
#            smooth_param_eim <- GaussianSmoothParam(window = smooth_window_eim)
#            smooth_param_eic <- GaussianSmoothParam(window = smooth_window_eic)
#          },
#          "loess" = {
#            smooth_param_eim <- LOESSSmoothParam(window = smooth_window_eim)
#            smooth_param_eic <- LOESSSmoothParam(window = smooth_window_eic)
#          },
#          stop("Unknown smooth method")
#   )
# 
#   if (use_cmpd_id && ("cmpd" %in% colnames(info))) {
#     if (!"precursor_info" %in% colnames(query_data)) {
#       stop("Precursor info is not available for non DDA data!")
#     }
#     info$target_frame <- sapply(info[, "cmpd"], function(cmpd) {
#       query_data$precursor_info$Parent[query_data$precursor_info$Id == cmpd]
#     })
#   } else {
#     info$target_frame <- sapply(info[, "rt"], function(rt) {
#       query_data$frame_info$Id[which.min(abs(query_data$frame_info$Time - rt))]
#     })
#   }
# 
#   peak_quality_tmp <- as.list(c(rep(FALSE, 2), rep(NA, 4)))
#   names(peak_quality_tmp) <- c("eim_peak", "eic_peak", "eim_sn", "eic_sn", "eim_baseline", "eic_baseline")
# 
#   res <- apply(info[, c("mz", intensity_col, mobility_col, "mz_tol", "target_frame")], 1, function(dr) {
#     mz <- dr["mz"]
#     mobility <- dr[mobility_col]
#     mz_tol <- dr["mz_tol"]
#     target_frame <- dr["target_frame"]
# 
#     peak_quality <- peak_quality_tmp
# 
#     target_idx <- which(query_data$frame_info$Id == target_frame)
#     frame_start <- max(0, target_idx - frame_range)
#     extract_frames <- query_data$frame_info$Id[frame_start:min(num_frames, target_idx + frame_range)]
#     num_extract_frames <- length(extract_frames)
# 
#     # get eim from the target frame
#     eim <- .get_eims2(query_data$all_frames[as.character(target_frame)], query_data$all_mobility,
#                       mz, mz_tol, mobility, mobility_range)
#     # eim <- .get_eims(D, target_frame, mz, mz_tol, mobility, mobility_range, all_columns)
# 
#     if (sum(eim[, 1]) == 0) {
#       return(NULL)
#     }
# 
#     if (is.null(eim)) {
#       return(NULL)
#     }
#     colnames(eim)[1] <- "intensity"
#     eim <- .interpolate_data(eim, c(-1, 1) * mobility_range + mobility, interpolate_method)
# 
#     eim$intensity_smooth <- splus2R::ifelse1(smooth_method == "gaussian",
#                                              .smooth_gaussian(data = eim$intensity, window = smooth_param_eim@window),
#                                              .smooth_loess(data = eim$intensity, degree = smooth_param_eim@degree, window = smooth_param_eim@window)
#     )
# 
#     # if (.get_smooth_sd(eim$intensity, eim$intensity_smooth) > 0.35) {
#     #   return(NULL)
#     # }
# 
#     # deterim the eim apex for optimized eic-eim profile extraction
#     ref_index <- which.min(abs(eim$k0 - mobility))
#     apex_eim <- eim$k0[.find_apex(eim$intensity, eim$intensity_smooth,
#                                   span = peak_span_eim,
#                                   min_points = min_points,
#                                   min_intensity = min_intensity,
#                                   n_skip = n_skip,
#                                   ref_index = ref_index,
#                                   find_roi = FALSE)]
#     if (length(apex_eim) == 0) {
#       return(NULL)
#     }
# 
#     # extract eic-eim profile around the target frame and eim apex
#     eims <- .get_eims2(query_data$all_frames[as.character(extract_frames)], query_data$all_mobility,
#                        mz, mz_tol, apex_eim, mobility_range)
#     # eims <- .get_eims(D, extract_frames, mz, mz_tol, apex_eim, mobility_range, all_columns)
#     # colnames(eims)[1:num_extract_frames] <- extract_frames
#     eic_rt <- query_data$frame_info[as.character(extract_frames), "Time"]
#     im_data <- new("IMData",
#                    eic_rt = eic_rt,
#                    eic_mz = attributes(eims)$mz,
#                    eim_mobility = eims[, num_extract_frames + 1],
#                    profile_data = eims[, 1:num_extract_frames])
#     im_data <- SetData(im_data, "eic", apex_eim + c(-1, 1) * mobility_intgration_range)
#     im_data <- SmoothData(im_data, "eic", smooth_param_eic)
#     apex_eic <- .find_apex(im_data@eic$raw, im_data@eic$smooth,
#                            span = peak_span_eic,
#                            min_points = min_points,
#                            min_intensity = min_intensity,
#                            n_skip = n_skip,
#                            ref_index = frame_range + 1,
#                            force_max = TRUE)
#     if (is.null(apex_eic)) {
#       return(NULL)
#     }
# 
#     if (.get_smooth_sd(im_data@eic$raw, im_data@eic$smooth) > 0.35) {
#       if (skip_invalid_eic_peaks) {
#         return(NULL)
#       }
#     } else {
#       peak_quality$eic_peak <- TRUE
#     }
# 
#     eic_quality <- .get_peak_quality(im_data@eic$raw,
#                                      apex_eic,
#                                      min_intensity = 0,
#                                      min_points = min_points,
#                                      n_skip = n_skip,
#                                      snthreshold = snthreshold,
#                                      skip_invalid_peaks = skip_invalid_eic_peaks)
# 
#     if (is.null(eic_quality) && skip_invalid_eic_peaks) {
#       return(NULL)
#     }
# 
#     peak_quality$eic_sn <- eic_quality$snr
#     peak_quality$eic_baseline <- eic_quality$baseline
# 
#     im_data <- SetData(im_data, "eim", seq(max(1, apex_eic - 5), min(apex_eic + 5, num_extract_frames)))
#     im_data <- InterpolateEIM(im_data, c(-1, 1) * mobility_range + mobility, interpolate_method)
#     im_data <- SmoothData(im_data, "eim", smooth_param_eim)
# 
#     if (.get_smooth_sd(im_data@eim$raw, im_data@eim$smooth) > 0.35) {
#       if (skip_invalid_eim_peaks) {
#         return(NULL)
#       }
#     } else {
#       peak_quality$eim_peak <- TRUE
#     }
# 
#     ref_index <- which.min(abs(im_data@eim_mobility_interpolate - mobility))
# 
#     apex_eim_index <- .find_apex(im_data@eim$interpolate, im_data@eim$smooth,
#                                  span = peak_span_eim,
#                                  min_points = min_points,
#                                  min_intensity = min_intensity,
#                                  n_skip = n_skip,
#                                  ref_index = ref_index,
#                                  find_roi = FALSE)
#     
#     eim_quality <- .get_peak_quality(im_data@eim$interpolate,
#                                      apex_eim_index,
#                                      min_intensity = 0,
#                                      min_points = min_points, 
#                                      n_skip = n_skip,
#                                      snthreshold = snthreshold,
#                                      skip_invalid_peaks = skip_invalid_eim_peaks)
# 
#     if (is.null(eim_quality) && skip_invalid_eim_peaks) {
#       return(NULL)
#     }
# 
#     peak_quality$eim_sn <- eim_quality$snr
#     peak_quality$eim_baseline <- eim_quality$baseline
# 
#     apex_eim <- im_data@eim_mobility_interpolate[apex_eim_index]
#     if (length(apex_eim) == 0) {
#       return(NULL)
#     }
#     if (filter_outlier_peaks) {
#       if (abs(mobility - apex_eim) > allowed_mobility_shift ||
#         abs(query_data$frame_info[as.character(target_frame), "Time"] - eic_rt[apex_eic]) > allowed_rt_shift) {
#         return(NULL)
#       }
#     }
#     im_data <- SetData(im_data, "eic", apex_eim + c(-1, 1) * mobility_intgration_range)
#     im_data <- SetInfo(im_data, apex_eic, apex_eim, unname(dr[intensity_col]), frame_integration_range)
#     im_data@peak_quality <- peak_quality
#     if (remove_profile) {
#       im_data <- DropProfile(im_data)
#     }
#     return(im_data)
#   })
#   names(res) <- rownames(info)
# 
#   return(res)
# }
.extract_im_data <- function(
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
    allowed_mobility_shift = 0.015,
    allowed_rt_shift = 10,
    ...
) {
  # browser()
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
    
    # extract eic-eim profile around the target frame and eim apex
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
  names(res) <- rownames(info)
  return(res)
}
.get_peaks <- function(
  spectra_file,
  im_data_file,
  smp_idx,
  ...
) {
  im_data <- readRDS(im_data_file)
  im_data <- im_data[!sapply(im_data, is.null)]
  peaks <- do.call(rbind, lapply(im_data, GetPeakInfo))
  # peak_list <- lapply(im_data, slot, "info")
  #
  # peaks <- do.call(rbind, peak_list)
  # info <- readRDS(spectra_file)@info
  # colnames(info)[match(c("mz", "rt"), colnames(info))] <- c("target_mz", "target_rt")
  # peaks <- cbind(info[match(rownames(peaks), rownames(info)),], peaks)
  peaks$ccs <- .mobility2ccs(peaks$mobility, peaks$mz)
  peaks$spec_idx <- rownames(peaks)
  peaks$smp_idx <- smp_idx
  return(peaks)
}

.dereplicate_peaks <- function(
  peaks,
  spectra_data,
  precursor_bins,
  order_column = "area",
  mz_tol = 20,
  rt_tol = 20,
  mobility_tol = 0.015,
  match_msms = FALSE,
  msms_cutoff = 0.8,
  method_match = c("direct", "bootstrap"),
  method_score = c("dp", "msdial", "combined"),
  weight_intensity = 1,
  weight_mz = 0,
  ms2_range = NULL,
  int_threshold_abs = 0,
  int_threshold_rel = 0,
  include_precursor = TRUE,
  mz_tol_ms1 = 20,
  mz_tol_ms2 = 35,
  res_define_at = 200,
  ...
) {

  peaks <- peaks[order(peaks[, order_column], decreasing = TRUE), , drop = FALSE]
  precursor_bins <- precursor_bins[rownames(peaks)]
  if (match_msms) {
    match_param <- SpectraTools::MatchParam(ppm = mz_tol_ms2,
                                            cutoff = 0,
                                            methodMatch = method_match,
                                            methodScore = method_score,
                                            includePrecursor = include_precursor,
                                            ppmPrecursorFilter = mz_tol_ms1,
                                            ms2range = ms2_range,
                                            intensityExpNormed = FALSE,
                                            intensityLibNormed = FALSE,
                                            tuneLibSpectra = TRUE,
                                            thrIntensityAbs = int_threshold_abs,
                                            thrIntensityRel = int_threshold_rel,
                                            weightIntensity = weight_intensity,
                                            weightMZ = weight_mz,
                                            useMS1ResDefine = TRUE,
                                            resDefineAt = res_define_at)
  }

  if (mz_tol > 1) {
    peaks$mz_tol <- sapply(peaks$mz, function(mz) .ppm2dalton(mz, mz_tol, res_define_at))
  } else {
    peaks$mz_tol <- mz_tol
  }
  while (nrow(peaks) > 0) {
    mz <- peaks[1, "mz"]
    rt <- peaks[1, "rt"]
    mobility <- peaks[1, "mobility"]
    mz_tol <- peaks[1, "mz_tol"]
    is_match <- abs(peaks$mz - mz) <= mz_tol &
      abs(peaks$rt - rt) <= rt_tol &
      abs(peaks$mobility - mobility) <= mobility_tol
    idx_remove <- which(is_match)
    target_cads <- rownames(peaks)[idx_remove]
    if (match_msms && length(idx_remove) > 1) {
      msms_score <- c(1, .match_spectra(spectra_data[precursor_bins[[target_cads[1]]]],
                                        spectra_data[unlist(precursor_bins[target_cads[-1]])],
                                        match_param))
      is_keep <- msms_score >= msms_cutoff
      idx_remove <- idx_remove[is_keep]
      target_cads <- target_cads[is_keep]
    }
    # peak_remove <- peaks[idx_remove, , drop = FALSE]
    peaks <- peaks[-idx_remove, , drop = FALSE]
    if (length(idx_remove) > 1) {
      precursor_bins[[target_cads[1]]] <- unname(do.call(c, precursor_bins[target_cads]))
      precursor_bins <- precursor_bins[-match(target_cads[-1], names(precursor_bins))]
    }
  }
  return(precursor_bins)
}

.correct_rt_landmarks <- function(
  peak_ref,
  peak_query,
  spectra_data,
  rt_tol = 30,
  cutoff = 0.8,
  method_match = c("direct", "bootstrap"),
  method_score = c("dp", "msdial", "combined"),
  weight_intensity = 1,
  weight_mz = 0,
  int_threshold_abs = 0,
  int_threshold_rel = 0,
  include_precursor = TRUE,
  mz_tol_ms1 = 20,
  mz_tol_ms2 = 35,
  ms2_range = NULL,
  rt_range = c(0, 720),
  res_define_at = 200,
  data_file = NULL,
  ...
) {
  match_param <- SpectraTools::MatchParam(ppm = mz_tol_ms2,
                                          cutoff = 0,
                                          methodMatch = method_match,
                                          methodScore = method_score,
                                          includePrecursor = include_precursor,
                                          ppmPrecursorFilter = mz_tol_ms1,
                                          ms2range = ms2_range,
                                          intensityExpNormed = FALSE,
                                          intensityLibNormed = FALSE,
                                          tuneLibSpectra = TRUE,
                                          thrIntensityAbs = int_threshold_abs,
                                          thrIntensityRel = int_threshold_rel,
                                          weightIntensity = weight_intensity,
                                          weightMZ = weight_mz,
                                          useMS1ResDefine = TRUE,
                                          resDefineAt = res_define_at)
  # peak_ref <- peaks[peaks$smp_idx == ref_index, , drop = FALSE]
  # peak_query <- peaks[peaks$smp_idx == query_index, , drop = FALSE]
  peak_ref$rt_correct <- peak_ref$rt

  peak_cads <- .find_peak_cads(peak_ref, peak_query, rt_tol)
  peak_cads2 <- .match_peak_cads(peak_cads, spectra_data, match_param, cutoff)
  land_marks <- .finalize_peak_cads(peak_cads2)

  dt_rt <- data.frame("ref" = c(rt_range[1], peak_ref[land_marks[, "ref"], "rt_correct"], rt_range[2]),
                      "query" = c(rt_range[1], peak_query[land_marks[, "query"], "rt"], rt_range[2]))
  dt_rt <- dt_rt[order(dt_rt$ref), , drop = FALSE]

  model <- loess(ref ~ query, data = dt_rt, span = 0.1, degree = 1L)
  peak_query$rt_correct <- predict(model, peak_query$rt)
  peak_query$rt_correct[peak_query$rt_correct < rt_range[1]] <- rt_range[1]
  peak_query$rt_correct[peak_query$rt_correct > rt_range[2]] <- rt_range[2]

  res <- list("rt" = peak_query$rt_correct, "land_marks" = land_marks, 'peak_cads' = peak_cads)
  if (!is.null(data_file)) {
    saveRDS(res, file = data_file, version = 2)
    return(res)
  } else {
    return(res)
  }
}

.align_peaks_raw <- function(
  ref_index,
  query_index,
  peaks,
  query_land_marks,
  query_peak_cads,
  rt_tol = 15,
  data_file,
  ...
) {
  pooled <- peaks[peaks$smp_idx == ref_index, , drop = FALSE]

  peak_groups <- apply(peaks[peaks$smp_idx %in% c(ref_index, query_index), , drop = FALSE], 1, function(dr) {
    x <- dr["peak_idx"]
    names(x) <- dr["smp_idx"]
    list(x)
  })
  peak_groups <- lapply(peak_groups, unlist)
  num_smps <- length(unique(peaks$smp_idx))
  for (idx in seq_along(query_index)) {
    smp_idx <- query_index[idx]
    land_marks <- query_land_marks[[idx]]
    peak_cads <- query_peak_cads[[idx]]
    peak_query <- peaks[peaks$smp_idx == smp_idx, , drop = FALSE]
    nms_pool <- rownames(pooled)
    nms_query <- rownames(peak_query)

    peak_cads <- peak_cads[!names(peak_cads) %in% land_marks[, 'ref']]
    peak_cads <- lapply(peak_cads, function(cad_info) {
      cad_idx <- setdiff(cad_info$query, land_marks[, 'query'])
      return(cad_info[cad_idx, , drop = FALSE])
    })
    peak_cads <- peak_cads[sapply(peak_cads, nrow) > 0]
    matched_peaks <- .finalize_peak_cads(peak_cads, rt_tol)

    overlaps <- rbind(land_marks, matched_peaks)
    for (nr_overlap in seq(nrow(overlaps))) {
      peak_groups[[overlaps[nr_overlap, "ref"]]] <- c(peak_groups[[overlaps[nr_overlap, "ref"]]],
                                                      peak_groups[[overlaps[nr_overlap, "query"]]])
    }
    peak_groups <- peak_groups[-match(overlaps[, "query"], names(peak_groups))]
    pooled <- rbind(pooled, peak_query[-match(overlaps[, "query"], nms_query), , drop = FALSE])
  }

  peak_groups <- peak_groups[match(pooled$peak_idx, names(peak_groups))]

  res <- list("pooled" = pooled, "peak_groups" = peak_groups)

  if (!is.null(data_file)) {
    saveRDS(res, file = data_file, version = 2)
    return(res)
  } else {
    return(res)
  }
}

.align_peaks_pool <- function(
  ref_index,
  ref_names,
  peak_align_data,
  rt_tol,
  data_file,
  ...
) {
  ref_groups <- sapply(ref_names, function(ref_name) {
    unlist(lapply(peak_align_data, function(data) {
      data$peak_groups[[ref_name]]
    }))
  })

  if (length(peak_align_data) > 1) {
    ref_peaks <- peak_align_data[[1]]$pooled[peak_align_data[[1]]$pooled$smp_idx == ref_index, , drop = FALSE]
    pooled <- peak_align_data[[1]]$pooled[peak_align_data[[1]]$pooled$smp_idx != ref_index, , drop = FALSE]
    peak_groups <- peak_align_data[[1]]$peak_groups
    peak_groups <- peak_groups[!names(peak_groups) %in% ref_names]
    for (idx in seq(2, length(peak_align_data))) {
      peaks_query <- peak_align_data[[idx]]$pooled[peak_align_data[[idx]]$pooled$smp_idx != ref_index, , drop = FALSE]
      peak_groups_query <- peak_align_data[[idx]]$peak_groups
      peak_groups_query <- peak_groups_query[!names(peak_groups_query) %in% ref_names]
      peak_groups <- c(peak_groups, peak_groups_query)
      peak_cads <- .find_peak_cads(pooled, peaks_query, rt_tol = rt_tol)
      matched_peaks <- .finalize_peak_cads(peak_cads, rt_tol)
      for (nr_matched in seq(nrow(matched_peaks))) {
        peak_groups[[matched_peaks[nr_matched, "ref"]]] <- c(peak_groups[[matched_peaks[nr_matched, "ref"]]],
                                                             peak_groups[[matched_peaks[nr_matched, "query"]]])
      }
      peak_groups <- peak_groups[-match(matched_peaks[, "query"], names(peak_groups))]
      pooled <- rbind(pooled, peaks_query[-match(matched_peaks[, "query"], rownames(peaks_query)), , drop = FALSE])
    }
    res <- list("peak_groups" = c(ref_groups, peak_groups),
                "peaks" = rbind(ref_peaks, pooled))
  } else {
    res <- list("peak_groups" = peak_align_data[[1]]$peak_groups,
                "peaks" = peak_align_data[[1]]$pooled)
  }

  return(res)
}

.align_peaks <- function(
  ref_index,
  peaks,
  spectra_files,
  mz_tol = 10,
  strict_rt_constrains = TRUE,
  rt_tol_landmark = 30,
  rt_tol_match = 15,
  ccs_tol = 1.5,
  cutoff = 0.8,
  min_num_fragments = 3,
  method_match = c("direct", "bootstrap"),
  method_score = c("dp", "msdial", "combined"),
  weight_intensity = 1,
  weight_mz = 0,
  int_threshold_abs = 0,
  int_threshold_rel = 0,
  include_precursor = TRUE,
  mz_tol_ms1 = 20,
  mz_tol_ms2 = 35,
  ms1_range = NULL,
  ms2_range = NULL,
  rt_range = c(0, 720),
  res_define_at = 200,
  ...
) {
  if (strict_rt_constrains) {
    peaks <- peaks[peaks$rt >= rt_range[1] & peaks$rt <= rt_range[2], , drop = FALSE]
  }
  peaks$ccs_tol <- .percentage2value(peaks$ccs, ccs_tol)
  if (mz_tol > 1) {
    peaks$mz_tol <- sapply(peaks$mz, function(mz) .ppm2dalton(mz, mz_tol, res_define_at))
  } else {
    peaks$mz_tol <- mz_tol
  }
  match_param <- SpectraTools::MatchParam(ppm = mz_tol_ms2,
                                          cutoff = 0,
                                          methodMatch = method_match,
                                          methodScore = method_score,
                                          includePrecursor = include_precursor,
                                          ppmPrecursorFilter = mz_tol_ms1,
                                          ms2range = ms2_range,
                                          intensityExpNormed = FALSE,
                                          intensityLibNormed = FALSE,
                                          tuneLibSpectra = TRUE,
                                          thrIntensityAbs = int_threshold_abs,
                                          thrIntensityRel = int_threshold_rel,
                                          weightIntensity = weight_intensity,
                                          weightMZ = weight_mz,
                                          useMS1ResDefine = TRUE,
                                          resDefineAt = res_define_at)
  pooled <- peaks[peaks$smp_idx == ref_index, , drop = FALSE]
  pooled$rt_align <- pooled$rt

  spectra_all <- .get_peak_spectra(peaks, spectra_files, min_num_fragments, update_info = TRUE)

  # peak_groups <- as.list(data.frame(apply(peaks, 1, function(dr) {
  #   dr[c("peak_idx", "smp_idx")]
  # }), check.names = FALSE))
  peak_groups <- apply(peaks, 1, function(dr) {
    x <- dr["peak_idx"]
    names(x) <- dr["smp_idx"]
    list(x)
  })
  peak_groups <- lapply(peak_groups, unlist)
  num_smps <- length(unique(peaks$smp_idx))
  rev_models <- vector("list", num_smps)

  for (smp_idx in seq(num_smps)[-ref_index]) {
    peak_query <- peaks[peaks$smp_idx == smp_idx, , drop = FALSE]
    nms_pool <- rownames(pooled)
    nms_query <- rownames(peak_query)

    peak_cads <- .find_peak_cads(pooled, peak_query, rt_tol_landmark)
    peak_cads2 <- .match_peak_cads(peak_cads, spectra_all, match_param, cutoff)
    land_marks <- .finalize_peak_cads(peak_cads2)

    dt_rt <- data.frame("ref" = c(rt_range[1], pooled[land_marks[, "ref"], "rt_align"], rt_range[2]),
                        "query" = c(rt_range[1], peak_query[land_marks[, "query"], "rt"], rt_range[2]))
    dt_rt <- dt_rt[order(dt_rt$ref), , drop = FALSE]
    model <- loess(ref ~ query, data = dt_rt, span = 0.1, degree = 1L)
    peak_query$rt_align <- predict(model, peak_query$rt)
    peak_query$rt_align[peak_query$rt_align < rt_range[1]] <- rt_range[1]
    peak_query$rt_align[peak_query$rt_align > rt_range[2]] <- rt_range[2]

    rev_models[[smp_idx]] <- loess(query ~ ref, data = dt_rt[order(dt_rt$query), , drop = FALSE], span = 0.1, degree = 1L)

    peak_cads <- peak_cads[!names(peak_cads) %in% land_marks[, 'ref']]
    peak_cads <- lapply(peak_cads, function(cad_info) {
      cad_idx <- setdiff(cad_info$query, land_marks[, 'query'])
      return(cad_info[cad_idx, , drop = FALSE])
    })
    peak_cads <- peak_cads[sapply(peak_cads, nrow) > 0]
    matched_peaks <- .finalize_peak_cads(peak_cads, rt_tol_match)

    overlaps <- rbind(land_marks, matched_peaks)
    for (nr_overlap in seq(nrow(overlaps))) {
      peak_groups[[overlaps[nr_overlap, "ref"]]] <- c(peak_groups[[overlaps[nr_overlap, "ref"]]],
                                                      peak_groups[[overlaps[nr_overlap, "query"]]])
    }
    peak_groups <- peak_groups[-match(overlaps[, "query"], names(peak_groups))]
    pooled <- rbind(pooled, peak_query[-match(overlaps[, "query"], nms_query), , drop = FALSE])
  }
  peak_groups <- peak_groups[match(pooled$peak_idx, names(peak_groups))] ### lmd 20220617
  return(list("features" = pooled, "peak_groups" = peak_groups, "rev_models" = rev_models))
}

.group_peaks_density <- function(
  peaks,
  bw = 30,
  mz_bin_size = 0.01,
  mobility_bin_size = 0.015,
  max_features = 50,
  res_define_at = 200,
  plot_density = FALSE,
  plot_dir = '.',
  ...
) {
  # modified from xcms
  mz_start <- min(peaks$mz)
  mz_end <- max(peaks$mz) + ifelse(mz_bin_size > 1, .ppm2dalton(max(peaks$mz), mz_bin_size, res_define_at), mz_bin_size)
  mz_bin_size_half <- mz_bin_size / 2

  if (mz_bin_size > 1) {
    mass <- vector("numeric", length = (mz_end - mz_start) / .ppm2dalton(res_define_at, mz_bin_size_half, res_define_at))
    mass_1 <- mass[1] <- mz_start
    mass_idx <- 1
    while(mass_1 <= mz_end) {
      mass_1 <- mass_1 + .ppm2dalton(mass_1, mz_bin_size_half, res_define_at)
      mass_idx <- mass_idx + 1
      mass[mass_idx] <- mass_1
    }
    mass <- mass[mass > 0]
  } else {
    mass <- seq(min(peaks$mz), max(peaks$mz) + mz_bin_size, by = mz_bin_size/2)
  }
  mobility <- seq(min(peaks$mobility), max(peaks$mobility) + mobility_bin_size, by = mobility_bin_size/2)

  peaks <- peaks[order(peaks[, 'mobility']), , drop = FALSE]

  mobility_pos <- find_greater_equal_than(peaks$mobility, mobility)
  res_list <- vector('list', length(mobility) - 2)
  names(res_list) <- paste0('#', seq_along(res_list))
  res_mass_list_tmp <- vector('list', length(mass) - 2)
  names(res_mass_list_tmp) <- paste0('#', seq_along(res_mass_list_tmp))
  mz_idxs <- seq(length(res_mass_list_tmp))

  dens_from = min(peaks$rt_align) - 3 * bw
  dens_to = max(peaks$rt_align) + 3 * bw
  dens_num = max(512, 2 * 2^(ceiling(log2(diff(range(peaks$rt_align)) / (bw / 2)))))

  if (plot_density) {
    pdf(file.path(plot_dir, 'group_density.pdf'), width = 8, height = 4)
  }

  for (idx_mobility in seq(length(mobility) - 2)) {
    idx_mobility_start <- mobility_pos[idx_mobility]
    idx_mobility_end <- mobility_pos[idx_mobility + 2] - 1
    if (idx_mobility_end - idx_mobility_start < 0) {
      next
    }
    peaks_mz <- peaks[idx_mobility_start:idx_mobility_end, ]
    peaks_mz <- peaks_mz[order(peaks_mz$mz), , drop = FALSE]
    mass_pos <- find_greater_equal_than(peaks_mz$mz, mass)
    res_mass_list <- res_mass_list_tmp
    for (idx_mz in mz_idxs) {
      idx_mz_start <- mass_pos[idx_mz]
      idx_mz_end <- mass_pos[idx_mz + 2] - 1
      if (idx_mz_end - idx_mz_start < 0) {
        next
      }
      # res_mass_list[[idx_mz]] <- peaks_mz[idx_mz_start:idx_mz_end, , drop=FALSE]
      res_mass_list[[idx_mz]] <- .group_density(peaks_mz[idx_mz_start:idx_mz_end, ],
                                                bw = bw,
                                                dens_from = dens_from,
                                                dens_to = dens_to,
                                                dens_num = dens_num,
                                                max_features = max_features,
                                                plot_density = plot_density)
    }
    res_mass_list <- res_mass_list[!sapply(res_mass_list, is.null)]
    res_list[[idx_mobility]] <- res_mass_list
  }
  if (plot_density) {
    dev.off()
  }
  res_list <- res_list[!sapply(res_list, is.null)]
  res <- do.call(rbind, unlist(res_list, recursive = FALSE))

  if (nrow(res)) {
    ## Remove groups that overlap with more "well-behaved" groups
    uorder <- order(res[, "npeaks"])

    uindex <- rect_unique(
      as.matrix(res[, c("mzmin", "mzmax", "mobilitymin", "mobilitymax", "rtmin", "rtmax"),
                      drop = FALSE]), uorder-1, c(0, 0, 0))
    res <- res[which(uindex == 1), , drop = FALSE]
    # peaks2 <- peaks[!peaks$peak_idx %in% do.call(c, res_u[, 'peakidx']), , drop = FALSE]
  }

  group_index <- res[, 'peakidx']
  res <- res[, 1:(ncol(res)-1), drop = FALSE]
  rownames(res) <- names(group_index) <- .gen_indexes(res)
  smp_idx <- peaks$smp_idx
  names(smp_idx) <- rownames(peaks)
  # group_index <- lapply(group_index, function(peak_idx) {
  #   names(peak_idx) <- smp_idx[peak_idx]
  #   peak_idx
  # })
  ext_data <- {
    ft_idx <- unname(do.call(c, mapply(rep, names(group_index), sapply(group_index, length))))
    apply(peaks[do.call(c, group_index), c("ccs", "target_intensity", "height", "height_fit", "area")], 2, function(dc) tapply(dc, ft_idx, median))
  }
  res <- cbind(res, ext_data[rownames(res), , drop = FALSE])
  return(list('features' = res, 'peak_groups' = group_index))
}

.match_between_runs <- function(
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
  im_data <- .extract_im_data(
    info = info,
    tims_data_file = tims_data_file,
    precursor_bin_file = NULL,
    order_column = "area",
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

.get_matched_peaks <- function(
  smp_index,
  features,
  match_between_run_file,
  rt_model = NULL,
  data_file = NULL,
  ...
) {
  fill_gaps_data <- readRDS(match_between_run_file)
  fill_gaps_data <- fill_gaps_data[!sapply(fill_gaps_data, is.null)]
  # peak_quality <- data.frame(t(sapply(fill_gaps_data, slot, "peak_quality")), stringsAsFactors = FALSE)
  # res <- data.frame(t(sapply(fill_gaps_data, slot, "info")), stringsAsFactors = FALSE)
  res <- do.call(rbind, lapply(fill_gaps_data, GetPeakInfo))
  res$ccs <- .mobility2ccs(res$mobility, res$mz)
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

.finalize_data <- function(
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
  ...
) {
  # browser()
  features <- features[, -match(c("npeaks", "area"), colnames(features)), drop = FALSE]

  sample_classes <- .sample_groups_to_class(sample_groups)

  spectra_all <- .get_peak_spectra(peaks, spectra_files, update_info = TRUE)

  filled_groups <- sapply(unique(filled_peaks$feature_idx), function(feature_idx) {
    rownames(filled_peaks[filled_peaks$feature_idx == feature_idx, "smp_idx", drop = FALSE])
  }, simplify = FALSE)
  all_groups <- sapply(rownames(features), function(nm) {
    c(peak_groups[[nm]], filled_groups[[nm]])
  }, simplify = FALSE)

  cols_keep <- intersect(colnames(peaks), colnames(filled_peaks))
  all_peaks <- rbind(peaks[, cols_keep], filled_peaks[, cols_keep])
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

  spectra_assign_nms <- sapply(rownames(features), function(nm) {
    nms_spec <- intersect(peak_groups[[nm]], rownames(group_sample_peaks[[nm]]))
    idx_max <- which.max(spectra_all[nms_spec]@info[[col_max]])
    nms_spec[idx_max]
  })
  peak_groups <- peak_groups[is_keep]
  features$ccs <- .mobility2ccs(features$mobility, features$mz)
  features$name <- .get_peak_name(features[, c("mz", "rt", "ccs")])
  features <- features[, c("name", colnames(features)[1:(ncol(features) - 1)])]

  smp_idxs <- seq(nrow(sample_groups))
  all_peaks <- rbind(peaks[, cols_keep], filled_peaks[, cols_keep])
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
  spec <- SpectraTools::UpdateNames(spec, rownames(features))
  return(list("features" = cbind(features, quant_data), "spectra" = spec, "peak_groups" = peak_groups))
}

.fill_peaks <- function(
  info,
  tims_data_file,
  mz_tol = 20,
  frame_integration_range = 5,
  mobility_intgration_range = 0.015,
  res_define_at = 200,
  data_file = NULL,
  ...
) {
  # browser()
  query_data <- readRDS(tims_data_file)
  num_frames <- nrow(query_data$ms1_frame_info)

  if (mz_tol > 1) {
    info$mz_tol <- sapply(info$mz, function(mz) .ppm2dalton(mz, mz_tol, res_define_at))
  } else {
    info$mz_tol <- mz_tol
  }

  mobility_col <- ifelse("k0" %in% colnames(info), "k0", "mobility")

  info$target_frame <- sapply(info[, "rt"], function(rt) {
    query_data$ms1_frame_info$Id[which.min(abs(query_data$ms1_frame_info$Time - rt))]
  })

  res <- apply(info[, c("mz", mobility_col, "mz_tol", "target_frame")], 1, function(dr) {
    mz <- dr["mz"]
    mobility <- dr[mobility_col]
    mz_tol <- dr["mz_tol"]
    target_frame <- dr["target_frame"]

    target_idx <- which(query_data$ms1_frame_info$Id == target_frame)
    frame_start <- max(0, target_idx - frame_integration_range)
    extract_frames <- query_data$ms1_frame_info$Id[frame_start:min(num_frames, target_idx + frame_integration_range)]
    num_extract_frames <- length(extract_frames)

    # get eim from the target frame
    eim <- .get_eims2(query_data$all_frames[as.character(extract_frames)], query_data$all_mobility,
                      mz, mz_tol, mobility, mobility_intgration_range)
    sum(eim[, seq(length(extract_frames))])
  })
  if (!is.null(data_file)) {
    saveRDS(res, file = data_file, version = 2)
  } else {
    return(res)
  }
}
