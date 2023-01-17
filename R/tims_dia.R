setGeneric("QueryTimsData_DIA",
           function(object, param, ...)
             standardGeneric("QueryTimsData_DIA"))

setGeneric("ExtractIMMSMS",
           function(object, param, ...)
             standardGeneric("ExtractIMMSMS"))

setGeneric("FinalizeFeatures_DIA",
           function(object, param, ...)
             standardGeneric("FinalizeFeatures_DIA"))

####### query ms2 #############
#' @export
setMethod(
  "QueryTimsData_DIA",
  signature = c("TimsData", "QueryTimsDataParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    
    message("Querying tims data...")
    
    param_list <- as.list(param)
    par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)
    
    files <- .tmp_files(object@files, object@experiment@tmp_dir, 'tims_data_ms2')
    names(files) <- object@files
    
    tims_data_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(object@files[idxs], function(data_file) {
        c("tims_data_file" = data_file,
          param_list)
      })
      .parallel_parser(".query_tims_data_ms2", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))
    names(tims_data_files) <- object@files
    object@tmp_data_files$tims_data_files_ms2 <- tims_data_files
    
    setwd(wd0)
    return(object)
  })

.query_tims_data_ms2 <- function(
    tims_data_file,
    opentims_thread = 2L,
    data_file = NULL,
    ...
) {
  opentimsr::opentims_set_threads(opentims_thread)
  
  all_columns <- .set_opentims()
  D <- opentimsr::OpenTIMS(tims_data_file)
  
  frame_id <- opentimsr::table2df(D, 'Frames')$Frames
  
  if('PasefFrameMsMsInfo' %in% opentimsr::tables_names(D)) {
    msms_info <- opentimsr::table2df(D, 'PasefFrameMsMsInfo')$PasefFrameMsMsInfo
  } else {
    msms_info <- merge(opentimsr::table2df(D, 'DiaFrameMsMsInfo')$DiaFrameMsMsInfo,
                       opentimsr::table2df(D, 'DiaFrameMsMsWindows')$DiaFrameMsMsWindows,
                       by = 'WindowGroup')
  }
  
  frame_info <- frame_id[frame_id$MsMsType != 0, c("Id", "Time")]
  all_frames <- lapply(frame_info$Id, function(frame) {
    # opentimsr::query(D, frames = frame, columns = c("scan", "mz", "intensity"))
    opentimsr::query(D, frames = frame, columns = all_columns)
  })
  names(all_frames) <- frame_info$Id
  all_mobility <- .query_scan_mobility(D, frame_info$Id)
  res <-   list("frame_info" = frame_info,
                "msms_info" = msms_info,
                "all_frames" = all_frames,
                "all_mobility" = all_mobility)
  
  opentimsr::CloseTIMS(D)
  if (!is.null(data_file)) {
    saveRDS(res, file = data_file, version = 2)
  } else {
    return(res)
  }
}

####### extract ms2 #############
#' @export
setMethod(
  "ExtractIMMSMS",
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
        c(list("precursor_info" = object@filled_peaks[object@filled_peaks$smp_idx == smp_idx, , drop = FALSE],
               "tims_data_file" = unname(object@tmp_data_files$tims_data_files_ms2[object@files[smp_idx]])),
          param_list)
      })
      .parallel_parser(".extract_im_msms",
                       arg_list, files[idxs], 
                       object@experiment@BPPARAM, 
                       save_in_analysis = TRUE)
    }, simplify = FALSE))
    
    names(spec_files) <- object@files
    object@tmp_data_files$spectra_files <- spec_files
    
    setwd(wd0)
    return(object)
  })


.extract_im_msms <- function(
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
    smooth_method = 'gaussian',
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
    spec_profile <- target_spec[target_spec$scan >= max(scan_range$ScanNumBegin, scan_from) &
                                  target_spec$scan <= min(scan_range$ScanNumEnd, scan_to) &
                                  target_spec$mz <= .ppm2dalton(mz, mz_tol, res_define_at) + mz,
                                c('mz', 'intensity', 'inv_ion_mobility', 'scan')]
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
  names(res) <- row.names(precursor_info)
  res <- res[!sapply(res, is.null)]
  if (!is.null(data_file)) {
    saveRDS(res, file = data_file, version = 2)
  } else {
    return(res)
  }
}


##### finalize feature #######
#' @export
setMethod(
  "FinalizeFeatures_DIA",
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
                       "bpparam" = splus2R::ifelse1(par_analysis, object@experiment@BPPARAM, NULL)),
                  param_list)
    
    finalize_data_file <- .analysis_parser(".finalize_data_dia", arg_list, files)
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


.finalize_data_dia <- function(
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
  
  spectra_assign_nms <- sapply(rownames(features), function(nm) {
    # nms_spec <- intersect(peak_groups[[nm]], rownames(group_sample_peaks[[nm]]))
    nms_spec <- rownames(group_sample_peaks[[nm]])
    idx_max <- which.max(spectra_all[nms_spec]@info[[col_max]])
    # nms_spec[idx_max]
    if(length(idx_max) == 0){
      return(NA)
    }else{
      return(nms_spec[idx_max])
    }
  }, simplify = TRUE)
  peak_groups <- peak_groups[is_keep]
  features$ccs <- .mobility2ccs(features$mobility, features$mz)
  features$name <- .get_peak_name(features[, c("mz", "rt", "ccs")])
  features <- features[, c("name", colnames(features)[1:(ncol(features) - 1)])]
  
  smp_idxs <- seq(nrow(sample_groups))
  # all_peaks <- rbind(peaks[, cols_keep], filled_peaks[, cols_keep])
  all_peaks <- filled_peaks
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
  # spec <- SpectraTools::UpdateNames(spec, rownames(features))
  spec <- SpectraTools::UpdateNames(spec, rownames(features)[which(row.names(features) %in% spec@info$feature_idx)])
  # browser()
  return(list("features" = cbind(features, quant_data), "spectra" = spec, "peak_groups" = peak_groups))
}

