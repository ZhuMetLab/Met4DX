setGeneric("QueryTimsData_MS2",
           function(object, param, ...)
             standardGeneric("QueryTimsData_MS2"))

setGeneric("GenerateMS2",
           function(object, param, ...)
             standardGeneric("GenerateMS2"))


#' @export
setMethod(
  "QueryTimsData_MS2",
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
      .parallel_parser(".query_tims_data_ms2_dda", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))
    names(tims_data_files) <- object@files
    object@tmp_data_files$tims_data_files_ms2 <- tims_data_files
    
    setwd(wd0)
    return(object)
  })


.query_tims_data_ms2_dda <- function(
    tims_data_file,
    opentims_thread = 2L,
    data_file = NULL,
    ...
) {
  opentimsr::opentims_set_threads(opentims_thread)
  
  all_columns <- .set_opentims()
  D <- opentimsr::OpenTIMS(tims_data_file)
  
  frame_id <- opentimsr::table2df(D, 'Frames')$Frames
  # browser()
  if('PasefFrameMsMsInfo' %in% opentimsr::tables_names(D)) {
    msms_info <- opentimsr::table2df(D, 'PasefFrameMsMsInfo')$PasefFrameMsMsInfo
    precursor <- opentimsr::table2df(D, 'Precursors')$Precursors
  } else {
    msms_info <- merge(opentimsr::table2df(D, 'DiaFrameMsMsInfo')$DiaFrameMsMsInfo,
                       opentimsr::table2df(D, 'DiaFrameMsMsWindows')$DiaFrameMsMsWindows,
                       by = 'WindowGroup')
  }
  
  frame_info <- frame_id[,c("Id", "Time")]
  all_frames <- lapply(frame_info$Id, function(frame) {
    # opentimsr::query(D, frames = frame, columns = c("scan", "mz", "intensity"))
    opentimsr::query(D, frames = frame, columns = all_columns)
  })
  names(all_frames) <- frame_info$Id
  all_mobility <- .query_scan_mobility(D, frame_info$Id)
  res <-   list("frame_info" = frame_info,
                "msms_info" = msms_info,
                "precursor" = precursor,
                "all_frames" = all_frames,
                "all_mobility" = all_mobility)
  
  opentimsr::CloseTIMS(D)
  if (!is.null(data_file)) {
    saveRDS(res, file = data_file, version = 2)
  } else {
    return(res)
  }
}



#' @export
setMethod(
  "GenerateMS2",
  signature = c("TimsData", "ReadSpectraParam"),
  function(object, param) {
    wd0 <- getwd()
    setwd(object@experiment@wd)
    message("Generating spectra data...")
    object@.processHistory <- c(object@.processHistory, param)
    param_list <- as.list(param)
    par_idx <- .gen_parallel_indexes(length(object@files), object@experiment@BPPARAM$workers)
    
    files <- .tmp_files(object@files, object@experiment@tmp_dir, 'spectra')
    names(files) <- object@files
    
    spec_files <- do.call(c, apply(par_idx, 1, function(dr) {
      idxs <- seq(dr[1], dr[2])
      cat('processing files: ', dr[1], '-', dr[2], '\n')
      arg_list <- lapply(object@files[idxs], function(data_file) {
        # browser()
        c("data_file" = unname(object@tmp_data_files$tims_data_files_ms2[data_file]),
          "res_define_at" = object@experiment@res_define_at,
          param_list)
      })
      .parallel_parser(".generate_ms2", arg_list, files[idxs], object@experiment@BPPARAM)
    }, simplify = FALSE))
    
    names(spec_files) <- object@files
    object@tmp_data_files$spectra_files <- spec_files
    
    setwd(wd0)
    return(object)
  })


.generate_ms2 <- function(
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
  # include_precursor = TRUE
  # mz_tol = 20
  # denoise = TRUE
  # int_threshold_abs = 30
  # int_threshold_rel = 0.01
  # ms2range = NULL
  # res_define_at = 200
  # data_file <- tims_data@tmp_data_files$tims_data_files_ms2[1]
  # browser()
  res <- readRDS(data_file)
  intensity_from <- match.arg(intensity_from)
  # param <- SpectraTools::ParseSpectraParam(
  #   type = "mgf",
  #   denoise = FALSE,
  #   ms2range = ms2range,
  #   includePrecursor = include_precursor,
  #   ppmPrecursorFilter = mz_tol,
  #   thrIntensityAbs = int_threshold_abs,
  #   thrIntensityRel = int_threshold_rel,
  #   labelKeep = c("PEPMASS", "RTINSECONDS", "TITLE", "RAWSCANS"),
  #   labelName = c("precursor_info", "rt", "info", "raw_scans"),
  #   resDefineAt = res_define_at
  # )
  # spectra <- SpectraTools::ParseSpectra(param, data_file)
  # mz_info_names <- c('mz', 'intensity')
  # info <- t(apply(spectra@info, 1, function(dr) {
  #   .info <- dr["info"]
  #   one_over_k0 <- regmatches(.info, regexpr('(?<=1/K0=)\\d+\\.\\d+(?=,)', .info, perl = TRUE))
  #   cmpd <- regmatches(.info, regexpr('(?<=Cmpd\\s)\\d+(?=,)', .info, perl = TRUE))
  #   scans <- strsplit(dr["raw_scans"], "-|,")[[1]]
  #   mz_info <- head(strsplit(dr["precursor_info"], "\t| ")[[1]], 2)
  #   names(mz_info) <- mz_info_names[seq_along(mz_info)]
  #   c(dr,
  #     mz_info,
  #     "k0" = unname(one_over_k0),
  #     "cmpd" = unname(cmpd),
  #     "ms1_scans" = scans[1],
  #     "ms2_scans" = paste0(scans[-1], collapse = ","))
  # }))
  res_spec <- lapply(unique(res$msms_info$Precursor), function(id){
    ms2_frame <- as.character(res$msms_info[res$msms_info$Precursor == id, ]$Frame)
    
    ScanNumBegin <- as.integer(res$msms_info[res$msms_info$Precursor == id, ]$ScanNumBegin[1])
    ScanNumEnd <- as.integer(res$msms_info[res$msms_info$Precursor == id, ]$ScanNumEnd[1])
    
    all_ion <- lapply(ms2_frame, function(fr){
      temp_frame <- res$all_frames[[fr]]
      temp_frame <- temp_frame[temp_frame$scan > ScanNumBegin & temp_frame$scan < ScanNumEnd, ]
      return(temp_frame)
    })
    
    all_ion <- do.call(rbind, all_ion)
    
    all_ion <- .combine_spec(all_ion)
    return(all_ion)
  })
  
  names(res_spec) <- unique(res$msms_info$Precursor)
  # View(res$msms_info)
  
  info <- lapply(unique(res$msms_info$Precursor), function(cmpd){
    # cmpd <- 17
    idx <- which(res$msms_info$Precursor == cmpd)
    temp_info <- res$msms_info[idx, ]
    ms2_scans <- paste(temp_info$Frame, collapse = ',')
    # mz <- temp_info$IsolationMz[1]
    intensity <- sum(res_spec[[cmpd]]$intensity)
    ms1_scans <- res$precursor$Parent[res$precursor$Id == cmpd]
    mz <- res$precursor$LargestPeakMz[res$precursor$Id == cmpd]
    targeted_scan <- as.integer(res$precursor$ScanNumber[res$precursor$Id == cmpd])
    k0 <- unname(res$all_mobility[as.character(targeted_scan)])
    rt <- res$frame_info$Time[res$frame_info$Id == ms1_scans]
    
    return(data.frame(mz = mz,
                      rt = rt,
                      intensity = intensity,
                      k0 = k0, 
                      cmpd = cmpd,
                      ms1_scans = ms1_scans,
                      ms2_scans = ms1_scans))
  })
  info <- do.call(rbind, info)
  row.names(info) <- info$cmpd
  spectra <- SpectraTools::SpectraData(info = info, 
                                       spectra = res_spec)
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


.combine_spec <- function(all_ion){
  all_ion <- all_ion[all_ion$intensity >= 10, ]
  all_ion <- all_ion[order(all_ion$intensity, decreasing = TRUE), ]
  
  idx <- which(abs(all_ion$tof - all_ion$tof[1]) <= 3)
  # ion_res <- c(ion_res, all_ion[idx, ])
  
  
  ion_res <- list()
  ion_res[[1]] <- all_ion[idx, ]
  all_ion <- all_ion[-idx, ]
  while(nrow(all_ion) > 0){
    idx <- which(abs(all_ion$tof - all_ion$tof[1]) <= 3)
    # ion_res <- c(ion_res, data.frame(all_ion[idx, ]))
    ion_res[[length(ion_res) + 1]] <- all_ion[idx, ]
    all_ion <- all_ion[-idx, ]
  }
  
  ion_res <- lapply(ion_res, .weight_mz_intensity)
  ion_res <- do.call(rbind, ion_res)
  ion_res <- ion_res[order(ion_res$mz, decreasing = FALSE), ]
  return(ion_res)
}


.weight_mz_intensity <- function(ion_res){
  sum_intensity <- sum(ion_res$intensity)
  new_mz <- sum(ion_res$intensity / sum_intensity * ion_res$mz)
  return(data.frame(mz = as.numeric(new_mz), 
                    intensity = sum_intensity))
}

