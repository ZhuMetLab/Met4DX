# project: Met4DX
# File name: Utils-TimsData.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/2/22 10:00
# Copyright (c) 2022- ZhuMSLab ALL right reserved

.get_eims2 <- function(frame_list, scan_mobility, mz, mz_tol, mobility, mobility_tol) {
  # browser()
  target_mobility <- scan_mobility[abs(scan_mobility - mobility) <= mobility_tol]
  target_scans <- as.numeric(names(target_mobility))
  temp <- rep(0, length(target_mobility))
  names(temp) <- target_scans
  eims <- lapply(frame_list, function(frame) {
    eim <- .get_frame_data(frame, target_scans, temp, mz, mz_tol)
  })
  mzs <- sapply(eims, function(eim) {
    ifelse(is.null(eim), mz, attributes(eim)$mz)
  })
  if (all(sapply(eims, is.null))) {
    return(NULL)
  }

  res <- data.frame(do.call(cbind, eims), check.names = FALSE)
  res$k0 <- target_mobility
  attr(res, "mz") <- unname(mzs)
  return(res)
}

.get_frame_data <- function(frame, target_scans, temp, mz, mz_tol) {
  mz_diff <- abs(frame$mz - mz)
  is_keep <- (mz_diff <= mz_tol) & (frame$scan %in% target_scans)
  frame <- frame[is_keep, , drop = FALSE]
  if (nrow(frame) == 0 | all(frame$intensity == 0)) {
    attr(temp, "mz") <- mz
    return(temp)
  }

  frame$mz_diff <- mz_diff[is_keep]
  if (nrow(frame) == 0) {
    return(NULL)
  }
  dups <- frame$scan[duplicated(frame$scan)]
  if (length(dups) > 0) {
    idx_remove <- lapply(dups, function(scan) {
      idx <- which(frame$scan == scan)
      return(idx[-which.min(frame$mz_diff[idx])])
    })
    frame <- frame[-unlist(idx_remove), , drop = FALSE]
  }
  temp[as.character(frame$scan)] <- frame$intensity
  attr(temp, "mz") <- .get_weighted_mz(frame$mz, frame$intensity)
  return(temp)
}

.get_eims <- function(D, frames, mz, mz_tol, mobility, mobility_tol, all_columns) {
  dt_frames <- opentimsr::query(D, frames = frames, columns = all_columns)
  scan_mobility <- .get_scan_mobility(dt_frames)
  mz_diff <- abs(dt_frames$mz - mz)
  is_keep <- (mz_diff <= mz_tol) & abs(dt_frames$inv_ion_mobility - mobility) <= mobility_tol
  dt_frames <- dt_frames[is_keep, , drop = FALSE]
  if (nrow(dt_frames) == 0) {
    return(NULL)
  }
  dt_frames$mz_diff <- mz_diff[is_keep]

  scan_rg <- range(dt_frames$scan)
  scans <- seq(scan_rg[1], scan_rg[2])
  tmp <- rep(0, length(scans))
  names(tmp) <- scans

  res <- sapply(frames, function(frame) {
    dt <- dt_frames[dt_frames$frame == frame, , drop = FALSE]
    dups <- dt$scan[duplicated(dt$scan)]
    if (length(dups) > 0) {
      idx_remove <- lapply(dups, function(scan) {
        idx <- which(dt$scan == scan)
        return(idx[-which.min(dt$mz_diff[idx])])
      })
      dt <- dt[-unlist(idx_remove), -1, drop = FALSE]
    }
    tmp[as.character(dt$scan)] <- dt$intensity

    tmp <- c(tmp, 'mz' = .get_weighted_mz(dt$mz, dt$intensity))
  })
  res <- data.frame(res)
  colnames(res) <- frames
  mzs <- as.numeric(res['mz',])
  mzs[is.nan(mzs)] <- mz
  res <- res[1:(nrow(res) - 1), , drop = FALSE]
  res$k0 <- scan_mobility[rownames(res)]
  attr(res, "mz") <- mzs
  return(res)
}

.get_scan_mobility <- function(dt) {
  scan_mobility <- unique(dt[, c("scan", "inv_ion_mobility")])
  scan_mobility <- scan_mobility[order(scan_mobility[, 1]), , drop = FALSE]
  mobility <- scan_mobility[, 2]
  names(mobility) <- scan_mobility[, 1]
  return(mobility)
}

.query_scan_mobility <- function(D, ms1_frames) {
  tmp <- query(D, frames = ms1_frames, columns = c("scan", "inv_ion_mobility"))
  scans <- sort(unique(tmp$scan))
  res <- sapply(scans, function(scan) {
    tmp[match(scan, tmp$scan), "inv_ion_mobility"]
  })
  names(res) <- scans
  return(res)
}

.clust_peaks <- function(distance_matrix, cutoff) {
  clust <- hclust(as.dist(distance_matrix))
  h_cut <- cutree(clust, h = cutoff)
  h_cut_names <- names(h_cut)
  res <- lapply(unique(h_cut), function(v) {
    h_cut_names[h_cut == v]
  })
  return(res)
}

.get_peak_distance <- function(scores, weights) {
  return(sqrt(rowSums((1 - scores)^2 * weights)))
}

.match_spectra <- function(spec_ref, spec_cad, param) {
  msms_match <- SpectraTools::MatchSpectra(spec_ref, spec_cad, param)
  msms_score <- splus2R::ifelse1(is.null(msms_match), 0, msms_match@info$scoreForward)
  return(msms_score)
}

.get_combinations <- function(n) {
  if (n == 1) {
    res <- NA
  } else {
    res <- vector("numeric")
    for (i in seq(n - 1)) {
      res <- rbind(res, cbind(i, seq(i + 1, n)))
    }
    colnames(res) <- NULL
  }
  return(res)
}

.find_apex <- function(data_raw,
                       data_smooth = NULL,
                       min_points = 4,
                       min_intensity = 0,
                       span = 7,
                       ignore_threshold = 0.01,
                       n_skip = 0,
                       ref_index = NULL,
                       find_roi = TRUE,
                       force_max = FALSE) {
  if (is.null(data_smooth)) {
    data_smooth <- data_raw
  } else {
    data_smooth <- round(data_smooth, 10)
  }

  if (find_roi) {
    roi <- .find_roi(data_raw,
                     min_points = min_points,
                     min_intensity = min_intensity,
                     n_skip = n_skip,
                     ref = ref_index)
    if (is.null(roi)) {
      return(NULL)
    }
  } else {
    roi <- matrix(c(1, length(data_raw)), nrow = 1)
  }
  if (nrow(roi) > 0) {
    apex <- which(ggpmisc:::find_peaks(data_smooth, span = span, ignore_threshold = ignore_threshold))
    apex <- apex[apex >= roi[, 1] & apex <= roi[, 2]]
    if (length(apex) > 1 && !is.null(ref_index)) {
      apex <- apex[which.min(abs(apex - ref_index))]
    } else if (length(apex) == 0 && force_max) {
      apex <- which.max(data_smooth)
    }
  } else {
    apex <- NULL
  }
  return(apex)
}

.find_roi <- function(data, min_intensity = 0, min_points = 4, n_skip = 0, ref = NULL) {
  is_roi <- get_continuous_points_above_thr_idx(data, i_start = 0, num = min_points, thr = min_intensity, n_skip = n_skip)
  if (all(is_roi == 0) || length(is_roi) == 0) {
    return(NULL)
  }
  dev_1 <- diff(is_roi)
  scan_min <- which(dev_1 == 1)
  scan_max <- which(dev_1 == -1)
  if (is_roi[1] == 1) {
    scan_min <- c(0, scan_min)
  }
  if (is_roi[length(is_roi)] == 1) {
    scan_max <- c(scan_max, length(is_roi))
  }
  if (length(scan_min) == 0) {
    scan_min <- 1
  } else {
    scan_min <- scan_min + 1
  }
  if (length(scan_max) == 0) {
    scan_max <- length(is_roi)
  }
  roi <- cbind(scan_min, scan_max)
  if (!is.null(ref)) {
    is_match <- (roi[, 1] <= ref) & (roi[, 2] >= ref)
    roi <- splus2R::ifelse1(any(is_match), roi[is_match, , drop = FALSE], NULL)
  }
  return(roi)
}

.get_peak_spectra <- function(peaks, spectra_files, min_num_fragments = NULL, update_info = FALSE) {
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
    spectra_data[[idx]]@spectra
  }))
  idxxxxx <- fastmatch::fmatch(rownames(info), names(spec))
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

.find_peak_cads <- function(peak_ref, peak_query, rt_tol = 30) {
  nms_ref <- rownames(peak_ref)
  nms_query <- rownames(peak_query)
  cads <- apply(peak_query, 1, function(dr) {
    rt <- as.numeric(dr["rt"])
    mz <- as.numeric(dr["mz"])
    ccs <- as.numeric(dr["ccs"])
    rt_diff <- abs(peak_ref$rt_align - rt)
    names(rt_diff) <- nms_ref
    cad_idx <- nms_ref[abs(peak_ref$mz - mz) <= peak_ref$mz_tol &
                         abs(peak_ref$ccs - ccs) <= peak_ref$ccs_tol &
                         rt_diff <= rt_tol]
    if (length(cad_idx) == 0) {
      return(NA)
    }
    res <- data.frame("ref" = cad_idx, "query" = unname(dr["peak_idx"]), "rt_diff" = unname(rt_diff[cad_idx]))
    return(res)
  })
  cad_info <- do.call(rbind, cads[!is.na(cads)])
  matched_ref <- unique(cad_info$ref)
  res <- lapply(matched_ref, function(ref_idx) {
    cad <- cad_info[cad_info$ref == ref_idx, , drop = FALSE]
    rownames(cad) <- cad$query
    return(cad)
  })
  names(res) <- matched_ref
  return(res)
}

.match_peak_cads <- function(peak_cads, spectra_data, match_param, cutoff = 0.8) {
  cad_nms <- intersect(names(peak_cads), names(spectra_data))
  cads <- lapply(cad_nms, function(nm) {
    info <- peak_cads[[nm]]
    spec_ref <- spectra_data[nm]
    spec_cad <- spectra_data[info$query]
    if (length(spec_cad) * length(spec_ref) == 0) {
      return(NA)
    }

    info[names(spec_cad), 'msms_score'] <- .match_spectra(spec_ref, spec_cad, match_param)
    info[names(spec_cad), 'ccs_diff'] <- spec_ref@info$ccs-spec_cad@info$ccs

    spec_cad_nms <- names(spec_cad)
    nm_keep <- spec_cad_nms[info[names(spec_cad), 'msms_score'] >= cutoff]
    cad_idx <- info[nm_keep, "query"]
    rt_diff <- info[nm_keep, "rt_diff"]
    ccs_diff <- abs(info[nm_keep, "ccs_diff"])
    len_cad <- length(cad_idx)
    if (len_cad == 0) {
      return(NA)
    } else if (len_cad > 1) {
      # cad_idx <- cad_idx[which.min(rt_diff)]
      cad_idx <- cad_idx[which.min(ccs_diff)]
    }
    return(info[cad_idx, , drop = FALSE])
  })
  names(cads) <- cad_nms
  cads <- cads[!is.na(cads)]
  return(cads)
}

.finalize_peak_cads <- function(peak_cads, rt_tol = NULL) {
  idx_multiple_match <- which(sapply(peak_cads, nrow) > 1)
  peak_cads[idx_multiple_match] <- lapply(peak_cads[idx_multiple_match], function(info) {
    # return(info[which.min(info$rt_diff), , drop = FALSE])
    return(info[which.min(abs(info$ccs_diff)), , drop = FALSE])
  })
  res <- do.call(rbind, peak_cads)
  if (!is.null(rt_tol)) {
    res <- res[res$rt_diff <= rt_tol, , drop = FALSE]
  }
  res <- t(sapply(unique(res[, "query"]), function(idx) {
    is_keep <- res[, "query"] == idx
    idx_r <- res[is_keep, "ref"]
    if (length(idx_r) > 1) {
      # idx_r <- idx_r[which.min(res[is_keep, "rt_diff"])]
      idx_r <- idx_r[which.min(abs(res[is_keep, "ccs_diff"]))]
    }
    c(idx_r, idx)
  }))
  colnames(res) <- c("ref", "query")
  rownames(res) <- NULL
  return(res)
}

.mobility_to_ccs <- function(mobility, exact_mass,
                             slope = 200.27492809353222,
                             charge = 1) {
  gamma <- (1 / charge) * (exact_mass / (exact_mass + 28.0061))^0.5
  ccs <- slope * (1 / gamma) * mobility
  return(ccs)
}

.mobility2ccs <- function(k0, mz, t = 305, z = 1) {
  factor <- 18509.863216340458
  # em_n2 <- 28.0134
  em_n2 <- 28.00615
  ccs <- factor *
    z *
    k0 *
    sqrt((mz + em_n2) / (t * em_n2 * mz))
  return(ccs)
}

.match_peaks <- function(peak_ref, peak_query,
                         spectra_data = NULL,
                         mz_tol = 10,
                         rt_tol = 30,
                         ccs_tol = 1.5,
                         cutoff = 0.8,
                         res_define_at = 200,
                         match_param = NULL) {
  if (mz_tol > 1) {
    peak_ref$mz_tol <- sapply(peak_ref$mz, function(mz) .ppm2dalton(mz, mz_tol, res_define_at))
  } else {
    peak_ref$mz_tol <- mz_tol
  }

  match_msms <- !is.null(spectra_data)
  nms_ref <- rownames(peak_ref)
  nms_query <- rownames(peak_query)

  if (match_msms) {
    nms_ref <- nms_ref[nms_ref %in% names(spectra_data)]
  }

  cads <- sapply(nms_ref, function(nr) {
    rt_diff <- abs(peak_query$rt - peak_ref[nr, "rt_align"])
    names(rt_diff) <- nms_query
    cad_idx <- nms_query[abs(peak_query$mz - peak_ref[nr, "mz"]) <= peak_ref[nr, "mz_tol"] &
                           rt_diff <= rt_tol &
                           abs(peak_query$ccs - peak_ref[nr, "ccs"]) / peak_ref[nr, "ccs"] * 100 <= ccs_tol]
    len_cad <- length(cad_idx)

    if (match_msms) {
      if (len_cad == 0) {
        return(NA)
      }

      spec_ref <- spectra_data[nr]
      if (length(spec_ref) == 0) {
        return(NA)
      }

      spec_cad <- spectra_data[cad_idx]
      if (length(spec_cad) == 0) {
        return(NA)
      }

      msms_score <- .match_spectra(spec_ref, spec_cad, match_param)

      is_keep <- msms_score > cutoff
      cad_idx <- cad_idx[is_keep]
      len_cad <- length(cad_idx)
    }

    if (len_cad == 0) {
      return(NA)
    }
    if (len_cad > 1) {
      cad_idx <- cad_idx[which.min(rt_diff[cad_idx])]
    }
    return(cad_idx)
  })
  res <- cbind("idx_ref" = nms_ref[!is.na(cads)], "idx_query" = na.omit(cads))
  res <- t(sapply(unique(res[, "idx_query"]), function(idx) {
    idx_r <- res[res[, "idx_query"] == idx, "idx_ref"]
    if (length(idx_r) > 1) {
      idx_r <- idx_r[which.min(abs(peak_ref[idx_r, "rt"] - peak_query[idx, "rt"]))]
    }
    c(idx_r, idx)
  }))
  colnames(res) <- c("ref", "query")
  rownames(res) <- NULL
  return(res)
}

.get_peak_name <- function(info) {
  info <- round(info)
  res <- apply(mapply(paste0, c("M", "T", "C"), info), 1, paste0, collapse = "")
  res <- .make_unique_names(res)
  return(res)
}

.make_unique_names <- function (name) {
  duplicated <- TRUE
  i <- 2
  while (duplicated) {
    idxDuplicated <- which(duplicated(name))
    if (length(idxDuplicated) > 0) {
      if (i == 2) {
        name[idxDuplicated] <- paste0(name[idxDuplicated], "_", i)
      }
      else {
        name[idxDuplicated] <- gsub(paste0("_", i - 1), paste0("_", i), name[idxDuplicated])
      }
      i <- i + 1
    }
    else {
      duplicated <- FALSE
    }
  }
  return(name)
}

.get_smooth_sd <- function(raw, smoothed) {
  sd(raw - smoothed)/diff(range(smoothed))
}

.estimate_noise <- function(x, trim = 0.05, min_points = 20) {
  # modified from xcms
  gz <- which(x > 0)
  if (length(gz) < min_points)
    return(mean(x))

  mean(x[gz], trim=trim)
}

.trimm <- function(x, trim=c(0.05,0.95)) {
  a <- sort(x[x > 0])
  num <- length(a)
  quant <- round((num * trim[1]) + 1):round(num * trim[2])
  return(a[quant])
}

.estimate_local_noise <- function(d, td, ftd, noiserange, num_points, threshold, num_continuous = 4) {
  # from xcms
  if (length(d) < num_points) {

    ## noiserange[2] is full d-range
    drange <- which(td %in% ftd)
    n1 <- d[-drange] ## region outside the detected ROI (wide)
    n1.cp <- get_continuous_points_above_thr_idx(n1, thr = threshold, num = num_continuous, i_start = 0, n_skip = 0) ## continousPtsAboveThreshold (probably peak) are subtracted from data for local noise estimation
    n1 <- n1[!n1.cp]
    if (length(n1) > 1)  {
      baseline1 <- mean(n1)
      sdnoise1 <- sd(n1)
    } else
      baseline1 <- sdnoise1 <- 1

    ## noiserange[1]
    d1 <- drange[1]
    d2 <- drange[length(drange)]
    nrange2 <- c(max(1,d1 - noiserange[1]) : d1, d2 : min(length(d),d2 + noiserange[1]))
    n2 <- d[nrange2] ## region outside the detected ROI (narrow)
    n2.cp <- get_continuous_points_above_thr_idx(n2, thr=threshold, num=num_continuous, i_start = 0, n_skip = 0) ## continousPtsAboveThreshold (probably peak) are subtracted from data for local noise estimation
    n2 <- n2[!n2.cp]
    if (length(n2) > 1)  {
      baseline2 <- mean(n2)
      sdnoise2 <- sd(n2)
    } else
      baseline2 <- sdnoise2 <- 1

  } else {
    trimmed <- .trimm(d, c(0.05,0.95))
    baseline1 <- baseline2 <- mean(trimmed)
    sdnoise1 <- sdnoise2 <- sd(trimmed)
  }

  c(min(baseline1, baseline2), min(sdnoise1, sdnoise2))
}

.get_peak_quality <- function(
  data,
  apex_index,
  min_intensity = 0,
  min_points = 4,
  n_skip = 0,
  snthreshold = 3,
  skip_invalid_peaks = TRUE
) {
  baseline <- snr <- NA
  roi <- .find_roi(data, min_intensity = 0,
                   min_points = min_points, n_skip = n_skip,
                   ref = apex_index)


  if (is.null(roi)) {
    if (skip_invalid_peaks) {
      return(NULL)
    }
  } else {
    num_data <- length(data)

    # estimate snr of the detected eim peak (modified from xcms)
    noise <- .estimate_noise(data)
    td <- roi[1, 1]:roi[1, 2]
    ftd <- max(1, roi[1, 1] - 10):min(roi[1, 2] + 10, num_data)
    local_noise <- .estimate_local_noise(data, td, ftd,
                                         noiserange = c(12, 42),
                                         num_points = 1000,
                                         threshold = noise,
                                         num_continuous = 4)
    sdnoise <- max(1, local_noise[2])
    baseline <- max(1, min(local_noise[1], noise))
    snr <- (data[apex_index] - baseline) / sdnoise
    if (snr < snthreshold) {
      if (skip_invalid_peaks) {
        return(NULL)
      }
      snr <- NA
    }
  }
  return(list('snr' = snr, 'baseline' = baseline))
}

.group_density <- function(
  peaks,
  bw,
  dens_from,
  dens_to,
  dens_num,
  max_features,
  plot_density = FALSE
) {
  den <- density(peaks$rt_align,
                 bw=bw,
                 from = dens_from,
                 to = dens_to,
                 n = dens_num)
  maxden <- max(den$y)
  deny <- den$y
  col_nms <- c("mz", "mzmin", "mzmax",
               "mobility", "mobilitymin", "mobilitymax",
               "rt", "rtmin", "rtmax",
               "npeaks")
  res_mat <- matrix(nrow = 0, ncol = length(col_nms),
                    dimnames = list(character(), col_nms))

  pk_nms = rownames(peaks)

  res_idx <- list()
  # ana_set <- c()
  while (deny[maxy <- which.max(deny)] > maxden / 20 && nrow(res_mat) < max_features) {
    grange <- find_dense_min(deny, maxy -  1)
    gidx <- which(peaks[,"rt_align"] >= den$x[grange[1]] &
                    peaks[,"rt_align"] <= den$x[grange[2]] &
                    !is.na(pk_nms))
    deny[grange[1]:grange[2]] <- 0

    if (length(gidx)) {
      res_mat <- rbind(res_mat,
                       c(median(peaks[gidx, "mz"]),
                         range(peaks[gidx, "mz"]),
                         median(peaks[gidx, "mobility"]),
                         range(peaks[gidx, "mobility"]),
                         median(peaks[gidx, "rt_align"]),
                         range(peaks[gidx, "rt_align"]),
                         length(gidx))
      )
      res_idx <- c(res_idx, list(pk_nms[gidx]))
      pk_nms[gidx] <- NA
    }
  }
  res <- as.data.frame(res_mat)
  res$peakidx <- res_idx

  if (plot_density) {
    plot(den, main = paste("mobility: ",
                           round(min(peaks[,"mobility"]), 2), "-",
                           round(max(peaks[,"mobility"]), 2), "| mz: ",
                           round(min(peaks[,"mz"]), 2), "-",
                           round(max(peaks[,"mz"]), 2))
    )

    for (j in seq_len(nrow(res))) {
      points(peaks[res_idx[[j]], "rt_align"], peaks[res_idx[[j]], "area"] / max(peaks[res_idx[[j]], "area"]) * maxden, col = j, pch=20)
      abline(v = res[j, 8:9], lty = "dashed", col = j)
    }
  }
  return(res)
}
