# project: Met4DX
# File name: Functions-Params.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/1/17 9:50
# Copyright (c) 2022- ZhuMSLab ALL right reserved

## GaussianSmoothParam
#' params for gaussian smoothing
#' @export
GaussianSmoothParam <- function(
  window = getOption("smoother.window"),
  alpha = getOption("smoother.gaussianwindow.alpha")
) {
  new("GaussianSmoothParam",
      window = window,
      alpha = alpha
  )
}

## LOESSSmoothParam
#' params for LOESS smoothing
#' @export
LOESSSmoothParam <- function(
  span = 0.3,
  degree = 1L,
  window = 10L
) {
  new("LOESSSmoothParam",
      span = span,
      degree = degree,
      window = window
  )
}

## FillPeakParam
#' params for filling peaks
#' @export
FillPeakParam <- function(
  mz_tol = 20,
  frame_integration_range = 5,
  mobility_tol = 0.015,
  mobility_intgration_range = 0.015,
  rerun = FALSE
) {
  new("FillPeakParam",
      mz_tol = mz_tol,
      frame_integration_range = frame_integration_range,
      mobility_intgration_range = mobility_intgration_range,
      rerun = rerun
  )
}

## FinalizeFeatureParam
#‘ parameters for finalize data
#' @export
FinalizeFeatureParam <- function(
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
  rerun = FALSE
) {
  quant_method <- match.arg(quant_method)
  col_max <- match.arg(col_max)
  col_quant <- match.arg(col_quant)

  new("FinalizeFeatureParam",
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
      rerun = rerun
  )
}

## MatchBetweenRunParam
#' parameters for gap filling
#' @export
MatchBetweenRunParam <- function(
  mz_tol = 20,
  frame_range = 30,
  frame_integration_range = 5,
  mobility_range = 0.1,
  mobility_intgration_range = 0.015,
  min_points = 4,
  min_intensity = 0,
  n_skip = 0,
  interpolate_method = NULL,
  smooth_method = c("gaussian", "loess"),
  peak_span_eim = 21,
  peak_span_eic = 7,
  smooth_window_eim = 16,
  smooth_window_eic = 8,
  keep_profile = FALSE,
  filter_outlier_peaks = TRUE,
  allowed_mobility_shift = 0.015,
  allowed_rt_shift = 10,
  rerun = FALSE
) {
  if (!is.null(interpolate_method)) {
    interpolate_method <- match.arg(interpolate_method, c("linear", "scans"))
  }

  smooth_method <- match.arg(smooth_method)

  new("MatchBetweenRunParam",
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
      peak_span_eim = peak_span_eim,
      peak_span_eic = peak_span_eic,
      smooth_window_eim = smooth_window_eim,
      smooth_window_eic = smooth_window_eic,
      keep_profile = keep_profile,
      filter_outlier_peaks = filter_outlier_peaks,
      allowed_mobility_shift = allowed_mobility_shift,
      allowed_rt_shift = allowed_rt_shift,
      rerun = rerun
  )

}

## GroupDensityParam
#' parameters for density peak grouping
#' @export
GroupDensityParam <- function(
  bw = 5,
  mz_bin_size = 0.015,
  mobility_bin_size = 0.015,
  max_features = 50,
  plot_density = FALSE,
  rerun = FALSE
) {
  new("GroupDensityParam",
      bw = bw,
      mz_bin_size = mz_bin_size,
      mobility_bin_size = mobility_bin_size,
      max_features = max_features,
      plot_density = plot_density,
      rerun = rerun
  )
}

## DereplicatePeaksParam
#' Parameters for dereplicate MS1 peaks
#' @param mz_tol \code{numeric} m/z tolerance for binning precursor ions
#' @param rt_tol \code{numeric} retention time tolerance for binning precursor ions
#' @param mobility_tol \code{numeric} ion mobility tolerance for binning precursor ions
#' @param rerun \code{logical} if rerun the process
#'
#' @export
DereplicatePeaksParam <- function(
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
  int_threshold_abs = 0,
  int_threshold_rel = 0,
  include_precursor = TRUE,
  mz_tol_ms1 = 20,
  mz_tol_ms2 = 35,
  rerun = FALSE
) {

  method_match <- match.arg(method_match)
  method_score <- match.arg(method_score)
  new("DereplicatePeaksParam",
      order_column = order_column,
      mz_tol = mz_tol,
      rt_tol = rt_tol,
      mobility_tol = mobility_tol,
      match_msms = match_msms,
      msms_cutoff = msms_cutoff,
      method_match = method_match,
      method_score = method_score,
      weight_intensity = weight_intensity,
      weight_mz = weight_mz,
      int_threshold_abs = int_threshold_abs,
      int_threshold_rel = int_threshold_rel,
      include_precursor = include_precursor,
      mz_tol_ms1 = mz_tol_ms1,
      mz_tol_ms2 = mz_tol_ms2,
      rerun = rerun
  )
}

## AlignPeakParamRaw
#' Parameters for aligning MS1 peaks
#' @param mz_tol \code{numeric} m/z tolerance for binning precursor ions
#' @param strict_rt_constrains \code{logical} if remove peaks out of defined rt range in the experiment
#' @param rt_tol_landmark \code{numeric} retention time tolerance for finding landmarks
#' @param rt_tol_match \code{numeric} retention time tolerance for finding matched peaks
#' @param ccs_tol \code{numeric} CCS tolerance for binning precursor ions in percentage
#' @param cutoff \code{numeric} cutoff value for defing a precursor ion bin
#' @param weight_rt \code{numeric} weight for retention time to calculate precursor ion distance
#' @param weight_mobility \code{numeric} weight for ion mobility to calculate precursor ion distance
#' @param weight_msms \code{numeric} weight for MSMS intensity to calculate precursor ion distance
#' @param method_match \code{character} method for matching MSMS spectra
#' @param method_score \code{character} method for scoring MSMS spectra
#' @param weight_mz \code{numeric} weight of m/z when scoring
#' @param weight_intensity \code{numeric} weight of intensity when scoring
#' @param int_threshold_abs \code{numeric} minimal absolute intensity to be kept in MSMS spectra
#' @param int_threshold_rel \code{numeric} minimal relative intensity to be kept in MSMS spectra
#' @param include_precursor \code{logical} if include precursor ion in MSMS spectra when matching
#' @param mz_tol_ms1 \code{numeric} m/z tolerance for defining a precursor ion
#' @param mz_tol_ms2 \code{numeric} m/z tolerance to define a matched MSMS spectra in msms matching
#' @param rerun \code{logical} if rerun the process
#'
#' @export
AlignPeakParamRaw <- function(
  mz_tol = 10,
  strict_rt_constrains = TRUE,
  rt_tol_landmark = 30,
  rt_tol_match = 15,
  ccs_tol = 1.5,
  cutoff = 0.8,
  min_num_fragments = 4,
  method_match = c("direct", "bootstrap"),
  method_score = c("dp", "msdial", "combined"),
  weight_intensity = 1,
  weight_mz = 0,
  int_threshold_abs = 0,
  int_threshold_rel = 0,
  include_precursor = TRUE,
  mz_tol_ms1 = 20,
  mz_tol_ms2 = 35,
  rerun = FALSE
) {

  method_match <- match.arg(method_match)
  method_score <- match.arg(method_score)

  new("AlignPeakParamRaw",
      mz_tol = mz_tol,
      strict_rt_constrains = strict_rt_constrains,
      rt_tol_landmark = rt_tol_landmark,
      rt_tol_match = rt_tol_match,
      ccs_tol = ccs_tol,
      cutoff = cutoff,
      min_num_fragments = min_num_fragments,
      method_match = method_match,
      method_score = method_score,
      weight_intensity = weight_intensity,
      weight_mz = weight_mz,
      int_threshold_abs = int_threshold_abs,
      int_threshold_rel = int_threshold_rel,
      include_precursor = include_precursor,
      mz_tol_ms1 = mz_tol_ms1,
      mz_tol_ms2 = mz_tol_ms2,
      rerun = rerun
  )

}

## AlignPeakParam
#' Parameters for aligning MS1 peaks
#' @param mz_tol \code{numeric} m/z tolerance for binning precursor ions
#' @param strict_rt_constrains \code{logical} if remove peaks out of defined rt range in the experiment
#' @param rt_tol_landmark \code{numeric} retention time tolerance for finding landmarks
#' @param rt_tol_match \code{numeric} retention time tolerance for finding matched peaks
#' @param ccs_tol \code{numeric} CCS tolerance for binning precursor ions in percentage
#' @param cutoff \code{numeric} cutoff value for defing a precursor ion bin
#' @param weight_rt \code{numeric} weight for retention time to calculate precursor ion distance
#' @param weight_mobility \code{numeric} weight for ion mobility to calculate precursor ion distance
#' @param weight_msms \code{numeric} weight for MSMS intensity to calculate precursor ion distance
#' @param method_match \code{character} method for matching MSMS spectra
#' @param method_score \code{character} method for scoring MSMS spectra
#' @param weight_mz \code{numeric} weight of m/z when scoring
#' @param weight_intensity \code{numeric} weight of intensity when scoring
#' @param int_threshold_abs \code{numeric} minimal absolute intensity to be kept in MSMS spectra
#' @param int_threshold_rel \code{numeric} minimal relative intensity to be kept in MSMS spectra
#' @param include_precursor \code{logical} if include precursor ion in MSMS spectra when matching
#' @param mz_tol_ms1 \code{numeric} m/z tolerance for defining a precursor ion
#' @param mz_tol_ms2 \code{numeric} m/z tolerance to define a matched MSMS spectra in msms matching
#' @param rerun \code{logical} if rerun the process
#'
#' @export
AlignPeakParam <- function(
  mz_tol = 10,
  strict_rt_constrains = TRUE,
  rt_tol_landmark = 30,
  rt_tol_match = 15,
  ccs_tol = 1.5,
  cutoff = 0.8,
  min_num_fragments = 4,
  method_match = c("direct", "bootstrap"),
  method_score = c("dp", "msdial", "combined"),
  weight_intensity = 1,
  weight_mz = 0,
  int_threshold_abs = 0,
  int_threshold_rel = 0,
  include_precursor = TRUE,
  mz_tol_ms1 = 20,
  mz_tol_ms2 = 35,
  rerun = FALSE
) {

  method_match <- match.arg(method_match)
  method_score <- match.arg(method_score)

  new("AlignPeakParam",
      mz_tol = mz_tol,
      strict_rt_constrains = strict_rt_constrains,
      rt_tol_landmark = rt_tol_landmark,
      rt_tol_match = rt_tol_match,
      ccs_tol = ccs_tol,
      cutoff = cutoff,
      min_num_fragments = min_num_fragments,
      method_match = method_match,
      method_score = method_score,
      weight_intensity = weight_intensity,
      weight_mz = weight_mz,
      int_threshold_abs = int_threshold_abs,
      int_threshold_rel = int_threshold_rel,
      include_precursor = include_precursor,
      mz_tol_ms1 = mz_tol_ms1,
      mz_tol_ms2 = mz_tol_ms2,
      rerun = rerun
  )

}

## CorrectRTLandmarksParam
#' Parameters for RT correrction of MS1 peaks
#' @param mz_tol \code{numeric} m/z tolerance for binning precursor ions
#' @param strict_rt_constrains \code{logical} if remove peaks out of defined rt range in the experiment
#' @param rt_tol_landmark \code{numeric} retention time tolerance for finding landmarks
#' @param rt_tol_match \code{numeric} retention time tolerance for finding matched peaks
#' @param ccs_tol \code{numeric} CCS tolerance for binning precursor ions in percentage
#' @param cutoff \code{numeric} cutoff value for defing a precursor ion bin
#' @param weight_rt \code{numeric} weight for retention time to calculate precursor ion distance
#' @param weight_mobility \code{numeric} weight for ion mobility to calculate precursor ion distance
#' @param weight_msms \code{numeric} weight for MSMS intensity to calculate precursor ion distance
#' @param method_match \code{character} method for matching MSMS spectra
#' @param method_score \code{character} method for scoring MSMS spectra
#' @param weight_mz \code{numeric} weight of m/z when scoring
#' @param weight_intensity \code{numeric} weight of intensity when scoring
#' @param int_threshold_abs \code{numeric} minimal absolute intensity to be kept in MSMS spectra
#' @param int_threshold_rel \code{numeric} minimal relative intensity to be kept in MSMS spectra
#' @param include_precursor \code{logical} if include precursor ion in MSMS spectra when matching
#' @param mz_tol_ms1 \code{numeric} m/z tolerance for defining a precursor ion
#' @param mz_tol_ms2 \code{numeric} m/z tolerance to define a matched MSMS spectra in msms matching
#' @param rerun \code{logical} if rerun the process
#'
#' @export
CorrectRTLandmarksParam <- function(
  mz_tol = 10,
  strict_rt_constrains = TRUE,
  rt_tol_landmark = 30,
  rt_tol_match = 15,
  ccs_tol = 1.5,
  cutoff = 0.8,
  min_num_fragments = 4,
  method_match = c("direct", "bootstrap"),
  method_score = c("dp", "msdial", "combined"),
  weight_intensity = 1,
  weight_mz = 0,
  int_threshold_abs = 0,
  int_threshold_rel = 0,
  include_precursor = TRUE,
  mz_tol_ms1 = 20,
  mz_tol_ms2 = 35,
  rerun = FALSE
) {

  method_match <- match.arg(method_match)
  method_score <- match.arg(method_score)

  new("CorrectRTLandmarksParam",
      mz_tol = mz_tol,
      strict_rt_constrains = strict_rt_constrains,
      rt_tol_landmark = rt_tol_landmark,
      rt_tol_match = rt_tol_match,
      ccs_tol = ccs_tol,
      cutoff = cutoff,
      min_num_fragments = min_num_fragments,
      method_match = method_match,
      method_score = method_score,
      weight_intensity = weight_intensity,
      weight_mz = weight_mz,
      int_threshold_abs = int_threshold_abs,
      int_threshold_rel = int_threshold_rel,
      include_precursor = include_precursor,
      mz_tol_ms1 = mz_tol_ms1,
      mz_tol_ms2 = mz_tol_ms2,
      rerun = rerun
  )

}

## ExtractIMDataParam
#' parameters for extracting target ion EIC/EIMs
#' @export
ExtractIMDataParam <- function(
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
  smooth_method = "gaussian",
  snthreshold = 3,
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
  new("ExtractIMDataParam",
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

## QueryTimsDataParam
#' parameters for querying TIMS data
#' @export
QueryTimsDataParam <- function(
  opentims_thread = 2,
  rerun = FALSE
) {
  new("QueryTimsDataParam",
      opentims_thread = opentims_thread,
      rerun = rerun
  )
}

## ReadSpectraParam
#' Parameters for reading spectra
#' @param intensity_from \code{character} where the precursor intensities were obtained
#' @param include_precursor \code{logical} if include precursor ion in MSMS spectra
#' @param mz_tol \code{numeric} m/z tolerance for defining a precursor ion
#' @param denoise \code{logical} if denoising MSMS spectra
#' @param int_threshold_abs \code{numeric} minimal absolute intensity to be kept in MSMS spectra
#' @param int_threshold_rel \code{numeric} minimal relative intensity to be kept in MSMS spectra
#' @param rerun \code{logical} if rerun the process
#'
#' @export
ReadSpectraParam <- function(
  intensity_from = c("pepmass", "ms2_intensity"),
  include_precursor = TRUE,
  mz_tol = 20,
  denoise = TRUE,
  int_threshold_abs = 30,
  int_threshold_rel = 0.01,
  rerun = FALSE
) {
  intensity_from <- match.arg(intensity_from)

  new("ReadSpectraParam",
      intensity_from = intensity_from,
      include_precursor = include_precursor,
      mz_tol = mz_tol,
      denoise = denoise,
      int_threshold_abs = int_threshold_abs,
      int_threshold_rel = int_threshold_rel,
      rerun = rerun
  )
}

## BinPrecursorParam
#' Parameters for binning precursor ions
#' @param mz_tol \code{numeric} m/z tolerance for binning precursor ions
#' @param rt_tol \code{numeric} retention time tolerance for binning precursor ions
#' @param mobility_tol \code{numeric} ion mobility tolerance for binning precursor ions
#' @param distance_cutoff \code{numeric} cutoff value for defing a precursor ion bin
#' @param weight_rt \code{numeric} weight for retention time to calculate precursor ion distance
#' @param weight_mobility \code{numeric} weight for ion mobility to calculate precursor ion distance
#' @param weight_msms \code{numeric} weight for MSMS intensity to calculate precursor ion distance
#' @param method_match \code{character} method for matching MSMS spectra
#' @param method_score \code{character} method for scoring MSMS spectra
#' @param weight_mz \code{numeric} weight of m/z when scoring
#' @param weight_intensity \code{numeric} weight of intensity when scoring
#' @param int_threshold_abs \code{numeric} minimal absolute intensity to be kept in MSMS spectra
#' @param int_threshold_rel \code{numeric} minimal relative intensity to be kept in MSMS spectra
#' @param include_precursor \code{logical} if include precursor ion in MSMS spectra when matching
#' @param mz_tol_ms1 \code{numeric} m/z tolerance for defining a precursor ion
#' @param mz_tol_ms2 \code{numeric} m/z tolerance to define a matched MSMS spectra in msms matching
#' @param rerun \code{logical} if rerun the process
#'
#' @export
BinPrecursorParam <- function(
  mz_tol = 20,
  rt_tol = c(10, 20),
  mobility_tol = c(0.015, 0.03),
  distance_cutoff = 1,
  weight_rt = 1,
  weight_mobility = 1,
  weight_msms = 1,
  method_match = c("direct", "bootstrap"),
  method_score = c("dp", "msdial", "combined"),
  weight_intensity = 1,
  weight_mz = 0,
  int_threshold_abs = 0,
  int_threshold_rel = 0,
  include_precursor = TRUE,
  mz_tol_ms1 = 20,
  mz_tol_ms2 = 35,
  rerun = FALSE
) {

  method_match <- match.arg(method_match)
  method_score <- match.arg(method_score)

  new("BinPrecursorParam",
      mz_tol = mz_tol,
      rt_tol = rt_tol,
      mobility_tol = mobility_tol,
      distance_cutoff = distance_cutoff,
      weight_rt = weight_rt,
      weight_mobility = weight_mobility,
      weight_msms = weight_msms,
      method_match = method_match,
      method_score = method_score,
      weight_intensity = weight_intensity,
      weight_mz = weight_mz,
      int_threshold_abs = int_threshold_abs,
      int_threshold_rel = int_threshold_rel,
      include_precursor = include_precursor,
      mz_tol_ms1 = mz_tol_ms1,
      mz_tol_ms2 = mz_tol_ms2,
      rerun = rerun
  )

}

## Functions related to the Param class and sub-classes.##
#' @description Extract all slot values and put them into a list, names being
#'     the slot names. If a slot \code{addParams} exist its content will be
#'     appended to the returned list.
#'
#' @param x A Param class.
#'
#' @author Johannes Rainer, modified by YD Yin
#'
#' @noRd
.param2list <- function(x) {
  ## Get all slot names, skip those matching the provided pattern.
  sNames <- slotNames(x)
  skipSome <- grep(sNames, pattern = "^\\.")
  if (length(skipSome) > 0) {
    sNames <- sNames[-skipSome]
  }
  ## handle a slot called "addParams" differently: this is thougth to contain
  ## ... arguments thus we have to skip this one too.
  if (any(sNames == "addParams")) {
    sNames <- sNames[sNames != "addParams"]
    addP <- x@addParams
  } else {
    addP <- list()
  }
  if (length(sNames) > 0) {
    resL <- vector("list", length(sNames))

    for (i in seq_along(sNames)) {
      resL[i] <- list(slot(x, name = sNames[i]))
    }
    names(resL) <- sNames
    resL <- c(resL, addP)
    return(resL)
  }else {
    return(list())
  }
}





#' Searching spectra Parameters
#'
#' Parameters for searching experimentall related SpectraData from refrence
#' SpectraData
#' @param ppm \code{numeric} ppm tolerance for searching precursor mz
#' @param scoreRT \code{logical} if comparing RT
#' @param toleranceRT \code{numeric (2)} c(Tmin, Tmax), RT range for trapezoidal score.
#' \itemize{
#'    \item Tmin: Topline for trapezoidal score
#'    \item Tmax: Baseline for trapezoidal score
#'  }
#' @param scoreCCS \code{logical} if comparing CCS
#' @param toleranceCCS \code{numeric (2)} c(Cmin, Cmax), CCS range for trapezoidal
#'  score.
#' \itemize{
#'    \item Cmin: Topline for trapezoidal score
#'    \item Cmax: Baseline for trapezoidal score
#'  }
#' @param typeCCS \code{character} CCS tolerance type, either "percentage" or
#'  "absolute"
#' @param adductIncluded \code{character} adduct types to be included in MS1 match
#' @param adductExclude \code{character} adduct types to be excluded in MS1 match
#' @param adductFile \code{character} file path for user provided adduct table
#' @param classIncluded \code{character} compound classes to be included in MS1 match
#' @param classExclude \code{character} compound classes to be excluded in MS1 match
#' @param useMS1ResDefine \code{logical} if use resDefineAt when calculating precursor m/z tolerance
#' @param updateRefMZ \code{logical} if update reference precursor m/z with it's corresponding
#' adduct type
#' @param resDefineAt \code{numeric} m/z (Da) value for resolution definition
#' 
#' @export
SearchParam <- function(ppm = 20,
                        scoreRT = TRUE,
                        toleranceRT = c(30, 60),
                        scoreCCS = TRUE,
                        toleranceCCS = c(2, 4),
                        typeCCS = c("percentage", "absolute"),
                        adductIncluded = NULL,
                        adductExcluded = NULL,
                        adductFile = NULL,
                        classIncluded = NULL,
                        classExcluded = NULL,
                        useMS1ResDefine = TRUE,
                        updateRefMZ = TRUE,
                        resDefineAt = 200) {
  typeCCS = match.arg(typeCCS)
  return(new("SearchParam",
             ppm = ppm,
             scoreRT = scoreRT,
             toleranceRT = toleranceRT,
             scoreCCS = scoreCCS,
             toleranceCCS = toleranceCCS,
             typeCCS = typeCCS,
             adductIncluded = adductIncluded,
             adductExcluded = adductExcluded,
             adductFile = adductFile,
             classIncluded = classIncluded,
             classExcluded = classExcluded,
             useMS1ResDefine = useMS1ResDefine,
             updateRefMZ = updateRefMZ,
             resDefineAt = resDefineAt))
}



#' Matching parameter setup
#'
#' @param ppm \code{numeric} ppm tolerance for MS2 m/z match
#' @param cutoff \code{numeric} threshold for a avaliable match
#' @param methodMatch \code{character} method for matching with librarial spectra
#'  (eithor "direct" or "bootstrapping")
#' @param methodScore \code{character} method for scoring the MSMS match ("dp",
#'  "msdial", "nist", "bs", "combined" are avaliable)
#' @param weightMZ \code{numeric} weight of m/z when scoring
#' @param weightIntensity \code{numeric} weight of intensity when scoring
#' @param includePrecursor \code{logical} if consider the precursor fragments when
#'  matching
#' @param ms2range \code{numeric} mass range setup when acquiring MSMS data
#' @param thrIntensityAbs \code{numeric} absolute intensity threshold to be removed
#' @param thrIntensityRel \code{numeric} relative intensity threshold to be removed
#' @param intensityExpNormed \code{logical} if the spectral intensity is normalized
#'  in experiment spectra
#' @param intensityLibNormed \code{logical} if the spectral intensity is normalized
#'  in library spectra

#' @param tuneLibSpectra \code{logical} if apply thrIntensityAbs or thrIntensityRel
#'  to reference
#' @param useMS1ResDefine \code{logical} if use resDefineAt when calculating precursor m/z tolerance
#' @param resDefineAt \code{numeric} m/z threshold for using ppm tolerance for MS1 or
#'  MS2 match, for smaller m/z values, using the tolerance of ppm.ms1 * res.defineat
#'  to mathing experimental and librarial fragments.
#' @return an \code{ParamMatch} object
#' @rdname MatchParam
#' @export
MatchParam <- function(
    ppm = 35,
    cutoff = 0.8,
    methodMatch = c('direct', 'bootstrap'),
    methodScore = c('dp', 'msdial', 'nist', 'bootstrap', 'combined',
                    'ratio', 'gnps', 'bonanza', 'hybrid'),
    weightMZ = 0,
    weightIntensity = 1,
    includePrecursor = TRUE,
    ppmPrecursorFilter = 20,
    ms2range = NULL,
    thrIntensityAbs = 0,
    thrIntensityRel = 0.01,
    intensityExpNormed = TRUE,
    intensityLibNormed = TRUE,
    tuneLibSpectra = FALSE,
    useMS1ResDefine = TRUE,
    resDefineAt = 200,
    normIntensity = TRUE,
    intensityNormedMethod = c('maximum', 'bonanza', 'gnps', 'hybrid')
) {
  
  methodMatch <- match.arg(methodMatch)
  methodScore <- match.arg(methodScore)
  intensityNormedMethod <- match.arg(intensityNormedMethod)
  
  return(new("MatchParam",
             ppm = ppm,
             cutoff = cutoff,
             methodMatch = methodMatch,
             methodScore = methodScore,
             weightMZ = weightMZ,
             weightIntensity = weightIntensity,
             includePrecursor = includePrecursor,
             ppmPrecursorFilter = ppmPrecursorFilter,
             ms2range = ms2range,
             thrIntensityAbs = thrIntensityAbs,
             thrIntensityRel = thrIntensityRel,
             intensityExpNormed = intensityExpNormed,
             intensityLibNormed = intensityLibNormed,
             tuneLibSpectra = tuneLibSpectra,
             useMS1ResDefine = useMS1ResDefine,
             resDefineAt = resDefineAt,
             normIntensity = normIntensity,
             intensityNormedMethod = intensityNormedMethod))
}




#' Combine score parameter setup
#'
#' @param cutoff \code{numeric} threshold for a avaliable match of combined score
#' @param weightRT \code{numeric} weight of RT when combining the scors
#' @param weightCCS \code{numeric} weight of CCS when combining the scors
#' @param weightMSMS \code{numeric} weight of MSMS when combining the scors
#' @param scoreMSMS \code{character} MSMS score used for score combination
#'  (only "reverse" or "forward" is supported)
#' @return a \code{CombineParam} object
#'
#' @export
CombineParam <- function(
    cutoff = 0.6,
    weightRT = 0.2,
    weightCCS = 0.4,
    weightMSMS = 0.4,
    scoreMSMS = c("reverse", "forward")
) {
  scoreMSMS <- match.arg(scoreMSMS)
  return(new("CombineParam",
             cutoff = cutoff,
             weightRT = weightRT,
             weightCCS = weightCCS,
             weightMSMS = weightMSMS,
             scoreMSMS = scoreMSMS))
}




