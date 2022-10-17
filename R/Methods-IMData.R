# project: Met4DX
# File name: Methods-IMData.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/1/26 16:28
# Copyright (c) 2022- ZhuMSLab ALL right reserved

#' @export
setMethod(
  "InterpolateEIM",
  signature = c("IMData", "numeric", "nullOrCharacter"),
  function(object, mobility_range, interpolate_method) {
    data <- data.frame("k0" = object@eim_mobility, "intensity" = object@eim$raw)
    rownames(data) <- rownames(object@profile_data)
    data <- .interpolate_data(data, mobility_range, interpolate_method)
    object@eim$interpolate <- data$intensity
    object@eim_mobility_interpolate <- data$k0
    return(object)
  }
)

#' @export
setMethod(
  "SetInfo",
  signature = c("IMData", "numeric", "numeric", "numeric"),
  function(object, apex_eic, apex_eim, target_intensity, frame_integration_range) {
    eic_start <- max(apex_eic - frame_integration_range, 1)
    eic_end <- min(apex_eic + frame_integration_range, length(object@eic$raw))
    object@info <- c("mz" = .get_weighted_mz(object@eic_mz[eic_start:eic_end],
                                             object@eic$raw[eic_start:eic_end]),
                     "rt" = object@eic_rt[apex_eic],
                     "rtmin" = object@eic_rt[eic_start],
                     "rtmax" = object@eic_rt[eic_end],
                     "mobility" = apex_eim,
                     "target_intensity" = target_intensity,
                     "height" = object@eic$raw[apex_eic],
                     "height_fit" = object@eic$smooth[apex_eic],
                     "area" = sum(object@eic$raw[eic_start:eic_end]) * median(diff(object@eic_rt)))
    return(object);
  }
)

#' @export
setMethod(
  "SetData",
  signature = c("IMData", "character", "numeric"),
  function(object, data_slot, data_range = NULL) {
    switch(data_slot,
           "eic" = {
             is_inrange <- object@eim_mobility >= data_range[1] & object@eim_mobility <= data_range[2]
             object@eic[["raw"]] <- unname(colSums(object@profile_data[is_inrange, , drop = FALSE]))
           },
           "eim" = {
             if (is.integer(data_range) && length(data_range) > 2) {
               is_inrange <- seq_along(object@eic_rt) %in% data_range
             } else {
               is_inrange <- object@eic_rt >= data_range[1] & object@eic_rt <= data_range[2]
             }
             object@eim[["raw"]] <- unname(rowSums(object@profile_data[, is_inrange, drop = FALSE]))
           },
           stop("Undefined data slot: ", data_slot)
    )
    return(object)
  }
)

#' @export
setMethod(
  "SmoothData",
  signature = c("IMData", "character", "GaussianSmoothParam"),
  function(object, data_slot, param) {
    data <- slot(object, data_slot)
    if ("interpolate" %in% names(data)) {
      values <- data$interpolate
    } else {
      values <- data$raw
    }
    slot(object, data_slot)[["smooth"]] <- .smooth_gaussian(data = unname(values),
                                                             param@window,
                                                             param@alpha)
    return(object)
  }
)

#' @export
setMethod(
  "SmoothData",
  signature = c("IMData", "character", "LOESSSmoothParam"),
  function(object, data_slot, param) {
    data <- slot(object, data_slot)
    if ("interpolate" %in% names(data)) {
      values <- data$interpolate
    } else {
      values <- data$raw
    }
    slot(object, data_slot)[["smooth"]] <- .smooth_loess(data = unname(slot(object, data_slot)[["raw"]]),
                                                         param@span,
                                                         param@degree,
                                                         param@window)
    return(object)
  }
)

#' @export
setMethod(
  "DropProfile",
  signature = c("IMData"),
  function(object) {
    object@profile_data <- data.frame()
    return(object)
  }
)

#' @export
setMethod(
  "GetPeakInfo",
  signature = c("IMData"),
  function(object) {
    return(as.data.frame(c(as.list(object@info), object@peak_quality)))
})

#' @export
setMethod(
  "plot",
  signature = c("IMData"),
  function(x, main = NULL, check_peak = FALSE) {
    if (nrow(x@profile_data) == 0) {
      par(mfrow = c(1, 2))
      plot(x@eic_rt, x@eic$raw, type = 'b',
           xlab = 'Retention time (s)', ylab = 'Intensity',
           main = paste(main, "EIC", sep = "-"))
      lines(x@eic_rt, x@eic$smooth, type = 'l', col = 'red')
      if (check_peak) {
        abline(v = x@info["rt"], lty = 2, col = 'blue')
      }
      plot(x@eim_mobility_interpolate, x@eim$interpolate, type = 'b',
           xlab = 'Mobility (cm^2/s)', ylab = 'Intensity',
           main = paste(main, "EIM", sep = " - "))
      lines(x@eim_mobility_interpolate, x@eim$smooth, type = 'l', col = 'red')
    } else {
      stop("Not defined for IMData with profile data yet!")
    }
})