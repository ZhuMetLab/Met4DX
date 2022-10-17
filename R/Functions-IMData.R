# project: Met4DX
# File name: Functions-IMData.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/1/26 16:30
# Copyright (c) 2022- ZhuMSLab ALL right reserved

.smooth_gaussian <- function(data, window = 0.1, alpha = 2.5, ...) {
  tryCatch({
    res <- smoother::smth.gaussian(x = data, window = window, alpha = alpha)
    res[is.na(res)] <- 0
    names(res) <- names(data)
    return(res)
  }, error = function(e) {
    warning("Smoothing with gaussian failed, use raw values instead")
    return(data)
  })
}

.smooth_loess <- function(data, span = 0.75, degree = 2, window = NULL, ...) {
  if (!is.null(window)) {
    len_data <- length(data)
    span <- max(4/len_data, window / len_data)
  }
  df <- data.frame("index" = seq_along(data), "data" = data)
  model <- loess(data~index, data = df, span = span, degree = degree)
  res <- predict(model, df$index)
  return(res)
}

.check_peak <- function(data) {

}