# project: Met4DX
# File name: Methods-Experiment.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/1/17 13:04
# Copyright (c) 2022- ZhuMSLab ALL right reserved

#' @export
setMethod(
  "TimsData",
  signature = "Experiment",
  function(experiment) {
    # browser()
    wd0 <- getwd()
    setwd(experiment@wd)

    # list raw data files
    raw_files <- .get_tims_files(experiment@wd)
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