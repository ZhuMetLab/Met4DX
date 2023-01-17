# library(Met4DX)
wd <- 'E:/02_4d_metdna/25_review1_data_processing/20230110_generateMS2_test/data_processing/'
setwd(wd)

### read the spectra ####
exp <- Experiment(wd = wd,
                  nSlaves = 6, # Number of threads to be used
                  rt_range = c(0, 725), # RT rang for LC separation in seconds
                  lc_column = 'HILIC', # LC column, either 'HILIC' or 'RP'
                  ion_mode = 'positive')
tims_data <- TimsData(exp)


param <- QueryTimsDataParam(rerun = F)
tims_data <- QueryTimsData(tims_data, param)

param <- QueryTimsDataParam(rerun = F)
tims_data <- QueryTimsData_MS2(tims_data, param)


param <- ReadSpectraParam(intensity_from = "ms2_intensity", # pseudo presursor intensity when fragmentation
                          rerun = TRUE)
tims_data <- GenerateMS2(tims_data, param)



### MS2 spectral dereplication ####
param <- BinPrecursorParam(rerun = T)
tims_data <- BinPrecursors(tims_data, param)



### bottom-up assembly peak detection #####
param <- ExtractIMDataParam(rerun = T,
                            smooth_method = 'loess',
                            snthreshold = 3,
                            smooth_window_eim = 15,
                            order_column = "intensity",
                            peak_span_eim = 13,
                            peak_span_eic = 11,
                            keep_profile = FALSE)
tims_data <- ExtractIMData(tims_data, param)
param <- DereplicatePeaksParam(rerun = T, match_msms = FALSE)
tims_data <- DereplicatePeaks(tims_data, param)

### RT alignment #####
param <- CorrectRTLandmarksParam(rerun = T, ccs_tol = 2)
tims_data <- CorrectRT(tims_data, param, ref_sample = "nist_urine_pos_1_p1-a1_1_4048.d")

### peak grouping ####
param <- GroupDensityParam(rerun = T, plot_density = FALSE, mz_bin_size = 0.015)
dereplication_param <- DereplicatePeaksParam(mz_tol = 0.015/2,
                                             mobility_tol = 0.015/2,
                                             rt_tol = 5,
                                             order_column = 'area',
                                             rerun = FALSE)
tims_data <- GroupPeaks(tims_data, param, dereplication_param)

### match sample between runs ####
param <- MatchBetweenRunParam(rerun = T,
                              peak_span_eim = 13,
                              peak_span_eic = 11)
tims_data <- MatchBetweenRuns(tims_data, param)

### finalize feature tbale ####
param <- FinalizeFeatureParam(rerun = T,
                              col_max = 'target_intensity')
tims_data <- FinalizeFeatures(tims_data, param)

### peak filling ####
param <- FillPeakParam(rerun = T)
tims_data <- FillPeaks(tims_data, param)


# ttt <- list.files('./results/tmp/spectra/',full.names = T)
# xxx <- sapply(ttt, function(t){
#   k <- readRDS(t)
#   nrow(k@info)
# })
# 
# 
# ttt <- list.files('./results/tmp/im_data/smooth_loess/', pattern = 'peaks$',full.names = T)
# xxx <- sapply(ttt, function(t){
#   k <- readRDS(t)
#   nrow(k)
# })




