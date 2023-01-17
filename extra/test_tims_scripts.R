wd <- 'E:/02_4d_metdna/23_ms1_peak_detection/02_ioi_tims/'
setwd(wd)
# library(Met4DX)

exp <- Experiment(wd = wd, nSlaves = 6, rt_range = c(0, 725), lc_column = 'HILIC', ion_mode = 'positive')
tims_data <- TimsData(exp)
param <- ReadSpectraParam(intensity_from = "ms2_intensity",
                          rerun = F)
tims_data <- ReadSpectraData(tims_data, param)


param <- QueryTimsDataParam(rerun = F)
tims_data <- QueryTimsData(tims_data, param)


param <- MatchBetweenRunParam(rerun = FALSE, 
                              peak_span_eim = 13,
                              peak_span_eic = 11, min_points = 4, n_skip = 0, 
                              allowed_mobility_shift = 0.02, allowed_rt_shift = 30, 
                              frame_range = 30, 
                              smooth_method = 'loess')
tims_data <- MatchBetweenRuns_IOI(tims_data, param, ion_of_interest_path = './20221223_sample_library - Copy.csv')



param <- FinalizeFeatureIOIParam(rerun = F,
                                 col_max = 'target_intensity')
tims_data <- FinalizeFeatures_IOI(tims_data, param)



param <- FillPeakParam(rerun = F)
tims_data <- FillPeaks_IOI(tims_data, param)






