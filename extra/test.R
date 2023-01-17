# library(Met4DX)
# wd <- 'E:/02_4d_metdna/23_ms1_peak_detection/03_6560_data/00_run_met4dx/03_commercial_plasma/'
wd <- 'E:/02_4d_metdna/23_ms1_peak_detection/03_6560_data/00_run_met4dx/02_mouse_liver/'
setwd(wd)

exp <- Experiment(wd = wd, nSlaves = 6, rt_range = c(0, 1440), lc_column = 'HILIC', ion_mode = 'positive')
tims_data <- ImsData(exp)

param <- QueryImsDataParam(rerun = T)
tims_data <- QueryImsData(tims_data, param)

param <- MatchBetweenRunParam(smooth_method = 'loess',
                            mz_tol = 20,
                            frame_range = 40,
                            frame_integration_range = 5,
                            mobility_range = 2,
                            mobility_intgration_range = 0.5,
                            min_points = 3,
                            min_intensity = 0,
                            n_skip = 0,
                            interpolate_method = NULL,
                            peak_span_eim = 5,
                            peak_span_eic = 5,
                            smooth_window_eim = 8,
                            smooth_window_eic = 4,
                            keep_profile = FALSE,
                            allowed_mobility_shift = 3,
                            allowed_rt_shift = 30,  
                            rerun = T)
tims_data <- MatchBetweenRuns_DTIM_IOI(tims_data, param,
                                       calibration_table_path = './calibration_table.csv',
                                       ion_of_interest_path = './sample_library.csv')


param <- QueryImsDataParam(rerun = T)
tims_data <- QueryImsmsData(tims_data, param)


# saveRDS(tims_data,
#         file = './tims_data',
#         version = 2)

### extract IMMSMS #######

# tims_data_file = './results/tmp/tims_data_ms2/S1_MA-d3-c3_SR.mzML'
# precursor_info = tims_data@filled_peaks[which(tims_data@filled_peaks$smp_idx == '1'), ]
# precursor_info <- precursor_info[which(is.number(precursor_info$mz)), ]
# mz_tol = 20
# frame_range = 30
# frame_integration_range = 5
# mobility_range = 1
# mobility_intgration_range = 0.5
# min_points = 2
# min_intensity = 0
# n_skip = 2
# interpolate_method = NULL
# peak_span_eim = 5
# peak_span_eic = 5
# snthreshold = 3
# smooth_window_eim = 16
# smooth_window_eic = 8
# keep_profile = FALSE
# res_define_at = 200
# use_cmpd_id = TRUE
# smooth_method = 'loess'
# skip_invalid_eic_peaks = FALSE
# skip_invalid_eim_peaks = FALSE
# filter_outlier_peaks = FALSE
# allowed_mobility_shift = 2
# allowed_rt_shift = 10
# data_file = NULL

# tims_data <- readRDS('./tims_data')

st <- Sys.time()
param <- ExtractIMMSMSParam(rerun = T, 
                            mz_tol = 20,
                            frame_range = 30,
                            frame_integration_range = 5,
                            mobility_range = 1,
                            mobility_intgration_range = 0.5,
                            min_points = 2,
                            min_intensity = 0,
                            n_skip = 2,
                            interpolate_method = NULL,
                            peak_span_eim = 5,
                            peak_span_eic = 5,
                            snthreshold = 3,
                            smooth_window_eim = 16,
                            smooth_window_eic = 8)
tims_data <- ExtractIMMSMS_DTIM(tims_data, param)
et <- Sys.time()
et - st





# tims_data <- readRDS('./tims_data')
# object <- tims_data
# "peaks" = object@peaks
# "filled_peaks" = object@filled_peaks
# "features" = object@features
# "sample_groups" = object@sample_groups
# "spectra_files" = object@tmp_data_files$spectra_files
# "peak_groups" = object@peak_groups
# min_fraction = 0.5
# min_num_samples = 1
# valid_eic_peak = FALSE
# valid_eim_peak = FALSE
# snthreshold_eic = NULL
# snthreshold_eim = NULL
# quant_method = "max"
# col_max = "area"
# col_quant = "area"
# res_define_at = 200
# 'tfix' = object@tmp_data_files$calibration_table$tfix[1]
# 'beta' = object@tmp_data_files$calibration_table$beta[1]
# bpparam = NULL
# xxx <- readRDS('./results/tmp/spectra/S1_MA-d3-c3_SR.mzML')
# tims_data <- readRDS('./tims_data')

param <- FinalizeFeatureParam(rerun = T,
                              min_fraction = 0.5,
                              min_num_samples = 1,
                              quant_method = "max",
                              col_max = "area",
                              col_quant = "area", 
                              valid_eic_peak = FALSE, 
                              valid_eim_peak = FALSE, 
                              snthreshold_eic = NULL,
                              snthreshold_eim = NULL)
tims_data <- FinalizeFeatures_DTIM(tims_data, param)

# saveRDS(tims_data,
#         file = './tims_data',
#         version = 2)
# 
# tims_data <- readRDS('./tims_data')

param <- FillPeakParam(rerun = T, 
                       mobility_intgration_range = 0.5)
tims_data <- FillPeaks_IOI(tims_data, param)



