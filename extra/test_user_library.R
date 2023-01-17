wd <- 'E:/02_4d_metdna/25_review1_data_processing/20230111_test_user_library/'
setwd(wd)

### read the spectra ####
exp <- Experiment(wd = wd,
                  nSlaves = 6, # Number of threads to be used
                  rt_range = c(0, 725), # RT rang for LC separation in seconds
                  lc_column = 'HILIC', # LC column, either 'HILIC' or 'RP'
                  ion_mode = 'positive')
tims_data <- TimsData(exp)
param <- ReadSpectraParam(intensity_from = "ms2_intensity", # pseudo presursor intensity when fragmentation
                          rerun = FALSE)
tims_data <- ReadSpectraData(tims_data, param)

### MS2 spectral dereplication ####
param <- BinPrecursorParam(rerun = FALSE)
tims_data <- BinPrecursors(tims_data, param)


### Query MS1 data frame ####
param <- QueryTimsDataParam(rerun = FALSE)
tims_data <- QueryTimsData(tims_data, param)


### bottom-up assembly peak detection #####
param <- ExtractIMDataParam(rerun = FALSE,
                            smooth_method = 'loess',
                            snthreshold = 3,
                            smooth_window_eim = 15,
                            order_column = "intensity",
                            peak_span_eim = 13,
                            peak_span_eic = 11,
                            keep_profile = FALSE)
tims_data <- ExtractIMData(tims_data, param)
param <- DereplicatePeaksParam(rerun = FALSE, match_msms = FALSE)
tims_data <- DereplicatePeaks(tims_data, param)

### RT alignment #####
param <- CorrectRTLandmarksParam(rerun = FALSE, ccs_tol = 2)
tims_data <- CorrectRT(tims_data, param, ref_sample = "nist_urine_pos_1_p1-a1_1_4048.d")

### peak grouping ####
param <- GroupDensityParam(rerun = FALSE, plot_density = FALSE, mz_bin_size = 0.015)
dereplication_param <- DereplicatePeaksParam(mz_tol = 0.015/2,
                                             mobility_tol = 0.015/2,
                                             rt_tol = 5,
                                             order_column = 'area',
                                             rerun = FALSE)
tims_data <- GroupPeaks(tims_data, param, dereplication_param)

### match sample between runs ####
param <- MatchBetweenRunParam(rerun = FALSE,
                              peak_span_eim = 13,
                              peak_span_eic = 11)
tims_data <- MatchBetweenRuns(tims_data, param)

### finalize feature tbale ####
param <- FinalizeFeatureParam(rerun = FALSE,
                              col_max = 'target_intensity')
tims_data <- FinalizeFeatures(tims_data, param)

### peak filling ####
param <- FillPeakParam(rerun = FALSE)
tims_data <- FillPeaks(tims_data, param)

### metabolite identification ####
param <- SearchParam(typeCCS = 'percentage',
                     toleranceCCS = c(3,6),
                     toleranceRT = c(30, 90))
match_para <- MatchParam(methodMatch = 'direct',
                         methodScore = 'dp',
                         intensityNormedMethod = 'maximum',
                         cutoff = 0.8)
combine_para <- CombineParam(scoreMSMS = 'reverse')
tims_data <- IdentifyPeaks_input_library(tims_data,
                                         param,
                                         match_para,
                                         combine_para,
                                         rt_exp_file = NULL,
                                         rt_ref_file = NULL, 
                                         lib_file = 'E:/02_4d_metdna/25_review1_data_processing/20230111_test_user_library/20220808_test_level12.msp',
                                         level3_lib_file = './test_level3_db.msp',
                                         demo_mode = FALSE)




