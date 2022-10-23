library(Met4DX)
wd <- '.'
setwd(wd)

### read the spectra ####
exp <- Experiment(wd = wd, nSlaves = 6, rt_range = c(0, 725), lc_column = 'HILIC', ion_mode = 'negative')
tims_data <- TimsData(exp)
param <- ReadSpectraParam(intensity_from = "ms2_intensity",
                          rerun = F)
tims_data <- ReadSpectraData(tims_data, param)

### MS2 spectral dereplication ####
param <- BinPrecursorParam(rerun = F)
tims_data <- BinPrecursors(tims_data, param)


### Query MS1 data frame ####
param <- QueryTimsDataParam(rerun = F)
tims_data <- QueryTimsData(tims_data, param)


### bottom-up assembly peak detection #####
param <- ExtractIMDataParam(rerun = F,
                            smooth_method = 'loess',
                            snthreshold = 3,
                            smooth_window_eim = 15,
                            order_column = "intensity",
                            peak_span_eim = 13,
                            peak_span_eic = 11,
                            keep_profile = T)
tims_data <- ExtractIMData(tims_data, param)
param <- DereplicatePeaksParam(rerun = F, match_msms = F)
tims_data <- DereplicatePeaks(tims_data, param)

### RT alignment #####
param <- CorrectRTLandmarksParam(rerun = F, ccs_tol = 2)
tims_data <- CorrectRT(tims_data, param, ref_sample = "nist_urine_neg_1_p1-a1_1_4095.d")

### peak grouping ####
param <- GroupDensityParam(rerun = F, plot_density = F, mz_bin_size = 0.015)
dereplication_param <- DereplicatePeaksParam(mz_tol = 0.015/2, mobility_tol = 0.015/2, rt_tol = 5, rerun = F, order_column = 'area')
tims_data <- GroupPeaks(tims_data, param, dereplication_param)

### match sample between runs ####
param <- MatchBetweenRunParam(rerun = F,
                              peak_span_eim = 13,
                              peak_span_eic = 11)
tims_data <- MatchBetweenRuns(tims_data, param)

### finalize feature tbale ####
param <- FinalizeFeatureParam(rerun = F,
                              col_max = 'target_intensity')
tims_data <- FinalizeFeatures(tims_data, param)

### peak filling ####
param <- FillPeakParam(rerun = F)
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
tims_data <- IdentifyPeaks(tims_data,
                           param,
                           match_para,
                           combine_para,
                           rt_exp_file = './rt.csv', # if RT calibration file is provided
                           demo_mode = TRUE)