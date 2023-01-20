library(Met4DX)
wd <- '.' # the wd should be the same as this 
setwd(wd)

### read the spectra ####
exp <- Experiment(wd = wd,
                  nSlaves = 6, # Number of threads to be used
                  rt_range = c(0, 725), # RT rang for LC separation in seconds
                  lc_column = 'HILIC', # LC column, either 'HILIC' or 'RP'
                  ion_mode = 'negative')
tims_data <- TimsData(exp)
param <- ReadSpectraParam(intensity_from = "ms2_intensity", # MS2 spectral intensity
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
                            snthreshold = 3, # the S/N threshold during peak detection
                            smooth_window_eim = 15,
                            order_column = "intensity",
                            peak_span_eim = 13, # minimum points to find a EIM. The whole mobility range was divided into 1000 portions, and user could set the peak span according to the observed peak width of EIM.
                            peak_span_eic = 11, # minimum points to find a EIC, the unit is point. User can set the peak span according to the observed peak width and cycle time of EIC.
                            keep_profile = FALSE)
tims_data <- ExtractIMData(tims_data, param)
param <- DereplicatePeaksParam(rerun = FALSE, match_msms = FALSE)
tims_data <- DereplicatePeaks(tims_data, param)

### RT alignment #####
param <- CorrectRTLandmarksParam(rerun = FALSE, ccs_tol = 2)
tims_data <- CorrectRT(tims_data, param, ref_sample = "nist_urine_neg_1_p1-a1_1_4095.d") # reference sample for landmark-based RT alignment. we recommend to used the QC sample or the middle sample in the injection order

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
                     toleranceCCS = c(3,6), # cca match tolerance (%)
                     toleranceRT = c(30, 90)) # RT match tolerance (second)
match_para <- MatchParam(methodMatch = 'direct',
                         methodScore = 'dp',
                         intensityNormedMethod = 'maximum',
                         cutoff = 0.8) # the cutoff of MS2 spectral match
combine_para <- CombineParam(scoreMSMS = 'reverse')
tims_data <- IdentifyPeaks(tims_data,
                           param,
                           match_para,
                           combine_para,
                           rt_exp_file = './rt_neg.csv', # if RT calibration file is provided
                           demo_mode = TRUE) # turn on the demo mode to use the metabolite library in Met4DX