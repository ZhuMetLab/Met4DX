wd <- 'E:/02_4d_metdna/23_ms1_peak_detection/04_tims_dia/'
setwd(wd)
# library(Met4DX)

exp <- Experiment(wd = wd, nSlaves = 6, rt_range = c(0, 725), lc_column = 'HILIC', ion_mode = 'positive')
tims_data <- TimsData(exp)

param <- QueryTimsDataParam(rerun = F)
tims_data <- QueryTimsData(tims_data, param)


param <- MatchBetweenRunParam(rerun = FALSE, 
                              peak_span_eim = 13,
                              peak_span_eic = 7, min_points = 3, n_skip = 0, 
                              allowed_mobility_shift = 0.02, allowed_rt_shift = 30, 
                              frame_range = 30, 
                              smooth_method = 'loess')
tims_data <- MatchBetweenRuns_IOI(tims_data, param, ion_of_interest_path = './20221223_sample_library - Copy.csv')


##### query ms2 spectra ######
param <- QueryTimsDataParam(rerun = F)
tims_data <- QueryTimsData_DIA(tims_data, param)


#### extract IMMS ############
param <- ExtractIMMSMSParam(rerun = FALSE, mobility_intgration_range = 0.01, 
                            min_points = 2,
                            n_skip = 2)
tims_data <- ExtractIMMSMS(tims_data, param)

# xxx <- readRDS('./results/tmp/spectra/urine_pos_dia1_P1-A-1_1_16950.d')
#### finalize data ###########
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
tims_data <- FinalizeFeatures_DIA(tims_data, param)

### fill gaps ################
param <- FillPeakParam(rerun = F)
tims_data <- FillPeaks_IOI(tims_data, param)


param <- SearchParam(typeCCS = 'percentage',
                     toleranceCCS = c(3,6),
                     toleranceRT = c(30, 90))
match_para <- MatchParam(methodMatch = 'direct',
                         methodScore = 'dp',
                         intensityNormedMethod = 'maximum',
                         cutoff = 0.8)
combine_para <- CombineParam(scoreMSMS = 'reverse')

