wd <- '.'
setwd(wd)
library(Met4DX)
### set an experiment ####
exp <- Experiment(wd = wd, 
                  nSlaves = 6, 
                  rt_range = c(0, 725), 
                  lc_column = 'HILIC', 
                  ion_mode = 'positive')
tims_data <- TimsData(exp)
#### read spectra ####
param <- ReadSpectraParam(intensity_from = "ms2_intensity",
                          rerun = F)
tims_data <- ReadSpectraData(tims_data, param)

##### query ms1 raw data #####
param <- QueryTimsDataParam(rerun = F)
tims_data <- QueryTimsData(tims_data, param)

##### peak detection #####
param <- MatchBetweenRunParam(rerun = T, 
                              peak_span_eim = 11, # minimum points to find a EIM. The whole mobility range was divided into 1000 portions, and user could set the peak span according to the observed peak width of EIM.
                              peak_span_eic = 11, # minimum points to find a EIC, the unit is point. User can set the peak span according to the observed peak width and cycle time of EIC.
                              min_points = 3, # minimal continous points to form a EIC
                              n_skip = 0, 
                              allowed_mobility_shift = 0.015, # mobility tolerance to search precursor in the raw data (V*s/cm^2)
                              allowed_rt_shift = 15, # RT tolerance to search precursor in the raw data (second)
                              smooth_window_eim = 16,
                              smooth_window_eic = 8,
                              frame_range = 30, # frame range to search precursor
                              smooth_method = 'loess', 
                              frame_integration_range = 5,  # frame range for signal integration
                              mobility_intgration_range = 0.015)  # mobility range for signal integration 
tims_data <- MatchBetweenRuns_IOI(tims_data, param, ion_of_interest_path = './precursor_ion_list_pos.csv')


#### finalize raw data ####
param <- FinalizeFeatureIOIParam(rerun = T,
                                 col_max = 'target_intensity', 
                                 valid_eic_peak = FALSE)
tims_data <- FinalizeFeatures_IOI(tims_data, param)


#### gap filling ####
param <- FillPeakParam(rerun = T)
tims_data <- FillPeaks_IOI(tims_data, param)


### metabolite identification ####
param <- SearchParam(typeCCS = 'percentage',
                     toleranceCCS = c(3,6), # ccs match tolerance (%)
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
                           rt_exp_file = './rt.csv', # if RT calibration file is provided
                           demo_mode = TRUE) # turn on the demo mode to use the metabolite library in Met4DX

