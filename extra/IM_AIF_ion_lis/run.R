library(Met4DX)
wd <- '.'
setwd(wd)

#### set an experiment ####
exp <- Experiment(wd = wd, nSlaves = 6, 
                  rt_range = c(0, 725), # RT range
                  lc_column = 'HILIC', #  LC column, either 'HILIC' or 'RP'
                  ion_mode = 'positive')
tims_data <- ImsData(exp)
##### query raw data ####3
param <- QueryImsDataParam(rerun = T)
tims_data <- QueryImsData(tims_data, param)

##### peak detection #####
param <- MatchBetweenRunParam(smooth_method = 'loess',
                              frame_range = 20, # frame range to search precursor
                              frame_integration_range = 5, # frame tolerance to signal integration
                              mobility_range = 2,  # drift time tolerance to search precursor (ms)
                              mobility_intgration_range = 0.5, # drift range to signal integration (ms)
                              min_points = 4, # minimal continous points to form a EIC
                              min_intensity = 0, 
                              n_skip = 0,
                              interpolate_method = NULL,
                              peak_span_eim = 9, # minimum points to find a EIM. The whole mobility range was divided into 1000 portions, and user could set the peak span according to the observed peak width of EIM.
                              peak_span_eic = 9, # minimum points to find a EIC, the unit is point. User can set the peak span according to the observed peak width and cycle time of EIC.
                              smooth_window_eim = 8,
                              smooth_window_eic = 4,
                              keep_profile = FALSE,
                              allowed_mobility_shift = 2, # mobility tolerance to search precursor in the raw data (ms)
                              allowed_rt_shift = 30,  # RT tolerance to search precursor in the raw data (second)
                              rerun = T)
tims_data <- MatchBetweenRuns_DTIM_IOI(tims_data, param,
                                       calibration_table_path = './calibration_table.csv', # path to ccs calibration coefficient tableb
                                       ion_of_interest_path = './precursor_ion_list_pos.csv') # path to the precursor ion list

##### query MS2 data #####
param <- QueryImsDataParam(rerun = T)
tims_data <- QueryImsmsData(tims_data, param)

#### extract ms2 raw data #####
st <- Sys.time()
param <- ExtractIMMSMSParam(rerun = T, 
                            mz_tol = 20,
                            mobility_intgration_range = 0.5, # drift range to signal integration (ms)
                            min_points = 6) 
tims_data <- ExtractIMMSMS_DTIM(tims_data, param)
et <- Sys.time()
et - st

#### finalize data #####
param <- FinalizeFeatureParam(rerun = T,
                              min_fraction = 0.5,
                              quant_method = "max",
                              col_max = "area",
                              col_quant = "area", 
                              valid_eic_peak = FALSE, 
                              valid_eim_peak = TRUE, 
                              snthreshold_eic = NULL,
                              snthreshold_eim = NULL)
tims_data <- FinalizeFeatures_DTIM(tims_data, param)


##### gap filling ####
param <- FillPeakParam(rerun = T, 
                       mobility_intgration_range = 0.5)
tims_data <- FillPeaks_IOI(tims_data, param)



