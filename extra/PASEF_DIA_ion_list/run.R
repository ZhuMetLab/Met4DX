wd <- '.'
setwd(wd)
library(Met4DX)
##### set the library ######
exp <- Experiment(wd = wd, nSlaves = 6, rt_range = c(0, 725), lc_column = 'HILIC', ion_mode = 'positive')
tims_data <- TimsData(exp)

##### query MS1 data ######
param <- QueryTimsDataParam(rerun = F)
tims_data <- QueryTimsData(tims_data, param)

##### Peak detection with ion list  ######
param <- MatchBetweenRunParam(rerun = FALSE, 
                              peak_span_eim = 13,# minimum points to find a EIM. The whole mobility range was divided into 1000 portions, and user could set the peak span according to the observed peak width of EIM.
                              peak_span_eic = 7, # minimum points to find a EIC, the unit is point. User can set the peak span according to the observed peak width and cycle time of EIC.
                              min_points = 3, # minimal continous points to form a EIC
                              n_skip = 0, 
                              allowed_mobility_shift = 0.02, # mobility tolerance to search precursor in the raw data (V*s/cm^2)
                              allowed_rt_shift = 30, # RT tolerance to search precursor in the raw data (second)
                              frame_range = 30, # frame ranges to search a EIM and EIC peak
                              smooth_method = 'loess')
tims_data <- MatchBetweenRuns_IOI(tims_data, param, ion_of_interest_path = './precursor_ion_list_pos.csv') # path to the precursor ion list


##### query ms2 spectra ######
param <- QueryTimsDataParam(rerun = F)
tims_data <- QueryTimsData_DIA(tims_data, param)


#### extract IMMS ############
param <- ExtractIMMSMSParam(rerun = FALSE, 
                            mobility_intgration_range = 0.01, # mobility ranges to integrate MS2 fragment signals
                            min_points = 2, # minimal continous points to form a EIC
                            n_skip = 2)
tims_data <- ExtractIMMSMS(tims_data, param)

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

