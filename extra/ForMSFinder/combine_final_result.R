######### 1. load the function #############
#   generateMsFinderFormulaDB --------------------------------------------------
setGeneric(name = 'generateMsFinderFormulaDB', 
           def = function(
    raw_data,
    dir_path = '.', 
    file_name = 'MsfinderFormulaDB-VS10.efd'
           ){
             idx <- match(raw_data$id, lib_meta$id)
             
             unique_formula <- unique(lib_meta$formula[idx])
             exact_mass <- sapply(unique_formula, ImmsTools::Calcu_EM)
             # exact_mass <- metfrag_candidate$MonoisotopicMass[idx]
             
             result <- data.frame(exact_mass=exact_mass,
                                  formula=unique_formula, 
                                  pubchem_cid='N/A', 
                                  records=1, 
                                  HMDB='N/A', 
                                  KNApSAcK='N/A', 
                                  ChEBI='N/A', 
                                  DrugBank='N/A', 
                                  SMPDB='N/A', 
                                  YMDB='N/A',
                                  T3DB='N/A',
                                  FoodDB='N/A', 
                                  NANPDB='N/A', 
                                  STOFF='N/A', 
                                  BMDB='N/A', 
                                  LipidMAPS='N/A', 
                                  Urine='N/A', 
                                  Saliva ='N/A', 
                                  Feces='N/A', 
                                  ECMDB='N/A', 
                                  CSF='N/A', 
                                  Serum='N/A', 
                                  PubChem=TRUE, 
                                  PlantCyc='N/A', 
                                  UNPD='N/A',
                                  MINE='N/A', 
                                  stringsAsFactors = F)
             
             colnames(result) <- c("Exact mass", "Formula", "PubChem CID", "Records", "HMDB", "KNApSAcK", "ChEBI",
                                   "DrugBank", "SMPDB", "YMDB", "T3DB", "FooDB", "NANPDB", "STOFF", "BMDB", 
                                   "LipidMAPS", "Urine", "Saliva", "Feces", "ECMDB", "CSF", "Serum",  "PubChem",
                                   "PlantCyc", "UNPD", "MINE")
             
             readr::write_tsv(result, 
                              path = file.path(dir_path, file_name), 
                              col_names = TRUE)
             
             
             
             
             
             return(result)
           }
)
#   generateUserDefineDatabaseMsFinder -----------------------------------------
setGeneric('generateUserDefinedDatabaseMsFinder', 
           def = function(
    raw_data,
    lib_meta,
    is.output = TRUE,
    dir_path = '.', 
    file_name = 'defined_DB_L0608.txt'
           ) {
             
             idx <- match(raw_data$labid, lib_meta$labid)
             # inchikey <- lib_meta$inchikey[idx]
             # short_inchikey <- lib_meta$inchikey1[idx]
             # exact_mass <- lib_meta$exact_mass[idx]
             # formula <- lib_meta$formula[idx]
             # exact_mass <- unname(exact_mass)
             
             # formula_mass_result <- pbapply::pblapply(lib_meta$smiles[idx], function(x){
             #   temp <- generateFormulaMass(smiles = x)
             #   temp
             # }) %>% bind_rows()
             
             final_result <- data.frame(title=lib_meta$labid[idx], 
                                        inchikey=lib_meta$inchikey[idx],
                                        short_inchikey=lib_meta$inchikey1[idx],
                                        pubchem_cid=0,
                                        # formula_mass_result,
                                        exact_mass = lib_meta$exact_mass[idx],
                                        formula = lib_meta$formula[idx],
                                        smiles=lib_meta$smiles[idx],
                                        database_id=lib_meta$labid[idx],
                                        stringsAsFactors = F)
             
             colnames(final_result) <- c('Title', 'InChIKey', 'Short InChIKey', 'PubChem CID', 
                                         'Exact mass', 'Formula', 'SMILES', 'Database ID')
             
             
             if (is.output) {
               readr::write_tsv(final_result, 
                                path = file.path(dir_path, file_name), 
                                col_names = TRUE)
             }
             
             return(final_result)
           }
)
#   runMSFinderPredict ---------------------------------------------------------
setGeneric(name = 'runMsFinderPredict',
           def = function(
    bat_file = 'test.bat',
    msfinder = 'I:/software/Tools/MSFinder/MSFINDER_ver_3.24/MsfinderConsoleApp.exe',
    input_folder = 'F:/MetIMMS_MSFinder/test8',
    output_folder = 'F:/MetIMMS_MSFinder/test8/result',
    method_file = 'F:/MetIMMS_MSFinder/test8/MsfinderConsoleApp-Param.txt'
           ){
             
             cmdl <- paste(msfinder, 'predict', '-i', input_folder, '-o', output_folder, '-m', method_file, sep=' ')
             readr::write_lines(x = cmdl, path = bat_file, append = TRUE)
             
             cat('MSFinder commend has been writen, please run bat file\n')
           }
)


#   generateMsFinderConsoleAppPara ---------------------------------------------
setGeneric(name = 'generateMsFinderConsoleAppPara', 
           def = function(
    file_para = 'MsfinderConsoleApp-Param-L0621.txt',
    template = 'F:/MetIMMS_MSFinder/test8/MsfinderConsoleApp-Param.txt',
    UserDefinedDbFilePath = 'F:/MetIMMS_MSFinder/test8/L0621_db.txt'
           ){
             template <- readr::read_lines(file = template)
             template[54] <- paste0('UserDefinedDbFilePath=', UserDefinedDbFilePath)
             readr::write_lines(template, path = file_para)
           }
)


#   get_the_combined_result_id -------------------------------------------------
get_the_combined_result_id <- function(file_to_match_result, 
                                       file_to_spec_searched){
  cat('\n', '\n', 'reading 4D match result')
  match_result <- readr::read_csv(file_to_match_result)
  cat('\n', '\n', 'reading MS-FINDER result')
  # read the MSFINDER result and assigned
  mf_result <- list.files('./msfinder/', 
                          pattern = '^Structure result-', 
                          recursive = T, 
                          full.names = T)
  
  mf_result <- lapply(mf_result, read.delim2)
  mf_result <- do.call(rbind, mf_result)
  test <- sapply(mf_result$Database, function(i){
    temp <- strsplit(i, ',')[[1]][1]
    temp <- str_replace(temp, pattern = 'Database ID=', 
                        replacement = '')
    return(temp)
  })
  
  names(test) <- NULL
  mf_result$met4dx_id <- test
  
  all_feature <- sapply(mf_result$File.name, function(file){
    strsplit(file, '_\\[')[[1]][1]
  })
  
  mf_result$feature <- all_feature
  rm(all_feature, test)
  
  mf_result <- mf_result %>% 
    filter(Rank <= 3)
  
  mf_result$final_adduct <- mf_result$Precursor.type
  
  # filter the id in higher level #
  cat('\n', '\n', 'filtering MS-FINDER result')
  match_result_filter <- match_result[which(!is.na(match_result$labids)), ]
  match_result_filter_col <- match_result_filter[, c(c('name', 'mz', 'rt', 'ccs'), colnames(match_result_filter)[(ncol(match_result_filter)-21):ncol(match_result_filter)])]
  match_result_filter_col <- match_result_filter_col[, -18]
  match_id <- lapply(seq(nrow(match_result_filter)), function(i){
    strsplit(match_result_filter$labids[i], ';')[[1]]
  })
  # match_id <- as.character(na.omit(unlist(match_id)))
  match_id <- unlist(match_id)
  mf_result <- mf_result[which(!mf_result$met4dx_id %in% match_id), ]
  
  
  # generate the final table # 
  cat('\n', '\n', 'generating final result')
  
  all_unique_feature <- unique(mf_result$feature)
  spec_match_for_msfinder <- readRDS(file_to_spec_searched)
  names(spec_match_for_msfinder) <- sapply(seq(length(spec_match_for_msfinder)), function(i){
    spec_match_for_msfinder[[i]]$dataExp$info$name
  })
  res <- lapply(all_unique_feature, function(ft){
    ids <- mf_result$met4dx_id[which(mf_result$feature == ft)]
    idx <- match(ft, names(spec_match_for_msfinder))
    
    temp_res <- spec_match_for_msfinder[[idx]]$dataRef$info
    check_in <- which(temp_res$labid %in% ids)
    temp_res <- temp_res[check_in, ]
    colnames(temp_res)[1] <- 'labids'
    colnames(temp_res)[3] <- 'compound_name'
    colnames(temp_res)[2] <- 'database'
    temp_res$level <- 3
    colnames(temp_res)[4] <- 'confidence_level'
    colnames(temp_res)[18] <- 'adducts'
    colnames(temp_res)[7] <- 'theo_mz'
    temp_res$scoreReverse <- NA
    temp_res$scoreForward <- NA
    temp_res$score <- NA
    
    temp_res$errorMZ <- round(temp_res$errorMZ, 0)
    temp_res$errorRT <- round(temp_res$errorRT, 0)
    temp_res$errorCCS <- round(temp_res$errorCCS, 3)
    temp_res$scoreRT <- round(temp_res$scoreRT, 2)
    temp_res$scoreCCS <- round(temp_res$scoreCCS, 2)
    
    temp_feature <- spec_match_for_msfinder[[idx]]$dataExp$info[, c('name', 'mz', 'rt', 'ccs')]
    
    temp_res_combined <- lapply(seq(ncol(temp_res)), function(ccc){
      return(paste0(temp_res[, ccc], collapse = ';'))
    })
    temp_res_combined <- as.data.frame(do.call(cbind, temp_res_combined))
    colnames(temp_res_combined) <- colnames(temp_res)
    
    fine_res <- cbind(temp_feature, temp_res_combined)
    idx_col_match <- match(colnames(match_result_filter_col), colnames(fine_res))
    return(fine_res[, idx_col_match])
  })
  res <- do.call(rbind, res)
  
  final_res <- rbind(match_result_filter_col, res)
  
  cat('\n', '\n', 'outputting final result: ScoreCombine.csv')
  write.csv(final_res, './ScoreCombine.csv', row.names = FALSE) # the final result was in ScoreCombine.csv file
  return(NULL)
}
########## 2. load the candidate search result from Met4DX ######
# create a folder and set it as the working directory, containing the "spec_searched", "result3_ScoreCombine_refined_level.csv" and "MsfinderConsoleApp-Param.txt"
setwd('E:/combined_multidimensional_match_reuslt/')
# Install required packages
if(!require(devtools)){
  install.packages("devtools")
}
devtools::install_github('justinZZW/ImmsTools')
install.packages('tidyverse')
library(tidyverse)
library(ImmsTools)

# create the folder call MSFINDER for the generated spectra 
dir.create('msfinder')

# read the intermediate data of candidate search result 
spec_searched <- readRDS('./spec_searched')

########## 3. prepare the msp file and run the MSFINDER ######
# set the directory of the MSFINDER parameter template and MSFINDER console
dir_to_template <- './MsfinderConsoleApp-Param.txt' # path to the parameter set template of MSFINDER, download from our github, putted in the file folder of wd
dir_to_msfinder <- 'E:/MSFINDER_ver_3.24/MsfinderConsoleApp.exe' # path to the "MsfinderConsoleApp.exe" on your windows computer
pol <- 'Positive' # Positive  or Negative mode
# generate the MS2 spectrum to run MS-FINDER
temp <- pbapply::pblapply(seq(length(spec_searched)), function(i){
  # browser()
  cat('\n', i, '\n\n')
  # i <- 1
  cat('Prepare data for retrieve candidates...\n\n')
  temp_data <- spec_searched[[i]]
  # temp_msms <- temp_data$dataExp
  temp_name <- temp_data$dataExp$info$name
  root_path <- file.path(getwd(), '/msfinder', temp_name)
  
  # preparing data ---------------------------------------------
  
  for (adduct_form in unique(temp_data$dataRef$info$adduct)) {
    cat('Start processing', adduct_form, '\n\n')
    
    
    temp_path <- file.path(root_path, paste0(temp_name, '_', adduct_form))
    dir.create(temp_path, recursive = TRUE, showWarnings = FALSE)
    
    
    ImmsTools::GenerateMSP(file_name = file.path(temp_path,
                                                 paste0(temp_name, '_', adduct_form, '.msp')),
                           cmp_name = temp_name,
                           precusormz = as.numeric(temp_data$dataExp$info$mz),
                           ce = '30',
                           adduct = adduct_form,
                           polarity = pol,
                           spec = temp_data$dataExp$spectra[[1]])
    
  
    cat('Generate user defined db\n')
    temp <- generateUserDefinedDatabaseMsFinder(raw_data = temp_data$dataRef$info,
                                                lib_meta = temp_data$dataRef$info,
                                                dir_path = temp_path,
                                                file_name = paste(temp_name, 
                                                                  adduct_form, 
                                                                  'db.txt', 
                                                                  sep = '_'))
    cat('Generate consoleApp param\n')
    generateMsFinderConsoleAppPara(file_para = file.path(temp_path, 
                                                         paste0('MsfinderConsoleApp-Param-', 
                                                                temp_name,
                                                                '.txt')),
                                   template = dir_to_template,
                                   UserDefinedDbFilePath = file.path(temp_path, 
                                                                     paste(temp_name, 
                                                                           adduct_form, 
                                                                           'db.txt', 
                                                                           sep = '_')))
    
    
    runMsFinderPredict(bat_file = './msfinder/run_msfinder.bat',
                       msfinder = dir_to_msfinder,
                       input_folder = temp_path,
                       output_folder = file.path(temp_path, 'result'),
                       method_file = file.path(temp_path, 
                                               paste0('MsfinderConsoleApp-Param-', 
                                                      temp_name, 
                                                      '.txt')))
    
    cat('\n\n')
    
  }
})
# run MS-FINDER through the bat file
shell.exec(file = '.\\msfinder\\run_msfinder.bat')


########## 4. integrate the result from multidimensional match and MSFINDER ####
final_result <- get_the_combined_result_id(file_to_match_result = './result3_ScoreCombine_refined_level.csv', 
                                           file_to_spec_searched = './spec_searched')
# the final result was in ScoreCombine.csv file
################################################################################

