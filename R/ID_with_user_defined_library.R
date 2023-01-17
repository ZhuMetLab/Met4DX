setGeneric("IdentifyPeaks_input_library",
           function(object, search_param, match_param, combine_param, ...)
             standardGeneric("IdentifyPeaks_input_library"))


#' @export
setMethod(
  "IdentifyPeaks_input_library",
  signature = c("TimsData", "SearchParam", "MatchParam", "CombineParam"),
  function(object,
           search_param, match_param, combine_param,
           lib_file = NULL,
           rt_exp_file=NULL, rt_ref_file=NULL, 
           level3_lib_file = NULL, level3_lib_info = NULL,
           demo_mode = FALSE
  ) {
    # browser()
    wd0 <- getwd()
    object <- tims_data
    setwd(object@experiment@wd)
    object@.processHistory <- c(object@.processHistory, param)
    
    message("Identify peaks...")
    
    files <- .tmp_files(object@files, object@experiment@tmp_dir, "identify_peaks")
    names(files) <- object@files
    
    if (!demo_mode && is.null(lib_file)) {
      stop("Either lib_file must be specified!")
    }
    
    pkg <- getPackageName()
    if (pkg == ".GlobalEnv") {
      pkg <- "Met4DX"
    }
    # 
    # if (demo_mode) {
    #   
    #   lib_data <- SpectraTools:::ReadSpectraLib(system.file(package = pkg, "library",
    #                                                         paste0('tims', object@experiment@ion_mode, '.lib')))
    # } else {
    #   lib_data <- SpectraTools::ParseSpectra(SpectraTools::ParseSpectraParam('msp',
    #                                                                          labelKeep = NULL,
    #                                                                          labelName = NULL,
    #                                                                          autoRename = TRUE,
    #                                                                          resDefineAt = 200,
    #                                                                          thrIntensityAbs = 0,
    #                                                                          ppmPrecursorFilter = 20),
    #                                          lib_file)
    #   
    # }
    # browser()
    lib_data <- SpectraTools::ParseSpectra(SpectraTools::ParseSpectraParam('msp',
                                                                           labelKeep = NULL,
                                                                           labelName = NULL,
                                                                           autoRename = TRUE,
                                                                           resDefineAt = 200,
                                                                           thrIntensityAbs = 0,
                                                                           ppmPrecursorFilter = 20),
                                           lib_file)
    lib_data <- SpectraTools::setRT(lib_data, "OTHER")
    
    exp_data <- SpectraTools::SpectraData(object@features, object@spectra)
    # browser()
    # rt_ref <- rt_exp <- NULL
    # if (search_param@scoreRT) {
    #   if (is.null(rt_ref_file)) {
    #     rt_ref_file <- system.file('rt_calibration',
    #                                paste0(object@experiment@lc_column, '_', object@experiment@ion_mode, '.csv'),
    #                                package = pkg)
    #   }
    #   rt_ref <- read.csv(rt_ref_file, stringsAsFactors = FALSE)
    #   rt_exp <- read.csv(rt_exp_file, stringsAsFactors = FALSE)
    # }
    if (search_param@scoreCCS) {
      adduct_table <- read.csv(system.file('adducts',
                                           paste0("adducts_", object@experiment@ion_mode, ".csv"),
                                           package = pkg), stringsAsFactors = FALSE)
    } else {
      adduct_table <- NULL
    }
    lib_data@ccsInfo$`[M-H]-` <- as.numeric(lib_data@ccsInfo$`[M-H]-`)
    lib_data@ccsInfo$`[M+Na-2H]-` <- as.numeric(lib_data@ccsInfo$`[M+Na-2H]-`)
    lib_data@ccsInfo$`[M+HCOO]-` <- as.numeric(lib_data@ccsInfo$`[M+HCOO]-`)
    # browser()
    
    spec_searched <- SpectraTools::SearchSpectra(exp_data, lib_data, search_param,
                                                 rtcalExp = NULL, rtcalRef = NULL,
                                                 adductTable = adduct_table)
    # score_match <- BiocParallel::bplapply(spec_searched, function(specData) {
    score_match <- lapply(spec_searched, function(specData) { 
      # browser()
      # cat(names(specData))
      # specData <- spec_searched[["#989"]]
      dataExp <- specData$dataExp
      dataRef <- specData$dataRef
      dataRef <- SpectraTools::SpectraData(info = dataRef@info,
                                           spectra = dataRef@spectra)
      SpectraTools::MatchSpectra(dataExp, dataRef, match_param)
    })
    # browser()
    score_match <- score_match[!sapply(score_match, is.null)]
    # browser()
    cat("Plotting MSMS match figures ...\n")
    dirPlot <- file.path(object@experiment@tmp_dir, "MSMSMatchPlot")
    PlotMatchResult(score_match, expinfo = exp_data@info, dirPlot, addname = FALSE)
    scTable <- do.call(rbind, lapply(score_match, SpectraTools::GenOutputScore,
                                     match_param@cutoff, type = "metabolites"))
    
    pkTable <- MergeResTable(exp_data@info, scTable)
    pkTable <- add_level_class(pkTable, lib_data)
    write.csv(pkTable, file.path(object@experiment@res_dir, "result1_MSMSmatch.csv"),
              row.names = FALSE)
    
    cat("Finalizing scores ...\n")
    # score_match <- BiocParallel::bplapply(score_match, function(sc) {
    score_match <- lapply(score_match, function(sc) {
      SpectraTools::CombineScore(sc, combine_param)
    })
    score_match <- score_match[!sapply(score_match, is.null)]
    # browser()
    scTable <- do.call(rbind, lapply(score_match, SpectraTools::GenOutputScore,
                                     type = "metabolites"))
    pkTable <- MergeResTable(exp_data@info, scTable)
    pkTable <- add_level_class(pkTable, lib_data)
    write.csv(pkTable, file.path(object@experiment@res_dir, "result2_ScoreCombine.csv"),
              row.names = FALSE)
    
    ### the final table : remain the highest level for each feature #####
    cat("Refine confidence level ...\n")
    scTable <- remain_the_highest_level_res(scTable, 
                                            lib_data)
    # browser()
    pkTable <- MergeResTable(exp_data@info, scTable)
    pkTable <- add_level_class(pkTable, lib_data)
    write.csv(pkTable, file.path(object@experiment@res_dir, "result3_ScoreCombine_refined_level.csv"),
              row.names = FALSE)
    
    cat("Generate ms1 match result for MS-FIDNER ...\n")
    # browser()
    idx_not_id <- which(!row.names(exp_data@info) %in% row.names(scTable))
    
    new_exp_info <- exp_data@info[idx_not_id, ]
    new_exp_spec <- exp_data@spectra[idx_not_id]
    names(new_exp_spec) <- row.names(new_exp_info)
    new_expe_data <- SpectraTools::SpectraData(new_exp_info, 
                                               new_exp_spec)
    new_adduct_table <- adduct_table[1, ]
    
    # if (demo_mode){
    #   level3_lib_file <- system.file('library',
    #                                  'msfinderpos',
    #                                  package = pkg)
    #   
    #   level3_lib_info <- system.file('library',
    #                                  'ms1info',
    #                                  package = pkg)
    # }else{
    #   level3_lib_file <- level3_lib_file
    #   level3_lib_info <- level3_lib_info
    # }
    
    # level3_db <- gen_lib(lib_file = level3_lib_file, 
    #                      info_file = level3_lib_info, 
    #                      col_lib = c('Ion_mode', 'ExactMass', 'level'))
    # level3_lib_file
    
    level3_db <- SpectraTools::ParseSpectra(SpectraTools::ParseSpectraParam('msp',
                                                                           labelKeep = NULL,
                                                                           labelName = NULL,
                                                                           autoRename = TRUE,
                                                                           resDefineAt = 200,
                                                                           thrIntensityAbs = 0,
                                                                           ppmPrecursorFilter = 20),
                                           level3_lib_file)
    # level3_db@info$`num peaks` <- NULL
    level3_db@ccsInfo$`[M-H]-` <- as.numeric(level3_db@ccsInfo$`[M-H]-`)
    level3_db <- SpectraTools::setRT(level3_db, "OTHER")
    
    new_spec_searched <- SpectraTools::SearchSpectra(new_expe_data, level3_db, search_param,
                                                     rtcalExp = NULL, rtcalRef = NULL,
                                                     adductTable = new_adduct_table)
    
    # browser()
    # ms1_info <- readRDS(level3_lib_info)
    for(i in seq(length(new_spec_searched))){
      # cat(i, '\t')
      idx <- match(new_spec_searched[[i]]$dataRef@info$labid,
                   level3_db@info$labid)
      new_spec_searched[[i]]$dataRef@info$exact_mass <- level3_db@info$mz[idx]
    }
    
    # if(tims_data@experiment@ion_mode == 'positive'){
    #   for(i in seq(length(new_spec_searched))){
    #     # cat(i, '\t')
    #     idx <- match(new_spec_searched[[i]]$dataRef@info$labid,
    #                  ms1_info$id)
    #     new_spec_searched[[i]]$dataRef@info$exact_mass <- ms1_info$mass_cal[idx]
    #   }
    # }else{
    #   for(i in seq(length(new_spec_searched))){
    #     # cat(i, '\t')
    #     idx <- match(new_spec_searched[[i]]$dataRef@info$labid,
    #                  ms1_info$id)
    #     new_spec_searched[[i]]$dataRef@info$exact_mass <- ms1_info$mass_cal[idx]
    #   }
    # }
    # 
    rm(i, idx)
    
    spec_searched <- new_spec_searched
    spec_searched <- lapply(seq(length(spec_searched)), function(i){
      res <- list()
      res$dataExp <- list()
      res$dataExp$info <- data.frame(spec_searched[[i]]$dataExp@info)
      res$dataExp$spectra <- spec_searched[[i]]$dataExp@spectra
      res$dataRef <- list()
      res$dataRef$info <- spec_searched[[i]]$dataRef@info
      return(res)
    })
    # saveRDS(spec_searched, './spec_searched', version = 2)
    
    saveRDS(spec_searched, 
            file.path(object@experiment@res_dir, "spec_searched"), 
            version = 2)
    
    setwd(wd0)
    return(object)
  })



PlotMatchResult <- function(scoreMatch, expinfo, dirPlot,
                            addname = FALSE, plotPNG = FALSE) {
  if (!dir.exists(dirPlot)) {
    dir.create(dirPlot)
  }
  
  # BiocParallel::bplapply(names(scoreMatch), function(nm) {
  lapply(names(scoreMatch), function(nm) {
    # cat(nm, '\t')
    require(ggplot2)
    pkname <- expinfo[nm, "name"]
    matchScore <- scoreMatch[[nm]]
    if (!is.null(matchScore)) {
      filePlot = file.path(dirPlot, pkname)
      SpectraTools::PlotMirror(matchScore, pkname, addname = addname,
                               plotPNG = plotPNG,
                               plotPDF = TRUE,
                               filePlot = filePlot,
                               direction = 'both')
    }
  })
  invisible()
}

MergeResTable <- function(info, res) {
  if (!is.data.frame(res)) {
    res <- SpectraTools:::.Col2Numeric(res)
  }
  tmp <- data.frame(matrix(ncol = ncol(res), nrow = nrow(info)))
  colnames(tmp) <- colnames(res)
  rownames(tmp) <- rownames(info)
  tmp[, sapply(res, is.character)] <- ""
  tmp[, sapply(res, is.numeric)] <- 0
  tmp[rownames(res), ] <- res
  info <- cbind(info, tmp)
  return(info)
}


add_level_class <- function(pkTable, 
                            lib_data){
  # parse the labid and match with the lib_data_info
  extra_info <- lapply(pkTable$labids, function(ids){
    ids <- strsplit(ids, ';')[[1]]
    idx <- match(ids, lib_data@info$labid)
    name <- paste0(lib_data@info$name[idx], collapse = ';')
    database <- paste0(lib_data@info$raw_id[idx], collapse = ';')
    level <- paste0(lib_data@info$level[idx], collapse = ';')
    formula <- paste0(lib_data@info$formula[idx], collapse = ';')
    exact_mass <- paste0(lib_data@info$mz[idx], collapse = ';')
    smiles <- paste0(lib_data@info$smiles[idx], collapse = ';')
    inchi <- paste0(lib_data@info$inchi[idx], collapse = ';')
    inchikey <- paste0(lib_data@info$inchikey[idx], collapse = ';')
    kingdom <- paste0(lib_data@info$kingdom[idx], collapse = ';')
    superclass <- paste0(lib_data@info$superclass[idx], collapse = ';')
    class <- paste0(lib_data@info$class[idx], collapse = ';')
    subclass <- paste0(lib_data@info$subclass[idx], collapse = ';')
    
    return(data.frame(compound_name = name, 
                      database = database,
                      formula = formula, 
                      exact_mass = exact_mass, 
                      smiles = smiles, 
                      inchi = inchi, 
                      inchikey = inchikey,
                      confidence_level = level, 
                      kingdom = kingdom, 
                      superclass = superclass, 
                      class = class, 
                      subclass = subclass))
  })
  extra_info <- do.call(rbind, extra_info)
  return(cbind(pkTable, extra_info))
}

remain_the_highest_level_res <- function(scTable, 
                                         lib_data){
  res_refined <- lapply(seq(nrow(scTable)), function(t){
    temp_row <- scTable[t, ]
    id <- strsplit(temp_row$labids, ';')[[1]]
    n_id <- length(id)
    
    if(n_id <= 1){
      return(temp_row)
    }else{
      # browser()
      level <- as.numeric(lib_data@info$level[match(id,
                                                    lib_data@info$labid)])
      highest_level <- min(level)
      final_idx <- which(level == highest_level)
      
      final_res <- apply(temp_row, 2, function(v){
        remain_v <- strsplit(v, ';')[[1]][final_idx]
        return(paste0(remain_v, collapse = ';'))
      })
      temp_res <- as.data.frame(t(final_res))
      row.names(temp_res) <- row.names(temp_row)
      return(temp_res)
    }
  })
  return(do.call(rbind, res_refined))
}

