gen_lib <- function(lib_file, info_file, col_lib = c('Ion_mode', 'ExactMass', 'level')) {
  lib_data <- readRDS(lib_file)
  lib_info <- lib_data@info
  lib_spec <- lib_data@spectra
  rm(lib_data)
  info_data <- readRDS(info_file)
  adduct_cols <- c('[M+H]+', '[M+Na]+', '[M+NH4]+', '[M-H2O+H]+', '[M-H]-', '[M+Na-2H]-', '[M+HCOO]-')
  colnames(info_data)[19:25] <- adduct_cols
  colnames(info_data)[1] <- 'labid'
  rt_cols <- c('rt_rp', 'rt_hilic')

  col_info <- c("labid", "name", 'raw_id', 'formula', 'formula_cal', 'mass_cal', 'monoisotope_mass',
                'smiles', 'inchi', 'inchikey', 'inchikey1',
                'kingdom', 'superclass', 'class', 'subclass')
  info <- cbind(lib_info[, col_lib], info_data[match(lib_info$met4dx_id, info_data$labid), col_info])
  info <- unique(info)
  info$level <- gsub("level", "", info$level)
  info <- SpectraTools:::.Col2Numeric(info)
  info$mz <- round(info$mass_cal, 4)
  info$Ion_mode <- sapply(info$Ion_mode, function(x) ifelse(x == "P", "positive", "negative"))
  info <- info[, c('labid', 'raw_id', "name", 'level', 'Ion_mode',
                   'formula', 'mz',
                   'smiles', 'inchi', 'inchikey', 'inchikey1',
                   'kingdom', 'superclass', 'class', 'subclass')]
  colnames(info) <- tolower(colnames(info))
  rownames(info) <- info$labid

  info_rt <- info_data[match(info$labid, info_data$labid), rt_cols]
  colnames(info_rt) <- toupper(gsub('rt_', '', rt_cols))
  info_ccs <- info_data[match(info$labid, info_data$labid), adduct_cols]
  # lib_data <-  SpectraTools::SpectraData(info=info, rtInfo = info_rt, ccsInfo = info_ccs, spectra=lib_spec)
  lib_data <-  SpectraTools::SpectraData(info=info, rtInfo = info_rt, ccsInfo = info_ccs, spectra=lib_spec)
  lib_data <- SpectraTools::UpdateNames(lib_data, .gen_indexes(lib_data@info))
  return(lib_data)
}