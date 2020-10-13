#' Get st parameter value for RamclustR
#'
#' @param XObj A lcms_dataset
#' @return recomended st parameter value for RamclustR
#' @export
#'
getRamSt <- function(XObj) {
  featInfo <- featureDefinitions(XObj)
  histo <- hist((featInfo$rtmax-featInfo$rtmin)/2)
  st <- round(median(featInfo$rtmax-featInfo$rtmin)/2,
              digits = 2)
  abline(v=st)
  return(st)
}

#TODO
# ramclustR(ms = paste0(output_dir_node8,"/xdataImputed.csv"),
#           featdelim = "_",
#           st = st,
#           sr = sr,
#           ExpDes = Experiment,
#           normalize = "TIC",
#           sampNameCol = 1,
#           fftempdir = output_dir_node9)
#
# RC <- do.findmain(RC,
#                   nls = adducts_list,
#                   mode = "positive",
#                   mzabs.error = 0.005,
#                   ppm.error = 5,
#                   writeMat = FALSE)
