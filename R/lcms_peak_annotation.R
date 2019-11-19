#' Peak annotation
#'
#' `MAIT` Package allows to annotate the peaks of the peak table provided by `XCMS`.
#' This is done in three different stages:
#' * Annotation with `CAMERA`
#' * Annotation using predefined biotransformation
#' * Annotation using the Human Metabolome Database (HMDB)
#'
#' *Note: The dataset must be converted from an object of the XCMS Package
#' to an object of the MAIT Package.
#'
#' @inheritParams MAIT::peakAnnotation
#' @return A MAIT-class object with [xsAnnotate-class] in the rawData slot.
#' @export
#' @examples
#' \dontrun{
#  file_name <-  system.file("extdata","peak_table_MAIT.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name)
#' peak_table_ANN <- lcms_peak_annotation(MAIT.object = peak_table)
#' lcms_raw_data(peak_table_ANN)
#' }

lcms_peak_annotation <- function (MAIT.object = NULL,
                             corrWithSamp = 0.7,
                             perfwhm = 0.6,
                             sigma = 6,
                             adductTable = NULL,
                             printSpectraTable = TRUE,
                             corrBetSamp = 0.75,
                             pval = 0.05,
                             calcIso = TRUE,
                             calcCiS = TRUE,
                             calcCaS = TRUE,
                             graphMethod = "hcs",
                             annotateAdducts = TRUE){

        if (is.null(MAIT.object)) {
          stop("No MAIT object was given")
        }
        if (length(MAIT.object@RawData@data) != 1) {
          stop("Raw data is not correct. Try running again signalProcessing function")
        }
        parameters <- list(corrWithSamp, corrBetSamp, perfwhm, sigma,
                           adductTable, pval, calcIso, calcCiS, calcCaS, graphMethod,
                           annotateAdducts)
        names(parameters) <- c("corrWithSamp", "corrBetSamp", "perfwhm",
                               "sigma", "adductTable", "peakAnnotation pvalue", "calcIso",
                               "calcCiS", "calcCaS", "graphMethod", "annotateAdducts")
        MAIT.object@RawData@parameters@peakAnnotation <- parameters
        lcms_writeParameterTable(parameters(MAIT.object), folder = MAIT.object@PhenoData@resultsPath)
        if (is.null(adductTable)) {
          cat("WARNING: No input adduct/fragment table was given. Selecting default MAIT table for positive polarity...",
              fill = TRUE)
          cat("Set adductTable equal to negAdducts to use the default MAIT table for negative polarity",
              fill = TRUE)
          peakAnnEnv <- new.env()
          data(MAITtables, envir = peakAnnEnv)
          adducts <- get(x = "posAdducts", envir = peakAnnEnv)
        }
        else {
          if (adductTable == "negAdducts") {
            cat("adductTable has been set to negAdducts. The default MAIT adducts table for negative polarization is selected...",
                fill = TRUE)
            peakAnnEnv <- new.env()
            data(MAITtables, envir = peakAnnEnv)
            adducts <- get(x = "negAdducts", envir = peakAnnEnv)
          }
          else {
            adducts <- read.csv2(paste(adductTable, ".csv",
                                       sep = ""), dec = ".", header = TRUE, sep = ",")
          }
        }
        resultsPath <- MAIT.object@PhenoData@resultsPath
        xsa <- xsAnnotate(rawData(MAIT.object)$xcmsSet)
        xsaF <- groupFWHM(xsa, perfwhm = perfwhm, sigma = sigma)
        cat("Spectrum build after retention time done", fill = TRUE)
        xsaFA <- findIsotopes(xsaF)
        cat("Isotope annotation done", fill = TRUE)
        xsaFA <- groupCorr(xsaFA, cor_eic_th = corrWithSamp, cor_exp_th = corrBetSamp,
                           calcIso = calcIso, calcCiS = calcCiS, calcCaS = calcCaS,
                           pval = pval, graphMethod = graphMethod)
        cat("Spectrum number increased after correlation done",
            fill = TRUE)
        if (annotateAdducts == TRUE) {
          xsaFA <- findAdducts(xsaFA, rules = adducts, polarity = "positive")
          cat("Adduct/fragment annotation done", fill = TRUE)
        }
        peakList <- getPeaklist(xsaFA)
        peakList <- peakList[order(as.numeric(peakList$pcgroup)),
                             ]
        peakList[, 4] <- peakList[, 4]/60
        peakList[, 5] <- peakList[, 5]/60
        peakList[, 6] <- peakList[, 6]/60
        peakList[, 1:6] <- round(peakList[, 1:6], 4)
        if (printSpectraTable == TRUE) {
          tab <- cbind(peakList[, 1], peakList[, 4], peakList[,
                                                              dim(peakList)[2]])
          rownames(tab) <- rep("", dim(tab)[1])
          colnames(tab) <- c("mass", "rt", "spectra Index")
          if (!file.exists(paste(resultsPath, "Tables", sep = "/"))) {
            if (resultsPath == "") {
              dir.create("Tables")
            }
            else {
              dir.create(paste(resultsPath, "Tables", sep = "/"))
            }
          }
          else {
            cat(" ", fill = TRUE)
            warning(paste("Warning: Folder", paste(resultsPath,
                                                   "Tables", sep = "/"), "already exists. Possible file overwritting."))
          }
          if (resultsPath == "") {
            writeExcelTable(file = tab, file.name = "Tables/Spectra")
          }
          else {
            writeExcelTable(file = tab, file.name = paste(resultsPath,
                                                          "Tables/Spectra", sep = "/"))
          }
        }
        xsaFA <- list(xsaFA)
        names(xsaFA) <- "xsaFA"
        MAIT.object@RawData@data <- xsaFA
        return(MAIT.object)

}





#' Raw data extractor from a MAIT object
#'
#' Function rawData extracts the raw data used to build the MAIT-class object
#' @param MAIT.object A MAIT-class object
#'
#' @return A list containing either a xcmsSet or a xsAnnotate object.
#' @export
#' @examples
#' \dontrun{
#' peak_table_MAIT <- to_MAIT(dataDir = dataDir,
#'         project = project,
#'         preproc_params = preproc_params,
#'         peakTable = peakTable)
#'
#' peak_table_ANN <- peak_annotation(MAIT.object = peak_table_MAIT)
#' lcms_raw_data(peak_table_ANN)
#' }
lcms_raw_data <- function (MAIT.object) {
  MAIT::rawData(MAIT.object)
}

NULL

#' Extract significant features from a MAIT object
#'
#' Function spectralSigFeatures takes a MAIT-class object and obtains which of
#' the variables are significant given a p-value threshold. The parameters of
#' the significant features can ve printed to an output table (TRUE by default).
#' Depending on the number of classes in the data, the function chooses between
#' using ANOVA tests through function spectralAnova, or T-Student tests by using
#' function spectralTStudent.
#' @inheritParams MAIT::spectralSigFeatures
#'
#' @return A MAIT-class object containing the significant features of the scores
#'   slot of MAIT-class object used as an input.
#' @export
#'
#' @examples
#' \dontrun{
#' peak_table_MAIT <- to_MAIT(dataDir = dataDir,
#'         project = project,
#'         preproc_params = preproc_params,
#'         peakTable = peakTable)
#'
#' peak_table_ANN <- peak_annotation(MAIT.object = peak_table_MAIT)
#' peak_table_sig_ANN <- spectralSigFeatures(MAIT.object = peak_table_neg_ANN,
#' pvalue=0.05,
#' p.adj="none",
#' scale=FALSE)
#' }
spectral_sig_features <- function(MAIT.object = NULL,
                                  pvalue = 0.05,
                                  p.adj = "none",
                                  printCSVfile = FALSE,
                                  scale = FALSE,
                                  parametric = TRUE,
                                  var.equal = FALSE,
                                  test.fun = NULL,
                                  jitter = FALSE,
                                  jitter.factor = 1,
                                  jitter.amount = 0,
                                  namefun = NULL){
MAIT::spectralSigFeatures(MAIT.object = MAIT.object,
                          pvalue = pvalue,
                          p.adj = p.adj,
                          printCSVfile = printCSVfile,
                          scale = scale,
                          parametric = parametric,
                          var.equal = var.equal,
                          test.fun = test.fun,
                          jitter = jitter,
                          jitter.factor = jitter.factor,
                          jitter.amount = jitter.amount,
                          namefun = NULL)
}
NULL

#' Build a table of the information related to the significant features
#' contained in a MAIT object
#'
#' Function sigPeaksTable takes an MAIT-class object containing significant
#' feature information and builds a table with the information related to these
#' features.
#' @inheritParams MAIT::sigPeaksTable
#' @return A table containing:
#' First column (mz): Peak mass
#'
#' Second column(mzmin): Minimum peak mass of the peak group.
#'
#' Third column(mzmax): Maximum peak mass of the peak group.
#'
#' Fourth column(rt): Peak retention time (in minutes).
#'
#' Fifth column(rtmin): Minimum peak retention time of the peak group.
#'
#' Sixth column(rtmax): Maximum peak retention time of the peak group.
#'
#' Seventh column(npeaks): Number of samples where the peak has been detected.
#'
#' The columns from the nineth to the column labeled "isotopes" contain number
#' of class samples where the peak has been detected and the intensities of the
#' peak among samples.
#'
#' The isotopes column shows if the peak has been identified as a possible
#' isotope.
#'
#' The adduct column shows which kind of adduct could the peak be.
#'
#' The column labeled pcgroup contains the spectral ID of the peak.
#'
#' The P.adjust column contains the corrected peak p-value using post-hoc
#' methods.
#'
#' The p column shows the peak p-value with no multiple test correction.
#'
#' The Fisher column shows the Fisher test results for the peak. Each of the
#' letters separated by the character "_" corresponds to a class value. Classes
#' having the same letters are indistinguible whereas those having different
#' letters are statistically different clases.
#'
#' The last columns contain the mean and median values for each feature
#' @export
#'
#' @examples
#' \dontrun{
#' peak_table_MAIT <- to_MAIT(dataDir = dataDir,
#'         project = project,
#'         preproc_params = preproc_params,
#'         peakTable = peakTable)
#'
#' peak_table_ANN <- peak_annotation(MAIT.object = peak_table_MAIT)
#' peak_table_sig_ANN <- spectralSigFeatures(MAIT.object = peak_table_neg_ANN,
#' pvalue=0.05,
#' p.adj="none",
#' scale=FALSE)
#' sign_table <- sig_peaks_table(peak_table_sig_ANN)
#' }
sig_peaks_table <- function (MAIT.object=NULL,
                             printCSVfile=FALSE,
                             extendedTable = TRUE,
                             printAnnotation=TRUE){
  MAIT::sigPeaksTable(MAIT.object=MAIT.object,
                      printCSVfile=printCSVfile,
                      extendedTable = extendedTable,
                      printAnnotation=printAnnotation)
}

NULL

#' Metabolite identifier
#'
#' Takes a MAIT object and performs the metabolite search for the significant
#' features
#'
#' @inheritParams MAIT::identifyMetabolites
#'
#' @return An output table is stored in the folder (working
#'   directory)/Tables/SearchTable.csv if printCSVfile is set to TRUE. More info
#'   at metaboliteTable
#' @export
#'
#' @examples
#' \dontrun{
#' peak_table_MAIT <- to_MAIT(dataDir = dataDir,
#'         project = project,
#'         preproc_params = preproc_params,
#'         peakTable = peakTable)
#'
#' peak_table_ANN <- peak_annotation(MAIT.object = peak_table_MAIT)
#' peak_table_sig_ANN <- spectralSigFeatures(MAIT.object = peak_table_neg_ANN,
#' pvalue=0.05,
#' p.adj="none",
#' scale=FALSE)
#' sign_table <- sig_peaks_table(peak_table_sig_ANN)
#' Identif <- identifyMetabolites(MAIT.object = peak_table_sig_ANN,
#'                               peakTolerance = 0.005)
#' }
identify_metabolites <- function (MAIT.object=NULL,
                                  peakTolerance=0.005,
                                  database=NULL,
                                  polarity="positive",
                                  printCSVfile=TRUE){
  MAIT::identifyMetabolites(MAIT.object=MAIT.object,
                            peakTolerance=peakTolerance,
                            database=database,
                            polarity=polarity,
                            printCSVfile=printCSVfile)
}
