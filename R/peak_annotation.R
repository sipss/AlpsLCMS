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
#'
#' @return A MAIT-class object with [xsAnnotate-class] in the rawData slot.
#' @export
#' @examples
#' \dontrun{
#' peak_table_MAIT <- to_MAIT(dataDir = dataDir,
#'         project = project,
#'         preproc_params = preproc_params,
#'         peakTable = peakTable)
#'
#' peak_table_ANN <- peak_annotation(MAIT.object = peak_table_MAIT)
#' raw_data(peak_table_ANN)
#' }

peak_annotation <- function (MAIT.object = NULL,
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
                             annotateAdducts = TRUE)
{
MAIT::peakAnnotation(MAIT.object = MAIT.object,
corrWithSamp = corrWithSamp,
perfwhm = perfwhm,
sigma = sigma,
adductTable = adductTable,
printSpectraTable = printSpectraTable,
corrBetSamp = corrBetSamp,
pval = pval,
calcIso = calcIso,
calcCiS = calcCiS,
calcCaS = calcCaS,
graphMethod = graphMethod,
annotateAdducts = annotateAdducts)
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
#' raw_data(peak_table_ANN)
#' }
raw_data <- function (MAIT.object) {
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
