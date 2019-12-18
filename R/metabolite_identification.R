#' Write parameter table
#'
#' Write a parameter Table for MAIT functionalities.
#'
#' @param listParameters Parameters.
#' @param folder A directory to store a .csv file where the  parameter table is stored.
#' @return A template with MAIT parameters.
#' @family metabolite identification functions
#' @family import/export functions
#' @export
#' @examples
#' data_dir <-  system.file("extdata", "rearrange_mait", package = "NIHSlcms")
#' mait_params_path <- system.file("extdata", "results_project", package = "NIHSlcms")
#' opt_result_path <-  system.file("extdata", package = "NIHSlcms")
#' preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#' parameters <- c(dataDir = data_dir, preproc_params)
#'
#' quiet <- function(x) {
#'              base::sink(base::tempfile())
#'              base::on.exit(base::sink())
#'              base::invisible(base::force(x))}
#'
#' quiet(getClassDef("MAIT","MAIT"))
#' MAIT.object = signature(MAIT.Object = "MAIT")
#' MAIT.object <- methods::new("MAIT")
#'
#' MAIT.object@RawData@parameters@sampleProcessing <- parameters
#' lcms_write_parameter_table(MAIT::parameters(MAIT.object), folder = mait_params_path)
lcms_write_parameter_table <- function(listParameters, folder){
  outputTable <- as.matrix(c(unlist(listParameters@sampleProcessing),
                             unlist(listParameters@peakAnnotation),
                             unlist(listParameters@peakAggregation),
                             unlist(listParameters@sigFeatures),
                             unlist(listParameters@biotransformations),
                             unlist(listParameters@identifyMetabolites),
                             unlist(listParameters@classification),
                             unlist(listParameters@plotPCA),
                             unlist(listParameters@plotHeatmap)))
  colnames(outputTable) <- c("Value")
  utils::write.csv(file=paste(folder,"MAITparameters.csv", sep="/"), outputTable)
  return(outputTable)
}

#' Convert xcms to MAIT
#'
#' Function to create a MAIT object using data from xcms.
#'
#' `MAIT` Package allows to annotate the peaks of the peak table provided by `XCMS`.
#' *Note: The dataset must be converted from an object of the XCMS Package
#' to an object of the MAIT Package.
#'
#' @param data_dir directory with LC-MS datafiles.
#' @param project_dir path to the the project directory.
#' @param project name of the.
#' @param preproc_params params from XCMS.
#' @param peak_table XCMS procesed LC-MS data.
#' @return A MAIT object.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family import/export functions
#' @export
#' @examples
#' file_name <-  system.file("extdata", "peak_table_imputed.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name)
#' opt_result_path <-  system.file("extdata", package = "NIHSlcms")
#' preproc_params <- lcms_read_ipo_to_xcms(opt_result_path)
#'
#' data_dir <-  system.file("extdata", "rearrange_mait", package = "NIHSlcms")
#' project_dir <-  system.file("extdata", package = "NIHSlcms")
#' project <- "project"
#'
#' peak_table_mait <- lcms_to_mait(data_dir = data_dir,
#'                                 project_dir = project_dir,
#'                                 project = project,
#'                                 preproc_params = preproc_params,
#'                                 peak_table = peak_table)
#' print(peak_table_mait)
lcms_to_mait <- function (data_dir = NULL, project_dir = NULL, project = NULL, preproc_params = NULL,  peak_table = NULL){

  quiet <- function(x) {
    base::sink(base::tempfile())
    base::on.exit(base::sink())
    base::invisible(base::force(x))
  }

  cat("Building a MAIT object from xcms data.","\n")
  parameters <- c(dataDir = data_dir, preproc_params)
  methods::getClassDef("MAIT","MAIT")
  MAIT.object = methods::signature(MAIT.Object = "MAIT")
  MAIT.object <- methods::new("MAIT")




  if (is.null(data_dir)) {
    stop("No input directory was given")
  }
  if (is.null(preproc_params)) {
    stop("No parameters for preprocessing were included")
  }
  if (is.null(peak_table)) {
    stop("No peak was included")
  }
  #print(MAIT::resultsPath(MAIT.object))
  MAIT.object@RawData@parameters@sampleProcessing <- parameters


  class <- list.files(data_dir)
  classNum <- vector(length = length(class))
  fileList <- list.files(path = paste(data_dir, list.files(path = data_dir),
                                      sep = "/"), full.names = TRUE)
  for (i in 1:length(class)) {
    classNum[i] <- length(list.files(paste(data_dir, list.files(data_dir)[i],
                                           sep = "/")))
  }

  classes <- rep(class, classNum)

  if (length(list.files(data_dir)) == 1) {
    warning("Warning: Input data only has one class!")
  }

  if (is.null(project)) {
    warning("Warning: Project name is empty!")
  }

  if (!is.null(project_dir)) {
    resultsPath <- paste(project_dir, paste("results", project, sep = "_"), sep ="/")
    if (any(base::dir.exists(resultsPath))) {
      #cat("There are already directories / files in the folder. Not saving new ones.")
      cat("\n")
    }else{
      base::dir.create(resultsPath)
      peak_table <- base::suppressMessages(methods::as(peak_table, "xcmsSet"))
      peak_table <-  base::suppressMessages(
        base::suppressWarnings(quiet(xcms::fillPeaks(peak_table))
        )
      )
    }
  }else{
    resultsPath <- "results"
    base::dir.create(resultsPath)
  }

  lcms_write_parameter_table(base::suppressMessages(base::suppressWarnings(
    quiet(MAIT::parameters(MAIT.object)))),
    folder = resultsPath)

  peak_table <- base::suppressMessages(methods::as(peak_table, "xcmsSet"))
  peak_table <-  base::suppressMessages(
    base::suppressWarnings(quiet(xcms::fillPeaks(peak_table))
    )
  )

  fPeaks <- list(peak_table)
  names(fPeaks) <- "xcmsSet"
  MAIT.object@RawData@data <- fPeaks
  MAIT.object@PhenoData@classes <- class
  MAIT.object@PhenoData@classNum <- classNum
  MAIT.object@PhenoData@resultsPath <- resultsPath
  MAIT <- MAIT.object
}

#' Peak Annotation
#'
#' The `MAIT` Package allows to annotate the peaks of the peak table provided by `XCMS`.
#' This is done in three different stages:
#' * Annotation with `CAMERA`.
#' * Annotation using predefined biotransformation.
#' * Annotation using the Human Metabolome Database (HMDB).
#'
#' *Note: The dataset must be converted from an object of the XCMS Package
#' to an object of the MAIT Package.
#'
#' @inheritParams MAIT::peakAnnotation
#' @return A MAIT-class object with xsAnnotate-class in the rawData slot.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family import/export functions
#' @export
#' @examples
#' file_name <-  system.file("extdata","peak_table_mait.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name)
#' peak_table_ann <- lcms_peak_annotation(MAIT.object = peak_table)
#' lcms_raw_data(peak_table_ann)
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
MAITtables <- NULL
  quiet <- function(x) {
    base::sink(base::tempfile())
    base::on.exit(base::sink())
    base::invisible(base::force(x))
  }


  writeExcelTable <- function(file,file.name){
    utils::write.csv(file,paste(file.name,".csv",sep=""))
  }

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
  lcms_write_parameter_table(parameters(MAIT.object), folder = MAIT.object@PhenoData@resultsPath)
  if (is.null(adductTable)) {
    cat("WARNING: No input adduct/fragment table was given. Selecting default MAIT table for positive polarity...",
        fill = TRUE)
    cat("Set adductTable equal to negAdducts to use the default MAIT table for negative polarity",
        fill = TRUE)
    peakAnnEnv <- new.env()
    utils::data(MAITtables, envir = peakAnnEnv)
    adducts <- get(x = "posAdducts", envir = peakAnnEnv)
  }
  else {
    if (adductTable == "negAdducts") {
      cat("adductTable has been set to negAdducts. The default MAIT adducts table for negative polarization is selected...",
          fill = TRUE)
      peakAnnEnv <- new.env()
      utils::data(MAITtables, envir = peakAnnEnv)
      adducts <- get(x = "negAdducts", envir = peakAnnEnv)
    }
    else {
      adducts <- utils::read.csv2(paste(adductTable, ".csv",
                                 sep = ""), dec = ".", header = TRUE, sep = ",")
    }
  }
  resultsPath <- MAIT.object@PhenoData@resultsPath
  quiet(xsa <- CAMERA::xsAnnotate(MAIT::rawData(MAIT.object)$xcmsSet))
  quiet(xsaF <- CAMERA::groupFWHM(xsa, perfwhm = perfwhm, sigma = sigma))
  cat("Spectrum build after retention time done", fill = TRUE)
  quiet(xsaFA <- CAMERA::findIsotopes(xsaF))
  cat("Isotope annotation done", fill = TRUE)
  quiet(xsaFA <- CAMERA::groupCorr(xsaFA, cor_eic_th = corrWithSamp, cor_exp_th = corrBetSamp,
                           calcIso = calcIso, calcCiS = calcCiS, calcCaS = calcCaS,
                           pval = pval, graphMethod = graphMethod))
  cat("Spectrum number increased after correlation done",
      fill = TRUE)
  if (annotateAdducts == TRUE) {
    quiet(xsaFA <- CAMERA::findAdducts(xsaFA, rules = adducts, polarity = "positive"))
    cat("Adduct/fragment annotation done", fill = TRUE)
  }
  peakList <- CAMERA::getPeaklist(xsaFA)
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
    if (!file.exists(paste(resultsPath, "tables", sep = "/"))) {
      if (resultsPath == "") {
        dir.create("Tables")
      }
      else {
        dir.create(paste(resultsPath, "tables", sep = "/"))
      }
    }
    else {
      cat(" ", fill = TRUE)
      #warning(paste("Warning: Folder", paste(resultsPath,
      #                                       "tables", sep = "/"), "already exists. Possible file overwritting."))
    }
    if (resultsPath == "") {
      writeExcelTable(file = tab, file.name = "tables/spectra")
    }
    else {
      writeExcelTable(file = tab, file.name = paste(resultsPath,
                                                    "tables/spectra", sep = "/"))
    }
  }
  xsaFA <- list(xsaFA)
  names(xsaFA) <- "xsaFA"
  MAIT.object@RawData@data <- xsaFA
  return(MAIT.object)

}


#' Raw data extractor from a MAIT object
#'
#' Function lcms_raw_data extracts the raw data used to build the MAIT-class object.
#'
#' @param MAIT.object A MAIT-class object.
#' @return A list containing either a xcmsSet or a xsAnnotate object.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family import/export functions
#' @export
#' @examples
#' file_name <-  system.file("extdata", "peak_table_mait.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name)
#' peak_table_ann <- lcms_peak_annotation(MAIT.object = peak_table)
#'
#' lcms_raw_data(peak_table_ann)
lcms_raw_data <- function (MAIT.object) {
  MAIT::rawData(MAIT.object)
}



#' Extract significant features from a MAIT object
#'
#' Function lcms_spectral_sig_features takes a MAIT-class object and obtains which of
#' the variables are significant given a p-value threshold. The parameters of
#' the significant features can ve printed to an output table (TRUE by default).
#' Depending on the number of classes in the data, the function chooses between
#' using ANOVA test or T-Student test.
#'
#' @inheritParams MAIT::spectralSigFeatures
#' @return A MAIT-class object containing the significant features of the scores
#'   slot of MAIT-class object used as an input.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family import/export functions
#' @export
#' @examples
#' file_name <-  system.file("extdata", "peak_table_ann.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name)
#' peak_table_sig_ann <- lcms_spectral_sig_features(MAIT.object = peak_table,
#'                                           pvalue=0.05,
#'                                           p.adj="none",
#'                                           scale=FALSE)
#' print(peak_table_sig_ann)
lcms_spectral_sig_features <- function(MAIT.object = NULL,
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


  if (length(MAIT::classes(MAIT.object)) == 1) {
    if (is.na(MAIT::classes(MAIT.object)) == TRUE) {
      stop("No class information saved in the MAIT object.")
    }
  }
  if (MAIT::method(MAIT.object) == "") {
    cat("Skipping peak aggregation step...")
    MAIT.object <- lcms_peak_aggregation(MAIT.object, scale = scale)
  }
  if (is.null(test.fun)) {
    if (length(MAIT::classes(MAIT.object)) == 2) {
      if (parametric == TRUE) {
        if (var.equal == TRUE) {
          out <- lcms_spectral_tstudent(pvalue = pvalue, p.adj = p.adj,
                                        MAIT.object = MAIT.object, printCSVfile = printCSVfile)
        }
        else {
          out <- lcms_spectral_welch(pvalue = pvalue, p.adj = p.adj,
                                     MAIT.object = MAIT.object, printCSVfile = printCSVfile)
        }
      }
      else {
        out <- lcms_spectral_wilcox(pvalue = pvalue, p.adj = p.adj,
                                    MAIT.object = MAIT.object, printCSVfile = printCSVfile,
                                    jitter = jitter, jitter.factor = jitter.factor,
                                    jitter.amount = jitter.amount)
      }
    }
    else {
      if (parametric == TRUE) {
        out <- lcms_spectral_anova(pvalue = pvalue, p.adj = p.adj,
                                   MAIT.object = MAIT.object, printCSVfile = printCSVfile)
      }
      else {
        out <- lcms_spectral_kruskal(pvalue = pvalue, p.adj = p.adj,
                                     MAIT.object = MAIT.object, printCSVfile = printCSVfile)
      }
    }
  }
  else {
    out <- lcms_spectral_fun(pvalue = pvalue, p.adj = p.adj, MAIT.object = MAIT.object,
                             printCSVfile = printCSVfile, test.fun = test.fun,
                             namefun = namefun)
  }
  if (length(MAIT::featureSigID(out)) == 0) {
    warning("No significative features found with the selected parameters.")
  }
  else {
    aux <- lcms_sig_peaks_table(out, printCSVfile = printCSVfile)
  }
  return(out)
}

#' Peak aggregation
#'
#' lcms_peak_aggregation function applies a peak aggregation technique to the data of a MAIT-class object.
#'
#' @inheritParams MAIT::peakAggregation
#' @return An MAIT-class object.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @keywords internal
#' @noRd

lcms_peak_aggregation<-function(MAIT.object=NULL,
                                method="None",
                                clases=NULL,
                                samples=NULL,
                                PCAscale=FALSE,
                                PCAcenter=FALSE,
                                scale=FALSE,
                                signVariables=NULL,
                                RemoveOnePeakSpectra=FALSE,
                                printCSVfile=TRUE
){
  if (is.null(MAIT.object)){
    stop("No input MAIT object was given")
  }

  if (is.null(method)){
    stop("No input peak aggregation method was given")
  }

  if(method=="NMF"|method=="PCA"|method=="Mean"|method=="Single"){
    stop("Your MAIT version do not support using peak aggregation measures. Make sure that pagR package is installed and that you are using the correct MAIT version.")
  }

  parameters <- list(method,
                     PCAscale,
                     PCAcenter,
                     scale,
                     signVariables,
                     RemoveOnePeakSpectra)

  names(parameters) <- c("peakAggregation method",
                         "peakAggregation PCAscale",
                         "peakAggregation PCAcenter",
                         "peakAggregation scale",
                         "peakAggregation signVariables",
                         "peakAggregation RemoveOnePeakSpectra")
  MAIT.object@RawData@parameters@peakAggregation <- parameters
  lcms_write_parameter_table(parameters(MAIT.object),folder=MAIT.object@PhenoData@resultsPath)
  classes <- MAIT::classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath#resultsPath(MAIT.object)

  MAIT.object@FeatureInfo@peakAgMethod <- method


  aux <- MAIT::getScoresTable(MAIT.object=MAIT.object,getSpectra=TRUE,getExtendedTable=FALSE)
  data <- aux$scores
  if(scale==TRUE && method!="Single"){

    data<-data/rowMeans(data)
    data[is.nan(as.matrix(data))]<-0

  }
  idGroup <- aux$spectraID

  removeOnePeakSpectra <- function(data,
                                   idGroup){
    data <- as.data.frame(data)
    index<-idGroup[which(diff(idGroup)!=0)]
    spectra <- matrix(nrow=1,ncol=ncol(data))
    peaks <- vector(length=1)

    if(idGroup[length(idGroup)]!=idGroup[length(idGroup)-1]){
      index<-c(index,idGroup[length(idGroup)])
    }

    for (i in c(1:length(index))){
      spec <- as.matrix(data[which(index[i]==idGroup),])
      if(dim(spec)[1]>1){
        spectra <- rbind(spectra,spec)
        peaks <- c(peaks,rep(index[i],length(which(index[i]==idGroup))))
      }
    }
    spectra <- spectra[-1,]
    peaks <- peaks[-1]
    out <- list(spectra,peaks)
    names(out) <- c("spectra","idGroup")
    return(out)
  }


  if(RemoveOnePeakSpectra==TRUE){
    dataWithNoOnePeakSpectra <- removeOnePeakSpectra(data=data,
                                                     idGroup=idGroup)
    data <- dataWithNoOnePeakSpectra$spectra
    idGroup <- dataWithNoOnePeakSpectra$idGroup
  }

  if(!is.null(samples)){
    data <- data[,samples]
  }


  data[is.na(as.matrix(data))]<-0
  colnames(data) <- NULL
  rownames(data) <- NULL

  if(is.null(signVariables)){

    MAIT.object@FeatureData@scores <- data
    MAIT.object@FeatureData@featureID <- idGroup

    MAIT.object@FeatureData@scores <- as.matrix(data)
    MAIT.object@FeatureData@featureID <- idGroup
    MAIT.object@FeatureData@models <- list(rep(1,dim(data)[1]))



  }else{
    index <- vector(length=1)
    noneIdGroup <- vector(length=1)
    for(i in c(1:length(signVariables))){
      index<-c(index,which(signVariables[i]==idGroup))
      noneIdGroup<-c(noneIdGroup,rep(signVariables[i],length(which(signVariables[i]==idGroup))))
    }
    index <- sort(index)
    index <- index[-1]
    noneIdGroup <- noneIdGroup[-1]


    MAIT.object@FeatureData@scores <- as.matrix(data[index,])
    MAIT.object@FeatureData@featureID <- noneIdGroup
    MAIT.object@FeatureData@models <- list(rep(1,dim(data)[1]))
    MAIT.object@FeatureInfo@peakAgMethod <- "None"
  }

  if(printCSVfile==TRUE){
    if(!file.exists(paste(resultsPath,"tables",sep="/"))){

      dir.create(paste(resultsPath,"tables",sep="/"))
    }else{
      cat(" " ,fill=TRUE)
      #         warning(paste("Folder",paste(resultsPath,"tables",sep="/"),"already exists. Possible file overwritting.",sep=" "),fill=TRUE)
    }

    tabl <- MAIT::scores(MAIT.object)
    if(length(MAIT::rawData(MAIT.object))!=0){
      if(is.null(samples)){
        colnames(tabl) <- xcms::sampnames(MAIT::rawData(MAIT.object)[[1]]@xcmsSet)

      }else{
        colnames(tabl) <- xcms::sampnames(MAIT::rawData(MAIT.object)[[1]]@xcmsSet)[samples]
      }

    }else{
      colnames(tabl) <- colnames(MAIT::scores(MAIT.object))
    }
    rownames(tabl) <- paste("S",1:dim(tabl)[1],sep="")
    if (scale==TRUE){
      norm <- "scaled"
    }else{
      norm <- "notScaled"
    }
    utils::write.table(file=paste(resultsPath,paste("tables/dataSet_",norm,".csv",sep=""),sep="/"),x=tabl,row.names=TRUE,col.names=NA,sep=",")
  }

  return(MAIT.object)

}

#' Significant feature information
#'
#' Function lcms_sig_peaks_table takes an MAIT-class object containing significant
#' feature information and builds a table with the information related to these
#' features.
#' @inheritParams MAIT::sigPeaksTable
#' @return A table containing:
#' First column (mz): Peak mass
#' Second column (mzmin): Minimum peak mass of the peak group.
#' Third column (mzmax): Maximum peak mass of the peak group.
#' Fourth column (rt): Peak retention time (in minutes).
#' Fifth column (rtmin): Minimum peak retention time of the peak group.
#' Sixth column (rtmax): Maximum peak retention time of the peak group.
#' Seventh column (npeaks): Number of samples where the peak has been detected.
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
#' The p column shows the peak p-value with no multiple test correction.
#'
#' The Fisher column shows the Fisher test results for the peak. Each of the
#' letters separated by the character "_" corresponds to a class value. Classes
#' having the same letters are indistinguible whereas those having different
#' letters are statistically different clases.
#'
#' The last columns contain the mean and median values for each feature.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family import/export functions
#' @export
#' @examples
#' file_name <-  system.file("extdata", "peak_table_sig_ann.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name)
#' sig_table <- lcms_sig_peaks_table(peak_table,  printCSVfile=FALSE)
#' str(sig_table)
lcms_sig_peaks_table<-function(
  MAIT.object=NULL,
  printCSVfile=FALSE,
  extendedTable=TRUE,
  printAnnotation=TRUE){

  #median = NULL
  if (is.null(MAIT.object)) {
    stop("No MAIT object was given")
  }

  if(length(MAIT::featureSigID(MAIT.object))==0){
    stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and spectralSigFeatures were launched")
  }

  dat <- MAIT::getScoresTable(MAIT.object,getSpectra=TRUE,getExtendedTable=TRUE)

  if(printAnnotation==TRUE){extendedTable<-TRUE}

  if(extendedTable==TRUE){
    peakList <- dat$extendedTable
  }else{
    peakList <- dat$scores
  }
  spectraID <- dat$spectraID


  data <- MAIT::scores(MAIT.object)
  index <- MAIT::featureSigID(MAIT.object)
  classes <- MAIT::classes(MAIT.object)
  classNum <- MAIT::classNum(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath#resultsPath(MAIT.object)
  Fisher <- MAIT::LSDResults(MAIT.object)
  TTs <- MAIT::pvalues(MAIT.object)

  sigPeaksTable <- matrix(nrow=1,ncol=ncol(peakList))
  colnames(sigPeaksTable) <- colnames(peakList)


  if (length(index)!=0){
    if (MAIT::method(MAIT.object)!="None"){
      sigPeaksTable <- peakList[spectraID%in%index,]
    }else{
      sigPeaksTable <- peakList[spectraID%in%unique(spectraID[index]),]

    }

    p <- matrix(nrow=1,ncol=dim(sigPeaksTable)[1])
    fisher <- matrix(nrow=1,ncol=dim(sigPeaksTable)[1])
    if(class(MAIT::classes(MAIT.object))!="logical"|class(MAIT::classNum(MAIT.object))!="logical"){
      if(length(MAIT::classes(MAIT.object))>2){
        for(i in c(1:dim(sigPeaksTable)[1])){
          fisher[i] <- Fisher[as.numeric(sigPeaksTable[i,dim(sigPeaksTable)[2]])]
        }
        fisher <- as.vector(fisher)
        fisherNames <- paste(classes[1],classes[2],sep="_")
        for (i in c(3:length(classNum))){
          fisherNames <- paste(fisherNames,classes[i],sep="_")
        }
        NamesFisher <- paste("Fisher",as.character(fisherNames),sep=" ")
        names(fisher) <- NamesFisher

      }else{
        for(i in c(1:dim(sigPeaksTable)[1])){
          fisher[i] <- NA
        }
      }
    }


    p<-MAIT::pvalues(MAIT.object)
    p <- p[spectraID%in%unique(spectraID[index])]

    p <- matrix(p,ncol=1)
    if(MAIT.object@FeatureData@pvaluesCorrection==""){MAIT.object@FeatureData@pvaluesCorrection<-"none"}
    P.adjust <- stats::p.adjust(p,MAIT.object@FeatureData@pvaluesCorrection)

    pAux <- matrix(ncol=1,nrow=dim(sigPeaksTable)[1])
    pAux.adjust <- matrix(ncol=1,nrow=dim(sigPeaksTable)[1])


    if (MAIT::method(MAIT.object)!="None"){
      for (i in c(1:dim(sigPeaksTable)[1])){
        if(sum(sigPeaksTable$pcgroup==index)==length(index)){
          pAux[i,] <- unique(p[which(index%in%index[i])])
          pAux.adjust[i,] <- unique(P.adjust[which(index%in%index[i])])
        }else{
          pAux[i,] <- p[which(index%in%as.numeric(sigPeaksTable$pcgroup)[i])]
          pAux.adjust[i,] <- P.adjust[which(index%in%as.numeric(sigPeaksTable$pcgroup)[i])]
        }
      }

    }else{
      pAux <- p
      pAux.adjust <- P.adjust
    }

    sigPeaksTable <- cbind(sigPeaksTable,pAux.adjust,pAux)

    sigPeaksTable <- cbind(sigPeaksTable,matrix(fisher,ncol=1))
    if(length(MAIT::classes(MAIT.object))>2){
      colnames(sigPeaksTable)[dim(sigPeaksTable)[2]] <- NamesFisher
    }else{
      colnames(sigPeaksTable)[dim(sigPeaksTable)[2]] <- "Fisher.Test"
    }
    colnames(sigPeaksTable)[dim(sigPeaksTable)[2]-2] <- "P.adjust"
    colnames(sigPeaksTable)[dim(sigPeaksTable)[2]-1] <- "p"



    if(length(MAIT::rawData(MAIT.object))==0&is.null(sigPeaksTable$adduct)==TRUE){
      adduct <- character(length=dim(peakList)[1])
      peakList <- data.frame(peakList,adduct)
      adduct <- character(length=dim(sigPeaksTable)[1])

      sigPeaksTable <- data.frame(sigPeaksTable,adduct)
      sigPeaksTable$adduct<-as.character(sigPeaksTable$adduct)
      peakList$adduct<-as.character(peakList$adduct)
    }


    if(sum(is.na(MAIT.object@FeatureInfo@biotransformations))!=1&dim(MAIT.object@FeatureInfo@biotransformations)[1]!=0){
      aux<-as.data.frame(MAIT.object@FeatureInfo@biotransformations)
      aux[,1]<-as.numeric(unlist(MAIT.object@FeatureInfo@biotransformations[,1]))
      aux[,2]<-as.numeric(unlist(MAIT.object@FeatureInfo@biotransformations[,2]))
      biotrans <- aux
      for(i in c(1:dim(biotrans)[1])){
        if(sigPeaksTable[as.numeric(biotrans[i,7]),]$adduct==""){
          sigPeaksTable[as.numeric(biotrans[i,7]),]$adduct <- as.character(biotrans[i,3])
        }else{
          sigPeaksTable[as.numeric(biotrans[i,7]),]$adduct <- paste(sigPeaksTable[as.numeric(biotrans[i,7]),]$adduct,biotrans[i,3],sep=";")
        }
      }
    }

    sigPeaksTable$adduct <- unlist(sigPeaksTable$adduct)
    means <- vector("list",length=length(index))
    medians <- vector("list",length=length(index))

    if(class(MAIT::classes(MAIT.object))!="logical"|class(MAIT::classNum(MAIT.object))!="logical"){
      Fgroups <- as.factor(rep(MAIT::classes(MAIT.object),MAIT::classNum(MAIT.object)))
      for(i in c(1:dim(sigPeaksTable)[1])){
        if(length(MAIT::rawData(MAIT.object))==0){
          ind <- c(3:(dim(sigPeaksTable)[2]-5))
        }else{
          ind <- c((8+length(MAIT::classes(MAIT.object))):(dim(sigPeaksTable)[2]-6))
        }
        means[[i]] <- stats::aggregate(as.numeric(sigPeaksTable[i,ind])~Fgroups,FUN=base::mean)
        medians[[i]] <- stats::aggregate(as.numeric(sigPeaksTable[i,ind])~Fgroups,FUN=stats::median)
      }

      allMeans <- merge(means[[1]],means[[2]],by="Fgroups")
      if(length(MAIT::rawData(MAIT.object))!=0){
        colnames(allMeans)[2:dim(allMeans)[2]] <- rownames(dat$scores[spectraID%in%unique(spectraID[index]),])[c(1,2)]
      }else{
        colnames(allMeans)[2:dim(allMeans)[2]] <- 1:(dim(allMeans)[2]-1)
      }
      if(length(index)>2){
        for(i in c(3:length(means))){
          allMeans <- merge(allMeans,means[[i]],by="Fgroups")
          if(length(MAIT::rawData(MAIT.object))!=0){
            colnames(allMeans)[i+1] <-  rownames(dat$scores[spectraID%in%unique(spectraID[index]),])[i]
          }else{
            colnames(allMeans)[i+1]  <- as.numeric(colnames(allMeans)[i])+1
          }
        }
      }
      rownames(allMeans)<-paste("Mean Class",allMeans[,1])

      temp<-as.data.frame(t(allMeans[,-1]))
      for(i in c(1:length(MAIT::classes(MAIT.object)))){
        temp[,i]<-as.numeric(as.character(temp[,i]))
      }

      if(length(MAIT::rawData(MAIT.object))!=0){rownames(temp)<-NULL}
      allMedians <- merge(medians[[1]],medians[[2]],by="Fgroups")
      if(length(MAIT::rawData(MAIT.object))!=0){
        colnames(allMedians)[2:dim(allMedians)[2]] <- rownames(dat$scores[spectraID%in%unique(spectraID[index]),])[c(1,2)]
      }else{
        colnames(allMedians)[2:dim(allMedians)[2]] <- 1:(dim(allMedians)[2]-1)
      }
      if(length(index)>3){
        for(i in c(3:length(medians))){
          allMedians <- merge(allMedians,medians[[i]],by="Fgroups")
          if(length(MAIT::rawData(MAIT.object))!=0){
            colnames(allMedians)[i+1] <-  rownames(dat$scores[spectraID%in%unique(spectraID[index]),])[i]
          }else{
            colnames(allMedians)[i+1]  <- as.numeric(colnames(allMedians)[i])+1

          }
        }
      }
      rownames(allMedians)<-paste("Median Class",allMedians[,1])

      tempMed<-as.data.frame(t(allMedians[,-1]))
      for(i in c(1:length(MAIT::classes(MAIT.object)))){
        tempMed[,i]<-as.numeric(as.character(tempMed[,i]))
      }
      sigPeaksTable <- cbind(sigPeaksTable,temp,tempMed)
      if(length(MAIT::rawData(MAIT.object))!=0){rownames(tempMed)<-NULL}
    }else{
      temp <- matrix(rep(NA,dim(sigPeaksTable)[1]),ncol=1)
      tempMed <- matrix(rep(NA,dim(sigPeaksTable)[1]),ncol=1)
      sigPeaksTable <- cbind(sigPeaksTable,temp,tempMed)
      colnames(sigPeaksTable)[c(dim(sigPeaksTable)[2]-1,dim(sigPeaksTable)[2])] <- c("Class Mean", "Class Median")
    }

    if(printCSVfile==TRUE){
      if(!file.exists(paste(resultsPath,"tables",sep="/"))){
        dir.create(paste(resultsPath,"tables",sep="/"))
      }else{
        cat(" " ,fill=TRUE)
        #warning(paste("Folder",paste(resultsPath,"tables",sep="/"),"already exists. Possible file overwritting.",sep=" "))
      }
      utils::write.csv(x=sigPeaksTable,file=paste(resultsPath,"tables/significantFeatures.csv",sep="/"),row.names=FALSE)
    }

  }else{
    stop("There are no significant features for this pvalue")
  }
  return(sigPeaksTable)
}

#' Metabolite search
#'
#' Takes a MAIT object and performs the metabolite search for the significant
#' features
#'
#' @inheritParams MAIT::identifyMetabolites
#' @return An output table is stored in the folder (working
#'   directory)/tables/SearchTable.csv if printCSVfile is set to TRUE. More info
#'   at metaboliteTable.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family import/export functions
#' @export
#' @examples
#' file_name <- system.file("extdata", "peak_table_sig_ann.rds", package = "NIHSlcms")
#' peak_table <- readRDS(file_name)
#' metabololite_table <- lcms_identify_metabolites(MAIT.object = peak_table,
#'                               peakTolerance = 0.005)
lcms_identify_metabolites <- function(MAIT.object=NULL,
                                      peakTolerance=0.005,
                                      database=NULL,
                                      polarity="positive",
                                      printCSVfile=TRUE){
  MAITtables <- NULL

  if (is.null(MAIT.object)) {
    stop("No input MAIT object file was given")
  }

  if(length(MAIT::featureSigID(MAIT.object))==0){
    stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and spectralSignFeatures were launched")
  }
  if(length(MAIT.object@FeatureData@masses)==0 & length(MAIT::rawData(MAIT.object))==0){
    stop("No peak masses found in the MAIT object")
  }

  if (is.null(database)) {
    identMetEnv<-new.env()
    utils::data(MAITtables,envir=identMetEnv)
    Database<-get("Database",envir=identMetEnv)
    dataBase <- Database
    cat("WARNING: No input database table was given. Selecting default MAIT database...",fill=TRUE)
  }else{
    dataBase<-utils::read.csv(paste(database,".csv",sep=""),sep=",",header=TRUE)
    dataBase <- as.matrix(dataBase[order(dataBase[,1]),])
  }
  parameters <- list(peakTolerance,
                     database,
                     polarity)
  names(parameters) <- c("peakTolerance",
                         "database",
                         "polarity")

  MAIT.object@RawData@parameters@identifyMetabolites <- parameters
  lcms_write_parameter_table(parameters(MAIT.object),folder= MAIT.object@PhenoData@resultsPath)

  signSpectra <-MAIT::featureSigID(MAIT.object)
  sigPeaksTable <- sigPeaksTable(MAIT.object,printCSVfile=FALSE)
  resultsPath <- MAIT.object@PhenoData@resultsPath

  if(length(MAIT::rawData(MAIT.object))==0){
    Search <- matrix(nrow=1,ncol=8)
    colnames(Search) <- c("Query Mass","Database Mass (neutral mass)","rt","Adduct","Name","spectra","Biofluid","ENTRY")
  }else{
    Search <- matrix(nrow=1,ncol=9)
    colnames(Search) <- c("Query Mass","Database Mass (neutral mass)","rt","Isotope","Adduct","Name","spectra","Biofluid","ENTRY")
  }
  H.mass <- 1.00794

  temp <- MAIT::getScoresTable(MAIT.object,getExtendedTable=TRUE)
  peakList <- temp$extendedTable
  spec <- temp$spectraID

  signPeaklist <- peakList[spec%in%signSpectra,]

  aux <- signPeaklist[c(-grep("[M+1]",signPeaklist$isotope,fixed=TRUE),-grep("[M+2]",signPeaklist$isotope,fixed=TRUE)),]
  aux <- aux[-which(aux$isotope==""),]

  if(length(MAIT::rawData(MAIT.object))==0){

    peaksIsTP <- matrix(nrow=1,ncol=dim(aux)[2]+2*length(MAIT::classes(MAIT.object)))
    colnames(peaksIsTP) <- c(colnames(aux)[3:(dim(aux)[2]-1)],"p.adj","p",colnames(sigPeaksTable)[dim(sigPeaksTable)[2]-1-length(MAIT::classes(MAIT.object))*2],colnames(sigPeaksTable)[(dim(sigPeaksTable)[2]-2*length(MAIT::classes(MAIT.object))+1):dim(sigPeaksTable)[2]])

  }else{

    peaksIsTP <- matrix(nrow=1,ncol=dim(aux)[2]-7+2*length(MAIT::classes(MAIT.object)))
    colnames(peaksIsTP) <- c(colnames(aux)[8:(dim(aux)[2]-3)],"p.adj","p",colnames(sigPeaksTable)[dim(sigPeaksTable)[2]-length(MAIT::classes(MAIT.object))*2],colnames(sigPeaksTable)[(dim(sigPeaksTable)[2]-2*length(MAIT::classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
  }
  auxAd <- signPeaklist[which(signPeaklist$adduct!=""),]
  auxAd <- auxAd[order(as.numeric(auxAd$pcgroup)),]

  #cat("Metabolite identification initiated",fill=TRUE)
  cat("Metabolite identification initiated")
  cat('\nMetabolite identification in progress ')
  lp <- -1
  for (i in (1:dim(sigPeaksTable)[1])){
    spectra<-sigPeaksTable$pcgroup[i]
    index<-i
    if(polarity=="positive"){
      neutralPeak <- sigPeaksTable$mz[i]-H.mass
    }

    if(polarity=="negative"){
      neutralPeak <- sigPeaksTable$mz[i]+H.mass
    }

    SearchPeak <- lcms_SearchCand(candidate=neutralPeak,dataBase=dataBase,peakTolerance=peakTolerance)
    if(SearchPeak$unk==0){
      if (length(SearchPeak$SearchCand)==5){
        ref <- 1
        SearchPeak$SearchCand <- matrix(SearchPeak$SearchCand,nrow=1)
      }else{
        ref <- dim(SearchPeak$SearchCand)[1]
      }
      for (k in c(1:ref)){

        if(length(MAIT::rawData(MAIT.object))==0){

          Search <- rbind(Search,cbind(round(sigPeaksTable$mz[index],5),as.numeric(SearchPeak$SearchCand[k,4]),round(sigPeaksTable$rt[index],2),sigPeaksTable$adduct[index],as.character(SearchPeak$SearchCand[k,2]),spectra,as.character(SearchPeak$SearchCand[k,5]),as.character(SearchPeak$SearchCand[k,1])))

          if(length(MAIT::classes(MAIT.object))>=2){

            added <- cbind(sigPeaksTable[index,3:(dim(sigPeaksTable)[2]-5-length(MAIT::classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,grep("Fisher.",colnames(sigPeaksTable))],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(MAIT::classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
          }else{
            added <- cbind(sigPeaksTable[index,3:(dim(sigPeaksTable)[2]-5-length(MAIT::classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable$Fisher.Test[index],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(MAIT::classes(MAIT.object))+1):dim(sigPeaksTable)[2]])

          }
          colnames(added) <- colnames(peaksIsTP)
          peaksIsTP <- rbind(peaksIsTP,added)

        }else{
          Search <- rbind(Search,cbind(round(sigPeaksTable$mz[index],5),as.numeric(SearchPeak$SearchCand[k,4]),round(sigPeaksTable$rt[index],2),sigPeaksTable$isotopes[index],sigPeaksTable$adduct[index],as.character(SearchPeak$SearchCand[k,2]),spectra,as.character(SearchPeak$SearchCand[k,5]),as.character(SearchPeak$SearchCand[k,1])))
          added <- cbind(sigPeaksTable[index,8:(dim(sigPeaksTable)[2]-6-length(MAIT::classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,dim(sigPeaksTable)[2]-length(MAIT::classes(MAIT.object))*2],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(MAIT::classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
          colnames(added) <- colnames(peaksIsTP)
          peaksIsTP <- rbind(peaksIsTP,added)
        }
      }

    }else{
      ref <- 1
      if(length(MAIT::rawData(MAIT.object))==0){

        Search <- rbind(Search,cbind(round(sigPeaksTable$mz[index],5),"Unknown",round(sigPeaksTable$rt[index],2),sigPeaksTable$adduct[index],"Unknown",spectra,as.character(SearchPeak$SearchCand[5]),as.character(SearchPeak$SearchCand[1])))
        #			added <- cbind(round(sigPeaksTable[index,3:(dim(sigPeaksTable)[2]-5-length(classes(MAIT.object))*2)],0),sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,dim(sigPeaksTable)[2]-1-length(classes(MAIT.object))*2],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
        added <- cbind(sigPeaksTable[index,3:(dim(sigPeaksTable)[2]-5-length(MAIT::classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,grep("Fisher.",colnames(sigPeaksTable))],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(MAIT::classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
        colnames(added) <- colnames(peaksIsTP)
        peaksIsTP <- rbind(peaksIsTP,added)


      }else{

        Search <- rbind(Search,cbind(round(sigPeaksTable$mz[index],5),"Unknown",round(sigPeaksTable$rt[index],2),sigPeaksTable$isotopes[index],sigPeaksTable$adduct[index],"Unknown",spectra,as.character(SearchPeak$SearchCand[5]),as.character(SearchPeak$SearchCand[1])))

        if(length(MAIT::classes(MAIT.object))>=2){

          added <- cbind(sigPeaksTable[index,8:(dim(sigPeaksTable)[2]-6-length(MAIT::classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,grep("Fisher.",colnames(sigPeaksTable))],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(MAIT::classes(MAIT.object))+1):dim(sigPeaksTable)[2]])

        }else{
          added <- cbind(sigPeaksTable[index,8:(dim(sigPeaksTable)[2]-6-length(MAIT::classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable$Fisher.Test[index],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(MAIT::classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
        }
        colnames(added) <- colnames(peaksIsTP)
        peaksIsTP <- rbind(peaksIsTP,added)

      }
    }
    perc <-round(i/dim(sigPeaksTable)[1]*100)

    #if ((perc %% 10 == 0) && (perc != lp)) { cat(perc,' '); lp <- perc }

  }
  Search <- Search[-1,]
  peaksIsTP <- peaksIsTP[-1,]
  peaksIsTP <- as.matrix(peaksIsTP)
  cat('\nMetabolite identification finished')
  row.names(peaksIsTP) <- 1:dim(peaksIsTP)[1]
  row.names(Search) <- 1:dim(Search)[1]
  peaksIsTP.df <- as.data.frame(peaksIsTP)
  Search.df <- as.data.frame(Search)

  ID <- c(1:dim(Search)[1])
  metaboliteTable <- matrix(nrow=dim(Search)[1])
  rownames(metaboliteTable) <- rep("",dim(metaboliteTable)[1])
  rownames(Search) <- rep("",dim(metaboliteTable)[1])
  rownames(peaksIsTP) <- rep("",dim(peaksIsTP)[1])
  metaboliteTable <- merge(Search.df,peaksIsTP.df,by="row.names")[,-1]
  metaboliteTable[,1] <- round(as.numeric(as.character(metaboliteTable[,1])),4)
  if(length(MAIT::rawData(MAIT.object))==0){

    ind1 <- 9:(dim(metaboliteTable)[2]-2*length(MAIT::classes(MAIT.object))-3)
    ind2<-(dim(metaboliteTable)[2]-2*length(MAIT::classes(MAIT.object))-2):dim(metaboliteTable)[2]
    tem<-metaboliteTable[,ind2]
    metaboliteTable[,(dim(metaboliteTable)[2]+1-length(ind1)):dim(metaboliteTable)[2]] <- metaboliteTable[,ind1]
    colnames(metaboliteTable)[(dim(metaboliteTable)[2]+1-length(ind1)):dim(metaboliteTable)[2]] <- colnames(metaboliteTable)[ind1]
    metaboliteTable[,(dim(metaboliteTable)[2]+1-length(ind1)-length(ind2)):(dim(metaboliteTable)[2]-length(ind1))] <- tem
    colnames(metaboliteTable)[(dim(metaboliteTable)[2]+1-length(ind1)-length(ind2)):(dim(metaboliteTable)[2]-length(ind1))] <- colnames(tem)

  }else{

    ind1 <- 10:(dim(metaboliteTable)[2]-2*length(MAIT::classes(MAIT.object))-3)
    ind2<-(dim(metaboliteTable)[2]-2*length(MAIT::classes(MAIT.object))-2):dim(metaboliteTable)[2]
    tem<-metaboliteTable[,ind2]
    metaboliteTable[,(dim(metaboliteTable)[2]+1-length(ind1)):dim(metaboliteTable)[2]] <- metaboliteTable[,ind1]
    colnames(metaboliteTable)[(dim(metaboliteTable)[2]+1-length(ind1)):dim(metaboliteTable)[2]] <- colnames(metaboliteTable)[ind1]
    metaboliteTable[,(dim(metaboliteTable)[2]+1-length(ind1)-length(ind2)):(dim(metaboliteTable)[2]-length(ind1))] <- tem
    colnames(metaboliteTable)[(dim(metaboliteTable)[2]+1-length(ind1)-length(ind2)):(dim(metaboliteTable)[2]-length(ind1))] <- colnames(tem)
  }

  if(printCSVfile==TRUE){

    if(!file.exists(paste(resultsPath,"tables",sep="/"))){
      dir.create(paste(resultsPath,"tables",sep="/"))
      utils::write.table(metaboliteTable,paste(paste( MAIT.object@PhenoData@resultsPath,"tables","metaboliteTable",sep="/"),".csv",sep=""),col.names=NA,row.names=TRUE,sep=",")

    }else{

      cat(" " ,fill=TRUE)
      #cat(paste("Warning: Folder",paste(resultsPath,"tables",sep="/"),"already exists. Possible file overwritting.",sep=" "),fill=TRUE)
      #warning(paste("Folder",paste(resultsPath,"tables",sep="/"),"already exists. Possible file overwritting.",sep=" "))
      utils::write.table(metaboliteTable,paste(paste( MAIT.object@PhenoData@resultsPath,"tables","metaboliteTable",sep="/"),".csv",sep=""),col.names=NA,row.names=TRUE,sep=",")

    }
  }

  MAIT.object@FeatureInfo@metaboliteTable <- metaboliteTable
  return(MAIT.object)

}


#' Search metabolite candidetes
#'
#' Function SearchCand looks up for a peak into a database.
#'
#' @inheritParams MAIT::SearchCand
#' @return A matrix containing all the possible hits for that peak candidateiteTable
#' @family metabolite identification functions.
#' @noRd
lcms_SearchCand <- function(candidate,
                            dataBase,
                            peakTolerance){

  dataBase <- as.matrix(dataBase[order(as.numeric(as.character(dataBase[,4]))),])
  aux <- c(as.numeric(dataBase[,4]),(candidate-peakTolerance))
  #	lowIndex<- which(order(as.numeric(aux))==length(aux))
  lowIndex <- which(sort(aux)%in%(candidate-peakTolerance))
  aux <- c(as.numeric(dataBase[,4]),candidate+peakTolerance)
  highIndex <- which(sort(aux)%in%(candidate+peakTolerance))
  #highIndex<- which(order(as.numeric(aux))==length(aux))
  if(lowIndex!=highIndex){

    unk <- 0
    if(!is.matrix(dataBase[lowIndex:(highIndex-1),])){
      dataBase[lowIndex:(highIndex-1),] <- matrix(dataBase[lowIndex:(highIndex-1),],nrow=1)
    }
    out <- list(dataBase[lowIndex:(highIndex-1),],unk)
    names(out) <- c("SearchCand","unk")
  }else{
    unk <- 1
    out <- list(as.matrix(rep("unknown",5),nrow=1),unk)
    names(out) <- c("SearchCand","unk")

  }
  return(out)
}


#' Extract Significant Features From A MAIT Object For Two Classes
#'
#' Function slcms_spectral_tstudent takes a MAIT-class object and obtains which of the variables are significant
#' given a p-value threshold when there only are two classes in the raw data.
#' The parameters of the significant features can be printed to an output table (TRUE by default).
#' @inheritParams MAIT::spectralTStudent
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family statistical functions
#' @return A MAIT-class object containing the significant features of the scores slot of MAIT-class object used as an input.
#' @keywords internal
#' @noRd

lcms_spectral_tstudent <- function(MAIT.object=NULL,
                                   pvalue=0.05,
                                   p.adj="none",
                                   printCSVfile=TRUE){

  if (is.null(MAIT.object)) {
    stop("No input MAIT object was given")
  }
  parameters <- list(pvalue,
                     p.adj)
  names(parameters) <- c("T-Student pvalue",
                         "T-Student p.adj")
  MAIT.object@RawData@parameters@sigFeatures <- parameters
  lcms_write_parameter_table(parameters(MAIT.object),folder=MAIT.object@PhenoData@resultsPath)

  data <- MAIT::scores(MAIT.object)
  clases <- MAIT::classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  #xsaFA <- rawData(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath#resultsPath(MAIT.object)

  auxs<-MAIT::getScoresTable(MAIT.object=MAIT.object,getExtendedTable=TRUE)
  peakList <- auxs$extendedTable
  spec <- auxs$spectraID

  classes <- c(rep(clases[1],classNum[1]),rep(clases[2],classNum[2]))
  numbers <- classNum
  Bonferroni <- matrix(ncol=as.numeric(peakList$pcgroup[dim(peakList)[1]]),nrow=1)
  names(numbers) <- clases
  TTs <- matrix(ncol=1,nrow=dim(data)[1])
  colnames(TTs) <- paste(clases[1],"&",clases[2],sep="")
  Tresults <- matrix(ncol=1,nrow=spec[dim(peakList)[1]])
  group <- as.factor(c(rep(clases[1],numbers[1]),rep(clases[2],numbers[2])))

  for (i in c(1:(dim(data)[1]))){
    numbers <- classNum
    names(numbers) <- clases
    lmdata <- as.vector(t(data[i,]))
    model <- stats::lm(lmdata~group)
    if (length(which(diff(data)!=0))){
      ttest <- stats::t.test(as.vector(t(data[i,]))~group,var.equal=TRUE)
      TTs[i] <- ttest$p.value
      Tresults[i] <- ttest$statistic
    }else{
      TTs[i] <- 1
    }
  }
  MAIT.object@FeatureData@pvalues <- TTs
  p.corr <- stats::p.adjust(TTs,method=p.adj)
  if (p.adj!="none"){
    index <- which(p.corr<=pvalue)

  }else{
    index<-which(TTs<=pvalue)
  }

  MAIT.object@FeatureData@pvaluesCorrection <- p.adj
  MAIT.object@FeatureData@featureSigID<-index

  return(MAIT.object)
}

#' Extract Significant Features From A MAIT Object For Two Classes
#'
#' Function lcms_spectral_welch takes an MAIT-class object and obtains which of the variables are significant given a p-value
#' threshold following a Welch test. The parameters of the significant features can ve printed to an output table (TRUE by default).
#' @inheritParams MAIT::spectralWelch
#' @return A MAIT-class object containing the significant features of the scores slot of MAIT-class object used as an input.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family statistical functions
#' @keywords internal
#' @noRd
lcms_spectral_welch <- function(MAIT.object=NULL,
                                pvalue=0.05,
                                p.adj="none",
                                printCSVfile=TRUE){

  if (is.null(MAIT.object)) {
    stop("No input MAIT object was given")
  }

  parameters <- list(pvalue,
                     p.adj)
  names(parameters) <- c("Welch pvalue",
                         "Welch p.adj")

  MAIT.object@RawData@parameters@sigFeatures <- parameters
  lcms_write_parameter_table(parameters(MAIT.object),folder=MAIT.object@PhenoData@resultsPath)


  data <- MAIT::scores(MAIT.object)
  clases <- MAIT::classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  #xsaFA <- rawData(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath
  #	peakList <- getPeaklist(MAIT.object)
  #	peakList <- peakList[order(as.numeric(peakList$pcgroup)),]

  auxs<-MAIT::getScoresTable(MAIT.object=MAIT.object,getExtendedTable=TRUE)
  peakList <- auxs$extendedTable
  spec <- auxs$spectraID

  classes <- c(rep(clases[1],classNum[1]),rep(clases[2],classNum[2]))
  numbers <- classNum
  Bonferroni <- matrix(ncol=as.numeric(peakList$pcgroup[dim(peakList)[1]]),nrow=1)
  names(numbers) <- clases
  TTs <- matrix(ncol=1,nrow=dim(data)[1])
  colnames(TTs) <- paste(clases[1],"&",clases[2],sep="")
  #	Tresults <- matrix(ncol=1,nrow=as.numeric(peakList$pcgroup[dim(peakList)[1]]))
  Tresults <- matrix(ncol=1,nrow=spec[dim(peakList)[1]])
  group <- as.factor(c(rep(clases[1],numbers[1]),rep(clases[2],numbers[2])))

  for (i in c(1:(dim(data)[1]))){
    numbers <- classNum
    names(numbers) <- clases
    lmdata <- as.vector(t(data[i,]))
    model <- stats::lm(lmdata~group)
    if (length(which(diff(data)!=0))){
      ttest <- stats::t.test(as.vector(t(data[i,]))~group,var.equal=FALSE)
      TTs[i] <- ttest$p.value
      Tresults[i] <- ttest$statistic
    }else{
      TTs[i] <- 1
    }
  }

  MAIT.object@FeatureData@pvalues <- TTs
  p.corr <- stats::p.adjust(TTs,method=p.adj)

  if (p.adj!="none"){
    index <- which(p.corr<=pvalue)
  }else{
    index<-which(TTs<=pvalue)
  }
  MAIT.object@FeatureData@pvaluesCorrection <- p.adj
  MAIT.object@FeatureData@featureSigID<-index
  return(MAIT.object)
}


#' Extract Significant Features From A MAIT Object For Two Classes
#'
#' Function lcms_spectral_wilcox takes an MAIT-class object and obtains which of the variables are significant
#' given a p-value threshold following a Mann-Witney-Wilcoxon test.
#' The parameters of the significant features can ve printed to an output table (TRUE by default).
#' @inheritParams MAIT::spectralWilcox
#' @return A MAIT-class object containing the significant features of the scores slot of MAIT-class object used as an input.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family statistical functions
#' @keywords internal
#' @noRd
lcms_spectral_wilcox <- function(MAIT.object=NULL,
                                 pvalue=0.05,
                                 p.adj="none",
                                 printCSVfile=TRUE,
                                 jitter=FALSE,
                                 jitter.factor=1,
                                 jitter.amount=0){
  if (is.null(MAIT.object)) {
    stop("No input MAIT object was given")
  }
  parameters <- list(pvalue,
                     p.adj)
  names(parameters) <- c("Mann-Whitney pvalue",
                         "Mann-Whitney p.adj")

  MAIT.object@RawData@parameters@sigFeatures <- parameters
  lcms_write_parameter_table(parameters(MAIT.object),folder = MAIT.object@PhenoData@resultsPath)

  data <- MAIT::scores(MAIT.object)
  clases <- MAIT::classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  #xsaFA <- rawData(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath#resultsPath(MAIT.object)


  auxs<-MAIT::getScoresTable(MAIT.object=MAIT.object,getExtendedTable=TRUE)
  peakList <- auxs$extendedTable
  spec <- auxs$spectraID

  classes <- c(rep(clases[1],classNum[1]),rep(clases[2],classNum[2]))
  numbers <- classNum
  #Bonferroni <- matrix(ncol=as.numeric(peakList$pcgroup[dim(peakList)[1]]),nrow=1)
  names(numbers) <- clases
  TTs <- matrix(ncol=1,nrow=dim(data)[1])
  colnames(TTs) <- paste(clases[1],"&",clases[2],sep="")
  Tresults <- matrix(ncol=1,nrow=spec[dim(peakList)[1]])
  group <- as.factor(c(rep(clases[1],numbers[1]),rep(clases[2],numbers[2])))

  for (i in c(1:(dim(data)[1]))){
    numbers <- classNum
    names(numbers) <- clases
    lmdata <- as.vector(t(data[i,]))
    model <- stats::lm(lmdata~group)
    if (length(which(diff(data)!=0))){
      if(jitter==TRUE){
        test<-stats::wilcox.test(jitter(data[i,],factor = jitter.factor, amount = jitter.amount)~group)
      }else{
        test<-stats::wilcox.test(data[i,]~group)
      }
      TTs[i] <- test$p.value
    }else{
      TTs[i] <- 1
    }
  }
  MAIT.object@FeatureData@pvalues <- TTs
  p.corr <- stats::p.adjust(TTs,method=p.adj)

  if (p.adj!="none"){

    index <- which(p.corr<=pvalue)
  }else{

    index<-which(TTs<=pvalue)
  }
  MAIT.object@FeatureData@pvaluesCorrection <- p.adj
  MAIT.object@FeatureData@featureSigID<-index

  return(MAIT.object)
}


#' Extract Significant Features From A MAIT Object
#'
#' Function lcms_spectral_anova takes an MAIT-class object and obtains which of the variables are significant given a p-value threshold.
#' The parameters of the significant features can ve printed to an output table (TRUE by default).
#' @inheritParams MAIT::spectralAnova
#' @return A MAIT-class object containing the significant features of the scores slot of MAIT-class object used as an input.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family statistical functions
#' @keywords internal
#' @noRd
lcms_spectral_anova <- function (pvalue=0.05,
                                 p.adj="none",
                                 MAIT.object=NULL,
                                 printCSVfile=TRUE)
{

  if (is.null(MAIT.object)) {
    stop("No input MAIT object was given")
  }

  parameters <- list(pvalue,
                     p.adj)
  names(parameters) <- c("Anova p-value",
                         "Anova p.adj")
  MAIT.object@RawData@parameters@sigFeatures <- parameters
  lcms_write_parameter_table(parameters(MAIT.object),folder=MAIT.object@PhenoData@resultsPath)

  data <- MAIT::scores(MAIT.object)
  clases <- MAIT::classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath#resultsPath(MAIT.object)
  auxs<-MAIT::getScoresTable(MAIT.object=MAIT.object,getExtendedTable=TRUE)
  peakList <- auxs$extendedTable

  Fgroups <- matrix(nrow=1)
  Fgroups <- rep(clases,classNum)
  Fgroups <- as.factor(Fgroups[order(Fgroups)])
  #        peakList <- peakList[order(as.numeric(peakList$pcgroup)),]

  TTs <- matrix(ncol=1,nrow=dim(data)[1])
  Fisher <- matrix(ncol=1,nrow=dim(data)[1])

  Tresults <- matrix(ncol=1,nrow=as.numeric(peakList$pcgroup[dim(peakList)[1]]))


  FisherLSD <- function(data,classes,index,DFerror,MSerror,numClasses){
    LSD <- agricolae::LSD.test(y=data[index,],trt=classes,DFerror=DFerror,MSerror=MSerror);
    LSD$groups$trt<- rownames(LSD$groups)
    LSD$groups$M<- LSD$groups$groups
    LSD <- LSD$groups[order(LSD$groups$trt),]
    groups <- paste(LSD$M[1],LSD$M[2],sep=" ")
    for (i in c(3:(numClasses))){
      groups <- paste(groups,LSD$M[i],sep=" ")
    }
    out <- list(LSD$trt,groups,LSD$means)
    names(out) <- c("trt","group","means")
    return(out)
  }


  for (i in c(1:dim(data)[1])){
    lmdata <- as.vector(t(data[i,]))
    numbers <- classNum
    names(numbers) <- clases
    model <- stats::lm(lmdata~Fgroups)
    an <- stats::anova(model)
    Fisher[i] <- FisherLSD(data=data,DFerror=an$Df[2],MSerror=an$Mean[2],index=i,classes=Fgroups,numClasses=length(classNum))[[2]]
    TTs[i] <- an$Pr[1]
  }
  MAIT.object@FeatureData@pvalues <- TTs
  MAIT.object@FeatureData@LSDResults <- Fisher

  p.corr <- stats::p.adjust(TTs,method=p.adj)

  if (p.adj!="none"){
    index <- which(p.corr<=pvalue)
  }else{
    index<-which(TTs<=pvalue)
  }
  MAIT.object@FeatureData@pvaluesCorrection <- p.adj
  MAIT.object@FeatureData@featureSigID<-index
  return(MAIT.object)
}


#' Extract Significant Features From A MAIT Object using an user-defined test function
#'
#' Function lcms_spectral_kruskal takes an MAIT-class object and obtains which of the variables are significant given a
#' p-value threshold following a Kruskal-Wallis test.
#' The parameters of the significant features can ve printed to an output table (TRUE by default).
#' @inheritParams MAIT::spectralKruskal
#' @return A MAIT-class object containing the significant features of the scores slot of MAIT-class object used as an input.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family statistical functions
#' @keywords internal
#' @noRd
lcms_spectral_kruskal <- function (pvalue=0.05,
                                   p.adj="none",
                                   MAIT.object=NULL,
                                   printCSVfile=TRUE)
{
  if (is.null(MAIT.object)) {
    stop("No input MAIT object was given")
  }

  parameters <- list(pvalue,
                     p.adj)
  names(parameters) <- c("Kruskal-Wallis p-value",
                         "Kruskal-Wallis  p.adj")

  MAIT.object@RawData@parameters@sigFeatures <- parameters
  lcms_write_parameter_table(parameters(MAIT.object),folder=MAIT.object@PhenoData@resultsPath)

  data <- MAIT::scores(MAIT.object)
  clases <- MAIT::classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath #resultsPath(MAIT.object)
  auxs<-MAIT::getScoresTable(MAIT.object=MAIT.object,getExtendedTable=TRUE)
  peakList <- auxs$extendedTable

  Fgroups <- matrix(nrow=1)
  Fgroups <- rep(clases,classNum)
  Fgroups <- as.factor(Fgroups[order(Fgroups)])

  TTs <- matrix(ncol=1,nrow=dim(data)[1])
  Fisher <- matrix(ncol=1,nrow=dim(data)[1])

  Tresults <- matrix(ncol=1,nrow=as.numeric(peakList$pcgroup[dim(peakList)[1]]))

  for (i in c(1:dim(data)[1])){

    lmdata <- as.vector(data[i,])
    numbers <- classNum
    names(numbers) <- clases
    #	   model <- lm(lmdata~Fgroups)
    an <- stats::kruskal.test(x=lmdata,g=Fgroups)
    TTs[i] <- an$p.value
  }

  MAIT.object@FeatureData@pvalues <- TTs
  #    MAIT.object@FeatureData@LSDResults <- Fisher
  p.corr <- stats::p.adjust(TTs,method=p.adj)

  if (p.adj!="none"){
    index <- which(p.corr<=pvalue)
  }else{
    index<-which(TTs<=pvalue)
  }

  MAIT.object@FeatureData@pvaluesCorrection <- p.adj
  MAIT.object@FeatureData@featureSigID<-index

  return(MAIT.object)
}

#' Extract Significant Features From A MAIT Object using an user-defined test function.
#'
#' Function lcms_spectral_fun takes an MAIT-class object and obtains which of the variables are significant given
#' a p-value threshold following a user-defined statistical test.
#' The parameters of the significant features can ve printed to an output table (TRUE by default).
#'
#' @inheritParams MAIT::spectralFUN
#' @return A MAIT-class object containing the significant features of the scores slot of MAIT-class object used as an input.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family statistical functions
#' @keywords internal
#' @noRd
lcms_spectral_fun <- function (pvalue=0.05,
                               p.adj="none",
                               MAIT.object=NULL,
                               printCSVfile=TRUE,
                               test.fun=NULL,
                               namefun=NULL)
{

  if (is.null(test.fun)) {
    stop("No input test function object was given")
  }

  if (is.null(MAIT.object)) {
    stop("No input MAIT object was given")
  }

  if (!is.null(namefun)&class(namefun)!="character"){
    stop("The name of the user-defined test function is not a character string")
  }

  parameters <- list(pvalue,
                     p.adj)
  if(is.null(namefun)){
    names(parameters) <- c("user-defined test p-value",
                           "user-defined test p.adj")
  }else{
    names(parameters) <- c(paste(namefun,"p-value",sep=" "),
                           paste(namefun,"p-value","p.adj",sep=" "))
  }
  MAIT.object@RawData@parameters@sigFeatures <- parameters
  lcms_write_parameter_table(parameters(MAIT.object),folder=MAIT.object@PhenoData@resultsPath)

  data <- MAIT::scores(MAIT.object)
  clases <- MAIT::classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  xsaFA <- MAIT.object@RawData@data
  resultsPath <-MAIT.object@PhenoData@resultsPath
  peakList <- CAMERA::getPeaklist(MAIT.object)

  Fgroups <- matrix(nrow=1)
  Fgroups <- rep(clases,classNum)
  Fgroups <- as.factor(Fgroups[order(Fgroups)])
  peakList <- peakList[order(as.numeric(peakList$pcgroup)),]

  TTs <- matrix(ncol=1,nrow=dim(data)[1])
  Fisher <- matrix(ncol=1,nrow=dim(data)[1])

  Tresults <- matrix(ncol=1,nrow=as.numeric(peakList$pcgroup[dim(peakList)[1]]))
  #	for (i in c(1:dim(data)[1])){
  TTs <- apply(X=data,MARGIN=1,FUN=test.fun,group=Fgroups)

  #	}

  MAIT.object@FeatureData@pvalues <- TTs
  p.corr <- stats::p.adjust(TTs,method=p.adj)

  if (p.adj!="none"){
    index <- which(p.corr<=pvalue)
  }else{
    index<-which(TTs<=pvalue)
  }

  MAIT.object@FeatureData@pvaluesCorrection <- p.adj
  MAIT.object@FeatureData@featureSigID<-index

  return(MAIT.object)
}


#' Boxplots for the significant peak table features
#'
#' It performs boxplots for any of significant features found in a peak table. All the
#' boxplots are stored in a directory /Boxplots/Boxplot_spectra_.
#'
#' @param MAIT.object MAIT object where it is found an annotated peak table.
#' @param treatment_col Treatment for the samples.
#' @return BoxPlots are stored in folders associated to their corresponding spectra. No explicit plot is produced by the device.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family statistical functions
#' @family visualization functions
#' @export
#' @examples
#' file_name_1 <-  system.file("extdata","peak_table_sig_ann.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name_1)
#' file_name_2 <-  system.file("extdata","dataset_pos_rt_rs.rds", package = "NIHSlcms")
#' dataset <-  base::readRDS(file_name_2)
#' treatment_col <- scales::hue_pal()(length(unique(dataset$treatment)))
#' lcms_peak_table_boxplots(peak_table,
#'                          treatment_col = treatment_col)
#'
lcms_peak_table_boxplots <- function (MAIT.object = NULL, treatment_col) {
 treatment <- NULL
 peaks <- NULL

   if (is.null(treatment_col)) {
    stop("No input treatment column was given")
  }
  if (is.null(MAIT.object)) {
    stop("No input MAIT object file was given")
  }
  if (length(MAIT::featureSigID(MAIT.object)) == 0) {
    stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and peakAggregation were launched")
  }
  data <- MAIT::scores(MAIT.object)
  index <- MAIT::featureSigID(MAIT.object)
  class <- MAIT::classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath
  clases <- matrix(nrow = 1)
  for (i in c(1:length(class))) {
    clases <- c(clases, rep(class[i], classNum[i]))
  }
  clases <- clases[-1]
  clases <- as.factor(clases)
  aux <- t(data[index, ])
  if (!file.exists(paste(resultsPath, "Boxplots", sep = "/"))) {
    dir.create(paste(resultsPath, "Boxplots", sep = "/"))
  }else {
    cat(" ", fill = TRUE)
    #cat(paste("Warning: Folder", paste(resultsPath, "Boxplots",
    #                                   sep = "/"), "already exists. Possible file overwritting.",
    #          sep = " "), fill = TRUE)
  }
  for (i in c(1:length(index))) {
    peaks_df <- data.frame(peaks = aux[,i],
                           treatment = clases,
                           stringsAsFactors = FALSE)
    ggplot2::ggplot(peaks_df) +
      ggplot2::geom_boxplot(ggplot2::aes(x = treatment, y = peaks, fill = treatment)) +
      ggplot2::scale_fill_manual("Treatment", values = treatment_col) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::ggtitle(paste("Boxplot of the Spectra", index[i]))
    base::suppressWarnings(
      base::suppressMessages(ggplot2::ggsave(paste(paste(resultsPath, "Boxplots/Boxplot_spectra_",
                                                         sep = "/"), index[i], ".png", sep = ""))))
  }
}

#' Principal Component Analysis (PCA)
#'
#' It performs PCA using on the annotated peak table obtained from a MAIT object.
#'
#' @param MAIT.object MAIT object where it is found an annotated peak table.
#' @param treatment_col Treatment for the samples.
#' @param Log Set to TRUE if the data should be plotted using the logarithm of the intensity.
#' @param center to TRUE if the data should be centered around its mean. See scale.
#' @param scale Set to TRUE if the data should be scaled.
#' @return A MAIT.object. Additionaly: Three different PCA scoreplots are printed in three png files.
#' One using PC1 vs PC2, another with PC1 vs PC3 and the last one with PC2 vs PC3.
#' The files will be stored in the directory /PCA_Scoreplots.
#' @family metabolite identification functions
#' @family dataset_peak_table functions
#' @family statistical functions
#' @family visualization functions
#' @export
#' @examples
#' file_name_1 <-  system.file("extdata","peak_table_sig_ann.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name_1)
#' file_name_2 <-  system.file("extdata","dataset_pos_rt_rs.rds", package = "NIHSlcms")
#' dataset <-  base::readRDS(file_name_2)
#' treatment_col <- scales::hue_pal()(length(unique(dataset$treatment)))
#'
#' lcms_peak_table_pca(peak_table,
#'                    treatment_col = treatment_col,
#'                    Log = FALSE,
#'                    center = TRUE, scale = FALSE)

lcms_peak_table_pca <- function (MAIT.object = NULL,treatment_col, Log = FALSE, center = TRUE, scale = TRUE)
{

  prcomp <- NULL
  pc <- NULL
  PC1 <- NULL
  PC2 <- NULL
  PC3 <- NULL
  treatment <- NULL
  feature <- NULL
  loadings <- NULL
  PC <- NULL

  if (is.null(MAIT.object)) {
    stop("No input MAIT object file was given")
  }
  if (is.null(treatment_col)) {
    stop("No input treatment column was given")
  }
  if (length(MAIT::featureSigID(MAIT.object)) == 0) {
    stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and peakAggregation were launched")
  }
  parameters <- list(Log, center, scale)
  names(parameters) <- c("PCA data logarithm", "PCA data centered",
                         "PCA data scaled")
  MAIT.object@RawData@parameters@plotPCA <- parameters
  lcms_write_parameter_table(parameters(MAIT.object), folder = MAIT.object@PhenoData@resultsPath)
  data <- MAIT::scores(MAIT.object)
  clases <- MAIT::classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  xsaFA <- MAIT.object@RawData@data
  resultsPath <- MAIT.object@PhenoData@resultsPath
  index <- MAIT::featureSigID(MAIT.object)
  cols <- matrix(nrow = 1)
  textCols <- matrix(nrow = 1)
  for (i in 1:length(clases)) {
    cols <- c(cols, rep(i, classNum[i]))
  }
  cols <- as.character(cols[-1])
  textCols <- 1:length(clases)
  if (Log == FALSE) {
    data <- (scale(t(data[index, ]), center = center, scale = scale))
  }else {
    data <- (scale(t(log10(data[index, ] + 1)), center = center,
                   scale = scale))
  }
  if (!file.exists(paste(resultsPath, "PCA", sep = "/"))) {
    dir.create(paste(resultsPath, "PCA", sep = "/"))
  }else {
    cat(" ", fill = TRUE)
    #cat(paste("Warning: Folder", paste(resultsPath, "pca_results",
    #                                   sep = "/"), "already exists. Possible file overwritting.",
    #         sep = " "), fill = TRUE)
  }


  model <- prcomp(data)
  var_expl <-100 *((model$sdev ^ 2)/sum(model$sdev ^ 2))
  model_df <- data.frame(pc = seq(1, dim(model$x)[1], by = 1), model$x,
                         var_expl = var_expl,
                         treatment = clases,
                         stringsAsFactors = FALSE)


  var_plot <- ggplot2::ggplot(model_df) +
    ggplot2::geom_point(ggplot2::aes(x = pc , y = var_expl), color = "blue") +
    ggplot2::scale_x_continuous("Principal Component") +
    ggplot2::scale_y_continuous("Percentage of Variance (%)") +
    ggplot2::ggtitle(paste("Percentage of Variance per Principal Component"))
  print(var_plot)
  base::suppressWarnings(
    base::suppressMessages(ggplot2::ggsave(paste(paste(resultsPath,
                                                       "PCA/Percentage_of_Variance_per_PC.png",
                                                       sep = "/")))))

  PC1_PC2_plot  <- ggplot2::ggplot(model_df) +
    ggplot2::geom_point(ggplot2::aes(x = PC1, y = PC2, color = treatment)) +
    ggplot2::geom_hline(yintercept = 0, linetype="dashed", color = "black") +
    ggplot2::geom_vline(xintercept = 0, linetype="dashed", color = "black") +
    ggplot2::scale_color_manual("Treatment", values = treatment_col) +
    ggplot2::scale_x_continuous(paste0("PC1 (", round(model_df$var_expl[1], 2), " %)")) +
    ggplot2::scale_y_continuous(paste0("PC2 (", round(model_df$var_expl[2], 2), " %)"))  +
    ggplot2::ggtitle(paste("PCA Scoreplot", "PC1", "vs", "PC2"))
  print(PC1_PC2_plot)
  base::suppressWarnings(
    base::suppressMessages(ggplot2::ggsave(paste(paste(resultsPath,"PCA/Scoreplot_PC12.png", sep = "/")),
                                           plot = PC1_PC2_plot)))




  PC1_PC3_plot  <- ggplot2::ggplot(model_df) +
    ggplot2::geom_point(ggplot2::aes(x = PC1, y = PC3, color = treatment)) +
    ggplot2::geom_hline(yintercept = 0, linetype="dashed", color = "black") +
    ggplot2::geom_vline(xintercept = 0, linetype="dashed", color = "black") +
    ggplot2::scale_color_manual("Treatment", values = treatment_col) +
    ggplot2::scale_x_continuous(paste0("PC1 (", round(model_df$var_expl[1], 2), " %)")) +
    ggplot2::scale_y_continuous(paste0("PC3 (", round(model_df$var_expl[3], 2), " %)"))  +
    ggplot2::ggtitle(paste("PCA Scoreplot", "PC1", "vs", "PC3"))
  print(PC1_PC3_plot)
  base::suppressWarnings(
    base::suppressMessages(ggplot2::ggsave(paste(paste(resultsPath,
                                                       "PCA/Scoreplot_PC13.png", sep = "/")))))


  PC2_PC3_plot  <-ggplot2::ggplot(model_df) +
    ggplot2::geom_point(ggplot2::aes(x = PC2, y = PC3, color = treatment)) +
    ggplot2::geom_hline(yintercept = 0, linetype="dashed", color = "black") +
    ggplot2::geom_vline(xintercept = 0, linetype="dashed", color = "black") +
    ggplot2::scale_color_manual("Treatment", values = treatment_col) +
    ggplot2::scale_x_continuous(paste0("PC1 (", round(model_df$var_expl[2], 2), " %)")) +
    ggplot2::scale_y_continuous(paste0("PC2 (", round(model_df$var_expl[3], 2), " %)"))  +
    ggplot2::ggtitle(paste("PCA Scoreplot", "PC2", "vs", "PC3"))
  print(PC2_PC3_plot)
  base::suppressWarnings(
    base::suppressMessages(ggplot2::ggsave(paste(paste(resultsPath,
                                                       "PCA/Scoreplot_PC23.png", sep = "/")))))

  loadings_df <- data.frame(model$rotation[, 1:3],
                            feature = 1:dim(model$rotation[, 1:3])[1]) %>%
    tidyr::gather(key = "PC", value = "loadings", -feature)


  loadings_plot  <- ggplot2::ggplot(loadings_df) +
    ggplot2::geom_line(ggplot2::aes(x = feature, y = loadings, color = PC)) +
    ggplot2::ggtitle(paste("Loadings of the PCA model (3 PCs)")) +
    ggplot2::scale_x_continuous("Signicant Feature Index") +
    ggplot2::scale_y_continuous("Loadings (adim.)")
  base::suppressWarnings(
    base::suppressMessages(ggplot2::ggsave(paste(paste(resultsPath,
                                                       "PCA/Loadings.png", sep = "/")))))
  print(loadings_plot)

  MAIT.object@FeatureData@pcaModel <- list(model)
  return(MAIT.object)
}

