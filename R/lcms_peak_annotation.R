#' Peak Annotation
#'
#' The `MAIT` Package allows to annotate the peaks of the peak table provided by `XCMS`.
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
#' file_name <-  system.file("extdata","peak_table_mait.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name)
#' peak_table_ann <- lcms_peak_annotation(MAIT.object = peak_table)
#' lcms_raw_data(peak_table_ann)
#'
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

        quiet <- function(x) {
          base::sink(base::tempfile())
          base::on.exit(base::sink())
          base::invisible(base::force(x))
        }


        writeExcelTable <- function(file,file.name){
                   write.csv(file,paste(file.name,".csv",sep=""))
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
        quiet(xsa <- xsAnnotate(rawData(MAIT.object)$xcmsSet))
        quiet(xsaF <- groupFWHM(xsa, perfwhm = perfwhm, sigma = sigma))
        cat("Spectrum build after retention time done", fill = TRUE)
        quiet(xsaFA <- findIsotopes(xsaF))
        cat("Isotope annotation done", fill = TRUE)
        quiet(xsaFA <- groupCorr(xsaFA, cor_eic_th = corrWithSamp, cor_exp_th = corrBetSamp,
                           calcIso = calcIso, calcCiS = calcCiS, calcCaS = calcCaS,
                           pval = pval, graphMethod = graphMethod))
        cat("Spectrum number increased after correlation done",
            fill = TRUE)
        if (annotateAdducts == TRUE) {
          quiet(xsaFA <- findAdducts(xsaFA, rules = adducts, polarity = "positive"))
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
            warning(paste("Warning: Folder", paste(resultsPath,
                                                   "tables", sep = "/"), "already exists. Possible file overwritting."))
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
#' Function rawData extracts the raw data used to build the MAIT-class object
#' @param MAIT.object A MAIT-class object
#' @return A list containing either a xcmsSet or a xsAnnotate object.
#' @export
#' @examples
#' file_name <-  system.file("extdata", "peak_table_mait.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name)
#' peak_table_ann <- lcms_peak_annotation(MAIT.object = peak_table)
#'
#' lcms_raw_data(peak_table_ann)
#'
lcms_raw_data <- function (MAIT.object) {
  MAIT::rawData(MAIT.object)
}



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
#' @examples
#' \dontrun{
#' file_name <-  system.file("extdata", "peak_table_ANN.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name)
#' peak_table@PhenoData@classes <- c("plasma","plasma..ringer...palm")
#' peak_table_sig_ANN <- lcms_spectral_sig_features(MAIT.object = peak_table,
#'                                           pvalue=0.05,
#'                                           p.adj="none",
#'                                           scale=FALSE)
#' print(peak_table_sig_ANN)
#' }
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


  if (length(classes(MAIT.object)) == 1) {
    if (is.na(classes(MAIT.object)) == TRUE) {
      stop("No class information saved in the MAIT object.")
    }
  }
  if (method(MAIT.object) == "") {
    cat("Skipping peak aggregation step...")
    MAIT.object <- lcms_peakAggregation(MAIT.object, scale = scale)
  }
  if (is.null(test.fun)) {
    if (length(classes(MAIT.object)) == 2) {
      if (parametric == TRUE) {
        if (var.equal == TRUE) {
          out <- lcms_spectralTStudent(pvalue = pvalue, p.adj = p.adj,
                                  MAIT.object = MAIT.object, printCSVfile = printCSVfile)
        }
        else {
          out <- lcms_spectralWelch(pvalue = pvalue, p.adj = p.adj,
                               MAIT.object = MAIT.object, printCSVfile = printCSVfile)
        }
      }
      else {
        out <- lcms_spectralWilcox(pvalue = pvalue, p.adj = p.adj,
                              MAIT.object = MAIT.object, printCSVfile = printCSVfile,
                              jitter = jitter, jitter.factor = jitter.factor,
                              jitter.amount = jitter.amount)
      }
    }
    else {
      if (parametric == TRUE) {
        out <- lcms_spectralAnova(pvalue = pvalue, p.adj = p.adj,
                             MAIT.object = MAIT.object, printCSVfile = printCSVfile)
      }
      else {
        out <- lcms_spectralKruskal(pvalue = pvalue, p.adj = p.adj,
                               MAIT.object = MAIT.object, printCSVfile = printCSVfile)
      }
    }
  }
  else {
    out <- lcms_spectralFUN(pvalue = pvalue, p.adj = p.adj, MAIT.object = MAIT.object,
                       printCSVfile = printCSVfile, test.fun = test.fun,
                       namefun = namefun)
  }
  if (length(featureSigID(out)) == 0) {
    warning("No significative features found with the selected parameters.")
  }
  else {
    aux <- lcms_sig_peaks_table(out, printCSVfile = printCSVfile)
  }
  return(out)
}


#' lcms_peakAggregation function applies a peak aggregation technique to the data of a MAIT-class object.
#'
#' @inheritParams MAIT::peakAggregation
#' @return An MAIT-class object.
#' @keywords internal
#' @noRd

lcms_peakAggregation<-function(MAIT.object=NULL,
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
  lcms_writeParameterTable(parameters(MAIT.object),folder=MAIT.object@PhenoData@resultsPath)
  classes <- classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath#resultsPath(MAIT.object)

  MAIT.object@FeatureInfo@peakAgMethod <- method


  aux <- getScoresTable(MAIT.object=MAIT.object,getSpectra=TRUE,getExtendedTable=FALSE)
  data <- aux$scores
  if(scale==TRUE && method!="Single"){

    data<-data/rowMeans(data)
    data[is.nan(as.matrix(data))]<-0

  }
  idGroup <- aux$spectraID

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
    if(!file.exists(paste(resultsPath,"Tables",sep="/"))){

      dir.create(paste(resultsPath,"Tables",sep="/"))
    }else{
      cat(" " ,fill=TRUE)
      #         warning(paste("Folder",paste(resultsPath,"Tables",sep="/"),"already exists. Possible file overwritting.",sep=" "),fill=TRUE)
    }

    tabl <- scores(MAIT.object)
    if(length(rawData(MAIT.object))!=0){
      if(is.null(samples)){
        colnames(tabl) <- sampnames(rawData(MAIT.object)[[1]]@xcmsSet)

      }else{
        colnames(tabl) <- sampnames(rawData(MAIT.object)[[1]]@xcmsSet)[samples]
      }

    }else{
      colnames(tabl) <- colnames(scores(MAIT.object))
    }
    rownames(tabl) <- paste("S",1:dim(tabl)[1],sep="")
    if (scale==TRUE){
      norm <- "scaled"
    }else{
      norm <- "notScaled"
    }
    write.table(file=paste(resultsPath,paste("Tables/dataSet_",norm,".csv",sep=""),sep="/"),x=tabl,row.names=TRUE,col.names=NA,sep=",")
  }

  return(MAIT.object)

}



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
#' file_name <-  system.file("extdata", "peak_table_sig_ANN.rds", package = "NIHSlcms")
#' peak_table <- base::readRDS(file_name)
#' #peak_table@PhenoData@classes <- c("plasma","plasma..ringer...palm")
#' sig_table <- lcms_sig_peaks_table(peak_table,  printCSVfile=FALSE)
#' str(sig_table)
#' }
lcms_sig_peaks_table<-function(
                              MAIT.object=NULL,
                              printCSVfile=FALSE,
                              extendedTable=TRUE,
                              printAnnotation=TRUE){

  if (is.null(MAIT.object)) {
    stop("No MAIT object was given")
  }

  if(length(featureSigID(MAIT.object))==0){
    stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and spectralSigFeatures were launched")
  }

  dat <- getScoresTable(MAIT.object,getSpectra=TRUE,getExtendedTable=TRUE)

  if(printAnnotation==TRUE){extendedTable<-TRUE}

  if(extendedTable==TRUE){
    peakList <- dat$extendedTable
  }else{
    peakList <- dat$scores
  }
  spectraID <- dat$spectraID

  data <- scores(MAIT.object)
  index <- featureSigID(MAIT.object)
  classes <- classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath#resultsPath(MAIT.object)
  Fisher <- LSDResults(MAIT.object)
  TTs <- pvalues(MAIT.object)

  sigPeaksTable <- matrix(nrow=1,ncol=ncol(peakList))
  colnames(sigPeaksTable) <- colnames(peakList)

  if (length(index)!=0){
    if (method(MAIT.object)!="None"){
      sigPeaksTable <- peakList[spectraID%in%index,]
    }else{
      sigPeaksTable <- peakList[spectraID%in%unique(spectraID[index]),]

    }

    p <- matrix(nrow=1,ncol=dim(sigPeaksTable)[1])
    fisher <- matrix(nrow=1,ncol=dim(sigPeaksTable)[1])
    if(class(classes(MAIT.object))!="logical"|class(classNum(MAIT.object))!="logical"){
      if(length(classes(MAIT.object))>2){
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

    p<-pvalues(MAIT.object)
    p <- p[spectraID%in%unique(spectraID[index])]

    p <- matrix(p,ncol=1)
    if(MAIT.object@FeatureData@pvaluesCorrection==""){MAIT.object@FeatureData@pvaluesCorrection<-"none"}
    P.adjust <- p.adjust(p,MAIT.object@FeatureData@pvaluesCorrection)

    pAux <- matrix(ncol=1,nrow=dim(sigPeaksTable)[1])
    pAux.adjust <- matrix(ncol=1,nrow=dim(sigPeaksTable)[1])


    if (method(MAIT.object)!="None"){
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
    if(length(classes(MAIT.object))>2){
      colnames(sigPeaksTable)[dim(sigPeaksTable)[2]] <- NamesFisher
    }else{
      colnames(sigPeaksTable)[dim(sigPeaksTable)[2]] <- "Fisher.Test"
    }
    colnames(sigPeaksTable)[dim(sigPeaksTable)[2]-2] <- "P.adjust"
    colnames(sigPeaksTable)[dim(sigPeaksTable)[2]-1] <- "p"


    if(length(rawData(MAIT.object))==0&is.null(sigPeaksTable$adduct)==TRUE){
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

    if(class(classes(MAIT.object))!="logical"|class(classNum(MAIT.object))!="logical"){
      Fgroups <- as.factor(rep(classes(MAIT.object),classNum(MAIT.object)))
      for(i in c(1:dim(sigPeaksTable)[1])){
        if(length(rawData(MAIT.object))==0){
          ind <- c(3:(dim(sigPeaksTable)[2]-5))
        }else{
          ind <- c((8+length(classes(MAIT.object))):(dim(sigPeaksTable)[2]-6))
        }
        means[[i]] <- aggregate(as.numeric(sigPeaksTable[i,ind])~Fgroups,FUN=mean)
        medians[[i]] <- aggregate(as.numeric(sigPeaksTable[i,ind])~Fgroups,FUN=median)
      }

      allMeans <- merge(means[[1]],means[[2]],by="Fgroups")
      if(length(rawData(MAIT.object))!=0){
        colnames(allMeans)[2:dim(allMeans)[2]] <- rownames(dat$scores[spectraID%in%unique(spectraID[index]),])[c(1,2)]
      }else{
        colnames(allMeans)[2:dim(allMeans)[2]] <- 1:(dim(allMeans)[2]-1)
      }
      if(length(index)>2){
        for(i in c(3:length(means))){
          allMeans <- merge(allMeans,means[[i]],by="Fgroups")
          if(length(rawData(MAIT.object))!=0){
            colnames(allMeans)[i+1] <-  rownames(dat$scores[spectraID%in%unique(spectraID[index]),])[i]
          }else{
            colnames(allMeans)[i+1]  <- as.numeric(colnames(allMeans)[i])+1
          }
        }
      }
      rownames(allMeans)<-paste("Mean Class",allMeans[,1])

      temp<-as.data.frame(t(allMeans[,-1]))
      for(i in c(1:length(classes(MAIT.object)))){
        temp[,i]<-as.numeric(as.character(temp[,i]))
      }

      if(length(rawData(MAIT.object))!=0){rownames(temp)<-NULL}
      allMedians <- merge(medians[[1]],medians[[2]],by="Fgroups")
      if(length(rawData(MAIT.object))!=0){
        colnames(allMedians)[2:dim(allMedians)[2]] <- rownames(dat$scores[spectraID%in%unique(spectraID[index]),])[c(1,2)]
      }else{
        colnames(allMedians)[2:dim(allMedians)[2]] <- 1:(dim(allMedians)[2]-1)
      }
      if(length(index)>3){
        for(i in c(3:length(medians))){
          allMedians <- merge(allMedians,medians[[i]],by="Fgroups")
          if(length(rawData(MAIT.object))!=0){
            colnames(allMedians)[i+1] <-  rownames(dat$scores[spectraID%in%unique(spectraID[index]),])[i]
          }else{
            colnames(allMedians)[i+1]  <- as.numeric(colnames(allMedians)[i])+1

          }
        }
      }
      rownames(allMedians)<-paste("Median Class",allMedians[,1])

      tempMed<-as.data.frame(t(allMedians[,-1]))
      for(i in c(1:length(classes(MAIT.object)))){
        tempMed[,i]<-as.numeric(as.character(tempMed[,i]))
      }
      sigPeaksTable <- cbind(sigPeaksTable,temp,tempMed)
      if(length(rawData(MAIT.object))!=0){rownames(tempMed)<-NULL}
    }else{
      temp <- matrix(rep(NA,dim(sigPeaksTable)[1]),ncol=1)
      tempMed <- matrix(rep(NA,dim(sigPeaksTable)[1]),ncol=1)
      sigPeaksTable <- cbind(sigPeaksTable,temp,tempMed)
      colnames(sigPeaksTable)[c(dim(sigPeaksTable)[2]-1,dim(sigPeaksTable)[2])] <- c("Class Mean", "Class Median")
    }

    if(printCSVfile==TRUE){
      if(!file.exists(paste(resultsPath,"Tables",sep="/"))){
        dir.create(paste(resultsPath,"Tables",sep="/"))
      }else{
        cat(" " ,fill=TRUE)
        warning(paste("Folder",paste(resultsPath,"Tables",sep="/"),"already exists. Possible file overwritting.",sep=" "))
      }
      write.csv(x=sigPeaksTable,file=paste(resultsPath,"Tables/significantFeatures.csv",sep="/"),row.names=FALSE)
    }

  }else{
    stop("There are no significant features for this pvalue")
  }
  return(sigPeaksTable)
}

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
#' @examples
#' \dontrun{
#' file_name <- system.file("extdata", "peak_table_sig_ANN.rds", package = "NIHSlcms")
#' peak_table <- readRDS(file_name)
#' metabololite_table <- lcms_identify_metabolites(MAIT.object = peak_table,
#'                               peakTolerance = 0.005)
#' }
  lcms_identify_metabolites <- function(MAIT.object=NULL,
                                  peakTolerance=0.005,
                                  database=NULL,
                                  polarity="positive",
                                  printCSVfile=TRUE){

    if (is.null(MAIT.object)) {
      stop("No input MAIT object file was given")
    }

    if(length(featureSigID(MAIT.object))==0){
      stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and spectralSignFeatures were launched")
    }
    if(length(MAIT.object@FeatureData@masses)==0 & length(rawData(MAIT.object))==0){
      stop("No peak masses found in the MAIT object")
    }

    if (is.null(database)) {
      identMetEnv<-new.env()
      data(MAITtables,envir=identMetEnv)
      Database<-get("Database",envir=identMetEnv)
      dataBase <- Database
      cat("WARNING: No input database table was given. Selecting default MAIT database...",fill=TRUE)
    }else{
      dataBase<-read.csv(paste(database,".csv",sep=""),sep=",",header=TRUE)
      dataBase <- as.matrix(dataBase[order(dataBase[,1]),])
    }
    parameters <- list(peakTolerance,
                       database,
                       polarity)
    names(parameters) <- c("peakTolerance",
                           "database",
                           "polarity")

    MAIT.object@RawData@parameters@identifyMetabolites <- parameters
    lcms_writeParameterTable(parameters(MAIT.object),folder= MAIT.object@PhenoData@resultsPath)

    signSpectra <-featureSigID(MAIT.object)
    sigPeaksTable <- sigPeaksTable(MAIT.object,printCSVfile=FALSE)
    resultsPath <- MAIT.object@PhenoData@resultsPath

    if(length(rawData(MAIT.object))==0){
      Search <- matrix(nrow=1,ncol=8)
      colnames(Search) <- c("Query Mass","Database Mass (neutral mass)","rt","Adduct","Name","spectra","Biofluid","ENTRY")
    }else{
      Search <- matrix(nrow=1,ncol=9)
      colnames(Search) <- c("Query Mass","Database Mass (neutral mass)","rt","Isotope","Adduct","Name","spectra","Biofluid","ENTRY")
    }
    H.mass <- 1.00794

    temp <- getScoresTable(MAIT.object,getExtendedTable=TRUE)
    peakList <- temp$extendedTable
    spec <- temp$spectraID

    signPeaklist <- peakList[spec%in%signSpectra,]

    aux <- signPeaklist[c(-grep("[M+1]",signPeaklist$isotope,fixed=TRUE),-grep("[M+2]",signPeaklist$isotope,fixed=TRUE)),]
    aux <- aux[-which(aux$isotope==""),]

    if(length(rawData(MAIT.object))==0){

      peaksIsTP <- matrix(nrow=1,ncol=dim(aux)[2]+2*length(classes(MAIT.object)))
      colnames(peaksIsTP) <- c(colnames(aux)[3:(dim(aux)[2]-1)],"p.adj","p",colnames(sigPeaksTable)[dim(sigPeaksTable)[2]-1-length(classes(MAIT.object))*2],colnames(sigPeaksTable)[(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])

    }else{

      peaksIsTP <- matrix(nrow=1,ncol=dim(aux)[2]-7+2*length(classes(MAIT.object)))
      colnames(peaksIsTP) <- c(colnames(aux)[8:(dim(aux)[2]-3)],"p.adj","p",colnames(sigPeaksTable)[dim(sigPeaksTable)[2]-length(classes(MAIT.object))*2],colnames(sigPeaksTable)[(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
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

          if(length(rawData(MAIT.object))==0){

            Search <- rbind(Search,cbind(round(sigPeaksTable$mz[index],5),as.numeric(SearchPeak$SearchCand[k,4]),round(sigPeaksTable$rt[index],2),sigPeaksTable$adduct[index],as.character(SearchPeak$SearchCand[k,2]),spectra,as.character(SearchPeak$SearchCand[k,5]),as.character(SearchPeak$SearchCand[k,1])))

            if(length(classes(MAIT.object))>=2){

              added <- cbind(sigPeaksTable[index,3:(dim(sigPeaksTable)[2]-5-length(classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,grep("Fisher.",colnames(sigPeaksTable))],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
            }else{
              added <- cbind(sigPeaksTable[index,3:(dim(sigPeaksTable)[2]-5-length(classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable$Fisher.Test[index],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])

            }
            colnames(added) <- colnames(peaksIsTP)
            peaksIsTP <- rbind(peaksIsTP,added)

          }else{
            Search <- rbind(Search,cbind(round(sigPeaksTable$mz[index],5),as.numeric(SearchPeak$SearchCand[k,4]),round(sigPeaksTable$rt[index],2),sigPeaksTable$isotopes[index],sigPeaksTable$adduct[index],as.character(SearchPeak$SearchCand[k,2]),spectra,as.character(SearchPeak$SearchCand[k,5]),as.character(SearchPeak$SearchCand[k,1])))
            added <- cbind(sigPeaksTable[index,8:(dim(sigPeaksTable)[2]-6-length(classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,dim(sigPeaksTable)[2]-length(classes(MAIT.object))*2],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
            colnames(added) <- colnames(peaksIsTP)
            peaksIsTP <- rbind(peaksIsTP,added)
          }
        }

      }else{
        ref <- 1
        if(length(rawData(MAIT.object))==0){

          Search <- rbind(Search,cbind(round(sigPeaksTable$mz[index],5),"Unknown",round(sigPeaksTable$rt[index],2),sigPeaksTable$adduct[index],"Unknown",spectra,as.character(SearchPeak$SearchCand[5]),as.character(SearchPeak$SearchCand[1])))
          #			added <- cbind(round(sigPeaksTable[index,3:(dim(sigPeaksTable)[2]-5-length(classes(MAIT.object))*2)],0),sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,dim(sigPeaksTable)[2]-1-length(classes(MAIT.object))*2],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
          added <- cbind(sigPeaksTable[index,3:(dim(sigPeaksTable)[2]-5-length(classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,grep("Fisher.",colnames(sigPeaksTable))],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
          colnames(added) <- colnames(peaksIsTP)
          peaksIsTP <- rbind(peaksIsTP,added)


        }else{

          Search <- rbind(Search,cbind(round(sigPeaksTable$mz[index],5),"Unknown",round(sigPeaksTable$rt[index],2),sigPeaksTable$isotopes[index],sigPeaksTable$adduct[index],"Unknown",spectra,as.character(SearchPeak$SearchCand[5]),as.character(SearchPeak$SearchCand[1])))

          if(length(classes(MAIT.object))>=2){

            added <- cbind(sigPeaksTable[index,8:(dim(sigPeaksTable)[2]-6-length(classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,grep("Fisher.",colnames(sigPeaksTable))],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])

          }else{
            added <- cbind(sigPeaksTable[index,8:(dim(sigPeaksTable)[2]-6-length(classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable$Fisher.Test[index],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
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
    if(length(rawData(MAIT.object))==0){

      ind1 <- 9:(dim(metaboliteTable)[2]-2*length(classes(MAIT.object))-3)
      ind2<-(dim(metaboliteTable)[2]-2*length(classes(MAIT.object))-2):dim(metaboliteTable)[2]
      tem<-metaboliteTable[,ind2]
      metaboliteTable[,(dim(metaboliteTable)[2]+1-length(ind1)):dim(metaboliteTable)[2]] <- metaboliteTable[,ind1]
      colnames(metaboliteTable)[(dim(metaboliteTable)[2]+1-length(ind1)):dim(metaboliteTable)[2]] <- colnames(metaboliteTable)[ind1]
      metaboliteTable[,(dim(metaboliteTable)[2]+1-length(ind1)-length(ind2)):(dim(metaboliteTable)[2]-length(ind1))] <- tem
      colnames(metaboliteTable)[(dim(metaboliteTable)[2]+1-length(ind1)-length(ind2)):(dim(metaboliteTable)[2]-length(ind1))] <- colnames(tem)

    }else{

      ind1 <- 10:(dim(metaboliteTable)[2]-2*length(classes(MAIT.object))-3)
      ind2<-(dim(metaboliteTable)[2]-2*length(classes(MAIT.object))-2):dim(metaboliteTable)[2]
      tem<-metaboliteTable[,ind2]
      metaboliteTable[,(dim(metaboliteTable)[2]+1-length(ind1)):dim(metaboliteTable)[2]] <- metaboliteTable[,ind1]
      colnames(metaboliteTable)[(dim(metaboliteTable)[2]+1-length(ind1)):dim(metaboliteTable)[2]] <- colnames(metaboliteTable)[ind1]
      metaboliteTable[,(dim(metaboliteTable)[2]+1-length(ind1)-length(ind2)):(dim(metaboliteTable)[2]-length(ind1))] <- tem
      colnames(metaboliteTable)[(dim(metaboliteTable)[2]+1-length(ind1)-length(ind2)):(dim(metaboliteTable)[2]-length(ind1))] <- colnames(tem)
    }

    if(printCSVfile==TRUE){

      if(!file.exists(paste(resultsPath,"Tables",sep="/"))){
        dir.create(paste(resultsPath,"Tables",sep="/"))
        write.table(metaboliteTable,paste(paste( MAIT.object@PhenoData@resultsPath,"Tables","metaboliteTable",sep="/"),".csv",sep=""),col.names=NA,row.names=TRUE,sep=",")

      }else{

        cat(" " ,fill=TRUE)
        cat(paste("Warning: Folder",paste(resultsPath,"Tables",sep="/"),"already exists. Possible file overwritting.",sep=" "),fill=TRUE)
        #warning(paste("Folder",paste(resultsPath,"Tables",sep="/"),"already exists. Possible file overwritting.",sep=" "))
        write.table(metaboliteTable,paste(paste( MAIT.object@PhenoData@resultsPath,"Tables","metaboliteTable",sep="/"),".csv",sep=""),col.names=NA,row.names=TRUE,sep=",")

      }
    }

    MAIT.object@FeatureInfo@metaboliteTable <- metaboliteTable
    return(MAIT.object)

  }


  #' Metabolite identification
  #'
  #' Function SearchCand looks up for a peak into a database
  #' @inheritParams MAIT::SearchCand
  #' @return A matrix containing all the possible hits for that peak candidateiteTable
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
