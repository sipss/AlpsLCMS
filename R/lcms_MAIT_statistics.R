#' Extract Significant Features From A MAIT Object For Two Classes
#'
#' Function spectralTStudent takes a MAIT-class object and obtains which of the variables are significant
#' given a p-value threshold when there only are two classes in the raw data.
#' The parameters of the significant features can be printed to an output table (TRUE by default).
#' @inheritParams MAIT::spectralTStudent
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

  data <- scores(MAIT.object)
  clases <- classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  #xsaFA <- rawData(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath#resultsPath(MAIT.object)

  auxs<-getScoresTable(MAIT.object=MAIT.object,getExtendedTable=TRUE)
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
    model <- lm(lmdata~group)
    if (length(which(diff(data)!=0))){
      ttest <- t.test(as.vector(t(data[i,]))~group,var.equal=TRUE)
      TTs[i] <- ttest$p.value
      Tresults[i] <- ttest$statistic
    }else{
      TTs[i] <- 1
    }
  }
  MAIT.object@FeatureData@pvalues <- TTs
  p.corr <- p.adjust(TTs,method=p.adj)
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
#' Function spectralWelch takes an MAIT-class object and obtains which of the variables are significant given a p-value
#' threshold following a Welch test. The parameters of the significant features can ve printed to an output table (TRUE by default).
#' @inheritParams MAIT::spectralWelch
#' @return A MAIT-class object containing the significant features of the scores slot of MAIT-class object used as an input.
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


  data <- scores(MAIT.object)
  clases <- classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  #xsaFA <- rawData(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath
  #	peakList <- getPeaklist(MAIT.object)
  #	peakList <- peakList[order(as.numeric(peakList$pcgroup)),]

  auxs<-getScoresTable(MAIT.object=MAIT.object,getExtendedTable=TRUE)
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
    model <- lm(lmdata~group)
    if (length(which(diff(data)!=0))){
      ttest <- t.test(as.vector(t(data[i,]))~group,var.equal=FALSE)
      TTs[i] <- ttest$p.value
      Tresults[i] <- ttest$statistic
    }else{
      TTs[i] <- 1
    }
  }

  MAIT.object@FeatureData@pvalues <- TTs
  p.corr <- p.adjust(TTs,method=p.adj)

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
#' Function spectralWilcox takes an MAIT-class object and obtains which of the variables are significant
#' given a p-value threshold following a Mann-Witney-Wilcoxon test.
#' The parameters of the significant features can ve printed to an output table (TRUE by default).
#' @inheritParams MAIT::spectralWilcox
#' @return A MAIT-class object containing the significant features of the scores slot of MAIT-class object used as an input.
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

  data <- scores(MAIT.object)
  clases <- classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  #xsaFA <- rawData(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath#resultsPath(MAIT.object)


  auxs<-getScoresTable(MAIT.object=MAIT.object,getExtendedTable=TRUE)
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
    model <- lm(lmdata~group)
    if (length(which(diff(data)!=0))){
      if(jitter==TRUE){
        test<-wilcox.test(jitter(data[i,],factor = jitter.factor, amount = jitter.amount)~group)
      }else{
        test<-wilcox.test(data[i,]~group)
      }
      TTs[i] <- test$p.value
    }else{
      TTs[i] <- 1
    }
  }
  MAIT.object@FeatureData@pvalues <- TTs
  p.corr <- p.adjust(TTs,method=p.adj)

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
#' Function spectralAnova takes an MAIT-class object and obtains which of the variables are significant given a p-value threshold.
#' The parameters of the significant features can ve printed to an output table (TRUE by default).
#' @inheritParams MAIT::spectralAnova
#' @return A MAIT-class object containing the significant features of the scores slot of MAIT-class object used as an input.
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

  data <- scores(MAIT.object)
  clases <- classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath#resultsPath(MAIT.object)
  auxs<-getScoresTable(MAIT.object=MAIT.object,getExtendedTable=TRUE)
  peakList <- auxs$extendedTable

  Fgroups <- matrix(nrow=1)
  Fgroups <- rep(clases,classNum)
  Fgroups <- as.factor(Fgroups[order(Fgroups)])
  #        peakList <- peakList[order(as.numeric(peakList$pcgroup)),]

  TTs <- matrix(ncol=1,nrow=dim(data)[1])
  Fisher <- matrix(ncol=1,nrow=dim(data)[1])

  Tresults <- matrix(ncol=1,nrow=as.numeric(peakList$pcgroup[dim(peakList)[1]]))

  for (i in c(1:dim(data)[1])){
    lmdata <- as.vector(t(data[i,]))
    numbers <- classNum
    names(numbers) <- clases
    model <- lm(lmdata~Fgroups)
    an <- anova(model)
    Fisher[i] <- FisherLSD(data=data,DFerror=an$Df[2],MSerror=an$Mean[2],index=i,classes=Fgroups,numClasses=length(classNum))[[2]]
    TTs[i] <- an$Pr[1]
  }
  MAIT.object@FeatureData@pvalues <- TTs
  MAIT.object@FeatureData@LSDResults <- Fisher

  p.corr <- p.adjust(TTs,method=p.adj)

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
#' Function spectralKruskal takes an MAIT-class object and obtains which of the variables are significant given a
#' p-value threshold following a Kruskal-Wallis test.
#' The parameters of the significant features can ve printed to an output table (TRUE by default).
#' @inheritParams MAIT::spectralKruskal
#' @return A MAIT-class object containing the significant features of the scores slot of MAIT-class object used as an input.
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

  data <- scores(MAIT.object)
  clases <- classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  resultsPath <- MAIT.object@PhenoData@resultsPath #resultsPath(MAIT.object)
  auxs<-getScoresTable(MAIT.object=MAIT.object,getExtendedTable=TRUE)
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
    an <- kruskal.test(x=lmdata,g=Fgroups)
    TTs[i] <- an$p.value
  }

  MAIT.object@FeatureData@pvalues <- TTs
  #    MAIT.object@FeatureData@LSDResults <- Fisher
  p.corr <- p.adjust(TTs,method=p.adj)

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
#' Function spectralFUN takes an MAIT-class object and obtains which of the variables are significant given
#' a p-value threshold following a user-defined statistical test.
#' The parameters of the significant features can ve printed to an output table (TRUE by default).
#'
#' @inheritParams MAIT::spectralFUN
#' @return A MAIT-class object containing the significant features of the scores slot of MAIT-class object used as an input.
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

  data <- scores(MAIT.object)
  clases <- classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  xsaFA <- MAIT.object@RawData@data
  resultsPath <-MAIT.object@PhenoData@resultsPath
  peakList <- getPeaklist(MAIT.object)

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
  p.corr <- p.adjust(TTs,method=p.adj)

  if (p.adj!="none"){
    index <- which(p.corr<=pvalue)
  }else{
    index<-which(TTs<=pvalue)
  }

  MAIT.object@FeatureData@pvaluesCorrection <- p.adj
  MAIT.object@FeatureData@featureSigID<-index

  return(MAIT.object)
}
