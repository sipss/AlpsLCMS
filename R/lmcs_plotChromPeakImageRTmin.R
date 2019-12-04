#' Image of Chromatographic Peaks by sample
#'
#' It plots the an image of the chromatographic peaks for each sample. This function is useful if
#' you are interested in knowking the effect of the retention time correction on the chromatographic axis.
#'
#' @param lcms_dataset An lcms_dataset
#' @examples
#' \dontrun{
#' file_name <-  system.file("extdata", "peakdet_align.rds", package = "NIHSlcms")
#' lcms_dataset <- base::readRDS(file_name)
#' chr_peak_image <- lmcs_plotChromPeakImageRTmin(lcms_dataset, binSize = 5, xlim = NULL, log = FALSE,
#'                                                xlab = "retention time (min)", yaxt = par("yaxt"))
#'                                                title(main ="Detected Peaks (Aligned)")
#' }
#' @return An image with the detected chromatographic peak, for each sample
#' @export



lmcs_plotChromPeakImageRTmin <- function (x, binSize = 30, xlim = NULL, log = FALSE, xlab = "retention time",
                                     yaxt = par("yaxt"), main = "Chromatographic peak counts",
                                     ...)
{
  min2sec <- 60
  if (!is(x, "XCMSnExp"))
    stop("'x' is supposed to be an 'XCMSnExp' object, but I got a ",
         class(x))
  if (is.null(xlim))
    xlim <- c(floor(min(rtime(x))), ceiling(max(rtime(x))))
  brks <- seq(xlim[1], xlim[2], by = binSize)
  if (brks[length(brks)] < xlim[2])
    brks <- c(brks, brks[length(brks)] + binSize)
  pks <- xcms::chromPeaks(x, rt = xlim)
  if (nrow(pks)) {
    rts <- split(pks[, "rt"], pks[, "sample"])
    cnts <- lapply(rts, function(z) {
      hst <- hist(z, breaks = brks, plot = FALSE)
      hst$counts
    })
    n_samples <- length(fileNames(x))
    sample_idxs <- 1:n_samples
    sample_idxs <- sample_idxs[!(as.character(sample_idxs) %in%
                                   names(rts))]
    if (length(sample_idxs)) {
      all_cnts <- vector("list", n_samples)
      all_cnts[as.numeric(names(cnts))] <- cnts
      zeros <- rep(0, (length(brks) - 1))
      all_cnts[sample_idxs] <- list(zeros)
      cnts <- all_cnts
    }
    cnts <- t(do.call(rbind, cnts))
    if (log)
      cnts <- log2(cnts)
    image(z = cnts, x = (brks - (brks[2] - brks[1])/2) / min2sec, xaxs = "r",
          xlab = xlab, yaxt = "n", ...)
      sample_labels <- stringr::str_extract(basename(fileNames(x)), "\\w+")
      axis(side = 2, at = seq(0, 1, length.out = n_samples),
           labels = FALSE)
      text(y = seq(0, 1, length.out = n_samples), par("usr")[1],
           cex = 0.6, labels = sample_labels,
           srt = 60, pos = 2, xpd = TRUE)
  }
}
