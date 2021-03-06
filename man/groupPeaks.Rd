% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pre_processing.R
\name{groupPeaks}
\alias{groupPeaks}
\title{Peak correspondence}
\usage{
groupPeaks(...)
}
\arguments{
\item{...}{
  Arguments passed on to \code{\link[xcms:groupChromPeaks]{xcms::groupChromPeaks}}
  \describe{
    \item{\code{}}{}
  }}
}
\value{
XCMSnExp object
}
\description{
\code{groupPeaks} is a wrapper of the \code{xcms::groupChromPeaks} from \code{xcms}
package. It performs the \strong{correspondence} (grouping) of chromatographic
peaks of the samples.
}
\seealso{
\code{\link[xcms]{featureDefinitions}} and
\code{\link[xcms]{featureValues,XCMSnExp-method}} for methods to access peak
grouping results.

\code{\link[xcms]{featureChromatograms}} to extract ion chromatograms for each
feature.

\code{\link[xcms]{group}} for the \emph{old} peak grouping methods.

Other peak grouping methods: 
\code{\link[xcms]{groupChromPeaks-density}},
\code{\link[xcms]{groupChromPeaks-mzClust}},
\code{\link[xcms]{groupChromPeaks-nearest}}
}
\author{
Johannes Rainer
}
