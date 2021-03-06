% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_optimization.R
\name{lcms_default_retcorgroup_params}
\alias{lcms_default_retcorgroup_params}
\title{Default parameters for optimization of retention time correction and grouping parameters}
\usage{
lcms_default_retcorgroup_params(
  profStep = c(0.7, 1),
  gapExtend = c(2.1, 2.7),
  optimize = TRUE
)
}
\arguments{
\item{profStep}{set the m/z step for generating profile (matrix) data from raw mass spectral data.
Smaller steps yield more precision at the cost of greater memory usage.}

\item{gapExtend}{numeric(1) defining the penalty for gap enlargement.
The default value for gapExtend depends on the value of distFun, for distFun = "cor"
and distFun = "cor_opt" it is 2.4, for distFun = "cov" 11.7, for distFun = "euc" 1.8
and for distFun = "prd" 7.8.}

\item{optimize}{by default is TRUE. If FALSE, the function does not optimize
the parameters.}
}
\value{
A parameters template for retention time correction and grouping optimization
}
\description{
The function creates default parameters for optimizing
retention time correction and grouping algorithms.
Parameter optimization on the XCMS preprocessing
algorithms is performed using the IPO Package.
Use this function to create the template within \code{lcms_corgroup_optimization} function.
}
\examples{
\dontrun{
default_retcorgroup_params <- lcms_default_retcorgroup_params(optimize = TRUE)
print(default_retcorgroup_params)
}
}
\concept{optimization functions.}
