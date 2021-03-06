% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_management.R
\name{lcms_meta_add}
\alias{lcms_meta_add}
\title{Add metadata}
\usage{
lcms_meta_add(object, metadata, by = "sampleNames")
}
\arguments{
\item{object}{A lcms_dataset.}

\item{metadata}{A data frame to be merged.}

\item{by}{A column present both in \code{metadata} and in \code{Biobase::pData(object)}.}
}
\value{
A lcms_dataset with the added metadata.
}
\description{
Add metadata to MSnExp object.
}
\examples{
\dontrun{
dataset <- lcms_dataset_load(system.file
                                  ("extdata","dataset.rds",
                                  package = "AlpsLCMS"))

metadata <- lcms_meta_read(system.file("extdata",
                                       "metadata.xlsx",
                                       package = "AlpsLCMS"))

dataset_metadata <- lcms_meta_add(dataset,
                               metadata,
                               by = "sampleNames")
print(dataset_metadata)
}
}
\seealso{
Other metadata functions: 
\code{\link{lcms_meta_export}()},
\code{\link{lcms_meta_read}()},
\code{\link{phData}()}

Other dataset functions: 
\code{\link{IPO_group_peaks}()},
\code{\link{find_peaks_cwp}()},
\code{\link{lcms_dataset_load}()},
\code{\link{lcms_dataset_save}()},
\code{\link{lcms_fill_chrom_peaks}()},
\code{\link{lcms_filter_mz}()},
\code{\link{lcms_filter_polarity}()},
\code{\link{lcms_filter_rt_min}()},
\code{\link{lcms_filter_sample_type}()},
\code{\link{lcms_meta_export}()},
\code{\link{lcms_plot_chrom_peak_image}()},
\code{\link{lcms_plot_tics}()},
\code{\link{lcms_rearrange_datafiles_by_class}()},
\code{\link{lcms_retention_time_alignment_plot}()},
\code{\link{lcms_tics}()},
\code{\link{obiAdjust_Rtime}()},
\code{\link{peakpicking_optimization}()},
\code{\link{phData}()}
}
\concept{dataset functions}
\concept{metadata functions}
