% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pre_processing.R
\name{lcms_filter_rt_min}
\alias{lcms_filter_rt_min}
\title{Filter by retention time}
\usage{
lcms_filter_rt_min(dataset, rt = c(4, 14))
}
\arguments{
\item{dataset}{A lcms_dataset filtered by retention time.}

\item{rt}{Range of the retention time to keep in minutes.}
}
\value{
A filtered lcms_dataset with the selected retention time range.
}
\description{
This function converts seconds into minutes to cut and keep
a range of the retention time in minutes.
}
\examples{
\dontrun{
rtime_range = c(5,10)
dataset <- lcms_dataset_load(system.file("extdata",
                                         "dataset_pos.rds",
                                         package = "AlpsLCMS"))
dataset_pos_rt <-lcms_filter_rt_min(dataset, rt = rtime_range)

print(dataset_pos_rt)
}
}
\seealso{
Other dataset functions: 
\code{\link{IPO_group_peaks}()},
\code{\link{find_peaks_cwp}()},
\code{\link{lcms_dataset_load}()},
\code{\link{lcms_dataset_save}()},
\code{\link{lcms_fill_chrom_peaks}()},
\code{\link{lcms_filter_mz}()},
\code{\link{lcms_filter_polarity}()},
\code{\link{lcms_filter_sample_type}()},
\code{\link{lcms_meta_add}()},
\code{\link{lcms_meta_export}()},
\code{\link{lcms_plot_chrom_peak_image}()},
\code{\link{lcms_plot_tics}()},
\code{\link{lcms_rearrange_datafiles_by_class}()},
\code{\link{lcms_retention_time_alignment_plot}()},
\code{\link{lcms_tics}()},
\code{\link{obiAdjust_Rtime}()},
\code{\link{peakpicking_optimization}()},
\code{\link{phData}()}

Other filtering functions: 
\code{\link{lcms_filter_mz}()},
\code{\link{lcms_filter_polarity}()},
\code{\link{lcms_filter_sample_type}()}
}
\concept{dataset functions}
\concept{filtering functions}
