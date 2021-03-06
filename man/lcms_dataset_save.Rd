% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_management.R
\name{lcms_dataset_save}
\alias{lcms_dataset_save}
\title{Save lcms_datasets}
\usage{
lcms_dataset_save(dataset, file_name, ...)
}
\arguments{
\item{dataset}{A lcms_dataset.}

\item{file_name}{The file name to save to.}

\item{...}{Additional arguments passed to \link{saveRDS}.}
}
\description{
Function to save a lcms_dataset object.
}
\examples{
\dontrun{
dataset <- lcms_dataset_load(system.file("extdata", "dataset.rds", package = "AlpsLCMS"))
file_name <- "dataset.rds"
lcms_dataset_save(dataset, file_name)
print(dataset)
}
}
\seealso{
Other dataset functions: 
\code{\link{IPO_group_peaks}()},
\code{\link{find_peaks_cwp}()},
\code{\link{lcms_dataset_load}()},
\code{\link{lcms_fill_chrom_peaks}()},
\code{\link{lcms_filter_mz}()},
\code{\link{lcms_filter_polarity}()},
\code{\link{lcms_filter_rt_min}()},
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

Other import/export functions: 
\code{\link{lcms_convert_ipo_to_xcms}()},
\code{\link{lcms_dataset_load}()},
\code{\link{lcms_identify_metabolites}()},
\code{\link{lcms_meta_export}()},
\code{\link{lcms_meta_read}()},
\code{\link{lcms_peak_annotation}()},
\code{\link{lcms_raw_data}()},
\code{\link{lcms_read_ipo_to_xcms}()},
\code{\link{lcms_read_samples}()},
\code{\link{lcms_rearrange_datafiles_by_class}()},
\code{\link{lcms_sig_peaks_table}()},
\code{\link{lcms_spectral_sig_features}()},
\code{\link{lcms_to_mait}()},
\code{\link{lcms_write_opt_params}()},
\code{\link{lcms_write_parameter_table}()},
\code{\link{phData}()}
}
\concept{dataset functions}
\concept{import/export functions}
