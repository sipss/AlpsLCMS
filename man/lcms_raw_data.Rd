% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metabolite_identification.R
\name{lcms_raw_data}
\alias{lcms_raw_data}
\title{Raw data extractor from a MAIT object}
\usage{
lcms_raw_data(MAIT.object)
}
\arguments{
\item{MAIT.object}{A MAIT-class object.}
}
\value{
A list containing either a xcmsSet or a xsAnnotate object.
}
\description{
Function lcms_raw_data extracts the raw data used to build the MAIT-class object.
}
\examples{
\dontrun{
file_name <-  system.file("extdata", "peak_table_mait.rds", package = "AlpsLCMS")
peak_table <- base::readRDS(file_name)
peak_table_ann <- lcms_peak_annotation(MAIT.object = peak_table)

lcms_raw_data(peak_table_ann)
}
}
\seealso{
Other metabolite identification functions: 
\code{\link{lcms_identify_metabolites}()},
\code{\link{lcms_peak_annotation}()},
\code{\link{lcms_peak_table_boxplots}()},
\code{\link{lcms_peak_table_pca}()},
\code{\link{lcms_sig_peaks_table}()},
\code{\link{lcms_spectral_sig_features}()},
\code{\link{lcms_to_mait}()},
\code{\link{lcms_write_parameter_table}()}

Other dataset_peak_table functions: 
\code{\link{lcms_dataset_load}()},
\code{\link{lcms_identify_metabolites}()},
\code{\link{lcms_peak_annotation}()},
\code{\link{lcms_peak_table_boxplots}()},
\code{\link{lcms_peak_table_pca}()},
\code{\link{lcms_plot_tics}()},
\code{\link{lcms_sig_peaks_table}()},
\code{\link{lcms_spectral_sig_features}()},
\code{\link{lcms_tics}()},
\code{\link{lcms_to_mait}()}

Other import/export functions: 
\code{\link{lcms_convert_ipo_to_xcms}()},
\code{\link{lcms_dataset_load}()},
\code{\link{lcms_dataset_save}()},
\code{\link{lcms_identify_metabolites}()},
\code{\link{lcms_meta_export}()},
\code{\link{lcms_meta_read}()},
\code{\link{lcms_peak_annotation}()},
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
\concept{dataset_peak_table functions}
\concept{import/export functions}
\concept{metabolite identification functions}
