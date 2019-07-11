#' Read LCMS samples on mzXML format
#'
#' @inheritDotParams MSnbase::readMSData
#'
#' @inherit MSnbase::readMSData return
#' @export
#'
lcms_read_samples <- MSnbase::readMSData


#' Title
#'
#' @param msnexp
#' @param metadata
#' @param by
#'
#' @return
#' @export
lcms_meta_add <- function(msnexp, metadata, by = "sampleNames") {
  phenotype_data <- Biobase::pData(msnexp)
  phenotype_data$sampleNames <- as.character(phenotype_data$sampleNames)
  phenotype_data_extra <- dplyr::left_join(phenotype_data, metadata, by = by)
  Biobase::pData(msnexp) <- phenotype_data_extra
  msnexp
}

convert_samples_linux <- function(samples, options, workdir = NULL) {
  if (is.null(workdir)) {
    stop("workdir must be the working directory where RawConverter was installed. See ?install_RawConverter_Linux")
  }
  wineprefix <- file.path(workdir, "wineprefix")
  rawconverterx86_dir <- file.path(workdir, "RawConverter_x86")
  rawconverter_exe <- file.path(rawconverterx86_dir, "Release", "RawConverter.exe")
  Sys.setenv(WINEARCH = "win32", WINEPREFIX = wineprefix, WINEDEBUG = "-all")
  furrr::future_map(
    samples, function(samp) {
      system2("wine", args = c(rawconverter_exe, samp, options),
              stdout = FALSE)
    },
    .progress = show_progress_bar(),
    .options = furrr::future_options(globals = character(), packages = character())
  )
}

convert_samples_windows <- function(samples, options, raw_converter_exe) {
  furrr::future_map(
    samples, function(samp) {
      system2(command = raw_converter_exe, args = c(samp, options))
    }
  )
}

setup_wineprefix <- function(wineprefix, force = FALSE) {
  Sys.setenv(WINEARCH = "win32", WINEPREFIX = wineprefix, WINEDEBUG = "-all")
  if (!dir.exists(wineprefix) || force) {
    for (what in c("dotnet40")) {
      message("Installing ", what, " through winetricks in ", wineprefix, ".\nPlease wait and follow the instructions")
      message("Some expected wine errors may appear. Ignore them")
      system2("winetricks", what)
      message(what, " installation completed")
    }
  }
}

setup_msfilereader <- function(msfilereader_dir, msfilereader_zip, force = FALSE) {
  if (!file.exists(msfilereader_zip)) {
    message("Downloading MSFileReader 2.2.262")
    utils::download.file("http://fields.scripps.edu/rawconv/download/MSFileReader%202.2.62.zip",
                         destfile = msfilereader_zip)
  }
  message("Unzipping MSFileReader_2.2.62")
  utils::unzip(msfilereader_zip, exdir = msfilereader_dir)
  message("MSFileReader_2.2.62 unzipping completed")

  msfilereader_exe <- file.path(msfilereader_dir, "MSFileReader.exe")
  message("Installing MSFileReader_2.2.62, please follow installer instructions...")
  if (.Platform$OS.type != "windows") {
    system2("wine", msfilereader_exe)
  } else {
    system2(msfilereader_exe)
  }
  message("MSFileReader installation completed")
}

setup_RawConverter <- function(rawconverterx86_dir, rawconverterx86_zip, force = FALSE) {
  rawconverter_exe <- file.path(rawconverterx86_dir, "Release", "RawConverter.exe")

  if (file.exists(rawconverter_exe) && !force) {
    message("RawConverter is installed")
    return(invisible(NULL))
  }

  if (!file.exists(rawconverterx86_zip)) {
    message("Downloading RawConverter_x86.zip")
    utils::download.file("http://fields.scripps.edu/rawconv/download/RawConverter_x86.zip")
  }

  message("Unzipping RawConverter_x86")
  utils::unzip(rawconverterx86_zip, exdir = rawconverterx86_dir)
  message("RawConverter is installed")
  return(invisible(NULL))
}

#' Download and Install RawConverter
#'
#' This function requires wine and winetricks to be installed
#'
#' @param workdir Working directory to install RawConverter
#' @param force logical. Set to `TRUE` to install even if RawConverter already exists
#'
#' @return the working directory
#' @export
#'
install_RawConverter_Linux <- function(workdir, force = FALSE) {
  workdir <- normalizePath(workdir)
  dir.create(workdir, recursive = TRUE)

  wineprefix <- file.path(workdir, "wineprefix")
  setup_wineprefix(wineprefix)

  msfilereader_zip <- file.path(workdir, "MSFileReader 2.2.62.zip")
  msfilereader_dir <- file.path(workdir, "MSFileReader_2.2.62")
  setup_msfilereader(msfilereader_dir = msfilereader_dir,
                     msfilereader_zip = msfilereader_zip)

  rawconverterx86_dir <- file.path(workdir, "RawConverter_x86")
  rawconverterx86_zip <- file.path(workdir, "RawConverter_x86.zip")
  setup_RawConverter(rawconverterx86_dir = rawconverterx86_dir,
                     rawconverterx86_zip = rawconverterx86_zip)
  workdir
}


#' Convert RAW samples to mzxml samples
#'
#' @param samples ThermoFisher raw files
#'
#' @param rawconverter On Linux, the workdir for `install_RawConverter_Linux`. On
#'  Windows, the path to RawConverter. See the RawConverter section.
#'
#' @section RawConverter:
#' Raw samples from ThermoFisher are digitally encrypted and cannot be read by
#' any software, allegedly to prevent sample forgery in pharmaceutical applications.
#'
#' This function uses a tool named RawConverter available at:
#' http://fields.scripps.edu/rawconv/
#'
#' RawConverter uses the ThermoFisher tool MSFileReader that takes care of the
#' decryption.
#'
#' Installing that tool may not be straightforward on all OS, so we provide some
#' helpers.
#'
#' On Linux systems:
#'   You will need wine and winetricks installed. Then call `install_RawConverter_Linux("/your/installation/directory")`
#' On Windows systems:
#'   You will need to install RawConverter yourself.
#'
#' On OSX, the Linux solution may work, assuming you can install wine and winetricks.
#'
#'
#' @export
#'
lcms_raw_to_mzxml <- function(samples, options = c("--mzxml", "--select_mono_prec"), rawconverter = NULL) {
  samples <- normalizePath(samples)
  if (.Platform$OS.type == "windows") {
    if (is.null(rawconverter)) {
      stop("rawconverter should be the path to RawConverter.exe")
    }
    convert_samples_windows(samples = samples, options = options, raw_converter_exe = rawconverter)
  } else if (.Platform$OS.type == "unix") {
    if (is.null(rawconverter)) {
      stop('rawconverter should be the path to the workdir used at install_RawConverter_Linux(workdir = "/any/empty/directory")')
    }
    rawconverter <- normalizePath(rawconverter)
    convert_samples_linux(samples, options, workdir = rawconverter)
  } else {
    stop("Sorry, lcms_raw_to_mzxml is not supported on your OS. See ?lcms_raw_to_mzxml for details")
  }
}
