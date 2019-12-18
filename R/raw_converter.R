convert_samples <- function(samples, options, workdir = NULL) {
  if (is.null(workdir)) {
    stop("workdir must be the working directory where RawConverter was installed. See ?install_RawConverter")
  }
  rawconverter_exe <- file.path(workdir, "RawConverter_x86", "Release", "RawConverter.exe")
  if (!file.exists(rawconverter_exe)) {
    stop("RawConverter.exe was not found. Please see ?install_RawConverter")
  }

  if (.Platform$OS.type != "windows") {
    wineprefix <- file.path(workdir, "wineprefix")
    setup_wineprefix_env(wineprefix)
    retcode <- furrr::future_map_int(
      samples, function(samp) {
        system2(command = "wine", args = c(rawconverter_exe, samp, options),
                stdout = FALSE)
      },
      .progress = show_progress_bar(),
      .options = furrr::future_options(globals = character(), packages = character())
    )
  } else {
    retcode <- furrr::future_map_int(
      samples, function(samp) {
        system2(command = rawconverter_exe, args = c(samp, options),
                stdout = FALSE)
      },
      .progress = show_progress_bar(),
      .options = furrr::future_options(globals = character(), packages = character())
    )
  }
  if (any(retcode != 0)) {
    number_of_errors <- sum(retcode != 0)
    if (number_of_errors > 10) {
      print_only <- 5
      samp <- glue::glue("- ",
                         glue::glue_collapse(utils::head(samples[retcode != 0], n = print_only),
                                             sep = "\n - ", last = "\n -"),
                         "- ... [ and ", number_of_errors - print_only, " more samples]")
    } else {
      samp <- glue::glue("- ", glue::glue_collapse(samples[retcode != 0], sep = "\n - ", last = "\n -"))
    }
    stop("The conversion of the following samples did not end successfully:\n", samp)
  }
  gsub(pattern = "\\.mzXML$", ".raw", samples)
}

setup_wineprefix_env <- function(wineprefix) {
  if (.Platform$OS.type == "windows") {
    return(NULL)
  }
  Sys.setenv(WINEARCH = "win32", WINEPREFIX = wineprefix, WINEDEBUG = "-all")
  NULL
}

setup_wineprefix <- function(wineprefix, force = FALSE) {
  if (.Platform$OS.type == "windows") {
    return(NULL)
  }
  setup_wineprefix_env(wineprefix)
  if (!dir.exists(wineprefix) || force) {
    if (Sys.info()[["sysname"]] == "Linux" && !capabilities()["X11"]) {
      stop("Wine setup needs an (interactive) X Server and R reports you lack X11 capabilities().\n",
           "Please run this from an X desktop session")
    }
    for (what in c("dotnet40")) {
      message("Installing ", what, " through winetricks in ", wineprefix, ".\nPlease wait and follow the instructions")
      message("Some expected wine errors may appear. Ignore them")
      retcode <- system2("winetricks", what)
      if (retcode != 0) {
        stop("winetricks did not end successfully. ", what, " was not installed properly")
      }
      message(what, " installation completed")
    }
  }
  return(NULL)
}

setup_msfilereader <- function(msfilereader_dir, msfilereader_zip, force = FALSE) {
  if (.Platform$OS.type != "windows") {
    wineprefix <- Sys.getenv("WINEPREFIX")
    XRawfile2_dll <- file.path(wineprefix, "drive_c", "Program Files",
                               "Thermo", "MSFileReader", "XRawfile2.dll")
    if (file.exists(XRawfile2_dll) && !force) {
      message("MSFileReader was already installed")
      return(NULL)
    }
  } else {
    XRawfile2_dll <- file.path("C:", "Program Files",
                               "Thermo", "MSFileReader", "XRawfile2.dll")
    if (file.exists(XRawfile2_dll) && !force) {
      message("MSFileReader was already installed")
      return(NULL)
    }
  }

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
    if (!capabilities()["X11"]) {
      stop("Wine setup needs an (interactive) X Server and R reports you lack X11 capabilities().\n",
           "Please run this from an interactive X desktop session (e.g. install Quartz on OSX)")
    }
    retcode <- system2("wine", msfilereader_exe)
    if (retcode != 0) {
      stop("wine did not end successfully. ", msfilereader_exe, " was not installed properly")
    }
  } else {
    retcode <- system2(msfilereader_exe)
    if (retcode != 0) {
      stop(msfilereader_exe, " did not end successfully.")
    }
  }
  message("MSFileReader installation completed")
  NULL
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
#' On non-windows systems, this function requires wine and winetricks to be installed
#'
#' @param rawconverter Working directory to install RawConverter.
#' @param force logical. Set to `TRUE` to install even if RawConverter already exists.
#' @return the working directory
#' @export
#'
install_RawConverter <- function(rawconverter, force = FALSE) {
  workdir <- fs::path_abs(rawconverter)
  dir.create(workdir, recursive = TRUE)

  wineprefix <- file.path(workdir, "wineprefix")
  setup_wineprefix(wineprefix, force = force)

  msfilereader_zip <- file.path(workdir, "MSFileReader 2.2.62.zip")
  msfilereader_dir <- file.path(workdir, "MSFileReader_2.2.62")
  setup_msfilereader(msfilereader_dir = msfilereader_dir,
                     msfilereader_zip = msfilereader_zip,
                     force = force)

  rawconverterx86_dir <- file.path(workdir, "RawConverter_x86")
  rawconverterx86_zip <- file.path(workdir, "RawConverter_x86.zip")
  setup_RawConverter(rawconverterx86_dir = rawconverterx86_dir,
                     rawconverterx86_zip = rawconverterx86_zip,
                     force = force)
  workdir
}

#' Raw to mzxml
#'
#' Convert RAW samples to mzxml samples.
#'
#' @param samples ThermoFisher .raw files.
#'
#' @param options Optios for the conversion.
#'
#' @param rawconverter The working directory given to [install_RawConverter].
#'
#' @section RawConverter:
#'
#' Raw samples from ThermoFisher are not easily parsed. All implementations go
#' through the proprietary dll from Thermo.
#'
#' This function uses a tool named RawConverter available at:
#' http://fields.scripps.edu/rawconv/
#'
#' RawConverter uses the ThermoFisher tool MSFileReader that takes care of the
#' decryption.
#'
#' Installing that tool may not be straightforward on all OS, so some
#' helpers are provided to the user.
#'
#' On non-Windows systems:
#'   You will need wine and winetricks installed. Then call `install_RawConverter("/your/installation/directory")`
#'
#' @export
#'
lcms_raw_to_mzxml <- function(samples, options = c("--mzxml", "--select_mono_prec"), rawconverter = NULL) {
  samples <- fs::path_abs(samples)
  if (is.null(rawconverter)) {
    stop('rawconverter should be the path to the workdir used at install_RawConverter(rawconverter = "/any/empty/directory")')
  }
  rawconverter <- fs::path_abs(rawconverter)
  convert_samples(samples, options, workdir = rawconverter)
}

#' RAW converter
#'
#' The function lists and converts samples in ".raw" format to ".mzXML" format.
#' It is the first step to create a `lcms_dataset`. It requires the previous
#' download and installation of the *RawConverter* application.
#'
#' @param sample_path Directory in which the samples are.
#' @param file_format Format of the LC-MS files (e.g. file_format = "raw").
#' @param rawconverter_path Directory in where the RawConverter application is located.
#' @return A list of the LC-MS files in a readable format to create the `lcms_dataset`.
#' @export
#' @examples
#' sample_path <- system.file("extdata", package = "NIHSlcms")
#' samples_mzxml <- lcms_list_mzxml_samples(sample_path,
#'                                     file_format = "mzXML")
#' print(samples_mzxml)
lcms_list_mzxml_samples <- function(sample_path,
                                    file_format = "mzXML",
                                    rawconverter_path = NULL){

  if (file_format == "raw") {
    #install_RawConverter(rawconverter)
    samples_raw <- fs::dir_ls(sample_path, glob = "*.raw")
    future::plan(multiprocess)
    lcms_raw_to_mzxml(samples = samples_raw, rawconverter = rawconverter_path)
    future::plan(sequential)
    sample_names_mzxml <- fs::path_abs(fs::dir_ls(sample_path, glob = "*.mzXML"))
  } else if(file_format == "mzXML") {
    sample_names_mzxml<- fs::path_abs(fs::dir_ls(sample_path, glob = "*.mzXML"))
  } else {
    stop("Not allowed file format. Use only either *.raw or *.mzMXML files")
  }
  sample_names_mzxml
}


