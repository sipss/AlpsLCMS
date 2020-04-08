# library(NIHSlcms)
# library(dplyr)
#
# path <- "../../datasets/dont_sync/20180319_Lumos_Test_Samples/For_Sergio/"
# samples <- fs::dir_ls(path, glob = "*.raw")
# rawconverter <- '../../lcms/rawconverter/'
# # install_RawConverter_Linux(rawconverter)
# plan(multiprocess)
# cosa <- lcms_raw_to_mzxml(samples = samples, rawconverter = rawconverter)
# plan(sequential)
#
# # Load samples
# samples <- normalizePath(fs::dir_ls("../../datasets/dont_sync/20180319_Lumos_Test_Samples/For_Sergio/", glob = "*.mzXML"))
#
# not_dia_samples <- samples[!grepl(pattern = "_dia.mzXML$", x = samples)]
#
# dataset_not_dia <- lcms_read_samples(not_dia_samples, mode = "onDisk")
#
# # Load and tidy metadata:
# basic_meta <- tibble::tibble(sampleNames = basename(not_dia_samples)) %>%
#   dplyr::mutate(sampleNumber = as.numeric(stringr::str_extract(.data$sampleNames, "[0-9]+")),
#                 DIA = grepl("dia", .data$sampleNames))
#
# metadata <- readxl::read_excel(
#   "../../datasets/dont_sync/20180319_Lumos_Test_Samples/20180319_Lumos_Test_Samples_key.xlsx",
#   range = "A3:C21") %>%
#   dplyr::rename(sampleNumber = `sample number`)
#
# full_metadata <- left_join(basic_meta, metadata, by = "sampleNumber")
#
# # Append tidy metadata to dataset:
# dataset_not_dia <- lcms_meta_add(dataset_not_dia, full_metadata, by = "sampleNames")
#
# saveRDS(object = dataset_not_dia, file = "dataset_not_dia.RDS")
#
# dataset_not_dia <- readRDS("dataset_not_dia.RDS")
# # Filter DIA/NonDIA samples
#
#
#
# ## Get the base peak chromatograms. This reads data from the files.
# bpis <- xcms::chromatogram(dataset_not_dia, aggregationFun = "max")
# saveRDS(bpis, "basepeaks.rds")
#
# bpi <- bpis[1,1]
#
# plot(bpi, xlim = c(500, 520))
#
#
#
# # Total Ion Current
# tc_pol <- S4Vectors::split(MSnbase::polarity(dataset_not_dia), f = MSnbase::fromFile(dataset_not_dia))
# tc_rt <- S4Vectors::split(MSnbase::rtime(dataset_not_dia), f = MSnbase::fromFile(dataset_not_dia))
# tc <- S4Vectors::split(MSnbase::tic(dataset_not_dia), f = MSnbase::fromFile(dataset_not_dia))
#
# plot(tc_rt[[1]][tc_pol[[1]] == 0], tc[[1]][tc_pol[[1]] == 0], type = "l")
#
#
# treatments <- unique(dataset_not_dia$treatment)
# treatment_colors <- RColorBrewer::brewer.pal(n = length(treatments), name = "Set1")
# names(treatment_colors) <- treatments
#
# boxplot(tc, col = treatment_colors[dataset_not_dia$treatment],
#         ylab = "intensity", main = "Total ion current")
#
# filterPolarity <- function(object, polarity.) {
#   if (missing(polarity.)) return(object)
#   polarity. <- as.numeric(polarity.)
#   object[MSnbase::polarity(dataset_not_dia) %in% polarity.]
# }
#
#
# dataset_not_dia_pos <- filterPolarity(dataset_not_dia, polarity. = 1)
# dataset_not_dia_neg <- filterPolarity(dataset_not_dia, polarity. = 0)
#
# # Peak detection
#
# ## Define the rt and m/z range of the peak area
# rtr <- c(300, 400)
# mzr <- c(334.9, 335.1)
# ## extract the chromatogram
# chr_raw <- MSnbase::chromatogram(dataset_not_dia_pos, mz = mzr, rt = rtr)
# plot(chr_raw, col = treatment_colors[chr_raw$treatment])
#
# dataset_not_dia %>%
#   filterPolarity(polarity. = 0) %>%
#   MSnbase::filterFile(1) %>%
#   MSnbase::filterRt(rt = c(100, 150)) %>%
#   MSnbase::filterMz(mz = c(100, 300)) %>%
#   MSnbase::plot(type = "XIC")
#
#
# cwp <- xcms::CentWaveParam(peakwidth = c(6, 80), noise = 5000)
# peakdet_pos <- xcms::findChromPeaks(dataset_not_dia_pos, param = cwp)
#
