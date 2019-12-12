#' Metabolite Candidates
#'
#' Post processing after obtaining the list of potential metabolite candidates.
#' These metabolites can be filtered by biofluid and are sorted by p-value in ascending order.
#' Unknown metabolites can be removed.
#' After that, the results the results are plotted using a word cloud. The lower the p-value, the higher the font size
#' of the metabolite.
#' Additionally we plot the a bar plot corresponding to the number of metabolites candidates associated to a single compound.
#'
#' @param path_result The folder where the .csv file where the annotated table of metabolites is stored
#' @param biofluid_type Filters the results by biofludid type.
#' @param significance Filters the results by significance.
#' @param metabolite_rm Logical. Remove "Unknown" metabolites.
#' @param  max_words  Maximum number of words allowed in the word cloud.
#' @return Two plots: A word cloud with the potential identified metabolites
#' and a bar plot with the number of potential cantidates by detected compound.
#' @family data visualization functions
#' @export
#' @examples
#' \dontrun{
#' path_result <-  system.file("extdata","results_project","tables","metaboliteTable.csv", package = "NIHSlcms")
#' plots <- lcms_plot_metabolites(path_result, biofluid_type = "Any",
#'                                    significance = 0.05,
#'                                    metabolite_rm = FALSE,
#'                                    max_words = 250)
#' plots$cloud
#'
#' plots$freq
#' }

lcms_plot_metabolites <- function(path_result, biofluid_type = "Any",
                                      significance = 0.05,
                                      metabolite_rm = TRUE,
                                      max_words = 250){

#Load Measured Data
suppressWarnings(Measured <- readr::read_csv(file = path_result,
                    col_names = TRUE,
                    col_types = readr::cols_only(spectra = readr::col_number(),
                                          rt = readr::col_double(),
                                          Name = readr::col_factor(),
                                          Biofluid = readr::col_character(),
                                          p.adj = readr::col_double()),
                    skip_empty_rows = TRUE))
attr(Measured, "spec") <- NULL

# Some Checks

if (!is.logical(metabolite_rm)){
  stop("The variable Metabolite_rm must be of class logical")
}


Biofluid_names <- Measured %>%
                    purrr::pluck("Biofluid") %>%
                    stringr::str_split( ";")  %>%
                    base::unlist() %>%
                    stringr::str_replace("^ ", "") %>%
                    base::unique() %>%
                    base::sort()


if (!is.character(biofluid_type)){
  stop("The variable Biofluid_type must be of class character")
} else if (!(biofluid_type %in% Biofluid_names)){
    if (biofluid_type != "Any"){
      stop(paste0("Your selected Biofluid is not among the following:",
                  "\n",str_c(Biofluid_names, collapse = ", "), "."))
    }

}


if (!is.null(significance)){
  if(is.numeric(significance)){
    if((significance > 1) | (significance < 0)){
      stop("Significance must be set to a value in the range [0, 1]")
    }
  }else {
    stop("Significance must be a variable with class either numeric or NULL")
  }
}

# Modify Column_names and drop na

Measured <- Measured %>%
              dplyr::rename(Spectra = spectra,
              Retention_Time = rt,
              Metabolite = Name,
              Pvalue_Adj = p.adj) %>%
              tidyr::drop_na()


# Remove unknown Metabolites if needed

if(metabolite_rm == TRUE){
  Measured <- Measured %>%
    dplyr::filter(Metabolite != "Unknown")
}


# Filter by Biofluid

if(biofluid_type == "Any"){
  Measured <- Measured %>%
                dplyr::mutate(Selected_Biofluid = "Any")

  max_area_wordcloud <- 15#6

} else{
  Measured <- Measured %>%
                dplyr::filter_at(vars(Biofluid), all_vars(str_detect(Biofluid, biofluid_type))) %>%
                dplyr::mutate(Selected_Biofluid = biofluid_type)

  max_area_wordcloud <- 20#11
}

# Select and arrange the data

Measured <- Measured %>%
              dplyr::select(Retention_Time, Spectra,  Metabolite, Selected_Biofluid, Biofluid, dplyr::everything()) %>%
              dplyr::arrange(Retention_Time, Spectra)

#Filter by p-value

if(!is.null(significance)){
  Measured <- Measured %>%
              dplyr::filter(Pvalue_Adj <= significance)
}


# plots
  #word_cloud
plot_df <- Measured %>%
  dplyr::select(Retention_Time, Spectra, Metabolite, Pvalue_Adj) %>%
  dplyr::mutate(Spectra = as.factor(Spectra)) %>%
  dplyr::mutate(Angle = 0 * sample(c(0, 1),#90
                             dplyr::n(), replace = TRUE,
                             prob = c(60, 40))) %>%
  dplyr::top_n(max_words,-log10(Pvalue_Adj))


#
text_title_cloud <- paste0("Metabolite Word Cloud. Biofluid: ",
                    base::unique(Measured$Selected_Biofluid), ". ","Pvalue < 0.05.")

set.seed(42)
cloud <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    label = Metabolite, size =  -log10(Pvalue_Adj),
    colour = Retention_Time,
    angle = Angle
  )
) +
ggwordcloud::geom_text_wordcloud_area() +
  ggplot2::scale_size_area(max_size = max_area_wordcloud) +#8 guay export, 4 plot
  ggplot2::theme_minimal() +
  ggplot2::scale_colour_continuous(type = "viridis") +
  ggplot2::ggtitle(text_title_cloud)

# plots
  #Barplot:  Spectra - Metabolite Candidate

freq <- plot_df %>%
  dplyr::group_by(Spectra) %>%
  dplyr::count() %>%
  ggplot2::ggplot(ggplot2::aes(x = stats::reorder(Spectra, n),
                          y = n, fill = stats::reorder(Spectra, n),
                          color = stats::reorder(Spectra, n))) +
  ggplot2::geom_bar(stat="identity", position = "stack", alpha = 0.7, size = 1) +
  ggplot2::scale_x_discrete("Spectrum Number") +
  ggplot2::scale_y_continuous("Number of Metabolite Candidates") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 17),
        axis.text.x = ggplot2::element_text(angle = 90),
        axis.text = ggplot2::element_text(size = 8, color = "black"),
        axis.title = ggplot2::element_text(size = 14, color = "black"),
        axis.line = ggplot2::element_line(color = "black",
                                 size = 1, linetype = "solid"),
        legend.position = "none") +
  ggplot2::ggtitle("Metabolite Candidates for Spectrum")

plot_list <- list(cloud = cloud, freq = freq)
plot_list

}


