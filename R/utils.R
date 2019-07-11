show_progress_bar <- function(...) {
  all(...) && interactive() && is.null(getOption("knitr.in.progress"))
}
