#' filter an experiment by its polarity
#'
#' This function was merged in
#' the MSnbase package on 2018-12-10 as proposed at:
#'
#' https://github.com/lgatto/MSnbase/issues/388
#'
#' Once MSnbase is released (Bioconductor 3.9, maybe once R 3.6 is out),
#' remove this function.
#'
#' @param object An MSnExp object
#' @param polarity. The polarity to keep
#'
#' @export
filterPolarity <- function(object, polarity.) {
  if (missing(polarity.)) return(object)
  polarity. <- as.numeric(polarity.)
  object[MSnbase::polarity(object) %in% polarity.]
}
