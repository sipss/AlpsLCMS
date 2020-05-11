#' HMDB metabolites in positive ionization
#'
#' The function compares (M + H) mass of metabolites contained in the The Human
#' Metabolome Database (HMDB) with features from the positive ionized feature
#' table.
#'
#' @param feature_table A dataframe containing feature in rows and samples in
#'   columns. A column called "mz" is required with the correspoinding mass of each feature.
#'
#' @return The same feature table with an extra column with matched metabolites
#' @export
#'
#' @examples
#' \dontrun{
#' FT <- assignation_pos_HMDB(feature_table)
#' head(FT)
#' }

assignation_pos_HMDB <- function(feature_table){

  HMDB_db <- NULL
  utils::data("HMDB_db", package = "NIHSlcms", envir = environment())

  output_assignation_list <- HMDB_db[NULL,]
  MZs <- rownames(feature_table)
  FT <- feature_table
  FT$HMDB <- NA
  FT$Assignation <- NA

  for (i in 1: nrow(feature_table)) {
    k = feature_table[i, "mz"]
    lower_mz_edge <- k - 0.002
    higher_mz_edge <- k + 0.002

    ind <- intersect(which(HMDB_db$M.plus.H < higher_mz_edge),
                     which(HMDB_db$M.plus.H > lower_mz_edge))
    # assignation_list <- as.data.frame(HMDB_db[ind,])
    if (length(ind) == 0){

      FT[i,"HMDB"] <- NA
      FT[i,"Assignation"] <- NA


    } else {

      FT[i,"HMDB"] <- paste(as.character(HMDB_db[ind, "ENTRY"]), collapse="--- ")
      FT[i,"Assignation"] <- paste(as.character(HMDB_db[ind, "NAME"]), collapse="---")
    }

    output_assignation_list <- rbind(output_assignation_list, HMDB_db[ind,])

  }

  return(FT)
}


#' HMDB metabolites in negative ionization
#'
#' The function compares (M - H) mass of metabolites contained in the The Human
#' Metabolome Database (HMDB) with features from the negative ionized feature
#' table.
#'
#' @param feature_table A dataframe containing feature in rows and samples in
#'   columns. First and second columns need to be the feature ID and mass.
#'
#' @return The same feature table with an extra column with matched metabolites
#' @export
#'
#' @examples
#' \dontrun{
#' FT <- assignation_neg_HMDB(feature_table)
#' head(FT)
#' }

assignation_neg_HMDB <- function(feature_table){

  HMDB_db <- NULL
  utils::data("HMDB_db", package = "NIHSlcms", envir = environment())

  output_assignation_list <- HMDB_db[NULL,]
  MZs <- rownames(feature_table)
  FT <- feature_table
  FT$HMDB <- NA
  FT$Assignation <- NA

  for (i in 1: nrow(feature_table)) {
    k = feature_table[i, 2]
    lower_mz_edge <- k - 0.002
    higher_mz_edge <- k + 0.002

    ind <- intersect(which(HMDB_db$M.minus.H < higher_mz_edge),
                     which(HMDB_db$M.minus.H > lower_mz_edge))
    # assignation_list <- as.data.frame(HMDB_db[ind,])
    if (length(ind) == 0){

      FT[i,"HMDB"] <- NA
      FT[i,"Assignation"] <- NA


    } else {

      FT[i,"HMDB"] <- paste(as.character(HMDB_db[ind, "ENTRY"]), collapse="--- ")
      FT[i,"Assignation"] <- paste(as.character(HMDB_db[ind, "NAME"]), collapse="---")
    }

    output_assignation_list <- rbind(output_assignation_list, HMDB_db[ind,])

  }

  return(FT)
}

#' The Human Metabolome DataBase template: Table with metabolites along with
#' their ionized masses
#'
#' @name HMDB_db
#' @docType data
#' @references \url{hmdb.ca}
#' @keywords data
NULL
