#' Known metabolites in the dataset
#'
#' The function compares retention time and mass of metabolites contained in a
#' given template with features in the feature table.
#'
#' @param feature_table_complete A dataframe containing feature in rows and
#'   samples in columns. First, second and third columns need to be the feature
#'   ID, mass and retention time (min). column needs
#' @param metabolites A dataframe with metabolites containing a "MZ" (mass) and
#'   "RT" (retention time in min) columns.
#'
#' @return The same feature table with an extra column with matched metabolites
#' @export
#'
#' @examples
#' FT <- known_metabolites(feature_table_complete, metabolites)
#' head(FT)
#'
known_metabolites <- function(feature_table_complete, metabolites){

  output_assignation_list <- metabolites[NULL,]
  MZs <- rownames(feature_table_complete)
  FT <- feature_table_complete
  FT$Assignation <- NA

  for (i in 1: nrow(feature_table_complete)) {
    k = feature_table_complete[i, 2]
    j = feature_table_complete[i, 3]
    lower_mz_edge <- k - 0.002
    higher_mz_edge <- k + 0.002
    lower_rt_edge <- j - 0.4
    higher_rt_edge <- j + 0.4

    ind <- intersect(which(metabolites$MZ < higher_mz_edge & metabolites$RT < higher_rt_edge),
                     which(metabolites$MZ > lower_mz_edge & metabolites$RT > lower_rt_edge))
    # assignation_list <- as.data.frame(metabolites[ind,])
    if (length(ind) == 0){

      FT[i,"Assignation"] <- NA

    } else {

      FT[i,"Assignation"] <- paste(as.character(metabolites[ind, 1]), collapse="--")
    }

    output_assignation_list <- rbind(output_assignation_list, metabolites[ind,])

  }

  # output_assignation_list <- output_assignation_list[order(output_assignation_list$MZ),]
  return(FT)
}


