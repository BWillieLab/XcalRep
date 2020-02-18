#' Get Unique Features
#'
#' Get unique features in input dataset
#'
#'
#' @param df HR-pQCT dataset (data.frame). Following column names are expected: value, site, section, timePoint, phantom, parameter, scanDate
#' @name getUniqueFeatures
#' @return list of unique.features (unique.features), number of variables per unique feature (N), and names of variables (variable)
#'
getUniqueFeatures <- function(df){

  # prep data for assay class
  # accepted.variables = c("value", "site", "section", "timePoint", "phantom", "parameter", "scanDate")
  available.variables <- colnames(df)

  unique.features <- list()
  N <- list()
  N[["n.entries"]] <- nrow(df)
  for (i in 1:length(available.variables)){
    if (available.variables[i] != "value"){
      n.var.name <- paste("n.", available.variables[i], sep = "")
        unique.features[[available.variables[i]]] <- unique(df[, available.variables[i]])
        N[[n.var.name]] <- length(unique.features[[available.variables[i]]])
    }
  }

  uf.output <- list(
    unique.features = unique.features,
    N = N,
    variables= available.variables
  )
  return(uf.output)
}


#' Delete Assay
#'
#' Delete assay from calibration object
#'
#' If specified assay was current assay, new current assay is assigned.
#'
#' @param object Calibration Object
#' @param which.assay A character specifying which assay to delete.
#' @param verbose Default is TRUE. Logical specifying whether to print progress.
#' @name deleteAssay
#' @return calibration object with specified assay deleted.
#'
deleteAssay <- function(object, which.assay, verbose = TRUE){

  existing.assays <- getAssay(object, which.assays = "all")

  if (!(which.assay %in% existing.assays)){
    stop(paste("'", which.assay, "' does not exist. No assays were deleted.", sep = ""))
  }

  # delete assay
  object@assays <- object@assays[!(existing.assays %in% which.assay)]

  # if deleted assay was default, set new default assay
  current.assay <- getAssay(object)
  existing.assays <- getAssay(object, which.assays = "all")

  if (current.assay == which.assay){
    new.default.assay <- existing.assays[length(existing.assays)]
    object <- set.default.assay(object, which.assay = new.default.assay)
  }

  if (verbose){
    cat("\n")
    cat(paste("'", which.assay, "' was successfully deleted.", sep = ""))
    cat("\n")
  }

  return(object)
}


#' Get Assay
#'
#' Get name of assays
#'
#' @param object Calibration Object
#' @param which.assays A character specifying which assay(s) to return. Accepted arguments are "default" or "all"
#' \itemize{
#' \item "default" - current assay
#' \item "all" - all existing assays
#' }
#' @name getAssay
#' @return name(s) of assay(s); character
#'
getAssay <- function(object, which.assays = "default"){
  # verbose = TRUE,

  # GIGO handling
  if (!(which.assays %in% c("default", "all"))) stop("'which.assays' must be specified as 'all' or 'default'")

  # get assay names
  if (which.assays == "default"){
    which.assay <- object@current.assay
  } else if (which.assays == "all"){
    which.assay <- names(object@assays)
  }

  return(which.assay)
}

#' Set Default Assay
#'
#' Set default (current) assay. Analyses are performed on default assay only.
#'
#' @param object calibration object
#' @param which.assay specifies which assay is set to default (character)
#' @param verbose print progress
#' @name setDefaultAssay
#' @return calibration object with default assay set to specified assay
#'
setDefaultAssay <- function(object, which.assay, verbose = TRUE){

  stopifnot(exists("which.assay"))
  if (class(which.assay) != "character") {
    stop(paste("'which.assay' is a ", class(which.assay), ". Must be a 'character'", sep = ""))
  }
  # stopifnot(class(which.assay) == "character")
  stopifnot(length(which.assay) == 1)
  stopifnot(which.assay %in% names(object@assays))

  object@current.assay <- which.assay

  if (verbose){
    cat("\n")
    cat(paste("Default assay set to '", which.assay, "'", sep  = ""))
  }

  return(object)
}

#' Get Features
#'
#' Get list of unique features from specified assay in calibration object
#'
#' @param object Calibration Object
#' @param which.assay Character specifying which assay to retrieve features from.
#' @name getFeatures
#' @return list of unique features
#'
getFeatures <- function(object, which.assay = NULL){
  # verbose = TRUE,

  if (is.null(which.assay)) which.assay <- getAssay(object)

  unique.features <- object@assays[[which.assay]]@feature.types

  return(unique.features)
}

#' Get Unique Feature Counts
#'
#' Get unique feature countsfor each ufeature present in specified assay of calibration object
#'
#' @param object Calibration Object
#' @param which.assay Specifies which assay to retrieve feature counts from (character)
#' @name getUniqueFeatureCount
#' @return data frame
#'
getUniqueFeatureCount <- function(object, which.assay = NULL){
  if (is.null(which.assay)) which.assay <- getAssay(object)
  N.features <- object@assays[[which.assay]]@N
  return(as.data.frame(N.features))
}


#' Get dataset names
#'
#' Get names of datasets stored in specified assay in calibration object
#'
#' @param object Calibration Object
#' @param which.assay Specifies which assay to retrieve feature counts from (character)
#' @name getDatasets
#' @return character
#'
getDatasets <- function(object, which.assay = NULL) {
  if (is.null(which.assay)) which.assay <- getAssay(object)
  existing.data <- names(co@assays[[which.assay]]@data)
  return(existing.data)
}


#' Get Analyses
#'
#' Get names of available analyses for specified assay in Calibration Object.
#'
#' @param object Calibration Object
#' @param which.assay Specifies which assay to retrieve feature counts from (character)
#' @name getAnalyses
#' @seealso \code{\link{svpAnalysis}}, \code{\link{mvpAnalysis}}
#' @return Character
#'
getAnalyses <- function(object, which.assay = NULL) {
  if (is.null(which.assay)) which.assay <- getAssay(object)
  existing.analyses <- names(object@assays[[which.assay]]@analysis)
  return(existing.analyses)
}

#' Delete Analysis
#'
#' Delete analysis from specified assay of calibration object
#'
#' @param object calibration object
#' @param which.assay specifies assay
#' @param which.analysis specifies analysis
#' @param verbose print progress
#' @name deleteAnalysis
#' @return character
#'
deleteAnalysis <- function(object, which.assay, which.analysis, verbose = TRUE){

  # GIGO handling
  if (is.null(which.assay)) which.assay <- getAssay(object)
  existing.analyses <- getAnalyses(object, which.assay = which.assay)
  if (!(which.analysis %in% c(existing.analyses, "all"))) stop("user-specified 'which.analysis' does not exist")

  # delete analysis
  if (which.analysis %in% existing.analyses){
    object@assays[[which.assay]]@analysis <- object@assays[[which.assay]]@analysis[!(existing.analyses %in% which.analysis)]

    if (verbose){
      cat("\n")
      cat(paste("'", which.analysis, "' was successfully deleted.", sep = ""))
      cat("\n")
    }

  } else if (which.analysis == "all"){
    object@assays[[which.assay]]@analysis <- list()

    if (verbose){
      cat("\n")
      cat(paste("All analyses were successfully deleted.", sep = ""))
      cat("\n")
    }
  }
  return(object)
}

#' Clone assay object
#'
#' Create duplicate of assay within Calibration Object.
#'
#' @param object Calibration Object
#' @param which.assay Specifies which assay to clone (character)
#' @param cloned.assay.name Name of clone. If not specified, "-copy" is appended to name of cloned Assay.
#' @param set.clone.as.default Logical specifying whether to set clone as default assay
#' @param verbose Logical to report progress.
#' @name cloneAssay
#' @return calibration object
#'
cloneAssay <- function(object, which.assay = NULL, cloned.assay.name = NULL, set.clone.as.default = FALSE, verbose = T) {

  #GIGO handling

  # ensure assay exists
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% getAssay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  }

  # ensure assay is specified
  if (is.null(which.assay)) which.assay <- getAssay(object)

  # ensure new assay name is compatible
  if (is.null(cloned.assay.name)) cloned.assay.name <- paste(getAssay(object), "-copy", sep = "")

  stopifnot(class(cloned.assay.name) == "character")
  if ((cloned.assay.name %in% getAssay(object, which.assay = "all"))){
    stop (paste(cloned.assay.name, "already exists. Specify different 'cloned.assay.name'", sep = ""))
  }

  # clone assay
  object@assays[[cloned.assay.name]] <- object@assays[[which.assay]]

  # set cloned assay to default (optional)
  if (set.clone.as.default) object <- setDefaultAssay(object, which.assay = cloned.assay.name)

  if (verbose){
    cat("\n")
    cat(paste("'", which.assay, "' was successfully cloned", sep = ""))
    cat("\n")
  }

  # return
  return(object)
}

#' Convert multi-tiered list of data.frames to single-tired list of data.frames
#'
#' conversion to single-tiered list of data.frames only supported for 3 tiers currently.
#'
#' @param a list of data.frames
#' @name get.df.from.list
#' @return list
#'
get.df.from.list <- function (a){
  # convert multi-layered list of data.frames (up to 3 tiers) to single-layered list of data.frames

  a.class <- class(a)

  if ("list" %in% a.class){
  out.df.list <- list()
  for (i in 1:length(a)){
    if ("data.frame" %in% class(a[[i]])){
      out.df.list[[names(a)[i]]] <- a[[i]]
      next
    } else if (class(a[i]) == "list"){
      cur.name.i <- names(a)[i]
      for (j in 1:length(a[[i]])){
        if ("data.frame" %in% class(a[[i]][[j]])){
          cur.name.j <- cur.name.i
          out.df.list[[cur.name.j]] <- a[[i]][[j]]
          next
        } else if (class(a[[i]][[j]]) == "list"){
          for (k in 1:length(a[[i]][[j]])){
            if ("data.frame" %in% class(a[[i]][[j]][[k]])){
              cur.name.k <- cur.name.i
              out.df.list[[cur.name.k]] <- a[[i]][[j]][[k]]
              next
            }
          }
        }
      }
    }
  }

  } else if  ("data.frame" %in% a.class){
    out.df.list <- list(a)
  }
  return(out.df.list)

}


#' Get Results Table
#'
#' Returns table as data frame or interactive data table
#'
#' Calibration Object (input option 1) or data frame (input option 2) are accepted as inputs (Calibration Object takes precedent), and tables are reformatted and returned as data.frame or data.table, as specified by 'format' parameter.
#' If Calibration Object is provided, all result tables for specified analysis are retrieved.
#'
#' @param object Calibration Object (input option 1)
#' @param results.table A data frame (input option 2)
#' @param which.assay A character specifying which assay to use. Only specify for input option 1.
#' @param which.results A character specifying function to retrieve results. Only specify for input option 1. One of: (
#' \itemize{
#' \item svp.analysis - returns replicate and rms statistics
#' \item mvp.analysis - results replicate and rms statistics
#' \item fit.calibration - returns calibration equations
#' }
#' @param format Output table format. One of:
#' \itemize{
#' \item df - Data frame. Preferred if user intends to do further data analysis
#' \item dt - Data table. Recommended for visualization and interactive table format. Data.table has interactive option to save table as csv or excel, or copy contents into clipboard.
#' }
#' @param max.page.length Default is 10. Maximal number of entries shown per page. Relevant only for datatable-formated tables (i.e., format = dt)
#' @name getResults
#' @return list of table(s)
#'
getResults <- function(object = NULL, which.results = NULL, which.data = NULL, results.table = NULL, which.assay = NULL, format = "df", max.page.length = 10) {

  if(!(format %in% (c("dt", "df")))) stop("format incorrectly specified")

  # initiate results table list
  tbl.rt.list <- list()

  # if option 1 (Calibration Object provided)
  if (!is.null(object) & !is.null(which.results)){

    # GIGO handling
    # ensure assay is specified
    if (is.null(which.assay)) which.assay <- getAssay(object)

    analysis.flag <- F
    calibration.flag <- F
    # retrive results from analysis output?
    if (any(grepl("vp", which.results))) {
      analysis.flag <- T

      # get existing dataset
      existing.data <- getDatasets(object)

      # handle which data was selected
      if (!is.null(which.data)){
        if (!(which.data %in% existing.data)) stop ("Specified data does not exist")
      } else {
        if (("calibrated" %in% existing.data) & (any(grepl("svp", which.results)))) {
          which.data <- "calibrated"
          warning("returning results for 'calibrated' dataset")
        } else if (("uncalibrated" %in% existing.data) & (any(grepl("svp", which.results)))) {
          which.data <- "uncalibrated"
          warning("returning results for 'uncalibrated' dataset")
        } else if (any(grepl("mvp", which.results))) {

        } else {
          stop("unaccounted for condition encountered. troubleshooting requried")
        }
      }

      # get existing analyses
      existing.analysis <- names(object@assays[[which.assay]]@analysis)

      # check if criteria are fulfilled
      match.ind <- grepl(which.results, existing.analysis)
      if (sum(match.ind) == 0) stop ("queried results do not exist")
      existing.analysis <- existing.analysis[match.ind]

      # check if further filtering is required
      if (length(existing.analysis) > 1) {
        match.ind <- grepl(paste(".", which.data, sep = ""), existing.analysis, fixed = T)
        if (sum(match.ind) == 0) stop ("queried results do not exist")
        existing.analysis <- existing.analysis[match.ind]
      }

      # ensure single result was selected
      if (length(existing.analysis) != 1) stop("Cannot resolve set of queried results)")

      # retrieve results
      mt.df.list <- object@assays[[which.assay]]@analysis[[existing.analysis]]

      if (any(grepl("mvp", which.results))) {
        mt.df.list[["replicate.statistics"]] <- as.data.frame(dplyr::select( mt.df.list[["replicate.statistics"]], -c("value")))
        mt.df.list[["rms.statistics"]][["results"]] <- as.data.frame( dplyr::select( mt.df.list[["rms.statistics"]][["results"]], -c("std.value", "cv.value")))
        mt.df.list[["rms.statistics.no.outliers"]][["results"]] <-  as.data.frame(dplyr::select( mt.df.list[["rms.statistics.no.outliers"]][["results"]], -c("std.value", "cv.value")))
      }


    } else if (any(grepl("calibration", which.results))){
      calibration.flag <- T

      existing.calibration <- names(object@assays[[which.assay]]@calibration)
      match.ind <- grepl("calibration.eq", existing.calibration)

      if (sum(match.ind) > 1) warning(paste("Multiple ", which.results , " exist", sep = ""))
      if (sum(match.ind) == 0) stop(paste(which.results , " results do not exist", sep = ""))

      mt.df.list <- object@assays[[which.assay]]@calibration[[which(match.ind)]]

      if (class(mt.df.list) == "data.frame"){
        mt.df.list <- list(mt.df.list)
        names(mt.df.list) <- existing.calibration[which(match.ind)]
      }
    }

    st.df.list <- get.df.from.list(mt.df.list)

    # name check
    if ((which.results == "fit.calibration") & (length(st.df.list) == 1)){
      names(st.df.list) <- "calibration.equations"
    }

    if (format == "dt"){
      tbl.rt.list <- lapply(st.df.list, datatable, filter="top",
                            width = "100%",
                            height=  "auto",
                            extensions = c('Buttons'),
                            options = list(pageLength = max.page.length,
                                           autoWidth = TRUE,
                                           dom = 'Bfrtip',
                                           buttons = c('copy', 'csv', 'excel')))
      return(tbl.rt.list)
    } else if (format == "df"){
      return(st.df.list)
    }
  } else if (!is.null(results.table)){
    if (format == "dt"){

      tbl.rt.list[[1]] <- datatable(results.table, filter="top",
                                    width = "100%",
                                    height=  "auto",
                                    extensions = 'Buttons',
                                    options = list(pageLength = max.page.length,
                                                   autoWidth = TRUE,
                                                   dom = 'Bfrtip',
                                                   buttons = c('copy', 'csv', 'excel')))
      return(tbl.rt.list)
    } else if (format == "df"){
      return(st.df.list)
    }
  }
}

#' Define replicate sets
#'
#' This function groups data frame entries by specified features, and entries within each grouping as assigned to the same replicate set. Values within a replicate set are expected to be compatible, such that they can be pooled as means and variances prior to computing precision errors.
#'
#' For example, if a phantom was scanned in triplicate, resulting in multiple parameter measurements across several phantom sections, a single replicate set would consist of 3 measurements, all coming from the same section of a given phantom (replicate.strata = c("section", "phantom")). If the same replicate scans were repeated at a later date, or different site, this would be considered a separate replicate set (replicate.strata = c("section", "phantom", "timePoint", "site")).
#'
#' @param df data frame
#' @param replicate.strata Feature grouping used to define a replicate set. Provided as a vector of characters, all which are expected to be column headers in the input data frame.
#' @name defineReplicateSet
#' @return Vector of values (same length as number of entires in input data frame), where the same values are assigned to common replicate sets.
#'
defineReplicateSet <- function(df, replicate.strata){
  # e.g., replicate.strata <- c("scanID", "parameter", "phantom", "site", "timePoint")

  # GIGO HANDLING
  if (!all(replicate.strata %in% colnames(df))) stop("Specified replicate.strata do not exist in data frame")

  # df$replicateSet <- NA
  all.scan.rep <-apply(df[ , replicate.strata], 1, paste, collapse = "-")
  recode.val <- seq(1, length(unique(as.character(all.scan.rep))))
  names(recode.val) <- unique(all.scan.rep)

  replicateSet <- as.numeric(as.character(recode.val[as.character(all.scan.rep)]))

  return(replicateSet)
}

#' Filter features
#'
#' This function filters a dataframe by a specified feature type(s)
#'
#' @param data data frame
#' @param feature Character indicating which feature (i.e., data.frame column) to consider for filtering
#' @param which.feature.type Character vector indicating which feature types to include. If NULL, all are included.
#' @name filterFeatures
#' @return filtered data frame
#'
filterFeatures <- function(data, feature, which.feature.type){
  # filter parameters
  data <- as.data.frame(data)
  u.feat <- unique(as.vector(data[ , feature]))
  if (is.null(which.feature.type)){
    which.feature.type <- u.feat
  } else {
    which.valid <- u.feat %in% which.feature.type
    which.feature.type <- u.feat[which.valid]
  }
  data <- data[data[ , feature] %in% which.feature.type, ]

  return(data)
}

#' Get Calibration information
#'
#' Reports calibration information for specified Assay
#'
#' @param object Calibration Object
#' @param which.assay Character indicating which assay to retrive calibration information for. If unspecified, current Assay is reported.
#' @name getCalibrationInfo
#' @return Calibration information is printed.
#'
getCalibrationInfo <- function(object, which.assay = NULL){

  # ensure assay is specified
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% getAssay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  } else {
    which.assay <- getAssay(object)
  }


  existing.calibrations <- names(object@assays[[which.assay]]@calibration)

  if (length(existing.calibrations) == 0){
    cat(paste("No calibration exists for", which.assay))
  } else {
    cat("Existing Calibration")
    cat(paste("\n\tAssay: ", which.assay, sep = ""))
    cat(paste("\n\tReference Site: ", object@assays[[which.assay]]@calibration[["reference.site"]], sep = ""))
    cat(paste("\n\tReference Time: ", object@assays[[which.assay]]@calibration[["reference.time"]], sep = ""))
    cat(paste("\n\tPhantom: ", object@assays[[which.assay]]@calibration[["phantom"]], sep = ""))
    cat(paste("\n\tParameters: ", paste(object@assays[[which.assay]]@calibration[["parameter"]], collapse = ", "), sep = ""))
  }

}
