#' Get Unique Features
#'
#' Get unique features in input dataset
#'
#'
#' @param df HR-pQCT dataset (data.frame). Following column names are expected: value, site, section, timePoint, phantom, parameter, scanDate
#' @name get.unique.features
#' @return list of unique.features (unique.features), number of variables per unique feature (N), and names of variables (variable)
#'
get.unique.features <- function(df){

  # prep data for assay class
  accepted.variables = c("value", "site", "section", "timePoint", "phantom", "parameter", "scanDate")
  available.variables <- colnames(df)

  match.ind <- accepted.variables  %in%  available.variables
  variables <- accepted.variables[match.ind]

  unique.features <- list()
  N <- list()
  N[["n.entries"]] <- nrow(df)
  for (i in 1:length(accepted.variables)){
    if (accepted.variables[i] != "value"){
      n.var.name <- paste("n.", accepted.variables[i], sep = "")
      if (accepted.variables[i] %in% variables){
        unique.features[[accepted.variables[i]]] <- unique(df[, accepted.variables[i]])
        N[[n.var.name]] <- length(unique.features[[accepted.variables[i]]])
      } else if ((accepted.variables[i] == "parameter") |
                 (accepted.variables[i] == "timePoint") |
                 (accepted.variables[i] == "sections")) {
        unique.features[[accepted.variables[i]]] <- NULL
        N[[n.var.name]] <- 1
      }

    }
  }

  uf.output <- list(
    unique.features = unique.features,
    N = N,
    variables= variables
  )
  return(uf.output)
}


#' Delete Assay
#'
#' Delete assay from calibration object
#'
#' @param object calibration object
#' @param which.assay which assay to delete (character)
#' @param verbose print progress (Deafult = TRUE)
#' @name deleteAssay
#' @return calibration object with specified assay deleted. If specified assay was current assay, new current assay is assigned.
#'
deleteAssay <- function(object, which.assay, verbose = TRUE){

  existing.assays <- get.assay(object, verbose = FALSE, which.assays = "all")

  if (!(which.assay %in% existing.assays)){
    stop(paste("'", which.assay, "' does not exist. No assays were deleted.", sep = ""))
  }

  # delete assay
  object@assays <- object@assays[!(existing.assays %in% which.assay)]

  # if deleted assay was default, set new default assay
  current.assay <- get.assay(object, verbose = FALSE)
  existing.assays <- get.assay(object, verbose = FALSE, which.assays = "all")

  if (current.assay == which.assay){
    new.default.assay <- existing.assays[length(existing.assays)]
    object <- set.default.assay(object, which.assay = new.default.assay)
  }

  if (verbose){
    cat("\n============================\n")
    cat(paste("'", which.assay, "' was successfully deleted.", sep = ""))
    cat("\n")
  }

  return(object)
}


#' Get Assay
#'
#' get name of assays (current or all)
#'
#' @param object calibration object
#' @param verbose print progress (Default = TRUE)
#' @param which.assays which assay(s) to return. Accepted arguments are "default" or "all"
#' @name get.assay
#' @return name(s) of assay(s); character
#'
get.assay <- function(object, verbose = TRUE, which.assays = "default"){

  # GIGO handling
  if (!(which.assays %in% c("default", "all"))) stop("'which.assays' must be specified as 'all' or 'default'")

  # get assay names
  if (which.assays == "default"){
    which.assay <- object@current.assay
  } else if (which.assays == "all"){
    which.assay <- names(object@assays)
  }

  if (verbose){
    cat("\n============================\n")
    if (which.assays == "default"){
      cat("Default Assay:\n")
      cat(paste(which.assay))
    } else if (which.assays == "all"){
      cat("Existing Assay(s):\n")
      cat(paste(which.assay, collapse = ", "))
    }
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
#' @name set.default.assay
#' @return calibration object with default assay set to specified assay
#'
set.default.assay <- function(object, which.assay, verbose = TRUE){

  stopifnot(exists("which.assay"))
  if (class(which.assay) != "character") {
    stop(paste("'which.assay' is a ", class(which.assay), ". Must be a 'character'", sep = ""))
  }
  # stopifnot(class(which.assay) == "character")
  stopifnot(length(which.assay) == 1)
  stopifnot(which.assay %in% names(object@assays))

  object@current.assay <- which.assay

  if (verbose){
    cat("\n============================\n")
    cat(paste("Default assay set to '", which.assay, "'", sep  = ""))
  }
  return(object)
}

#' Get Features
#'
#' Get list of unique features from specified assay in calibration object
#'
#' @param object calibration object
#' @param verbose print progress
#' @param which.assay specifies which assay to retrieve features from (character)
#' @name get.features
#' @return list of unique features
#'
get.features <- function(object, verbose = TRUE, which.assay = NULL){

  if (is.null(which.assay)) which.assay <- get.assay(object, verbose = FALSE)

  unique.features <- object@assays[[which.assay]]@unique.features

  if (verbose){
    cat("\n============================\n")
    cat(paste("Assay: ", which.assay, sep = ""))
    cat("\n")

    available.variables <- names(unique.features)
    for (i in 1:length(available.variables)){
      if (available.variables[i] != "value") {

        cat(paste(available.variables[i], "(s)", sep = ""))
        cat(": \n")
        cat( paste(unique.features[[i]], collapse = ", "))
        cat("\n")
      }

    }
  }
  return(unique.features)
}

#' Get Unique Feature Counts
#'
#' Get counts of variables for each unique feature present in specified assay of calibration object
#'
#' @param object calibration object
#' @param verbose print progress
#' @param which.assay specifies which assay to retrieve feature counts from (character)
#' @name get.unique.feature.count
#' @return list
#'
get.unique.feature.count <- function(object, verbose = TRUE, which.assay = NULL){

  if (is.null(which.assay)) which.assay <- get.assay(object, verbose = FALSE)

  N.features <- object@assays[[which.assay]]@N

  if (verbose){
    cat("\n============================\n")
    cat(paste("Assay: ", which.assay, sep = ""))
    cat("\n")

    list.names <- names(N.features)
    for (i in 1:length(list.names)){
      cat(paste(list.names[i], ": ", sep = ""))
      cat(paste(N.features[[i]], collapse = ", "))
      cat("\n")
    }
  }

  return(N.features)
}

#' Get Analyses
#'
#' Get names of avaialble analyses for specified assay in calibration object
#'
#' @param object calibration object
#' @param verbose print progress
#' @param which.assay specifies which assay to retrieve feature counts from (character)
#' @name get.analyses
#' @return character
#'
get.analyses <- function(object, verbose = TRUE, which.assay = NULL) {
  if (is.null(which.assay)) which.assay <- get.assay(object, verbose = FALSE)

  existing.analyses <- names(object@assays[[which.assay]]@analysis)

  if (verbose){
    cat("\n============================\n")
    if (length(existing.analyses) == 0){
      cat(paste("No analyses exist for '", which.assay, "'\n", sep = ""))
    } else {
      cat(paste("Existing analyses for '", which.assay, "': ", sep = ""))
      cat(paste(existing.analyses, collapse = ", "))
      cat("\n")
    }
  }
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
  if (is.null(which.assay)) which.assay <- get.assay(object, verbose = FALSE)
  existing.analyses <- get.analyses(object, verbose = FALSE, which.assay = which.assay)
  if (!(which.analysis %in% c(existing.analyses, "all"))) stop("user-specified 'which.analysis' does not exist")

  # delete analysis
  if (which.analysis %in% existing.analyses){
    object@assays[[which.assay]]@analysis <- object@assays[[which.assay]]@analysis[!(existing.analyses %in% which.analysis)]

    if (verbose){
      cat("\n============================\n")
      cat(paste("'", which.analysis, "' was successfully deleted.", sep = ""))
      cat("\n")
    }

  } else if (which.analysis == "all"){
    object@assays[[which.assay]]@analysis <- list()

    if (verbose){
      cat("\n============================\n")
      cat(paste("All analyses were successfully deleted.", sep = ""))
      cat("\n")
    }
  }
  return(object)
}

#' Clone assay object
#'
#' Create duplicate of assay within calibration object
#'
#' @param object calibration object
#' @param which.assay specifies which assay to clone (character)
#' @param cloned.assay.name name of clone. If not specified, "-copy" is appended to name of cloned assay.
#' @param set.clone.as.default logical specifying whether to set clone as default assay
#' @name clone.assay
#' @return calibration object
#'
clone.assay <- function(object, which.assay = NULL, cloned.assay.name = NULL, set.clone.as.default = FALSE) {


  #GIGO handling

  # ensure assay exists
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% get.assay(object, verbose = FALSE, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  }

  # ensure assay is specified
  if (is.null(which.assay)) which.assay <- get.assay(object, verbose = FALSE)

  # ensure new assay name is compatible
  if (is.null(cloned.assay.name)) cloned.assay.name <- paste(get.assay(object, verbose = FALSE), "-copy", sep = "")

  stopifnot(class(cloned.assay.name) == "character")
  if ((cloned.assay.name %in% get.assay(object, verbose = FALSE, which.assay = "all"))){
    stop (paste(cloned.assay.name, "already exists. Specify different 'cloned.assay.name'", sep = ""))
  }

  # clone assay
  object@assays[[cloned.assay.name]] <- object@assays[[which.assay]]

  # set cloned assay to default (optional)
  if (set.clone.as.default) object <- set.default.assay(object, which.assay = cloned.assay.name)

  # return
  return(object)
}


