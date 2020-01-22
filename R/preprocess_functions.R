

#' Filter dataset
#'
#' Filter HR-pQCT data to only include features of interest
#'
#' @param df HR-pQCT dataset (data.frame) containing value column along with descriptors site, section, timePoint, phantom, parameter and/or scanDate
#' @param analyze.these specifies which features to include in analysis
#' @name filterData
#' @return data.frame
#'
filterData <- function(df, analyze.these){

  # GIGO handles
  stopifnot(class(df) == "data.frame")
  stopifnot(class(analyze.these) == "list")
  stopifnot(all(names(analyze.these) %in% colnames(df)))

  # get features
  u.specifications <- names(analyze.these)

  # copy data
  df.filtered <- df

  # iterate through each feature
  for (i in 1:length(u.specifications)){

    # get current feature
    current.feature <- u.specifications[i]

    # subset data & and force class
    df.subset <- as.character(df.filtered[ ,current.feature])
    analyze.these[[current.feature]] <- as.character(analyze.these[[current.feature]])

    # ensure data are same class
    if (class(df.subset) != class(analyze.these[[current.feature]])) class(df.subset) <- class(analyze.these[[current.feature]])
    stopifnot(class(df.subset) == class(analyze.these[[current.feature]]))

    # get indices of entries to include
    include.ind <- which(df.subset %in% analyze.these[[current.feature]])

    # filter data
    df.filtered <- df.filtered[include.ind, ]

  }

  return(df.filtered)

}

#' Preprocess data
#'
#' Preprocess HR-pQCT data for subsequent analysis
#'
#' @param object calibration object
#' @param analyze.which list specifying which features to include in analysis (generate using analyzeWhich() function)
#' @param new.assay.name name of new preprocessed assay
#' @param which.assay specifies which assay to preprocess. Set to current.assay if which.assay is unspecified.
#' @name preprocess.data
#' @return calibration object
#'
preprocess.data <- function(object, analyze.which, new.assay.name, which.assay = NULL) {


  #GIGO handling
  if (is.null(which.assay)) which.assay <- get.assay(object, verbose = FALSE)
  stopifnot(exists("new.assay.name"))
  stopifnot(class(analyze.which) == "list")
  stopifnot(class(new.assay.name) == "character")
  if (new.assay.name %in% get.assay(object, verbose = FALSE, which.assays = "all")) {
    stop(paste("'", new.assay.name, "' already exists.", sep = ""))
  }

  # get data
  df <- object@assays[[which.assay]]@data[["uncalibrated.data"]]

  # filter data
  df <- filterData(df, analyze.which)

  # get unique features
  uf.output <- get.unique.features(df)

  # create new assay object
  as <- new("assay",
            data = list(uncalibrated = df[ ,uf.output$variables]),
            variables = as.character(uf.output$variables),
            unique.features = uf.output$unique.features,
            N = uf.output$N,
            description = paste("Source: ", which.assay, sep = ""))
  # as <- list(as)

  # assign new assay to existing list of assays
  existing.as <- object@assays
  existing.assay.names <- names(existing.as)



  # REVISIT THIS, ENSURE PROPER NAME HANDLING
  # if (is.null(new.assay.name) & !("input" %in% existing.assay.names)) new.assay.name <- "input"

  existing.as[length(existing.as)+1] <- as
  names(existing.as)[length(existing.as)] <- new.assay.name

  object@assays <- existing.as

  object <- set.default.assay(object, new.assay.name)

  return(object)
}



#' Specify features to analyze
#'
#' Preprocess HR-pQCT data for subsequent analysis.
#' Inclusion parameters are considered first, followed by omission parameters.
#' For inclusion parameters, "all", c("parameter_names") or NULL are supported inputs.
#' For omission parameters, c("parameter_names") or NULL are supported inputs.
#'
#' @param object calibration object
#' @param include.sites sites to include in analysis
#' @param include.sections sections to include in analysis
#' @param include.times timePoints to include in analysis
#' @param include.phantoms phantoms to include in analysis
#' @param include.parameters parameters to include in analysis
#' @param include.scanDates scanDates to include in analysis
#' @param omit.sites sites to omit from analysis
#' @param omit.sections phantom sections to omit from analysis
#' @param omit.times timePoints to omit from analysis
#' @param omit.phantoms phantoms to omit from analysis
#' @param omit.parameters parameters to omit from analysis
#' @param omit.scanDates scanDates to omit from analysis
#' @param which.assay specifies assay
#' @name analyzeWhich
#' @return list
#'
analyzeWhich<-function(object, include.sites = "all", include.sections = "all", include.times = "all",
                       include.phantoms = "all", include.parameters = NULL, include.scanDates = "all",
                       omit.sites = NULL, omit.sections = NULL, omit.times = NULL, omit.phantoms = NULL,
                       omit.parameters = NULL, omit.scanDates = NULL, which.assay = NULL){

  # inclusion parameters are processed first, followed by omission parameters.
  # options: inclusion parameters: "all", c("parameter names"), NULL
  #          omission parameters: c("parameter names"), NULL

  # Inclusion list
  inclusion.list <- list(site = include.sites,
                         section = include.sections,
                         timePoint = include.times,
                         phantom = include.phantoms,
                         parameter = include.parameters,
                         scanDate = include.scanDates)

  # Omission list
  omission.list <- list(site = omit.sites,
                        section = omit.sections,
                        timePoint = omit.times,
                        phantom = omit.phantoms,
                        parameter = omit.parameters,
                        scanDate = omit.scanDates)


  # Available features list
  if (is.null(which.assay)) which.assay <- get.assay(object, verbose = FALSE)
  available.features <- get.features(object, verbose = FALSE, which.assay = which.assay)
  available.feature.names <- names(available.features)

  # ensure features are non-factor class
  for (i in 1:length(available.features)){
    if (class( available.features[[i]] ) == "factor") {
      available.features[[i]] <- as.character(available.features[[i]])
    }
  }

  # check that inclusion and omission list is available in queried calibration object
  inclusion.availability.check <- NULL
  omission.availability.check <- NULL
  for (i in 1:length(inclusion.list)){
    if ((names(inclusion.list)[i] %in% available.feature.names) & (!is.null(inclusion.list[[i]]))){
      inclusion.availability.check[i] <- TRUE
    } else  inclusion.availability.check[i] <- FALSE
    if ((names(omission.list)[i] %in% available.feature.names) & (!is.null(omission.list[[i]]))){
      omission.availability.check[i] <- TRUE
    } else  omission.availability.check[i] <- FALSE

  }

  # filter out nonsense variables
  inclusion.list <- inclusion.list[inclusion.availability.check]
  omission.list <- omission.list[omission.availability.check]

  # retain available features (pass only those features that are available, drop all else)
  analyze.these <- list()
  if (all(as.vector(unlist(inclusion.list)) == "all")) {
    for (i in 1:length(inclusion.list)){
      if (length(available.features[[names(inclusion.list)[i]]]) > 0){

        analyze.these[[names(inclusion.list)[i]]] <- (available.features[[names(inclusion.list)[i]]])

      }
    }
  } else {
    for (i in 1:length(inclusion.list)){
      if (inclusion.list[[i]] == "all"){
        analyze.these[[names(inclusion.list)[i]]] <- available.features[[names(inclusion.list)[i]]]
      } else {
        match.ind <- inclusion.list[[i]] %in%  available.features[[names(inclusion.list)[i]]]
        if (sum(match.ind) == 0) stop("Must include atleast one feature in analysis")
        analyze.these[[names(inclusion.list)[i]]] <- inclusion.list[[i]][match.ind]
      }
    }
  }

  # omit any remaining features that were specified in omission list
  if (length(omission.list) > 0){
    for (i in 1:length(omission.list)){
      match.ind <- analyze.these[[names(omission.list)[i]]] %in%  omission.list[[i]]
      if (!(sum(match.ind) < length(match.ind))) stop("Cannot omit all data")
      if (sum(match.ind) > 0){
        analyze.these[[names(omission.list)[i]]] <- analyze.these[[names(omission.list)[i]]][!match.ind]
      }

    }

  }
  return(analyze.these)
}


#' Omit outliers
#'
#' Credit: aL3xa, StackOVerflow (https://stackoverflow.com/a/4788102)
#'
#' @param x values
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds.
#' @name omit.outliers
#' @return values
#'
omit.outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}




