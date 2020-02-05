
#' Create Calibration Object
#'
#' The Calibration object is the engine that drives each calibration and precision analysis.
#' It stores all information associated with the dataset, including the original input,
#' transformations, analyses, results and plots. All that is needed to construct a Calibration
#' object is raw data formatted as a data.frame. See documentation for details
#'
#'
#' @param df HR-pQCT dataset (data.frame) containing value column along with descriptors site, section, timePoint, phantom, parameter and/or scanDate
#' @param assay.name assay name
#' @param project.name project name
#' @name createCalibrationObject
#' @return Calibration Object
#'
createCalibrationObject <- function(df, assay.name = NULL, project.name = NULL){

  if (is.null(assay.name)) assay.name <- "input"

  if (!("value" %in% colnames(df))) stop("values were not specified")
  df$value <- as.numeric(as.character(df$value))
  try(df$site <- as.character(df$site), silent = TRUE)
  try(df$scanID <- as.character(df$scanID), silent = TRUE)
  try(df$section <- as.character(df$section), silent = TRUE)
  try(df$timePoint <- as.numeric(as.character(df$timePoint)), silent = TRUE)
  if (!is.numeric(df$timePoint))  try(df$timePoint <- as.character(df$timePoint), silent = TRUE)
  try(df$phantom <- as.character(df$phantom), silent = TRUE)
  try(df$parameter <- as.character(df$parameter), silent = TRUE)

  uf.output <- getUniqueFeatures(df)

  # handle missing features
  expected.features <- c("value", "site", "scanID", "section", "timePoint", "phantom", "parameter", "scanDate")
  available.features <- uf.output[["variables"]]
  missing.features <- expected.features[!(expected.features %in% available.features)]

  if ("timePoint" %in% missing.features) {
    cat("\n")
    warning("timePoint feature was not specified in input dataset; timePoint feature created and set to 'baseline'.
            If multiple timePoints expected, ensure that input is correctly specified")
    cat("\n")
    df$timePoint <- 'baseline'
  }
  if ("phantom" %in% missing.features) {
    cat("\n")
    warning("phantom feature was not specified in input dataset; phantom feature created and set to 'unnamed.phantom'")
    cat("\n")
    df$phantom <- "unnamed.phantom"
  }
  if ("site" %in% missing.features) {
    cat("\n")
    warning("site feature was not specified in input dataset; site feature created and set to 'unnamed.site'")
    cat("\n")
    df$site <- "unnamed.site"
  }
  if ("parameter" %in% missing.features) {
    cat("\n")
    warning("parameter feature was not specified in input dataset; parameter feature created and set to 'unnamed.parameter'")
    cat("\n")
    df$parameter <- "unnamed.parameter"
  }
  if ("section" %in% missing.features) {
    cat("\n")
    warning("section feature was not specified in input dataset; parameter feature created and set to 'unnamed.section'.
            Note that it is impossible to cross-calibrate instruments with a single section type.")
    cat("\n")
    df$section <- "unnamed.section"
  }
  if ("scanID" %in% missing.features) {
    cat("\n")
    warning("scanID feature was not specified in input dataset; parameter feature created and set to 'unnamed.scan.id'")
    cat("\n")
    df$scanID <- "unnamed.scan.id"
  }
  if ("scanDate" %in% missing.features) {
    cat("\n")
    warning("scanDate feature was not specified in input dataset; scanDate feature created and set to NA")
    cat("\n")
    df$scanDate <- NA
  }

  uf.output <- getUniqueFeatures(df)

  # define assay class
  as <- new("assay",
            data = list(uncalibrated.data = df[ ,uf.output$variables]),
            features = as.character(uf.output$variables),
            feature.types = uf.output$unique.features,
            N = uf.output$N,
            description = assay.name)

  # store assay in list
  as <- list(as)
  names(as)[1] <- assay.name


  # prep meta data
  meta.data <- data.frame(site = as.character(df$site),
                          time = df$timePoint,
                          scanID = as.character(df$scanID),
                          section = as.character(df$section),
                          parameter = as.character(df$parameter),
                          phantom = as.character(df$phantom),
                          scanDate = df$scanDate)


  # define calibration class
  if (is.null(project.name)) project.name <- "calibration_project"
  co <- new("calibration",
            assays = as,
            metadata = meta.data,
            current.assay = assay.name,
            project.name = project.name)
  return(co)

}


