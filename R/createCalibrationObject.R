
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

  df$value <- as.numeric(as.character(df$value))
  try(df$site <- as.character(df$site), silent = TRUE)
  try(df$section <- as.character(df$section), silent = TRUE)
  try(df$timePoint <- as.numeric(as.character(df$timePoint)), silent = TRUE)
  if (!is.numeric(df$timePoint))  try(df$timePoint <- as.character(df$value), silent = TRUE)
  try(df$phantom <- as.character(df$phantom), silent = TRUE)
  try(df$parameter <- as.character(df$parameter), silent = TRUE)

  uf.output <- get.unique.features(df)

  # define assay class
  as <- new("assay",
            data = list(uncalibrated.data = df[ ,uf.output$variables]),
            variables = as.character(uf.output$variables),
            unique.features = uf.output$unique.features,
            N = uf.output$N,
            description = assay.name)
  as <- list(as)
  names(as)[1] <- assay.name


  # prep meta data
  meta.data <- data.frame(site = as.character(df$site),
                          time = df$timePoint,
                          scanID = as.character(df$scanNumber),
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


