#' Calibration Class
#'
#' This class defines a calibration object
#'
#' @slot assays list of assay objects
#' @slot metadata input dataset (data.frame)
#' @slot current.assay species current assay on which analyses are called
#' @slot project.name project name
#' @slot misc miscilaneous features
#' @name calibration-class
#'
calibration <- setClass("calibration",
                        slots = c(assays = "list",
                                  metadata = "data.frame",
                                  current.assay = "character",
                                  project.name = "character",
                                  misc = "list")
)
#' Assay Class
#'
#' This class defines an assay object
#'
#' @slot data list of datasets
#' @slot analysis list of analyses
#' @slot calibration calibration status (pre- or post-calibration)
#' @slot plots list of plots generated
#' @slot features features available for analysis in dataset
#' @slot feature.types list of features provided in dataset
#' @slot N list of number of available variables for each feature
#' @slot description description
#' @name assay-class
#'
assay <- setClass("assay", slots = c(data = "list",
                                     analysis = "list",
                                     calibration  = "list",
                                     plots = "list",
                                     features = "character",
                                     feature.types = "list",
                                     N = "list",
                                     description = "character"),
                  prototype = list(
                    features = NA_character_,
                    description = NA_character_
                  ))

# devtools::document()
