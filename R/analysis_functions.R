#' Single variant precision/reproducibility analysis
#'
#' Short-term single-center precision analysis of repeated phantom/in vivo measurements.
#' Parameter-specific precision errors are calculated as root-mean-square standard deviation (rms.std) and root-mean-square coefficient of variance (rms.cv).
#' Results are saved in calibration object as 'svp.analysis'.
#'
#' @param object Calibration Object
#' @param which.data A character specifying which data to analyze. One of:
#' \itemize{
#' \item uncalibrated (Default)
#' \item calibrated
#' }
#' @param replicate.strata Vector specifying how a set of scan replicates are defined. Default replicate definition is stratification by site, timePoint, section, parameter, phantom. That is, all scans grouped by the specified features, and pooled within each strata prior to pooling across strata to compute parameter-specific rms-statistics.
#' @param which.assay A character specifying assay for analysis. If unspecified, set to current assay.
#' @param n.signif Number of significant digits reported.
#' @param verbose Logical specifying where to report progress.
#' @name svp.analysis
#' @return Calibration Object
#'
svp.analysis <- function(object, which.data = "uncalibrated",  replicate.strata = NULL, n.signif = 3, which.assay = NULL, verbose = T){

  # calibrated data: logical (true/false); if null, checks for calibrated data first.

  #GIGO handling
  if (!is.null(which.assay)) stopifnot(class(which.assay) == "character")
  if (is.null(which.assay)) which.assay <- get.assay(object)
  stopifnot((class(n.signif) == "numeric") & (length(n.signif) == 1))
  if (!(which.data %in% names(object@assays[[which.assay]]@data))) stop(paste("'", which.data , "' data does not exist", sep = ""))

# get grouping features
  ufeatures <- get.features(object = object, which.assay = which.assay)
  ufeatures.names <- names(ufeatures)

  if (is.null(replicate.strata)){
      group.by <- c("site", "timePoint", "section", "parameter", "phantom")
  } else if (sum(replicate.strata %in% ufeatures.names)> 0){
    group.by <- ufeatures.names[ufeatures.names %in% replicate.strata]
  } else {
    stop("'replicate.strata' is incorrectly specified")
  }

  # get data
  df <- object@assays[[which.assay]]@data[[which.data]]

  # ensure values are numeric
  df$value <- as.numeric(as.vector(df$value))

  # calculate precision errors
  df.stats <- df %>%
    dplyr::group_by(.dots = group.by) %>%
    dplyr::summarize(mean.value = mean(value, na.rm = T),
              median.value = median(value, na.rm = T),
              std.value = sd(value, na.rm = T),
              cv.value = (sd(value, na.rm = T)/ mean(value, na.rm = T)),
              n.scans = length(value))

  # flag outliers (CHECK)
  df.stats <- df.stats %>%
    dplyr::group_by(.dots = c("parameter")) %>%
    dplyr::mutate(outlier.flag = omit.outliers(cv.value))
  df.stats$outlier.flag <- is.na(df.stats$outlier.flag)

    # outliers.as.na <- as.numeric(omit.outliers(as.matrix(df.stats[ , "cv.value"])))
    # df.stats$outlier.flag <- F
    # df.stats$outlier.flag <- is.na(outliers.as.na)
    n.outliers <- sum(df.stats$outlier.flag)
    if (verbose){
      if (n.outliers > 0) warning(paste("\nWarning: ", n.outliers, " outliers detected \n", sep = ""))
    }

  # parameter-level precision errors
    # check if parameter is specified - if not, assume single parameter
    if (("parameter" %in% ufeatures.names) & (length(unique(df.stats$parameter)) > 1)){
      df.stats_pooled <- df.stats %>%
        dplyr::group_by(parameter) %>%
        dplyr::summarize(rms.std = signif(sqrt(mean(std.value^2)), n.signif),
                         median.std = signif(median(std.value), n.signif),
                         iqr.std = paste(signif(quantile(std.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                         rms.cv = signif(sqrt(mean(cv.value^2)), n.signif),
                         median.cv = signif(median(cv.value), n.signif),
                         iqr.cv = paste(signif(quantile(cv.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                         n.scans = sum(n.scans),
                         n.independent = round(length(cv.value)/ length(unique(section))))

      df.stats_pooled.no.outliers <- df.stats %>%
        dplyr::filter(outlier.flag == FALSE) %>%
        dplyr::group_by(parameter) %>%
        dplyr::summarize(rms.std = signif(sqrt(mean(std.value^2)), n.signif),
                         median.std = signif(median(std.value), n.signif),
                         iqr.std = paste(signif(quantile(std.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                         rms.cv = signif(sqrt(mean(cv.value^2)), n.signif),
                         median.cv = signif(median(cv.value), n.signif),
                         iqr.cv = paste(signif(quantile(cv.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                         n.scans = sum(n.scans),
                         n.independent = round(length(cv.value)/ length(unique(section))))
    } else {
      df.stats_pooled <- df.stats %>%
        dplyr::summarize(rms.std = signif(sqrt(mean(std.value^2)), n.signif),
                         median.std = signif(median(std.value), n.signif),
                         iqr.std = paste(signif(quantile(std.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                         rms.cv = signif(sqrt(mean(cv.value^2)), n.signif),
                         median.cv = signif(median(cv.value), n.signif),
                         iqr.cv = paste(signif(quantile(cv.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                         n.scans = sum(n.scans),
                         n.independent = round(length(cv.value)/ length(unique(section))))

      df.stats_pooled.no.outliers <- df.stats %>%
        dplyr::filter(outlier.flag == FALSE) %>%
        dplyr::summarize(rms.std = signif(sqrt(mean(std.value^2)), n.signif),
                         median.std = signif(median(std.value), n.signif),
                         iqr.std = paste(signif(quantile(std.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                         rms.cv = signif(sqrt(mean(cv.value^2)), n.signif),
                         median.cv = signif(median(cv.value), n.signif),
                         iqr.cv = paste(signif(quantile(cv.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                         n.scans = sum(n.scans),
                         n.independent = round(length(cv.value)/ length(unique(section))))
    }


    # get extra features and append to df.replicate
    input.features <- colnames(df)
    output.features <- colnames(df.stats)

    # find extra features
    which.extra <- input.features[!(input.features %in% output.features)]

    df.extra.list <- list()
    df.extra <- df.stats[ ,group.by]
    for (i in 1:length(which.extra)){

      if (is.factor(df[, which.extra[i]])) df[, which.extra[i]] <- as.character(df[, which.extra[i]])
      if (class(df[, which.extra[i]]) == "Date") df[, which.extra[i]] <- as.character(df[, which.extra[i]])

      df.extra.list[[which.extra[i]]] <- df %>%
        dplyr::group_by(.dots = c(group.by))  %>%
        dplyr::summarize(mean.value = list(unique(get(which.extra[i]))))

      colnames(df.extra.list[[which.extra[i]]])[dim(df.extra.list[[which.extra[i]]])[2]] <- which.extra[i]

      if (all(unlist( lapply( df.extra.list[[which.extra[i]]][[which.extra[i]]], length)) == 1)){
        df.extra.list[[which.extra[i]]][[which.extra[i]]] <-  unlist(df.extra.list[[which.extra[i]]][[which.extra[i]]])
      }
      df.extra <- merge(df.extra,  df.extra.list[[which.extra[i]]], by =  group.by)
    }

    # enforce significant figures
    df.stats_sigfigs <- df.stats
    is.num <- sapply(df.stats_sigfigs, is.numeric)
    df.stats_sigfigs[is.num] <- lapply(df.stats_sigfigs[is.num], signif, n.signif)

    df.stats_sigfigs <- merge(df.stats_sigfigs,  df.extra, by =  group.by)

  # store results
    replicate.statistics = list(
      results = df.stats_sigfigs,
      grouping.variable = group.by)

  results.list <- list(
    replicate.statistics = replicate.statistics,
    rms.statistics = list(
      results = df.stats_pooled,
      grouping.variable = "parameter",
      notes = "all data included"),
    rms.statistics.no.outliers = list(
      results = df.stats_pooled.no.outliers,
      grouping.variable = "parameter",
      notes = "outliers omitted"
    )
  )

  # assign names
  analysis.name <- paste("svp.analysis", ".", which.data, sep = "")
  existing.analyses <- object@assays[[which.assay]]@analysis
  existing.analysis.names <- names(existing.analyses)

  if (verbose){
    if (analysis.name %in% existing.analysis.names){
      warning(paste("Pre-existing '", analysis.name, "' Analysis in '", which.assay , "' Assay was overwritten", sep = ""))
    } else {
      cat("\n")
      cat(paste("analysis saved as '", analysis.name, "'", sep = ""))
      cat("\n")}
  }


  existing.analyses[[analysis.name]] <- results.list
  object@assays[[which.assay]]@analysis <- existing.analyses

# return calibration object
  return(object)
}





#' Multi variant precision/reproducibility analysis
#'
#' Longitudinal (i.e., multiple time points) and multicenter (i.e., multiple instruments) precision analysis of
#' repeated phantom measurements. Precisions are calculated as root-mean-square standard deviation (rms.std) and
#' root-mean-square coefficient of variance (rms.cv). Results are saved in Calibration Object as 'mvp.analysis'.
#'
#' @param object Calibration Object
#' @param which.data A character specifying which data to analyze. One of:
#' \itemize{
#' \item all (Default)
#' \item uncalibrated
#' \item calibrated
#' }
#' @param replicate.strata Vector specifying how a set of scan replicates are defined. Default replicate definition is stratification by site, timePoint, section, parameter, phantom. That is, all scans grouped by the specified features, and pooled within each strata prior to pooling across strata to compute parameter-specific rms-statistics.
#' @param which.assay A character specifying assay for analysis. If unspecified, set to current assay.
#' @param n.signif number of significant digits to report
#' @param verbose Logical specifying where to report progress.
#' @name mvp.analysis
#' @return Calibration Object
#'
mvp.analysis <- function(object, which.data = "all", replicate.strata = NULL, which.assay = NULL, n.signif = 3, verbose = T){

  # @param var2plot precision error to plot (cv.value or std.value)
  # @param show.data.table logical specifying whether to generate datatable summary of results
  # @param trim.plot logical specifying whether to trim outliers from plots
  # @param plot.flag logical specifying whether to plot results

  # calibrated data: logical (true/false); if null, checks for calibrated data first.

  #GIGO handling
  if (!is.null(which.assay)) stopifnot(class(which.assay) == "character")
  if (is.null(which.assay)) which.assay <- get.assay(object)
  stopifnot((class(n.signif) == "numeric") & (length(n.signif) == 1))

  # get grouping features
  ufeatures <- get.features(object = object, which.assay = which.assay)
  ufeatures.names <- names(ufeatures)

  if (is.null(replicate.strata)){
    replicate.strata <- c("site", "timePoint", "section", "parameter", "phantom")
  } else if (sum(replicate.strata %in% ufeatures.names)> 0){
    replicate.strata <- ufeatures.names[ufeatures.names %in% replicate.strata]
  } else {
    stop("'replicate.strata' is incorrectly specified")
  }

  existing.datasets <- object@assays[[which.assay]]@data
  existing.data.names <- names(existing.datasets)
  if (which.data == "all"){
    df <- NULL
    for (i in 1:length(existing.data.names)){
      df.cur <- existing.datasets[[existing.data.names[i]]]
      df.cur$which.data <- existing.data.names[i]
      df <- bind_rows(df, df.cur)
    }
  } else {
    df <- object@assays[[which.assay]]@data[[which.data]]
    df$which.data <- which.data
  }

  # set datatype (OMIT?)
  # if (var2plot == "cv"){
  #   val <- "CV"
  #   var2plot <- "cv.value"
  #   pooled.par <- "rms.cv"
  # } else if (var2plot == "std"){
  #   val <- "STD"
  #   var2plot <- "std.value"
  #   pooled.par <- "rms.std"
  # }

  # ensure values are numeric
  df$value <- as.numeric(as.vector(df$value))

  # reorder datasets
  u.data <- unique(df$which.data)
  if (("uncalibrated" %in% u.data) & ("uncalibrated" %in% u.data)){
    df$which.data <- factor(df$which.data)
    match.uncal <- match("uncalibrated", levels(df$which.data))
    match.cal <- match("calibrated", levels(df$which.data))
    df$which.data <- factor(df$which.data, levels(df$which.data)[c(match.uncal, match.cal)])
  } else {df$which.data <- factor(df$which.data)}

  # unique timepoints and sites
  # u.site <- unique(df$site)
  # u.time <- unique(df$timePoint)
  group.by <- c(replicate.strata, "which.data")

  # single center single time
  df.unpooled <- NULL

  # define baseline
  u.time <- unique(df$timePoint)
  if ("baseline" %in% u.time){
    t.base <- "baseline"
  } else {
    t.base <- min(u.time)
    if (class(t.base) != "numeric") stop ("Could not identify baseline timePoint")
  }
  t.followup <- u.time[!(u.time %in% t.base)]

  df.ssst <- NULL
  try({
    group.by.ssst <- group.by
    df.ssst <- df %>%
      dplyr::group_by(.dots = group.by.ssst) %>%
      dplyr::summarize(mean.value = mean(value, na.rm = T),
                median.value = median(value, na.rm = T),
                std.value = sd(value, na.rm = T),
                cv.value = (sd(value, na.rm = T)/ mean(value, na.rm = T)),
                n.scans = length(value))

    outliers.as.na <- as.numeric(omit.outliers(as.matrix(df.ssst[ , "cv.value"])))
    df.ssst$outlier.flag <- F
    df.ssst$outlier.flag <- is.na(outliers.as.na)

    df.ssst$t.span <- gsub(" ", "",paste("t", as.character(t.base), sep = ""))
    df.ssst$s.span <- "single"
  }, silent = TRUE)

  # single center long time
  df.sslt <- NULL
  try({
    group.by.sslt <- group.by[!(group.by %in% "timePoint")]
    for (i in 1:length(t.followup)){
      df.sslt.cur <- df %>%
        dplyr::filter(timePoint %in% c(t.base, t.followup[i])) %>%
        dplyr::group_by(.dots = group.by.sslt) %>%
        dplyr::summarize(mean.value = mean(value, na.rm = T),
                         median.value = median(value, na.rm = T),
                         std.value = sd(value, na.rm = T),
                         cv.value = (sd(value, na.rm = T)/ mean(value, na.rm = T)),
                         n.scans = length(value))

      # flag outliers
      df.sslt.cur <- df.sslt.cur %>%
        dplyr::group_by(.dots = c("parameter", "which.data")) %>%
        dplyr::mutate(outlier.flag = omit.outliers(cv.value))
      df.sslt.cur$outlier.flag <- is.na(df.sslt.cur$outlier.flag)

      # specify calculation type
      df.sslt.cur$t.span <- gsub(" ", "",paste("t", as.character(t.followup[i]), sep = ""))
      df.sslt.cur$s.span <- "single"
      df.sslt <- bind_rows(df.sslt, df.sslt.cur)
    }

  }, silent = TRUE)

  # multi center single time
  df.msst <- NULL
  try({
    group.by.msst <- group.by[!(group.by %in% "site")]
    df.msst <- df %>%
      dplyr::group_by(.dots = group.by.msst) %>%
      dplyr::summarize(mean.value = mean(value, na.rm = T),
                median.value = median(value, na.rm = T),
                std.value = sd(value, na.rm = T),
                cv.value = (sd(value, na.rm = T)/ mean(value, na.rm = T)),
                n.scans = length(value))

    # flag outliers
    df.msst <- df.msst %>%
      dplyr::group_by(.dots = c("parameter", "which.data")) %>%
      dplyr::mutate(outlier.flag = omit.outliers(cv.value))
    df.msst$outlier.flag <- is.na(df.msst$outlier.flag)

    # specify calculation type
    df.msst$t.span <- gsub(" ", "",paste("t", as.character(t.base), sep = ""))
    df.msst$s.span <- "multi"
  }, silent = TRUE)

  # multi center long time
  df.mslt <- NULL
  try({
    group.by.mslt <- group.by[!(group.by %in% c("site", "timePoint"))]
    for (i in 1:length(t.followup)){
    df.mslt.cur <- df %>%
      dplyr::filter(timePoint %in% c(t.base, t.followup[i])) %>%
      dplyr::group_by(.dots = group.by.mslt) %>%
      dplyr::summarize(mean.value = mean(value, na.rm = T),
                median.value = median(value, na.rm = T),
                std.value = sd(value, na.rm = T),
                cv.value = (sd(value, na.rm = T)/ mean(value, na.rm = T)),
                n.scans = length(value))

    # flag outliers
    df.mslt.cur <- df.mslt.cur %>%
      dplyr::group_by(.dots = c("parameter", "which.data")) %>%
      dplyr::mutate(outlier.flag = omit.outliers(cv.value))
    df.mslt.cur$outlier.flag <- is.na(df.mslt.cur$outlier.flag)

    # specify calculation type
    df.mslt.cur$t.span <- gsub(" ", "",paste("t", as.character(t.followup[i]), sep = ""))
    df.mslt.cur$s.span <- "multi"
    df.mslt <- bind_rows(df.mslt, df.mslt.cur)
    }
  }, silent = TRUE)

  # combine dataframes
  df.unpooled <- bind_rows(df.ssst,
                           df.sslt,
                           df.msst,
                           df.mslt)

  # cast precision types as factors, and ensure correct ordering.
  df.unpooled$precision.type <- paste(df.unpooled$t.span, ".", df.unpooled$s.span, sep = "")
  p.types.ordered <- unique(df.unpooled$precision.type )
  df.unpooled$precision.type <- factor(df.unpooled$precision.type, levels = p.types.ordered)


  # flag outliers
  n.outliers <- sum(df.unpooled$outlier.flag )
  if (verbose){
    if (n.outliers > 0) warning(paste("\nWarning: ", n.outliers, " outliers detected \n", sep = ""))
  }


  # estimate rms-statistics
  df.pooled <- df.unpooled %>%
    dplyr::group_by(parameter, precision.type, which.data) %>%
    dplyr::summarize(rms.std = signif(sqrt(mean(std.value^2)), n.signif),
                     median.std = signif(median(std.value), n.signif),
                     iqr.std = paste(signif(quantile(std.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                     rms.cv = signif(sqrt(mean(cv.value^2)), n.signif),
                     median.cv = signif(median(cv.value), n.signif),
                     iqr.cv = paste(signif(quantile(cv.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                     n.scans = sum(n.scans),
                     n.independent = length(cv.value))

  df.pooled.no.outliers <- df.unpooled %>%
    dplyr::filter(outlier.flag == FALSE) %>%
    dplyr::group_by(parameter, precision.type, which.data) %>%
    dplyr::summarize(rms.std = signif(sqrt(mean(std.value^2)), n.signif),
                     median.std = signif(median(std.value), n.signif),
                     iqr.std = paste(signif(quantile(std.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                     rms.cv = signif(sqrt(mean(cv.value^2)), n.signif),
                     median.cv = signif(median(cv.value), n.signif),
                     iqr.cv = paste(signif(quantile(cv.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                     n.scans = sum(n.scans),
                     n.independent = length(cv.value))

  results.list <- list(
    replicate.statistics = df.unpooled,
    rms.statistics = list(
      results = df.pooled,
      grouping.variable = c("parameter", "precision.type", which.data),
      notes = "all data included"),
    rms.statistics.no.outliers = list(
      results = df.pooled.no.outliers,
      grouping.variable = c("parameter", "precision.type", which.data),
      notes = "outliers omitted"
    )
  )


  # assign names
  analysis.name <- paste("mvp.analysis", ".", which.data, sep = "")
  existing.analyses <- object@assays[[which.assay]]@analysis
  existing.analysis.names <- names(existing.analyses)

  if (verbose){
    if (analysis.name %in% existing.analysis.names){
      warning(paste("Pre-existing '", analysis.name, "' Analysis in '", which.assay , "' Assay was overwritten", sep = ""))
    } else {
      cat("\n")
      cat(paste("Analysis saved as '", analysis.name, "'", sep = ""))
      cat("\n")}
  }


  existing.analyses[[analysis.name]] <- results.list
  object@assays[[which.assay]]@analysis <- existing.analyses

  # return calibration object
  return(object)


}




