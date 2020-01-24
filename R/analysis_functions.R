#' Single variant precision/reproducibility analysis
#'
#' Short-term single-center reproducibility/precision analysis of repeated phantom/in vivo measurements.
#' Precisions are calculated as root-mean-square standard deviation (rms.std) and root-ean-square coefficient of variance (rms.cv).
#' Results are saved in calibration object as 'svp.analysis'.
#'
#' @param object Calibration Object
#' @param which.data A character that specifies which data to analysis. One of  (uncalibrated or calibrated)
#' \itemize{
#' \item uncalibrated
#' \item calibrated
#' }
#' @param stratify.by features to stratify analysis by. If left unspecified, default grouping is used (recommended), i.e., site, timePoint, section, parameter, phantom
#' @param n.signif number of significant digits to report
#' @param which.assay specifies assay for analysis
#' @param omit.outliers logical specifying whether to omit outliers from analysis. If TRUE, outliers will be saved as separate dataset within current Assay.
#' @name svp.analysis
#' @return calibration object
#'
svp.analysis <- function(object, which.data,  stratify.by = NULL, n.signif = 3, which.assay = NULL, omit.outliers = FALSE){

  # calibrated data: logical (true/false); if null, checks for calibrated data first.

  #GIGO handling
  if (!is.null(which.assay)) stopifnot(class(which.assay) == "character")
  if (is.null(which.assay)) which.assay <- get.assay(object)
  stopifnot((class(n.signif) == "numeric") & (length(n.signif) == 1))
  if (!(which.data %in% names(object@assays[[which.assay]]@data))) stop(paste("'", which.data , "' data does not exist", sep = ""))

# get grouping features
  ufeatures <- get.features(object = object, which.assay = which.assay)
  ufeatures.names <- names(ufeatures)

  if (is.null(stratify.by)){
      group.by <- c("site", "timePoint", "section", "parameter", "phantom")
  } else if (sum(stratify.by %in% ufeatures.names)> 0){
    group.by <- ufeatures.names[ufeatures.names %in% stratify.by]
  } else {
    stop("'stratify.by' is incorrectly specified")
  }

  # get data
  # if (is.null(calibrated.data))
  df <- object@assays[[which.assay]]@data[[which.data]]

  # ensure values are numeric
  df$value <- as.numeric(as.vector(df$value))

  # calculate precision errors
  df.stats <- df %>%
    group_by(.dots = group.by) %>%
    summarize(mean.value = mean(value, na.rm = T),
              median.value = median(value, na.rm = T),
              std.value = sd(value, na.rm = T),
              cv.value = (sd(value, na.rm = T)/ mean(value, na.rm = T)),
              n.scans = length(value))

  # omit outliers and store in df.outliers for downstream evaluation
  # outliers removed according to cv values (not std values)
  if (omit.outliers){
    outliers.as.na <- as.numeric(omit.outliers(as.matrix(df.stats[ , "cv.value"])))
    df.outliers <- df.stats[is.na(outliers.as.na), ]
    df.stats <- df.stats[!is.na(outliers.as.na), ]
    n.outliers <- nrow(df.outliers)
    if (n.outliers > 0) warning(paste("\nWarning: ", n.outliers, " outliers omitted \n", sep = ""))
  }


  # parameter-level precision errors
  # check if parameter is specified - if not, assume single parameter
  if (("parameter" %in% ufeatures.names) & (length(unique(df.stats$parameter)) > 1)){
    df.stats_pooled <- df.stats %>%
      group_by(parameter) %>%
      summarize(rms.std = signif(sqrt(mean(std.value^2)), n.signif),
                rms.cv = signif(sqrt(mean(cv.value^2)), n.signif),
                n.scans = sum(n.scans),
                n.indepedent = round(length(cv.value)/ length(unique(section))))
  } else {
    # cat("\nAll data pooled, single parameter input assumed \n")
    df.stats_pooled <- df.stats %>%
      summarize(rms.std = signif(sqrt(mean(std.value^2)), n.signif),
                rms.cv = signif(sqrt(mean(cv.value^2)), n.signif),
                n.scans = sum(n.scans),
                n.indepedent = round(length(cv.value)/ length(unique(section))))
  }


  # enforce significant figures
  df.stats_sigfigs <- df.stats
  is.num <- sapply(df.stats_sigfigs, is.numeric)
  df.stats_sigfigs[is.num] <- lapply(df.stats_sigfigs[is.num], signif, n.signif)

  # print results in data table
  # if (show.data.table){
  #   print(datatable(df.stats_pooled, filter = 'top'))
  #   print(datatable(df.stats_sigfigs, filter = 'top'))
  # }

  # save results to calibration object
  if (exists("df.outliers")){
    replicate.statistics = list(
      results = df.stats_sigfigs,
      outliers = df.outliers,
      grouping.variable = group.by)
  } else {
    replicate.statistics = list(
      results = df.stats_sigfigs,
      grouping.variable = group.by)
  }

  results.list <- list(
    replicate.statistics = replicate.statistics,
    rms.statistics = list(
      results = df.stats_pooled,
      grouping.variable = "parameter")
  )

  # save analysis
  if (is.null(stratify.by)){
    analysis.grouping <- ""
  } else {
    analysis.grouping <- paste(".", stratify.by, collapse = "_") # append stratification gorup to end of analysis name
  }

  analysis.name <- paste("svp.analysis", ".", which.data, analysis.grouping, sep = "")
  existing.analyses <- object@assays[[which.assay]]@analysis
  existing.analysis.names <- names(existing.analyses)

  if (analysis.name %in% existing.analysis.names){
    warning(paste("Pre-existing '", analysis.name, "' Analysis in '", which.assay , "' Assay was overwritten", sep = ""))
  } else {
    # cat("\n============================\n")
    cat("\n")
    cat(paste("analysis saved as '", analysis.name, "'", sep = ""))
    cat("\n")}

  existing.analyses[[analysis.name]] <- results.list
  object@assays[[which.assay]]@analysis <- existing.analyses


  return(object)
}

#' Multi variant precision/reproducibility analysis
#'
#' Longitudinal (i.e., multiple time points) and multicenter (i.e., multiple instruments) reproducibility/precision analysis of
#' repeated phantom/in vivo measurements. #' Precisions are calculated as root-mean-square standard deviation (rms.std) and
#' root-mean-square coefficient of variance (rms.cv). Results are saved in calibration object as 'mvp.analysis'.
#'
#' @param object calibration object
#' @param which.data specifies data to analysis (uncalibrated or calibrated)
#' @param var2plot precision error to plot (cv.value or std.value)
#' @param show.data.table logical specifying whether to generate datatable summary of results
#' @param n.signif number of significant digits to report
#' @param which.assay specifies assay for analysis
#' @param trim.plot logical specifying whether to trim outliers from plots
#' @param plot.flag logical specifying whether to plot results
#' @name mvp.analysis
#' @return calibration object
#'
mvp.analysis <- function(object, which.data = "all",var2plot = "cv.value",show.data.table = TRUE, n.signif = 3, which.assay = NULL, trim.plot = FALSE, plot.flag = TRUE)  {

  # calibrated data: logical (true/false); if null, checks for calibrated data first.

  #GIGO handling
  if (!is.null(which.assay)) stopifnot(class(which.assay) == "character")
  if (is.null(which.assay)) which.assay <- get.assay(object)
  stopifnot(class(show.data.table) == "logical")
  stopifnot((class(n.signif) == "numeric") & (length(n.signif) == 1))


  # grouping fe
  ufeatures <- names(get.features(object = object, which.assay = which.assay))
  ufeatures <- ufeatures[!(ufeatures %in% "scanDate")]
  group.by <- ufeatures

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

  # set datatype
  if (var2plot == "cv.value"){
    val <- "CV"
    pooled.par <- "rms.cv"
  } else if (var2plot == "std.value"){
    val <- "STD"
    pooled.par <- "rms.std"
  }

  if (!("parameter" %in% colnames(df))) df$parameter <- "parameter"

  # ensure values are numeric
  df$value <- as.numeric(as.vector(df$value))

  # reorder datasets
  if (("uncalibrated" %in% u.data) & ("uncalibrated" %in% u.data)){
    df$which.data <- factor(df$which.data)
    match.uncal <- match("uncalibrated", levels(df$which.data))
    match.cal <- match("calibrated", levels(df$which.data))
    df$which.data <- factor(df$which.data, levels(df$which.data)[c(match.uncal, match.cal)])
  } else {df$which.data <- factor(df$which.data)}

  # unique timepoints and sites
  u.data <- unique(df$which.data)
  # u.site <- unique(df$site)
  # u.time <- unique(df$timePoint)
  group.by <- c(group.by, "which.data")

  # single center single time
  df.unpooled <- NULL

  df.ssst <- NULL
  try({
    group.by.ssst <- group.by
    df.ssst <- df %>%
      group_by(.dots = group.by.ssst) %>%
      summarize(mean.value = mean(value),
                median.value = median(value),
                std.value = sd(value),
                cv.value = (sd(value)/ mean(value)),
                n.replicates = length(value))
  }, silent = TRUE)

  # single center long time
  df.sslt <- NULL
  try({
    group.by.sslt <- group.by[!(group.by %in% "timePoint")]
    df.sslt <- df %>%
      group_by(.dots = group.by.sslt) %>%
      summarize(mean.value = mean(value),
                median.value = median(value),
                std.value = sd(value),
                cv.value = (sd(value)/ mean(value)),
                n.replicates = length(value))
  }, silent = TRUE)

  # multi center single time
  df.msst <- NULL
  try({
    group.by.msst <- group.by[!(group.by %in% "site")]
    df.msst <- df %>%
      group_by(.dots = group.by.msst) %>%
      summarize(mean.value = mean(value),
                median.value = median(value),
                std.value = sd(value),
                cv.value = (sd(value)/ mean(value)),
                n.replicates = length(value))
  }, silent = TRUE)

  # multi center long time
  df.mslt <- NULL
  try({
    group.by.mslt <- group.by[!(group.by %in% c("site", "timePoint"))]
    df.mslt <- df %>%
      group_by(.dots = group.by.mslt) %>%
      summarize(mean.value = mean(value),
                median.value = median(value),
                std.value = sd(value),
                cv.value = (sd(value)/ mean(value)),
                n.replicates = length(value))
  }, silent = TRUE)

  # assign precision types
  df.ssst$precision.type <- "shortTerm.singleSite"
  df.sslt$precision.type <- "longTerm.singleSite"
  df.msst$precision.type <- "shortTerm.multiSite"
  df.mslt$precision.type <- "longTerm.multiSite"

  # combine dataframes
  df.unpooled <- bind_rows(df.ssst,
                           df.sslt,
                           df.msst,
                           df.mslt)

  # reorder datasets
  df.unpooled$precision.type <- factor(df.unpooled$precision.type)
  match.ssst <- match("shortTerm.singleSite", levels(df.unpooled$precision.type))
  match.sslt <- match("longTerm.singleSite", levels(df.unpooled$precision.type))
  match.msst <- match("shortTerm.multiSite", levels(df.unpooled$precision.type))
  match.mslt <- match("longTerm.multiSite", levels(df.unpooled$precision.type))
  df.unpooled$precision.type <- factor(df.unpooled$precision.type,
                                       levels(df.unpooled$precision.type)[c(match.ssst, match.sslt, match.msst, match.mslt)])

  # pooled results
  df.pooled <- NULL

  df.ssst_pooled <- NULL
  try({
    df.ssst_pooled <- df.ssst %>%
      group_by(parameter, which.data) %>%
      summarize(rms.std = signif(sqrt(mean(std.value^2)), n.signif),
                median.std = signif(median(std.value), n.signif),
                iqr.std = signif(IQR(std.value), n.signif),
                lower.upper.quartile.std = paste(signif(quantile(std.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                rms.cv = signif(sqrt(mean(cv.value^2)), n.signif),
                median.cv = signif(median(cv.value), n.signif),
                iqr.cv = signif(IQR(cv.value), n.signif),
                lower.upper.quartile.cv = paste(signif(quantile(cv.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                n.replicates = sum(n.replicates),
                n.indepedent = length(cv.value))
    df.ssst_pooled$precision.type <- "shortTerm.singleSite"
  }, silent = TRUE)

  df.sslt_pooled <- NULL
  try({
    df.sslt_pooled <- df.sslt %>%
      group_by(parameter, which.data) %>%
      summarize(rms.std = signif(sqrt(mean(std.value^2)), n.signif),
                median.std = signif(median(std.value), n.signif),
                iqr.std = signif(IQR(std.value), n.signif),
                lower.upper.quartile.std = paste(signif(quantile(std.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                rms.cv = signif(sqrt(mean(cv.value^2)), n.signif),
                median.cv = signif(median(cv.value), n.signif),
                iqr.cv = signif(IQR(cv.value), n.signif),
                lower.upper.quartile.cv = paste(signif(quantile(cv.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                n.replicates = sum(n.replicates),
                n.indepedent = length(cv.value))
    df.sslt_pooled$precision.type <- "longTerm.singleSite"
  }, silent = TRUE)

  df.msst_pooled <- NULL
  try({
    df.msst_pooled <- df.msst %>%
      group_by(parameter, which.data) %>%
      summarize(rms.std = signif(sqrt(mean(std.value^2)), n.signif),
                median.std = signif(median(std.value), n.signif),
                iqr.std = signif(IQR(std.value), n.signif),
                lower.upper.quartile.std = paste(signif(quantile(std.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                rms.cv = signif(sqrt(mean(cv.value^2)), n.signif),
                median.cv = signif(median(cv.value), n.signif),
                iqr.cv = signif(IQR(cv.value), n.signif),
                lower.upper.quartile.cv = paste(signif(quantile(cv.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                n.replicates = sum(n.replicates),
                n.indepedent = length(cv.value))
    df.msst_pooled$precision.type <- "shortTerm.multiSite"
  }, silent = TRUE)

  df.mslt_pooled <- NULL
  try({
    df.mslt_pooled <- df.mslt %>%
      group_by(parameter, which.data) %>%
      summarize(rms.std = signif(sqrt(mean(std.value^2)), n.signif),
                median.std = signif(median(std.value), n.signif),
                iqr.std = signif(IQR(std.value), n.signif),
                lower.upper.quartile.std = paste(signif(quantile(std.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                rms.cv = signif(sqrt(mean(cv.value^2)), n.signif),
                median.cv = signif(median(cv.value), n.signif),
                iqr.cv = signif(IQR(cv.value), n.signif),
                lower.upper.quartile.cv = paste(signif(quantile(cv.value, probs = c(0.25, 0.75)), n.signif), collapse = ", "),
                n.replicates = sum(n.replicates),
                n.indepedent = length(cv.value))
    df.mslt_pooled$precision.type <- "longTerm.multiSite"
  }, silent = TRUE)

  df.pooled <- bind_rows(df.ssst_pooled,
                         df.sslt_pooled,
                         df.msst_pooled,
                         df.mslt_pooled)

  u.par <- unique(df.pooled$parameter)

  plt.mvp.list <- list()
  for (i in 1:length(u.par)){

    current.parameter <- u.par[i]

    df.pooled_subset <- df.pooled %>% filter(parameter == current.parameter)
    df.unpooled_subset <- df.unpooled %>% filter(parameter == current.parameter)

    if (trim.plot){

      # df.temp <- df.unpooled_subset[, var2plot]
      # ylim1 = boxplot.stats(as.matrix(df.temp))$stats[c(1, 5)]

      # ylim1 =

      plt.mvp <- df.unpooled_subset %>%
        ggplot(aes(x = precision.type, y = get(var2plot), fill = which.data)) +
        geom_boxplot(outlier.shape = NA) +
        ggtitle(paste(current.parameter, sep = "")) +
        theme_bw() +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black")) +
        ylab(val)
      # +
      # ylim(0, ylim1[2]*1.2)

    } else{

      jitter.width <- 0.1

      plt.mvp <- df.unpooled_subset %>%
        ggplot() +
        geom_boxplot(aes(x = precision.type, y = get(var2plot), fill = which.data), outlier.shape = NA) +
        geom_point(aes(x = precision.type, y =  get(var2plot), fill = which.data, colour = which.data),
                   position = position_jitterdodge(jitter.width = jitter.width, jitter.height = 0, seed = 1), show.legend = F ) +
        geom_point(aes(x = precision.type, y =  get(var2plot), fill = which.data),
                   shape = 1, colour = "black",
                   position = position_jitterdodge(jitter.width = jitter.width, jitter.height = 0, seed = 1), show.legend = F ) +
        geom_point(data = df.pooled_subset, aes(x = precision.type, y =get(pooled.par), fill = which.data),
                   colour = "black",
                   position = position_jitterdodge(jitter.width = 0, jitter.height = 0, seed = 1),
                   size = 4, show.legend = F ) +
        ggtitle(paste(current.parameter, sep = "")) +
        scale_fill_discrete(which.data) +
        theme_bw() +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black")) +
        ylab(val)

    }

    plt.mvp.list[[current.parameter]] <- plt.mvp
    if (plot.flag) print(plt.mvp)

  }


  # save plot handle
  plt.name <- paste("mvp.boxplot.", val, ".", which.data, sep = "")
  existing.plots <- object@assays[[which.assay]]@plots
  existing.plot.names <- names(existing.plots)

  if (plt.name %in% existing.plot.names){
    warning(paste("Pre-existing'", plt.name, "' in '", which.assay , "' was overwritten", sep = ""))
  } else {
    cat("\n============================\n")
    cat(paste("plot(s) saved as '", plt.name, "'", sep = ""))
    cat("\n")}

  existing.plots[[plt.name]] <- plt.mvp.list
  object@assays[[which.assay]]@plots <- existing.plots

  # enforce significant figures
  df.pooled_sigfigs <- df.pooled
  is.num <- sapply(df.pooled_sigfigs, is.numeric)
  df.pooled_sigfigs[is.num] <- lapply(df.pooled_sigfigs[is.num], signif, n.signif)

  df.unpooled_sigfigs <- df.unpooled
  is.num <- sapply(df.unpooled_sigfigs, is.numeric)
  df.unpooled_sigfigs[is.num] <- lapply(df.unpooled_sigfigs[is.num], signif, n.signif)

  # print results in data table
  if (show.data.table){
    print(datatable(df.unpooled_sigfigs, filter = 'top'))
    print(datatable(df.pooled_sigfigs, filter = 'top'))
  }


  # save analysis
  results.list <- list(
    unpooled.statistics = df.unpooled_sigfigs,
    pooled.statistics = df.pooled_sigfigs
  )

  analysis.name <- paste("mvp.analysis.", val, ".", which.data, sep = "")
  existing.analyses <- object@assays[[which.assay]]@analysis
  existing.analysis.names <- names(existing.analyses)

  if (analysis.name %in% existing.analysis.names){
    warning(paste("Pre-existing'", analysis.name, "' in '", which.assay , "' was overwritten", sep = ""))
  } else {
    cat("\n============================\n")
    cat(paste("analysis saved as '", analysis.name, "'", sep = ""))
    cat("\n")}

  existing.analyses[[analysis.name]] <- results.list
  object@assays[[which.assay]]@analysis <- existing.analyses


  return(object)
}




