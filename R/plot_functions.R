
#' Short-term reproducibility boxplot
#'
#' Generates boxplot summarizing results from svp.analysis.
#'
#' @param object calibration object
#' @param var2plot precision error to plot (cv.value or std.value)
#' @param stratify.by features to stratify analysis by.
#' @param which.analysis which analysis to plot (character)
#' @param which.assay which assay to plot
#' @param show.outliers logical specifying whether to omit outliers from analysis
#' @param save.interactive.pooled.plot
#' @param jitter.width width of jitter in boxplot
#' @name svp.boxplot
#' @return calibration object
#'
svp.boxplot <- function(object, var2plot = "cv.value", stratify.by = NULL, which.analysis = NULL, which.assay = NULL, show.outliers = FALSE, save.interactive.pooled.plot = FALSE, jitter.width = 0.1)
{

  # TODO
  # 1. update interactive plots
  # 2. show.outliers (merge df.outliers with df.results)

  #GIGO handling
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% get.assay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  }
  if (!is.null(which.analysis)){
    stopifnot(class(which.analysis) == "character")
    if (!(which.analysis %in% get.analyses(object))){
      stop ("'which.analysis' does not exist")
    }
  }
  if (!(var2plot %in% c("cv.value", "std.value"))) {
    stop ("'var2plot' is incorrectly defined. Must be 'cv.value' or 'std.value'")
  }

  # if null arguments (set to current assay, and first available "precisionAnalysis" in list)
  if (is.null(which.assay)) which.assay <- get.assay(object)
  if (is.null(which.analysis)){
    match.ind <- grepl("svp.analysis", get.analyses(object))
    if (sum(match.ind) == 0) stop("No precision statistics available to plot")
    if (sum(match.ind) > 0){
      which.analysis <-  get.analyses(object)[match.ind[1]]
    }
  }

  ufeatures <- names(get.features(object = object, which.assay = which.assay))
  if (!is.null(stratify.by)){
    stopifnot(class(stratify.by) == "character")
    if (length(stratify.by) != 1) stop("'stratify.by' must specify one feature")
    if (!(stratify.by %in% ufeatures)) stop ("'stratify.by' does not exist")

  }

  # define group.by features
  if (is.null(stratify.by)) {
    group.by <- "pooled"
  } else if (sum(stratify.by %in% ufeatures)> 0){
    group.by <- ufeatures[ufeatures %in% stratify.by]
  }

  # get data
  df.replicate <- object@assays[[which.assay]]@analysis[[which.analysis]][["replicate.statistics"]][["results"]]
  df.pooled <- object@assays[[which.assay]]@analysis[[which.analysis]][["pooled.statistics"]][["results"]]

  # set datatype
  if (var2plot == "cv.value"){
    val <- "CV"
    pooled.par <- "rms.cv"
  } else if (var2plot == "std.value"){
    val <- "STD"
    pooled.par <- "rms.std"
  }

  # if (save.interactive.plot){
  # if more than one site and timepoint, include labels( for interactive plots)
  # REVISIT THIS
  additional_labels_of_interest <- c("site", "timePoint")
  N.counts <- unlist(object@assays[[current.assay]]@N)
  N.counts <- N.counts[N.counts > 1]
  pattern <- paste(additional_labels_of_interest, collapse = "|")
  labels_to_add <- sapply(strsplit(names(N.counts)[grepl(pattern, names(N.counts))], "[.]"), "[[", 2)
  # }
  suppressWarnings({
    if (group.by == "pooled"){


      plt.precision <- df.replicate %>%
        ggplot() +
        geom_boxplot(aes(x = parameter, y = get(var2plot), fill = parameter), outlier.shape = NA) +
        geom_point(aes(x = parameter, y = get(var2plot), fill = parameter, colour = parameter),
                   position = position_jitter(w = jitter.width, h = 0, seed = 1), show.legend = F ) +
        geom_point(aes(x = parameter, y = get(var2plot),
                       text = paste(
                         "Mean: ", mean.value,
                         "\nStd: ", std.value,
                         "\nCV: ", cv.value,
                         "\nn: ", n.replicates)),
                   shape = 1, colour = "black",
                   position = position_jitter(w = jitter.width, h = 0, seed = 1), show.legend = F ) +
        geom_point(data = df.pooled, aes(x = parameter, y = get(pooled.par)),
                   colour = "black",
                   size = 4, show.legend = F ) +
        ggtitle(paste(val, ": pooled", sep = "")) +
        theme_bw() +
        theme(legend.position = "none",
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"))


    } else {

      df.replicate[,group.by] <- as.factor(as.matrix(df.replicate[,group.by]))

      plt.precision <- df.replicate %>%
        ggplot() +
        geom_boxplot(aes(x = parameter, y =  get(var2plot), fill = get(group.by)), outlier.shape = NA) +
        geom_point(aes(x = parameter, y = get(var2plot), fill = get(group.by), colour = get(group.by)),
                   position = position_jitterdodge(jitter.width = jitter.width, jitter.height = 0, seed = 1), show.legend = F )  +
        geom_point(aes(x = parameter, y = get(var2plot), fill = get(group.by)),
                   shape = 1, colour = "black",
                   position = position_jitterdodge(jitter.width = jitter.width, jitter.height = 0, seed = 1), show.legend = F ) +
        ggtitle(paste(val, ": ", group.by, "-specific", sep = "")) +
        scale_fill_discrete(group.by) +
        theme_bw() +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"))
    }
  })

  if (val == "CV"){
    # define y limits
    ylim.scale.factor <- 1.2
    ylim1 <- max(as.matrix(df.replicate[ , var2plot]))
    plt.precision <- plt.precision + ylim(0, ylim1*ylim.scale.factor)
  }

  plt.precision <- plt.precision + ylab(val)


  if (val == "STD") plt.precision <- plt.precision + facet_wrap(~parameter, scales="free")

  # save plot handle
  plt.name <- paste("svp.boxplot.", val, ".", which.analysis,".", group.by, sep = "")
  existing.plots <- object@assays[[which.assay]]@plots
  existing.plot.names <- names(existing.plots)

  if (plt.name %in% existing.plot.names){
    warning(paste("Pre-existing'", plt.name, "' in '", which.assay , "' was overwritten", sep = ""))
  } else {
    cat("\n============================\n")
    cat(paste("plot(s) saved as '", plt.name, "'", sep = ""))
    cat("\n")}


  existing.plots[[plt.name]] <- plt.precision
  object@assays[[which.assay]]@plots <- existing.plots

  # generate interactive html plot
  if (save.interactive.pooled.plot & (group.by == "pooled")){
    # file.name <- paste(which.analysis , "_", group.by, "_", format(Sys.time(), "_%b%d_%Hh_%Mm_%Y"), ".html" , sep = "")
    file.name <- paste(plt.name, "_", format(Sys.time(), "_%b%d_%Hh_%Mm_%Y"), ".html" , sep = "")
    p <- ggplotly(plt.precision, tooltip = "text")
    htmlwidgets::saveWidget(p, file.path(getwd(), file.name))
  }

  # print plot
  print(plt.precision)

  return(object)
}




#' Short-term reproducibility violin plot
#'
#' Generates violin plot summarizing results from svp.analysis.
#'
#' @param object calibration object
#' @param var2plot precision error to plot (cv.value or std.value)
#' @param stratify.by features to stratify analysis by.
#' @param which.analysis which analysis to plot (character)
#' @param which.assay which assay to plot
#' @param show.outliers logical specifying whether to omit outliers from analysis
#' @param save.interactive.pooled.plot
#' @param jitter.width width of jitter in boxplot
#' @name svp.boxplot
#' @return calibration object
#'
svp.violinplot=function(object, var2plot = "cv.value", stratify.by = NULL, which.analysis = NULL, which.assay = NULL, show.outliers = FALSE, save.interactive.pooled.plot = FALSE, jitter.width = 0.1) {


  # TODO
  # 1. update interactive plots
  # 2. show.outliers (merge df.outliers with df.results)

  #GIGO handling
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% get.assay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  }
  if (!is.null(which.analysis)){
    stopifnot(class(which.analysis) == "character")
    if (!(which.analysis %in% get.analyses(object))){
      stop ("'which.analysis' does not exist")
    }
  }
  if (!(var2plot %in% c("cv.value", "std.value"))) {
    stop ("'var2plot' is incorrectly defined. Must be 'cv.value' or 'std.value'")
  }

  # if null arguments (set to current assay, and first available "precisionAnalysis" in list)
  if (is.null(which.assay)) which.assay <- get.assay(object)
  if (is.null(which.analysis)){
    match.ind <- grepl("svp.analysis", get.analyses(object))
    if (sum(match.ind) == 0) stop("No precision statistics available to plot")
    if (sum(match.ind) > 0){
      which.analysis <-  get.analyses(object)[match.ind[1]]
    }
  }

  ufeatures <- names(get.features(object = object, which.assay = which.assay))
  if (!is.null(stratify.by)){
    stopifnot(class(stratify.by) == "character")
    if (length(stratify.by) != 1) stop("'stratify.by' must specify one feature")
    if (!(stratify.by %in% ufeatures)) stop ("'stratify.by' does not exist")

  }

  # define group.by features
  if (is.null(stratify.by)) {
    group.by <- "pooled"
  } else if (sum(stratify.by %in% ufeatures)> 0){
    group.by <- ufeatures[ufeatures %in% stratify.by]
  }

  # get data
  df.replicate <- object@assays[[which.assay]]@analysis[[which.analysis]][["replicate.statistics"]][["results"]]
  df.pooled <- object@assays[[which.assay]]@analysis[[which.analysis]][["pooled.statistics"]][["results"]]

  # set datatype
  if (var2plot == "cv.value"){
    val <- "CV"
    pooled.par <- "rms.cv"
  } else if (var2plot == "std.value"){
    val <- "STD"
    pooled.par <- "rms.std"
  }



  # if (save.interactive.plot){
  # if more than one site and timepoint, include labels( for interactive plots)
  # REVISIT THIS
  additional_labels_of_interest <- c("site", "timePoint")
  N.counts <- unlist(object@assays[[current.assay]]@N)
  N.counts <- N.counts[N.counts > 1]
  pattern <- paste(additional_labels_of_interest, collapse = "|")
  labels_to_add <- sapply(strsplit(names(N.counts)[grepl(pattern, names(N.counts))], "[.]"), "[[", 2)
  # }
  suppressWarnings({
    if (group.by == "pooled"){


      plt.precision <- df.replicate %>%
        ggplot() +
        geom_violin(aes(x = parameter, y = get(var2plot), fill = parameter), outlier.shape = NA) +
        geom_point(aes(x = parameter, y = get(var2plot), fill = parameter, colour = parameter),
                   position = position_jitter(w = jitter.width, h = 0, seed = 1), show.legend = F ) +
        geom_point(aes(x = parameter, y = get(var2plot),
                       text = paste(
                         "Mean: ", mean.value,
                         "\nStd: ", std.value,
                         "\nCV: ", cv.value,
                         "\nn: ", n.replicates)),
                   shape = 1, colour = "black",
                   position = position_jitter(w = jitter.width, h = 0, seed = 1), show.legend = F ) +
        geom_point(data = df.pooled, aes(x = parameter, y = get(pooled.par)),
                   colour = "black",
                   size = 4, show.legend = F ) +
        ggtitle(paste(val, ": pooled", sep = "")) +
        theme_bw() +
        theme(legend.position = "none",
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"))


    } else {

      df.replicate[,group.by] <- as.factor(as.matrix(df.replicate[,group.by]))

      plt.precision <- df.replicate %>%
        ggplot() +
        geom_violin(aes(x = parameter, y =  get(var2plot), fill = get(group.by)), outlier.shape = NA) +
        geom_point(aes(x = parameter, y = get(var2plot), fill = get(group.by), colour = get(group.by)),
                   position = position_jitterdodge(jitter.width = jitter.width, jitter.height = 0, seed = 1), show.legend = F )  +
        geom_point(aes(x = parameter, y = get(var2plot), fill = get(group.by)),
                   shape = 1, colour = "black",
                   position = position_jitterdodge(jitter.width = jitter.width, jitter.height = 0, seed = 1), show.legend = F ) +
        ggtitle(paste(val, ": ", group.by, "-specific", sep = "")) +
        scale_fill_discrete(group.by) +
        theme_bw() +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"))
    }
  })

  if (val == "CV"){
    # define y limits
    ylim.scale.factor <- 1.2
    ylim1 <- max(as.matrix(df.replicate[ , var2plot]))
    plt.precision <- plt.precision + ylim(0, ylim1*ylim.scale.factor)
  }

  plt.precision <- plt.precision + ylab(val)

  if (val == "STD") plt.precision <- plt.precision + facet_wrap(~parameter, scales="free")

  # save plot handle
  plt.name <- paste("svp.violinplot.", val, ".", which.analysis,".", group.by, sep = "")
  existing.plots <- object@assays[[which.assay]]@plots
  existing.plot.names <- names(existing.plots)

  if (plt.name %in% existing.plot.names){
    warning(paste("Pre-existing'", plt.name, "' in '", which.assay , "' was overwritten", sep = ""))
  } else {
    cat("\n============================\n")
    cat(paste("plot(s) saved as '", plt.name, "'", sep = ""))
    cat("\n")}


  existing.plots[[plt.name]] <- plt.precision
  object@assays[[which.assay]]@plots <- existing.plots

  # generate interactive html plot
  if (save.interactive.pooled.plot & (group.by == "pooled")){
    # file.name <- paste(which.analysis , "_", group.by, "_", format(Sys.time(), "_%b%d_%Hh_%Mm_%Y"), ".html" , sep = "")
    file.name <- paste(plt.name, "_", format(Sys.time(), "_%b%d_%Hh_%Mm_%Y"), ".html" , sep = "")
    p <- ggplotly(plt.precision, tooltip = "text")
    htmlwidgets::saveWidget(p, file.path(getwd(), file.name))
  }

  # print plot
  print(plt.precision)

  return(object)
}


#' Diagnostic plot
#'
#' Generates diagnostic curves for evaluation of calibration. Require that data have been calibrated.
#'
#' @param object calibration object
#' @param which.assay which assay to plot
#' @name diagnostic.plot
#' @return calibration object
#'
diagnostic.plot <- function(object, which.assay = NULL) {


  #GIGO handling

  # ensure assay exists
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% get.assay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  }

  # ensure assay is specified
  if (is.null(which.assay)) which.assay <- get.assay(object)

  existing.data <- object@assays[[which.assay]]@data
  existing.calibration <- object@assays[[which.assay]]@calibration
  if (!("calibrated" %in% names(existing.data))) stop("calibrated data does not exist")
  if (!("uncalibrated" %in% names(existing.data))) stop("uncalibrated data does not exist")
  if (!("calibration.fit" %in% names(existing.calibration))) stop("calibration fit does not exist")

  # get data
  df.uncal <- object@assays[[which.assay]]@data[["uncalibrated"]]
  df.cal <- object@assays[[which.assay]]@data[["calibrated"]]
  df.fit <- object@assays[[which.assay]]@calibration[["calibration.fit"]]

  which.uncal.val <- which(colnames(df.uncal) == "value")
  which.cal.val <- which(colnames(df.cal) == "value")

  colnames(df.uncal)[which.uncal.val] <- "x"
  colnames(df.cal)[which.cal.val] <- "y"

  # add unique identifier to each entry to allow merge
  if (nrow(df.uncal) != nrow(df.cal)) stop("uncalibrated and calibrated data are different sizes and incompatible")
  df.uncal$id <- seq(1, nrow(df.uncal))
  df.cal$id <- df.uncal$id

  # join datasets
  df.merge <- suppressMessages({join(df.uncal, df.cal)})
  df.merge <- df.merge[, c(names(df.uncal)[(seq(2, ncol(df.uncal)-2))], "x", "y")]

  # ensure grouping variables exist
  if (!("site" %in% colnames(df.merge))) stop("site data is missing")
  if ("parameter" %in% colnames(df.merge)) {
    u.par <- unique(df.merge$parameter)
  } else {
    df.merge$parameter <- "parameter"
    u.par <- unique(df.merge$parameter)
  }

  if ("timePoint" %in% colnames(df.merge)) {
    u.time <- unique(df.merge$timePoint)
  } else {
    df.merge$timePoint <- 0
    u.time <- unique(df.merge$timePoint)
  }

  plt.calibration.list <- list()

  for (i in 1:length(u.par)){

    # get current parameter
    current.parameter <- u.par[i]
    current.time <- u.time[i]

    # subset data
    df.merge.sub <- df.merge %>%
      filter(parameter == current.parameter)

    # if unique reference time, get reference time means
    reference.site <- unique(df.fit$reference.site)
    reference.time <- unique(df.fit$reference.time)
    if (length(reference.time) == 1){
      df.reference <- df.merge %>%
        filter(parameter == current.parameter,
               timePoint == reference.time,
               site == reference.site) %>%
        group_by(section) %>%
        summarize(value.mean = mean(x),
                  value.median = median(x))
      reference.means <- df.reference$value.mean
    } else {reference.means <- NA}

    # plot data
    plt.calibration <- df.merge.sub %>%
      ggplot(aes(x, y, colour = site, by = timePoint)) +
      geom_point() +
      geom_smooth(method='lm',formula=y~x) + ggtitle("calibration") +
      geom_abline(slope=1, intercept=0, linetype = "dashed") +
      xlab(paste("uncalibrated ", current.parameter, sep = "")) +
      ylab(paste("calibrated ", current.parameter, sep = "")) +
      ggtitle(paste("Calibration Curves: ", current.parameter, sep = "")) +
      geom_hline(yintercept = reference.means, linetype = "dashed", alpha=.5) +
      facet_wrap(~timePoint)

    plt.name.current <- paste("pre.post.calibration.diagnostic.", current.parameter, sep = "")

    plt.calibration.list[[plt.name.current]] <- plt.calibration
    print(plt.calibration)
  }


  # save plot handle
  plt.name <- paste("diagnostic.plot", sep = "")
  existing.plots <- object@assays[[which.assay]]@plots
  existing.plot.names <- names(existing.plots)

  if (plt.name %in% existing.plot.names){
    warning(paste("Pre-existing'", plt.name, "' in '", which.assay , "' was overwritten", sep = ""))
  } else {
    cat("\n============================\n")
    cat(paste("'", plt.name, "' created and saved", sep = ""))
    cat("\n")}

  existing.plots[[plt.name]] <- plt.calibration.list
  object@assays[[which.assay]]@plots <- existing.plots

  return(object)
}




