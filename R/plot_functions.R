
#' Display table in markdown
#'
#' create asthetically-pleasing tables in R markdown.
#'
#' @param input.table calibration object
#' @param head.flag Logical specifying whether to return first 10 entries of tabke. Only specified for data.frame input.
#' @name show.table
#'
show.table <- function(input.table, head.flag = F){

  # get class of input.table
  tbl.class <- class(input.table)

  if (any(grepl("data.frame", tbl.class))) {
    tbl.format <- "df"
  } else if (any(grepl("datatables", tbl.class))) {
    input.table <- list(input.table)
    tbl.format <- "dt"
  } else if (grepl("list", tbl.class)) {
    tbl.class <- class(input.table[[1]])
    if (any(grepl("datatables", tbl.class))){
      tbl.format <- "dt"
    } else {
      stop ("input table is not data.frame or data.table")
    }
  } else {
    stop ("input table is not data.frame or data.table")
  }

  if (tbl.format == "df"){
    if (head.flag){
      head(input.table)  %>% kable %>% kable_styling()
    } else {
      input.table %>% kable %>% kable_styling()
    }
  } else if (tbl.format == "dt"){

    htmltools::tagList(
      input.table
    )

  }


}



#' Plot Calibration Curves
#'
#' Pair-wise calibration curves for each parameter at each timePoint are plotted with 95% confidence interval and dashed line reference indicating x=y.
#'
#' Currently intercepts and slopes are fitted for all calibration curves, independent of whether intercepts were fit for calibration equation using fit.equation. This will be updated to accomodate no intercept calibration curve plots in future releasaed of XcalRep.
#'
#' @param object Calibration Object.
#' @param which.assay Character specifying which assay to get calibration curves from.
#' @param which.parameter Vector specifying which parameter calibration curves to plot. If unspecified, all parameters plotted.
#' @param which.time Vector specifying which times to plot calibration curves for. If unspecified, all timePoints plotted.
#' @param return.plt.handle Logical specifying whether to return list of plot handles. If TRUE, list of plot handles is returned. If FALSE, plots are printed without returning handle.
#' @name calibration.plot
#' @return plot handle (see return.plt.handle argument)
#' @seealso fit.calibration
calibration.plot <- function(object, which.assay = NULL, which.parameter = NULL, which.time = NULL, return.plt.handle = F) {

  #GIGO handling
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% get.assay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  }
  if (is.null(which.assay)) which.assay <- get.assay(object)

  # check that calibraiton curves exist
  if ((!"calibration.curve.plots" %in% names(object@assays[[which.assay]]@calibration))) stop ("calibration curves do not exist")
  calibration.curve.plt <- object@assays[[which.assay]]@calibration[["calibration.curve.plots"]]

  # get specified plots
  which.entries <- names(calibration.curve.plt)

  if (is.null(which.parameter)) which.parameter <- "all"
  if (is.null(which.time)) which.time <- "all"

  # filter entries based on parameter and time specifications
  if (which.parameter != "all"){
    match.ind <- grepl(paste(which.parameter, collapse = "|"), which.entries)
    if (sum(match.ind) == 0) stop("'which.parameter' is incorrectly specified")
    which.entries <- which.entries[match.ind]
  }

  if (which.time != "all"){
    if (tolower(which.time) == tolower("baseline")){
      df <- object@assays[[which.assay]]@data[["uncalibrated"]]
      which.time <- min(as.matrix((df %>% dplyr::select(timePoint))))
    }
    which.time <- gsub(" ", "", paste("t", as.character(which.time), ""))
    match.ind <- grepl(paste(which.time, collapse = "|"), which.entries)
    if (sum(match.ind) == 0) stop("'which.time' is incorrectly specified")
    which.entries <- which.entries[match.ind]
  }

if (length(which.entries) < 1) stop ("not calibration curves exist for specified parameters and times")


  if (return.plt.handle){

    return(calibration.curve.plt[which.entries])

  } else {
    for (i in 1:length(which.entries)){
      print(calibration.curve.plt[[which.entries[i]]])
    }
  }
}





#' Short-term reproducibility boxplot
#'
#' Generates boxplot summarizing results generated during svp.analysis.
#'
#' @param object Calibration Object.
#' @param which.data Character specifying which dataset to plot svp analysis for. One of:
#' \itemize{
#' \item uncalibrated - Default
#' \item calibrated
#' }
#' @param group.by Vector of features to stratify analysis by. If specified, rms-statistic is not overlaid on boxplot.
#' @param var2plot Character specifying which precision error metric to plot. One of:
#' \itemize{
#' \item cv - coefficient of variance
#' \item std - standard deviation
#' }
#' @param which.assay Character specifying which assay to plot.
#' @param omit.outliers Logical specifying whether to omit outliers from analysis
#' @param jitter.width Numerical specifying width of jitter in boxplot. Default = 0.1.
#' @name svp.plot
#' @seealso svp.analysis
#' @return ggplot object
#'
svp.plot <- function(object, which.data = "uncalibrated", group.by = NULL,  var2plot = "cv", which.assay = NULL,
                        omit.outliers = FALSE, jitter.width = 0.1) {

  #GIGO handling
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% get.assay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  } else {
    which.assay <- get.assay(object)
  }

  if (!is.null(which.data)){
    stopifnot(class(which.data) == "character")
    if (!(which.data %in% get.datasets(object))){
      stop ("'which.data' does not exist")
    }
  }

  if (!(var2plot %in% c("cv", "std"))) {
    stop ("'var2plot' is incorrectly defined. Must be 'cv' or 'std'")
  }

  # specify which analysis to plot
  which.analysis <-  get.analyses(object)
    match.ind <- grepl("svp.analysis", which.analysis)

    if (sum(match.ind) == 0) stop("No precision statistics available to plot. Must run svp.analysis first.")
    if (sum(match.ind) > 0){
      which.analysis <- which.analysis[match.ind]
      match.ind <- grepl(paste(".", which.data, sep = ""), get.analyses(object), fixed = T)
      which.analysis <-  get.analyses(object)[match.ind]
    }

    if (length(which.analysis) > 1) stop("Error encountered specifying analysis to plot. Troubleshooting required.")

  ufeatures <- names(get.features(object = object, which.assay = which.assay))
  if (!is.null(group.by)){
    stopifnot(class(group.by) == "character")
    if (length(group.by) != 1) stop("'group.by' must specify one feature")
    if (!(group.by %in% ufeatures)) stop ("'group.by' does not exist")

  }

  # define group.by features
  if (is.null(group.by)) {
    group.by <- "pooled"
  } else if (sum(group.by %in% ufeatures)> 0){
    group.by <- ufeatures[ufeatures %in% group.by]
  }

  # get data
  df.replicate <- object@assays[[which.assay]]@analysis[[which.analysis]][["replicate.statistics"]][["results"]]

  if (omit.outliers == T){
    df.pooled <- object@assays[[which.assay]]@analysis[[which.analysis]][["rms.statistics.no.outliers"]][["results"]]
  } else if (omit.outliers == F){
    df.pooled <- object@assays[[which.assay]]@analysis[[which.analysis]][["rms.statistics"]][["results"]]
  }

  # set datatype
  if (var2plot == "cv"){
    val <- "CV"
    var2plot <- "cv.value"
    pooled.par <- "rms.cv"
  } else if (var2plot == "std"){
    val <- "STD"
    var2plot <- "std.value"
    pooled.par <- "rms.std"
  }

  # construct filter criteria based on outlier flag
  if (omit.outliers){
    outlier.filter <- c(F)
  } else {
    outlier.filter <- c(T, F)
  }

  suppressWarnings({
    if (group.by == "pooled"){

      plt.precision <- df.replicate %>%
        dplyr::filter(outlier.flag %in% outlier.filter) %>%
        ggplot() +
        geom_boxplot(aes(x = parameter, y = get(var2plot), fill = parameter), outlier.shape = NA) +
        geom_point(aes(x = parameter, y = get(var2plot), fill = parameter, colour = parameter),
                   position = position_jitter(w = jitter.width, h = 0, seed = 1), show.legend = F ) +
        geom_point(aes(x = parameter, y = get(var2plot)),
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
        dplyr::filter(outlier.flag %in% outlier.filter) %>%
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

  plt.precision <- plt.precision + ylab(val)


  if (val == "STD") plt.precision <- plt.precision + facet_wrap(~parameter, scales="free")

    return(plt.precision)

}


#' Diagnostic plot
#'
#' Generates diagnostic curves for evaluation of calibration. Require that data have been calibrated.
#'
#' @param object Calibration Object
#' @param which.assay Character specifying which assay to use.
#' @param which.parameter Vector of parameters to plot diagnostics for. If unspecified, all parameters are plotted.
#' @param highlight.site Character specifying site to highlight in diagnostic plot.
#' @param return.plt.handle Logical specifying whether to return list of plot handles. If TRUE, list of plot handles is returned. If FALSE, plots are printed without returning handle.
#' @name diagnostic.plot
#' @return plot handle (see return.plt.handle argument)
#' @seealso fit.calibration, calibrate.data
diagnostic.plot <- function(object, which.assay = NULL, which.parameter = "all", highlight.site = NULL, return.plt.handle = F) {

  #GIGO handling

  # ensure assay exists
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% get.assay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  } else {
    which.assay <- get.assay(object)
  }

  existing.data <- object@assays[[which.assay]]@data
  existing.calibration <- object@assays[[which.assay]]@calibration
  if (!("calibrated" %in% names(existing.data))) stop("calibrated data does not exist")
  if (!("uncalibrated" %in% names(existing.data))) stop("uncalibrated data does not exist")
  if (!("fit.calibration" %in% names(existing.calibration))) stop("calibration fit does not exist")

  # get data
  df.uncal <- object@assays[[which.assay]]@data[["uncalibrated"]]
  df.cal <- object@assays[[which.assay]]@data[["calibrated"]]
  df.fit <- object@assays[[which.assay]]@calibration[["fit.calibration"]]

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
  # df.merge <- df.merge[, c(names(df.uncal)[(seq(2, ncol(df.uncal)-2))], "x", "y")]

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

  # get parameter list to plot
  if (which.parameter == "all"){
    which.parameter <- u.par
  } else {
    if (sum(which.parameter %in% u.par) == 0) stop("specified parameters do not exist")
    which.parameter <- which.parameter[which.parameter %in% u.par]
  }

  u.par <- which.parameter

  for (i in 1:length(u.par)){

    # get current parameter
    current.parameter <- u.par[i]
    current.time <- u.time[i]

    # subset data
    df.merge.sub <- df.merge %>%
      dplyr::filter(parameter == current.parameter)

    # if unique reference time, get reference time means
    reference.site <- unique(df.fit$reference.site)
    reference.time <- unique(df.fit$reference.time)
    if (length(reference.time) == 1){
      df.reference <- df.merge %>%
        dplyr::filter(parameter == current.parameter,
               timePoint == reference.time,
               site == reference.site) %>%
        dplyr::group_by(section) %>%
        dplyr::summarize(value.mean = mean(x),
                  value.median = median(x))
      reference.means <- df.reference$value.mean
    } else {reference.means <- NA}

    # plot data


    if (!is.null(highlight.site)){
      col1 <- "grey20"
      col2 <- "tomato"
      df.merge.sub$highlight <- col1
      df.merge.sub$highlight[grepl(paste(highlight.site, collapse = "|"), df.merge.sub$site)] <- col2

      df.sub1 <- df.merge.sub[df.merge.sub$highlight ==  col1, ]
      df.sub2 <- df.merge.sub[df.merge.sub$highlight ==  col2, ]

      plt.calibration <- df.merge.sub %>%
        ggplot(aes(x, y, color = site, by = timePoint)) +
        geom_point(color = df.merge.sub$highlight) +
        geom_smooth(data = df.sub1, aes(x, y,  group = interaction(site, timePoint)), method='lm',formula=y~x, color =  col1, alpha=.5) +
        geom_smooth(data = df.sub2, aes(x, y, group = interaction(site, timePoint)), method='lm',formula=y~x, color =  col2) +
        ggtitle("calibration") +
        geom_abline(slope=1, intercept=0, linetype = "dashed") +
        xlab(paste("uncalibrated ", current.parameter, sep = "")) +
        ylab(paste("calibrated ", current.parameter, sep = "")) +
        ggtitle(paste(current.parameter, ": pre- vs. post-calibration", sep = "")) +
        geom_hline(yintercept = reference.means, linetype = "dashed", alpha=.5) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        facet_wrap(~timePoint)
    } else {
      plt.calibration <- df.merge.sub %>%
        ggplot(aes(x, y, colour = site, by = timePoint)) +
        geom_point() +
        geom_smooth(method='lm',formula=y~x) + ggtitle("calibration") +
        geom_abline(slope=1, intercept=0, linetype = "dashed") +
        xlab(paste("uncalibrated ", current.parameter, sep = "")) +
        ylab(paste("calibrated ", current.parameter, sep = "")) +
        ggtitle(paste(current.parameter, ": pre- vs. post-calibration", sep = "")) +
        geom_hline(yintercept = reference.means, linetype = "dashed", alpha=.5) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        facet_wrap(~timePoint)
    }


    plt.name.current <- paste("pre.post.calibration.diagnostic.", current.parameter, sep = "")

    plt.calibration.list[[plt.name.current]] <- plt.calibration

      # print(plt.calibration)


  }

    if (return.plt.handle){

      return(plt.calibration.list)

    } else {
      for (i in 1:length(plt.calibration.list)){
        print(plt.calibration.list[[i]])
      }
    }

}


#' Multi-variant precision boxplots
#'
#' Generates boxplot summarizing results generated during svp.analysis.
#'
#' @param object Calibration Object.
#' @param which.data Character specifying which dataset to plot svp analysis for. One of:
#' \itemize{
#' \item uncalibrated - Default
#' \item calibrated
#' }
#' @param group.by Vector of features to stratify analysis by. If specified, rms-statistic is not overlaid on boxplot.
#' @param var2plot Character specifying which precision error metric to plot. One of:
#' \itemize{
#' \item cv - coefficient of variance
#' \item std - standard deviation
#' }
#' @param which.assay Character specifying which assay to plot.
#' @param omit.outliers Logical specifying whether to omit outliers from analysis
#' @param jitter.width Numerical specifying width of jitter in boxplot. Default = 0.1.
#' @name svp.plot
#' @seealso svp.analysis
#' @return ggplot object
#'
mvp.plot <- function(object, which.data = "uncalibrated", group.by = NULL,  var2plot = "cv", which.assay = NULL,
                     omit.outliers = FALSE, jitter.width = 0.1) {



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
