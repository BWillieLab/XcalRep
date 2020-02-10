#' Compare Calibration Equations
#'
#' Generates boxplot summarizing results generated during svp.analysis. Additionally, Kruskal-Wallis rank sum test is performed (using group.by variable for data stratification) and p-values are shown for each comparison group.
#'
#' @param object Calibration Object.
#' @param which.assay.1 First Assay for calibraiton comparison. Assumed to have fitted calibration equations.
#' @param which.assay.2 Second Assay for calibraiton comparison. Assumed to have fitted calibration equations.
#' @param which.parameter Character specifying which parameters to plot. If unspecified, one parameter is expected in the dataset.
#' @param which.plot Character specifying which plot to generate. One of:
#' \itemize{
#' \item scatter - Default. Pairwise scatter plot between slopes and intercepts in both assays. If scatter, spearman rho is calculated and reported.
#' \item box - Boxplot of slopes and intercepts for each assay. If box,  Kruskal-Wallis rank sum test is conducted and p-value is reported.
#' }
#' @param color.option Character indicating which color map option to use (from viridis palette). One of:
#' \itemize{
#' \item "magma" (or "A")
#' \item "inferno" (or "B")
#' \item "plasma" (or "C")
#' \item "viridis" (or "D")
#' \item "cividis" (or "E") - Default
#' }
#' @param color.begin Hue in [0,1] at which the viridis color map begins. Default is 0.
#' @param color.end Hue in [0,1] at which the viridis color map ends. Default is 1.
#' @param point.size Numeric specifying size of scatter plot points. Default is 2.
#' @param point.alpha Numeric [0,1] specifying degree of transparency for points. Default is 1.
#' @param box.alpha Numeric [0,1] specifying degree of transparency for box plots. Default is 0.5.
#' @param slope.scale Range of values [a,b] specifying min (a) and max(a) limits of slope-related plots. Specify as vector c(a,b).
#' @param intercept.scale Range of values [a,b] specifying min (a) and max(a) limits of intercept-related plots. Specify as vector c(a,b).
#' @param show.reference.line Logical indicating whether reference x=y line is shown as dashed line. Only specified if which.plot is scatter. Default is true.
#' @name compareCalibrationPlot
#' @import viridis
#' @seealso \code{\link{calibrationPlot}}, \code{\link{fitCalibration}}
#' @return ggplot object
#'
compareCalibrationPlot <- function(object, which.assay.1, which.assay.2, which.parameter = NULL, which.plot = "scatter",
                                color.option = "cividis",color.begin = 0, color.end = 1, box.alpha = 0.5, point.size = 2, point.alpha = 1,
                                slope.scale = NULL, intercept.scale = NULL, show.reference.line = T){

  # GIGO Handling
  if (!(which.assay.1 %in% getAssay(object, which.assay = "all"))) stop (which.assay.1,  " does not exist")
  if (!(which.assay.2 %in% getAssay(object, which.assay = "all"))) stop (which.assay.2,  " does not exist")


  # retrieve calibration results (as datatables)
  calibration.results.1 <-getResults(object = co,
                                     which.results = "calibration",
                                     which.assay = which.assay.1,
                                     format = 'df')

  calibration.results.2 <-getResults(object = co,
                                     which.results = "calibration",
                                     which.assay = which.assay.2,
                                     format = 'df')

  res.1 <- calibration.results.1[["calibration.equations"]]
  res.2 <- calibration.results.2[["calibration.equations"]]


  if (!is.null(which.parameter)){
    if (!(which.parameter %in% as.character(unique(res.1$parameter)))) stop (which.parameter,  " does not exist in " ,which.assay.1)
    if (!(which.parameter %in% as.character(unique(res.2$parameter)))) stop (which.parameter,  " does not exist in " ,which.assay.2)
  } else {
    which.parameter <- unique(as.character(unique(res.1$parameter)), as.character(unique(res.2$parameter)))
  }

  if (length(which.parameter) > 1) stop("Multiple parameters detected. Only one must be specified.")

  res.1 <- filter.features(res.1, "parameter", which.parameter)
  res.2 <- filter.features(res.2, "parameter", which.parameter)

  which.var <- c("slope", "intercept", "phantom")

  rename.list.1 <- as.list(which.var)
  names(rename.list.1) <- paste(which.assay.1,".", rename.list.1, sep = "")
  cal.eq.1 <- dplyr::rename(res.1, !!!rename.list.1)

  rename.list.2 <- as.list(which.var)
  names(rename.list.2) <- paste(which.assay.2,".", rename.list.2, sep = "")
  cal.eq.2 <- dplyr::rename(res.2, !!!rename.list.2)

  cal.eq <- join(cal.eq.1, cal.eq.2, by = c("parameter", "reference.site", "calibration.site", "reference.time", "calibration.time"))
  cal.eq <- cal.eq[ ,unique(colnames(cal.eq))]



  if (which.plot == "scatter"){
    hex.col <- viridis_pal(option = color.option, begin = color.begin, end = color.end)(2)

    slope.rho <- signif(cor.test(x=cal.eq[,names(rename.list.1)[[1]]],
                        y=cal.eq[, names(rename.list.2)[[1]]],
                        method = 'spearman', exact = F)[["estimate"]][["rho"]], 3)

    plt.slope.zoom.out <- cal.eq %>%
      ggplot(aes(x = get(names(rename.list.1)[[1]]), y =  get(names(rename.list.2)[[1]]))) +
      geom_point(size = point.size, alpha = point.alpha, color = hex.col[1]) +
      xlab(paste(which.assay.1,  "slope")) +
      ylab(paste(which.assay.2,  "slope")) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill=NA)) +
      ggtitle(paste("Slope Comparison\nrho = ", slope.rho, sep = ""))

    if (show.reference.line){
      plt.slope.zoom.out <- plt.slope.zoom.out + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    }

    if (!is.null(slope.scale)){
      plt.slope.zoom.out <-  plt.slope.zoom.out +
        xlim(slope.scale[1], slope.scale[2])+
        ylim(slope.scale[1],slope.scale[2])
    }

    intercept.rho <- signif(cor.test(x=cal.eq[,names(rename.list.1)[[2]]],
                                   y=cal.eq[, names(rename.list.2)[[2]]],
                                   method = 'spearman', exact = F)[["estimate"]][["rho"]], 3)

    plt.intr.zoom.auto <- cal.eq %>%
      ggplot(aes(x =  get(names(rename.list.1)[[2]]), y =  get(names(rename.list.2)[[2]]))) +
      geom_point(size = point.size, alpha = point.alpha, color = hex.col[2]) +
      xlab(paste(which.assay.1,  "intercept"))  +
      ylab(paste(which.assay.2,  "intercept")) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill=NA)) +
      ggtitle(paste("Intercept Comparison\nrho = ", intercept.rho, sep = ""))

    if (show.reference.line){
      plt.intr.zoom.auto <- plt.intr.zoom.auto + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    }

    if (!is.null(intercept.scale)) {
      plt.intr.zoom.auto <- plt.intr.zoom.auto +
        xlim(intercept.scale[1], intercept.scale[2])+
        ylim(intercept.scale[1],intercept.scale[2])
    }

    return(ggpubr::ggarrange(plt.slope.zoom.out, plt.intr.zoom.auto, ncol = 2, nrow = 1))
    # print(plt.intr.zoom.auto)


  } else if (which.plot == "box"){

    df.stat.sub <-  cal.eq %>%pivot_longer(cols = c(names(rename.list.1)[[2]], names(rename.list.2)[[2]]))
    kruskal.aov.p <- signif(kruskal.test(value ~ name, data = df.stat.sub)[["p.value"]], 3)

    plt.box.int <- cal.eq %>%
      pivot_longer(cols = c(names(rename.list.1)[[2]], names(rename.list.2)[[2]])) %>%
      ggplot(aes(x = name, y = value, fill = name)) +
      ylab("Intecept") + xlab("") +
      geom_boxplot(alpha = box.alpha)  + ggtitle(paste("Intercepts\np = ", kruskal.aov.p, sep = "")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 35, hjust = 1),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_fill_viridis_d(aesthetics = "fill", option =color.option, begin = color.begin, end = color.end)

    df.stat.sub <-  cal.eq %>%pivot_longer(cols = c(names(rename.list.1)[[1]], names(rename.list.2)[[1]]))
    kruskal.aov.p <- signif(kruskal.test(value ~ name, data = df.stat.sub)[["p.value"]], 3)

    plt.box.slope <- cal.eq %>%
      pivot_longer(cols = c(names(rename.list.1)[[1]], names(rename.list.2)[[1]])) %>%
      ggplot(aes(x = name, y = value, fill = name)) +
      ylab("Slope") + xlab("") +
      geom_boxplot(alpha = box.alpha) + ggtitle(paste("Slopes\np = ", kruskal.aov.p, sep = "")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 35, hjust = 1),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_fill_viridis_d(aesthetics = "fill", option = color.option, begin = color.begin, end = color.end)

    if (!is.null(slope.scale)) plt.box.slope <- plt.box.slope + ylim(slope.scale[1],slope.scale[2])
    if (!is.null(intercept.scale)) plt.box.int <- plt.box.int + ylim(intercept.scale[1],intercept.scale[2])

    return(ggpubr::ggarrange(plt.box.slope, plt.box.int, ncol = 2, nrow = 1))

  }

}



#' mean-variance relationship
#'
#' Evalute relationship between replicate means and replicate std and cvs
#'
#' @param object Calibration Object
#' @param which.assay Character specifying which assay to identify reference site for.
#' @param which.data Character specifying which data to identify reference site for.
#' @param which.plot Character specifying which type of plot to generate. One of:
#' \itemize{
#' \item "scatter" -  Default. scatter plot of measure values vs. std and cv values
#' \item "box" - boxplot of std and cv values stratified by phantom section
#' }
#' @param which.parameter Character specfiying which parameters to plot. If unspecified, all parameters are plotted.
#' @param outliers Logical specifying whether to include outliers. Default is true.
#' @param return.plt.handle Logical specifying whether to return plot handle. If TRUE, list of plot handles is returned. If FALSE, plots are printed and not returned.
#' @param color.option Character indicating which color map option to use (from viridis palette). One of:
#' \itemize{
#' \item "magma" (or "A")
#' \item "inferno" (or "B")
#' \item "plasma" (or "C")
#' \item "viridis" (or "D")
#' \item "cividis" (or "E") - Default
#' }
#' @param color.begin Hue in [0,1] at which the viridis color map begins
#' @param color.end Hue in [0,1] at which the viridis color map ends
#' @param point.size Numeric specifying size of scatter plot points. Only implemented if which.plot is "scatter". Default is 2.
#' @param point.alpha Numeric [0,1] specifying transparency of points. Only implemented if which.plot is "scatter". Default is 1.
#' @param error.squared Logical indicating whether to evaluate CV and STD (FALSE) or CV^2 and STD^2 (TRUE). Default is false
#' @name meanVarPlot
#' @import gridExtra, ggpubr, viridis
#' @return List of ggplot handles
#'
meanVarPlot <- function(object,  which.assay = NULL, which.data = "uncalibrated", which.plot = "scatter",which.parameter = NULL,
                          outliers = F, return.plt.handle = F,  color.option = "cividis", color.begin = 0, color.end = 1, point.size = 2, point.alpha = 1, error.squared = F) {

  # ensure assay is specified
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% getAssay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  } else {
    which.assay <- getAssay(object)
  }

  # ensure correct plot specified
  if (!(which.plot %in% c("scatter", "box"))) stop("'which.plot' is incorrectly specified. Must be one of 'scatter' or 'box'")

  # ensure data is specified
  if (!is.null(which.data)) {
    stopifnot(class(which.data) == "character")
    if (!(which.data %in% getDatasets(object, which.assay = which.assay))){
      stop (paste(which.data, " does not exist", sep = ""))
    }
  }

  # get data
  df <- object@assays[[which.assay]]@data[[which.data]]

  u.par <- as.vector(unique(df$parameter))

  # ensure paameters are correctly  specified
  if (!is.null(which.parameter)) {
    stopifnot(class(which.parameter) == "character")
    if (!(which.parameter %in% u.par)){
      stop ("'which.parameter' does not exist")
    }
  } else {
    which.parameter <- u.par
  }

  # ensure values are numeric
  df$value <- as.numeric(as.vector(df$value))

  # calculate precision errors
  group.by <- c("site", "timePoint", "section", "parameter", "phantom")
  df.stats <- df %>%
    dplyr::group_by(.dots = group.by) %>%
    dplyr::summarize(mean.value = mean(value, na.rm = T),
                     median.value = median(value, na.rm = T),
                     std.value = sd(value, na.rm = T),
                     cv.value = (sd(value, na.rm = T)/ mean(value, na.rm = T)),
                     n.scans = length(value))

  # square values
  if (error.squared){
    df.stats$std.value <- (df.stats$std.value)^2
    df.stats$cv.value <- (df.stats$cv.value)^2
    plt.labels <- c("CV^2", "STD^2")
  } else {
    plt.labels <- c("CV", "STD")
    }

  # flag outliers
  if (!outliers){
    df.stats <- df.stats %>%
      dplyr::group_by(.dots = c("parameter")) %>%
      dplyr::mutate(outlier.flag = omitOutliers(cv.value))
    df.stats$outlier.flag <- is.na(df.stats$outlier.flag)
  } else {
    df.stats$outlier.flag <- F
  }

  plt.list.cont <- list()
  plt.list.disc.cv <- list()
  plt.list.disc.sd <- list()
  plt.list.disc.val <- list()

  plt.list <- list()

  for (i in 1:length(which.parameter)){

    df.sum <- df.stats %>% dplyr::filter(parameter == which.parameter[i], !outlier.flag)

    mean.std <- mean(df.sum$std.value)
    mean.cv <- mean(df.sum$cv.value)
    scale.ratio <- (mean.std)/(mean.cv)

    df.rms <- df.sum %>%
      summarize(
        rms.std = mean(std.value),
        rms.cv = mean(cv.value)
      )

    hex.col <- viridis_pal(option = color.option, begin = color.begin, end = color.end)(2)

    if (which.plot == "scatter"){
      plt.list.cont[[which.parameter[i]]] <- df.stats %>%
        dplyr::filter(parameter == which.parameter[i], !outlier.flag) %>%
        ggplot(aes(x = mean.value)) +
        geom_smooth(aes(x = mean.value, y = cv.value), method = "lm", color = hex.col[1], fill =hex.col[1]) +
        geom_smooth(aes(x = mean.value, y = mean.cv*std.value/mean.std), method = "lm", color = hex.col[2], fill =hex.col[2]) +
        geom_point(aes(x = mean.value, y = cv.value,  color = hex.col[1]), size = point.size, alpha = point.alpha) +
        geom_point(aes(x = mean.value, y = cv.value), shape = 1, colour = "black", size = point.size, alpha = point.alpha) +
        geom_point(aes(x = mean.value, y = mean.cv*std.value/mean.std, color = hex.col[2]), size = point.size , alpha = point.alpha)  +
        geom_point(aes(x = mean.value, y = mean.cv*std.value/mean.std), shape = 1, colour = "black", size = point.size, alpha = point.alpha) +
        scale_y_continuous(name = paste("Replicate ", plt.labels[1], sep  = ""), sec.axis = sec_axis(~.*scale.ratio, name = paste("Replicate ", plt.labels[2], sep  = ""))) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        xlab("Replicate Mean") +
        ggtitle(which.parameter[i]) +
        scale_color_manual(name = "Error Type", values = hex.col, labels = plt.labels)

      plt.list[[which.parameter[i]]] <-  plt.list.cont[[which.parameter[i]]]
    } else if (which.plot == "box"){

      df.stat.sub <- df.stats %>% dplyr::filter(parameter == which.parameter[i], !outlier.flag)
      kruskal.aov.val.p <- signif(kruskal.test(mean.value ~ section, data = df.stat.sub)[["p.value"]], 3)
      kruskal.aov.cv.p <- signif(kruskal.test(cv.value ~ section, data = df.stat.sub)[["p.value"]], 3)
      kruskal.aov.sd.p <- signif(kruskal.test(std.value ~ section, data = df.stat.sub)[["p.value"]], 3)

      plt.list.disc.val[[which.parameter[i]]] <- df.stats %>%
        dplyr::filter(parameter == which.parameter[i], !outlier.flag) %>%
        ggplot(aes(x =  paste("Section ",as.character(section), sep = ""), y = mean.value, fill = paste("Section ",as.character(section), sep = ""))) +
        geom_boxplot() +
        theme_bw() +
        theme(legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        xlab("Phantom Section") +
        ylab("Replicate Mean")+
        ggtitle(paste(which.parameter[i], ": Mean \n(p=", kruskal.aov.val.p, ", Kruskal-Wallis Test)", sep = "")) +
        scale_fill_viridis(discrete = T, option = color.option, begin = color.begin , end = color.end)

      plt.list.disc.cv[[which.parameter[i]]] <- df.stats %>%
        dplyr::filter(parameter == which.parameter[i], !outlier.flag) %>%
        ggplot(aes(x =  paste("Section ",as.character(section), sep = ""), y = cv.value, fill = paste("Section ",as.character(section), sep = ""))) +
        geom_boxplot() +
        theme_bw() +
        theme(legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        xlab("Phantom Section") +
        ylab(paste("Replicate ", plt.labels[1], sep  = ""))+
        ggtitle(paste(which.parameter[i], ": ", plt.labels[1], "\n(p=", kruskal.aov.cv.p, ", Kruskal-Wallis Test)", sep = "")) +
        scale_fill_viridis(discrete = T, option = color.option, begin = color.begin , end = color.end)

      plt.list.disc.sd[[which.parameter[i]]] <- df.stats %>%
        dplyr::filter(parameter == which.parameter[i], !outlier.flag) %>%
        ggplot(aes(x =  paste("Section ",as.character(section), sep = ""), y = std.value, fill = paste("Section ",as.character(section), sep = ""))) +
        geom_boxplot() +
        theme_bw() +
        theme(legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        xlab("Phantom Section") +
        ylab(paste("Replicate ", plt.labels[2], sep  = ""))+
        ggtitle(paste(which.parameter[i], ": ", plt.labels[2], "\n(p=", kruskal.aov.sd.p, ", Kruskal-Wallis Test)", sep = "")) +
        scale_fill_viridis(discrete = T, option = color.option, begin = color.begin , end = color.end)

      plt.list[[which.parameter[i]]] <- (gridExtra::grid.arrange(plt.list.disc.val[[which.parameter[i]]],
                                                                 plt.list.disc.cv[[which.parameter[i]]],
                                                                 plt.list.disc.sd[[which.parameter[i]]], nrow = 1))
    }

    if (!return.plt.handle){
      if (which.plot == "scatter"){
        print(plt.list[[which.parameter[i]]])
      }else if (which.plot == "box"){
        ggpubr::as_ggplot(plt.list[[which.parameter[i]]])
      }
    }

  }

  if (return.plt.handle) return(plt.list)


}




#' Visualize site rankings according to mean-squared errors
#'
#' Site rankings are visualized
#'
#' @param object Calibration Object
#' @param var2plot Character indicating which ranking meausre to plot. On of:
#' \itemize{
#' \item "rank" - Default
#' \item "mse" - Normalized mean-squared errors (mse) are plotted.
#' }
#' @param which.phantom Specify which phantom to plot rankings for. If unspecified and multiple phantoms exist in dataset, function will return error requesting that specific phantom be specified.
#' @param which.parameter character indicating which parameters to include. If unspecified, all are included.
#' @param group.by Grouping variable
#' @param which.assay Character specifying which assay to identify reference site for.
#' @param which.data Character specifying which data to identify reference site for.
#' @param which.plot One of:
#' \itemize{
#' \item "box" - boxplot of site-specific parameter rankings
#' \item "tile" - Default. heatmap of site-specific parameter rankings
#' }
#' @param color.option Character indicating which color map option to use (from viridis palette). One of:
#' \itemize{
#' \item "magma" (or "A")
#' \item "inferno" (or "B")
#' \item "plasma" (or "C")
#' \item "viridis" (or "D")
#' \item "cividis" (or "E") - Default
#' }
#' @param show.tile.value Logical specifying whether to show site rankings, overlaid on heamap (if which.plot == "tile")
#' @name consistencyPlot
#' @seealso \code{\link{identifyReference}}, \code{\link{consistencyAnalysis}}
#' @return Character
#'
consistencyPlot <- function(object,  var2plot = "rank", which.phantom = NULL, which.parameter = NULL, group.by = NULL, which.assay = NULL, which.data = "uncalibrated", which.plot = "box", color.option = "cividis", show.tile.value = F) {

  # ensure assay is specified
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% getAssay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  } else {
    which.assay <- getAssay(object)
  }

  # ensure data is specified
  if (!is.null(which.data)) {
    stopifnot(class(which.data) == "character")
    if (!(which.data %in% getDatasets(object, which.assay = which.assay))){
      stop ("'which.data' does not exist")
    }
  }

  # get data
  df <- object@assays[[which.assay]]@data[[which.data]]

  # filter parameters and phantoms
  df <- filter.features(df, "parameter", which.parameter)
  df <- filter.features(df, "phantom", which.phantom)

  if (is.null(which.phantom)){
    if (length(unique(df$phantom)) > 1) stop ("Multiple phantom provided. Specify only one.")
  }

  # compute median values
  df.Med <- df %>%
    dplyr::group_by(phantom, parameter, section) %>%
    dplyr::summarize(median.value = median(value))

  # merge datasets
  df <- merge(df,df.Med, by = c("phantom", "parameter", "section"))

  # compute parameter-specific mean squared errors
  if (is.null(group.by)) group.by <- NA
  if (group.by %in% colnames(df)){
    grp.by.v1 <- c("phantom", "parameter", "site", group.by)
    grp.by.v2 <- c("phantom",  "site", group.by)
  } else {
    grp.by.v1 <- c("phantom", "parameter", "site")
    grp.by.v2 <- c("phantom", "site")
  }

  df.MSE <- df %>%
    dplyr::group_by_at((grp.by.v1)) %>%
    dplyr::summarize(mse = signif(sum((median.value - value)^2)/length(value), 3))

  # rank sites
  df.rank <- df.MSE %>%
    dplyr::group_by_at(c("phantom", "parameter")) %>%
    dplyr::mutate(mse.per = signif(mse / max(mse), 3)) %>%
    dplyr::mutate(site.rank = rank(mse.per))

  # rank best
  df.rank.best <- df.rank %>%
    dplyr::group_by_at(grp.by.v2) %>%
    dplyr::summarize(mean.rank = mean(site.rank)) %>%
    dplyr::arrange(desc(mean.rank))

  df.rank$site <- factor(df.rank$site, levels = as.vector(df.rank.best$site))



  # variable to plot
  if (var2plot == "rank"){
    fill.var <- "site.rank"
    legend.label = "Rank"
    # create tile plot
    c.breaks <- ceiling((seq(1, dim(df.rank.best)[1], by = dim(df.rank.best)[1]/4)))
  } else if (var2plot == "mse"){
    fill.var <- "mse.per"
    legend.label = "MSE (norm.)"
    # create tile plot
    c.breaks <- seq(0,1, by = 0.25)
  }



  if (show.tile.value){
    plt.tile <- ggplot(df.rank, aes(x = parameter, y = site, fill = get(fill.var))) +
      geom_tile() +
      geom_text(aes(label =get(fill.var))) +
      scale_fill_viridis_c(option = color.option, direction = -1, breaks = c.breaks) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 35, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      ylab("Site") +
      xlab("Parameter") +
      ggtitle("Parameter-Stratified Site Rankings") + labs(fill = legend.label)
  } else {
    plt.tile <- ggplot(df.rank, aes(x = parameter, y = site, fill = get(fill.var))) +
      geom_tile() +
      scale_fill_viridis_c(option = color.option, direction = -1, breaks = c.breaks) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 35, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      ylab("Site") +
      xlab("Parameter") +
      ggtitle("Parameter-Stratified Site Rankings") + labs(fill = legend.label)
  }

  # create boxplot
  plt.box <- ggplot(df.rank) +
    geom_boxplot(aes(x = site, y = as.numeric(get(fill.var)), fill = site)) +
    ylab(legend.label) +
    xlab("Site") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_fill_viridis_d(aesthetics = "fill", option = "cividis") +
    ggtitle("Overall Site Ranking")


  if (!is.na(group.by)){
    plt.box <- plt.box + facet_wrap(~get(group.by))
    plt.tile <- plt.tile + facet_wrap(~get(group.by))
  }

  if (which.plot == "box"){
    return(plt.box)
  } else if (which.plot == "tile"){
    return(plt.tile)
  }
}




#' Display table
#'
#' Create asthetically-pleasing tables in R markdown.
#'
#' @param input.table Data frame or data table.
#' @param head.flag Logical specifying whether to return first 10 entries of table. Only specified for data.frame input.
#' @param as.dt Logical specifying whether to cast data.frame as data.table for interactive exploration of data contents. Only works if input is data.frame.
#' @seealso \code{\link{getResults}}
#' @name showTable
#'
showTable <- function(input.table, head.flag = F, as.dt = F){

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
      if (as.dt){
        input.table <- datatable(head(input.table),
                                 filter="top",
                                 width = "100%",
                                 height=  "auto",
                                 extensions = c('Buttons'),
                                 options = list(pageLength = 10,
                                                autoWidth = TRUE,
                                                dom = 'Bfrtip',
                                                buttons = c('copy', 'csv', 'excel')))
        htmltools::tagList(
          input.table
        )
      } else {
        head(input.table)  %>% kable %>% kable_styling()
      }
    } else {
      if (as.dt){
        input.table <- datatable((input.table),
                                 filter="top",
                                 width = "100%",
                                 height=  "auto",
                                 extensions = c('Buttons'),
                                 options = list(pageLength = 10,
                                                autoWidth = TRUE,
                                                dom = 'Bfrtip',
                                                buttons = c('copy', 'csv', 'excel')))
        htmltools::tagList(
          input.table
        )
      } else {
        input.table %>% kable %>% kable_styling()
      }
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
#' @param color.option Character indicating which color map option to use (from viridis palette). One of:
#' \itemize{
#' \item "magma" (or "A")
#' \item "inferno" (or "B")
#' \item "plasma" (or "C")
#' \item "viridis" (or "D")
#' \item "cividis" (or "E") - Default
#' }
#' @param color.begin Hue in [0,1] at which the viridis color map begins. Default is 0.
#' @param color.end Hue in [0,1] at which the viridis color map ends. Default is 0.
#' @param point.size Numeric specifying size of scatter plot points. Only implemented if which.plot is "scatter". Default is 2.
#' @name calibrationPlot
#' @return plot handle (see return.plt.handle argument)
#' @import viridis
#' @seealso \code{\link{fitCalibration}}
calibrationPlot <- function(object, which.assay = NULL, which.parameter = NULL, which.time = NULL, return.plt.handle = F,
                             color.option = "cividis", color.begin = 0, color.end = 0, point.size = 2) {

  #GIGO handling
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% getAssay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  }
  if (is.null(which.assay)) which.assay <- getAssay(object)

  # check that calibraiton curves exist
  if ((!"calibration.curves" %in% names(object@assays[[which.assay]]@calibration))) stop ("calibration curves do not exist")
  calibration.curve.plt <- object@assays[[which.assay]]@calibration[["calibration.curves"]]

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

  # apply color options
  for (i in 1:length(calibration.curve.plt)){
    calibration.curve.plt[[i]] <- calibration.curve.plt[[i]] +
      scale_fill_viridis(discrete = T, option = color.option, begin = color.begin, end = color.end) +
      scale_colour_viridis(discrete = T, option = color.option, begin = color.begin, end = color.end) +
      geom_point(size = point.size)
  }

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
#' Generates boxplot summarizing results generated during svp.analysis. Additionally, Kruskal-Wallis rank sum test is performed (using group.by variable for data stratification) and p-values are shown for each comparison group.
#'
#' @param object Calibration Object.
#' @param which.data Character specifying which dataset to plot svp analysis for. One of:
#' \itemize{
#' \item uncalibrated - Default
#' \item calibrated
#' }
#' @param group.by Vector of features to stratify analysis by. If specified, rms-statistic is not overlaid on boxplot.
#' @param which.parameter Character specifying which parameters to plot. If unspecified, all parameters are plotted.
#' @param which.phantom Character specifying which phantoms to plot. If unspecified, all available phantoms are plotted.
#' @param var2plot Character specifying which precision error metric to plot. One of:
#' \itemize{
#' \item cv - coefficient of variance
#' \item std - standard deviation
#' }
#' @param which.assay Character specifying which assay to plot.
#' @param outliers Logical specifying whether to omit outliers from analysis
#' @param jitter.width Numerical specifying width of jitter in boxplot. Default is 0.1.
#' @param color.option Character indicating which color map option to use (from viridis palette). One of:
#' \itemize{
#' \item "magma" (or "A")
#' \item "inferno" (or "B")
#' \item "plasma" (or "C")
#' \item "viridis" (or "D")
#' \item "cividis" (or "E") - Default
#' }
#' @param color.begin Hue in [0,1] at which the viridis color map begins. Default is 0.
#' @param color.end Hue in [0,1] at which the viridis color map ends. Default is 1.
#' @param point.size Numeric specifying size of scatter plot points. Default is 2.
#' @param point.alpha Numeric [0,1] specifying degree of transparency for points. Default is 1.
#' @param box.alpha Numeric [0,1] specifying degree of transparency for box plots. Default is 0.5.
#' @param txt.size Numeric specifying size of text
#' @param xlab.direction Numeric specifying direciton (i.e., angle) of x axis labels
#' @param xlab.hjust Horizontal justification of x axis labels; 0 is left-justified, 0.5 is centered, and 1 is right-justified.
#' @param xlab.vjust Vertical justification of x axis labels;
#' @param show.p Logical specifying if p-values for Kruskal-Wallis rank sum test are shown. Default is true.
#' @name svpPlot
#' @import viridis
#' @seealso \code{\link{svpAnalysis}}
#' @return ggplot object
#'
svpPlot <- function(object, which.data = "uncalibrated", group.by = NULL, which.parameter = NULL, which.phantom = NULL, var2plot = "cv", which.assay = NULL,
                     outliers = TRUE, jitter.width = 0.1, color.option = "cividis", color.begin = 0, color.end = 1, point.size = 2, point.alpha = 1,
                    box.alpha = 0.5, txt.size = NULL, xlab.direction = NULL, xlab.hjust = NULL, xlab.vjust = NULL, show.p = T) {

  #GIGO handling
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% getAssay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  } else {
    which.assay <- getAssay(object)
  }

  if (!is.null(which.data)){
    stopifnot(class(which.data) == "character")
    if (!(which.data %in% getDatasets(object))){
      stop ("'which.data' does not exist")
    }
  }

  if (!(var2plot %in% c("cv", "std"))) {
    stop ("'var2plot' is incorrectly defined. Must be 'cv' or 'std'")
  }

  # specify which analysis to plot
  which.analysis <-  getAnalyses(object)
  match.ind <- grepl("svpAnalysis", which.analysis)

  if (sum(match.ind) == 0) stop("No precision statistics available to plot. Must run svp.analysis first.")
  if (sum(match.ind) > 0){
    which.analysis <- which.analysis[match.ind]
    match.ind <- grepl(paste(".", which.data, sep = ""), getAnalyses(object), fixed = T)
    which.analysis <-  getAnalyses(object)[match.ind]
  }

  if (length(which.analysis) > 1) stop("Error encountered specifying analysis to plot. Troubleshooting required.")

  ufeatures <- names(getFeatures(object = object, which.assay = which.assay))
  if (!is.null(group.by)){
    stopifnot(class(group.by) == "character")
    if (length(group.by) != 1) stop("'group.by' must specify one feature")
    if (!(group.by %in% ufeatures)) stop (paste(group.by, " does not exist", sep = ""))

  }

  # define group.by features
  if (is.null(group.by)) {
    group.by <- "pooled"
  } else if (sum(group.by %in% ufeatures)> 0){
    group.by <- ufeatures[ufeatures %in% group.by]
  }

  # get data
  df.replicate <- object@assays[[which.assay]]@analysis[[which.analysis]][["replicate.statistics"]][["results"]]

  # filter parameters & parameters
  u.par <- unique(df.replicate$parameter)
  u.phan <- unique(df.replicate$phantom)
  df.replicate <- filter.features(df.replicate, "parameter", which.parameter)
  df.replicate <- filter.features(df.replicate, "phantom", which.phantom)


  if (outliers == F){
    df.pooled <- object@assays[[which.assay]]@analysis[[which.analysis]][["rms.statistics.no.outliers"]][["results"]]
  } else if (outliers == T){
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
  if (!outliers){
    outlier.filter <- c(F)
  } else {
    outlier.filter <- c(T, F)
  }

  # ensure no NA entries exist for given variable of interest
  na.match <- is.na(df.replicate[ ,var2plot])

  if (sum(na.match) > 0) warning(paste(sum(na.match), "NA values omitted. Likely due to presence of single replicate values."))
  df.replicate <- df.replicate[!na.match, ]

  suppressWarnings({
    if (group.by == "pooled"){

      # if CV, run ANOVA (STD comparions across parameters are not approprirate)
      if ((val == "CV") & (length(unique(df.replicate$parameter)) > 1) & (show.p)){
        df.stat.sub <- df.replicate %>% dplyr::filter( outlier.flag %in% outlier.filter)
        kruskal.aov.p <- signif(kruskal.test(get(var2plot) ~ parameter, data = df.stat.sub)[["p.value"]], 3)
        x.label <- paste("Parameter", "\np=", as.vector(kruskal.aov.p), sep = "")
      } else {
        x.label <- "Parameter"
      }

      plt.precision <- df.replicate %>%
        dplyr::filter(outlier.flag %in% outlier.filter) %>%
        ggplot() +
        geom_boxplot(aes(x = parameter, y = get(var2plot), fill = parameter), outlier.shape = NA, alpha = box.alpha) +
        geom_point(aes(x = parameter, y = get(var2plot), fill = parameter, colour = parameter),
                   position = position_jitter(w = jitter.width, h = 0, seed = 1), show.legend = F, size = point.size , alpha = point.alpha) +
        geom_point(aes(x = parameter, y = get(var2plot)),
                   shape = 1, colour = "black",
                   position = position_jitter(w = jitter.width, h = 0, seed = 1), show.legend = F , size = point.size , alpha = point.alpha) +
        geom_point(data = df.pooled, aes(x = parameter, y = get(pooled.par)),
                   colour = "black",
                   size = point.size*2, show.legend = F ) +
        geom_point(data = df.pooled, aes(x = parameter, y = get(pooled.par)),
                   colour = "black", shape = 1,
                   size = point.size*2, show.legend = F ) +
        ggtitle(paste(val, ": pooled", sep = "")) +
        theme_bw() +
        theme(legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA)) +
        scale_fill_viridis(discrete = T, option = color.option, begin = color.begin, end = color.end) +
        scale_color_viridis(discrete = T, option = color.option, begin = color.begin, end = color.end) +
        labs(fill = group.by) +
        xlab(x.label)

    } else {

      df.replicate[,group.by] <- as.factor(as.matrix(df.replicate[,group.by]))
      u.par <- unique(df.replicate$parameter)
      kruskal.aov.p <- c()

      if (show.p){
        for (i in 1:length(u.par)){
          df.stat.sub <- df.replicate %>% dplyr::filter(parameter == u.par[i]) %>% dplyr::filter(outlier.flag %in% outlier.filter)

          u.grp <- as.character(unique(df.stat.sub[ ,group.by]))

          if (length(u.grp) > 1){
            kruskal.aov.p[i] <- signif(kruskal.test(get(var2plot) ~ get(group.by), data = df.stat.sub)[["p.value"]], 3)
            names(kruskal.aov.p)[i] <- u.par[i]
          } else {
            kruskal.aov.p[i] <- 1
            names(kruskal.aov.p)[i] <- u.par[i]
          }
        }
        x.label <-  paste(as.character(u.par), "\np=", as.vector(kruskal.aov.p), sep = "")
      } else {
        x.label <-  as.character(u.par)
      }

      plt.precision <- df.replicate %>%
        dplyr::filter(outlier.flag %in% outlier.filter) %>%
        ggplot() +
        geom_boxplot(aes(x = parameter, y =  get(var2plot), fill = get(group.by)), outlier.shape = NA, alpha = box.alpha) +
        geom_point(aes(x = parameter, y = get(var2plot), fill = get(group.by), colour = get(group.by)),
                   position = position_jitterdodge(jitter.width = jitter.width, jitter.height = 0, seed = 1), show.legend = F , size = point.size , alpha = point.alpha)  +
        geom_point(aes(x = parameter, y = get(var2plot), fill = get(group.by)),
                   shape = 1, colour = "black",
                   position = position_jitterdodge(jitter.width = jitter.width, jitter.height = 0, seed = 1), show.legend = F , size = point.size , alpha = point.alpha) +
        ggtitle(paste(val, ": ", group.by, "-specific", sep = "")) +
        scale_fill_discrete(group.by) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA)) +
        scale_fill_viridis(discrete = T, option = color.option, begin = color.begin, end = color.end) +
        scale_color_viridis(discrete = T, option = color.option, begin = color.begin, end = color.end) +
        labs(fill = group.by) +
        scale_x_discrete("Parameter", labels = x.label, breaks = u.par)
    }
    plt.precision <- plt.precision + ylab(val)
    if (val == "STD") plt.precision <- plt.precision + facet_wrap(~parameter, scales="free")

    # adjust x axis labels
    if (!is.null(txt.size)) plt.precision <- plt.precision + theme(text = element_text(size=txt.size))
    if (!is.null(xlab.direction)) plt.precision <- plt.precision + theme(axis.text.x = element_text(angle = xlab.direction))
    if (!is.null(xlab.hjust)) plt.precision <- plt.precision + theme(axis.text.x = element_text(hjust=xlab.hjust))
    if (!is.null(xlab.vjust)) plt.precision <- plt.precision + theme(axis.text.x = element_text(vjust = xlab.vjust))

  })
  return(plt.precision)
}


#' Diagnostic plot
#'
#' Generates diagnostic curves for evaluation of calibration. Require that data have been calibrated.
#'
#' Diagnostics are only run for data originating from the calibrating phantom. I.e.,if data from multiple phantoms are present within the dataset, only the data from the calibrating phantom is retained for diagnostic purposes.
#'
#' @param object Calibration Object
#' @param which.assay Character specifying which assay to use.
#' @param which.parameter Character indicating which parameters to plot diagnostics for. If unspecified, all parameters are plotted.
#' @param which.plot Character indicating which diagnostic plot to produce. One of:
#' \itemize{
#' \item "line" - Line plot used to compare pre and post calibrated values.
#' \item "bar" - Horizontal barplots used to compare pre and post calibrated values
#' \item "residual" - Histogram used to visualize pre- and post-calibration residuals (i.e., difference between overall mean and replicate mean)
#' }
#' @param color.option Character indicating which color map option to use (from viridis palette). One of:
#' \itemize{
#' \item "magma" (or "A")
#' \item "inferno" (or "B")
#' \item "plasma" (or "C")
#' \item "viridis" (or "D")
#' \item "cividis" (or "E") - Default
#' }
#' @param color.begin Hue in [0,1] at which the viridis color map begins
#' @param color.end Hue in [0,1] at which the viridis color map ends
#' @param highlight.site Character specifying site to highlight in diagnostic plot (valid for line plots)
#' @param bar.errors Logical indicating whether to show error bars for barplot (valid for bar plots)
#' @param fix.axis Logical indicating whether to fix x-axis. Default is True. (Valid for bar and residual plots)
#' @param return.plt.handle Logical specifying whether to return list of plot handles. If TRUE, list of plot handles is returned. If FALSE, plots are printed without returning handle.
#' @import viridis
#' @name diagnosticPlot
#' @return plot handle (see return.plt.handle argument)
#' @seealso \code{\link{fitCalibration}}, \code{\link{calibrateData}}
diagnosticPlot <- function(object, which.assay = NULL, which.parameter = NULL, which.plot = "bar", highlight.site = NULL,
                            color.option = "cividis", color.begin = 0, color.end = 1, bar.errors = F, fix.axis = T, return.plt.handle = F) {

  #GIGO handling

  # ensure assay exists
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% getAssay(object, which.assay = "all"))){
      stop (paste(which.assay, " does not exist", sep = ""))
    }
  } else {
    which.assay <- getAssay(object)
  }

  existing.data <- object@assays[[which.assay]]@data
  existing.calibration <- object@assays[[which.assay]]@calibration
  if (!("calibrated" %in% names(existing.data))) stop("calibrated data does not exist")
  if (!("uncalibrated" %in% names(existing.data))) stop("uncalibrated data does not exist")
  if (!("calibration.equations" %in% names(existing.calibration))) stop("calibration fit does not exist")

  # get data
  df.uncal <- object@assays[[which.assay]]@data[["uncalibrated"]]
  df.cal <- object@assays[[which.assay]]@data[["calibrated"]]
  df.fit <- object@assays[[which.assay]]@calibration[["calibration.equations"]]

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

  # only keep values that match calibration phantom (i.e., omit other phantoms from diagnostics)
  calibration.phantom <- as.character(unique(df.merge$cal.phantom))
  df.merge <- filter.features(df.merge, "cal.phantom", calibration.phantom)
  df.merge <- filter.features(df.merge, "phantom", calibration.phantom)
  na.match <- is.na(df.merge$cal.phantom)
  df.merge <- df.merge[!na.match, ]

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
  if (is.null(which.parameter)){
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


    if (which.plot == "residual"){
      #######################

      suppressWarnings({
      u.time <- unique(df.merge.sub$timePoint)
      df.merge.sub$timePoint <- factor(df.merge.sub$timePoint, levels = u.time[order(u.time)])
      df.merge.sub$section <- paste("Section ", df.merge.sub$section, sep = "")

      df.merge.sub.residual <- df.merge.sub %>%
        dplyr::group_by(section, timePoint, parameter, phantom) %>%
        dplyr::mutate(res.x = x - median(x),
                         res.y = y - median(y))

      x.min <- min(df.merge.sub.residual$res.x, df.merge.sub.residual$res.y)
      x.max <- max(df.merge.sub.residual$res.x, df.merge.sub.residual$res.y) *1.05

      plt.hist.pre <- df.merge.sub.residual %>%
        ggplot(aes(x = res.x)) +
        geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                       colour="black") +
        geom_density(alpha=.2, fill="#FF6666") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        facet_wrap(~section, scales = "free_x") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA))+
        scale_fill_viridis(discrete = T, option = color.option, begin = 0*color.begin, end = 0.5*color.end) +
        scale_color_viridis(discrete = T, option = color.option, begin = 0*color.begin, end = 0.5*color.end) +
        ggtitle(paste(current.parameter, ": Pre-Calibration Residuals", sep = "")) +
        xlab("Pre-Calibration Residuals (Replicate Mean - Pooled Median)") +
        ylab("Density")

      plt.hist.post <- df.merge.sub.residual %>%
        ggplot(aes(x = res.y)) +
        geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                       colour="black") +
        geom_density(alpha=.2, fill="#FF6666") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        facet_wrap(~section, scales = "free_x") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA))+
        scale_fill_viridis(discrete = T, option = color.option, begin = 0*color.begin, end = 0.5*color.end) +
        scale_color_viridis(discrete = T, option = color.option, begin = 0*color.begin, end = 0.5*color.end) +
        ggtitle(paste(current.parameter, ": Post-Calibration Residuals", sep = "")) +
        xlab("Post-Calibration Residuals (Replicate Mean - Pooled Median)") +
        ylab("Density")

      if (fix.axis)  {
        plt.hist.pre <-plt.hist.pre + xlim(x.min, x.max)
        plt.hist.post <-plt.hist.post + xlim(x.min, x.max)
      }

      plt.calibration.list[[u.par[i]]] <- list(plt.hist.pre, plt.hist.post)
})

    }else if (which.plot == "bar"){

      suppressWarnings({
      u.time <- unique(df.merge.sub$timePoint)
      df.merge.sub$timePoint <- factor(df.merge.sub$timePoint, levels = u.time[order(u.time)])

      df.merge.sub$section <- paste("Section ", df.merge.sub$section, sep = "")

      df.merge.sub.sum <- df.merge.sub %>%
        dplyr::group_by(site, timePoint, section, parameter, phantom) %>%
        dplyr::summarize(mean.x = mean(x, na.rm = T),
                         std.x = sd(x, na.rm = T),
                         mean.y = mean(y, na.rm = T),
                         std.y = sd(y, na.rm = T))

      x.min <- min(df.merge.sub.sum$mean.x, df.merge.sub.sum$mean.y)
      x.max <- max(df.merge.sub.sum$mean.x, df.merge.sub.sum$mean.y) * 1.1
      if (x.min > 0) x.min <- 0

      plt.bar.pre <- df.merge.sub.sum %>%
        ggplot(aes(x = site, y = mean.x, fill = timePoint, ymin=mean.x-std.x, ymax=mean.x+std.x)) +
        geom_bar(position = "dodge", stat = "identity", color = "black", alpha = 0.7) +
        geom_hline(data = aggregate(df.merge.sub.sum[, "mean.x"], by = df.merge.sub.sum[, c("section")], FUN =  median),
                   mapping = aes(yintercept = mean.x), linetype = "dashed", color = "red") +
        coord_flip() +
        facet_wrap(~section, scales = "free_y", ncol = 1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.position="bottom")+
        scale_fill_viridis(discrete = T, option = color.option, begin = 0*color.begin, end = 0.5*color.end) +
        scale_color_viridis(discrete = T, option = color.option, begin = 0*color.begin, end = 0.5*color.end) +
        ylab(paste(current.parameter, ": Pre-Calibration Replicate Means", sep = "")) +
        xlab("Site") +
        ggtitle(paste(current.parameter, ": Pre-Calibration", sep = ""))

      if (bar.errors)  plt.bar.pre <- plt.bar.pre + geom_errorbar(color = "black", width = 0.5, position = position_dodge(0.9))

      plt.bar.post <- df.merge.sub.sum %>%
        ggplot(aes(x = site, y = mean.y, fill = timePoint, ymin=mean.y-std.y, ymax=mean.y+std.y)) +
        geom_bar(position = "dodge", stat = "identity", color = "black", alpha = 0.7) +

        geom_hline(data = aggregate(df.merge.sub.sum[, "mean.y"], by = df.merge.sub.sum[, c("section")], FUN =  median),
                   mapping = aes(yintercept = mean.y), linetype = "dashed", color = "red") +
        coord_flip() +
        facet_wrap(~section, scales = "free_y", ncol = 1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.position="bottom")+
        scale_fill_viridis(discrete = T, option = color.option, begin = 0*color.begin, end = 0.5*color.end) +
        scale_color_viridis(discrete = T, option = color.option, begin = 0*color.begin, end = 0.5*color.end) +
        ylab(paste(current.parameter, ": Post-Calibration Replicate Means", sep = "")) +
        xlab("Site") +
        ggtitle(paste(current.parameter, ": Post-Calibration", sep = ""))

      if (bar.errors)  plt.bar.post <- plt.bar.post +  geom_errorbar(color = "black", width = 0.5, position = position_dodge(0.9))

      if (fix.axis)  {
        plt.bar.pre <-plt.bar.pre + ylim(x.min, x.max)
        plt.bar.post <-plt.bar.post + ylim(x.min, x.max)
      }

      plt.calibration.list[[u.par[i]]] <- list(plt.bar.pre, plt.bar.post)

      })

    } else if (which.plot == "line"){

      u.time <- unique(df.merge.sub$timePoint)
      u.time.Out <- paste("TimePoint ",u.time, sep = "")
      df.merge.sub$timePoint <- paste("TimePoint ", df.merge.sub$timePoint, sep = "")
      df.merge.sub$timePoint <- factor(df.merge.sub$timePoint, levels = u.time.Out[order(u.time)])

      if (!is.null(highlight.site)){
        col1 <- "grey20"
        col2 <- "tomato"
        df.merge.sub$highlight <- col1
        df.merge.sub$highlight[grepl(paste(highlight.site, collapse = "|"), df.merge.sub$site)] <- col2

        df.sub1 <- df.merge.sub[df.merge.sub$highlight ==  col1, ]
        df.sub2 <- df.merge.sub[df.merge.sub$highlight ==  col2, ]

        if (length(highlight.site) == 1){
          gtitle <- paste(current.parameter, ": Pre- vs. Post-Calibration (", highlight.site, ")", sep = "")
        } else {
          gtitle <- paste(current.parameter, ": Pre- vs. Post-calibration", sep = "")
        }

        plt.calibration <- df.merge.sub %>%
          ggplot(aes(x, y, color = site, by = timePoint)) +
          geom_point(color = df.merge.sub$highlight) +
          geom_smooth(data = df.sub1, aes(x, y,  group = interaction(site, timePoint)), method='lm',formula=y~x, color =  col1, alpha=.5) +
          geom_smooth(data = df.sub2, aes(x, y, group = interaction(site, timePoint)), method='lm',formula=y~x, color =  col2) +
          ggtitle("calibration") +
          geom_abline(slope=1, intercept=0, linetype = "dashed") +
          xlab(paste("uncalibrated ", current.parameter, sep = "")) +
          ylab(paste("calibrated ", current.parameter, sep = "")) +
          ggtitle(gtitle) +
          geom_hline(yintercept = reference.means, linetype = "dashed", alpha=.5) +
          theme_bw() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.border = element_rect(colour = "black", fill=NA))+
          facet_wrap(~timePoint)

      } else {
        plt.calibration <- df.merge.sub %>%
          ggplot(aes(x, y, colour = site, by = timePoint)) +
          geom_point() +
          geom_smooth(method='lm',formula=y~x) + ggtitle("calibration") +
          geom_abline(slope=1, intercept=0, linetype = "dashed") +
          xlab(paste("uncalibrated ", current.parameter, sep = "")) +
          ylab(paste("calibrated ", current.parameter, sep = "")) +
          ggtitle(paste(current.parameter, ": Pre- vs. Post-Calibration", sep = "")) +
          geom_hline(yintercept = reference.means, linetype = "dashed", alpha=.5) +
          theme_bw() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.border = element_rect(colour = "black", fill=NA)) +
          scale_fill_viridis(discrete = T, option = color.option, begin = color.begin, end = color.end) +
          scale_color_viridis(discrete = T, option = color.option, begin = color.begin, end = color.end) +
          facet_wrap(~timePoint)
      }
      plt.calibration.list[[u.par[i]]] <- plt.calibration

    }
  }

  if (return.plt.handle){
    return(plt.calibration.list)
  } else {

    suppressWarnings({
    for (i in 1:length(plt.calibration.list)){
      if (which.plot == "line"){
        print(plt.calibration.list[[i]])
      }else if ((which.plot == "bar") ){
        gridExtra::grid.arrange(plt.calibration.list[[i]][[1]], plt.calibration.list[[i]][[2]], ncol = 2)
      } else if ((which.plot == "residual")){
        gridExtra::grid.arrange(plt.calibration.list[[i]][[1]], plt.calibration.list[[i]][[2]], ncol = 1)
      }
    }
})
  }

}



#' Multi-variant reproducibility boxplots
#'
#' Generates boxplots summarizing results generated during mvp.analysis.
#'
#' @param object Calibration Object.
#' @param which.data Character specifying which dataset to plot svp analysis for. One of:
#' \itemize{
#' \item "all" - Default
#' \item "uncalibrated"
#' \item "calibrated"
#' }
#' @param group.by Vector of features to stratify analysis by. If specified, rms-statistic is not overlaid on boxplot.
#' @param which.parameter Character specifying which parameters to plot. If unspecified, all parameters are plotted.
#' @param which.phantom Character specifying which phantoms to plot. If unspecified, all phantoms are plotted.
#' @param which.precision Character specifying which precision errors to plot. One of:
#' \itemize{
#' \item "all" - Default. Short and longitudinal, single- and multi-site precision errors
#' \item "short" - Short-term precision errors only
#' \item "long" - Longitudinal precision errors only
#' \item "single" - Single-site precision errors only
#' \item "multi" - Multi-site precision errors only
#' }
#' @param which.analysis Character indicating which mvp.analysis to get results from. If unspecified, available analyses will be used.
#' @param var2plot Character specifying which precision error metric to plot. One of:
#' \itemize{
#' \item cv - coefficient of variance
#' \item std - standard deviation
#' }
#' @param which.assay Character specifying which assay to plot.
#' @param outliers Logical specifying whether to omit outliers from analysis
#' @param jitter.width Numerical specifying width of jitter in boxplot. Default is 0.1.
#' @param color.option Character indicating which color map option to use (from viridis palette). One of:
#' \itemize{
#' \item "magma" (or "A")
#' \item "inferno" (or "B")
#' \item "plasma" (or "C")
#' \item "viridis" (or "D")
#' \item "cividis" (or "E") - Default
#' }
#' @param color.begin Hue in [0,1] at which the viridis color map begins. Default is 0.
#' @param color.end Hue in [0,1] at which the viridis color map ends. Default is 1.
#' @param point.size Numeric specifying size of scatter plot points.
#' @param point.alpha Numeric [0,1] specifying degree of transparency for points. Default is 1.
#' @param box.alpha Numeric [0,1] specifying degree of transparency for box plots. Default is 0.5.
#' @param return.plt.handle Logical indicating whether list of ggplot handle is returned. Default is false.
#' @param combine.plots Logical indicating whether plots are combined. Default is false.
#' @param combine.ncols Number of columns when combining plots (Only if combine.plot is true). Default is 2.
#' @param show.rms.statistic Logical indicating to overlay RMS statistic (rms-std or rms-cv) on boxplots. Default is false.
#' @param show.legend.flag Logical indicating to show legend
#' @name mvpPlot
#' @import viridis
#' @seealso \code{\link{mvpAnalysis}}
#' @return list of ggplot objects (if return.plt.handle is true)
#'
mvpPlot <- function(object, which.data = "all", group.by = NULL, which.parameter = NULL, which.phantom = NULL, which.precision = "all", which.analysis = NULL, var2plot = "cv", which.assay = NULL,
                     outliers = T, jitter.width = 0.1, color.option = "cividis", color.begin = 0, color.end = 1, point.size = 2,
                point.alpha = 1, box.alpha = 0.5, return.plt.handle = F, combine.plots = F, combine.ncols = 2, show.rms.statistic = F, show.legend.flag = T) {

    #GIGO handling
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% getAssay(object, which.assay = "all"))){
      stop (paste(which.assay , " does not exist", sep = ""))
    }
  } else {
    which.assay <- getAssay(object)
  }

  # specify which variable to plot
  if (!(var2plot %in% c("cv", "std"))) {
    stop ("'var2plot' is incorrectly defined. Must be 'cv' or 'std'")
  }

  # specify which precision errors to plot
  if (!(which.precision %in% c("all", "short", "long", "multi", "single"))) {
    stop ("'which.precision' is incorrectly defined.")
  }


  # specify which analysis to plot
  if (is.null(which.analysis)){
    which.analysis <-  getAnalyses(object, which.assay = which.assay)
    match.ind <- grepl("mvpAnalysis", which.analysis)
    if (sum(match.ind) == 0) stop("No precision statistics available to plot. Must run mvp.analysis first.")
    if (sum(match.ind) > 0){
      which.analysis <- which.analysis[match.ind]
    }
    if (length(which.analysis) > 1) stop("Multiple mvp analyses exists. Specify 'which.analysis' to select which to plot.")
  } else {
    if (!(which.analysis %in% getAnalyses(object))) stop(paste(which.analysis, " does not exist", sep = ""))
  }

  # omit outliers if specified
  if (!outliers){
    df.pooled <- object@assays[[which.assay]]@analysis[[which.analysis]][["rms.statistics.no.outliers"]][["results"]]
    df.unpooled <- object@assays[[which.assay]]@analysis[[which.analysis]][["replicate.statistics"]]
    df.unpooled <- df.unpooled[!df.unpooled$outlier.flag, ]
  } else {
    df.pooled <- object@assays[[which.assay]]@analysis[[which.analysis]][["rms.statistics"]][["results"]]
    df.unpooled <- object@assays[[which.assay]]@analysis[[which.analysis]][["replicate.statistics"]]
    df.unpooled$outlier.flag <- F
  }


  #filter phantoms and parameter
  df.unpooled <- filter.features(df.unpooled, "phantom", which.phantom)
  df.pooled <- filter.features(df.pooled, "phantom", which.phantom)

  df.unpooled <- filter.features(df.unpooled, "parameter", which.parameter)
  df.pooled <- filter.features(df.pooled, "parameter", which.parameter)

  # filter data by precision type
  if (which.precision != "all"){
    df.pooled <- df.pooled[grepl(which.precision, df.pooled$precision.type), ]
    df.unpooled <- df.unpooled[grepl(which.precision, df.unpooled$precision.type), ]
  }
  # if (!is.null(which.parameter)){
  #   df.pooled <- df.pooled[grepl(which.parameter, df.pooled$parameter), ]
  #   df.unpooled <- df.unpooled[grepl(which.parameter, df.unpooled$parameter), ]
  # }

  # get unique parametes
  u.par <- as.character(unique(df.pooled$parameter))

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

  # define group.by features
  ufeatures <- names(getFeatures(object = object, which.assay = which.assay))
  if (!is.null(group.by)){
    stopifnot(class(group.by) == "character")
    if (length(group.by) != 1) stop("'group.by' must specify one feature")
    if (!(group.by %in% ufeatures)) stop (paste(group.by, " does not exist", sep = ""))
  }

  if (is.null(group.by)) {
    group.by <- "pooled"
  } else if (sum(group.by %in% ufeatures)> 0){
    group.by <- ufeatures[ufeatures %in% group.by]
  }

  # specify which data to evaluate
  if (which.data != "all"){
    u.data <- as.character(unique(df.unpooled$which.data))
    if(!(which.data %in% u.data)) {
      stop (paste(which.data, " does not exist", sep = ""))
    } else {
      df.pooled <- df.pooled[df.pooled$which.data %in% which.data,  ]
      df.unpooled <- df.unpooled[df.unpooled$which.data %in% which.data,  ]
    }
  }

  plt.mvp.list <- list()
  for (i in 1:length(u.par)){

    current.parameter <- u.par[i]
    df.pooled_subset <- df.pooled %>% dplyr::filter(parameter == current.parameter)
    df.unpooled_subset <- df.unpooled %>% dplyr::filter(parameter == current.parameter)

    u.prec.types <- as.character(unique(df.pooled_subset$precision.type))

    # if (group.by == "pooled"){

      # if CV, run ANOVA (STD comparions across parameters are not approprirate)
      if ((length(unique(df.pooled_subset$which.data)) > 1) & (group.by == "pooled")){

        show.legend.flag <- T
        x.tick <- c()
        for (j in 1:length(u.prec.types)){
          df.stat.sub <- df.pooled_subset %>% dplyr::filter( precision.type %in% u.prec.types[j])

          df.stat.unlist <- NULL
          for (k in 1:nrow(df.stat.sub)){
            values <- unlist(df.stat.sub[k, var2plot])
            df.stat.unlist <- bind_rows(df.stat.unlist, data.frame(which.data = df.stat.sub$which.data[k],
                                                                   values = values))
          }

          u.prec.types.split <- strsplit(u.prec.types[j], "[.]")
          # u.prec.types.split[[1]][1]

          kruskal.aov.p <- signif(kruskal.test(values ~ which.data, data = df.stat.unlist)[["p.value"]], 3)
          x.tick[j] <- paste(u.prec.types.split[[1]][1], "\n", u.prec.types.split[[1]][2], "\np=", as.vector(kruskal.aov.p), sep = "")
        }

      } else {

        show.legend.flag <- F
        x.tick <- c()

        for (j in 1:length(u.prec.types)){
          u.prec.types.split <- strsplit(u.prec.types[j], "[.]")
          x.tick[j] <- paste(u.prec.types.split[[1]][1], "\n", u.prec.types.split[[1]][2], sep = "")
        }
      }

      plt.mvp <- df.unpooled_subset %>%
        ggplot() +
        geom_boxplot(aes(x = precision.type, y = get(var2plot), fill = which.data), outlier.shape = NA, alpha = box.alpha) +
        geom_point(aes(x = precision.type, y =  get(var2plot), fill = which.data, colour = which.data),
                   position = position_jitterdodge(jitter.width = jitter.width, jitter.height = 0, seed = 1), show.legend = F , alpha = point.alpha,  size = point.size) +
        geom_point(aes(x = precision.type, y =  get(var2plot), fill = which.data),
                   shape = 1, colour = "black",
                   position = position_jitterdodge(jitter.width = jitter.width, jitter.height = 0, seed = 1), show.legend = F , alpha = point.alpha,  size = point.size) +
        ggtitle(paste(current.parameter, sep = "")) +
        theme_bw() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.border = element_rect(colour = "black", fill=NA)) +
        ylab(val) +
        xlab("Precision Type") +
        scale_x_discrete(breaks=u.prec.types,
                         labels=x.tick) +
          scale_fill_viridis(discrete = T, option = color.option, begin = color.begin, end = color.end, name = "Calibration Status") +
          scale_color_viridis(discrete = T, option = color.option, begin = color.begin, end = color.end)

      if (show.rms.statistic){
        plt.mvp <- plt.mvp +
          geom_point(data = df.pooled_subset, aes(x = precision.type, y =get(pooled.par), fill = which.data),
                                                colour = "black",
                                                position = position_jitterdodge(jitter.width = 0, jitter.height = 0, seed = 1),
                                                size = 2*point.size, show.legend = F ) +
          geom_point(data = df.pooled_subset, aes(x = precision.type, y =get(pooled.par), fill = which.data),
                     position = position_jitterdodge(jitter.width = 0, jitter.height = 0, seed = 1),
                     size = 2*point.size,  shape = 1, colour = "black", show.legend = F )
      }

      # }

    # else {
    #
    #   }

  if (!(group.by == "pooled")){
    plt.mvp <- plt.mvp + facet_wrap(~get(group.by))
  }

    # if only one strata, remove legend
    if (!show.legend.flag) plt.mvp <- plt.mvp + theme(legend.position = "none")

    plt.mvp.list[[current.parameter]] <- plt.mvp

    if ((!return.plt.handle) & (!combine.plots)) print(plt.mvp.list[[current.parameter]])

    } # end loop

  if ((!return.plt.handle) & (combine.plots)) {
    do.call("grid.arrange", c(plt.mvp.list, ncol=combine.ncols))
  }

  if (return.plt.handle) return(plt.mvp.list)

} # end function







