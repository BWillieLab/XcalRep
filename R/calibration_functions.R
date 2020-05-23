


#' Identify Reference Site
#'
#' Sites are ranked according to mean-squared error calculated between site-specific parameters and overall median parameter values (pooled across all timePoints and sites), and top ranking site is returned.
#'
#' @param object Calibration Object
#' @param which.assay Character specifying which assay to identify reference site for.
#' @param which.parameter Character specifying which parameters to consider when selecting reference site
#' @param which.phantom Character specifying which phantoms to select reference sites for. If unspecified, all phantoms are considered.
#' @param which.data Character specifying which data to identify reference site for.
#' @param pool.phantoms Logical specifying whether to pool ranking scores across all phantoms (true), or to rank phantoms separately (false, default).
#' @param return.scores Logical indicating whether list of scores are returned. If false (default), reference sites are returned only.
#' @name identifyReference
#' @seealso \code{\link{consistencyAnalysis}}, \code{\link{consistencyPlot}}
#' @return If return.scores is false, named vector of character(s) specfying reference site is returned. In case of tie, multiple sites are returned. If return.scores is true, list of site rankings is returned.
#'
identifyReference <- function(object,  which.assay = NULL, which.parameter= NULL, which.phantom= NULL, which.data = "uncalibrated", pool.phantoms = F, return.scores = F) {

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

  # filter parameters
  df <- filterFeatures(df, "parameter", which.parameter)
  df <- filterFeatures(df, "phantom", which.phantom)

  # compute median values
  df.Med <- df %>%
    dplyr::group_by(phantom, parameter, section) %>%
    dplyr::summarize(median.value = median(value, na.rm = T))

  # merge datasets
  df <- merge(df,df.Med, by = c("phantom", "parameter", "section"))

  # compute parameter-specific mean squared errors
  df.MSE <- df %>%
    dplyr::group_by(phantom, parameter, site) %>%
    dplyr::summarize(mse = signif(mean((median.value - value)^2, na.rm = T), 3))


  df.rank.list <- list()
  df.rank.best.list <- list()
  u.phan <- unique(df$phantom)
  reference.site <- c()
  for (i in 1:length(u.phan)){
    df.rank.list[[u.phan[i]]] <- df.MSE %>%
      dplyr::filter(phantom == u.phan[i]) %>%
      dplyr::group_by(phantom, parameter) %>%
      dplyr::mutate(mse.normalized = signif(mse / max(mse), 3)) %>%
      dplyr::mutate(site.rank = rank(mse.normalized))

    df.rank.best.list[[u.phan[i]]] <- df.rank.list[[u.phan[i]]]  %>%
      dplyr::group_by(phantom, site) %>%
      dplyr::summarize(mean.rank = mean(site.rank)) %>%
      dplyr::arrange(desc(mean.rank))

    reference.site[u.phan[i]] <- as.vector(df.rank.best.list[[u.phan[i]]]$site[which(min(df.rank.best.list[[u.phan[i]]]$mean.rank) == df.rank.best.list[[u.phan[i]]]$mean.rank)])
  }

  if (pool.phantoms){
    reference.site <- NULL
    merge.scores <- NULL
    for (i in 1:length(u.phan)){
      merge.scores <- bind_rows(merge.scores, df.rank.list[[u.phan[i]]])
    }
    df.rank.list[["pool"]] <- merge.scores

    df.rank.best.pooled <- merge.scores  %>%
      dplyr::group_by(site) %>%
      dplyr::summarize(mean.rank = mean(site.rank)) %>%
      dplyr::arrange(desc(mean.rank))

    ref.site <- as.vector(df.rank.best.pooled$site[which(min(df.rank.best.pooled$mean.rank) ==df.rank.best.pooled$mean.rank)])
    if (length(ref.site) > 1) {
      reference.site["pooled.score"] <- paste(ref.site, collapse = ", ")
    } else {
      reference.site["pooled.score"] <- ref.site
    }
  }

  if (return.scores){
    return(df.rank.list)
  } else {
    return(reference.site)
  }

}


#' Evaluate site consistency
#'
#' Sites are ranked according to mean-squared error calculated between site-specific parameters and overall median parameter values (pooled across all timePoints and sites).
#'
#' @param object Calibration Object
#' @param which.assay Character specifying which assay to identify reference site for.
#' @param which.data Character specifying which data to identify reference site for.
#' @name consistencyAnalysis
#' @seealso \code{\link{identifyReference}}, \code{\link{consistencyPlot}}
#' @return data frame
#'
consistencyAnalysis <- function(object,  which.assay = NULL, which.data = "uncalibrated") {

  # ensure assay is specified
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% getAssay(object, which.assay = "all"))){
      stop (paste(which.assay, " does not exist", sep = ""))
    }
  } else {
    which.assay <- getAssay(object)
  }

  # ensure data is specified
  if (!is.null(which.data)) {
    stopifnot(class(which.data) == "character")
    if (!(which.data %in% getDatasets(object, which.assay = which.assay))){
      stop(paste(which.data, " does not exist", sep = ""))
    }
  }

  # get data
  df <- object@assays[[which.assay]]@data[[which.data]]

  # compute median values
  df.Med <- df %>%
    dplyr::group_by(phantom, parameter, section) %>%
    dplyr::summarize(median.value = median(value))

  # merge datasets
  df <- merge(df,df.Med, by = c("phantom", "parameter", "section"))

  # compute parameter-specific mean squared errors
  df.MSE <- df %>%
    dplyr::group_by(phantom, parameter, site) %>%
    dplyr::summarize(mse = signif(sum((median.value - value)^2)/length(value), 3))

  # rank sites
  df.rank <- df.MSE %>%
    dplyr::group_by(phantom, parameter) %>%
    dplyr::mutate(mse.normalized = signif(mse / max(mse), 3)) %>%
    dplyr::mutate(site.rank = base::rank(mse.normalized))

  # rank best
  df.rank.best <- df.rank %>%
    dplyr::group_by(phantom, site) %>%
    dplyr::summarize(mean.rank = mean(site.rank)) %>%
    dplyr::arrange(desc(mean.rank))
  df.rank$site <- factor(df.rank$site, levels = as.vector(df.rank.best$site))

  return(df.rank)
}

#' Fit calibration curves
#'
#' Fit pairwise linear regression between calibration and reference sites, using uncalibrated dataset in specified Assay.
#'
#' @param object Calibration Object
#' @param reference.site Character specifying reference site (i.e., reference instrument to calibrated other instruments)
#' @param reference.time Character specifying reference time. Does not need to be specified if only one timePoint is available. One of:
#' \itemize{
#' \item baseline - Default. All calibration curves are generated using reference site measurements at baseline as the reference timpoint.
#' \item match - Calibration curves are constructed using matched time points. E.g., Site measurements are baseline are calibrated using reference site measurements at baseline, 6 month measurements are calibrated using references measurements at 6 months, etc...
#' }
#' @param which.phantom Character indicating which phantom to fit calibration curves for. If unspecified, assumed that only one phantom is provided in data. In this case, if multiple exist, and error will occur.
#' @param sig.intercept.only Logical specifying whether to include intercept in fitted regression equation.
#' \itemize{
#' \item TRUE - Only include intercept in model if significant (p < 0.05).
#' \item FALSE - Default. Include intercept in model regradless of significance.
#' }
#' @param which.parameter Character specifying which parameters to fit calibration curves for. If unspecified, all parameters are fit.
#' @param omit.parameter Character specifying which parameters to omit from calibration curve fitting. Overrides parameters specified by which.parameter argument.
#' @param which.assay Character specifying which assay to fit calibration curves for.
#' @param n.signif Number of significant figures to report.
#' @param verbose Logical specify whether to report progress.
#' @param which.center Specify which measure of central tendency to use when fitting calibration curves. One of:
#' \itemize{
#' \item "mean" - Deafult. mean replicate values.
#' \item "median" - median replicate values. Recommended if >2 replicates scanned, and outliers are present.
#' }
#' @name fitCalibration
#' @seealso \code{\link{calibrateData}}
#' @return Calibration Object
#'
fitCalibration <- function(object, reference.site = NULL, reference.time = "baseline",which.phantom = NULL,
                            sig.intercept.only = F, which.parameter = NULL, omit.parameter = NULL, which.assay = NULL, n.signif = 3, verbose = T, which.center = "mean") {

  # reference.time options:
  #   baseline    baseline time
  #   match       matched to calibration site time
  #   #           custom specified entry (must exist within df)

  # x.fold.val = T,
  # @param n.fold.val Logical indicating whether to perform n-fold calibration. Recommended to reduce influence of outliers.

  #GIGO handling

  # ensure assay exists
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% getAssay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  }

  # ensure assay is specified
  if (is.null(which.assay)) which.assay <- getAssay(object)

  # get data
  df <- object@assays[[which.assay]]@data[["uncalibrated"]]

  # ensure that only one phantom has been specified
  u.phan <- as.character(unique(df$phantom))
  if (is.null(which.phantom)){
    which.phantom <- u.phan
  } else {
    if (!(which.phantom %in% u.phan)) stop(paste(which.phantom, "does not exist"))
  }
  if (length(which.phantom)>1) stop("Multiple phantoms detected. Only one expected.")

  # filter phantoms and parameters
  df <- filterFeatures(df, "phantom", which.phantom)
  df <- filterFeatures(df, "parameter", which.parameter)
  which.parameter <- as.character(unique(df$parameter))
  if (!is.null(omit.parameter)){
    which.parameter <- which.parameter[!(which.parameter %in% omit.parameter)]
    df <- filterFeatures(df, "parameter", which.parameter)
  }

  # ensure sites are properly specified
  if (!("site" %in% colnames(df))) stop ("'site' feature does not exist")
  u.sites <- as.character(unique(df$site))
  if (length(u.sites) < 2) stop ("insufficient number of sites for calibration")
  if (!(reference.site %in% u.sites)) stop (paste(reference.site, " does not exist", sep = ""))
  reference.site.opt <- reference.site

  # ensure times are properly specified
  match.flag <- FALSE
  reference.time.opt <- reference.time

  # check that timepoint feature exists
  if ("timePoint" %in% colnames(df)) {

    # if specified reference.time does not exist, manage contingencies
    if (!(reference.time %in% as.character(unique(df$timePoint)))) {
      if (reference.time == "baseline") {
        reference.time <- min(as.matrix((df %>%
                                           dplyr::filter(site == reference.site) %>%
                                           dplyr::select(timePoint))))
      } else if (reference.time == "match") {
        match.flag <- TRUE
      } else {
        stop ("Specified 'reference.time' does not exist")
    }
  } else {
    df$timePoint <- 0
    reference.time <- 0
  }
  } else {
    stop("'timePoint' feature does not exist")
  }

  # ensure atleast 3 sections are available for cross-calibration
  if (!("section" %in% colnames(df))) stop ("'section' feature does not exist")
  u.sections <- unique(df$section)
  if (length(u.sections) < 3) stop ("Atleast 3 unique imaging phantom sections are required to perform cross-calibration")

  # check parameters
  if (!("parameter" %in% colnames(df))) df$parameter <- "parameter"
  u.parameter <- as.character(unique(df$parameter))

  # ensure data.frame entries are correct class
  try({df$timePoint <- as.numeric(as.character(df$timePoint))}, silent = TRUE)
  if (!is.numeric(df$timePoint)) df$timePoint <- as.character(df$timePoint)
  u.time <- unique(df$timePoint)
  df$value <- as.numeric(as.character(df$value))

  # get all calibration sites
  calibration.sites <- u.sites[u.sites != reference.site]

  # compute calibration curves for each calibration site using specified reference site
  calibrations <- NULL
  calibration.curve.plt <- list()

  for (i in 1:length(u.parameter)){
    for (j in 1:length(u.time)){

      current.parameter <- u.parameter[i]
      current.time <- u.time[j]

      # match if possible
      if (match.flag){
        reference.time <- current.time
        if (length(df %>% dplyr::select(timePoint) %>% dplyr::filter(site == reference.site, time = reference.time)) == 0){
          stop("Cannot match calibration and reference sites at all time points.")
        }
      }

      ref.data <-  df %>%
        dplyr::filter(parameter == current.parameter,
               timePoint == reference.time,
               site == reference.site) %>%
        dplyr::group_by(section) %>%
        dplyr::summarize(mean.val = mean(value, na.rm = T),
                         median.val = median(value, na.rm = T))

      cal.data <- df %>%
        dplyr::filter(parameter == current.parameter,
               timePoint == current.time,
               site %in%  calibration.sites) %>%
        dplyr::group_by(site, section) %>%
        dplyr::summarize(mean.val = mean(value, na.rm = T),
                         median.val = median(value, na.rm = T))


      #####
      repsets <- df %>%
        dplyr::filter(parameter == current.parameter,
                      timePoint == current.time,
                      site %in%  calibration.sites) %>%
        dplyr::group_by(site, parameter, timePoint,phantom) %>%
        dplyr::summarize(repSet = list(unique(scanID)),
                         n.rep = length(unique(scanID)))
      #####

      # generate plots
      match.section.ind <- match(cal.data$section, ref.data$section)

      df.merge2plot <- cal.data
      if (which.center == "mean"){
        colnames(df.merge2plot)[(colnames(df.merge2plot) == "mean.val")] <- "x"
        df.merge2plot$y <- ref.data$mean.val[match.section.ind]
      } else if (which.center == "median"){
        colnames(df.merge2plot)[(colnames(df.merge2plot) == "median.val")] <- "x"
        df.merge2plot$y <- ref.data$median.val[match.section.ind]
      }


      plt.name <- paste("calibrationCurve.", current.parameter, ".t", current.time, sep = "")

      plt.calibration <- df.merge2plot %>%
        ggplot(aes(x, y, colour = site)) +
        geom_point() +
        geom_smooth(method='lm',formula=y~x) + ggtitle("calibration") +
        stat_poly_eq(formula = y~x, aes(label = paste(..rr.label.., sep = "~~~")), rr.digits = 3,label.y = 0.05, label.x = 0.95,parse = TRUE) +
        geom_abline(slope=1, intercept=0, linetype = "dashed") +
        xlab("Calibration Site") +
        ylab("Reference Site") +
        ggtitle(paste("Cross Calibration: ", current.parameter, " (t=", current.time,")", sep = "")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none") +
        facet_wrap(~site)

      calibration.curve.plt[[plt.name]] <- plt.calibration

      for (k in 1:length(calibration.sites)){
        calibration.df <- NULL

        if (which.center == "mean"){
          if (length(dplyr::filter(cal.data, site == calibration.sites[k])$mean.val) == 0){next}
        } else if (which.center == "median"){
          if (length(dplyr::filter(cal.data, site == calibration.sites[k])$median.val) == 0){next}
        }

        ref.cal <- NULL
        try({
          if (which.center == "mean"){
            ref.cal <- data.frame(cal = dplyr::filter(cal.data, site == calibration.sites[k])$mean.val,
                                  ref = ref.data$mean.val)
          } else if (which.center == "median"){
            ref.cal <- data.frame(cal = dplyr::filter(cal.data, site == calibration.sites[k])$median.val,
                                  ref = ref.data$median.val)
          }
        }, silent = T)
        if (is.null(ref.cal)) stop("Reference and calibration pairs could not be matched for ", calibration.sites[k], " and " ,reference.site,
                                   " at t= " , current.time, " for ", current.parameter, " parameter. Consider omitting ", current.parameter,
                                   " prior to fitting calibration curves (see omit.parameter argument).")

        calibration.curve <- lm( ref ~ cal, data = ref.cal)
        calibration.summary <- summary(calibration.curve)


        ##### Xfold calibration
        # if (x.fold.val){
        #
        #
        # cur.repset.df <- as.data.frame(filter(repsets, site ==  calibration.sites[k]))
        # cur.repset <- unlist(cur.repset.df$repSet)
        #
        # n.comb <- combn(cur.repset, 2)
        #
        # calibration.curve.in <- list()
        # calibration.summary.in <- list()
        # effective.int.in <-c()
        # effective.slope.in <- c()
        # effective.r2.in <- c()
        #
        # for (f in 1:ncol(n.comb)){
        #
        #   in.site <- n.comb[ ,f]
        #   out.site <- cur.repset[!(cur.repset %in% in.site)]
        #
        #   in.cal <- df %>%
        #     dplyr::filter(parameter == current.parameter,
        #                   timePoint == current.time,
        #                   site %in% calibration.sites[k],
        #                   scanID %in% in.site) %>%
        #     dplyr::group_by(site, section) %>%
        #     dplyr::summarize(mean.val = mean(value),
        #                      median.val = median(value),
        #                      n.val = length(value))
        #
        #   out.cal <- df %>%
        #     dplyr::filter(parameter == current.parameter,
        #                   timePoint == current.time,
        #                   site %in% calibration.sites[k],
        #                   scanID %in% out.site) %>%
        #     dplyr::group_by(site, section) %>%
        #     dplyr::summarize(mean.val = mean(value),
        #                      median.val = median(value),
        #                      n.val = length(value))
        #
        #   if (which.center == "mean"){
        #     ref.cal.in <- data.frame(cal = dplyr::filter(in.cal, site == calibration.sites[k])$mean.val,
        #                              ref = ref.data$mean.val)
        #   } else if (which.center == "median"){
        #     ref.cal.in <- data.frame(cal = dplyr::filter(in.cal, site == calibration.sites[k])$median.val,
        #                              ref = ref.data$median.val)
        #   }
        #
        #   calibration.curve.in[[f]] <- lm( ref ~ cal, data = ref.cal.in)
        #   calibration.summary.in[[f]] <- summary(calibration.curve)
        #
        #   effective.int.in[f] <- calibration.curve.in[[f]][["coefficients"]][["(Intercept)"]]
        #   effective.slope.in[f] <- calibration.curve.in[[f]][["coefficients"]][["cal"]]
        #   effective.r2.in[f]  <- calibration.summary.in[[f]][["adj.r.squared"]]
        #
        #   out.cal$pred <- calibration.curve.in[[f]][["fitted.values"]]
        #
        # }
        #
        # if (abs(sd(effective.int.in)/mean(effective.int.in)) > 0.5) warning("intercept flag: ", calibration.sites[k])
        # if (abs(sd(effective.slope.in)/mean(effective.slope.in)) > 0.5) warning("slope flag: ", calibration.sites[k])
        #
        # calibration.curve[["coefficients"]][["(Intercept)"]] <- median(effective.int.in)
        # calibration.curve[["coefficients"]][["cal"]] <- median(effective.slope.in)
        # calibration.summary[["adj.r.squared"]] <- median(effective.r2.in)
        # # calibration.summary[["sigma"]]
        # }


        ####

        # intercept handling
        int.p <- calibration.summary[["coefficients"]][1,4]

        if ((sig.intercept.only & int.p <= 0.05) | (!sig.intercept.only)){
          effective.int <- calibration.curve[["coefficients"]][["(Intercept)"]]
          effective.slope <- calibration.curve[["coefficients"]][["cal"]]
          effective.r2  <- calibration.summary[["adj.r.squared"]]
          effective.residual.sem <- calibration.summary[["sigma"]]
          effective.p.int <- calibration.summary[["coefficients"]][1,4]
          effective.p.slope <- calibration.summary[["coefficients"]][2,4]
        } else if (sig.intercept.only & int.p > 0.05) {

          # recompute curve without intercept
          calibration.curve.noint <- lm( ref ~ cal - 1, data = ref.cal)
          calibration.noint.summary <- summary(calibration.curve.noint)

          effective.int <- 0
          effective.slope <- calibration.curve.noint[["coefficients"]][["cal"]]
          effective.r2  <- calibration.summary[["adj.r.squared"]]
          effective.residual.sem <- calibration.summary[["sigma"]]
          effective.p.int <- 1
          effective.p.slope <- calibration.noint.summary[["coefficients"]][1,4]
        } else {
          stop ("troubleshooting needed - unaccounted condition encountered")
        }

        # store results

        calibration.df <- data.frame(parameter = current.parameter,
                                     phantom = which.phantom,
                                     reference.site = reference.site,
                                     calibration.site = calibration.sites[k],
                                     reference.time = reference.time,
                                     calibration.time = current.time,
                                     intercept = effective.int,
                                     slope = effective.slope,
                                     r2 = effective.r2,
                                     residual.sem = effective.residual.sem,
                                     p.intercept = effective.p.int,
                                     p.slope = effective.p.slope)


        suppressWarnings({calibrations <- bind_rows(calibrations, calibration.df)})
      }

    }
  }

  # assign names to calibration variables
  variables2round <- c("intercept", "slope", "r2", "residual.sem", "p.intercept", "p.slope")
  calibrations[,variables2round] <- lapply(calibrations[, variables2round], signif, n.signif)

  # overwrite pre-existing calibration fits
  existing.calibrations <- names(object@assays[[which.assay]]@calibration)
  if (verbose){
    if (length(existing.calibrations) > 0) {
      warning(paste("Pre-existing'", which.assay, "' calibration curves were overwritten", sep = ""))
    } else {
      cat("\n")
      cat(paste("fitCalibration results created for ", which.phantom, sep = ""))
      cat("\n")
      }
  }

  # store calibration results
  object@assays[[which.assay]]@calibration <- list(calibration.equations = calibrations,
                                                   reference.site = reference.site.opt,
                                                   reference.time = reference.time.opt,
                                                   calibration.curves = calibration.curve.plt,
                                                   phantom = which.phantom,
                                                   parameter = which.parameter)

  return(object)
}


#' Calibrate Data
#'
#' Calibrate data set using fit calibration curves. Calibration curves must exist, having been generated by fit.calibration function.
#'
#' New column is added to calibrated dataset, name cal.phantom, specifying which phantom was used for calibraiton. That is, if the uncalibrated data contains datasets from multiple phantom-types, the 'phantom' column specifies which phantom the values originate from and the 'cal.phantom' column specifies which phantom was used for calibration.
#'
#' @param object Calibration Object
#' @param which.assay Character specifying which assay to calibate.
#' @param verbose Logical specifying whether progress is reported.
#' @name calibrateData
#' @return Calibration Object
#' @seealso \code{\link{fitCalibration}}
calibrateData <- function(object, which.assay = NULL, verbose = T) {

  #GIGO handling

  # ensure assay exists
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% getAssay(object, which.assay = "all"))){
      stop (paste(which.assay, " does not exist"), sep = "")
    }
  }

  # ensure assay is specified
  if (is.null(which.assay)) which.assay <- getAssay(object)

  # get data
  df <- object@assays[[which.assay]]@data[["uncalibrated"]]

  # get calibration
  calibration <- object@assays[[which.assay]]@calibration[["calibration.equations"]]
  cal.sub <- calibration %>% dplyr::select(parameter, phantom, calibration.site, calibration.time, intercept, slope, p.intercept)
  colnames(cal.sub) <- c("parameter", "cal.phantom", "site", "timePoint", "intercept", "slope", "p.intercept")

  # question - MERGE BY PHANTOM? for now, no.

  # join dataframes
  df.merge <- suppressMessages({join(df, cal.sub)})
  df.merge$slope[is.na(df.merge$slope)] <- 1
  df.merge$intercept[is.na(df.merge$intercept)] <- 0
  df.merge$p.intercept[is.na(df.merge$p.intercept)] <- 1

  df.merge$value.cal <- (df.merge$value * df.merge$slope) + df.merge$intercept

  # remove uncalibrated values
  df.merge <- dplyr::select(df.merge, -c("value"))

  # rename calibrated values
  colnames(df.merge)[colnames(df.merge) %in% "value.cal"] <- "value"

  df.final <- df.merge %>% dplyr::select(c(colnames(df), "cal.phantom"))

  existing.data <- object@assays[[which.assay]]@data
  if (verbose){
    if ("calibrated" %in% names(existing.data)){
      warning(paste("Pre-existing '", which.assay, "' calibrated data was overwritten", sep = ""))
    } else {
      cat("\n")
      cat(paste("data succesfully calibrated", sep = ""))
      cat("\n")}
  }

  object@assays[[which.assay]]@data[["calibrated"]] <- df.final

  return(object)
}


