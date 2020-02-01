


#' Identify Reference Site
#'
#' Sites are ranked according to mean-squared error calculated between site-specific parameters and overall median parameter values (pooled across all timePoints and sites), and top ranking site is returned.
#'
#' @param object Calibration Object
#' @param which.assay Character specifying which assay to identify reference site for.
#' @param which.data Character specifying which data to identify reference site for.
#' @name identify.reference
#' @seealso \code{\link{consistency.analysis}}, \code{\link{consistency.plot}}
#' @return Character
#'
identify.reference <- function(object,  which.assay = NULL, which.data = "uncalibrated") {

  # ensure assay is specified
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% get.assay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  } else {
    which.assay <- get.assay(object)
  }

  # ensure data is specified
  if (!is.null(which.data)) {
    stopifnot(class(which.data) == "character")
    if (!(which.data %in% get.datasets(object, which.assay = which.assay))){
      stop ("'which.data' does not exist")
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
    dplyr::mutate(site.rank = rank(mse.normalized))

  # rank best
  df.rank.best <- df.rank %>%
    dplyr::group_by(phantom, site) %>%
    dplyr::summarize(mean.rank = mean(site.rank)) %>%
    dplyr::arrange(desc(mean.rank))

  reference.site <- as.vector(df.rank.best$site[which(min(df.rank.best$mean.rank) == df.rank.best$mean.rank)])

  return(reference.site)
}


#' Evaluate site consistency
#'
#' Sites are ranked according to mean-squared error calculated between site-specific parameters and overall median parameter values (pooled across all timePoints and sites).
#'
#' @param object Calibration Object
#' @param which.assay Character specifying which assay to identify reference site for.
#' @param which.data Character specifying which data to identify reference site for.
#' @name consistency.analysis
#' @seealso \code{\link{identify.reference}}, \code{\link{consistency.plot}}
#' @return data.frame
#'
consistency.analysis <- function(object,  which.assay = NULL, which.data = "uncalibrated") {

  # ensure assay is specified
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% get.assay(object, which.assay = "all"))){
      stop ("'which.assay' does not exist")
    }
  } else {
    which.assay <- get.assay(object)
  }

  # ensure data is specified
  if (!is.null(which.data)) {
    stopifnot(class(which.data) == "character")
    if (!(which.data %in% get.datasets(object, which.assay = which.assay))){
      stop ("'which.data' does not exist")
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
#' @param sig.intercept.only Logical specifying whether to include intercept in fitted regression equation.
#' \itemize{
#' \item TRUE - Only include intercept in model if significant (p < 0.05).
#' \item FALSE - Default. Include intercept in model regradless of significance.
#' }
#' @param which.assay Character specifying which assay to fit calibration curves for.
#' @param n.signif Number of significant figures to report.
#' @param verbose Logical specify whether to report progress.
#' @name fit.calibration
#' @return Calibration Object
#'
fit.calibration <- function(object, reference.site = NULL, reference.time = "baseline",
                            sig.intercept.only = F, which.assay = NULL, n.signif = 3, verbose = T) {

  # reference.time options:
  #   baseline    baseline time
  #   match       matched to calibration site time
  #   #           custom specified entry (must exist within df)

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

  # get data
  df <- object@assays[[which.assay]]@data[["uncalibrated"]]

  # ensure sites are properly specified
  if (!("site" %in% colnames(df))) stop ("'site' feature does not exist")
  u.sites <- as.character(unique(df$site))
  if (length(u.sites) < 2) stop ("insufficient number of sites for calibration")
  if (!is.null(reference.site)){
    if (!(reference.site %in% u.sites)) stop ("specified 'reference.site' does not exist")
  }
  if (is.null(reference.site)) reference.site <- identify.referenceSite(object, which.assay)
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

      # cur.calibration <- calibrate.values(df, u.par[i], u.time[j], reference.site)
      ref.data <-  df %>%
        dplyr::filter(parameter == current.parameter,
               timePoint == reference.time,
               site == reference.site) %>%
        dplyr::group_by(section) %>%
        dplyr::summarize(mean.val = mean(value))

      cal.data <- df %>%
        dplyr::filter(parameter == current.parameter,
               timePoint == current.time,
               site %in%  calibration.sites) %>%
        dplyr::group_by(site, section) %>%
        dplyr::summarize(mean.val = mean(value))


      # generate plots
      match.section.ind <- match(cal.data$section, ref.data$section)

      df.merge2plot <- cal.data
      colnames(df.merge2plot)[(colnames(df.merge2plot) == "mean.val")] <- "x"
      df.merge2plot$y <- ref.data$mean.val[match.section.ind]

      plt.name <- paste("calibration.curve.", current.parameter, ".t", current.time, sep = "")

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

      # if (plot.flag) print(plt.calibration)

      for (k in 1:length(calibration.sites)){
        calibration.df <- NULL

        if (length(dplyr::filter(cal.data, site == calibration.sites[k])$mean.val) == 0){next}


        ref.cal <- data.frame(cal = dplyr::filter(cal.data, site == calibration.sites[k])$mean.val,
                              ref = ref.data$mean.val)

        calibration.curve <- lm( ref ~ cal, data = ref.cal)
        calibration.summary <- summary(calibration.curve)

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
      cat(paste("fit.calibration results created", sep = ""))
      cat("\n")
      }
  }

  # store calibration results
  object@assays[[which.assay]]@calibration <- list(fit.calibration = calibrations,
                                                   reference.site = reference.site.opt,
                                                   reference.time = reference.time.opt,
                                                   calibration.curve.plots = calibration.curve.plt)


  # save plot handle
  plt.name <- paste("calibration.curves.", reference.site.opt, ".", reference.time.opt, sep = "")
  existing.plots <- object@assays[[which.assay]]@plots
  existing.plot.names <- names(existing.plots)

  if (verbose){
    if (plt.name %in% existing.plot.names){
      warning(paste("Pre-existing'", plt.name, "' in '", which.assay , "' was overwritten", sep = ""))
    } else {
      cat("\n")
      cat(paste("calibration curves created and saved as '", plt.name, "'", sep = ""))
      cat("\n")
    }
  }

  existing.plots[[plt.name]] <- calibration.curve.plt
  object@assays[[which.assay]]@plots <- existing.plots


  return(object)
}


#' Calibrate Data
#'
#' Calibrate data set using fit calibration curves. Calibration curves must exist, having been generated by fit.calibration function.
#'
#' @param object Calibration Object
#' @param which.assay Character specifying which assay to calibate.
#' @param verbose Logical specifying whether progress is reported.
#' @name calibrate.data
#' @return Calibration Object
#' @seealso \code{\link{fit.calibration}}
calibrate.data <- function(object, which.assay = NULL, verbose = T) {

  #GIGO handling

  # ensure assay exists
  if (!is.null(which.assay)) {
    stopifnot(class(which.assay) == "character")
    if (!(which.assay %in% get.assay(object, which.assay = "all"))){
      stop ("Specified 'which.assay' does not exist")
    }
  }

  # ensure assay is specified
  if (is.null(which.assay)) which.assay <- get.assay(object)

  # get data
  df <- object@assays[[which.assay]]@data[["uncalibrated"]]

  # get calibration
  calibration <- object@assays[[which.assay]]@calibration[["fit.calibration"]]
  cal.sub <- calibration %>% select(parameter, calibration.site, calibration.time, intercept, slope, p.intercept)
  colnames(cal.sub) <- c("parameter", "site", "timePoint", "intercept", "slope", "p.intercept")

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

  df.final <- df.merge %>% dplyr::select(colnames(df))

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


