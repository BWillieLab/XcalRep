

#' Fit calibration curves
#'
#' Fit pairwise and
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
fit.calibration <- function(object, reference.site = NULL, reference.time = "baseline", which.assay = NULL, n.signif = 3, show.data.table = TRUE, plot.flag = TRUE) {


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
    if (!(reference.site %in% u.sites)) stop ("specfied 'reference.site' does not exist")
  }
  if (is.null(reference.site)) reference.site <- identify.referenceSite(object, which.assay)
  reference.site.opt <- reference.site

  # ensure times are properly specified
  match.flag <- FALSE
  reference.time.opt <- reference.time
  if ("timePoint" %in% colnames(df)) {
    # u.time <- as.character(unique(df$timePoint))
    if (!(reference.time %in% as.character(unique(df$timePoint)))) {
      if (reference.time == "baseline") {
        reference.time <- min(as.matrix((df %>%
                                           filter(site == reference.site) %>%
                                           select(timePoint))))
      } else if (reference.time == "match") match.flag <- TRUE
      else stop ("Specified 'reference.time' does not exist")
    }
  } else {
    df$timePoint <- 0
    reference.time <- 0
  }

  # ensure atleast 3 sections are available for cross-calibration
  if (!("section" %in% colnames(df))) stop ("'section' feature does not exist")
  u.sections <- unique(df$section)
  if (length(u.sections) < 3) stop ("Atleast 3 unique sections are required to perform cross-calibration")

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
        if (length(df %>% select(timePoint) %>% filter(site == reference.site, time = reference.time)) == 0){
          stop("Cannot match calibration and reference sites at all time points.")
        }
      }

      # cur.calibration <- calibrate.values(df, u.par[i], u.time[j], reference.site)
      ref.data <-  df %>%
        filter(parameter == current.parameter,
               timePoint == reference.time,
               site == reference.site) %>%
        group_by(section) %>%
        summarize(mean.val = mean(value))

      cal.data <- df %>%
        filter(parameter == current.parameter,
               timePoint == current.time,
               site %in%  calibration.sites) %>%
        group_by(site, section) %>%
        summarize(mean.val = mean(value))


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
        theme(legend.position = "none") +
        facet_wrap(~site)

      calibration.curve.plt[[plt.name]] <- plt.calibration

      if (plot.flag) print(plt.calibration)

      for (k in 1:length(calibration.sites)){
        calibration.df <- NULL

        if (length(filter(cal.data, site == calibration.sites[k])$mean.val) == 0){next}


        ref.cal <- data.frame(cal = filter(cal.data, site == calibration.sites[k])$mean.val,
                              ref = ref.data$mean.val)

        calibration.curve <- lm( ref ~ cal, data = ref.cal)
        calibration.summary <- summary(calibration.curve)

        # store results
        calibration.df <- data.frame(parameter = current.parameter,
                                     reference.site = reference.site,
                                     calibration.site = calibration.sites[k],
                                     reference.time = reference.time,
                                     calibration.time = current.time,
                                     intercept = calibration.curve[["coefficients"]][["(Intercept)"]],
                                     slope = calibration.curve[["coefficients"]][["cal"]],
                                     r2 = calibration.summary[["adj.r.squared"]],
                                     residual.sem = calibration.summary[["sigma"]],
                                     p.intercept = calibration.summary[["coefficients"]][1,4],
                                     p.slope = calibration.summary[["coefficients"]][2,4])


        suppressWarnings({calibrations <- bind_rows(calibrations, calibration.df)})
      }

    }
  }

  variables2round <- c("intercept", "slope", "r2", "residual.sem", "p.intercept", "p.slope")
  calibrations[,variables2round] <- lapply(calibrations[, variables2round], signif, n.signif)

  # print results in data table
  if (show.data.table) print(datatable(calibrations, filter = 'top'))

  # overwrite pre-existing calibration fits
  existing.calibrations <- names(object@assays[[which.assay]]@calibration)
  if (length(existing.calibrations) > 0) {
    warning(paste("Pre-existing'", which.assay, "' calibration curves were overwritten", sep = ""))
  } else {
    cat("\n============================\n")
    cat(paste("calibration.fit created", sep = ""))
    cat("\n")}
  object@assays[[which.assay]]@calibration <- list(calibration.fit = calibrations,
                                                   reference.site = reference.site.opt,
                                                   reference.time = reference.time.opt,
                                                   calibration.curve.plots = calibration.curve.plt)


  # save plot handle
  plt.name <- paste("calibration.curves.", reference.site.opt, "RefSite.", reference.time.opt, "RefTime", sep = "")
  existing.plots <- object@assays[[which.assay]]@plots
  existing.plot.names <- names(existing.plots)

  if (plt.name %in% existing.plot.names){
    warning(paste("Pre-existing'", plt.name, "' in '", which.assay , "' was overwritten", sep = ""))
  } else {
    cat("\n============================\n")
    cat(paste("calibration curves created and saved as '", plt.name, "'", sep = ""))
    cat("\n")}

  existing.plots[[plt.name]] <- calibration.curve.plt
  object@assays[[which.assay]]@plots <- existing.plots


  return(object)
}





calibrate.data <- function(object, which.assay = NULL, sig.intercept.only = TRUE) {

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

  # stopifnot(class(new.data.name) == "character")


  # get data
  df <- object@assays[[which.assay]]@data[["uncalibrated"]]

  # get calibration
  calibration <- object@assays[[which.assay]]@calibration[["calibration.fit"]]
  cal.sub <- calibration %>% select(parameter, calibration.site, calibration.time, intercept, slope, p.intercept)
  colnames(cal.sub) <- c("parameter", "site", "timePoint", "intercept", "slope", "p.intercept")

  # join dataframes
  df.merge <- suppressMessages({join(df, cal.sub)})
  df.merge$slope[is.na(df.merge$slope)] <- 1
  df.merge$intercept[is.na(df.merge$intercept)] <- 0
  df.merge$p.intercept[is.na(df.merge$p.intercept)] <- 1

  # specify effective intercept term
  if (sig.intercept.only){
    df.merge$effective.intercept <- df.merge$intercept
    df.merge$effective.intercept[df.merge$p.intercept >0.5] <- 0
  } else {df.merge$effective.intercept}

  # calibrate data

  df.merge$value.cal <- (df.merge$value * df.merge$slope) + df.merge$effective.intercept


  df.final <- df.merge %>% select(value.cal, colnames(df)[seq(2, ncol(df))])

  which.value.cal <- which(colnames(df.final) == "value.cal")
  colnames(df.final)[which.value.cal] <- "value"

  existing.data <- object@assays[[which.assay]]@data
  if ("calibrated" %in% names(existing.data)){
    warning(paste("Pre-existing '", which.assay, "' calibrated data was overwritten", sep = ""))
  } else {
    cat("\n============================\n")
    cat(paste("data succesfully calibrated", sep = ""))
    cat("\n")}

  object@assays[[which.assay]]@data[["calibrated"]] <- df.final

  return(object)
}


