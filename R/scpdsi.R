# Computation of the Self-calibrating Palmer Drought Severity Index (scPDSI).

#' @useDynLib scPDSIr
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  ops <- options()
  pdsi.ops <- list(
    # Calculating the conventional PDSI
    # Duration factors
    PDSI.p = 0.897,
    PDSI.q = 1/3,

    # Coefficients of climate characteristic
    PDSI.coe.K1.1 = 1.5,
    PDSI.coe.K1.2 = 2.8,
    PDSI.coe.K1.3 = 0.5,

    PDSI.coe.K2 = 17.67
  )

  toset <- !(names(pdsi.ops) %in% names(ops))
  if(any(toset)) {
    options(pdsi.ops[toset])
  }
}

#' Calculate the (sc)PDSI
#' @description Calculating the monthly Self-calibrating PDSI (scPDSI) using the precipitation
#'              and potential evapotranspiration.
#'
#'
#' @details
#'
#' The Palmer Drought Severity Index (PDSI), proposed by Palmer (1965), is a
#' widely used drought indicator to quantify the long-term drought conditions, for
#' an area at a certain time. The PDSI is a semi-physical based drought index calculated
#' using the precipitation and potential evapotranspiration data, based on a simple
#' two-layer bucket water balance model. Conventionally, the constants to calculate
#' the PDSI were firstly empirically derived by using the meteorological records in
#' Kansas and Iowa in middle US with a semi-arid climate conditions, therefore the
#' conventional PDSI usually could not satisfactorily represent the drought conditions
#' for other areas around the world, which also makes spatial comparisons of PDSI values
#' difficult.
#'
#' For this, Wells et al. (2004) proposed a self-calibrating Palmer Drought Severity
#' Index (scPDSI). The scPDSI could automatically adjust the empirical constants
#' in the PDSI computation with dynamically calculated values. Several works have
#' proved that the scPDSI performs better in spatially comparison than the conventional
#' PDSI. For more details please see the works of Wells et al. (2004).
#'
#'
#' @return
#' This function return a vector of scpdsi.
#'
#'
#' @references Palmer W., 1965. Meteorological drought. U.s.department of Commerce
#'             Weather Bureau Research Paper.
#'
#'             Wells, N., Goddard, S., Hayes, M. J., 2004. A Self-Calibrating Palmer
#'             Drought Severity Index. Journal of Climate, 17(12):2335-2351.
#'
#' @examples
#' 
#' library(scPDSIr)
#' scpdsi(indat,  start, end, cal_start, cal_end)
#' Where the dataframe "indat" contains the following variables:
#' 
#' P = Monthly precipitation series without NA [mm]. Can be a time series.
#'
#' PE = Monthly potential evapotranspiration corresponding to the precipitation
#'            series. Can be calculated by the Penman-Monteith or the Thonthwate
#'            equation [mm].
#' lat = Latitude  (degrees)
#' 
#' lon = Longitude (degrees)
#' 
#' The following are provided as separate arguments:
#' 
#' AWC Available Water Capacity of the soil layer [mm]. This will be retrieved automatically
#' from a 0.1 degree European grid and does not need to be supplied. However 
#' if using a single site this can be specified if required.
#'
#' start Integer. Start year of the PDSI to be calculated default 1.
#'
#' end Integer. End year of the PDSI to be calculated.
#'
#' cal_start Integer. Start year of the calibrate period. Default is start year.
#'
#' cal_end Integer. End year of the calibrate period. Default is end year.
#' 
#' @importFrom stats ts
#'
#' @export
pdsi <- function(indat, AWC = NULL, start = NULL, end = NULL, cal_start = NULL, cal_end = NULL,
                 sc = TRUE) {

  P <- indat$prec
  PE <- indat$pet
 
  ndat <- length(P)
  if(sum(is.na(P)) > 0 | sum(is.na(PE)) > 0)
  {
    psdi <- array(dim=ndat, NA)
    return(psdi)
  }
  
  # check P for low values
  t.str <- strptime(indat$prec.dates, "%Y-%m-%d")
  yr <- as.numeric(format(t.str, "%Y"))
  mn <- as.numeric(format(t.str, "%m"))
  dat <- data.frame("yr", "mn", "P")
  P_yr_sum <- aggregate(P~yr, dat , sum)
  P_yr_sum_mean <- mean(P_yr_sum$P, na.rm=TRUE)
  if(P_yr_sum_mean < 50)
  {
    psdi <- array(dim=ndat, NA)
    return(psdi)
  }
  lat <- indat$lat
  lon <- indat$lon
  freq <- 12

  if(is.null(start)) start <-  1;
  if(is.null(end)) end <- start + ceiling(length(P)/freq) - 1

  if(is.null(cal_start)) cal_start <- start
  if(is.null(cal_end)) cal_end <- end

  data("awc_dat", package="scPDSIr")
  if(is.null(AWC))
  {
    AWC <- get_awc(lat, lon)
    if(is.na(AWC)){
      AWC <- 100
    }
  }

  res <- C_pdsi(P, PE, AWC, start, end, cal_start, cal_end, sc,
                getOption("PDSI.coe.K1.1"),
                getOption("PDSI.coe.K1.2"),
                getOption("PDSI.coe.K1.3"),
                getOption("PDSI.coe.K2"),
                getOption("PDSI.p"),
                getOption("PDSI.q"))

  #names(res) <- c("inter.vars", "clim.coes", "calib.coes")

  inter.vars <- res[[1]]
  clim.coes <- res[[2]]
  calib.coes <- res[[3]]
  inter.vars[inter.vars == -999.] <- NA
# browser()
  
  # out <- list(call = match.call(expand.dots=FALSE),
  #             X = ts(inter.vars[, 14], start = start, frequency = freq),
  #             PHDI = ts(inter.vars[, 15], start = start, frequency = freq),
  #             WPLM = ts(inter.vars[, 16], start = start, frequency = freq),
  #             inter.vars = ts(inter.vars[, 3:13], start = start, frequency = freq))
  # 
  # colnames(out$inter.vars) <- c("P", "PE", "PR", "PRO", "PL", "d", "Z",
  #                               "Prob", "X1", "X2", "X3")
  # 
  # dim(calib.coes) <- c(2, 5)
  # colnames(calib.coes) <- c("m", "b", "p", "q", "K2")
  # rownames(calib.coes) <- c('wet', 'dry')
  # 
  # colnames(clim.coes) <- c("alpha", "beta", "gamma", "delta", "K1")
  # rownames(clim.coes) <- month.name
  # 
  # out$clim.coes <- clim.coes
  # out$calib.coes <- calib.coes
  # 
  # out$self.calib <- sc
  # out$range <- c(start, end)
  # out$range.ref <- c(cal_start, cal_end)
  # 
  # class(out) <- "pdsi"
  # out
  
  scpdsi <- inter.vars[, 14]
  # big_vals_index <- which(scpdsi > 10 | scpdsi < -10 | is.infinite(scpdsi))
  # if(sum(big_vals_index) > 0)
  # {
  #   filename <- sprintf("big_vals_log_thresh_0/bigvals_%f_%f.csv", lat,lon)
  #   output <- t(cbind(big_vals_index, scpdsi[big_vals_index]))
  #   write(output, file=filename, ncolumns=2, sep=',')
  #   # write coord logfile
  #   coord_out <- t(cbind(lat,lon))
  #   write(coord_out, "coord_log_thresh_0.csv",  ncolumns=2, sep=',', append=TRUE)
  #   
  # }
  scpdsi[is.infinite(scpdsi) & scpdsi > 0] <- 10
  scpdsi[is.infinite(scpdsi) & scpdsi < 0] <- -10
  scpdsi[scpdsi > 10] <- 10
  scpdsi[scpdsi < -10] <- -10
  
  return(scpdsi)
}
