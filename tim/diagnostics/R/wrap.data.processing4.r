#======================================================================
#
# DART software - Copyright © 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$
#
# ncdf.functions.r:get.var.ncdf.whole() maintains singleton array dimensions
# and returns a list with the storage order. The required storage order is 
#
# in netcdf: (analysisT, stations, levels, copy, nmembers, forecast_lead)
# in      R:  forecast_lead, nmembers, copy, levels, stations, analysisT 
#
# Most of the context and functions came from
# /ptmp/syha/R/work/JME_conus/Section_5/
#======================================================================

# data and plotting script ...
/glade/user/syha/R/work/JME_conus/Section_4.3/V.T.P4.Q.0.75boot.100.RData
/glade/user/syha/R/work/JME_conus/module_rel_bw.48h.R

rm(list=ls())

library(ncdf)
library(verification)

source("DART.functions.r")
source("ens.fcast.bootstrap.r")

out.dir <- "/Users/thoar/R"

quantiles     <- c(0, 0.25, 0.5, 0.75, 1.0)   ## set quantile array
quantiles     <- c(0.5)   ## set quantile array

add.obs.error <- FALSE   ## add observational error to each forecast
bias.dims     <- c("ensemble","level","station")
group.dims    <- c("level","station")
n.boot        <- 100     # number of bootstrap samples
n.boot        <- 1       # TJH just to test

prob.method   <- "quantile" ## or "fixed"
fixed.thres   <-  5         ## m/s fixed threshold

## by default, bs is alway calculates because it is pretty fast.

run.bs        <- TRUE
run.crps      <- TRUE
run.rank      <- FALSE
run.attr      <- TRUE

# files <- system('ls omb*.nc', intern = TRUE)
fname <- 'Original/omb_metar_BS_d01.nc'

##======================================================================
## Written to calculate one level at a time.  
## P = 1 for surface if using ncdf 1.7

varname   <- 'METAR_U_10_METER_WIND'
varname   <- 'T'

## these only used to constuct out.file ...

P         <- 1  
Q         <- 0.5
var.char1 <- "T"       ## U, V, T or WS

out.file <- sprintf('%s/V_%s_level_%d_quantile_%.2f_boot_%d.RData',out.dir,var.char1,P,Q,n.boot)

#source("/glade/home/syha/R/scripts/ens.dataprocessing4.R")

bob = ens.fcast.bootstrap(filename=fname, varname=varname, quantiles=quantiles,
            add.obs.error=add.obs.error, bias.dims=bias.dims, group.dims=group.dims,
            n.boot=n.boot, prob.method=prob.method,fixed.thres=fixed.thres,
            run.bs=run.bs, run.crps=run.crps, run.rank=run.rank, run.attr=run.attr, 
            out.file=out.file)

##======================================================================

var.char1 <- "T"
out.file  <- paste(out.dir,"V.", var.char1, ".P", P,".Q.", Q, "boot.", n.boot,  ".RData", sep = "") 

source("/glade/home/syha/R/scripts/ens.dataprocessing4.R")

##======================================================================

n.boot    <- 1
run.attr  <- TRUE
run.rank  <- TRUE
var.char1 <- "WS"
out.file  <- paste(out.dir,"V.", var.char1, ".P", P,".Q.", Q, "boot.", n.boot,  ".RData", sep = "") 

source("/glade/home/syha/R/scripts/ens.dataprocessing4.R")

##======================================================================

var.char1 <- "WS"
out.file  <- paste(out.dir,"V.", var.char1, ".P", P,".Q.", Q, "boot.", n.boot,  ".RData", sep = "") 

source("/glade/home/syha/R/scripts/ens.dataprocessing4.R")
