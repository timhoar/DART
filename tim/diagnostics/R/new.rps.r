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
#======================================================================
#
# obs   <- rnorm(10)
# pred  <- matrix(rnorm(100), ncol = 10 )
# q     <- quantile(obs, c(1/3, 2/3) )
# INF   <- 1000
# thres <- c(-INF, q, INF)
#
# rps2(obs, pred, thres)

rps2 <- function(obs, pred, thres){
   ### assume obs and pred are in form of vector and matrix, original values
   ne <- ncol(pred)
   mn <- min(c(obs,unlist(pred)), na.rm = TRUE)
   mx <- max(c(obs,unlist(pred)), na.rm = TRUE)

   if(mn < min(thres)|mx > max(thres)){stop("thresholds too narrow")}

   OBSbin  <- as.numeric(cut(obs, thres))
   PREDbin <- as.numeric(cut(as.matrix(pred), thres))
   PREDbin <- matrix(PREDbin, ncol = ne)

   #### turn into cumulative matrix
   nf   <- length(thres)-1
   PRED <- matrix(NA, nrow = nrow(pred), ncol = nf)

   for(i in 1:nf){
      f <- function(x, i){sum(x <= i)}
      PRED[,i] <- apply(PREDbin, 1, f, i)/ ne
   }

   O <- matrix(0, nrow = nrow(pred), ncol = nf)

   for(i in 1:nrow(O)){
      O[i,OBSbin[i]:nf] <- 1
   }

   ##########
   rps <- mean((PRED - O)^2)

   ### rps climo
   x         <- table(OBSbin)/length(OBSbin)
   pred.clim <- matrix(rep(x, length(OBSbin)), byrow = TRUE, ncol = nf)

   PRED.clim <- matrix(NA, nrow = nrow(pred), ncol = nf )

   for(i in 1:nf){
      PRED.clim[,i] <- apply(as.matrix(pred.clim[,1:i] ),  1, sum)
   }

   rps.clim <- mean((PRED.clim - O)^2)
   rps.ss   <- 1 - rps/rps.clim

   return(list(rps = rps, rps.clim = rps.clim, rps.ss = rps.ss))
}

