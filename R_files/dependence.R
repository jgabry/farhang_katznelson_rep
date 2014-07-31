library(coda)
library(MCMCpack)
library(energy)
library(rstan)
library(parallel)

setGeneric("dependence", function(object, ...) standardGeneric("dependence"))
setMethod("dependence", "array",
          function(object, lag = 1, biased = FALSE) {
            pars <- dim(object)[3]
            if(is.na(pars)) stop("'object' must be an array with 3 dimensions")
            chains <- ncol(object)
            sims <- nrow(object)
            arr <- array(apply(object, MARGIN = 2, function(x) {
              cbind(x[-(sims:(sims - lag + 1)),], x[-(1:lag),])
            }), dim = c(sims - lag, 2 * pars, chains))
            FUN <- if(biased) dcor else bcdcor
            apply(arr, MARGIN = 1, function(y) {
              FUN(t(y[1:pars,]), t(y[-c(1:pars),]))
            })
          }
)

setMethod("dependence", "stanfit",
          function(object, lag = 1:7, biased = FALSE,
                   plot = TRUE, pars = object@sim$pars_oi, f = 1/4, ...) {
            stopifnot(object@mode == 0)
            x <- extract(object, pars, permuted = FALSE, inc_warmup = TRUE)
            dcors <- sapply(lag, FUN = function(l) dependence(x, lag = l))
            warmup <- object@stan_args[[1]]$warmup / object@stan_args[[1]]$thin    
            if(plot && length(lag) == 1) {
              lag <- lag * object@stan_args[[1]]$thin
              plot(dcors, pch = ".", ylim = 0:1, las = 1,
                   xlab = "Iteration",
                   ylab = paste("Estimated distance correlation among", 
                                ncol(object), "chains"),           
                   main = paste(object@model_name, "model, lag =", lag), ...)
              lines(lowess(dcors, f = f), type = "l", col = 2, lty = 2)
              abline(v = warmup, col = 3, lty = 3)
              legend("topright", legend = c("Value", "Smoothed"),
                     col = 1:2, lty = c(NA_integer_, 2), 
                     pch = c(20L, NA_integer_), cex = 0.75)
              text(x = warmup, y = 1.02, labels = "Warmup  | Retained")
            }
            else if(plot) {
              plot(dcors[[1]], type = "n", ylim = 0:1, las = 1,
                   xlab = "Iteration", 
                   ylab = paste("Estimated distance correlation among", 
                                ncol(object), "chains"),
                   main = paste(object@model_name, "model"), ...)
              sapply(1:length(dcors), FUN = function(i) {
                lines(lowess(dcors[[i]], f = f), col = i + 1)
              })
              abline(v = warmup, col = 1, lty = 3)
              legend("topright", legend = lag * object@stan_args[[1]]$thin, 
                     title = "lag length", ncol = 2,
                     col = 1 + lag, lty = 1, cex = 0.75)
              text(x = warmup, y = 1.02, labels = "Warmup  | Retained")
            }
            return(dcors)
          }
)

setOldClass("mcmc.list")
setMethod("dependence", "mcmc.list",
          function(object, lag = 1:7, biased = FALSE, 
                   plot = TRUE, pars = colnames(object[[1]]), f = 1/4, ...) {
            x <- array(NA_real_, 
                       dim = c(nrow(object[[1]]), length(object), ncol(object[[1]])))
            dimnames(x)[[3]] <- pars
            for(i in seq_along(object)) x[,i,] <- as.matrix(object[[i]])[,pars]
            dcors <- sapply(lag, FUN = function(l) dependence(x, lag = l))
            burnin <- attributes(object[[1]])$mcpar[1] - 1L
            chains <- length(object)
            title <- attributes(object[[1]])$title
            thin <- attributes(object[[1]])$mcpar[3]
            if(plot && length(lag) == 1) {
              plot(dcors, pch = ".", ylim = 0:1, las = 1,
                   xlab = "Iteration",
                   ylab = paste("Estimated distance correlation among", 
                                chains, "chains"),           
                   main = paste(title, "lag =", lag * thin), ...)
              lines(lowess(dcors, f = f), type = "l", col = 2, lty = 2)
              legend("topright", legend = c("Value", "Smoothed"),
                     col = 1:2, lty = c(NA_integer_, 2), 
                     pch = c(20L, NA_integer_), cex = 0.75)
            }
            else if(plot) {
              plot(dcors[[1]], type = "n", ylim = 0:1, las = 1,
                   xlab = "Iteration", 
                   ylab = paste("Estimated distance correlation among", 
                                chains, "chains"),
                   main = paste(title), ...)
              sapply(1:length(dcors), FUN = function(i) {
                lines(lowess(dcors[[i]], f = f), col = i + 1)
              })
              legend("topright", legend = lag * thin, title = "lag length", ncol = 2,
                     col = 1 + lag, lty = 1, cex = 0.75)
            }
            return(dcors)
          }
)
