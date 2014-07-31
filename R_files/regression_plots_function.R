regression.plots <- function (ests,varname,regionname) {
  lbub <- apply(ests, 2, quantile,(probs = c(lb, ub)))
  b.range <- range(lbub[1, ], lbub[2, ]) + range(lbub[1, ], lbub[2, ])*.01
  
  plot(congnum, apply(ests, 2, median),
       ylim = b.range,
       pch = 20,
       cex.axis = .8,
       xlab = "Congress",
       ylab = "Estimates",
       main = list(varname, cex = .9),
       sub = regionname)
  
  segments(congnum, lbub[1, ], congnum, lbub[2, ])
  lines (c(72, 81), c(0, 0), lty = 3, lwd = 1.5, col = "turquoise4")
  cat(regionname, varname, "\n")
  print(cbind(congnum, apply(ests, 2, median), lbub[1, ], lbub[2, ]))
}