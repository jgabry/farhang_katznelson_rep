# make_gg_lines function
make_gg_lines <- function(Bnames, point.est = "median", conf.level = 0.9, titles = Bnames,
                          lines = TRUE, linesize = 1, ribbon = TRUE, ribbon_alpha = 0.1, 
                          points = FALSE, pt.size = 3, palette, ...) {
  require(ggplot2)
  
  # palette check
  user_palette <- !missing(palette)
  
  if (user_palette) {
    user_fills <- scale_fill_manual(values = palette)
    user_colors <- scale_color_manual(values = palette)
  }
  
  
  NS <- seq(1, 22, 3) 
  BS <- seq(2, 23, 3) 
  DS <- seq(3, 24, 3)
  lb <- (1 - conf.level)/2
  ub <- 1 - lb
  
  estimate <- mat.or.vec(24, length(Bnames))
  lower <- mat.or.vec(24, length(Bnames))
  upper <- mat.or.vec(24, length(Bnames))
  for(i in 1:length(Bnames)){
    estimate[,i] <- apply(get(Bnames[i]), 2, point.est)  
    lower[,i] <- apply(get(Bnames[i]), 2, quantile, probs = lb)
    upper[,i] <- apply(get(Bnames[i]), 2, quantile, probs = ub)
  }
  
  g_plots <- list()
  for(i in 1:length(Bnames)){
    df <- data.frame(
      congress = rep(73:80,3),
      region = rep(c("Deep South", "Border South", "Non-South"), each = 8),
      estimate = c(estimate[c(DS,BS,NS),i]),
      lower = c(lower[c(DS,BS,NS),i]),
      upper = c(upper[c(DS,BS,NS),i])) 
    
    #     assign(names[i], df)
    gg <- ggplot(df, aes(x = congress, y = estimate, ymin = lower, ymax = upper,
                         group = region, color = region, fill = region))
    
    if (lines) {
      gg <- gg + geom_line(size = linesize, ...)
    }
    if (ribbon) {
      gg <- gg + geom_ribbon(alpha = ribbon_alpha, size = 0, ...) 
    }
    if (points) {
      gg <- gg + geom_point(size = pt.size, ...)
    }
    
    gg <- gg + ggtitle(titles[i]) + scale_x_continuous(breaks = 73:80)
    
    if (user_palette) {
      g_plots[[i]] <- gg + user_fills + user_colors
    }
    else {
      g_plots[[i]] <- gg 
    }
    
  }
  
  return(g_plots)
}

