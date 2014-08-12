# define make_gg_densities function --------------------------------------
make_gg_densities <- function(ests, plot_title){
  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  
  congress.labs <- 73:80
  NS <- seq(1, 22, 3) 
  BS <- seq(2, 23, 3) 
  DS <- seq(3, 24, 3)
  colnames(ests)[DS] <- paste0("DS_",1:8)
  colnames(ests)[BS] <- paste0("BS_",1:8)
  colnames(ests)[NS] <- paste0("NS_",1:8)

  melt.df <- list()
  g_plots <- list()
  for (i in 1:8) {  
    melt.df[[i]] <- melt(ests[,paste0(c("DS_","BS_","NS_"),i)])
    
    g_plots[[i]] <- (ggplot(melt.df[[i]], aes(x = value, group = Var2, color = Var2)) 
                     + geom_density()
                     + scale_color_manual(
                       name = "Region",
                       labels = c("Deep South", "Border South", "Non-South"),
                       values = c("red","purple2", "blue"))
#                      + ggtitle(paste("Year = ", period.labs[i]))
                     + xlab(paste("Congress = ", congress.labs[i]))
                     + ylab("")
                     + theme(plot.title = element_text(size = 13), 
                             legend.position = "none"))
  }
  
  legend <- (g_plots[[1]] + theme(legend.position = "bottom", 
                                  legend.direction = "horizontal",
                                  legend.key.size = unit(0.75, "cm"),
                                  legend.title = element_text(size = 14)))
  
  g_legend<-function(gg){
    g_tab <- ggplot_gtable(ggplot_build(gg))
    leg <- which(sapply(g_tab$grobs, function(k) k$name) == "guide-box")
    legend <- g_tab$grobs[[leg]]
    return(legend)
  }
  my.legend <- g_legend(legend)
  
  g_grob <- arrangeGrob(g_plots[[1]],
                        g_plots[[2]],
                        g_plots[[3]],
                        g_plots[[4]],
                        g_plots[[5]],
                        g_plots[[6]],
                        g_plots[[7]],
                        g_plots[[8]],
                        nrow=3, 
                        main = textGrob(plot_title, 
                                        gp=gpar(fontsize=14,font=2)))
  
  g_grid <- grid.arrange(g_grob, my.legend, nrow = 2, heights = c(10,1))
}
