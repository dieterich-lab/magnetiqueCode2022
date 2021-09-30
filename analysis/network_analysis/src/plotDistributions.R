plotDistributions <- function(data, path=getwd(), alpha=0.3, width=7, height=7, 
                              xsize=5, xangle=90, ysize=10, yangle=0){
  
  library(reshape)
  mm <- melt(data = data)
  
  #
  pdf(file = paste0(path, "/density_plot.pdf"), width = width, height = height)
  p <- ggplot(mm, aes(x = value, fill = variable)) + 
    geom_density(alpha = alpha) + 
    theme(axis.text.x=element_text(size=xsize, angle=xangle), axis.text.y=element_text(size=ysize, angle = yangle)) +
    ggtitle("Density plot of samples")
  plot(p)
  dev.off()
  
  # #
  # pdf(file = paste0(path, "/cloud_point.pdf"), width = width, height = height)
  # p <- ggplot(mm, aes(x = 1:nrow(mm), y = value, group = variable, color = variable)) + geom_point() + 
  #   theme(axis.text.x=element_text(size=xsize, angle=xangle), axis.text.y=element_text(size=ysize, angle = yangle)) +
  #   ggtitle("Density plot of samples")
  # plot(p)
  # dev.off()
  
  #
  pdf(file = paste0(path, "/box_plot.pdf"), width = width, height = height)
  p <- ggplot(mm, aes(x = 1:nrow(mm), y = value, group = variable, color = variable)) + geom_boxplot() + 
    theme(axis.text.x=element_text(size=xsize, angle=xangle), axis.text.y=element_text(size=ysize, angle = yangle)) +
    ggtitle("Density plot of samples")
  plot(p)
  dev.off()
  
}