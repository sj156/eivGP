# plot configurations
my_theme <- 
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 17, face = "bold"),
    axis.text = element_text(size = 15),      
    legend.title = element_text(size = 17, face = "bold"),
    legend.text = element_text(size = 15))

# plot the scatter plots of y versus (x,c) and y versus obversed (x,u) (only for 1d inputs)
Plot.Data <- function (c, X, y, U, u.obs.idx) {
  data_df <- data.frame(c = c, X = X, y = y)
  ggfig1 <- ggplot(data_df, aes(x = c, y = X, color = y)) +
    geom_point(size = 5) +
    scale_x_continuous(breaks = 1:6) +
    scale_color_viridis_c(option = 'inferno', limits = c(min(y), max(y))) +
    labs(color = "y") +
    my_theme
  
  data_df <- data.frame(u = U, X = X, y = y)
  ggfig2 <- ggplot(data_df[u.obs.idx, ], aes(x = u, y = X, color = y)) +
    geom_point(size = 5) +
    #scale_color_distiller(palette = "RdYlBu", limits = c(min(y), max(y))) +
    scale_color_viridis_c(option = 'inferno', limits = c(min(y), max(y))) +
    labs(color = "y") +
    my_theme
  
  return(list(ggfig1 = ggfig1, ggfig2 = ggfig2))
}



