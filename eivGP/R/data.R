# Generated from _main.Rmd: do not edit by hand

#' Annotated Branin Function for Nominal Case Simulation
#' 
#' Acts as the oracle for the simulation study in Section 5.2.
#' 
#' @param u1,u2 Latent continuous coordinates.
#' @return A list with response 'f' and category 'c'.
#' @details 
#' Maps U to 5 regions: 
#' Region 1-2: Separated by linear boundaries.
#' Region 3-5: Partitioned via Voronoi diagram based on Branin's global minima.
#' @export
Branin <- function (x1, x2, 
                    x1.ub = 2, x1.lb = -2, x2.ub = 2, x2.lb = -2, 
                    a = 1, b = 5.1/(4*pi^2), c = 5/pi, s = 10, r = 6, t = 1/(8*pi)) {
  # rescale 
  x1 <- matrix((x1 - x1.lb)/(x1.ub - x1.lb)*(10 - (-5)) + (-5), ncol = 1)  #[-5, 10]
  x2 <- matrix((x2 - x2.lb)/(x2.ub - x2.lb)*(15 - 0) + 0, ncol = 1)  #[0, 15]
  x <- cbind(x1, x2)
  
  # f value
  f <- c(a*(x2 - b*x1^2 + c*x1 - r)^2 + s*(1-t)*cos(x1) + s)
  
  # linear boundary
  # line 1
  p11 <- c(10, 8)
  p12 <- c(-1, 15)
  slope1 <- (p12[2] - p11[2])/(p12[1] - p11[1])
  intp1 <- p11[2] - slope1*p11[1]
  # line 2
  p21 <- c(-5, 10)
  p22 <- c(0, 0)
  slope2 <- (p22[2] - p21[2])/(p22[1] - p21[1])
  intp2 <- p21[2] - slope2*p21[1]

  # classify
  d1 <- x2 - (intp1 + slope1*x1)
  d2 <- x2 - (intp2 + slope2*x1)
  c <- rep(0, length(f))
  c[d1 > 0] <- 1
  c[(d1 <= 0) & (d2 < 0)] <- 2
  
  # Voronoi diagram
  x.min1 <- c(-pi, 12.275)
  x.min2 <- c(pi, 2.275)
  x.min3 <- c(9.42478, 2.475)
  d11 <- apply(x, 1, function(x) sum((x-x.min1)^2))
  d22 <- apply(x, 1, function(x) sum((x-x.min2)^2))
  d33 <- apply(x, 1, function(x) sum((x-x.min3)^2))
  
  c[(d1 <= 0) & (d2 >= 0)] <- 2 + apply(cbind(d11,d22,d33)[(d1 <= 0) & (d2 >= 0), ,drop = FALSE], 1, which.min)
  
  # return
  return(list(f = f, c = c))
}

# plot configurations
# my_theme <- 
  # theme_minimal() +
  # theme(
    # plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    # axis.title = element_text(size = 17, face = "bold"),
    # axis.text = element_text(size = 15),      
    # legend.title = element_text(size = 17, face = "bold"),
    # legend.text = element_text(size = 15))



#' Plot Scatterplots of Observed and Latent Data
#' 
#' @param c Vector of qualitative inputs (category labels).
#' @param X Matrix/Vector of quantitative inputs (assumed 1D for visualization).
#' @param y Vector of observed responses.
#' @param U Vector of true latent continuous variables.
#' @param u.obs.idx Indices of the samples where the latent variable U is observed.
#' 
#' @return A list containing two ggplot2 objects:
#'   \item{ggfig1}{Scatterplot of y vs (x, c). Corresponds to Figure 5(a). 
#'     Shows the data as seen by a standard mixed-input model.}
#'   \item{ggfig2}{Scatterplot of y vs observed (x, u). Corresponds to Figure 5(b). 
#'     Shows the "expensive" ground truth structure used as anchors.}
#'     
Plot.Data <- function (c, X, y, U, u.obs.idx) {
  
  my_theme <- 
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 17, face = "bold"),
    axis.text = element_text(size = 15),      
    legend.title = element_text(size = 17, face = "bold"),
    legend.text = element_text(size = 15))
  
  
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
