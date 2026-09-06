library(ggplot2)
library(kernlab)  # For Gaussian process regression

set.seed(123)

# Generate synthetic data
x <- seq(0, 10, length.out = 20)
y <- sin(x) + rnorm(length(x), sd = 0.3)
data <- data.frame(x = x, y = y)

# Fit Gaussian process
gp <- gausspr(x = data$x, y = data$y, kernel = "rbfdot")

# Create prediction grid
x_new <- seq(-2, 12, length.out = 100)
pred <- predict(gp, newdata = as.data.frame(x_new), type = "sdev")
mean_pred <- pred$mean
sd_pred <- pred$sd

# Create plot
ggplot() +
  geom_ribbon(aes(x = x_new, ymin = mean_pred - 1.96*sd_pred, 
                  ymax = mean_pred + 1.96*sd_pred), 
              fill = "skyblue", alpha = 0.3) +
  geom_line(aes(x = x_new, y = mean_pred), color = "blue", linewidth = 1) +
  geom_point(aes(x = x, y = y), data = data, size = 3, color = "red") +
  geom_line(aes(x = x_new, y = sin(x_new)), color = "darkred", linetype = "dashed") +
  labs(title = "Gaussian Process Regression Example",
       x = "Input (x)", y = "Output (y)",
       caption = "Red points: Observed data | Blue line: GP mean | Shaded: 95% CI") +
  theme_minimal()

# Save plot
ggsave("gp_example.png", width = 8, height = 5, dpi = 300)
