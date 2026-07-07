# eivGP
For original code, see following link [https://github.com/ffpphh/eivGP](https://github.com/ffpphh/eivGP)  

Download through following command.

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install bdynets from GitHub
pak::pak("sj156/eivGP/eivGP", dependencies = TRUE)

## `install_github()` was deprecated in devtools 2.5.0. 
## devtools::install_github("sj156/bayesdir", subdir = "bayesdir", dependencies = TRUE) 
```

Very first raw version. For examples, see `script\5simulation_ordinal.Rmd`, where mixing report is nearby `script\5.1mixing_ordinal.Rmd`

1. The Examples may be too long to be shown here?
2. The Plots in mixing or just basic visualization, package it?
3. Is the notes good enough?
