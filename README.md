# JMtwostage

`JMtwostage` is an R package that provides a collection of functions for **two-stage joint modeling** of longitudinal and survival data using various methods. It is designed for flexible modeling and imputation in survival analysis with time-dependent covariates and missing data.

## ✨ Features

- Implements multiple approaches for handling missing data:
  - Complete Case Analysis
  - Last Observation Carried Forward (LOCF)
  - Multiple Imputation (MI)
  - Inverse Probability Weighting (IPW)
  
- Supports modeling of **time-dependent covariates** in survival analysis.
- Provides tools for:
  - Model fitting
  - Model evaluation
  - Imputation of longitudinal biomarkers
- Useful in joint modeling of longitudinal and time-to-event data.

## 📦 Installation

You can install the package directly from GitHub using the `remotes` package:

```r
# Install remotes if not already installed
install.packages("remotes")

# Install JMtwostage from GitHub
remotes::install_github("kumarbhrigu/JMtwostage")
