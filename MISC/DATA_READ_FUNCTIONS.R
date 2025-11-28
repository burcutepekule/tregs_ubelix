
# Function to calculate Cohen's d
cohens_d = function(x, y) {
  nx = length(x)
  ny = length(y)
  mx = mean(x, na.rm = TRUE)
  my = mean(y, na.rm = TRUE)
  sx = sd(x, na.rm = TRUE)
  sy = sd(y, na.rm = TRUE)
  
  # Pooled standard deviation
  pooled_sd = sqrt(((nx - 1) * sx^2 + (ny - 1) * sy^2) / (nx + ny - 2))
  
  # Cohen's d
  d = (my - mx) / pooled_sd
  
  if (is.na(pooled_sd) || pooled_sd == 0) {
    return(0)  # or NA, depending on how you want to interpret it
  }
  return(d)
}

calculate_js_divergence_simple = function(vec1, vec2, n_bins = 50) {
  # Create matching histograms with the same breaks
  breaks = seq(min(c(vec1, vec2)),
                max(c(vec1, vec2)),
                length.out = n_bins + 1)
  
  hist1 = hist(vec1, breaks = breaks, plot = FALSE)
  hist2 = hist(vec2, breaks = breaks, plot = FALSE)
  
  # Convert counts to probabilities
  p = hist1$counts / sum(hist1$counts)
  q = hist2$counts / sum(hist2$counts)
  
  # Calculate JS divergence
  js_div = JSD(rbind(p, q), unit = "log2")
  
  return(js_div)
}


calculate_js_divergence = function(vec1, vec2, n_bins = NULL, method = "sturges") {
  
  # If n_bins not specified, calculate it
  if (is.null(n_bins)) {
    n1 = length(vec1)
    n2 = length(vec2)
    n = n1 + n2  # Total sample size
    
    n_bins = switch(method,
                     # Sturges' rule (default in R's hist)
                     "sturges" = ceiling(log2(n) + 1),
                     
                     # Scott's rule (assumes normal-like distribution)
                     "scott" = {
                       sd_val = sd(c(vec1, vec2))
                       if (sd_val == 0 || is.na(sd_val)) {
                         ceiling(log2(n) + 1)  # Fallback to Sturges
                       } else {
                         h = 3.5 * sd_val / (n^(1/3))
                         bins = ceiling((max(c(vec1, vec2)) - min(c(vec1, vec2))) / h)
                         if (is.na(bins) || bins <= 0 || is.infinite(bins)) {
                           ceiling(log2(n) + 1)  # Fallback
                         } else {
                           bins
                         }
                       }
                     },
                     
                     # Freedman-Diaconis rule (robust to outliers)
                     "fd" = {
                       iqr_val = IQR(c(vec1, vec2))
                       if (iqr_val == 0 || is.na(iqr_val)) {
                         ceiling(log2(n) + 1)  # Fallback to Sturges if IQR is zero - this happens when distributions are very very different
                       } else {
                         h = 2 * iqr_val / (n^(1/3))
                         bins = ceiling((max(c(vec1, vec2)) - min(c(vec1, vec2))) / h)
                         if (is.na(bins) || bins <= 0 || is.infinite(bins)) {
                           ceiling(log2(n) + 1)  # Fallback
                         } else {
                           bins
                         }
                       }
                     },
                     
                     # Square root rule
                     "sqrt" = ceiling(sqrt(n)),
                     
                     # Rice rule
                     "rice" = ceiling(2 * n^(1/3)),
                     
                     # Default
                     ceiling(log2(n) + 1)
    )
    
    # Ensure reasonable bounds
    n_bins = max(10, min(n_bins, 200))
  }
  
  # Create matching histograms with the same breaks
  breaks = seq(min(c(vec1, vec2)),
                max(c(vec1, vec2)),
                length.out = n_bins + 1)
  
  hist1 = hist(vec1, breaks = breaks, plot = FALSE)
  hist2 = hist(vec2, breaks = breaks, plot = FALSE)
  
  # Convert counts to probabilities
  p = hist1$counts / sum(hist1$counts)
  q = hist2$counts / sum(hist2$counts)
  
  # Calculate JS divergence
  js_div = suppressMessages(JSD(rbind(p, q), unit = "log2"))
  
  return(c(js_div, n_bins))
}

steady_state_idx = function(x, k = 20, tail_frac = 0.25,
                            tol_abs = 0.05*(150-25),     # 2.5 units
                            tol_sd  = 0.02*(150-25),# 1.25 units
                            tol_slope = 0.005*(150-25))  # 0.125/step
{
  n      = length(x)
  tail_n = ceiling(n*tail_frac)
  x_asym = mean(tail(x, tail_n))
  
  m    = rollapply(x, k, mean, align = "right", fill = NA)
  sd   = rollapply(x, k, sd,   align = "right", fill = NA)
  sl   = rollapply(x, k, function(v) mean(diff(v)), align = "right", fill = NA)
  cand = which(abs(m-x_asym) <= tol_abs & sd <= tol_sd & abs(sl) <= tol_slope)
  
  if (length(cand) == 0) return(NA_integer_)
  cand[1]
}