# =========================
# MODULE: UTILITY FUNCTIONS=====================================
# Utility functions: conversions, normalization, data processing
# ==============================================================

#' @title Conversion barn → m2
#' @description Converts cross section from barn to square meters
#' @param b Cross section in barn
#' @return Cross section in m2
#' @details 1 barn = 10exp(-28) m2 = 10exp(-24) cm2
#' @examples
#' barn_to_m2(1)      # 1e-28 m2
#' barn_to_m2(100)    # 1e-26 m2
barn_to_m2 <- function(b) {
  if (!is.numeric(b)) {
    stop("Input must be numeric")
  }
  b * 1e-28
}

#' @title Conversion m2 → barn
#' @description Converts cross section from square meters to barn
#' @param m Cross section in m2
#' @return Cross section in barn
#' @details 1 m2 = 10e28 barn
#' @examples
#' m2_to_barn(1e-28)  # 1 barn
#' m2_to_barn(1e-24)  # 10000 barn
m2_to_barn <- function(m) {
  if (!is.numeric(m)) {
    stop("Input must be numeric")
  }
  m / 1e-28
}

#' @title Conversion barn → cm2
#' @description Converts cross section from barn to cm²
#' @param b Cross section in barn
#' @return Cross section in cm2
#' @details 1 barn = 10e-24 cm2
barn_to_cm2 <- function(b) {
  if (!is.numeric(b)) {
    stop("Input must be numeric")
  }
  b * 1e-24
}

#' @title Conversion cm2 → barn
#' @description Converts cross section from cm2 to barn
#' @param cm Cross section in cm2
#' @return Cross section in barn
#' @details 1 cm2 = 10e24 barn
cm2_to_barn <- function(cm) {
  if (!is.numeric(cm)) {
    stop("Input must be numeric")
  }
  cm / 1e-24
}

#' @title Spectrum Normalization
#' @description Normalizes a vector by its maximum
#' @param x Numeric vector to normalize
#' @param method Normalization method: "max" (default), "sum", "range"
#' @return Normalized vector
#' @details
#'   - "max": divides by maximum (result in [0, 1])
#'   - "sum": divides by sum (integral = 1)
#'   - "range": normalizes between 0 and 1
#' @examples
#' x <- c(2, 4, 6, 8, 10)
#' normalize(x)                    # [0.2, 0.4, 0.6, 0.8, 1.0]
#' normalize(x, method = "sum")    # [0.067, 0.133, 0.2, 0.267, 0.333]
#' normalize(x, method = "range")  # [0, 0.25, 0.5, 0.75, 1.0]
normalize <- function(x, method = "max") {
  if (!is.numeric(x)) {
    stop("Input must be numeric")
  }

  if (length(x) == 0) {
    return(x)
  }

  method <- match.arg(method, c("max", "sum", "range"))

  if (method == "max") {
    max_val <- max(x, na.rm = TRUE)
    if (max_val == 0) {
      warning("Maximum = 0, returning original vector")
      return(x)
    }
    return(x / max_val)

  } else if (method == "sum") {
    sum_val <- sum(x, na.rm = TRUE)
    if (sum_val == 0) {
      warning("Sum = 0, returning original vector")
      return(x)
    }
    return(x / sum_val)

  } else if (method == "range") {
    min_val <- min(x, na.rm = TRUE)
    max_val <- max(x, na.rm = TRUE)
    range_val <- max_val - min_val
    if (range_val == 0) {
      warning("Range = 0, all elements are identical")
      return(rep(0, length(x)))
    }
    return((x - min_val) / range_val)
  }
}

#' @title Create Energy Grid
#' @description Generates uniform or logarithmic energy grid
#' @param Emin Minimum energy
#' @param Emax Maximum energy
#' @param n Number of points (default: 100)
#' @param log If TRUE, logarithmic grid; if FALSE, linear grid
#' @return Energy vector
#' @details
#'   Log grid: useful for covering multiple orders of magnitude
#'   Linear grid: useful for studying a specific region
#' @examples
#' # Log grid to study from 0.1 eV to 10 MeV
#' E_log <- energy_grid(0.1, 1e7, n = 200, log = TRUE)
#'
#' # Linear grid for a resonance
#' E_lin <- energy_grid(2.8, 3.2, n = 100, log = FALSE)
energy_grid <- function(Emin, Emax, n = 100, log = TRUE) {
  if (Emin <= 0 && log) {
    stop("Emin must be > 0 for logarithmic grid")
  }

  if (Emin >= Emax) {
    stop("Emin must be < Emax")
  }

  if (n < 2) {
    stop("n must be >= 2")
  }

  if (log) {
    10^seq(log10(Emin), log10(Emax), length.out = n)
  } else {
    seq(Emin, Emax, length.out = n)
  }
}

#' @title Data Interpolation
#' @description Interpolates values from experimental data
#' @param x_exp Experimental x-values vector (must be sorted)
#' @param y_exp Experimental y-values vector
#' @param x_new New x-values for interpolation
#' @param method Interpolation method: "linear" or "spline"
#' @param extrap Allow extrapolation (default: FALSE)
#' @return Interpolated values at x_new points
#' @details
#'   - "linear": piecewise linear interpolation (fast)
#'   - "spline": cubic spline interpolation (smooth)
#'   If extrap = FALSE, values outside domain return NA
#' @examples
#' E_exp <- c(1, 2, 3, 4, 5)
#' sigma_exp <- c(2, 8, 6, 4, 2)
#' E_new <- seq(1, 5, by = 0.1)
#' sigma_interp <- interpolate_data(E_exp, sigma_exp, E_new)
interpolate_data <- function(x_exp, y_exp, x_new,
                             method = "linear",
                             extrap = FALSE) {

  if (length(x_exp) != length(y_exp)) {
    stop("x_exp and y_exp must have the same length")
  }

  if (length(x_exp) < 2) {
    stop("At least 2 points are required for interpolation")
  }

  method <- match.arg(method, c("linear", "spline"))

  if (method == "linear") {
    # Native R linear interpolation
    result <- approx(x_exp, y_exp, xout = x_new,
                     method = "linear",
                     rule = if(extrap) 2 else 1)$y

  } else if (method == "spline") {
    # Cubic spline interpolation
    result <- spline(x_exp, y_exp, xout = x_new,
                     method = "natural")$y

    # Handle extrapolation if not allowed
    if (!extrap) {
      out_of_domain <- x_new < min(x_exp) | x_new > max(x_exp)
      result[out_of_domain] <- NA
    }
  }

  return(result)
}

#' @title Smooth Noisy Data
#' @description Applies smoothing filter to data
#' @param x Vector to smooth
#' @param window_size Smoothing window size (odd number)
#' @param method Method: "moving_avg" or "savitzky_golay"
#' @return Smoothed vector
#' @details
#'   - "moving_avg": simple moving average
#'   - "savitzky_golay": polynomial filter (preserves peaks better)
smooth_data <- function(x, window_size = 5, method = "moving_avg") {
  if (!is.numeric(x)) {
    stop("Input must be numeric")
  }

  if (window_size %% 2 == 0) {
    window_size <- window_size + 1
    warning("window_size must be odd, adjusted to ", window_size)
  }

  method <- match.arg(method, c("moving_avg", "savitzky_golay"))

  n <- length(x)
  if (window_size > n) {
    stop("window_size cannot be larger than vector length")
  }

  half_window <- (window_size - 1) / 2
  x_smooth <- numeric(n)

  if (method == "moving_avg") {
    # Moving average
    for (i in 1:n) {
      idx_start <- max(1, i - half_window)
      idx_end <- min(n, i + half_window)
      x_smooth[i] <- mean(x[idx_start:idx_end], na.rm = TRUE)
    }

  } else if (method == "savitzky_golay") {
    # Simplified Savitzky-Golay filter (order 2)
    # Coefficients for window of 5: [-3, 12, 17, 12, -3] / 35
    if (window_size == 5) {
      coeff <- c(-3, 12, 17, 12, -3) / 35
    } else {
      # Fallback: moving average
      warning("Savitzky-Golay not implemented for window_size != 5, using moving_avg")
      return(smooth_data(x, window_size, method = "moving_avg"))
    }

    for (i in 1:n) {
      idx_start <- max(1, i - half_window)
      idx_end <- min(n, i + half_window)
      idx_range <- idx_start:idx_end

      if (length(idx_range) == window_size) {
        x_smooth[i] <- sum(x[idx_range] * coeff)
      } else {
        x_smooth[i] <- mean(x[idx_range], na.rm = TRUE)
      }
    }
  }

  return(x_smooth)
}

#' @title Numerical Derivative
#' @description Calculates derivative of sampled function
#' @param x Vector of x-values
#' @param y Vector of y-values
#' @param method Method: "forward", "backward", "central"
#' @return Vector of derivatives dy/dx
#' @details
#'   - "forward": (y[i+1] - y[i]) / (x[i+1] - x[i])
#'   - "backward": (y[i] - y[i-1]) / (x[i] - x[i-1])
#'   - "central": (y[i+1] - y[i-1]) / (x[i+1] - x[i-1]) (more accurate)
derivative <- function(x, y, method = "central") {
  if (length(x) != length(y)) {
    stop("x and y must have the same length")
  }

  n <- length(x)
  if (n < 2) {
    stop("At least 2 points are required")
  }

  method <- match.arg(method, c("forward", "backward", "central"))
  dy_dx <- numeric(n)

  if (method == "forward") {
    for (i in 1:(n-1)) {
      dy_dx[i] <- (y[i+1] - y[i]) / (x[i+1] - x[i])
    }
    dy_dx[n] <- dy_dx[n-1]  # Extrapolate last point

  } else if (method == "backward") {
    dy_dx[1] <- (y[2] - y[1]) / (x[2] - x[1])  # Forward for first
    for (i in 2:n) {
      dy_dx[i] <- (y[i] - y[i-1]) / (x[i] - x[i-1])
    }

  } else if (method == "central") {
    dy_dx[1] <- (y[2] - y[1]) / (x[2] - x[1])  # Forward for first
    for (i in 2:(n-1)) {
      dy_dx[i] <- (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
    }
    dy_dx[n] <- (y[n] - y[n-1]) / (x[n] - x[n-1])  # Backward for last
  }

  return(dy_dx)
}

#' @title Peak Finding
#' @description Detects local maxima in a signal
#' @param x Vector of values
#' @param threshold Relative threshold (0-1) for peak detection
#' @param min_distance Minimum distance between two peaks
#' @return Indices of detected peaks
#' @examples
#' E <- seq(0, 10, length.out = 500)
#' sigma <- breit_wigner(E, 3, 0.2, 10) +
#'          breit_wigner(E, 6, 0.3, 8)
#' peaks <- find_peaks(sigma, threshold = 0.3)
find_peaks <- function(x, threshold = 0.1, min_distance = 1) {
  if (!is.numeric(x)) {
    stop("Input must be numeric")
  }

  n <- length(x)
  if (n < 3) {
    return(integer(0))
  }

  # Normalize for threshold
  x_norm <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

  peaks <- c()

  for (i in 2:(n-1)) {
    # Check if it's a local maximum
    if (x[i] > x[i-1] && x[i] > x[i+1]) {
      # Check threshold
      if (x_norm[i] > threshold) {
        # Check minimum distance with already detected peaks
        if (length(peaks) == 0 || min(abs(i - peaks)) >= min_distance) {
          peaks <- c(peaks, i)
        }
      }
    }
  }

  return(peaks)
}

#' @title Spectrum Descriptive Statistics
#' @description Calculates statistics of a cross section
#' @param E Energy vector
#' @param sigma Cross section vector
#' @return List with: mean, median, max, max_E, integral, FWHM
spectrum_stats <- function(E, sigma) {
  if (length(E) != length(sigma)) {
    stop("E and sigma must have the same length")
  }

  # Basic statistics
  idx_max <- which.max(sigma)

  # Integral (area under curve)
  if (length(E) >= 2) {
    integral <- sum(diff(E) * (sigma[-1] + sigma[-length(sigma)]) / 2)
  } else {
    integral <- NA
  }

  # FWHM (Full Width at Half Maximum) - approximate
  half_max <- max(sigma) / 2
  idx_above <- which(sigma > half_max)
  if (length(idx_above) >= 2) {
    FWHM <- E[max(idx_above)] - E[min(idx_above)]
  } else {
    FWHM <- NA
  }

  list(
    mean = mean(sigma, na.rm = TRUE),
    median = median(sigma, na.rm = TRUE),
    max = max(sigma, na.rm = TRUE),
    max_E = E[idx_max],
    integral = integral,
    FWHM = FWHM,
    n_points = length(E)
  )
}
