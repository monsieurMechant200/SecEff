# ===================
# MODULE 4: VIZ_TOOLS========================================
# Advanced visualization functions for nuclear cross sections
# ===========================================================

#' @title Generic Cross-Section Plot
#'
#' @description
#' Produces a standard line plot of a nuclear cross section \eqn{\sigma(x)} as
#' a function of any variable (energy, angle, wavelength, etc.). The plot can be
#' displayed on a linear or logarithmic scale.
#'
#' @param x Numeric vector. The variable over which the cross section is defined.
#'          Typically energy (MeV), angle (degrees), or wavelength.
#' @param sigma Numeric vector of the same length as `x`. Cross-section values
#'          expressed in barns.
#' @param xlab Character. Label for the x-axis.
#' @param ylab Character. Label for the y-axis.
#' @param main Character. Title of the plot.
#' @param col Color for the curve.
#' @param log_scale Logical. If `TRUE`, the plot uses logarithmic scales on both axes.
#'
#' @return A base R plot is generated. No value is returned.
#'
#' @examples
#' E <- seq(0.01, 10, length.out = 200)
#' sigma <- exp(-E) * 5
#' plot_sigma(E, sigma, main = "Example Cross Section")
#'
#' @export
plot_sigma <- function(x, sigma, xlab = "Variable", ylab = "Cross section (barn)",
                       main = "Cross section", col = "darkblue",
                       log_scale = FALSE) {
  if (log_scale) {
    plot(x, sigma, type = "l", lwd = 2, log = "xy",
         xlab = xlab, ylab = ylab, main = main, col = col)
  } else {
    plot(x, sigma, type = "l", lwd = 2,
         xlab = xlab, ylab = ylab, main = main, col = col)
  }
  grid(col = "gray80")
}



#' @title Resonance Curve Plot
#'
#' @description
#' Plots a nuclear resonance peak \eqn{\sigma(E)} as a function of incident
#' particle energy. The function automatically marks the maximum resonance value.
#'
#' @param E Numeric vector of energies (MeV).
#' @param sigma Numeric vector of cross-section values (barns).
#' @param title Character. Title of the plot.
#'
#' @return A resonance plot with the maximum resonance point highlighted.
#'
#' @examples
#' E <- seq(0, 2, length.out = 300)
#' sigma <- dnorm(E, mean = 1, sd = 0.1) * 100
#' plot_resonance(E, sigma)
#'
#' @export
plot_resonance <- function(E, sigma, title = "Resonance curve") {
  plot(E, sigma, type = "l", col = "blue", lwd = 2,
       xlab = "Energy (MeV)", ylab = "σ(E) (barn)",
       main = title)
  grid(col = "gray80")

  # Mark the maximum
  idx_max <- which.max(sigma)
  points(E[idx_max], sigma[idx_max], pch = 19, col = "red", cex = 1.5)
  text(E[idx_max], sigma[idx_max],
       labels = sprintf("Max: %.2f barn", sigma[idx_max]),
       pos = 3, col = "red")
}



#' @title Compare Multiple Cross Sections
#'
#' @description
#' Displays several cross-section curves \eqn{\sigma_i(E)} on the same graph
#' for comparison. Supports linear or logarithmic axes.
#'
#' @param E Numeric vector of energies.
#' @param sigma_list List of numeric vectors. Each element is a cross-section curve.
#' @param labels Character vector. Names of each curve.
#' @param colors Optional vector of colors.
#' @param main Title of the plot.
#' @param log_scale If TRUE, uses logarithmic axes.
#'
#' @return A multi-curve comparison plot.
#'
#' @examples
#' E <- seq(0.01, 5, length.out = 200)
#' sigma1 <- 10 * exp(-E)
#' sigma2 <- 5 * exp(-0.5 * E)
#' plot_compare_sigma(E, list(sigma1, sigma2), labels = c("Material A", "Material B"))
#'
#' @export
plot_compare_sigma <- function(E, sigma_list, labels, colors = NULL,
                               main = "Cross-section comparison",
                               log_scale = FALSE) {

  if (is.null(colors)) colors <- rainbow(length(sigma_list))

  ylim <- range(unlist(sigma_list), na.rm = TRUE)

  if (log_scale) {
    plot(E, sigma_list[[1]], type = "l", lwd = 2, col = colors[1],
         xlab = "Energy (MeV)", ylab = "σ(E) (barn)",
         main = main, log = "xy", ylim = ylim)
  } else {
    plot(E, sigma_list[[1]], type = "l", lwd = 2, col = colors[1],
         xlab = "Energy (MeV)", ylab = "σ(E) (barn)",
         main = main, ylim = ylim)
  }

  for (i in 2:length(sigma_list)) {
    lines(E, sigma_list[[i]], lwd = 2, col = colors[i])
  }

  grid(col = "gray80")
  legend("topright", legend = labels, col = colors, lwd = 2, bg = "white")
}



#' @title Angular Differential Cross-Section Plot
#'
#' @description
#' Visualizes a differential cross section \eqn{\frac{d\sigma}{d\Omega}(\theta)} either
#' in polar coordinates or in standard Cartesian coordinates.
#'
#' @param theta Numeric vector of angles (degrees).
#' @param sigma_theta Numeric vector of differential cross-section values (barn/sr).
#' @param type Either `"polar"` or `"cartesian"`.
#' @param main Title of the plot.
#'
#' @return A polar or Cartesian angular distribution plot.
#'
#' @examples
#' theta <- seq(0, 360, by = 1)
#' sigma <- 3 + sin(theta * pi/180)
#' plot_angular_distribution(theta, sigma, type = "polar")
#'
#' @export
plot_angular_distribution <- function(theta, sigma_theta,
                                      type = "polar",
                                      main = "Angular distribution") {

  if (type == "polar") {
    theta_rad <- theta * pi / 180
    x <- sigma_theta * cos(theta_rad)
    y <- sigma_theta * sin(theta_rad)

    plot(x, y, type = "l", lwd = 2, col = "darkgreen",
         xlab = "σ·cos(θ)", ylab = "σ·sin(θ)",
         main = main, asp = 1)
    lines(c(0, 0), range(y), lty = 2, col = "gray50")
    lines(range(x), c(0, 0), lty = 2, col = "gray50")

  } else {
    plot(theta, sigma_theta, type = "l", lwd = 2, col = "darkgreen",
         xlab = "Angle (degrees)", ylab = "dσ/dΩ (barn/sr)",
         main = main)
    grid(col = "gray80")
  }
}



#' @title 2D Cross-Section Heatmap
#'
#' @description
#' Generates a color heatmap of a cross section depending simultaneously on
#' energy and angle: \eqn{\sigma(E, \theta)}.
#'
#' @param E Vector of energies (x-axis).
#' @param theta Vector of angles (y-axis).
#' @param sigma_matrix Matrix of cross-section values.
#' @param main Title of the plot.
#'
#' @return A 2D heatmap with contour lines.
#'
#' @export
plot_heatmap_sigma <- function(E, theta, sigma_matrix,
                               main = "2D cross section") {
  image(E, theta, sigma_matrix,
        xlab = "Energy (MeV)", ylab = "Angle (degrees)",
        main = main, col = heat.colors(50))
  contour(E, theta, sigma_matrix, add = TRUE, col = "black", lwd = 1)
}



#' @title Full Cross-Section Dashboard
#'
#' @description
#' Creates a 2×2 panel dashboard including:
#' - linear scale cross section
#' - logarithmic scale cross section
#' - experimental residuals (if provided)
#' - summary statistics
#'
#' Useful for nuclear data exploration and model validation.
#'
#' @param E Energy vector.
#' @param sigma Theoretical cross section.
#' @param sigma_exp Experimental cross section (optional).
#' @param E_exp Energies corresponding to `sigma_exp`.
#'
#' @return A 4-panel dashboard plot.
#'
#' @export
plot_dashboard <- function(E, sigma, sigma_exp = NULL, E_exp = NULL) {
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

  # 1. Linear scale
  plot_sigma(E, sigma, main = "Linear scale")
  if (!is.null(sigma_exp)) {
    points(E_exp, sigma_exp, pch = 19, col = "red", cex = 0.8)
  }

  # 2. Logarithmic scale
  plot_sigma(E, sigma, main = "Logarithmic scale", log_scale = TRUE)
  if (!is.null(sigma_exp)) {
    points(E_exp, sigma_exp, pch = 19, col = "red", cex = 0.8)
  }

  # 3. Residuals
  if (!is.null(sigma_exp)) {
    sigma_interp <- interpoler_sigma(E, sigma, E_exp)
    residuals <- sigma_exp - sigma_interp
    plot(E_exp, residuals, pch = 19, col = "purple",
         xlab = "Energy (MeV)", ylab = "Residuals (barn)",
         main = "Residuals exp - theory")
    abline(h = 0, lty = 2, col = "gray50")
    grid(col = "gray80")
  } else {
    plot(E, sigma, type = "n", main = "No experimental data")
    text(mean(E), mean(sigma), "Experimental data\nnot available")
  }

  # 4. Statistics panel
  plot.new()
  text(0.5, 0.8, "STATISTICS", cex = 1.5, font = 2)
  text(0.5, 0.6, sprintf("σ max: %.2e barn", max(sigma)), cex = 1.1)
  text(0.5, 0.5, sprintf("σ min: %.2e barn", min(sigma)), cex = 1.1)
  text(0.5, 0.4, sprintf("σ mean: %.2e barn", mean(sigma)), cex = 1.1)
  text(0.5, 0.3, sprintf("Points: %d", length(E)), cex = 1.1)

  par(mfrow = c(1, 1))
}



#' @title Time-Evolution Animation (Frame Export)
#'
#' @description
#' Generates a sequence of plots showing the evolution of a cross-section curve
#' over time. Useful for Monte-Carlo simulations, iterative models, and
#' parameter sweeps.
#'
#' If `output_dir` is provided, frames are saved as PNG files; otherwise, the
#' function displays frames interactively with a short delay.
#'
#' @param x_list List of x-vectors.
#' @param sigma_list List of cross-section vectors.
#' @param labels Character. Title base for the animation.
#' @param output_dir Directory in which to save frames (optional).
#'
#' @return PNG files or interactive animation.
#'
#' @export
plot_evolution <- function(x_list, sigma_list, labels,
                           output_dir = NULL) {

  n_frames <- length(sigma_list)

  for (i in 1:n_frames) {

    if (!is.null(output_dir)) {
      png(sprintf("%s/frame_%03d.png", output_dir, i),
          width = 800, height = 600)
    }

    plot(x_list[[i]], sigma_list[[i]], type = "l", lwd = 2, col = "blue",
         xlab = "Variable", ylab = "Cross section (barn)",
         main = sprintf("%s - Frame %d/%d", labels, i, n_frames),
         ylim = range(unlist(sigma_list)))
    grid(col = "gray80")

    if (!is.null(output_dir)) {
      dev.off()
    } else {
      Sys.sleep(0.2)
    }
  }
}
