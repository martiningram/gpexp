#' Compute an RBF Kernel
#'
#' This function computes the full RBF kernel for two vectors.
#'
#' @param v1 The first vector to compute distances for.
#' @param v2 The second vector to compute distances for.
#' @param l The length scale of the RBF kernel.
#' @param tau The standard deviation of the kernel.
#' @return A matrix containing the distances between each of the elements in v1
#' and v2.
full_rbf_kernel <- function(v1, v2, l, tau, jitter=1e-8) {
  # This is the non-squared distance matrix.
  if (identical(v1, v2)) {
    distances <- dist(v1)
  }
  else {
    distances <- pdist::pdist(v1, v2)
  }

  distances <- as.matrix(distances)
  kernel_matrix <- exp(-distances^2 / (2 * l^2))
  return(tau^2 * kernel_matrix + diag(rep(jitter, nrow(kernel_matrix)))
}

#' Compute a diagonal RBF kernel
#'
#' This function computes the diagonal RBF kernel for two vectors.
#'
#' @param v1 The first vector to compute distances for.
#' @param v2 The second vector to compute distances for.
#' @param l The length scale of the RBF kernel.
#' @param tau The standard deviation of the kernel.
#' @return A vector containing the diagonal elements of the RBF kernel matrix
#' between v1 and v2.
diag_rbf_kernel <- function(v1, v2, l, tau) {
  if (identical(v1, v2)) {
    return (tau^2 * rep(1, dim(v1)[1]))
  }
  else {
    # Make them conformable
    max_dim <- min(nrow(v1), nrow(v2))
    v1 <- v1[1:max_dim, ]
    v2 <- v2[1:max_dim, ]
    
    # Calculate the row-wise square distances
    row_diffs <- exp(-((v1 - v2) %*% (v1 - v2)) / (2 * l^2))
    return (tau^2 * row_diffs)
  }
}

#' Compute an ARD Kernel
#'
#' This function computes the ARD kernel for two data matrices.
#'
#' @param x1 Matrix of shape (N1, D).
#' @param x2 Matrix of shape (N2, D).
#' @param sigma_inv Diagonal matrix of shape (D, D) with the weights for each
#' dimension D.
#' @return The (N1, N2) kernel between each row of x1 and x2.
#' @export
ard_kernel <- function(x1, x2, sigma_inv) {

    # First, let's weight each dimension by its sigma.
    diag_sigma_inv <- diag(sigma_inv)
    x1_re_weighted <- sweep(x1, 2, diag_sigma_inv, 
                            function (x, y) x * sqrt(y))
    x2_re_weighted <- sweep(x2, 2, diag_sigma_inv, 
                            function (x, y) x * sqrt(y))
    
    if (identical(x1, x2)) {
        if (nrow(x1) == 1) {
            return(as.matrix(1))
        }
        distances <- as.matrix(dist(x1_re_weighted))^2
    } else {
        distances <- as.matrix(pdist::pdist(x1_re_weighted, x2_re_weighted))^2
    }
                            
    return(exp(-0.5 * distances))
                            
}


#' Compute an inverse using the Cholesky decomposition
#'
#' This function computes the inverse of a positive definite matrix using the
#' Cholesky decomposition.
#'
#' @param matrix_to_invert The matrix to find the inverse for.
#' @return The inverse of the matrix given.
cholesky_inverse <- function(matrix_to_invert) {

  return(chol2inv(chol(matrix_to_invert)))

}

predict_using_inverse <- function(x_new, x_train, y, inverse, kernel_fun) {

  k_star_star <- kernel_fun(x_new, x_new)
  new_with_train <- kernel_fun(x_new, x_train)
  times_inv <- new_with_train %*% inverse

  predicted_mean <- times_inv %*% y
  var_accounted_for <- times_inv %*% t(new_with_train)
  predicted_cov <- k_star_star - var_accounted_for

  return(list('mean' = predicted_mean, 'cov' = predicted_cov))

}

#' Fit an exact GP and predict a set of points
#'
#' This function fits and predicts a GP using the exact expressions for mean and
#' covariance.
#'
#' @param x_train The training covariates.
#' @param x_new The covariates to use for prediction,
#' @param sigma_noise The observation noise.
#' @param y The training observations.
#' @param kernel_fun The kernel to use.
#' @return A list containing entries `mean` and `cov` specifying the mean and
#' covariance at the points given by `x_new`.
#' @export
predict_points <- function(x_train, x_new, sigma_noise, y, kernel_fun,
                           mean_centre = TRUE, marginals_only = FALSE) {

  if (mean_centre) {
    # Standardise y first
    mean_y <- mean(y)
    y <- y - mean_y
  }

  # Compute the main inverse
  training_part <- kernel_fun(x_train, x_train)
  diag(training_part) <- diag(training_part) + sigma_noise^2
  inverse <- chol2inv(chol(training_part))

  if (marginals_only) {

    # We can do this in batches.
    batch_size <- 256
    assignments <- seq(nrow(x_new)) %/% batch_size
    predictions <- list('mean' = c(), 'cov' = c())

    for (cur_assignment in unique(assignments)) {
      cur_x_new <- as.matrix(x_new[assignments == cur_assignment, ])
      cur_predictions <- predict_using_inverse(cur_x_new, x_train, y, inverse,
                                               kernel_fun)
      # Pick out only the diagonal
      diag_cov <- diag(cur_predictions[['cov']])
      predictions[['cov']] <- c(predictions[['cov']], diag_cov)
      predictions[['mean']] <- c(predictions[['mean']],
                                 cur_predictions[['mean']])
    }

    # Turn the diagonal covariance matrix into a proper matrix
    new_size <- length(predictions[['cov']])
    predictions[['cov']] <- diag(predictions[['cov']], nrow=new_size,
                                 ncol=new_size)

  } else {

    predictions = predict_using_inverse(x_new, x_train, y, inverse, kernel_fun)

  }

  predicted_mean <- predictions[['mean']]
  predicted_cov <- predictions[['cov']]

  if (mean_centre) {
    # Add mean back on
    predicted_mean <- predicted_mean + mean_y
  }

  return(list('mean' = predicted_mean, 'cov' = predicted_cov))
}

#' Fits the hyperparameters of a single RBF kernel by mAP estimation
#'
#' This function optimises the hyperparameters, as well as the observation
#' noise, for an RBF kernel.
#'
#' @param x_train The training covariates.
#' @param y_train The training observations.
#' @param start_sigma The initial observation noise to use.
#' @param start_l The initial length scale to use.
#' @param start_tau The initial process standard deviation to use.
#' @return A list containing entries `sigma`, `l` and `tau` -- the optimised
#' hyperparameters of the RBF kernel.
fit_marginal_likelihood_rbf <- function(x_train, y_train, start_sigma = 10,
                                        start_l = 10, start_tau = 10) {

  y_train <- y_train - mean(y_train)

  to_optimize <- function(sigma, l, tau) {
    kernel <- full_rbf_kernel(x_train, x_train, tau = tau, l = l)
    diag(kernel) <- diag(kernel) + sigma^2

    signed_det <- determinant(kernel, logarithm = TRUE)
    first_part <- signed_det[['sign']] * signed_det[['modulus']]
    second_part <- t(y_train) %*% chol2inv(chol(kernel)) %*% y_train

    # We want to minimize, so it's the sum of these
    return(first_part + second_part)
  }

  vector_wrapper <- function(x) {
    return(to_optimize(sigma = x[1], l = x[2], tau = x[3]))
  }

  start_par <- c(start_sigma, start_l, start_tau)

  fit_result <- optim(start_par, vector_wrapper)

  params <- list('sigma' = fit_result$par[1],
                 'l' = fit_result$par[2],
                 'tau' = fit_result$par[3])

  return(params)
}

#' Optimises and fits an RBF Gaussian Process.
#'
#' This function optimises an RBF kernel using mAP estimation, trains the GP,
#' and predicts a set of new points.
#'
#' @param x_train The training covariates.
#' @param y_train The training observations.
#' @param x_new The covariates to predict for.
#' @param start_sigma The initial observation noise to use.
#' @param start_l The initial length scale to use.
#' @param start_tau The initial process standard deviation to use.
#' @return A list containing entries `predictions` and `hyperparameters`,
#' containing the `mean` and `cov` of the process in `predictions`, and the
#' kernel's hyperparameters in `hyperparameters`.
optimise_and_fit_rbf_gp <- function(x_train, y_train, x_new, start_sigma = 10,
                                    start_l = 10, start_tau = 10) {

  # Fit the kernel hyperparameters
  param_results <- fit_marginal_likelihood_rbf(x_train, y_train, start_sigma =
                                               start_sigma, start_tau =
                                               start_tau, start_l = start_l)

  # Curry the kernel
  kernel_fun <- function(x1, x2) full_rbf_kernel(x1, x2, 
                                                 l = param_results[['l']],
                                                 tau = param_results[['tau']])

  predictions <- predict_points(x_train, x_new, param_results[['sigma']],
                                y_train, kernel_fun)

  return(list('predictions' = predictions,
              'hyperparameters' = param_results))

}

#' @import ggplot2
#' @export
plot_gp <- function(x, mean, covariance, sigma, x_train = NULL, y_train = NULL,
                    save_to = NULL) {
  
  vars <- diag(covariance)
  
  plot_frame <- data.frame(x = x,
                           means = mean,
                           vars = diag(covariance),
                           vars_with_noise = diag(covariance) + sigma^2)
  
  p <- ggplot(data = plot_frame, aes(x = x, y = means)) +
    geom_point(aes(colour = 'Mean')) +
    geom_ribbon(aes(ymin = means - 2 * sqrt(vars),
                    ymax = means + 2 * sqrt(vars),
                    fill = 'Process noise'),
                alpha = 0.2) +
    geom_ribbon(aes(ymin = means - 2 * sqrt(vars_with_noise),
                    ymax = means + 2 * sqrt(vars_with_noise),
                    fill = 'Measurement noise'),
                alpha = 0.2) +
    geom_line(aes(colour = 'Mean')) +
    theme_classic()
  
  if (!(is.null(x_train)) & !(is.null(y_train))) {
    train_df <- data.frame(x = x_train,
                           y = y_train)
    
    p <- p +     
      geom_point(data = train_df, aes(x = x, y = y, 
                                      colour = 'Observations'), 
                            inherit.aes=FALSE)
  }

  if (!is.null(save_to)) {
    # Save the plot to the file specified
    ggsave(save_to, height=4, width=8, dpi=300)
  }
  
  return(p)
  
}
