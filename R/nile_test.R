#' @import ggplot2
run_nile_example <- function() {

  # Use Nile dataset
  data(Nile)

  # Get x and y
  x <- start(Nile)[[1]]:end(Nile)[[1]]
  y <- as.vector(Nile)
  y <- y - mean(y)
  y <- y / sd(y)

  # OK -- let's write some GP equations
  # Find the best parameters
  test_points <- seq(min(x) - 100, max(x) + 100, length.out = 200)

  fit <- optimise_and_fit_rbf_gp(as.matrix(x), y, as.matrix(test_points))
  param_results <- fit[['hyperparameters']]
  full_results <- fit[['predictions']]

  hyperparams <- fit_marginal_likelihood_rbf(as.matrix(x), y)
  kernel_fun <- function(x1, x2) full_rbf_kernel(x1, x2, l = hyperparams[['l']],
                                                 tau = hyperparams[['tau']])
  diag_kernel_fun <- function(x1, x2) diag_rbf_kernel(x1, x2, l = hyperparams[['l']],
                                                      tau = hyperparams[['tau']])
  inducing <- seq(min(x), max(x), length.out = 10)
  dic_results <- calculate_dic(hyperparams[['sigma']], as.matrix(inducing),
                               as.matrix(x), y, as.matrix(test_points),
                               kernel_fun)

  dic_vars <- diag(dic_results$cov)

  # Try instead to draw samples from this
  dic_frame <- data.frame(year = test_points, 
                          means = dic_results$mean, 
                          vars = dic_vars, 
                          vars_with_noise = diag(dic_results$cov) +
                            param_results[['sigma']]^2,
                          method = 'DIC')

  exact_frame <- data.frame(year = test_points,
                            means = full_results$mean,
                            vars = diag(full_results$cov),
                            vars_with_noise = diag(full_results$cov) +
                              param_results[['sigma']]^2,
                            method = 'Exact')

  inducing_frame <- data.frame(x = inducing, y = 0)

  combined <- rbind(dic_frame, exact_frame)

  train_df <- data.frame(year = x,
                         value = y,
                         dataset = 'train')

  p <- ggplot(data = combined, aes(x = year, y = means, fill = method)) +
    geom_point() +
    geom_ribbon(aes(ymin = means - 2 * sqrt(vars),
                    ymax = means + 2 * sqrt(vars),
                    colour = 'Process noise'),
                alpha = 0.2) +
    geom_ribbon(aes(ymin = means - 2 * sqrt(vars_with_noise),
                    ymax = means + 2 * sqrt(vars_with_noise),
                    colour = 'Measurement noise'),
                alpha = 0.2) +
    geom_line() +
    theme_classic() +
    geom_point(data = train_df, aes(x = year, y = value, 
                                    colour = 'Training data'), 
               inherit.aes=FALSE) +
    geom_point(data = inducing_frame, aes(x = x, y = y, 
                                          colour = 'Inducing points'),
               inherit.aes=FALSE) +
    ggtitle('DIC vs. Exact GP') +
    ylab('Nile level relative to mean')

  print(p)

  # ggsave('dic_vs_exact.png', dpi = 300, width = 10, height = 6)

}
