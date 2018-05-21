tol <- 1e-6

calculate_q <- function(kernel_fun, inducing, a, b, tol = 1e-6) {

  k_au <- kernel_fun(a, inducing)
  k_ub <- kernel_fun(inducing, b)
  k_uu <- kernel_fun(inducing, inducing) + diag(nrow(inducing)) * tol
  k_uu_inv <- cholesky_inverse(k_uu)

  return(k_au %*% k_uu_inv %*% k_ub)

}

calculate_dic <- function(sigma_noise, inducing, x_train, y_train, x_new,
                          kernel_fun) {

  k_uf <- kernel_fun(inducing, x_train)
  k_fu <- t(k_uf) # At least I think that's right...?
  k_uu <- kernel_fun(inducing, inducing)

  big_sigma <- cholesky_inverse(sigma_noise^(-2) * k_uf %*% k_fu + k_uu)

  k_star_u <- kernel_fun(x_new, inducing)

  k_star_u_times_sigma <- k_star_u %*% big_sigma 

  dic_mean <- sigma_noise^(-2) * k_star_u_times_sigma %*% k_uf %*% y_train
  dic_cov <- k_star_u %*% big_sigma %*% t(k_star_u)

  return(list('mean' = dic_mean, 'cov' = dic_cov))

}

calculate_dtc <- function(sigma_noise, inducing, x_train, y_train, x_new,
                          kernel_fun) {

  # TODO: This is code duplication with DIC. Maybe clean up.
  k_uf <- kernel_fun(inducing, x_train)
  k_fu <- t(k_uf)
  k_uu <- kernel_fun(inducing, inducing) + diag(nrow(inducing)) * 1e-6

  big_sigma <- cholesky_inverse(sigma_noise^(-2) * k_uf %*% k_fu + k_uu)
  k_star_u <- kernel_fun(x_new, inducing)
  k_star_u_times_sigma <- k_star_u %*% big_sigma 

  # Same as DIC
  dtc_mean <- sigma_noise^(-2) * k_star_u_times_sigma %*% k_uf %*% y_train

  q_star_star <- k_star_u %*% cholesky_inverse(k_uu) %*% t(k_star_u)

  k_star_star <- kernel_fun(x_new, x_new)
  dtc_cov <- k_star_star - q_star_star + k_star_u %*% big_sigma %*% t(k_star_u)

  return(list('mean' = dtc_mean, 'cov' = dtc_cov))

}

calculate_fitc <- function(sigma_noise, inducing, x_train, y_train, x_new,
                           kernel_fun, diag_kernel_fun) {

  k_uf <- kernel_fun(inducing, x_train)
  k_fu <- t(k_uf)
  k_uu <- kernel_fun(inducing, inducing)
  k_uu_inv <- cholesky_inverse(k_uu + diag(nrow(inducing)) * 1e-6)

  q_ff <- k_fu %*% k_uu_inv %*% k_uf

  diagonal_k_ff <- diag_kernel_fun(x_train, x_train)

  big_lambda <- -diag(q_ff) + diagonal_k_ff + sigma_noise^2
  big_lambda_inv <- 1. / big_lambda
  big_lambda_inv <- diag(length(big_lambda_inv)) * big_lambda_inv

  big_sigma_pre_inv <- k_uu + k_uf %*% big_lambda_inv %*% k_fu 
  big_sigma_pre_inv <- big_sigma_pre_inv + diag(nrow(big_sigma_pre_inv)) * 1e-6

  big_sigma <- cholesky_inverse(big_sigma_pre_inv)

  k_star_u <- kernel_fun(x_new, inducing)
  k_star_u_times_sigma <- k_star_u %*% big_sigma 

  k_star_star <- kernel_fun(x_new, x_new)
  q_star_star <- k_star_u %*% k_uu_inv %*% t(k_star_u)

  fitc_mean <- k_star_u_times_sigma %*% k_uf %*% big_lambda_inv %*% y_train
  fitc_cov <- k_star_star - q_star_star + k_star_u_times_sigma %*% t(k_star_u)

  return(list('mean' = fitc_mean, 'cov' = fitc_cov))
}

calculate_fitc_naive <- function(sigma_noise, inducing, x_train, y_train,
                                 x_new, kernel_fun, diag_kernel_fun) {
  # A naive, non-optimised implementation
  k_star_u <- kernel_fun(x_new, inducing)
  k_uu <- kernel_fun(inducing, inducing)
  k_uu_inv <- cholesky_inverse(k_uu + (diag(nrow(k_uu)) * 1e-6))
  k_fu <- kernel_fun(x_train, inducing)

  q_star_f <- k_star_u %*% k_uu_inv %*% t(k_fu)

  q_ff <- k_fu %*% k_uu_inv %*% t(k_fu)

  big_lambda <- diag(kernel_fun(x_train, x_train) - q_ff) + sigma_noise^2

  summed_inv <- cholesky_inverse(q_ff + diag(nrow(q_ff)) * big_lambda +
                                 (diag(nrow(q_ff)) * 1e-6))

  fitc_mean <- q_star_f %*% summed_inv %*% y_train
  k_star_star <- kernel_fun(x_new, x_new)
  fitc_cov <- k_star_star - (q_star_f %*% summed_inv %*% t(q_star_f))

  return(list('mean' = fitc_mean, 'cov' = fitc_cov))

}

calculate_fic_naive <- function(sigma_noise, inducing, x_train, y_train,
                                x_new, kernel_fun, diag_kernel_fun) {
  # A naive version of FIC which just does exact GP inference with the FIC
  # kernel

  new_kernel <- function(v1, v2) {
    # Calculate the SoR kernel
    sor_kernel <- calculate_q(kernel_fun, inducing, v1, v2)
    # Adapt it with the diagonal
    diag_kernel <- diag_kernel_fun(v1, v2)
    diag(sor_kernel) <- diag_kernel
    return(sor_kernel)
  }

  return(predict_points(x_train, x_new, sigma_noise, y_train, new_kernel))

}

calculate_fic <- function(sigma_noise, inducing, x_train, y_train, x_new,
                          kernel_fun, diag_kernel_fun) {

  # Calculate f, the diagonal correction matrix
  diagonal_kff <- diag_kernel_fun(x_train, x_train)

  # We may be able to optimise this as we only care about the diagonal:
  q_ff <- calculate_q(kernel_fun, inducing, x_train, x_train)
  f <- (diagonal_kff - diag(q_ff) + sigma_noise^2)
  f_inv <- 1. / f * diag(length(f)) # As a diagonal matrix

  # Calculate the main inverse using the woodbury version
  k_uu <- kernel_fun(inducing, inducing)
  k_fu <- kernel_fun(x_train, inducing)

  to_invert <- k_uu + t(k_fu) %*% f_inv %*% k_fu
  main_inv <- cholesky_inverse(to_invert + diag(nrow(to_invert)) * tol)
  wood <- f_inv - f_inv %*% k_fu %*% main_inv %*% t(k_fu) %*% f_inv

  fic_kernel <- function(v1, v2) {
    new_kernel <- calculate_q(kernel_fun, inducing, v1, v2)
    correction <- diag_kernel_fun(v1, v2)
    diag(new_kernel) <- correction
    return(new_kernel)
  }

  k_star_x <- fic_kernel(x_new, x_train)
  k_star_star <- fic_kernel(x_new, x_new)

  fic_mean <- k_star_x %*% wood %*% y_train
  fic_cov <- k_star_star - k_star_x %*% wood %*% t(k_star_x)

  return(list('mean' = fic_mean, 'cov' = fic_cov))

}
