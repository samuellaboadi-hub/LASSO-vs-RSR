## LASSO vs RSR on REAL data (mtcars)
## Outcome: mpg
## LASSO: main effects only
## RSR:   main + (quadratics + interactions) with rank weights


set.seed(123)

dat <- mtcars
y <- dat$mpg
X_main <- as.matrix(dat[, setdiff(names(dat), "mpg")])   # main effects only

n <- nrow(X_main)
p0 <- ncol(X_main)


## Helpers

soft_thresh <- function(z, lam) sign(z) * pmax(abs(z) - lam, 0)


## LASSO (Coordinate Descent) + CV

lasso_cd <- function(X, y, lam, max_iter = 2000, tol = 1e-6, standardize = TRUE) {
  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  
  # center y
  y_mean <- mean(y)
  y_c <- y - y_mean
  
  # standardize X
  x_means <- rep(0, p)
  x_sds   <- rep(1, p)
  Xs <- X
  
  if (standardize) {
    x_means <- colMeans(Xs)
    Xs <- sweep(Xs, 2, x_means, "-")
    x_sds <- apply(Xs, 2, sd)
    x_sds[x_sds == 0] <- 1
    Xs <- sweep(Xs, 2, x_sds, "/")
  }
  
  beta <- rep(0, p)
  f <- as.numeric(Xs %*% beta)
  
  for (iter in 1:max_iter) {
    beta_old <- beta
    
    for (j in 1:p) {
      rj <- y_c - (f - Xs[, j] * beta[j])
      zj <- sum(Xs[, j] * rj) / n
      aj <- sum(Xs[, j]^2) / n
      
      bj_new <- soft_thresh(zj, lam) / aj
      
      if (bj_new != beta[j]) {
        f <- f + Xs[, j] * (bj_new - beta[j])
        beta[j] <- bj_new
      }
    }
    
    if (max(abs(beta - beta_old)) < tol) break
  }
  
  # unstandardize
  beta_unscaled <- beta
  if (standardize) beta_unscaled <- beta / x_sds
  
  intercept <- y_mean - sum(beta_unscaled * x_means)
  
  list(intercept = intercept, beta = beta_unscaled, iter = iter)
}

predict_lasso <- function(fit, Xnew) {
  as.numeric(fit$intercept + as.matrix(Xnew) %*% fit$beta)
}

kfold_cv_lasso <- function(X, y, lambdas, K = 10, max_iter = 2000, tol = 1e-6) {
  n <- nrow(X)
  folds <- sample(rep(1:K, length.out = n))
  cv_mse <- numeric(length(lambdas))
  
  for (li in seq_along(lambdas)) {
    lam <- lambdas[li]
    mse_k <- numeric(K)
    
    for (k in 1:K) {
      te <- which(folds == k)
      tr <- setdiff(1:n, te)
      
      fit <- lasso_cd(X[tr, , drop = FALSE], y[tr], lam,
                      max_iter = max_iter, tol = tol, standardize = TRUE)
      yhat <- predict_lasso(fit, X[te, , drop = FALSE])
      mse_k[k] <- mean((y[te] - yhat)^2)
    }
    cv_mse[li] <- mean(mse_k)
  }
  
  best <- which.min(cv_mse)
  list(lambdas = lambdas, cv_mse = cv_mse, lambda_min = lambdas[best])
}

# lambda grid for LASSO
y_c <- y - mean(y)
Xs_tmp <- scale(X_main, center = TRUE, scale = TRUE)
lambda_max_lasso <- max(abs(colSums(Xs_tmp * y_c) / n))
lambdas_lasso <- exp(seq(log(lambda_max_lasso), log(lambda_max_lasso * 0.01), length.out = 60))

cv_lasso  <- kfold_cv_lasso(X_main, y, lambdas_lasso, K = 10, max_iter = 2000, tol = 1e-6)
fit_lasso <- lasso_cd(X_main, y, cv_lasso$lambda_min, max_iter = 4000, tol = 1e-8)

lasso_table <- data.frame(variable = colnames(X_main), beta = fit_lasso$beta)
lasso_table <- subset(lasso_table, abs(beta) > 1e-8)
lasso_table <- lasso_table[order(-abs(lasso_table$beta)), ]

cat("\n================ LASSO ================\n")
cat("Chosen lambda:", cv_lasso$lambda_min, "\n")
cat("Intercept:", fit_lasso$intercept, "\n")
cat("Selected variables:\n")
print(lasso_table, row.names = FALSE)


## Build ranked features for RSR
## Group 1: main effects
## Group 2: quadratics + interactions

main_names <- colnames(X_main)

X_quad <- X_main^2
quad_names <- paste0(main_names, "^2")

pair_idx <- combn(p0, 2)
X_int <- matrix(NA_real_, nrow = n, ncol = ncol(pair_idx))
int_names <- character(ncol(pair_idx))
for (m in 1:ncol(pair_idx)) {
  j <- pair_idx[1, m]
  k <- pair_idx[2, m]
  X_int[, m] <- X_main[, j] * X_main[, k]
  int_names[m] <- paste0(main_names[j], ":", main_names[k])
}

X_rsr <- cbind(X_main, X_quad, X_int)
feat_names <- c(main_names, quad_names, int_names)

group <- c(rep(1, p0), rep(2, ncol(X_quad) + ncol(X_int)))
stopifnot(length(group) == ncol(X_rsr))

## RSR (Coordinate Descent) + CV 

rsr_cd <- function(X, y, group, lambda, gamma = 1.0,
                   max_iter = 2000, tol = 1e-6, standardize = TRUE) {
  
  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  
  # center y
  y_mean <- mean(y)
  y_c <- y - y_mean
  
  # standardize X
  x_means <- rep(0, p)
  x_sds <- rep(1, p)
  Xs <- X
  if (standardize) {
    x_means <- colMeans(Xs)
    Xs <- sweep(Xs, 2, x_means, "-")
    x_sds <- apply(Xs, 2, sd)
    x_sds[x_sds == 0] <- 1
    Xs <- sweep(Xs, 2, x_sds, "/")
  }
  
  # group weights: w_k = p_k^(1 - 2*gamma)
  groups <- sort(unique(group))
  p_k <- sapply(groups, function(g) sum(group == g))
  w_k <- p_k^(1 - 2 * gamma)
  names(w_k) <- groups
  w_j <- w_k[as.character(group)]
  
  beta <- rep(0, p)
  f <- as.numeric(Xs %*% beta)
  
  for (iter in 1:max_iter) {
    beta_old <- beta
    
    for (j in 1:p) {
      rj <- y_c - (f - Xs[, j] * beta[j])
      zj <- sum(Xs[, j] * rj) / n
      aj <- sum(Xs[, j]^2) / n
      
      bj_new <- soft_thresh(zj, lambda * w_j[j]) / aj
      
      if (bj_new != beta[j]) {
        f <- f + Xs[, j] * (bj_new - beta[j])
        beta[j] <- bj_new
      }
    }
    
    if (max(abs(beta - beta_old)) < tol) break
  }
  
  beta_unscaled <- beta
  if (standardize) beta_unscaled <- beta / x_sds
  intercept <- y_mean - sum(beta_unscaled * x_means)
  
  list(intercept = intercept, beta = beta_unscaled, iter = iter, w_k = w_k)
}

predict_rsr <- function(fit, Xnew) {
  as.numeric(fit$intercept + as.matrix(Xnew) %*% fit$beta)
}

kfold_cv_rsr <- function(X, y, group, lambdas, K = 5, gamma = 1.0,
                         max_iter = 500, tol = 1e-5) {
  n <- nrow(X)
  folds <- sample(rep(1:K, length.out = n))
  cv_mse <- numeric(length(lambdas))
  
  for (li in seq_along(lambdas)) {
    lam <- lambdas[li]
    mse_k <- numeric(K)
    
    for (k in 1:K) {
      te <- which(folds == k)
      tr <- setdiff(1:n, te)
      
      fit <- rsr_cd(X[tr, , drop = FALSE], y[tr], group, lam, gamma = gamma,
                    max_iter = max_iter, tol = tol, standardize = TRUE)
      yhat <- predict_rsr(fit, X[te, , drop = FALSE])
      mse_k[k] <- mean((y[te] - yhat)^2)
    }
    cv_mse[li] <- mean(mse_k)
  }
  
  best <- which.min(cv_mse)
  list(lambdas = lambdas, cv_mse = cv_mse, lambda_min = lambdas[best])
}

# lambda grid for RSR (20 values)
Xs_tmp2 <- scale(X_rsr, center = TRUE, scale = TRUE)
lambda_max_rsr <- max(abs(colSums(Xs_tmp2 * y_c) / n))
lambdas_rsr <- exp(seq(log(lambda_max_rsr), log(lambda_max_rsr * 0.01), length.out = 20))

gamma <- 1.0

cv_rsr  <- kfold_cv_rsr(X_rsr, y, group, lambdas_rsr, K = 5, gamma = gamma,
                        max_iter = 500, tol = 1e-5)
fit_rsr <- rsr_cd(X_rsr, y, group, cv_rsr$lambda_min, gamma = gamma,
                  max_iter = 2000, tol = 1e-8)

rsr_table <- data.frame(feature = feat_names, group = group, beta = fit_rsr$beta)
rsr_table <- subset(rsr_table, abs(beta) > 1e-8)
rsr_table <- rsr_table[order(-abs(rsr_table$beta)), ]

cat("\n================ RSR ================\n")
cat("Chosen lambda:", cv_rsr$lambda_min, "\n")
cat("Intercept:", fit_rsr$intercept, "\n")
cat("Group weights w_k:\n")
print(fit_rsr$w_k)
cat("Selected features:\n")
print(rsr_table, row.names = FALSE)


## Comparison plots

par(mfrow = c(1, 2))
plot(log(cv_lasso$lambdas), cv_lasso$cv_mse, type = "b",
     xlab = "log(lambda)", ylab = "CV MSE", main = "LASSO")
abline(v = log(cv_lasso$lambda_min), lty = 2)

plot(log(cv_rsr$lambdas), cv_rsr$cv_mse, type = "b",
     xlab = "log(lambda)", ylab = "CV MSE", main = "RSR")
abline(v = log(cv_rsr$lambda_min), lty = 2)
par(mfrow = c(1, 1))

barplot(c(LASSO = min(cv_lasso$cv_mse), RSR = min(cv_rsr$cv_mse)),
        ylab = "Min CV MSE", main = "Predictive Performance")

barplot(c(LASSO = sum(abs(fit_lasso$beta) > 1e-8),
          RSR   = sum(abs(fit_rsr$beta) > 1e-8)),
        ylab = "# Non-zero coefficients", main = "Model Sparsity")



par(mfrow = c(1, 2))
yhat_lasso <- predict_lasso(fit_lasso, X_main)

yhat_rsr <- predict_rsr(fit_rsr, X_rsr)
plot(y, yhat_lasso,
     xlab = "Observed mpg",
     ylab = "Predicted mpg",
     main = "LASSO",
     pch = 19, col = "blue")
abline(0, 1, lty = 2, lwd = 2)

plot(y, yhat_rsr,
     xlab = "Observed mpg",
     ylab = "Predicted mpg",
     main = "RSR",
     pch = 19, col = "darkgreen")
abline(0, 1, lty = 2, lwd = 2)

par(mfrow = c(1, 1))

library(sparseR)

# Trying with the sparseR package
set.seed(123)

rsr_sparseR <- sparseR(
  mpg ~ .,
  data = mtcars,
  k = 1,          # include pairwise interactions
  family = "gaussian",
  seed = 123
)

coef_rsr <- coef(rsr_sparseR, at = "cvmin")

# Keep non-zero coefficients
selected_rsr <- coef_rsr[coef_rsr != 0]
selected_rsr

