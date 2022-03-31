cv_elasticnet <- function(x, y, alpha, lambda, nfolds = 5) {
  n <- nrow(x)
  ind <- sample(1:n)
  ni <- round(n/nfolds)
  group <- ceiling(seq_along(ind)/ni)
  group[group>nfolds] <- nfolds
  
  fold_list <- split(ind, factor(group))
  corrs <- c()
  
  for (i in 1:length(alpha)) {
    tmp <- c()
    for (j in 1:nfolds) {
      val_id <- fold_list[[j]]
      train_id <- setdiff(1:n, val_id)
      x_train <- x[train_id,]
      y_train <- y[train_id]
      x_val <- x[val_id, ]
      y_val <- y[val_id]
      glmnet_fit <- glmnet(x_train, y_train, alpha = alpha[i], lambda = lambda[i])
      glm_pred <- predict(glmnet_fit, x_val)[,1]
      tmp[j] <- cor(glm_pred, y_val)
    }
    tmp[is.na(tmp)] <- 0
    corrs[i] <- mean(tmp)
  }
  return(data.frame(alpha, lambda, corrs))
}
