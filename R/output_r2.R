fitting_lasso_return_r2 <- function(y, X, alpha = 1){

  require(glmnet)
  cv.lasso <- cv.glmnet(X, y, intercept = FALSE, alpha = alpha)
  r2 <- cv.lasso$glmnet.fit$dev.ratio[which(cv.lasso$glmnet.fit$lambda == cv.lasso$lambda.min)]
  coeff <- as.vector(coef(cv.lasso, s = cv.lasso$lambda.min))
  selected <- which(coeff!=0)
  res <- list(coeff = coeff[selected], selected = selected - 1, r2 = r2)
  return(res)
}

Output_r_2 <- function(gene.expression, percentage.cutoff = 0.1, num = 5000, ImputeAll = FALSE){

  xx <- gene.expression # p*n
  p <- nrow(xx)
  n <- ncol(xx)
  zero.rate <- apply(xx, 1, function(x){length(x[x == 0])})/n
  flag <- zero.rate <= percentage.cutoff

  logxx <- apply(xx, 2, function(y){log(y + 0.1)})
  data <- logxx[round(runif(num)*p), ]
  zero.matrix <- xx != 0
  zero.matrix <- apply(zero.matrix, 2, as.numeric)
  t_logxx <- t(logxx)
  t_zero.matrix <- t(zero.matrix)

  imputed <- logxx
  weights.list <- list(NULL)
  outlier.list <- NULL
  r_square_table <- list(NULL)

  for(j in 1:n){

    remain <- data[, -j]
    res <- fitting_lasso_return_r2(data[, j], remain)
    coeff <- res$coeff
    selected <- res$selected

    if(length(selected) < 3){

      weights.list[[j]] <- cbind(c(1:n)[-j][selected], rep(1, length(selected)))
      outlier.list <- c(outlier.list, j)

    } else {
      prior.weight <- calculate_weights(logxx[flag, j], logxx[flag, -j][, selected])

      kk <- summary(lm(logxx[flag, j]~logxx[flag, -j][, selected]))

      sub.selected <- selected[prior.weight >= 10^-4]
      sub.prior.weight <- prior.weight[prior.weight >= 10^-4]
      sub.prior.weight <- sub.prior.weight/sum(sub.prior.weight)

      weights.list[[j]] <- cbind(c(1:n)[-j][sub.selected], sub.prior.weight)

      Ymat <- t_logxx[-j, ][sub.selected, ]
      Yflagmat <- t_zero.matrix[-j, ][sub.selected, ]

      kk2 <- summary(lm(logxx[flag, j]~logxx[flag, -j][, sub.selected]))

      r_square_table[[j]] <- c(res$r2, kk$r.squared, kk$adj.r.squared, kk2$r.squared, kk2$adj.r.squared, length(selected), length(sub.selected))

      if(ImputeAll == TRUE){
        tt <- reweighting_sum_C(Ymat, Yflagmat, logxx[, j], zero.matrix[, j], sub.prior.weight, TRUE)
      } else {
        tt <- reweighting_sum_C(Ymat, Yflagmat, logxx[, j], zero.matrix[, j], sub.prior.weight, FALSE)
      }

      imputed[, j] <- tt
    }

  }

  res <- list(weights.list = weights.list, imputed = imputed, outlier.list = outlier.list, r2.table = r_square_table)
  return(res)
}



Output_r_2_gene_regression <- function(gene.expression, gene.num = 1000){

  xx <- gene.expression # p*n
  p <- nrow(xx)
  n <- ncol(xx)

  logxx <- apply(xx, 2, function(y){log(y + 0.1)})
  data <- t(logxx)
  zero.matrix <- xx != 0
  zero.matrix <- apply(zero.matrix, 2, as.numeric)
  t_logxx <- t(logxx)
  t_zero.matrix <- t(zero.matrix)

  imputed <- logxx
  weights.list <- list(NULL)
  outlier.list <- NULL
  r_square_table <- list(NULL)

  selected.nums <- sample(1:p, gene.num)

  for(j in 1:selected.nums){

    remain <- data[, -j]
    res <- try(fitting_lasso_return_r2(data[, j], remain))
    if(class(res) != "try-error"){
      coeff <- res$coeff
      selected <- res$selected

      if(length(selected) >= 1){
        kk <- summary(lm(data[, j]~remain[, selected]))
        r_square_table[[j]] <- c(res$r2, kk$r.squared, kk$adj.r.squared, length(selected))
      }
    }
  }

  res <- list(r2.table = r_square_table)
  return(res)
}



Output_lasso_neighbors <- function(gene.expression, percentage.cutoff = 0.1, num = 5000, ImputeAll = FALSE){

  xx <- gene.expression # p*n
  p <- nrow(xx)
  n <- ncol(xx)
  zero.rate <- apply(xx, 1, function(x){length(x[x == 0])})/n
  flag <- zero.rate <= percentage.cutoff

  logxx <- apply(xx, 2, function(y){log(y + 0.1)})
  data <- logxx[round(runif(num)*p), ]
  zero.matrix <- xx != 0
  zero.matrix <- apply(zero.matrix, 2, as.numeric)
  t_logxx <- t(logxx)
  t_zero.matrix <- t(zero.matrix)

  imputed <- logxx
  selection.list <- list(NULL)


  for(j in 1:n){

    remain <- data[, -j]
    res <- fitting_lasso_return_r2(data[, j], remain)
    coeff <- res$coeff
    selection.list[[j]] <- res$selected

  }

  res <- list(selection.list)
  return(res)
}


