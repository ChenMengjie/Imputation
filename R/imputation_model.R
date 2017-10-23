
Imputation1_cpp <- function(gene.expression, percentage.cutoff = 0.1, num = 10000,  minbool = FALSE){

  xx <- gene.expression # p*n
  p <- nrow(xx)
  n <- ncol(xx)
  if(p < num) num <- round(0.8*p)
  zero.rate <- apply(xx, 1, function(x){length(x[x == 0])})/n
  flag <-  zero.rate <= percentage.cutoff

  logxx <- apply(xx, 2, function(y){log(y + 0.1)})
  data <- logxx[round(runif(num)*p), ]
  zero.matrix <- xx != 0
  zero.matrix <- apply(zero.matrix, 2, as.numeric)
  selected_logxx <- logxx[flag, ]
  res_imp <- imputation_by_samples(data, selected_logxx, logxx, zero.matrix, n, p, minbool)
  outlier_flag <- apply(res_imp$sample_weights, 1, function(x){any(x==-1)})
  outliers <- c(1:n)[outlier_flag]

  nopredict <- logxx
  nopredict[gene.expression==0] <- res_imp$imputed[gene.expression==0]

  res <- list(predicted = res_imp$imputed, imputed = nopredict,
              sample_weights = res_imp$sample_weights, outliers=outliers)
  return(res)
}

Imputation3 <- function(gene.expression, percentage.cutoff = 0.1, num = 5000, percentage.samples = 0.8){

  xx <- gene.expression # p*n
  p <- nrow(xx)
  n <- ncol(xx)
  if(p < num) num <- round(0.8*p)
  zero.rate <- apply(xx, 1, function(x){length(x[x == 0])})/n
  flag <-  zero.rate <= percentage.cutoff

  logxx <- apply(xx, 2, function(y){log(y + 0.1)})
  data <- logxx[round(runif(num)*p), ]
  data.copy2 <- t(logxx[flag, round(runif(ceiling(percentage.samples*n))*n)])

  zero.matrix <- xx != 0
  zero.matrix <- apply(zero.matrix, 2, as.numeric)
  t_logxx <- t(logxx)
  t_zero.matrix <- t(zero.matrix)

  imputed <- logxx

  sample.weights.list <- matrix(0, nrow = n, ncol = n)
  gene.weights.list <- list(NULL)
  outlier.list <- NULL

  for(j in 1:n){

    remain <- data[, -j]
    res <- fitting_lasso(data[, j], remain)
    coeff <- res$coeff
    selected <- res$selected

    if(length(selected) < 3){

      sample.weights.list[j, c(1:n)[-j][selected]] <- rep(1, length(selected))
      outlier.list <- c(outlier.list, j)

    } else {

      prior.weight <- calculate_weights(logxx[flag, j], logxx[flag, -j][, selected])
      sub.selected <- selected[prior.weight >= 10^-4]
      sub.prior.weight <- prior.weight[prior.weight >= 10^-4]
      sub.prior.weight <- sub.prior.weight/sum(sub.prior.weight)
      sample.weights.list[j, c(1:n)[-j][sub.selected]] <-  sub.prior.weight

      Ymat <- t_logxx[-j, ][sub.selected, ]
      Yflagmat <- t_zero.matrix[-j, ][sub.selected, ]

      imputed[, j] <- reweighting_sum_C(Ymat, Yflagmat, logxx[, j], zero.matrix[, j], sub.prior.weight, TRUE)

    }

  }

  imputed.gene <- imputed
  which_flag <- which(flag)

  for(i in 1:ncol(data.copy2)){

      remain <- data.copy2[, -i]
      res <- cor(data.copy2[, i], remain)
      selected <- which(abs(res) > 0.8)
      len <- length(selected)
      if(len == 0 ){

        gene.weights.list[[which_flag[i]]] <- NULL

      } else if(len == 1){
        gene.weights.list[[which_flag[i]]] <-  cbind(selected, 1)
        imputed.gene[i, ] <- imputed[selected, ]
      } else {
        prior.weight <- calculate_weights(data.copy2[, i], data.copy2[, -i][, selected])
        sub.selected <- selected[prior.weight >= 10^-4]
        if(length(sub.selected) < 2){
          sub.prior.weight <- rep(1/len, len)
          aa <- which_flag[-i][selected]

        } else {
          sub.prior.weight <- prior.weight[prior.weight >= 10^-4]
          sub.prior.weight <- sub.prior.weight/sum(sub.prior.weight)
          aa <- which_flag[-i][sub.selected]

        }
        gene.weights.list[[which_flag[i]]] <-  cbind(aa, sub.prior.weight)

        imputed.from.other.genes <- apply(imputed[aa, ], 2, function(dd){sum(dd*sub.prior.weight)})
        imputed.gene[i, ] <- imputed.from.other.genes
        #imputed.gene[i, ] <- 0.5*imputed[i, ] + 0.5*imputed.from.other.genes

      }

    }


  res <- list(sample.weights.list = sample.weights.list, gene.weights.list = gene.weights.list,
              outlier.list = outlier.list, imputed = imputed, imputed.gene = imputed.gene)
  return(res)

}


Imputation3_cpp <- function(gene.expression, percentage.cutoff = 0.1, num = 10000, percentage.samples = 0.8, minbool = FALSE){

  xx <- gene.expression # p*n
  p <- nrow(xx)
  n <- ncol(xx)
  if(p < num) num <- round(0.8*p)
  zero.rate <- apply(xx, 1, function(x){length(x[x == 0])})/n
  flag <- zero.rate <= percentage.cutoff
  which_flag <- which(flag)
  logxx <- apply(xx, 2, function(y){log(y + 0.1)})
  data <- logxx[round(runif(num)*p), ]
  sample.flag <- round(runif(round(percentage.samples*n))*n)
  data.copy2 <- t(logxx[flag, sample.flag])

  zero.matrix <- xx != 0
  zero.matrix <- apply(zero.matrix, 2, as.numeric)
  selected_logxx <- logxx[flag, ]
  res_imp <- imputation_by_samples(data, selected_logxx, logxx, zero.matrix, n, p, minbool)
  outlier_flag <- apply(res_imp$sample_weights, 1, function(x){any(x == -1)})
  outliers <- c(1:n)[outlier_flag]
  imputed <- res_imp$imputed
  imputed.gene <- imputation_by_genes(imputed, data.copy2, which_flag-1)

  pair.flag <- imputed.gene$gene_weights_1!=0
  gene_weights <- cbind(imputed.gene$gene_weights_1[pair.flag], imputed.gene$gene_weights_2[pair.flag],
  imputed.gene$gene_weights_3[pair.flag])

  nopredict <- logxx
  nopredict[gene.expression==0] <- res_imp$imputed[gene.expression==0]

  nopredict2 <- logxx
  nopredict2[gene.expression==0] <- imputed.gene$imputed_gene[gene.expression==0]

  res <- list(sample_weights = res_imp$sample_weights, predicted = res_imp$imputed,
              imputed = nopredict, outliers = outliers, predicted.by.gene = imputed.gene$imputed_gene,
              imputed.by.gene = nopredict2, gene_weights = gene_weights)
  return(res)

}



reweighting_with_bulk <- function(prior, meanY, sdY = NULL, Yflag){
  k <- length(Y)
  zero_rate <- length(Yflag[Yflag == 0])/k
  if(zero_rate < 0.9){
    if(is.null(sdY)) sdY <- 1
    pN <- dnorm(0, meanY, sdY)
    non.dropout <- (1-zero_rate)*pN/(zero_rate + (1-zero_rate)*pN)
    zero_vec <- rep(1, k)
    zero_vec[Yflag == 0] <- non.dropout
    xx <- prior*zero_vec
    return(xx/sum(xx))
  } else {
    return(prior/sum(prior))
  }
}


Imputation4 <- function(gene.expression, bulk, sample.id.in.bulk, percentage.cutoff = 0.1, num = 10000,  minbool = FALSE){

  xx <- gene.expression # p*n
  p <- nrow(xx)
  n <- ncol(xx)
  if(p < num) num <- round(0.8*p)
  zero.rate <- apply(xx, 1, function(x){length(x[x == 0])})/n
  flag <-  zero.rate <= percentage.cutoff

  logxx <- apply(xx, 2, function(y){log(y + 0.1)})
  data <- logxx[round(runif(num)*p), ]
  zero.matrix <- xx != 0
  zero.matrix <- apply(zero.matrix, 2, as.numeric)
  t_logxx <- t(logxx)
  t_zero.matrix <- t(zero.matrix)

  imputed <- logxx
  weights.list <- list(NULL)
  outlier.list <- NULL

  for(j in 1:n){

    remain <- data[, -j]
    res <- fitting_lasso(data[, j], remain)
    coeff <- res$coeff
    selected <- res$selected

    if(length(selected) < 3){

      weights.list[[j]] <- cbind(c(1:n)[-j][selected], rep(1, length(selected)))
      outlier.list <- c(outlier.list, j)

    } else {
      prior.weight <- calculate_weights(logxx[flag, j], logxx[flag, -j][, selected])

      sub.selected <- selected[prior.weight >= 10^-4]
      sub.prior.weight <- prior.weight[prior.weight >= 10^-4]
      sub.prior.weight <- sub.prior.weight/sum(sub.prior.weight)

      weights.list[[j]] <- cbind(c(1:n)[-j][sub.selected], sub.prior.weight)

      Ymat <- t_logxx[-j, ][sub.selected, ]
      Yflagmat <- t_zero.matrix[-j, ][sub.selected, ]

      if(is.null(sample.id.in.bulk[[j]])){
        tt <- reweighting_sum_C(Ymat, Yflagmat, logxx[, j], zero.matrix[, j], sub.prior.weight, TRUE)
      } else {

        if(length(sample.id.in.bulk[[j]]) > 1){
          meanB <- apply(bulk[, sample.id.in.bulk[[j]]], 1, mean)
          sdB <- apply(bulk[, sample.id.in.bulk[[j]]], 1, sd)
        } else {
          meanB <- bulk[, sample.id.in.bulk[[j]]]
          sdB <- rep(1, p)
        }
        tt <- reweighting_with_bulk_C(Ymat, Yflagmat, meanY, sdY, prior_weight)

      }

      imputed[, j] <- tt
    }

  }

  res <- list(weights.list = weights.list, imputed = imputed, outlier.list = outlier.list)
  return(res)
}



