#
#       Adapt Variance partition to Brms
#
#.....................................................

## prepare element to feed to the variance partition functions adapted from Torsti

PrepbrmsVP = function(mod){
  
  ## the prepare function from brms give us access to all needed element
  Pred_brms = prepare_predictions(mod)
  
  # 1- we extract the fixed and random effect coefficients:
  # Fixed part 
  beta = Pred_brms$dpars$mu$fe$b
  colnames(beta)
  
  # random part
  lambda = do.call(cbind, Pred_brms$dpars$mu$re$r)
  colnames(lambda) # 
  
  
  # 2-  we prepare the data matrix with all environemntal and random covariates values
  ## envrionmental variable matrix
  
  X_Beta <- as.matrix(Pred_brms$dpars$mu$fe$X)
  colnames(X_Beta)[1] =  paste('(Intercept)')
  
  ## random variable matrix: 
  Eta <- do.call(cbind, Pred_brms$dpars$mu$re$Z)
  dim(Eta)
  colnames(Eta) = colnames(do.call(cbind, Pred_brms$dpars$mu$re$r))
  
  # 3- assemble them
  X <- cbind(X_Beta, Eta) ## combine environmental part and random effects
  
  brms_standat <- standata(mod)
  y_brms <- brms_standat$Y ## extract y
  
  params <- as.matrix(cbind(beta, lambda)) #converts to matrix
  colnames(params) = colnames(X)
  dimnames(params) = list('iterations' = 1:dim(params)[1], 
                          'parameters' = colnames(params))
  
  return(list(X = X, y = y_brms, params = params))
}


## VP functions adapted from Torsti .............................................

## Functions .........................

#builds matrix A of linear terms but combining effects from covariates, random effects (ID and year), and residuals
#represent specific contributions

## This functions now autmotically recognize the random effects!!
construct_lin_term_mat <- function(params, X, y) {
  
  ## first we need to differentiate if it a stan or brms model
  if(isTRUE(any(grepl("^r_", colnames(X))))){
    ## brms model
    
    ## Linear terms for the growth rate and density dependence
    X_cov <- X[, !grepl('Intercept', colnames(X), fixed = TRUE)]
    A_cov <- as.matrix(X_cov %*% Matrix::Diagonal(x = params[colnames(X_cov)]))
    colnames(A_cov) <- colnames(X_cov)
    
    # first we get the name of the random effect and how many they are
    names_X = unique(gsub("\\[.*","", colnames(X)))
    names_random = names_X[grep("^r_", names_X, invert = FALSE, ignore.case = TRUE)]
    
    ## Linear term for the ID random effect
    A_iid = list()
    
    for (rd in names_random) {
      X_iid <- X[, grep(rd, colnames(X), fixed = TRUE)]
      A_iid[[rd]] <- X_iid %*% params[colnames(X_iid)]
      colnames(A_iid[[rd]]) <- c(gsub(".*r_", "", rd))
    }
    
  }else{
    # stan model
    
    ## Linear terms for the growth rate and density dependence
    X_cov <- X[, !grepl('(Intercept)', colnames(X), fixed = TRUE)]
    A_cov <- as.matrix(X_cov %*% Matrix::Diagonal(x = params[colnames(X_cov)]))
    colnames(A_cov) <- colnames(X_cov)
    
    col_X = grep(":", colnames(X), value = TRUE)
    names_X = unique(gsub("\\:.*","", col_X))
    names_random = unique(gsub(".* ","", names_X))
    
    ## Linear term for the ID random effect
    A_iid = list()
    
    for (rd in names_random) {
      X_iid <- X[, grep(rd, colnames(X), fixed = TRUE)]
      A_iid[[rd]] <- X_iid %*% params[colnames(X_iid)]
      colnames(A_iid[[rd]]) <- c(gsub(".*  (.+) :.*", "\\1", rd))
    }
  }
  

  A_iid_all = do.call(cbind, A_iid)
  
  A <- cbind(A_cov, A_iid_all)
  
  ## Linear term for residuals (omiting the intercept above does not affect variance)
  #A <- cbind(A, residuals = y - rowSums(A))
  
  return(as.matrix(A))
}


## Covariance matrix
#calculate covariance matrix of linear terms
#measures relationship between linear terms
lin_term_cov <- function(params, X, y) {
  A <- construct_lin_term_mat(params, X, y)
  A <- A[,linear_terms]
  K <- cov(A)
  dimnames(K) <- list(colnames(A), colnames(A))
  return(K)
}


conditional_lin_term_cov <- function(params, G, g, X, y) {
  A <- construct_lin_term_mat(params, X, y)
  
  K <- t(sapply(G, function(x, A, g) { cov(A[g == x,])}, A, g))
  dim(K) <- c(length(G), rep(ncol(A), 2))
  dimnames(K) <- list(condition = G, colnames(A), colnames(A))
  
  return(K)
}

wit_btw_lin_term_cov <- function(params, H, X, y) {
  A <- construct_lin_term_mat(params, X, y)
  
  A_cm <- H %*% diag(1 / colSums(H)) %*% t(H) %*% A
  
  K_wit <- cov(A - A_cm)
  K_btw <- cov(A_cm)
  
  K <- abind::abind(K_wit, K_btw, rev.along = 0)
  dimnames(K) <- list(
    rows = colnames(A),
    cols = colnames(A),
    partition = c('within', 'between')
  )
  
  return(K)
}




vp_plot <- function(x, col = 'grey80', nm, main, las =1) {
  vp <- as.data.frame(aperm(x, c('iterations', 'variance')))
  colnames(vp) <- nm
  
  if (max(vp) <= 1 & min(vp) >= 0) {
    ylim <- 0:1
  } else {
    ylim <- NULL
  }
  
  vioplot(
    vp,
    col = col,
    las=las,
    drawRect = TRUE,
    pchMed = NA,
    ylim = ylim,
    main = main
  )
}


## VP functions adapted from Torsti .............................................

