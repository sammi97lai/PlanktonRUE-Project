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
  
  
  # if existing smooth term
  # if(length(Pred_brms$dpars$mu$sm)>0){
  #   # fixed part
  #   beta = cbind(beta, Pred_brms$dpars$mu$sm$fe$bs)
  #   colnames(beta)
  #   
  #   # smooth part
  #   smooth = Pred_brms$dpars$mu$sm$re$sDOY$s[[1]] # the list has only one element
  #   colnames(smooth)
  # }
  
  # random part
  lambda = do.call(cbind, Pred_brms$dpars$mu$re$r)
  colnames(lambda) # 
  
  
  # 2-  we prepare the data matrix with all environemntal and random covariates values
  ## envrionmental variable matrix
  
  X_Beta <- as.matrix(Pred_brms$dpars$mu$fe$X)
  colnames(X_Beta)[1] =  paste('(Intercept)')
  
  # if existing smooth term
  # if(length(Pred_brms$dpars$mu$sm)>0){
  #   # fixed part
  #   X_Beta = as.matrix(cbind(X_Beta, Pred_brms$dpars$mu$sm$fe$Xs))
  #   colnames(X_Beta)
  #   
  #   # smooth part
  #   S_mat = Pred_brms$dpars$mu$sm$re$sDOY$Zs[[1]] # the list has only one element
  #   colnames(S_mat) = colnames(smooth)
  # }
  
  
  ## random variable matrix: 
  Eta <- do.call(cbind, Pred_brms$dpars$mu$re$Z)
  dim(Eta)
  colnames(Eta) = colnames(do.call(cbind, Pred_brms$dpars$mu$re$r))
  
  # 3- assemble them
  X <- cbind(X_Beta, 
             # S_mat, 
             Eta) ## combine environmental part and random effects
  
  brms_standat <- standata(mod)
  y_brms <- brms_standat$Y ## extract y
  
  params <- as.matrix(cbind(beta, 
                            # smooth, 
                            lambda)) #converts to matrix
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
construct_lin_term_mat <- function(params, X, y, residuals = TRUE) {
  
  ## first we need to differentiate if it a stan or brms model
  if(isTRUE(any(grepl("^r_", colnames(X))))){
    ## brms model
    
    ## Linear terms for the growth rate and density dependence
    X_cov <- X[, !grepl('Intercept', colnames(X), fixed = TRUE) & !grepl('[', colnames(X), fixed = TRUE)] ## lower case here!
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
    
    # Same for smoothing if they exist
    # if(any(grep('^s_', names_X))){
    #   names_smooth = names_X[grep("^S_", names_X, invert = FALSE, ignore.case = TRUE)] # upper case
    #   
    #   As_iid = list()
    #   
    #   for (sm in names_smooth) {
    #     Xs_iid <- X[, grep(sm, colnames(X), fixed = TRUE)]
    #     As_iid[[sm]] <- Xs_iid %*% params[colnames(Xs_iid)]
    #     colnames(As_iid[[sm]]) <- c(gsub(".*s_", "", sm))  # lower case
    #     colnames(As_iid[[sm]]) <- c(gsub("_.", "", colnames(As_iid[[sm]])))
    #   }
    #   
    # }
  
    
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
  
  # Same for smoothing if they exist
  # if(any(grep('^s_', names_X))){
  #   As_iid_all = do.call(cbind, As_iid)
  #   
  #   A <- cbind(A_cov, As_iid_all, A_iid_all)
  # }
  
  
  ## Linear term for residuals (omiting the intercept above does not affect variance)
  if(isTRUE(residuals)){
    A <- cbind(A, residuals = y - rowSums(A))
    
  }else{
    A <- A
    
  }
  
  return(as.matrix(A))
}


## Covariance matrix
#calculate covariance matrix of linear terms
#measures relationship between linear terms
lin_term_cov <- function(params, X, y, residuals) {
  A <- construct_lin_term_mat(params, X, y, residuals)
  A <- A[,linear_terms]
  K <- cov(A)
  dimnames(K) <- list(colnames(A), colnames(A))
  return(K)
}


conditional_lin_term_cov <- function(params, G, g, X, y, residuals) {
  A <- construct_lin_term_mat(params, X, y, residuals)
  
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

