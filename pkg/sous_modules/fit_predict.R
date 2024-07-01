library(randtoolbox)
library(lhs)
L_case = function(case_str,quanti_matrix,prior_matrix,F_matrix,kernel_str,Y_vec,reg,opt_reg,pow_w_10){
  if(pow_w_10){
    if(case_str == "normal"){
      if(opt_reg){
        L = function(vec){
          p = length(quanti_matrix[1,])
          q = length(prior_matrix[1,])
          A_matrix = matrix(vec[(p+1):(length(vec)-1)],nrow = q)
          w_vec = vec[1:p]
          r = vec[length(vec)]
          return(L_normal_cpp(quanti_matrix,prior_matrix,F_matrix,kernel_str,A_matrix,10**w_vec,Y_vec,r))
        }
      } else{
        L = function(vec){
          p = length(quanti_matrix[1,])
          q = length(prior_matrix[1,])
          A_matrix = matrix(vec[(p+1):length(vec)],nrow = q)
          w_vec = vec[1:p]
          return(L_normal_cpp(quanti_matrix,prior_matrix,F_matrix,kernel_str,A_matrix,10**w_vec,Y_vec,reg))
        }
      }
    } else if(case_str == "quanti"){
      if(opt_reg){
        L = function(vec){
          vec1 = vec[1:(length(vec)-1)]
          r = vec[length(vec)]
          return(L_quanti_cpp(quanti_matrix,F_matrix,kernel_str,10**vec1,Y_vec,r))
        }
      } else{
        L = function(vec){
          return(L_quanti_cpp(quanti_matrix,F_matrix,kernel_str,10**vec,Y_vec,reg))
        }
      }
    } else if(case_str == "quali"){
      if(opt_reg){
        L = function(vec){
          q = length(prior_matrix[1,])
          A_matrix = matrix(vec[1:(length(vec)-1)],nrow = q)
          r = vec[length(vec)]
          return(L_quali_cpp(prior_matrix,F_matrix,kernel_str,A_matrix,Y_vec,r))
        }
      } else{
        L = function(vec){
          q = length(prior_matrix[1,])
          A_matrix = matrix(vec,nrow = q)
          return(L_quali_cpp(prior_matrix,F_matrix,kernel_str,A_matrix,Y_vec,reg))
        }
      }
    } else{
      stop("erreur L_case")
    }
  } else{
    if(case_str == "normal"){
      if(opt_reg){
        L = function(vec){
          p = length(quanti_matrix[1,])
          q = length(prior_matrix[1,])
          A_matrix = matrix(vec[(p+1):(length(vec)-1)],nrow = q)
          w_vec = vec[1:p]
          r = vec[length(vec)]
          return(L_normal_cpp(quanti_matrix,prior_matrix,F_matrix,kernel_str,A_matrix,w_vec,Y_vec,r))
        }
      } else{
        L = function(vec){
          p = length(quanti_matrix[1,])
          q = length(prior_matrix[1,])
          A_matrix = matrix(vec[(p+1):length(vec)],nrow = q)
          w_vec = vec[1:p]
          return(L_normal_cpp(quanti_matrix,prior_matrix,F_matrix,kernel_str,A_matrix,w_vec,Y_vec,reg))
        }
      }
    } else if(case_str == "quanti"){
      if(opt_reg){
        L = function(vec){
          vec1 = vec[1:(length(vec)-1)]
          r = vec[length(vec)]
          return(L_quanti_cpp(quanti_matrix,F_matrix,kernel_str,vec1,Y_vec,r))
        }
      } else{
        L = function(vec){
          return(L_quanti_cpp(quanti_matrix,F_matrix,kernel_str,vec,Y_vec,reg))
        }
      }
    } else if(case_str == "quali"){
      if(opt_reg){
        L = function(vec){
          q = length(prior_matrix[1,])
          A_matrix = matrix(vec[1:(length(vec)-1)],nrow = q)
          r = vec[length(vec)]
          return(L_quali_cpp(prior_matrix,F_matrix,kernel_str,A_matrix,Y_vec,r))
        }
      } else{
        L = function(vec){
          q = length(prior_matrix[1,])
          A_matrix = matrix(vec,nrow = q)
          return(L_quali_cpp(prior_matrix,F_matrix,kernel_str,A_matrix,Y_vec,reg))
        }
      }
    } else{
      stop("erreur L_case")
    }
  }
  return(L)
}

opt_vec_optim_bfgs = function(vec,L,control){
  result <- optim(vec,fn = L,control = control,method = "BFGS")
  return(result)
}

opt_vec_optim_lbfgs = function(vec,L,control,low,up){
  result <- optim(vec,fn = L,control = control,method = "L-BFGS-B", lower = low, upper = up)
  return(result)
}

opt_vec_bobyka = function(vec,L,control){
  result <- bobyqa(vec,fn = L,control = control)
  return(result)
}

opt_vec_cmaes = function(vec,L,control,low,up){
  result <- cma_es(par = vec,fn = L,control = control, lower = low, upper = up)
  return(result)
}


OPT = function(N_try,silent,case_str,quanti_matrix,prior_matrix,
               F_matrix,kernel_str,Y_vec,d,reg,amplitude_A,
               amplitude_w,amplitude_reg,opt_reg,optimiseur,
               max_it_optimisation,relative_tol_optimisation,pow_w_10,type_init){
  
  if(case_str == "normal"){
    size_v = length(quanti_matrix[1,]) + length(prior_matrix[1,])*d
    low = c(rep(amplitude_w[1],length(quanti_matrix[1,])),rep(amplitude_A[1],length(prior_matrix[1,])*d))
    up = c(rep(amplitude_w[2],length(quanti_matrix[1,])),rep(amplitude_A[2],length(prior_matrix[1,])*d))
  } else if(case_str == "quanti") {
    size_v = length(quanti_matrix[1,])
    low = rep(amplitude_w[1],size_v)
    up = rep(amplitude_w[2],size_v)
  } else if (case_str == "quali") {
    size_v = length(prior_matrix[1,])*d
    low = rep(amplitude_A[1],size_v)
    up = rep(amplitude_A[2],size_v)
  } else{
    stop("erreur OPT")
  }
  
  if(opt_reg){
    size_v = size_v + 1
    low = c(low,amplitude_reg[1])
    up = c(up,amplitude_reg[2])
  }
  
  if(type_init == "sobol"){
    vec_try = matrix(sobol(N_try, dim = size_v),ncol = size_v) %*% diag(up - low)
    for(i in 1:ncol(vec_try)){
      vec_try[,i] = vec_try[,i] + low[i]
    }
  } 
  
  
  if(type_init =="unif"){
    vec_try = matrix(runif(N_try * size_v, min = rep(low, each = N_try), max = rep(up, each = N_try)), nrow = N_try)
  }
  
  if(type_init == "lhs"){
    
    vec_try <- optimumLHS(n = N_try, k = size_v ) %*% diag(up - low)
    for(i in 1:ncol(vec_try)){
      vec_try[,i] = vec_try[,i] + low[i]
    }
  }
  if(opt_reg){
    vec_try[,size_v] = 10**vec_try[,size_v]
    low[size_v] = 10**low[size_v]
    up[size_v] = 10**up[size_v]
  }
  vec_values = rep(NA,N_try)
  M_points = matrix(NA,nrow = N_try,ncol = size_v)
  L = L_case(case_str,quanti_matrix,prior_matrix,F_matrix,kernel_str,Y_vec,reg,opt_reg,pow_w_10)
  if(optimiseur == "BFGS"){
    control = list()
    if(!is.na(max_it_optimisation)){
      control$maxit = max_it_optimisation
    }
    if(!is.na(relative_tol_optimisation)){
      control$reltol = relative_tol_optimisation
    }
    for(i in 1:N_try){
      try(
        expr = {
          op=opt_vec_optim_bfgs(vec = vec_try[i,],L = L,control)
          vec_values[i] = op$value
          M_points[i,] = op$par
        }, silent = silent
      )
    }
  } else if(optimiseur == "L-BFGS-B"){
    control = list()
    if(!is.na(max_it_optimisation)){
      control$maxit = max_it_optimisation
    }
    if(!is.na(relative_tol_optimisation)){
      control$reltol = relative_tol_optimisation
    }
    for(i in 1:N_try){
      try(
        expr = {
          op=opt_vec_optim_lbfgs(vec = vec_try[i,],L = L,control,low,up)
          vec_values[i] = op$value
          M_points[i,] = op$par
        }, silent = silent
      )
    }
  } else if (optimiseur == "cmaes"){
    control = list(maxit = 200 * length(vec_try[1,])^2)
    if(!is.na(max_it_optimisation)){
      control$maxit = max_it_optimisation
    }
    if(!is.na(relative_tol_optimisation)){
      warning("option relative_tol_optimisation non disponible avec cmaes")
    }
    for(i in 1:N_try){
      try(
        expr = {
          op =opt_vec_cmaes(vec = vec_try[i,],L = L,control,low,up)
          vec_values[i] = op$value
          M_points[i,] = op$par
        }, silent = silent
      )
    }
  } else{
    control = list()
    if(!is.na(max_it_optimisation)){
      control$maxeval = max_it_optimisation
    }
    if(!is.na(relative_tol_optimisation)){
      control$ftol_rel = relative_tol_optimisation
    }
    for(i in 1:N_try){
      try(
        expr = {
          op=opt_vec_bobyka(vec = vec_try[i,],L = L,control)
          vec_values[i] = op$value
          M_points[i,] = op$par
        }, silent = silent
      )
    }
  }
  indices = which(!is.na(vec_values))
  id = indices[which.min(vec_values[indices])]
  liks = -vec_values[indices]
  lik = max(liks)
  res = M_points[id,]
  if(length(res) == 0){
    stop("erreur de convergence, cela peut arriver si Ntry est trop petit, veuillez retenter")
  }
  return(list(res = res, lik = lik, liks = liks))
}

build_LMGP = function(fitted_vec,case_str,quanti_matrix,prior_matrix,kernel_str,F_matrix,Y_vec,reg,lik,opt_reg,pow_w_10){
  if(pow_w_10){
    if(opt_reg){
      reg = fitted_vec[length(fitted_vec)]
      fitted_vec = fitted_vec[1:(length(fitted_vec)-1)]
    }
    if(case_str == "normal"){
      p = length(quanti_matrix[1,])
      A_matrix = matrix(fitted_vec[(p+1):length(fitted_vec)],nrow = length(prior_matrix[1,]))
      w_vec = fitted_vec[1:p]
      post_matrix = prior_matrix %*% A_matrix
      inv_R_fit = invertMatrixChol(cpp_compute_R_normal(quanti_matrix,post_matrix,10**w_vec,kernel_str,reg))
    } else if(case_str == "quanti"){
      A_matrix = NA
      w_vec = fitted_vec
      post_matrix = NA
      inv_R_fit = invertMatrixChol(cpp_compute_R_quanti(quanti_matrix,10**w_vec,kernel_str,reg))
    } else if(case_str == "quali"){
      A_matrix = matrix(fitted_vec,nrow = length(prior_matrix[1,]))
      w_vec = NA
      post_matrix = prior_matrix %*% A_matrix
      inv_R_fit = invertMatrixChol(cpp_compute_R_quali(post_matrix,kernel_str,reg)) 
    } else{
      stop("erreur build_LMGP")
    }
    beta_fit = cpp_beta_hat(F_matrix,inv_R_fit,Y_vec)
    sigma2_fit = cpp_sigma_hat(F_matrix,inv_R_fit,beta_fit,Y_vec)
    inv_V = (1/sigma2_fit) * inv_R_fit
    v_y_f_beta = inv_V%*%(Y_vec - F_matrix %*% beta_fit)
    f_v = t(F_matrix) %*% inv_V
    inv_f_v_f = invertMatrixChol(f_v %*% F_matrix)
  } else{
    if(opt_reg){
      reg = fitted_vec[length(fitted_vec)]
      fitted_vec = fitted_vec[1:(length(fitted_vec)-1)]
    }
    if(case_str == "normal"){
      p = length(quanti_matrix[1,])
      A_matrix = matrix(fitted_vec[(p+1):length(fitted_vec)],nrow = length(prior_matrix[1,]))
      w_vec = fitted_vec[1:p]
      post_matrix = prior_matrix %*% A_matrix
      inv_R_fit = invertMatrixChol(cpp_compute_R_normal(quanti_matrix,post_matrix,w_vec,kernel_str,reg))
    } else if(case_str == "quanti"){
      A_matrix = NA
      w_vec = fitted_vec
      post_matrix = NA
      inv_R_fit = invertMatrixChol(cpp_compute_R_quanti(quanti_matrix,w_vec,kernel_str,reg))
    } else if(case_str == "quali"){
      A_matrix = matrix(fitted_vec,nrow = length(prior_matrix[1,]))
      w_vec = NA
      post_matrix = prior_matrix %*% A_matrix
      inv_R_fit = invertMatrixChol(cpp_compute_R_quali(post_matrix,kernel_str,reg)) 
    } else{
      stop("erreur build_LMGP")
    }
    beta_fit = cpp_beta_hat(F_matrix,inv_R_fit,Y_vec)
    sigma2_fit = cpp_sigma_hat(F_matrix,inv_R_fit,beta_fit,Y_vec)
    inv_V = (1/sigma2_fit) * inv_R_fit
    v_y_f_beta = inv_V%*%(Y_vec - F_matrix %*% beta_fit)
    f_v = t(F_matrix) %*% inv_V
    inv_f_v_f = invertMatrixChol(f_v %*% F_matrix)
  }
  return(list(case_str = case_str,kernel_str = kernel_str, quanti_matrix = quanti_matrix, post_matrix = post_matrix, w_vec = w_vec, A_matrix = A_matrix, beta = beta_fit, sigma2 = sigma2_fit, inv_V = inv_V, v_y_f_beta = v_y_f_beta, inv_f_v_f = inv_f_v_f, f_v = f_v, reg = reg, lik = lik))
}

LMGP_pred = function(fitted_LMGP,new_quanti_matrix,new_prior_matrix,new_F_matrix){
  if(fitted_LMGP$pow_w_10){
    if(fitted_LMGP$case_str == "normal"){
      new_post_matrix = new_prior_matrix %*% fitted_LMGP$A_matrix
      g = cpp_compute_g_normal(new_quanti_matrix,new_post_matrix,fitted_LMGP$quanti_matrix,fitted_LMGP$post_matrix,10**fitted_LMGP$w_vec,fitted_LMGP$kernel_str,fitted_LMGP$sigma2)
    } else if(fitted_LMGP$case_str == "quanti"){
      g = cpp_compute_g_quanti(new_quanti_matrix,fitted_LMGP$quanti_matrix,10**fitted_LMGP$w_vec,fitted_LMGP$kernel_str,fitted_LMGP$sigma2)
    } else if(fitted_LMGP$case_str == "quali"){
      new_post_matrix = new_prior_matrix %*% fitted_LMGP$A_matrix
      g = cpp_compute_g_quali(new_post_matrix,fitted_LMGP$post_matrix,fitted_LMGP$kernel_str,fitted_LMGP$sigma2)
    } else{
      stop("erreur predict_LMGP")
    }
  } else{
    if(fitted_LMGP$case_str == "normal"){
      new_post_matrix = new_prior_matrix %*% fitted_LMGP$A_matrix
      g = cpp_compute_g_normal(new_quanti_matrix,new_post_matrix,fitted_LMGP$quanti_matrix,fitted_LMGP$post_matrix,fitted_LMGP$w_vec,fitted_LMGP$kernel_str,fitted_LMGP$sigma2)
    } else if(fitted_LMGP$case_str == "quanti"){
      g = cpp_compute_g_quanti(new_quanti_matrix,fitted_LMGP$quanti_matrix,fitted_LMGP$w_vec,fitted_LMGP$kernel_str,fitted_LMGP$sigma2)
    } else if(fitted_LMGP$case_str == "quali"){
      new_post_matrix = new_prior_matrix %*% fitted_LMGP$A_matrix
      g = cpp_compute_g_quali(new_post_matrix,fitted_LMGP$post_matrix,fitted_LMGP$kernel_str,fitted_LMGP$sigma2)
    } else{
      stop("erreur predict_LMGP")
    }
  }
  E_y = new_F_matrix %*% fitted_LMGP$beta + t(g) %*% fitted_LMGP$v_y_f_beta
  V_y=cpp_var_gp(g,fitted_LMGP$inv_V,new_F_matrix,fitted_LMGP$f_v,fitted_LMGP$inv_f_v_f,fitted_LMGP$sigma2)
  V_y = V_y + fitted_LMGP$reg * fitted_LMGP$sigma2
  V_y[V_y < 0] = 0
  if(fitted_LMGP$scale){
    E_y = (E_y*fitted_LMGP$sigma_y) + fitted_LMGP$mu_y
    V_y = V_y * (fitted_LMGP$sigma_y**2)
  }
  return(list(E_y = E_y, V_y = V_y))
}
