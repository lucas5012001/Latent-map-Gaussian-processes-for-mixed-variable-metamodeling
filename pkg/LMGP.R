library(dplyr)
library(Rcpp)
library(nloptr)
library(cmaes)
source("LMGP/sous_modules/fit_predict.R")
source("LMGP/sous_modules/utils.R")

Rcpp::sourceCpp("LMGP/sous_modules/fit_predict.cpp")

LMGP.fit = function(data_frame_type, 
                    Y_vec, 
                    N_try = 12,
                    kernel_str = c("exp","mat"), 
                    F_matrix = matrix(1,nrow = length(Y_vec), ncol = 1), 
                    ncol_A = 2, 
                    reg = NA, 
                    silent = TRUE,
                    scale = TRUE, 
                    amplitude_A = c(-1,1), 
                    amplitude_w = c(-8,3),
                    amplitude_reg = c(-6,-1),
                    optimiseur = c("L-BFGS-B","BFGS","bobyka","cmaes"),
                    max_it_optimisation = NA,
                    relative_tol_optimisation = NA, 
                    pow_w_10 = TRUE,
                    type_init = "lhs"){
  
  if(optimiseur[1] == "optim"){
    optimiseur = "L-BFGS-B"
  }
  
  kernel_str = kernel_str[1]
  optimiseur = optimiseur[1]
  mu_y = mean(Y_vec)
  sigma_y = sd(Y_vec)
  if(scale){
    Y_vec = (Y_vec - mu_y)/sigma_y
  }
  df_list = cut_df(data_frame_type)
  quant_matrix = df_list[[1]]
  qual_matrix = df_list[[2]]
  case_str = df_list[[3]]
  rm(df_list)
  if(case_str == "normal"){
    quant_matrix = as.matrix(quant_matrix)
    encod = make_columns(qual_matrix)
    prior_matrix = one_hot(qual_matrix,encod)
    if(scale){
      quant_matrix = scale(quant_matrix)
      scales_quanti = attr("scaled:scale",x = quant_matrix)
      means_quanti = attr("scaled:center",x = quant_matrix)
      prior_matrix = scale(prior_matrix)
      scales_quali = attr("scaled:scale",x = prior_matrix)
      means_quali = attr("scaled:center",x = prior_matrix)
    }
    qual_matrix = NA
  } else if(case_str == "quanti"){
    quant_matrix = as.matrix(quant_matrix)
    if(scale){
      quant_matrix = scale(quant_matrix)
      scales_quanti = attr("scaled:scale",x = quant_matrix)
      means_quanti = attr("scaled:center",x = quant_matrix)
      scales_quali = NA
      means_quali = NA
    }
    prior_matrix = NA
    qual_matrix = NA
    encod = NA
  } else{
    quant_matrix = NA
    encod = make_columns(qual_matrix)
    prior_matrix = one_hot(qual_matrix,encod)
    if(scale){
      prior_matrix = scale(prior_matrix)
      scales_quali = attr("scaled:scale",x = prior_matrix)
      means_quali = attr("scaled:center",x = prior_matrix)
      scales_quanti = NA
      scales_quali = NA
    }
    qual_matrix = NA
  }
  opt_reg = is.na(reg)
  opt_list = OPT(N_try,
                 silent,
                 case_str,
                 quant_matrix,
                 prior_matrix,
                 F_matrix,
                 kernel_str,
                 Y_vec,
                 ncol_A,
                 reg,
                 amplitude_A,
                 amplitude_w,
                 amplitude_reg,
                 opt_reg,
                 optimiseur,
                 max_it_optimisation,
                 relative_tol_optimisation,
                 pow_w_10,
                 type_init)
  best_parameters = opt_list$res
  lik = opt_list$lik
  model = build_LMGP(best_parameters,
                     case_str,
                     quant_matrix,
                     prior_matrix,
                     kernel_str,
                     F_matrix,
                     Y_vec,
                     reg,
                     lik,
                     opt_reg,
                     pow_w_10)
  model$liks = opt_list$liks
  model$warning = "ok"
  model$pow_w_10 = pow_w_10
  if(model$reg < 0){
    ## On recommence sans regularisation.
    opt_reg = FALSE
    reg = 0
    opt_list = OPT(N_try,
                   silent,
                   case_str,
                   quant_matrix,
                   prior_matrix,
                   F_matrix,
                   kernel_str,
                   Y_vec,
                   ncol_A,
                   reg,
                   amplitude_A,
                   amplitude_w,
                   amplitude_reg,
                   opt_reg,
                   optimiseur,
                   max_it_optimisation,
                   relative_tol_optimisation,
                   pow_w_10,
                   type_init)
    best_parameters = opt_list$res
    lik = opt_list$lik
    model = build_LMGP(best_parameters,
                       case_str,
                       quant_matrix,
                       prior_matrix,
                       kernel_str,
                       F_matrix,
                       Y_vec,
                       reg,
                       lik,
                       opt_reg,
                       pow_w_10)
    model$liks = opt_list$liks
    model$warning = "reg inutile"
    model$pow_w_10 = pow_w_10
    ##
  }
  if(scale){
    model$scales_quali = scales_quali
    model$means_quali = means_quali
    model$scales_quanti = scales_quanti
    model$means_quanti = means_quanti
    model$mu_y = mu_y
    model$sigma_y = sigma_y
  }
  model$scale = scale
  model$encod = encod
  check_pred_train = as.vector(LMGP.predict(model,data_frame_type,F_matrix)[[1]])
  if(scale){
    check_pred_train = (check_pred_train - mu_y)/sigma_y
    MSE1 = sum((check_pred_train - Y_vec)**2)
    MSE2 = sum((0 - Y_vec)**2)
  } else{
    MSE1 = sum((check_pred_train - Y_vec)**2)
    MSE2 = sum((mu_y - Y_vec)**2)
  }
  if(MSE1 > MSE2){
    stop("erreur de convergence, la moyenne predit mieux que le modele")
  }
  model$MSE = MSE1
  return(model)
}

LMGP.predict = function(fitted_LMGP, new_data_frame_type, new_F_matrix = matrix(1,nrow = length(new_data_frame_type[,1]), ncol = 1)){
  new_df_list = cut_df(new_data_frame_type)
  new_quant_matrix = new_df_list[[1]]
  new_qual_matrix = new_df_list[[2]]
  rm(new_df_list)
  if(fitted_LMGP$case_str == "normal"){
    new_quant_matrix = as.matrix(new_quant_matrix)
    new_prior_matrix = one_hot(new_qual_matrix,fitted_LMGP$encod)
    if(fitted_LMGP$scale){
      new_quant_matrix = scale(new_quant_matrix,center = fitted_LMGP$means_quanti, scale = fitted_LMGP$scales_quanti)
      new_prior_matrix = scale(new_prior_matrix,center = fitted_LMGP$means_quali, scale = fitted_LMGP$scales_quali)
    }
    new_qual_matrix = NA
  } else if(fitted_LMGP$case_str == "quanti"){
    new_quant_matrix = as.matrix(new_quant_matrix)
    if(fitted_LMGP$scale){
      new_quant_matrix = scale(new_quant_matrix,center = fitted_LMGP$means_quanti, scale = fitted_LMGP$scales_quanti)
    }
    new_prior_matrix = NA
    new_qual_matrix = NA
  } else{
    new_quant_matrix = NA
    new_prior_matrix = one_hot(new_qual_matrix,fitted_LMGP$encod)
    if(fitted_LMGP$scale){
      new_prior_matrix = scale(new_prior_matrix,center = fitted_LMGP$means_quali, scale = fitted_LMGP$scales_quali) 
    }
    new_qual_matrix = NA
  }
  predictions = LMGP_pred(fitted_LMGP,new_quant_matrix,new_prior_matrix,new_F_matrix)
  return(predictions)
}
