cut_df = function(df){
  nb_col = ncol(df)
  quant_cols <- sapply(df, is.numeric)
  nb_quanti = sum(quant_cols)
  nb_quali = nb_col - nb_quanti
  if (nb_quanti == nb_col){
    case_str = 'quanti'
    list_df_case_str = list(df_quanti = df, df_quali = NA, case_str = case_str)
    return(list_df_case_str)
  }
  else if (nb_quali == nb_col){
    case_str = 'quali'
    list_df_case_str = list(df_quanti = NA, df_quali = df, case_str = case_str)
    return(list_df_case_str)
  }
  else {
    case_str = 'normal'
    df_quanti <- as.data.frame(df[, quant_cols])
    df_quali <- as.data.frame(df[, !quant_cols])
    list_df_case_str = list(df_quanti = df_quanti, df_quali = df_quali, case_str = case_str)
    return(list_df_case_str)
  }
}

make_columns = function(data){
  cols = character(0)
  list_col_mod = list()
  for (j in 1:ncol(data)){
    mod = unique(data[,j])
    cols = c(cols, paste0(colnames(data)[j], "_&_", mod))
    list_col_mod[[j]] = mod
  }
  return(list(cols, list_col_mod))
}

one_hot = function(data, encodage) {
  ni = nrow(data)
  nj = ncol(data)
  columns = encodage[[1]]
  list_col_mod = encodage[[2]]
  oh = matrix(0, nrow = ni, ncol = length(columns))
  k = 0
  for(j in seq_along(data)){
    mods = list_col_mod[[j]]
    for(i in seq_along(mods)){
      oh[data[,j] == mods[i], k + i] = 1
    }
    k=k+length(mods)
  }
  return(oh)
}

latent = function(model,num=FALSE){
  encod = model$encod
  if(num){
    for(i in 1:length(encod[[2]])){
      encod[[2]][[i]] = seq(1,length(encod[[2]][[i]]))
    }
  }
  nbvar = length(encod[[2]])
  compteur = rep(1,nbvar)
  bornes = sapply(encod[[2]], length)
  M = do.call(expand.grid, encod[[2]])
  M = as.data.frame(sapply(M,as.character))
  Z = one_hot(M,encod)
  return(return(list(labs = M,latent = Z%*%model$A_matrix)))
}


creer_encodage <- function(liste_variables){
  
  ind_var_quali <- lapply( liste_variables, function(sous_liste){
    
    if( sous_liste$type == "quali"){
      
      return(1)
    }
    
    else{ return(0)}
    
  })
  
  liste_var_quali <- liste_variables[which(ind_var_quali == 1)]
  liste_noms_encodage <-  lapply( liste_var_quali, function(sous_liste){
    
    vec_nom <- c()
    
    for(i in 1:sous_liste$max){
      
      vec_nom <- c(vec_nom, paste0(sous_liste$nom, "_", i))
      
    }
    
    liste_levels <- list(as.character(sous_liste$levels))
    
    return(vec_nom)
  })
  
  liste_levels_encodage <-  lapply(liste_var_quali, function(sous_liste){
    
    
    liste_levels <- list(as.character(sous_liste$levels))
    
    return(liste_levels)
  })
  for(i in 1:length(liste_levels_encodage)){
    liste_levels_encodage[[i]] = unlist(liste_levels_encodage[[i]])
  }
  
  
  return(list(unlist(liste_noms_encodage), liste_levels_encodage))
}



