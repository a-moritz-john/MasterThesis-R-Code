mxScores <- function(x, control) {
  if (OpenMx::imxHasDefinitionVariable(x)) {
    return(mxScores_df(x = x, control = control))
  } else {
    return(mxScores_standard(x = x, control = control))
  }
}


mxScores_standard <- function(x, control) {
  
  # Whole dataframe is being centered
  data_obs <- x$data$observed[, x$manifestVars, drop = FALSE]
  data_obs <- data_obs[rowSums(!is.na(data_obs)) > 0, , drop = FALSE]
  cd <- scale(x = data_obs, center = TRUE, scale = FALSE)
  
  # Get all parameter names
  full_param_names <- names(x$output$estimate)
  
  # Get missing data pattern 
  # unique_patterns <- control$scores_info$unique_patterns # resulted in out of bounds error
  # getting the info from the missing data pattern caused out of bounds errors
  pattern_string  <- apply(is.na(data_obs), 1, paste0, collapse = "")
  unique_patterns <- unique(pattern_string)
  
  # Create final score matrix
  scores <- matrix(0, nrow = nrow(data_obs), ncol = length(full_param_names))
  colnames(scores) <- full_param_names
  
  # Loop that iterates through all the missing data patterns
  for (pat in unique_patterns) {
    
    # Get row indices (individuals with this missing data pattern)
    row_idx <- which(pattern_string == pat)
    
    # Get column indices (variables that are observed for this pattern)
    col_idx <- which(!is.na(data_obs[row_idx[1], ]))
    
    # Extract the submatrix with only observed data
    data_obs_sub <- data_obs[row_idx, col_idx, drop = FALSE]
    cd_sub <- cd[row_idx, col_idx, drop = FALSE]
    
    
    p <- length(col_idx)# control$scores_info$p
    mean_structure <- any(x$M$free[, col_idx, drop = FALSE]) # control$scores_info$mean_structure
    p_star <- p * (p + 1) / 2 # control$scores_info$p_star
    p_star_seq <- seq_len(p_star)
    p_star_means <- p * (p + 3) / 2 # control$scores_info$p_star_means
    
    exp_cov <- OpenMx::mxGetExpected(model = x, component = "covariance")[col_idx, col_idx, drop = FALSE]
    exp_cov_inv <- solve(exp_cov)
    # data_obs <- x$data$observed[, x$manifestVars, drop = FALSE][, control$scores_info$not_missing, drop = FALSE] # Is not needed anymore since creation and subsetting already happened before
    N <- nrow(data_obs_sub)
    
    if (control$linear) {
      
      # Filter parameters by missing data pattern
      all_var_names <- colnames(data_obs)  
      
      # Get variable names for this pattern
      var_names <- all_var_names[col_idx] 
      
      # Map to positions in the full variable list
      var_indices <- match(var_names, all_var_names)
      
      # Create vech-style lower-triangle covariance param names
      cov_pairs <- expand.grid(i = var_indices, j = var_indices)
      cov_pairs <- cov_pairs[cov_pairs$i <= cov_pairs$j, ]
      cov_param_names <- paste0("covariance_model_S", cov_pairs$i, cov_pairs$j)
      
      # Create mean param names
      mean_param_names <- paste0("covariance_model_M1", var_indices)
      
      # Combine both
      param_names_needed <- c(cov_param_names, mean_param_names)
      
      # Get the indices of the relevant parameters
      active_param_indices <- which(full_param_names %in% param_names_needed)
      
      
      q <- length(active_param_indices)# control$scores_info$q # Still needs to be fixed -> has to be subsetted according to the specififc missing data pattern
      q_seq <- active_param_indices # control$scores_info$q_seq # Depends in line above
      p_unf <- NROW(x$A$values[col_idx, col_idx, drop = FALSE]) # control$scores_info$p_unf # 
      
      A_deriv <- control$scores_info$A_deriv # Probably needs to be subsetted -> in function before the whole A deriv is being created
      S_deriv <- control$scores_info$S_deriv # Same as above
      m_deriv <- control$scores_info$m_deriv # Same as above
      
      F_RAM <- x$F$values[col_idx, col_idx, drop = FALSE] # Moritz figures this out (does not work with latent  variables)
      m <- t(x$M$values[1, col_idx, drop = FALSE])
      A <- x$A$values[col_idx, col_idx, drop = FALSE]
      B <- solve(diag(x = 1, nrow = NROW(A)) - A)
      E <- B %*% x$S$values[col_idx, col_idx, drop = FALSE] %*% t(B)
      FB <- F_RAM %*% B
      
      jac <- matrix(0, nrow = p_star_means, ncol = q)
      
      for (k in seq_along(q_seq)) { 
        
        i <- q_seq[k] # index of the correct matrix derivative according to missing data pattern
        j <- seq_len(q)[k] # index for the jacobian column index
        
        symm <- FB %*% A_deriv[[i]][col_idx, col_idx] %*% E %*% t(F_RAM)
        jac[p_star_seq, j] <- lavaan::lav_matrix_vech(symm + t(symm) + FB %*% S_deriv[[i]][col_idx, col_idx] %*% t(FB))
      }
      
      for (k in seq_along(q_seq)) {
        
        i <- q_seq[k] # index of the correct matrix derivative according to missing data pattern
        j <- seq_len(q)[k] # index for the jacobian column index
        
        jac[(p_star+1):p_star_means, j] <- FB %*% A_deriv[[i]][col_idx, col_idx] %*% B %*% m +
          FB %*% m_deriv[[i]][col_idx,]
      }
      
      colnames(jac) <- param_names_needed
      
    } else {
      
      jac <- OpenMx::omxManifestModelByParameterJacobian(model = x)[q_seq, q_seq]
      
    }
    
    if (mean_structure == FALSE) {jac <- jac[p_star_seq, , drop = FALSE]}
    
    # Calculate weight matrix
    Dup <- lavaan::lav_matrix_duplication(n = p)
    V <- 0.5 * t(Dup) %*% kronecker(X = exp_cov_inv, Y = exp_cov_inv) %*% Dup
    if (mean_structure) {
      V_m_cov <- matrix(data = 0, nrow = p_star_means, ncol = p_star_means)
      V_m_cov[p_star_seq, p_star_seq] <- V
      V_m_cov[(p_star + 1):p_star_means, (p_star + 1):p_star_means] <- exp_cov_inv
      V <- V_m_cov
    }
    
    # Individual deviations from the sample moments
    # cd <- scale(x = data_obs, center = TRUE, scale = FALSE) # Not needed anymore
    if (p == 1) {
      mc <- matrix(apply(X = cd_sub, MARGIN = 1,
                         FUN = function(x) {lavaan::lav_matrix_vech(x %*% t(x))}))
    } else {
      mc <- t(apply(X = cd_sub, MARGIN = 1,
                    FUN = function(x) {lavaan::lav_matrix_vech(x %*% t(x))}))
    }
    vech_cov <- matrix(data = rep(x = lavaan::lav_matrix_vech(exp_cov), times = N),
                       byrow = TRUE, nrow = N, ncol = p_star)
    md <- mc - vech_cov
    if (mean_structure) {
      exp_means <- OpenMx::mxGetExpected(model = x, component = "means")[col_idx]
      means <- matrix(data = rep(x = exp_means, times = N), byrow = TRUE,
                      nrow = N, ncol = p)
      mean_dev <- data_obs_sub - means
      md <- as.matrix(cbind(md, mean_dev))
    }
    
    # Calculates scores
    scores_sub  <- md %*% V %*% jac
    
    # Insert scores into final matrix
    person_ids <- as.vector(row_idx) # take the actual positions and not the names
    param_cols <- colnames(scores_sub)
    
    scores[person_ids, param_cols] <- as.matrix(scores_sub)
    
  }
  
  return(scores)
  
}


mxScores_df <- function(x, control) {
  
  p <- control$scores_info$p
  mean_structure <- control$scores_info$mean_structure
  p_star <- control$scores_info$p_star
  p_star_seq <- seq_len(p_star)
  p_star_means <- control$scores_info$p_star_means
  p_star_p_means_seq <- (p_star + 1):p_star_means
  p_unf <- control$scores_info$p_unf
  
  data_obs <- x$data$observed[, x$manifestVars, drop = FALSE]
  N <- nrow(data_obs)
  Ident <- diag(x = 1, nrow = p_unf)
  Dup <- lavaan::lav_matrix_duplication(n = p)
  
  # Scores
  scores <- matrix(NA, nrow = N, ncol = control$scores_info$q)
  colnames(scores) <- names(x$output$estimate)
  
  # Individual sample moments
  cd <- scale(x = data_obs, center = TRUE, scale = FALSE)
  if (p == 1) {
    mc <- matrix(apply(X = cd, MARGIN = 1,
                       FUN = function(x) {lavaan::lav_matrix_vech(x %*% t(x))}))
  } else {
    mc <- t(apply(X = cd, MARGIN = 1,
                  FUN = function(x) {lavaan::lav_matrix_vech(x %*% t(x))}))
  }
  
  # Prepare empty object for the individual deviations from the sample moments
  if (mean_structure) {
    md <- matrix(NA, nrow = N, ncol = p_star_means)
  } else {
    md <- matrix(NA, nrow = N, ncol = p_star)
  }
  
  # Get definition variables
  df <- identify_definition_variables(x)
  df_labels <- df[[1]]
  df_nr <- length(df_labels)
  df_data <- x$data$observed[, df[[2]], drop = FALSE]
  df_indices <- seq_along(df_labels)
  
  ## RAM matrices
  F_RAM <- x$F$values
  A <- x$A$values
  S <- x$S$values
  m <- t(x$M$values)
  B <- solve(Ident - A)
  FB <- F_RAM %*% B
  E <- B %*% S %*% t(B)
  
  if (control$linear) {
    
    q <- control$scores_info$q
    q_seq <- control$scores_info$q_seq
    
    A_deriv <- control$scores_info$A_deriv
    S_deriv <- control$scores_info$S_deriv
    m_deriv <- control$scores_info$m_deriv
    
    jac <- matrix(0, nrow = p_star_means, ncol = q)
    
  }
  
  ## Assign definition variables to the corresponding RAM matrices
  RAM_df <- rep(NA, times = df_nr)
  RAM_df[which(df_labels %in% x$A$labels)] <- "A"
  RAM_df[which(df_labels %in% x$S$labels)] <- "S"
  RAM_df[which(df_labels %in% x$M$labels)] <- "M"
  
  # Get coordinates of the definition variables in the corresponding RAM matrices
  RAM_coord <- list()
  for (j in df_indices) {
    if (RAM_df[j] == "A") {
      RAM_coord[[j]] <- which(x$A$labels == df_labels[j], arr.ind = TRUE)
    }
    if (RAM_df[j] == "S") {
      RAM_coord[[j]] <- which(x$S$labels == df_labels[j], arr.ind = TRUE)
    }
    if (RAM_df[j] == "M") {
      RAM_coord[[j]] <- which(t(x$M$labels) == df_labels[j], arr.ind = TRUE)
    }
  }
  
  # Get groups of individuals with identical definition variables
  group <- transform(df_data,
                     group_ID = as.numeric(interaction(df_data,
                                                       drop = TRUE)))
  unique_groups <- unique(group$group_ID)
  
  # Loop over the all the different definition variable values
  for (i in unique_groups) { # Start loop with index i
    
    group_rows <- which(group$group_ID == i) # CHECK IF NEEDED
    group_n <- NROW(group_rows)
    df_values <- as.numeric(group[group$group_ID == i, ][1, ])
    
    # Update the definition variable values in the RAM matrices
    for (j in df_indices){
      
      if (RAM_df[j] == "A") {
        A[RAM_coord[[j]]] <- df_values[j]
      }
      
      if (RAM_df[j] == "S") {
        S[RAM_coord[[j]]] <- df_values[j]
      }
      
      if (RAM_df[j] == "M") {
        m[RAM_coord[[j]]] <- df_values[j]
      }
    } # end loop with index j
    
    # Update sample covariance
    B <- solve(Ident - A)
    FB <- F_RAM %*% B
    E <- B %*% S %*% t(B)
    exp_cov <- F_RAM %*% E %*% t(F_RAM)
    exp_cov_inv <- solve(exp_cov)  
    
    # Update Jacobian matrix
    if (control$linear) { # Analytic Jacobian matrix
      
      for (j in seq_len(q)) {
        symm <- FB %*% A_deriv[[j]] %*% E %*% t(F_RAM)
        jac[p_star_seq, j] <- lavaan::lav_matrix_vech(symm + t(symm) + FB %*% S_deriv[[j]] %*% t(FB))
      }
      
      if (mean_structure) {
        for (j in seq_len(q)) {
          jac[(p_star+1):p_star_means, j] <- FB %*% A_deriv[[j]] %*% B %*% m +
            FB %*% m_deriv[[j]]
        }
      }
      
      
    } else { # Numeric Jacobian matrix
      x <- OpenMx::omxSetParameters(model = x, labels = df$labels,
                                    values = df_values[df_indices])
      x <- suppressMessages(OpenMx::mxRun(model = x, useOptimizer = FALSE, silent=TRUE))
      jac <- OpenMx::omxManifestModelByParameterJacobian(model = x)
    }
    
    if (mean_structure == FALSE) {jac <- jac[p_star_seq, , drop = FALSE]}
    
    
    V <- 0.5 * t(Dup) %*% kronecker(X = exp_cov_inv, Y = exp_cov_inv) %*% Dup
    if (mean_structure) {
      V_m_cov <- matrix(data = 0, nrow = p_star_means, ncol = p_star_means)
      V_m_cov[p_star_seq, p_star_seq] <- V
      V_m_cov[p_star_p_means_seq, p_star_p_means_seq] <- exp_cov_inv
      V <- V_m_cov
    }
    
    # Individual deviations from the sample moments
    vech_cov <- matrix(data = rep(x = lavaan::lav_matrix_vech(exp_cov),
                                  times = group_n),
                       byrow = TRUE, nrow = group_n, ncol = p_star)
    md[group_rows, 1:p_star] <- mc[group_rows, ] - vech_cov
    if (mean_structure) {
      means <- matrix(data = rep(x = FB %*% m, times = group_n), byrow = TRUE,
                      nrow = group_n, ncol = p)
      means_dev <- as.matrix(data_obs[group_rows, ]) - means
      md[group_rows, (p_star+1):p_star_means] <- means_dev
    }
    
    scores[group_rows, ]  <- md[group_rows, ] %*% V %*% jac
  }
  
  scores
  
}


identify_definition_variables <- function(x) {
  definition_variables <- c()
  for (i in 1:length(x@matrices)) {
    definition_variables <- c(definition_variables,
                              sapply(x@matrices[[i]]$labels,
                                     OpenMx::imxIsDefinitionVariable))
  }
  definition_variables <- names(which(definition_variables))
  list(labels = definition_variables, # OpenMx labels
       data = sub(".*\\.", "", definition_variables)) # column names
}

