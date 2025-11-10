OpenMx_scores_input <- function(x, control) {
  
  # Determine missing data pattern
  data_obs <- x$data$observed[, x$manifestVars, drop = FALSE]
  data_obs <- data_obs[rowSums(!is.na(data_obs)) > 0, , drop = FALSE]
  
  # Step 1: Create missingness pattern string per row
  pattern_string <- apply(is.na(data_obs), 1, paste0, collapse = "")
  
  # Step 2: Identify unique patterns
  unique_patterns <- unique(pattern_string)
  
  row_groups <- list()
  col_groups <- list()
  
  for (pat in unique_patterns) {
    row_idx <- which(pattern_string == pat)
    col_idx <- which(!is.na(data_obs[row_idx[1], ]))  # assume all rows in pattern have same NA structure
    
    row_groups[[pat]] <- row_idx
    col_groups[[pat]] <- col_idx
  }
  
  names(row_groups) <- paste0("pattern_", seq_along(unique_patterns))
  names(col_groups) <- paste0("pattern_", seq_along(unique_patterns))
  
  
  p <- length(x$manifestVars)
  mean_structure <- any(x$M$free)
  p_star <- p * (p + 1) / 2
  p_star_means <- p * (p + 3) / 2
  
  # AB: give pseudo-labels to matrices if
  # unlabelled parameters are given
  candidate_param_id <- which(startsWith(x=names(x$output$estimate), prefix=x$name))
  if (length(candidate_param_id)>0) {
    for (k in candidate_param_id) {
      candidate_param_name <- names(x$output$estimate)[k]
      cplen <- nchar(x$name)
      candidate_matrix <- substr(candidate_param_name, cplen+2,cplen+2)
      candidate_pos <- as.integer(strsplit(substr(candidate_param_name, cplen+4, nchar(candidate_param_name)-1),",")[[1]])
      if (candidate_matrix=="A") {
        x$A$labels[candidate_pos[1], candidate_pos[2]]<-candidate_param_name        
      } else if (candidate_matrix=="S") {
        x$S$labels[candidate_pos[1], candidate_pos[2]]<-candidate_param_name        
      } else if (candidate_matrix == "M") {
        x$M$labels[candidate_pos[1], candidate_pos[2]]<-candidate_param_name
        # Changes made to the original function to fix bug
        # x$M$labels[candidate_pos]<-candidate_param_name        
      }
    }
  }
  
  if (control$linear | imxHasDefinitionVariable(x)) {
    
    param_names <- names(x$output$estimate)
    q <- length(param_names)
    
    if (q==0) {
      ui_stop("Error in OpenMx_scores_input() function. There are no free parameters in the model estimates. Model not run or converged?")
    }
    
    q_seq <- seq_len(q)
    p_unf <- NROW(x$A$values)
    
    Zero <- matrix(0, nrow = p_unf, ncol = p_unf)
    A_deriv <- lapply(q_seq, function(x) {Zero})
    S_deriv <- A_deriv
    zero <- matrix(0, nrow = p_unf, ncol = 1)
    m_deriv <- lapply(q_seq, function(x) {zero})
    
    for (i in q_seq) {
      A_deriv[[i]][which(x$A$labels == param_names[i], arr.ind = TRUE)] <- 1
    }
    
    for (i in q_seq) {
      # S_deriv[[i]][which(x$S$labels == param_names[i], arr.ind = TRUE)] <- 1
      # Changes made to the original function to fix bugs
      match_pos <- which(x$S$labels == param_names[i], arr.ind = TRUE)
      
      for (j in seq_len(nrow(match_pos))) {
        row <- match_pos[j, 1]
        col <- match_pos[j, 2]
        
        S_deriv[[i]][row, col] <- 1
        S_deriv[[i]][col, row] <- 1  
      }
    }
    
    for (i in q_seq) {
      m_deriv[[i]][which(x$M$labels == param_names[i])] <- 1
    }
    
    scores_info <- list(p = p, mean_structure = mean_structure, p_star = p_star,
                        p_star_means = p_star_means, q = q, q_seq = q_seq,
                        p_unf = p_unf, A_deriv = A_deriv, S_deriv = S_deriv,
                        m_deriv = m_deriv, row_groups = row_groups, col_groups = col_groups, unique_patterns = unique_patterns, pattern_string = pattern_string)
    
  } else {
    
    scores_info <- list(p = p, mean_structure = mean_structure,
                        p_star = p_star, p_star_means = p_star_means,
                        row_groups = row_groups,
                        col_groups = col_groups,
                        unique_patterns = unique_patterns, pattern_string = pattern_string)
  }
  
  return(scores_info)
  
}