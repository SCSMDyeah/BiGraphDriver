snf_struct <- function(matrices, num_iterations = 5) {
  n <- length(matrices)

  for (t in 1:num_iterations) {
    new_matrices <- list()
    for (i in 1:n) {
      temp_sum <- matrix(0, nrow = nrow(matrices[[i]]), ncol = ncol(matrices[[i]]))
      for (j in setdiff(1:n, i)) {
        temp_sum <- temp_sum + matrices[[j]]
      }
      new_matrix <- temp_sum / (n - 1)
      new_matrix <- new_matrix %*% t(new_matrix)
      diag(new_matrix) <- 1
      new_matrices[[i]] <- new_matrix
    }
    matrices <- new_matrices
  }
  
  final_matrix <- Reduce("+", matrices) / n
  return(final_matrix)
}

calculation_similarity = function(diff1, diff2, diff11, diff22, patMutMatrix, patMutMatrix1,mu) {
  pathways = read.table('.../data/BLCA.csv', sep = ',', header = TRUE)
  pat_pathways = list()
  for(i in 1:nrow(patMutMatrix)) {
    pathway_mut = list()
    mut_names = names(which(patMutMatrix[i,] != 0))
    for(j in 1:length(mut_names)) {
      ol = grepl(mut_names[j], pathways[,8])
      pathway_mut[[j]] = pathways[which(ol == TRUE), 2]
    }
    pat_pathways[[i]] = unique(unlist(pathway_mut))
  }
  names(pat_pathways) = row.names(patMutMatrix)

  common_samples <- Reduce(intersect, lapply(list(diff1, diff2, diff11, diff22,patMutMatrix,patMutMatrix1), rownames))
  common_genes <- Reduce(intersect, lapply(list(diff1, diff2, diff11, diff22,patMutMatrix,patMutMatrix1), colnames))
  mm <- diff1[common_samples, common_genes, drop = FALSE]
  oo <- diff2[common_samples, common_genes, drop = FALSE]
  mm_methyl <- diff11[common_samples, common_genes, drop = FALSE]
  oo_methyl <- diff22[common_samples, common_genes, drop = FALSE]
  patMutMatrix<- patMutMatrix[common_samples, common_genes, drop = FALSE]
  
  m_dist <- as.matrix(dist(mm))
  uu = matrix(rowMeans(m_dist), nrow = nrow(m_dist), ncol = nrow(m_dist))
  beta = (uu + t(uu) + as.matrix(dist(mm))) / 3
  E_sim = exp(-(m_dist) / (mu * beta))
  diag(E_sim) = 0
  E_sim = E_sim / rowSums(E_sim)
  diag(E_sim) = 1
  
  o_dist <- as.matrix(dist(oo))
  uo = matrix(rowMeans(o_dist), nrow = nrow(o_dist), ncol = nrow(o_dist))
  beta = (uo + t(uo) + as.matrix(dist(oo))) / 3
  o_sim = exp(-(o_dist) / (mu * beta))
  diag(o_sim) = 0
  o_sim = o_sim / rowSums(o_sim)
  diag(o_sim) = 1
  
  mm_methyl_dist <- as.matrix(dist(mm_methyl))
  uu_methyl = matrix(rowMeans(mm_methyl_dist), nrow = nrow(mm_methyl_dist), ncol = nrow(mm_methyl_dist))
  beta_methyl = (uu_methyl + t(uu_methyl) + as.matrix(dist(mm_methyl))) / 3
  E_sim_methyl = exp(-(mm_methyl_dist) / (mu * beta_methyl))
  diag(E_sim_methyl) = 0
  E_sim_methyl = E_sim_methyl / rowSums(E_sim_methyl)
  diag(E_sim_methyl) = 1
  
  o_methyl_dist <- as.matrix(dist(oo_methyl))
  uo_methyl = matrix(rowMeans(o_methyl_dist), nrow = nrow(o_methyl_dist), ncol = nrow(o_methyl_dist))
  beta_methyl = (uo_methyl + t(uo_methyl) + as.matrix(dist(oo_methyl))) / 3
  o_sim_methyl = exp(-(o_methyl_dist) / (mu * beta_methyl))
  diag(o_sim_methyl) = 0
  o_sim_methyl = o_sim_methyl / rowSums(o_sim_methyl)
  diag(o_sim_methyl) = 1
  
  Eo_sim <- snf_struct(list(E_sim, o_sim))
  Eo_sim_methyl <- snf_struct(list(E_sim_methyl, o_sim_methyl))
  final_sim_matrix <- snf_struct(list(Eo_sim, Eo_sim_methyl))
  
  pathway_sim = matrix(0, nrow = nrow(patMutMatrix), ncol = nrow(patMutMatrix), dimnames = list(rownames(patMutMatrix), rownames(patMutMatrix)))
  for (i in 1:nrow(pathway_sim)) {
    index_r = which(names(pat_pathways) == row.names(pathway_sim)[i])
    for (j in 1:ncol(pathway_sim)) {
      index_c = which(names(pat_pathways) == colnames(pathway_sim)[j])
      set_r = pat_pathways[[index_r]]
      set_c = pat_pathways[[index_c]]
      intersection_size = length(intersect(set_r, set_c))
      norm_r = sqrt(length(set_r))
      norm_c = sqrt(length(set_c))
      if (norm_r > 0 && norm_c > 0) {
        pathway_sim[i, j] = intersection_size / (norm_r * norm_c)
      } else {
        pathway_sim[i, j] = 0
      }
    }
  }
  diag(pathway_sim) = 0
 
  sim_pat = pathway_sim * final_sim_matrix
  sim_pat[is.na(sim_pat)] = 0
  sim_pat[rowSums(sim_pat) == 0, ] = 1
  
  final_rates = (rowMeans(mm) + rowMeans(mm_methyl)) + (sim_pat %*% (mm + mm_methyl)) / rowSums(sim_pat)
  ll1 = sort(colSums(final_rates), decreasing = TRUE)
  return(ll1)
}
