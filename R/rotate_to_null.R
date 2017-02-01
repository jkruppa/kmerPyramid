##' Experimental function for rotating the PCA coordinates
##'
##' Experimental function for rotating the PCA coordinates
##' @title Experimental function for rotating the PCA coordinates 
##' @param list 
##' @return List of rotated PCA coordinates 
##' @author Jochen Kruppa
##' @examples
##' 
##' ### edges
##' A <- c(1, 1, 1)
##' C <- c(1, 1, 2)
##' G <- c(1, 2, 1)
##' T <- c(2, 1, 1)
##' ### Further points
##' X <- matrix(rnorm(20*3), 20, 3)
##
##' rotate_to_null(A, C, G, T, X)
rotate_to_null <- function(list) {
  ## rotate to 0
  ## pca_rotated <- rotate_to_null(list[[1]]$pca_plot_list)
  ## ## start plotting
  ## pyramid_3d_draw_edges(pca_rotated$pca_edge_df)
  ## lines3d(pca_rotated$pca_plot_df$PC1,
  ##         pca_rotated$pca_plot_df$PC2,
  ##         pca_rotated$pca_plot_df$PC3,
  ##         col = color)
  ## draw the sphere
  ## spheres3d(list[[1]]$pca_plot_sphere_df$PC1,
  ##           list[[1]]$pca_plot_sphere_df$PC2,
  ##           list[[1]]$pca_plot_sphere_df$PC3,
  ##           alpha = 0.5,
  ##           col = color,
  ##           radius = list[[1]]$pca_plot_sphere_df$count_radius)
  ## assign the values
  A <- as.numeric(list$pca_edge_df[1,])
  C <- as.numeric(list$pca_edge_df[2,])
  G <- as.numeric(list$pca_edge_df[3,])
  T <- as.numeric(list$pca_edge_df[4,])
  X <- as.matrix(list$pca_plot_df)
  ## start function
  A = A - A
  C = C - A
  G = G - A
  T = T - A
  X = X - A
  ## Drehung um Y-Achse ###
  H1 = c(C[1], A[2], A[3])
  H2 = c(C[1], A[2], C[3])
  a = sqrt((A[1] - H2[1])^2 + (A[2] - H2[2])^2 + (A[3] - H2[3])^2)
  b = sqrt((H1[1] - H2[1])^2 + (H1[2] - H2[2])^2 + (H1[3] - H2[3])^2)
  c = sqrt((A[1] - H1[1])^2 + (A[2] - H1[2])^2 + (A[3] - H1[3])^2)
  alpha = acos((b^2 + c^2 - a^2)/(2*b*c))
  beta = acos((a^2 + c^2 - b^2)/(2*a*c))
  gamma = acos((a^2 + b^2 - c^2)/(2*a*b))
  R = diag(3)
  R[1,1] = cos(beta)
  R[3,3] = cos(beta)
  R[1,3] = sin(beta)
  R[3,1] = -sin(beta)
  A = R %*% A
  C = R %*% C
  G = R %*% G
  T = R %*% T
  X = t(R %*% t(X))
  ## Drehung um die Z-Achse ###
  H1 = c(C[1], A[2], A[3])
  a = sqrt((A[1] - C[1])^2 + (A[2] - C[2])^2 + (A[3] - C[3])^2)
  b = sqrt((H1[1] - C[1])^2 + (H1[2] - C[2])^2 + (H1[3] - C[3])^2)
  c = sqrt((A[1] - H1[1])^2 + (A[2] - H1[2])^2 + (A[3] - H1[3])^2)
  alpha = acos((b^2 + c^2 - a^2)/(2*b*c))
  beta = acos((a^2 + c^2 - b^2)/(2*a*c))
  gamma = acos((a^2 + b^2 - c^2)/(2*a*b))
  R = diag(3)
  R[1,1] = cos(beta)
  R[2,2] = cos(beta)
  R[1,2] = sin(beta)
  R[2,1] = -sin(beta)
  A = R %*% A
  C = R %*% C
  G = R %*% G
  T = R %*% T
  X = t(R %*% t(X))
  ## Drehung um die X-Achse ###
  a = sqrt((A[1] - G[1])^2 + (A[2] - G[2])^2 + (A[3] - G[3])^2)
  b = sqrt((G[1] - C[1])^2 + (G[2] - C[2])^2 + (G[3] - C[3])^2)
  c = sqrt((A[1] - C[1])^2 + (A[2] - C[2])^2 + (A[3] - C[3])^2)
  height = sqrt(2 * (a^2*b^2+b^2*c^2+c^2*a^2) - (a^4+b^4+c^4)) / (2*c)
  H1 = c(G[1], A[2], A[3])
  H2 = c(G[1], A[2], A[3] + height)
  a = sqrt((G[1] - H2[1])^2 + (G[2] - H2[2])^2 + (G[3] - H2[3])^2)
  b = sqrt((H1[1] - H2[1])^2 + (H1[2] - H2[2])^2 + (H1[3] - H2[3])^2)
  c = sqrt((G[1] - H1[1])^2 + (G[2] - H1[2])^2 + (G[3] - H1[3])^2)
  alpha = acos((b^2 + c^2 - a^2)/(2*b*c))
  beta = acos((a^2 + c^2 - b^2)/(2*a*c))
  gamma = acos((a^2 + b^2 - c^2)/(2*a*b))
  R = diag(3)
  R[2,2] = cos(alpha)
  R[3,3] = cos(alpha)
  R[3,2] = sin(alpha)
  R[2,3] = -sin(alpha)
  A = R %*% A
  C = R %*% C
  G = R %*% G
  T = R %*% T
  X = t(R %*% t(X))
  return(list(pca_plot_df = tibble(PC1 = X[,1],
                                   PC2 = X[,2],
                                   PC3 = X[,3]),
              pca_edge_df = setNames(tbl_df(rbind(t(A), t(C), t(G), t(T))),
                                     c("PC1", "PC2", "PC3"))))
  ## return(list(A=A, C=C, G=G, T=T, X=X))
}
