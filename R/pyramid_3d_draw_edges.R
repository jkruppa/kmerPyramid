##' Function to draw the edges of the 3D pyramid
##'
##' Function to draw the edges of the 3D pyramid
##' @title Function to draw the edges of the 3D pyramid
##' @param pca_edge_df Data frame created by the
##'   \code{\link{get_pca_plot_list}} function for drawing the edges
##'   of the 3D pyramid.
##' @return NULL
##' @author Jochen Kruppa
pyramid_3d_draw_edges <- function(pca_edge_df, alpha = 1){
  lines3d(pca_edge_df[c(1,2), c("PC1", "PC2", "PC3")], col = 1)
  lines3d(pca_edge_df[c(1,3), c("PC1", "PC2", "PC3")], col = 1)
  lines3d(pca_edge_df[c(1,4), c("PC1", "PC2", "PC3")], col = 1)
  lines3d(pca_edge_df[c(2,3), c("PC1", "PC2", "PC3")], col = 1)
  lines3d(pca_edge_df[c(2,4), c("PC1", "PC2", "PC3")], col = 1)
  lines3d(pca_edge_df[c(3,4), c("PC1", "PC2", "PC3")], col = 1)
  text3d(pca_edge_df[1, "PC1"], pca_edge_df[1, "PC2"], pca_edge_df[1, "PC3"], "A",
         color="black", alpha = alpha)
  text3d(pca_edge_df[2, "PC1"], pca_edge_df[2, "PC2"], pca_edge_df[2, "PC3"], "C",
         color="black", alpha = alpha)
  text3d(pca_edge_df[3, "PC1"], pca_edge_df[3, "PC2"], pca_edge_df[3, "PC3"], "G",
         color="black", alpha = alpha)
  text3d(pca_edge_df[4, "PC1"], pca_edge_df[4, "PC2"], pca_edge_df[4, "PC3"], "T",
         color="black", alpha = alpha)
}
