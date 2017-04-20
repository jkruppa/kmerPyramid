##' Calculate the PCA with pyramid edges
##'
##' Calculate the PCA with pyramid edges
##' @title Calculate the PCA with pyramid edges
##' @param df Data frame of k-mer freuencies. The single k-mers are
##'   the columns and the rows indicating different samples, sequence
##'   reads, or contigs.
##' @return A list containing the following elements: 
##' \describe{
##' \item{\emph{pca_plot_df}}{Data frame of the first three principle components from the PCA representing the frequencies of the k-mers}
##' \item{\emph{pca_edge_df}}{Data frame of the first three principle components from the PCA representing the edges of the 3D pyramid}
##' }
##' @author Jochen Kruppa
##' @export
get_pca_plot_list <- function(df){
  edges_mat <- setNames(data.frame(rbind(diag(1,4))), c("A", "C", "G", "T"))
  pca_df <- rbind(edges_mat, df[c("A", "C", "G", "T")])
  pca <- predict(prcomp(pca_df))[,1:3]
  return(list(pca_plot_df = tbl_df(pca[-c(1:4),]),
              pca_edge_df = tbl_df(pca[c(1:4),])))
}
