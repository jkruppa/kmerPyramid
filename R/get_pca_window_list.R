##' Get the PCA values for sliding window
##'
##' Get the PCA values for sliding window
##' @title Get the PCA values for sliding window
##' @param seqs DNAStringSet of sequences
##' @param window Size of the sliding window
##' @return A list per sequence containing the following elements:
##'   \describe{
##' \item{\emph{pca_plot_list}}{Includes a list with the
##'   first three principal components of the k-mer frequencies and the
##'   first three principal components of the edges of the 3D
##'   pyramid.}
##' \item{\emph{pca_plot_sphere_df}}{Returns a data frame
##'   with information on the frist three principal components and the
##'   count number of single windows. Further the plotting radius for
##'   the \code{\link{pyramid_3d_window}} function is returned.}
##' }
##' @author Jochen Kruppa
##' @export
##' @examples
##' data(viralExampleSeqs)
##' 
##' viral_window_list <- get_pca_window_list(viralExampleSeqs, window = 2)
get_pca_window_list <- function(seqs, window = 5) {
  p_load(plyr, dplyr, Biostrings, ShortRead, stringr)
  pca_list <- llply(seqs, function(x, win_size = window) {
    ## get the kmer sliding window distribution
    window_kmer_prob_df <- get_window_distr(x, window = win_size, prob = TRUE)
    pca_plot_list <- get_pca_plot_list(window_kmer_prob_df)
    ## get the numer of count for each k-mer
    window_kmer_count_df <- get_window_distr(x, window = win_size, prob = FALSE)
    window_kmer_count_df$ind <- laply(1:nrow(window_kmer_count_df), function(i) {
      str_c(rep(c("A", "C", "G", "T"), window_kmer_count_df[i,1:4]), collapse = "")
    })
    pca_sphere_df <- data.frame(pca_plot_list$pca_plot_df, window_kmer_count_df)
    ## ## get the radius for the difference sphere
    pca_sphere_df$count <- 1
    pca_sphere_df$sort_ind <- 1:nrow(pca_sphere_df)
    pca_sphere_count_df <- select(aggregate(. ~ ind, pca_sphere_df, sum), c(ind, count))
    pca_sphere_reduced_df <- pca_sphere_df[!duplicated(pca_sphere_df$ind),]
    pca_sphere_join_df <- left_join(pca_sphere_count_df, pca_sphere_reduced_df, by = c("ind" = "ind"))
    pca_sphere_join_df$count_radius <- pca_sphere_join_df$count.x / (max(pca_sphere_join_df$count.x) / 0.15)
    return(list(pca_plot_list = pca_plot_list,
                pca_plot_sphere_df = pca_sphere_join_df[order(pca_sphere_join_df$sort_ind),]))
  })
  return(pca_list)
}
