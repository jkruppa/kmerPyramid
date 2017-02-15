##' The acgt pyramid 3D plotting function
##'
##' The acgt pyramid 3D plotting function
##' @title ACGT pyramid 3D plotting function
##' @param df Data frame of k-mer freuencies. The single k-mers are
##'   the columns and the rows indicating different samples, sequence
##'   reads, or contigs.
##' @param color Single value 'black', if all points should be black
##'   or vector of length \emph{nrow(df)} if all points should be
##'   colored differently.
##' @param ids If \emph{identify = TRUE} the single points can be
##'   selected and the given id is shown. Must have the length
##'   \emph{nrow(df)}.
##' @param text Single value 'x', if all points should be printed as
##'   'x' or vector of length \emph{nrow(df)} if all points should be
##'   printed by a different letter differently.
##' @param cex Size of the text.
##' @param identify Set to TRUE, if points should be identified by
##'   their \emph{ids}.
##' @return NULL
##' @author Jochen Kruppa
##' @export
##' @examples
##' data(viralExampleSeqs)
##' 
##' kmer_distr <- get_kmer_distribution(viralExampleSeqs, k = 1)
##' 
##' pyramid_3d(kmer_distr,
##'            cex = 2,
##'            color = "blue")
##' 
##' ids <- names(viralExampleSeqs)
##' 
##' pyramid_3d(kmer_distr,
##'            ids = ids,
##'            cex = 2,
##'            color = "blue",
##'            identify = TRUE)
##' 
##' data(viralExampleCodingSeq)
##' 
##' kmer_distr <- get_kmer_distribution(viralExampleCodingSeq, k = 1)
##' text_ids <- ifelse(names(viralExampleCodingSeq) == "non_coding", "x", "o")
##' color_ids <- ifelse(names(viralExampleCodingSeq) == "non_coding", "black", "red")
##' 
##' pyramid_3d(kmer_distr,
##'            cex = 1,
##'            text = text_ids,
##'            color = color_ids)
##' 
##' ids <- names(viralExampleCodingSeq)
##' 
##' pyramid_3d(kmer_distr,
##'            ids = ids,
##'            cex = 1,
##'            text = text_ids,
##'            color = color_ids,
##'            identify = TRUE)
pyramid_3d <- function(df,
                       color = "black",
                       ids = NULL,
                       text = NULL,
                       cex = 1,
                       identify = FALSE)
{
  require(plyr); require(dplyr)
  require(rgl)
  ## get the matrix for the edges of the pyramid
  pca_plot_list <- get_pca_plot_list(df)
  ## rgl parameter
  rgl.open()
  material3d(back="culled", specular="black")
  bg3d(color = "white")
  par3d(family = "sans", cex = 2)
  ## start with the edges
  pyramid_3d_draw_edges(pca_plot_list$pca_edge_df)
  ## start plotting the points
  if(!is.null(text)) {
    texts3d(pca_plot_list$pca_plot_df,
            texts = text,
            cex = cex,
            col = color)
  } else {
    points3d(pca_plot_list$pca_plot_df,
             col = color,
             size = cex)
  }
  if(identify){
    bg3d(color = c("white", "black"))
    material3d(color = "black")      
    par3d(family = "sans", cex = 1)
    identify3d(pca_plot_list$pca_plot_df,
               label = ids)
  }
}
