##' The acgt pyramid 3D ploting function allows to plot the ACGT
##' distribution of a given sample into the 3D space by a principal
##' component analysis.
##'
##' The function allows to plot the ACGT distribution of a given
##' sample into the 3D space using a PCA. The PCA is done by the
##' function \emph{prcomp} with the default parameters. Further the
##' function is able to draw the points in different colors and text
##' symbols. Only points and letters are available. The identify
##' option allows to identify single points in the 3D plot.
##'
##' The sequences must be provided as a DNAStringSet object. See
##'   \url{https://web.stanford.edu/class/bios221/labs/biostrings/lab_1_biostrings.html}
##'   for more details and readDNAStringSet() for reading in fasta files.
##' @title ACGT pyramid 3D ploting function
##' @param df Data frame of k-mer frequencies. The single k-mers are
##'   the columns and the rows indicating different samples, sequence
##'   reads, or contigs. See teh function
##'   \code{link{get_kmer_distribution}} to generate the frequencies
##'   from a DNA sequence.
##' @param color Single value 'black', if all points should be black
##'   or vector of length \emph{nrow(df)} if all points should be
##'   colored differently.
##' @param ids If \emph{identify = TRUE} the single points can be
##'   selected and the given id is shown. Must have the length
##'   \emph{nrow(df)}.
##' @param text Single value 'x', if all points should be printed as
##'   'x' or vector of length \emph{nrow(df)} if all points should be
##'   printed by a different letter differently.
##' @param cex Size of the shown text.
##' @param identify Set to TRUE, if points should be identified by
##'   their \emph{ids}.
##' @return NULL
##' @author Jochen Kruppa
##' @export
##' @examples
##' ## Read in own DNA sequences by the package Biostrings (see Details for more information)
##' 
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
