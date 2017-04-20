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
##'   for more details and readDNAStringSet() for reading in fasta
##'   files.
##' @title ACGT pyramid 3D ploting function
##' @param df Data frame of k-mer frequencies. The single k-mers are
##'   the columns and the rows indicating different samples, sequence
##'   reads, or contigs. See teh function
##'   \code{link{get_kmer_distribution}} to generate the frequencies
##'   from a DNA sequence.
##' @param type Choose if single points [default] or lines between the
##'   points should be shown.
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
##' @param classify Sepcify the used classifier: 'kmeans' or 'hclust'
##' @param groups Number of assumed groups or 'k' [default = 2]
##' @return If classify is set a data.frame with predicated group
##'   labels.
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
##'
##' Use the build in classification
##'
##' pred_df <- pyramid_3d(kmer_distr,
##'                       cex = 1,
##'                       text = text_ids,
##'                       classify = "kmeans")
##'
##' pred_df <- pyramid_3d(kmer_distr,
##'                       cex = 1,
##'                       text = text_ids,
##'                       classify = "hclust")
pyramid_3d <- function(df,
                       type = "points",
                       color = "black",
                       ids = NULL,
                       text = NULL,
                       cex = 1,
                       identify = FALSE,
                       classify = NULL,
                       groups = 2)
{
  ## default color palette
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cbb_map <- hashmap(cbbPalette,
                     c("black", "orange", "lightblue", "green",
                       "yellow", "darkblue", "red", "violet"))
  p_load(plyr, dplyr, rgl)
  if(groups > 8) stop("More than eight groups are not supported!")
  ## get the matrix for the edges of the pyramid
  pca_plot_list <- get_pca_plot_list(df)
  ## get the classification
  if(!is.null(classify)) {
    if(!classify %in% c("hclust", "kmeans")) {
      stop("Classify must be 'hclust' or 'kmeans'")
    }
    pred <- switch(classify,
                   "kmeans" = {
                     kmeans(df, centers = groups)$cluster
                   },
                   "hclust" = {
                     cutree(hclust(dist(df), method = "average"), groups)
                   })
    ## set the new color variable
    color <- cbbPalette[as.numeric(as.factor(str_c(pred, text)))]
  }
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
    switch(type,
           "points" = {
             points3d(pca_plot_list$pca_plot_df,
                      col = color,
                      size = cex)
           },
           "lines" = {
             lines3d(pca_plot_list$pca_plot_df,
                     col = color)
           }, 
           stop("Choose type equal 'lines' or 'points'")
           )
  }
  if(identify){
    bg3d(color = c("white", "black"))
    material3d(color = "black")      
    par3d(family = "sans", cex = 1)
    identify3d(pca_plot_list$pca_plot_df,
               label = ids)
  }
  if(!is.null(classify)) {
    if(is.null(ids)) {
      return(data.frame(pred, color = cbb_map[[color]]))
    } else {
      return(data.frame(pred, resp = ids, color = cbb_map[[color]]))
    }
  }
}
