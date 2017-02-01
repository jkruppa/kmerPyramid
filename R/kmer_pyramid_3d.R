##' The acgt pyramid 3D plotting function
##'
##' The acgt pyramid 3D plotting function
##' @title acgt pyramid 3D plotting function
##' @param df 
##' @param color 
##' @param ids 
##' @param text 
##' @param cex 
##' @param show.identify 
##' @return NULL
##' @author Jochen Kruppa
##' @export
##' @examples
##' data(viralExampleSeqs)
##' 
##' kmer_distr <- getKmerDistribution(viralExampleSeqs, k = 1)
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
##' kmer_distr <- getKmerDistribution(viralExampleCodingSeq, k = 1)
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

##' Wrapper for getting the positions of the edges
##'
##' Wrapper for getting the positions of the edges
##' @title Wrapper for getting the positions of the edges
##' @param pca_edge_df 
##' @return NULL
##' @author Jochen Kruppa
pyramid_3d_draw_edges <- function(pca_edge_df){
  lines3d(pca_edge_df[c(1,2), c("PC1", "PC2", "PC3")], col = 1)
  lines3d(pca_edge_df[c(1,3), c("PC1", "PC2", "PC3")], col = 1)
  lines3d(pca_edge_df[c(1,4), c("PC1", "PC2", "PC3")], col = 1)
  lines3d(pca_edge_df[c(2,3), c("PC1", "PC2", "PC3")], col = 1)
  lines3d(pca_edge_df[c(2,4), c("PC1", "PC2", "PC3")], col = 1)
  lines3d(pca_edge_df[c(3,4), c("PC1", "PC2", "PC3")], col = 1)
  text3d(pca_edge_df[1, "PC1"], pca_edge_df[1, "PC2"], pca_edge_df[1, "PC3"], "A", color="black")
  text3d(pca_edge_df[2, "PC1"], pca_edge_df[2, "PC2"], pca_edge_df[2, "PC3"], "C", color="black")
  text3d(pca_edge_df[3, "PC1"], pca_edge_df[3, "PC2"], pca_edge_df[3, "PC3"], "G", color="black")
  text3d(pca_edge_df[4, "PC1"], pca_edge_df[4, "PC2"], pca_edge_df[4, "PC3"], "T", color="black")
}

##' Calculate the PCA with pyramid edges
##'
##' Calculate the PCA with pyramid edges
##' @title Calculate the PCA with pyramid edges
##' @param df 
##' @return NULL
##' @author Jochen Kruppa
get_pca_plot_list <- function(df){
  edges_mat <- setNames(data.frame(rbind(diag(1,4))), c("A", "C", "G", "T"))
  pca_df <- rbind(edges_mat, df[c("A", "C", "G", "T")])
  pca <- predict(prcomp(pca_df))[,1:3]
  return(list(pca_plot_df = tbl_df(pca[-c(1:4),]),
              pca_edge_df = tbl_df(pca[c(1:4),])))
}

##' Get the PCA values for sliding window
##'
##' Get the PCA values for sliding window
##' @title Get the PCA values for sliding window
##' @param seqs 
##' @param window 
##' @return NULL
##' @author Jochen Kruppa
##' @export
##' @examples
##' data(viralExampleSeqs)
##' 
##' viral_window_list <- get_pca_window_list(viralExampleSeqs, window = 2)
get_pca_window_list <- function(seqs, window) {
  require(plyr); library(dplyr)
  require(stringr)
  pca_list <- llply(seqs, function(x, window = 5) {
    ## get the kmer sliding window distribution
    window_kmer_prob_df <- get_window_distr(x, window = window, prob = TRUE)
    pca_plot_list <- get_pca_plot_list(window_kmer_prob_df)
    ## get the numer of count for each k-mer
    window_kmer_count_df <- get_window_distr(x, window = window, prob = FALSE)
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

##' Wrapper for the sliding window k-mer distribution
##'
##' Wrapper for the sliding window k-mer distribution
##' @title Wrapper for the sliding window k-mer distribution
##' @param seq 
##' @param window 
##' @param prob_flag 
##' @return NULL
##' @author Jochen Kruppa
get_window_distr <- function(seq, window, prob_flag = TRUE) {
  agct_win_distr <- letterFrequencyInSlidingView(unlist(seq), view.width = window,
                                                 letters = "ACGT", OR = 0, as.prob = prob_flag)
  return(tbl_df(agct_win_distr))
}


##' The k-mer window pyramid 3D plotting function
##'
##' The k-mer window pyramid 3D plotting function
##' @title The k-mer window pyramid 3D plotting function
##' @param list 
##' @param color 
##' @param difference 
##' @param identify 
##' @return NULL
##' @author Jochen Kruppa
##' @export
##' @examples
##' data(viralExampleSeqs)
##' 
##' viral_window_list <- get_pca_window_list(viralExampleSeqs, window = 2)
##' 
##' pyramid_3d_window(viral_window_list[1],
##'                   color = "red")
##' 
##' pyramid_3d_window(viral_window_list[1],
##'                   color = "red",
##'                   identify = TRUE)
##' 
##' pyramid_3d_window(viral_window_list[c(3,5)],
##'                   difference = TRUE)
##' 
##' pyramid_3d_window(viral_window_list[c(3,5)],
##'                   difference = TRUE,
##'                   identify = TRUE)
pyramid_3d_window <- function(list,
                              color = "black",
                              difference = FALSE,
                              identify = FALSE)
{
  library(rgl)
  if(!difference) {
    if(length(list) > 1) stop("Only one sample can be plotted! Choose difference = TRUE if two samples should be compared.")
    rgl.open()
    material3d(back="culled", specular="black")
    bg3d(color = "white")
    par3d(family = "sans", cex = 2)
    ## start plotting
    pyramid_3d_draw_edges(list[[1]]$pca_plot_list$pca_edge_df)
    lines3d(list[[1]]$pca_plot_list$pca_plot_df$PC1,
            list[[1]]$pca_plot_list$pca_plot_df$PC2,
            list[[1]]$pca_plot_list$pca_plot_df$PC3,
            col = color)
    ## draw the sphere
    spheres3d(list[[1]]$pca_plot_sphere_df$PC1,
              list[[1]]$pca_plot_sphere_df$PC2,
              list[[1]]$pca_plot_sphere_df$PC3,
              alpha = 0.5,
              col = color,
              radius = list[[1]]$pca_plot_sphere_df$count_radius)
    if(identify){
      bg3d(color = c("white", "black"))
      material3d(color = "black")      
      par3d(family = "sans", cex = 1)
      identify3d(list[[1]]$pca_plot_sphere_df$PC1,
                 list[[1]]$pca_plot_sphere_df$PC2,
                 list[[1]]$pca_plot_sphere_df$PC3,
                 label = list[[1]]$pca_plot_sphere_df$ind)
    }
  } else {
    if(length(list) != 2) stop("Only two samples can be compared!")
    sphere_plot_diff_df <- get_diff_values(list)
    rgl.open()
    material3d(back="culled", specular="black")
    bg3d(color = "white")
    par3d(family = "sans", cex = 2)
    pyramid_3d_draw_edges(list[[1]]$pca_plot_list$pca_edge_df)
    spheres3d(sphere_plot_diff_df$PC1,
              sphere_plot_diff_df$PC2,
              sphere_plot_diff_df$PC3,
              alpha = 0.5,
              col = ifelse(sphere_plot_diff_df$radius_sign == 1,
                           "blue",
                           "red"),
              radius = sphere_plot_diff_df$radius_abs)
    if(identify){
      bg3d(color = c("white", "black"))
      material3d(color = "black")      
      par3d(family = "sans", cex = 1)
      identify3d(sphere_plot_diff_df$PC1,
                 sphere_plot_diff_df$PC2,
                 sphere_plot_diff_df$PC3,
                 label = sphere_plot_diff_df$ind)
    }
  }
}

##' Get the difference in counts between to sequences
##'
##' Get the difference in counts between to sequences
##' @title Get the difference in counts between to sequences
##' @param list 
##' @return NULL
##' @author Jochen Kruppa
get_diff_values <- function(list){
  vir_1 <- list[[1]]$pca_plot_sphere_df
  vir_2 <- list[[2]]$pca_plot_sphere_df
  joined <- left_join(vir_1, vir_2, by = c("ind" = "ind"))
  radius_df <- tibble(PC1 = joined$PC1.x,
                      PC2 = joined$PC2.x,
                      PC3 = joined$PC3.x,
                      radius_diff = joined$count_radius.x - joined$count_radius.y,
                      radius_abs = abs(radius_diff),
                      radius_sign = sign(radius_diff),
                      ind = joined$ind)
  return(radius_df)
}
