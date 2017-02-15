##' Wrapper for the sliding window k-mer distribution
##'
##' Wrapper for the sliding window k-mer distribution
##' @title Wrapper for the sliding window k-mer distribution
##' @param seq DNAStringSet object 
##' @param window Size of the window
##' @param prob_flag Should the frequencies of the k-mers returned or the counts
##' @return Data frame of the k-mer distribution. 
##' @author Jochen Kruppa
get_window_distr <- function(seq, window, prob_flag = TRUE) {
  agct_win_distr <- letterFrequencyInSlidingView(unlist(seq), view.width = window,
                                                 letters = "ACGT", OR = 0, as.prob = prob_flag)
  return(tbl_df(agct_win_distr))
}
