##' This is a test function
##'
##' This is a test function
##' @title This is a test function
##' @param foo 
##' @param bar 
##' @return NULL
##' @author Jochen Kruppa
##' @export
kmerPyramid <- function(sampleKmerDistribution, virusKmerDistribution = NULL){
    if(is.null(virusKmerDistribution)){
        load("../data/virusKmerDistribution.rda")
    }
    
    
    print(foo)
}
