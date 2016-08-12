##' This is a test function
##'
##' This is a test function
##' @title This is a test function
##' @param sampleKmerDistribution 
##' @param virusKmerDistribution 
##' @param edges 
##' @param colorVirus 
##' @param colorSample 
##' @return NULL
##' @author Jochen Kruppa
##' @export
kmerPyramid <- function(sampleKmerDistribution, virusKmerDistribution = NULL, edges = FALSE,
                        colorVirus = "yellow", colorSample = NULL){
    require(plyr); require(dplyr)
    require(scatterplot3d)
    ## prepare the virus kmer distribution
    if(is.null(virusKmerDistribution)){
        load("../data/virusKmerDistribution.rda")
    }
    if(!all(c("A", "C", "T", "G") %in% names(virusKmerDistribution))){
        stop("The virus kmer distibution must include a data.frame with four columns named A, C, G, and T.")
    }
    kmerVirusDistDf <- select(tbl_df(virusKmerDistribution), A, C, G, T)
    ## get the sample kmer distribution
    if(!all(c("A", "C", "T", "G") %in% names(sampleKmerDistribution))){
        stop("The sample kmer distibution must include a data.frame with four columns named A, C, G, and T.")
    }
    kmerSampleDistDf <- select(tbl_df(sampleKmerDistribution), A, C, G, T)
    ## pca
    ## put everything together
    pcaDf <- rbind(kmerVirusDistDf, kmerSampleDistDf)
    lm <- max(pcaDf) - 0.1                             # scalling parameter for the plot 
    ## edges of the pyramide
    if(edges){
        lm <- 0.8                             # scalling parameter for the plot 
        edgesMat <- data.frame(rbind(diag(1,4)))
        names(edgesMat) <- c("A", "C", "G", "T")
        pcaDf <- rbind(pcaDf, edgesMat)
        ## run pca
        pca <- predict(prcomp(pcaDf))[,1:3]
        ## build plotting data frame
        plotDf <- tbl_df(data.frame(id = c(rep("virus", nrow(kmerVirusDistDf)),
                                           rep("sample", nrow(kmerSampleDistDf)),
                                           rep("edge", 4)),
                                    pca))
        ## get the 2D edge positions
        edgeDf <- filter(plotDf, id == "edge")
        s1 <- scatterplot3d.invisible(edgeDf[1, "PC1"], edgeDf[1, "PC2"], edgeDf[1, "PC3"],
                                      xlim = c(-lm, lm), ylim = c(-lm, lm), zlim = c(-lm, lm))
        e1 <- s1$xyz.convert(edgeDf[1, "PC1"], edgeDf[1, "PC2"], edgeDf[1, "PC3"])
        s2 <- scatterplot3d.invisible(edgeDf[2, "PC1"], edgeDf[2, "PC2"], edgeDf[2, "PC3"],
                                      xlim = c(-lm, lm), ylim = c(-lm, lm), zlim = c(-lm, lm))
        e2 <- s2$xyz.convert(edgeDf[2, "PC1"], edgeDf[2, "PC2"], edgeDf[2, "PC3"])
        s3 <- scatterplot3d.invisible(edgeDf[3, "PC1"], edgeDf[3, "PC2"], edgeDf[3, "PC3"],
                                      xlim = c(-lm, lm), ylim = c(-lm, lm), zlim = c(-lm, lm))
        e3 <- s3$xyz.convert(edgeDf[3, "PC1"], edgeDf[3, "PC2"], edgeDf[3, "PC3"])
        s4 <- scatterplot3d.invisible(edgeDf[4, "PC1"], edgeDf[4, "PC2"], edgeDf[4, "PC3"],
                                      xlim = c(-lm, lm), ylim = c(-lm, lm), zlim = c(-lm, lm))
        e4 <- s4$xyz.convert(edgeDf[4, "PC1"], edgeDf[4, "PC2"], edgeDf[4, "PC3"])
    } else {
        ## run pca
        pca <- predict(prcomp(pcaDf))[,1:3]
        ## build plotting data frame
        plotDf <- data.frame(id = c(rep("virus", nrow(kmerVirusDistDf)),
                                    rep("sample", nrow(kmerSampleDistDf))),
                             pca)
    }
    ## get the sample position
    sampleDf <- filter(plotDf, id == "sample")
    samplePosList <- llply(1:nrow(sampleDf), function(i) {
        s <- scatterplot3d.invisible(sampleDf[i, "PC1"], sampleDf[i, "PC2"], sampleDf[i, "PC3"],
                                     xlim = c(-lm, lm), ylim = c(-lm, lm), zlim = c(-lm, lm))
        p <- s$xyz.convert(sampleDf[i, "PC1"], sampleDf[i, "PC2"], sampleDf[i, "PC3"])
        return(p)
    })
    ## png("virusAGTCDist1.png", width = 20, height = 20, res = 300, unit = "cm")
    p <- suppressWarnings(scatterplot3d(select(filter(plotDf, id == "virus"), PC1, PC2, PC3),
                                        color = colorVirus,
                                        cex.symbols = 1, cex.lab = 1.5, cex.axis = 1.5,
                                        mar = c(3, 3, 0, 2),
                                        xlab = "PC1", ylab = "PC2", zlab = "PC3",
                                        xlim = c(-lm, lm), ylim = c(-lm, lm), zlim = c(-lm, lm)))
    if(edges){
        segments(e1$x, e1$y, e2$x, e2$y, lwd = 2, col = 1, lty = 3)
        segments(e1$x, e1$y, e3$x, e3$y, lwd = 2, col = 1, lty = 3)
        segments(e1$x, e1$y, e4$x, e4$y, lwd = 2, col = 1, lty = 3)
        segments(e2$x, e2$y, e3$x, e3$y, lwd = 2, col = 1, lty = 3)
        segments(e2$x, e2$y, e4$x, e4$y, lwd = 2, col = 1, lty = 3)
        segments(e3$x, e3$y, e4$x, e4$y, lwd = 2, col = 1, lty = 3)
        text(e1$x, e1$y, label = "A", pos = 1, cex = 2)
        text(e2$x, e2$y, label = "T", pos = 4, cex = 2)
        text(e3$x, e3$y, label = "C", pos = 4, cex = 2)
        text(e4$x, e4$y, label = "G", pos = 3, cex = 2)
    }
    if(is.null(colorSample)){
        colSample <- rainbow(length(samplePosList))
    }
    l_ply(seq_along(samplePosList), function(i) {
        text(samplePosList[[i]]$x, samplePosList[[i]]$y,
             label = "X", pos = 1, cex = 2, col = colSample[i])
    })
    legend(p$xyz.convert(lm, lm/2, lm/2),
           col = colSample, bg="white", lty=c(1,1), lwd=2, yjust=0,
           legend = c("A", "B", "C"), cex = 1.1)
    ## dev.off()
   
}



##' Small helper function to make a invisible scatterplot3d
##'
##' Small helper function to make a invisible scatterplot3d
##' @title Small helper function to make a invisible scatterplot3d
##' @param ... 
##' @return scatterplot3d object
##' @author Jochen Kruppa
scatterplot3d.invisible <- function(...){
    ff <- tempfile()
    png(filename = ff)
    res <- suppressWarnings(scatterplot3d(...))
    dev.off()
    unlink(ff)
    res
}

## sampleKmerDistribution <- virusKmerDistribution[1:3,]
## kmerPyramid(sampleKmerDistribution, edges = TRUE)
