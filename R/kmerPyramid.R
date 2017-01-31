##' This is a test function
##'
##' This is a test function
##' @title This is a test function
##' @param acgtDistributionList distribution of the acgt
##' @param colorList list of color values
##' @param ids names of samples or species
##' @param show.edges should the deges in the scatter plot shown? (default = FALSE) 
##' @param method should a scatterplot (scatter) or 3d plot (3d) be plotted (default = scatter)
##' @return NULL
##' @author Jochen Kruppa
##' @export
acgtPyramid <- function(acgtDistributionList,
                        colorList,
                        ids = NULL,
                        show.edges = FALSE,
                        text = NULL,
                        text.cex = 1,
                        show.identify = FALSE,
                        method = "scatter"){
    require(plyr); require(dplyr)
    require(scatterplot3d)
    require(rgl)
    kmerDistributionList <- acgtDistributionList ## rename befor finishing
    ## only select the ACTG information
    kmerPcaList <- llply(kmerDistributionList, function(x) {
        return(select(tbl_df(x), A, C, G, T))
    })
    ## pca
    pcaDf <- do.call(rbind, kmerPcaList)
    lm <- max(pcaDf) - 0.2                            # scalling parameter for the plot 
    ## edges of the pyramide
    idsList <- llply(names(kmerPcaList), function(x) rep(x, nrow(kmerDistributionList[[x]])))
    if(show.edges){
        lm <- 0.8                             # scalling parameter for the plot 
        edgesMat <- data.frame(rbind(diag(1,4)))
        names(edgesMat) <- c("A", "C", "G", "T")
        pcaDf <- rbind(pcaDf, edgesMat)
        ## run pca
        pca <- predict(prcomp(pcaDf))[,1:3]
        ## build plotting data frame
        plotDf <- tbl_df(data.frame(id = c(do.call(c, idsList),
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
        plotDf <- data.frame(id = do.call(c, idsList), pca)
    }
    if("sample" %in% names(kmerDistributionList)){
        ## get the sample position
        sampleDf <- filter(plotDf, id == "sample")
        samplePosList <- llply(1:nrow(sampleDf), function(i) {
            s <- scatterplot3d.invisible(sampleDf[i, "PC1"], sampleDf[i, "PC2"], sampleDf[i, "PC3"],
                                         xlim = c(-lm, lm), ylim = c(-lm, lm), zlim = c(-lm, lm))
            p <- s$xyz.convert(sampleDf[i, "PC1"], sampleDf[i, "PC2"], sampleDf[i, "PC3"])
            return(p)
        })
    }
    switch(method,
           "scatter" = {
               if("species" %in% names(kmerDistributionList)){
                   p <- suppressWarnings(scatterplot3d(select(filter(plotDf, id == "species"),
                                                              PC1, PC2, PC3),
                                                       color = colorList[["species"]],
                                                       cex.symbols = 1,
                                                       cex.lab = 1,
                                                       cex.axis = 1,
                                                       mar = c(3, 3, 0, 2),
                                                       pch = 1,
                                                       box = FALSE,
                                                       xlab = "PC1", ylab = "PC2", zlab = "PC3",
                                                       xlim = c(-lm, lm),
                                                       ylim = c(-lm, lm),
                                                       zlim = c(-lm, lm)))
               } else {
                   p <- suppressWarnings(scatterplot3d(select(filter(plotDf, id == "sample"),
                                                              PC1, PC2, PC3),
                                                       color = colorList[["sample"]],
                                                       cex.symbols = 1,
                                                       cex.lab = 1,
                                                       cex.axis = 1,
                                                       mar = c(3, 3, 0, 2),
                                                       pch = 1,
                                                       box = FALSE,
                                                       xlab = "PC1", ylab = "PC2", zlab = "PC3",
                                                       xlim = c(-lm, lm),
                                                       ylim = c(-lm, lm),
                                                       zlim = c(-lm, lm)))      
               }
           },
           "3d" = {
             rgl.open()
             bg3d(color = "white")
             par3d(family = "sans", cex = 2)
             if("species" %in% names(kmerDistributionList)){
               if(!is.null(text)) {
                 texts3d(select(filter(plotDf, id == "species"), PC1, PC2, PC3),
                         texts = text, cex = text.cex, col = "black")
               } else {
                 points3d(select(filter(plotDf, id == "species"), PC1, PC2, PC3),
                          col = colorList[["species"]])
               }
               if(show.identify){
                 identify3d(select(filter(plotDf, id == "species"), PC1, PC2, PC3),
                            label = ids)
               }
             } else {
               if(!is.null(text)) {
                 texts3d(select(filter(plotDf, id == "sample"), PC1, PC2, PC3),
                         texts = text, cex = text.cex, col = "black")
               } else {
                 points3d(select(filter(plotDf, id == "sample"), PC1, PC2, PC3),
                          col = colorList[["sample"]])
               }
               if(show.identify) {
                 identify3d(select(filter(plotDf, id == "sample"), PC1, PC2, PC3),
                            label = ids)
               }
             }
           })
    if(show.edges & method != "3d"){
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
    if(show.edges & method == "3d"){
      lines3d(edgeDf[c(1,2), c("PC1", "PC2", "PC3")], col = 1)
      lines3d(edgeDf[c(1,3), c("PC1", "PC2", "PC3")], col = 1)
      lines3d(edgeDf[c(1,4), c("PC1", "PC2", "PC3")], col = 1)
      lines3d(edgeDf[c(2,3), c("PC1", "PC2", "PC3")], col = 1)
      lines3d(edgeDf[c(2,4), c("PC1", "PC2", "PC3")], col = 1)
      lines3d(edgeDf[c(3,4), c("PC1", "PC2", "PC3")], col = 1)
      text3d(edgeDf[1, "PC1"], edgeDf[1, "PC2"], edgeDf[1, "PC3"], "A", color="black")
      text3d(edgeDf[2, "PC1"], edgeDf[2, "PC2"], edgeDf[2, "PC3"], "C", color="black")
      text3d(edgeDf[3, "PC1"], edgeDf[3, "PC2"], edgeDf[3, "PC3"], "G", color="black")
      text3d(edgeDf[4, "PC1"], edgeDf[4, "PC2"], edgeDf[4, "PC3"], "T", color="black")
    }
    if("sample" %in% names(kmerDistributionList) & method != "3d"){
      l_ply(seq_along(samplePosList), function(i) {
        text(samplePosList[[i]]$x, samplePosList[[i]]$y,
             label = "X", pos = 1, cex = 2, col = colorList[["sample"]][i])
      })
      legend("topleft", ##inset = 1,
             ##horiz = TRUE, xpd = TRUE,
             col = colorList[["sample"]], bg="white", lty=c(1,1), lwd=2, yjust=0.5,
             legend = ids, cex = 1.1)
      ## dev.off()
    } else {
      if(!is.null(ids) & method != "3d"){
        legend("topleft", ##inset = 1,
               ##horiz = TRUE, xpd = TRUE,
               col = unique(colorList[["species"]]), bg="white", lty=c(1,1), lwd=2, yjust=0.5,
               legend = unique(ids), cex = 1.1)
      }
    }
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

