# acgtPyramid

R package to visualize the acgt distribution between samples

<p align="center">
  <img src="img/kmerPlot.png" width="400">
  <img src="img/animated_pyramide.gif" width="400">
</p>


## Installation

Get the released version from CRAN:

```R
## not available yet
## install.packages("acgtPyramid")
```

Or the development version from github:

```R
# install.packages("devtools")
devtools::install_github("jkruppa/acgtPyramid")
```

## Examples

```R
load("../data/viralExampleSeqs.rda")
 
kmer_distr <- getKmerDistribution(viralExampleSeqs, k = 1)
 
pyramid_3d(kmer_distr,
           cex = 2,
           color = "blue")
 
ids <- names(viralExampleSeqs)
 
pyramid_3d(kmer_distr,
           ids = ids,
           cex = 2,
           color = "blue",
           identify = TRUE)
 
load("../data/viralExampleCodingSeq.rda")
 
kmer_distr <- getKmerDistribution(viralExampleCodingSeq, k = 1)
text_ids <- ifelse(names(viralExampleCodingSeq) == "non_coding", "x", "o")
color_ids <- ifelse(names(viralExampleCodingSeq) == "non_coding", "black", "red")
 
pyramid_3d(kmer_distr,
           cex = 1,
           text = text_ids,
           color = color_ids)
 
ids <- names(viralExampleCodingSeq)

pyramid_3d(kmer_distr,
           ids = ids,
           cex = 1,
           text = text_ids,
           color = color_ids,
           identify = TRUE)
```

```R
load("../data/viralExampleSeqs.rda")

viral_window_list <- get_pca_window_list(viralExampleSeqs, window = 2)

pyramid_3d_window(viral_window_list[1],
                  color = "red")

pyramid_3d_window(viral_window_list[1],
                  color = "red",
                  identify = TRUE)

pyramid_3d_window(viral_window_list[c(3,5)],
                  color = "red",
                  difference = TRUE)

pyramid_3d_window(viral_window_list[c(3,5)],
                  color = "red",
                  difference = TRUE,
                  identify = TRUE)
```
