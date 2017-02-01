#' k-mer distribution of virus DNA sequences
#'
#' A dataset containing the 1, 2, 3, 4, and 5-mer
#' distribution of the nucleic acid bases of 5785 viruses available as
#' embl files from NCBI.
#'
#' @format A data frame with 5000+ rows and 1367 variables:
#' \describe{
#'   \item{id}{Name of the virus (NCBI)}
#'   \item{type}{Type of the virus}
#'   \item{species}{Species of the virus}
#'   \item{A}{Percentage of nucleic acid base A in the 1-mers}
#'   \item{...}{...}
#'   \item{TTTTT}{Percentage of nucleic acid mer TTTTT in the 5-mers}
#' }
#' @source \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC99098/}
"virusKmerDistribution"


#' k-mer distribution of bacteria DNA sequences
#'
#' A dataset containing the 1, 2, 3, 4, and 5-mer
#' distribution of the nucleic acid bases of bacteria available as
#' embl files from NCBI.
#'
#' @format A data frame with 2000+ rows and 1365 variables:
#' \describe{
#'   \item{id}{Name of the bacteria (NCBI)}
#'   \item{type}{Type of the bacteria}
#'   \item{species}{Species of the bacteria}
#'   \item{A}{Percentage of nucleic acid base A in the 1-mers}
#'   \item{...}{...}
#'   \item{TTTTT}{Percentage of nucleic acid mer TTTTT in the 5-mers}
#' }
#' @source \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC99098/}
"bacteriaKmerDistribution"

#' DNAStringSet of ten viral sequences
#'
#' A DNAStringSet containing ten randomly choosen viral sequences for
#' example purpose.
#'
#' @format A DNAStringSet with 10 entries of different width
#' \describe{
#'   \item{width}{Length of the viral sequence}
#'   \item{seq}{DNA sequence of the viral sequence}
#'   \item{names}{NCBI genbank ID and description}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/}
"viralExampleSeqs"

#' DNAStringSet of coding and non coding of NC_015230.1 viral strain
#'
#' A DNAStringSet containing coding and non coding sequences of the
#' NC_015230.1 viral strain
#'
#' @format A DNAStringSet with 41 entries of different width
#' \describe{
#'   \item{width}{Length of the viral sequence}
#'   \item{seq}{DNA sequence of the viral sequence}
#'   \item{names}{NCBI genbank ID and dsecription if CDS otherwise names non_coding}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/}
"viralExampleCodingSeq"

