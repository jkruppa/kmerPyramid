##' Function to read in fastq files
##'
##' Function to read in fastq files
##' @title Function to read in fastq files
##' @param path path to the fastq file
##' @param nlines number of lines to read in 
##' @return
##' A list containing the following elements:
##' \describe{
##' \item{\emph{ids}}{The ids of the reads i.e. the names of the reads}
##' \item{\emph{reads}}{The read as character}
##' \item{\emph{phreds}}{The phreds scores per read as character}
##' }
##' @author Klaus Jung
##' @export
read.fastq = function(path, nlines=0) {
	f1 = paste(path, file, sep="")
	if (nlines==0) {
            print("Determining number of reads in file")
            nlines = countLines(f1)
	}
	nreads = nlines/4
	ids = as.character(rep(NA, nreads))
	reads = as.character(rep(NA, nreads))
	phreds = as.character(rep(NA, nreads))
	filein1 = file(f1, "r", blocking=FALSE)
	for (i in 1:nreads) {
		if (i%%round(nreads * 0.1)==0) print(paste(round(100 * i/nreads, 2), "% done"))
		readi = scan(filein1, what=character(), nlines=4, quiet=TRUE, sep="\n")
		ids[i] = readi[1]
		reads[i] = readi[2]
		phreds[i] = readi[4]
	}
	close(filein1)
	return(list(ids=ids, reads=reads, phreds=phreds))
}

##' Test title for later
##'
##' Test title for later
##' @title Test title for later
##' @param QS 
##' @param plot 
##' @param title 
##' @return NULL
##' @author Klaus Jung
##' @export
translate.phred = function(QS, plot=TRUE, title="") {
	qcode = unlist(strsplit("#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[]^_`abcdefghijklmnopqrstuvwxyz{|}~", ""))
	n = length(QS)
	QST = vector(mode="list", length=n)
	for (i in 1:n) {
		if (i%%round(n * 0.1)==0) print(paste(round(100 * i/n, 2), "% done"))
		QST[[i]] = match(unlist(strsplit(QS[i], "")), qcode)
	}
	if (plot==TRUE) {
		QST.lmax = max(lengths(QST))
		QST.distr = vector(mode="list", length=QST.lmax)
		for (j in 1:QST.lmax) {
		for (i in 1:n) {
			QST.distr[[j]][i] = QST[[i]][j]
		}}
		plot(1, 1, type="n", xlim=c(0, QST.lmax), ylim=c(0, 40), cex.lab=1.5, cex.axis=1.5, xlab="Position in read", ylab="Score", main=title)
		for (i in 1:QST.lmax) {
			B = boxplot(QST.distr[[i]], plot=FALSE, at=i, range=0, axes=FALSE)
			points(c(i, i), B$stats[c(1, 5)], type="l", col=8)
			points(c(i, i), B$stats[c(2, 4)], type="l", col=1)
			points(c(i, i), B$stats[c(3, 3)], type="l", col=2, lwd=3)
		}
		grid()
	}
	return(QST)
}
