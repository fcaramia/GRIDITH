\name{setWeights}
\alias{setWeights}
\alias{SETWEIGHTS}
\title{
    visualizes potential parameter settings for optimizeRotation
}
\description{
    visualizes potential parameter settings for \code{\link{optimizeRotation}}.
}
\usage{
    setWeights(matrix, mystart, weight.constant = 5000000, nmad = 3)
}
\arguments{
    \item{matrix}{
        input matrix. see details
    }
    \item{mystart}{
        data.frame: output from \code{\link{start.alphamax.f}}
    }
}
\details{
    matrix:
        a data.frame consisting of 7 columns:

        Chromosome: integer - Chromosome where segment is located

        Start.bp : interger - starting position of segment

        End.bp : integer - end position of segment

        a1: real - minor copy number estimate on arbitrary scale

        a2: real - major copy number estimate on arbitrary scale

        W: real - length of segment divided by the total genome length

        length: integer - length of segment  
}
\value{
    data.frame = modified input matrix, added rotation weights for \code{\link{optimizeRotation}}
}
\author{
    Ingrid Lonnstedt
}
\seealso{
   \code{\link{optimizeRotation}}
   \code{\link{start.alphamax.f}}
}
\examples{
    
    mystart = start.alphamax.f(matrix, colbychrom = TRUE, xlim = NULL, ylim = NULL,
                           dx.eq.dy = TRUE, force.diag = TRUE)

    mat = setWeights(matrix, mystart, weight.constant = 5000000, nmad = 3)

}
\keyword{setWeights}
