\name{plot.transformed}
\alias{plot.transformed}
\alias{PLOT.TRANSFORMED}
\title{
    Graphical representation of rotation starting values from start.alphamax.f.
}
\description{
    A longer description
}
\usage{
    plot.transformed(matrix, alpha, F, w=NULL, rhofunc = 'leastsq',
                            rhoargs = NULL, main = NULL, xlim = NULL, ylim = NULL,
                            maxCN = 24, color.by = 'eachrho')
}
\arguments{
    \item{matrix}{
        input matrix. see details
    }
    \item{alpha}{
        data.frame: output from \code{\link{start.alphamax.f}}
    }
     \item{xlim}{
        integer: lenght of x axis. If NULL will be calculated on basis of data scale
    }
    \item{ylim}
        integer: lenght of y axis. If NULL will be calculated on basis of data scale
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
    no return value
}
\author{
    Ingrid Lonnstedt
}
\seealso{
   \code{\link{start.alphamax.f}}
}
\examples{

    mystart = start.alphamax.f(mat, colbychrom = TRUE, xlim = NULL, ylim = NULL,
                           dx.eq.dy = TRUE, force.diag = TRUE)
    plot.transformed(mat, mystart$alpha, 
                 mystart$F, xlim = NULL, ylim = NULL)
}
\keyword{plot.transformed}
