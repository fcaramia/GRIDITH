\name{start.alphamax.f}
\alias{start.alphamax.f}
\alias{START.ALPHAMAX.F}
\title{
    input of starting values for the rotation function and visualization of input matrix
}
\description{
    Text or graphical input for starting values to be used in the rotation function optimizeRotation.
    Plotting of input matrix values on a transformed scale
    For graphical input vertical.cluster.line, horizontal.cluster.line, allelic.balance.cluster must be set to NULL
}
\usage{
    mystart = start.alphamax.f(matrix, maxCN, maxlines, colbychrom, xlim, ylim, dx.eq.dy,  vertical.cluster.line, horizontal.cluster.line, allelic.balance.cluster, force.diag)
}
\arguments{
    \item{matrix}{
        input matrix. see details
    }
    \item{colbychrom}{
        logical: colour minor and major CN intesity by chromosome. 
    }
    \item{maxCN}{
        integer: maximum CN value to plot 
    }
    \item{maxlines}{
        integer: maximum number of grid lines to plot
    }
    \item{xlim}{
        integer: lenght of x axis. If NULL will be calculated on basis of data
    }
    \item{ylim}
        integer: lenght of y axis. If NULL will be calculated on basis of data
    }
    \item{dx.eq.dy}{
        logical: scale axes equally. 
    }
    \item{vertical.cluster.line}{
        str: vertical line for cluster. see example 
    }
    \item{horizontal.cluster.line}{
        str: horizontal line for cluster. see example 
    }
    \item{allelic.balance.cluster}{
        str: point of allelic balance. see example 
    }
    \item{force.diag}{
        logical: force points on diagonal to stay on diagonal. 
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
    a data.frame containing input values for rotation. 
}
\author{
    Ingrid Lonnstedt
}
\seealso{
   \code{\link{plot.transformed}}
   \code{\link{optimizeRotation}}
}
\examples{
    v=list('x'=c(1.193036,1.243445),'y'=c(0.8479729,0.6519817))
    h= list('x'=c(0.9787967,1.1846342),'y'=c(0.9179698,0.8619723))
    a=list('x'=0.9493914,'y'=0.9249694)

    mystart = start.alphamax.f(mat, colbychrom = TRUE, xlim = NULL, ylim = NULL,
                           dx.eq.dy = TRUE, vertical.cluster.line=v, horizontal.cluster.line=h, allelic.balance.cluster=a,
                           force.diag = TRUE)
}
\keyword{start.alphamax.f}
