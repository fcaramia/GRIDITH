\name{setCutsUnscaled}
\alias{setCutsUnscaled}
\alias{SETCUTSUNSCALED}
\title{
    set grid limits
}
\description{
    set the upper, lower and right side limit of the grid plot rectangle within which a regular grid pattern is evident. The grid plot refers to a plot of minor versus major segment CN plot.
}
\usage{
    setCutsUnscaled(mat, a.cut0 = 3, a.cuttop = 6.5, a.cut = 7, xlim = c(0,10), ylim = c(0,10))
}
\arguments{
    \item{matrix}{
        input matrix. see details
    }
    \item{a.cut0}{
        
    }
    \item{a.cuttop}{
        
    }
    \item{a.cut}{
        
    }
    \item{xlim}{
        (Optional) Scale of x axis
    }
    \item{ylim}{
        (Optional) Scale of y axis
    }
}
\details{
    matrix:
        The ouptput matrix from \code{\link{optimizeRotation}}
}
\value{
    data.frame = modified input matrix, optimized values
}
\author{
    Ingrid Lonnstedt
}
\seealso{
   \code{\link{setWeights}}
   \code{\link{start.alphamax.f}}
   \code{\link{optimizeRotation}}
}
\examples{
    
    #LOAD DATA

    #EXAMPLE CLUSTER
    v=list('x'=c(1.193036,1.243445),'y'=c(0.8479729,0.6519817))
    h= list('x'=c(0.9787967,1.1846342),'y'=c(0.9179698,0.8619723))
    a=list('x'=0.9493914,'y'=0.9249694)

    mystart = start.alphamax.f(mat, colbychrom = TRUE, xlim = NULL, ylim = NULL, dx.eq.dy = TRUE, vertical.cluster.line=v, horizontal.cluster.line=h, allelic.balance.cluster=a, force.diag = TRUE)

    mat = optimizeRotation(mat, mystart, UB = c(1, Inf, 0, 0, Inf), eqfun = NULL)

    mat = setCutsUnscaled(mat, a.cut0 = 3, a.cuttop = 6.5, a.cut = 7, xlim = c(0,10), ylim = c(0,10))
}
\keyword{optimizeRotation}
\keyword{skewness}