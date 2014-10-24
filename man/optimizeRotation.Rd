\name{optimizeRotation}
\alias{optimizeRotation}
\alias{OPTIMIZEROTATION}
\title{
    normalize data for BAF bias (skewness)
}
\description{
    normalize SNP array segment CNs for BAF bias (seen as skewness in grid plots). It will optimize the fit of a rectangular lattice pattern to a grid plot by minimizing the summed distance from each grid point to its closest lattice point through Tukey's robust function for M-estimation. Note that the function can only find a locally optimized rotation. The starting values are important to find the best global fit.
}
\usage{
    setWeights(matrix, mystart, weight.constant = 5000000, nmad = 3)
}
\arguments{
    \item{matrix}{
        input matrix. see details
    }
    \item{LB}{
        (Optional) The lower bound on the parameters.
    }
    \item{UB}{
        (Optional) The upper bound on the parameters.
    }
    \item{eqfun}{
        (Optional) The equality constraint function returning the vector of evaluated equality constraints.
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
    data.frame = modified input matrix, optimized values
}
\author{
    Ingrid Lonnstedt
}
\seealso{
   \code{\link{setWeights}}
   \code{\link{start.alphamax.f}}
}
\examples{
    
    #LOAD DATA

    #EXAMPLE CLUSTER
    v=list('x'=c(1.193036,1.243445),'y'=c(0.8479729,0.6519817))
    h= list('x'=c(0.9787967,1.1846342),'y'=c(0.9179698,0.8619723))
    a=list('x'=0.9493914,'y'=0.9249694)

    mystart = start.alphamax.f(mat, colbychrom = TRUE, xlim = NULL, ylim = NULL, dx.eq.dy = TRUE, vertical.cluster.line=v, horizontal.cluster.line=h, allelic.balance.cluster=a, force.diag = TRUE)

    mat = optimizeRotation(mat, mystart, UB = c(1, Inf, 0, 0, Inf), eqfun = NULL)
}
\keyword{optimizeRotation}
\keyword{skewness}
