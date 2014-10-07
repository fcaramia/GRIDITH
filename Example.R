################################################################
################################################################
#
# Read example data
#
################################################################
################################################################

mat <- readRDS( 'segfile.RData')


################################################################
################################################################
#
# Rotation and ITH
#
################################################################
################################################################

#We clicked these positions (x, y):
#Two points along a vertical line: 
#  (1.193036,  0.8479729) and (1.243445, 0.6519817)
#Two points along a horizontal line: 
#  (0.9787967, 0.9179698) and (1.1846342, 0.8619723)
#An allelic balance cluster:
#  (0.9493914, 0.9249694)
#The algorithm only needs the clicks to be approximate.
#Make sure the preliminary rotation looks good before proceeding.
v=list('x'=c(1.193036,1.243445),'y'=c(0.8479729,0.6519817))
h= list('x'=c(0.9787967,1.1846342),'y'=c(0.9179698,0.8619723))
a=list('x'=0.9493914,'y'=0.9249694)

mystart = start.alphamax.f(mat, colbychrom = TRUE, xlim = NULL, ylim = NULL,
                           dx.eq.dy = TRUE, vertical.cluster.line=v, horizontal.cluster.line=h, allelic.balance.cluster=a,
                           force.diag = TRUE)


plot.transformed(mat, mystart, xlim = NULL, ylim = NULL)

mat = setWeights(mat, mystart, weight.constant = 5000000, nmad = 3)

#You may have to adjust weights above and rerun optimization below:
#Often it also helps to add the contraint eqfun = myeqfun
#For some reason solnp() complains about the myeqfun, but the
#results seem good anyway.
mat = optimizeRotation(mat, mystart, UB = c(1, Inf, 0, 0, Inf), 
                        eqfun = NULL)
#Fix cn.cut, cn.cut0, cn.cuttop until regions outside 
#the regular pattern of the main subclone are all grey:
mat = setCutsUnscaled(mat, a.cut0 = 3, a.cuttop = 6.5, a.cut = 7,
                      xlim = c(0,10), ylim = c(0,10))

#Choose cut so that red segments (bottom two panels) seem
#safely distant from the blue segments
mat = subclonalCNdist(mat, cut = 3)

mat = setType(mat, xlim = c(0,21), legend = TRUE)


################################################################
################################################################
################################################################
################################################################
