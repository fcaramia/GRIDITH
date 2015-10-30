################################################################
################################################################
#
# Read data and functions
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

#xlim and ylim can be omitted, but this zooming helped pick the cluster positions:

mystart = start.alphamax.f(mat, colbychrom = TRUE, xlim = c(0,2), ylim = c(0,1.5),
                           dx.eq.dy = TRUE, vertical.cluster.line=v, horizontal.cluster.line=h, allelic.balance.cluster=a,
                           force.diag = TRUE, printout = 'Plot1.png')


plot.transformed(mat, mystart, xlim = c(0,10), ylim = c(0,6), printout = 'Plot2.png')

mat = setWeights(mat, mystart, weight.constant = 5000000, nmad = 3, printout = 'Plot3.png')

#You may have to adjust weights above and rerun optimization below:
#Often it also helps to add the contraint eqfun = myeqfun
#For some reason solnp() complains about the myeqfun, but the
#results seem good anyway.
tmp = optimizeGrid(mat, mystart, UB = c(1, Inf, 0, 0, Inf), 
                        eqfun = NULL, printout = 'Plot4.png')
mat = tmp$mat
pars = tmp$pars

#Fix cn.cut, cn.cut0, cn.cuttop until regions outside 
#the regular pattern of the main subclone are all grey:
mat = setCutsUnscaled(mat, a.cut0 = 3, a.cuttop = 6.5, a.cut = 7,
                      xlim = c(0,10), ylim = c(0,8), printout = 'Plot5.png')

#Choose cut so that red segments (bottom two panels) seem
#safely distant from the blue segments
mat = subclonalCNdist(mat, cut = 3, printout = 'Plot6.png')

mat = setType(mat, xlim = c(0,21), legend = TRUE, printout = 'Plot7.png')


################################################################
################################################################
################################################################
################################################################
