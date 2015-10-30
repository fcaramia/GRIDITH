
################################################################
################################################################
# Script: Gridith beta version R code with example
# Author: Ingrid Lonnstedt
# Date:  25/06/2014
# R version: R 3.0.2
# 
################################################################
################################################################

################################################################
################################################################
#
# Functions to load
#
################################################################
################################################################

setType = function(mat, xlim = NULL, ylim = NULL, legend = TRUE,
                   printout = NULL, width = 500, height = 500){
  mat$type = 'A'
  mat$type[mat$D>mat$cut[1]] = 'B'
  mat$type[mat$a2unscaled>mat$a.cut[1] ] = 'D'
  mat$type[mat$a1unscaled<mat$a.cut0[1]] = 'C'
  mat$type[mat$a1unscaled>mat$a.cuttop[1]] = 'D'
  
  par(mfrow = c(1,1), mar = c(4, 4, 2, .1))
  if (is.null(xlim)) xlim = c(0, max(mat$a2unscaled))
  if (is.null(ylim)) ylim = c(0, max(mat$a1unscaled))
  cols = rep(col.transp('blue', .7), nrow(mat))
  cols[mat$type == 'D'] = col.transp('red', .5)
  cols[mat$type == 'B'] = col.transp('magenta', .5)
  cols[mat$type == 'C'] = col.transp('green', .5)
  plot(mat$a2unscaled, mat$a1unscaled, pch = 16, 
       xlim = xlim, ylim = ylim, xlab = 'Major array intensity, unscaled',
       ylab = 'Minor array intensity, unscaled', col = cols)
  abline(0,1)
  
  #Small grid
  maxCN = 24
  alpha = mat$alphamax[1]
  EEs = c(alpha*c(0:maxCN) + (1-alpha))
  for (j in 1:length(EEs)){
    lines(range(EEs), rep(EEs[j], 2), lty = 'solid', col = 'grey20', lwd = .5)
    lines(rep(EEs[j], 2), range(EEs), lty = 'solid', col = 'grey20', lwd = .5)
  }

  if (legend)
  {
    mtext(paste('A:', round(sum(mat$W[mat$type == 'A'])*100), '% of genome'), adj = 0, line = 1)
    mtext(paste('B:', round(sum(mat$W[mat$type == 'B'])*100), '% of genome'), adj = 0)
    mtext(paste('C:', round(sum(mat$W[mat$type == 'C'])*100), '% of genome'), adj = 0.5, line = 1)
    mtext(paste('D:', round(sum(mat$W[mat$type == 'D'])*100), '% of genome'), adj = 0.5)
    mat$ith = sum(mat$W[mat$type == 'B'])/(sum(mat$W[mat$type == 'A']) + 
      sum(mat$W[mat$type == 'B']))
    mtext(paste('ITH:', round(mat$ith[1]*100), '%'), adj = 1)
    legend('bottomright', legend = c('A','B','C','D'),
           pch = 16, col = c('blue','magenta','green','red'), bty = 'n')


  }
  print (paste('ith estimate: ',mat$ith[1]))
  
  if(!is.null(printout)) dev.print(png, file = printout, width= width, height = height)
  
    mat
}

subclonalCNdist = function(mat, cut = 3, w = 1-exp(-mat$length/500000),
                           printout = NULL, width = 700, height = 700){
  mymat = mat
  ix = which(mat$a1unscaled<= mat$a.cuttop[1] & mat$a1unscaled>=mat$a.cut0[1] & mat$a2unscaled<=mat$a.cut[1])
  mat = mat[ix,]
  alpha = mat$alphamax[1]
  mat$c1 = (mat$a1unscaled - (1 - alpha))/alpha
  mat$c2 = (mat$a2unscaled - (1 - alpha))/alpha
  if (length(w)==1) w = rep(w, nrow(mat)) 
  if (length(w)>nrow(mat)) w = w[ix]
  
  #Overlay lattice points into X
  maxCN = 24
  E2 = 0:maxCN
  E1 = matrix(E2, ncol = maxCN + 1, nrow = maxCN + 1)
  E2 = t(E1)
  E1 = c(E1)
  E2 = c(E2)
  
  x = cbind(mat$c1, mat$c2)
  alldist = apply(x, 1, dist, y = cbind(E1, E2))
  closest = apply(alldist, 2, which.min)
  e = cbind(E1, E2)[closest,]
  di = dist(x, e)
  xcoord = x[,2]-e[,2]
  ycoord = x[,1]-e[,1]
  X = cbind(ycoord, xcoord)
  
  #Colour scale for figures
  pal = brewer.pal(11, 'Spectral')
  pal = rev(pal[c(1:4, 9:11)])
  
  par(mfrow = c(2,2), mar = c(4, 4, 3, .1))
  #Multivariate t robust covariance qq=plots
  nu = 2
  D = mahalanobis.w(X, nu, w)  
  cols = (colors.by(D, pal, span = c(0,10))$mycols)[order(D)]
  qqplot(qexp(ppoints(length(D)), rate = .5), D, col = cols,
         main = 'Exponential QQ-plot',ylim = c(0,20),
         ylab = 'Mahalanobis squared distance',
         xlab = 'Exponential quantiles')
  abline(0, 1, col = 'gray')
  abline(h = cut, lty = 'dotted')
  #Overlay lattice points
  col1 = colors.by(D, pal, span = c(0,10))
  plot(X[,2], X[,1], xlim = c(-.75, .75), col = col1$mycols, pch = 16,
       ylim = c(-.75, .75), main = '', cex.main = 1,
       xlab = 'Distance to closest major integer CN',
       ylab = 'Distance to closest minor integer CN')
  mtext('Colouring by squared Mahalanobis distance to lattice point', cex = .8)
  abline(v=0, h=0)
  rect(xleft = -.5, ybottom = -.5, xright = .5, 
       ytop = .5, lty = 'dotted')
  #Colour scale
  tmp = c(0.1, .75)
  tmp = seq(tmp[1], tmp[2], length=length(col1$cols))
  points(tmp, rep(-.75, length(tmp)), col = col1$cols, pch = 16)
  
  #Subclonal CN segments given cutoff with 1-pexp(cut, .5) sign.level
  cols = rep('blue', nrow(mat))
  cols[D > cut] = 'red'
  cols = col.transp(cols)
  plot(X[,2], X[,1], xlim = c(-.75, .75), col = cols, pch = 16,
       ylim = c(-.75, .75), main = paste('Subclonal CNs with cutoff =', 
                                         cut), cex.main = 1,
       xlab = 'Distance to closest major integer CN',
       ylab = 'Distance to closest minor integer CN')
  abline(v=0, h=0)
  rect(xleft = -.5, ybottom = -.5, xright = .5, 
       ytop = .5, lty = 'dotted')
  legend('topleft', bty = 'n', pch = 16, col = c('blue','red'), 
         legend = c('Segment with integer CN','Segment with subclonal CN'))
  
  #Choose how extreme CNs can be judged to be subclonal: CN subclonality defined here
  xlim = c(0, quantile(mat$c2, probs = c(.9))*c(1.2))
  ylim = c(0, quantile(mat$c1, probs = c(.95))*c(1.2))
  plot.cn(x = mat$c2, y = mat$c1, mycols = cols, xlim = xlim, ylim = ylim,
          xlab = 'Major CN, unscaled', ylab = 'Minor CN, unscaled')
  mymat$D = NA
  mymat$D[ix] = D
  mymat$cut = cut
  
  if(!is.null(printout)) dev.print(png, file = printout, width= width, height = height)
  
  mymat
}
mahalanobis.w = function(X, nu, w){
  require(MASS)
  c = cov.trob(X, center = c(0,0), nu = nu)$cov
  D = numeric(nrow(X))
  for (j in 1:nrow(X)){
    D[j] = mahalanobis(X[j,], center = c(0,0), cov = c/w[j])
  } 
  D
}
colors.by = function(x, pal = brewer.pal(9,'YlOrRd'), span = range(na.omit(x)),
                     alpha = .5){
  
  cols = colorRampPalette(pal)(100)
  colbreaks = seq(min(span), max(span), length = length(cols))
  delta = colbreaks[2] - colbreaks[1]
  colind = pmax((x-colbreaks[1]) %/% delta + 1, 1)
  colind[colind>length(cols)] = length(cols)
  mycols = col.transp(cols[colind], alpha = alpha)
  list(mycols=mycols, cols = cols)
}


setCutsUnscaled = function(mat, a.cut= 2.2, a.cut0=1, a.cuttop=2.3,
                           xlim = NULL, ylim = NULL, printout = NULL,
                           width = 500, height = 500){
  
  alpha = mat$alphamax[1]
  maxCN = 24
  if (is.null(xlim)) xlim = c(0, max(mat$a2unscaled))
  if (is.null(ylim)) ylim = c(0, max(mat$a1unscaled))
  par(mfrow = c(1,1))
  plot(mat$a2unscaled, mat$a1unscaled, col = col.transp('blue', .7), pch = 16, 
       xlim = xlim, ylim = ylim, xlab = 'Major array intensity',
       ylab = 'Minor array intensity')
  ix = mat$a2unscaled>a.cut
  points(mat$a2unscaled[ix], mat$a1unscaled[ix], pch = 16, col = 'grey')
  ix = mat$a1unscaled<a.cut0
  points(mat$a2unscaled[ix], mat$a1unscaled[ix], pch = 16, col = 'grey')
  ix = mat$a1unscaled>a.cuttop
  points(mat$a2unscaled[ix], mat$a1unscaled[ix], pch = 16, col = 'grey')
  
  abline(0,1)
  #Small grid
  EEs = c(alpha*c(0:maxCN) + (1-alpha))
  for (j in 1:length(EEs)){
    lines(range(EEs), rep(EEs[j], 2), lty = 'solid', col = 'grey20', lwd = .5)
    lines(rep(EEs[j], 2), range(EEs), lty = 'solid', col = 'grey20', lwd = .5)
  }
  mat$a.cut = a.cut
  mat$a.cut0 = a.cut0
  mat$a.cuttop = a.cuttop

  if(!is.null(printout)) dev.print(png, file = printout, width= width, height = height)
  mat
 }

getdist = function(mat, alpha, F, maxCN = 24){
  x = cbind(mat$a1, mat$a2)
  Fm1 = solve(F)  
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  
  E1 = c(alpha*A1 + (1-alpha))
  E2 = c(alpha*A2 + (1-alpha))
  
  ### Positions of lattice points in observed scale and
  ### each segment's distance to it's lattice point -> weight
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  alldist = apply(x, 1, dist, y = cbind(e1, e2))
  closest = apply(alldist, 2, which.min)
  e = cbind(e1, e2)[closest,]
  
  di = dist(x, e)
  di
}
tukey = function(x, args, ...){
  k = args$c
  out =  k*k/6*(1-(1-(x/k)^2)^3)
  ix = which(abs(x)>k)
  out[ix] = k*k/6
  out
}
dist = function(x, y){
  if (is.null(dim(x))) x = matrix(x, ncol = 2)
  if (is.null(dim(y))) y = matrix(y, ncol = 2)
  sqrt((x[,1] - y[,1])**2 + (x[,2] - y[,2])**2)
}

col.transp <- function(col1, alpha=0.5){
  out = NULL
  for (j in 1:length(col1)){
    tmp=as.list(t(col2rgb(col1[j])))
    tmp$alpha=alpha*255
    tmp$maxColorValue=255
    out = c(out, do.call('rgb',tmp))
  }
  out
}

optimizeGrid = function(mat, mystart, UB = NULL, LB = NULL,
                            eqfun = NULL, printout = NULL, width = 640, height = 500){
  di = getdist(mat, mystart$alphamax, mystart$F) 
  mat$rhoargs = mad(di)*mat$nmad[1]
  x = cbind(mat$a1, mat$a2)
  
  require(Rsolnp)
  pars = solnp(pars = c(mystart$alphamax, mystart$F), fun = rho.alphamax.f, x = x,
               w = mat$rotationw, maxCN = 24, rhofunc = mat$rhofunc[1], 
               eqfun = eqfun, eqB = 0,
               rhoargs = list(c=mat$rhoargs), LB = LB, UB = UB)$pars
  
  alpha = pars[1]
  F = matrix(pars[-1], ncol = 2, nrow = 2)
  Fm1 = solve(F)
  a1p = mat$a1*Fm1[1,1] + mat$a2 * Fm1[1,2]
  a2p = mat$a1*Fm1[2,1] + mat$a2 * Fm1[2,2]
  xlim = c(0, quantile(a2p, probs = .95, na.rm = T)*1.2)
  ylim = c(0, quantile(a1p, probs = .95, na.rm = T)*1.2)
  mystart$F = F
  mystart$alpha = alpha
  plot.transformed(mat, mystart, w=mat$rotationw, 
                   rhoargs=list(c=mat$rhoargs[1]), 
                   rhofunc = mat$rhofunc[1], 
                   xlim = xlim, ylim = ylim)
  
  #All possible alphas and Fs
  maxCN = 24
  
  #All possible alphas
  EEs = alpha*c(0:maxCN) + (1-alpha) 
  alphas = NULL
  fs = NULL
  for (j in 1:length(EEs)){
    alphas = c(alphas, alpha/(EEs[j]+alpha))
    if (j == 1) fs = 1 else {
      fs = c(fs, fs[j-1]*alphas[j-1]/alphas[j])
    }
  }
  pars = list(alphas = alphas, fs = fs, rotation.pars = pars)

  mat$alphamax = pars$rotation.pars[1]
  F = matrix(pars$rotation.pars[-1], ncol = 2, nrow = 2)
  F = F*pars$fs[1]

  Fm1 = solve(F)
  a1r = mat$a1*Fm1[1,1] + mat$a2 * Fm1[1,2]
  a2r = mat$a1*Fm1[2,1] + mat$a2 * Fm1[2,2]
  mat$a1unscaled = pmin(a1r, a2r)
  mat$a2unscaled = pmax(a1r, a2r)
  
  if(!is.null(printout)) dev.print(png, file = printout, width= width, height = height)
  
  list(mat = mat, pars = pars)
}


myeqfun = function(pars, x, w, maxCN=24, rhofunc = 'tukey', 
                   rhoargs = NULL){
  F = matrix(pars[-1], ncol = 2, nrow = 2)
  F[1,1]+F[2,1] - F[1,2] - F[2,2]
}

setWeights = function(mat, mystart, weight.constant = 500000, rhofunc = 'tukey',
                      nmad = 3, printout = NULL,
                      width = 640, height = 500){
  mat$rotationw = 1-exp(-mat$length/weight.constant)
  mat$rhofunc = rhofunc
  mat$nmad = nmad
  evaluate.weights(mat, mystart$alphamax, mystart$F, 
                   nmad = nmad, w=mat$rotationw)
  title(paste('nmad = ', nmad, ', w = 1-exp(-segment length/', weight.constant, ')', 
              sep = ''), adj = 1, outer = T, line = -2, font.main = 1)
  mat$weight.constant = weight.constant
  print("Change weight constant until 1'000'000 length segs have almost top weight")
  print("Change nmad until Contrib to rho  gets constant at good Distance to lattice point")
  if(!is.null(printout)) dev.print(png, file = printout, width= width, height = height)
  mat
}


evaluate.weights = function(mat, alpha, F, rhofunc = 'tukey', 
                            nmad = 1, maxCN = 24, w=1){
  #Recalculate transformed lattice points from F and alpha
  x = cbind(mat$a1, mat$a2)
  Fm1 = solve(F)  
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  
  E1 = c(alpha*A1 + (1-alpha))
  E2 = c(alpha*A2 + (1-alpha))
  
  ### Positions of lattice points in observed scale and
  ### each segment's distance to it's lattice point -> weight
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  alldist = apply(x, 1, dist, y = cbind(e1, e2))
  closest = apply(alldist, 2, which.min)
  e = cbind(e1, e2)[closest,]
  
  di = dist(x, e)
  r = di*w
  rhoargs = list(k = nmad*mad(di), c = nmad*mad(di))
  eachrho = do.call(rhofunc, args = list(x = r, args = rhoargs))
  
  #Plot
  par(mfrow=c(2,2), mar = c(4, 4, 3, .5))
  plot(r, eachrho, log = 'xy', xlab = 'r = Distance * Weight', 
       ylab = 'Contribution to rho')
  plot(di, eachrho, log = 'xy', xlab = 'Distance to lattice point', 
       ylab = 'Contribution to rho')
  plot(mat$length, w, log = 'xy', xlab = 'Segment length', 
       ylab = 'Weight')
  plot(di, r, log = 'xy', ylab = 'r = Distance * Weight', 
       xlab = 'Distance to lattice point')
  title(main = 'Evaluation of rho function and weights', line = -1,
        outer = T)

}


plot.transformed = function(mat, start, w=NULL, rhofunc = 'leastsq',
                            rhoargs = NULL, main = NULL, xlim = NULL, ylim = NULL,
                            maxCN = 24, color.by = 'eachrho',
                            printout = NULL, width = 640, height = 500){
  require('RColorBrewer')
  require('graphics')
  F = start$F
  alpha = start$alpha
  #Recalculate transformed lattice points from F and alpha
  Fm1 = solve(start$F)  
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  
  E1 = c(alpha*A1 + (1-alpha))
  E2 = c(alpha*A2 + (1-alpha))
  
  ### Positions of lattice points in observed scale and
  ### each segment's distance to it's lattice point -> weight
  x = cbind(mat$a1, mat$a2)
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  alldist = apply(x, 1, dist, y = cbind(e1, e2))
  closest = apply(alldist, 2, which.min)
  e = cbind(e1, e2)[closest,]
  if (is.null(w)) w = 1
  r = dist(x, e)*w
  eachweight = do.call(rhofunc, args = list(x = r, args = rhoargs))
  
  #Colour palette
  #display.brewer.all()
  pal = brewer.pal(9, 'YlOrRd')
  #display.brewer.pal(8, 'Set1')
  cols = pal[-c(1:2)]
  cols = colorRampPalette(cols)(100)
  
  #Transformed positions of data
  a1p = x[,1]*Fm1[1,1] + x[,2] * Fm1[1,2]
  a2p = x[,1]*Fm1[2,1] + x[,2] * Fm1[2,2]
  
  par(mfrow = c(1,1))
  split.screen(rbind(c(0.2, 0.98, 0.2, 0.9),
                     c(0.1,0.2,0.2, 0.9), 
                     c(0.2, 0.98, 0.1, 0.2),
                     c(0.1, 0.2, 0.1, 0.2)))
  screen(1)
  par(mar = c(0,0,0,0))
  
  
  #Plot frame
  if (is.null(xlim)) xlim = c(0, max(c(a2p, E2), na.rm = T))
  if (is.null(ylim)) ylim = c(0, max(a1p, na.rm = T)*1.2)
  plot(a2p, a1p, pch = 16, 
       type = 'n', xlim = xlim, ylim = ylim, xaxt = 'n', yaxt = 'n')
  
  #Points
  colby = eachweight
  if (color.by[1] == 'log') colby = log(eachweight+1)
  if (is.numeric(color.by)) colby = color.by
  colbreaks = seq(min(colby), max(colby[a2p<xlim[2]]), length = length(cols))
  delta = colbreaks[2] - colbreaks[1]
  colind = (colby-colbreaks[1]) %/% delta + 1
  colind[colind>length(cols)] = length(cols)
  mycols = col.transp(cols[colind], .5)
  if (color.by == 'nocol') mycols = col.transp('grey', .3)
  points(a2p, a1p, pch = 16, col = mycols, cex = 1)
  
  #Gridlines
  EEs = unique(E1)
  for (j in 1:length(EEs)){
    lines(range(EEs), rep(EEs[j], 2), lty = 'dotted', col = 'grey20')
    lines(rep(EEs[j], 2), range(EEs), lty = 'dotted', col = 'grey20')
  }
  abline(0,1)
  
  #Colour scale
  if (color.by != 'nocol'){
    tmp = c(xlim[2]/2, xlim[2])
    tmp = seq(tmp[1], tmp[2], length=length(cols))
    points(tmp, rep(0, length(tmp)), col = cols, pch = 16)
  }
  
  #Single scale density y
  screen(2)
  par(mar = c(0,0,0,0))
  yden = density(a1p,kernel="gaussian", bw = .03)
  roughy = round(sum(diff(yden$y)**2)/diff(yden$x[1:2]), 3)
  yden$y = log(yden$y+1)
  plot(1, 1, ylim = ylim, xlim = range(yden$y), type = 'n', xaxt = 'n')
  polygon(yden$y, yden$x, col = col.transp(1, .5), border = NA)
  title(ylab = 'a1', line = -1, outer = T)
  
  #Single scale density x
  screen(3)
  par(mar = c(0,0,0,0))
  xden = density(a2p,kernel="gaussian", bw = .03)
  roughx = round(sum(diff(xden$y)**2)/diff(xden$x[1:2]), 3)
  xden$y = log(xden$y+1)
  plot(1, 1, xlim = xlim, ylim = range(xden$y), type = 'n', yaxt = 'n')
  polygon(xden$x, xden$y, col = col.transp(1, .5), border = NA)
  title(xlab = 'a2', line = -1, outer = T)
  
  screen(4)
  par(mar = c(0,0,0,0))
  plot(1,1, type = 'n', xaxt = 'n', yaxt = 'n', xlim = c(0,1), ylim = c(0,1), bty = 'n')
  text(0, 1, roughx, adj = c(0,1))
  text(1, 0, roughy, adj = c(1,0))
  
  title(main = main, outer = T, line = -1)
  close.screen(all.screens = TRUE)
  if(!is.null(printout)) dev.print(png, file = printout, width= width, height = height)
  
  c(roughx = roughx, roughy=roughy)
  
}

start.alphamax.f = function(mat, maxCN = 24, maxlines = 30, xlim = NULL, ylim = NULL,
                            colbychrom = FALSE, force.diag = TRUE, vertical.cluster.line=NULL,horizontal.cluster.line=NULL,allelic.balance.cluster=NULL,
                            dx.eq.dy = FALSE, printout = NULL, width = 640, height = 500, ...){
  ##dx.eq.dy: Scale axis equally (logic)?
  ##force.diag: Force points on the diagonal to stay on diagonal (logic)?
  #if (is.null(xlim)) xlim = c(0, quantile(mat$a2, probs = c(.9)))
  #if (is.null(ylim)) ylim = c(0, quantile(mat$a1, probs = c(.975)))
  cols = col.transp('grey', .3)
  if (colbychrom) 
  {
    cols = col.transp(as.numeric(mat$Chromosome), .3)
  }
  #plot unscaled values
  plot(mat$a2, mat$a1, pch = 16, col = cols, 
       xlim = xlim, ylim = ylim, xlab = 'Major intensity ratio a2', 
       ylab = 'Minor intensity ratio a1',
       main = 'Minor and major CN intensity')
  abline(0,1, col = 'red')
  
  posy = vertical.cluster.line
  if (is.null(vertical.cluster.line))
  {
    print('Mark two subsequent CN clusters along a "vertical" (to be) line:')
    posy = locator(2)
  }
  posl = posy
  #print(posy)
  #  print('Mark two distant points in one "vertical" line:')
  #  posl = locator(2)
  posx = horizontal.cluster.line
  if (is.null(horizontal.cluster.line))
  {
    
    print('Mark two subsequent CN clusters along a "horizontal" (to be) line:')
    posx = locator(2)
  }
  #print(posx)
  #  print('Mark two distant points in one "horizontal" line:')
  #  posk = locator(2)
  posk = posx


  posb = allelic.balance.cluster
  if (is.null(allelic.balance.cluster))
  {
    print('Mark an allelic balance cluster:')
    posb = locator(1)
  }
  #print(posb)
  dy = abs(diff(posy$y))
  lines(posx$x, posx$y)
  dx = abs(diff(posx$x))
  l = lm(x~y, data = posl)$coef[2]
  k = lm(y~x, data = posk)$coef[2]
  lines(posy$x, posy$y)
  lines(posl$x, posl$y, lty = 'dotted')
  lines(posk$x, posk$y, lty = 'dotted')
  points(posb$x, posb$y, pch = 4)
  if (dx.eq.dy){
    dx = (dx+dy)/2
    dy = dx
  }
  if (force.diag){
    k = (k+l)/2
    l = k
  }
  f11 = dy
  f22 = dx
  
  F = matrix(c(f11, l*f11, k*f22, f22), ncol = 2, nrow = 2)
  Fm1 = solve(F)
  
  
  ### Starting value positions of lattice points in E scale
  E1 = Fm1[1,1]*posy$y + Fm1[1,2]*posy$x
  deltay = max(E1) - min(E1)
  b = seq(max(E1), by = -deltay, to = 0)
  by = b[length(b)]
  ny = length(b)
  E1 = seq(by, by = deltay, length.out = maxlines)
  
  E2 = Fm1[2,1]*posx$y + Fm1[2,2]*posx$x
  deltax = max(E2) - min(E2)
  b = seq(max(E2), by = -deltax, to = 0)
  bx = b[length(b)]
  nx = length(b)
  E2 = seq(bx, by = deltax, length.out = maxlines)
  
  E1 = rep(E1, maxlines)
  E2 = rep(E2, each = maxlines)
  
  ### Starting value positions of lattice points in observed scale
  ### Before adjusted to common alpha
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  points(e2, e1)
  
  ### Estimate startvalue of maximum alpha which is <1
  balance.x = Fm1[2,1]*posb$y + Fm1[2,2]*posb$x
  balance.y = Fm1[1,1]*posb$y + Fm1[1,2]*posb$x
  dists = dist(x = c(balance.y, balance.x), y = cbind(E1, E2))
  closest = which.min(dists)
  bbx = seq(E2[closest], by = -deltax, to = 0)
  bby = seq(E1[closest], by = -deltay, to = 0)
  n = min(length(bbx), length(bby))
  nxs = which(bbx>=0)[(length(bbx)-n+1):length(bbx)] 
  nys = which(bby>=0)[(length(bby)-n+1):length(bby)] 
  
  alphax = (1)/(bx+nxs)
  alphay = 1/(by+nys)
  print(paste('Maximum % cells in main subclone (y):', alphay[1]))
  print(paste('Maximum % cells in main subclone (x):', alphax[1]))
  alphamax = mean(c(alphax[1], alphay[1]))
  
  ### Re-estimate starting value positions of Es and F from alpha
  f11 = dy/alphamax#(dx+dy)/2/alphamax#
  f22 = dx/alphamax#(dx+dy)/2/alphamax#
  F = matrix(c(f11, l*f11, k*f22, f22), ncol = 2, nrow = 2)
  Fm1 = solve(F)
  
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  
  E1 = c(alphamax*A1 + (1-alphamax))
  E2 = c(alphamax*A2 + (1-alphamax))
  
  ### Starting value positions of lattice points in observed scale
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  #  points(e2, e1, col = 'blue')
  
  ### Plot starting values on transformed scale
  par(mfrow = c(1,1))
  x = cbind(mat$a1, mat$a2)
  a1 = Fm1[1,1] * x[,1] + Fm1[1,2]*x[,2]
  a2 = Fm1[2,1] * x[,1] + Fm1[2,2]*x[,2]
  
  if(!is.null(printout)) dev.print(png, file = printout, width= width, height = height, ...)
  
  list(alphamax = alphamax, F=F)
}

leastsq = function(x, ...){
  x^2
}

rho.alphamax.f = function(pars, x, w, maxCN=24, rhofunc = 'tukey', 
                          rhoargs = NULL){
  alphamax = pars[1]
  F = matrix(pars[-1], ncol = 2, nrow = 2)
  
  #Position of latice points on observed scale
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  E1 = c(alphamax*A1 + (1-alphamax))
  E2 = c(alphamax*A2 + (1-alphamax))
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  
  #Assign observations to closest lattice point
  alldist = apply(x, 1, dist, y = cbind(e1, e2))
  closest = apply(alldist, 2, which.min)
  
  #Calculate residuals
  e = cbind(e1, e2)[closest,]
  r = dist(x, e)*w
  sum(do.call(rhofunc, args = list(x = r, args = rhoargs)))
}
plot.cn = function(x, y, main = NULL, xlim = NULL, xlab = NULL, ylab = NULL,
                   ylim = NULL, maxCN = 24, mycols = NULL, grid = T){
  E1 = E2 = 0:maxCN
  
  #Plot frame
  if (is.null(xlim)) xlim = c(0, max(c(x, E2), na.rm = T))
  if (is.null(ylim)) ylim = c(0, max(y, na.rm = T)*1.2)
  type = ifelse(is.null(mycols), 'p', 'n')
  plot(x, y, pch = 16, col = col.transp('blue', .3), type = type,
       xlim = xlim, ylim = ylim, xlab = '', ylab = '')
  if (is.vector(mycols)){
    points(x, y, pch = 16, col = mycols, cex = 1)    
  }
  if (is.list(mycols)) {
    points(x, y, pch = 16, col = mycols$mycols, cex = 1)
    
    #Colour scale
    tmp = c(xlim[2]/2, xlim[2])
    tmp = seq(tmp[1], tmp[2], length=length(mycols$cols))
    points(tmp, rep(ylim[1], length(tmp)), col = mycols$cols, pch = 16)
    
  }
  #Gridlines
  if (grid){
    EEs = unique(E1)
    for (j in 1:length(EEs)){
      lines(range(EEs), rep(EEs[j], 2), lty = 'dotted', col = 'grey20')
      lines(rep(EEs[j], 2), range(EEs), lty = 'dotted', col = 'grey20')
    }
  }
  abline(0,1)
  title(main = main, outer = T, line = -1)
  title(xlab = xlab, ylab = ylab, line = 2)
}





