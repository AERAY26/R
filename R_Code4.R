# Spatial statistics
# based on Modern Statistical Methods for Astronomy with R Applications
# Penn State Summer School   June 2019 


# I Setup R infrastructure

setwd('/Users/ericfeigelson/Desktop/Rdir')
install.packages('scatterplot3d')  ;  library(scatterplot3d)
install.packages('spatstat')  ;  library(spatstat)
install.packages('alphahull')  ;  library(alphahull)
install.packages('dixon')  ;  library(dixon)


# II Input EFIGI galaxy morphology catalog

efigi <- read.csv('EFIGI.csv', head=TRUE)
summary(efigi)
efigi[,5] <- as.numeric(as.character(efigi[,5]))
efigi[,9] <- as.numeric(as.character(efigi[,9]))
summary(efigi)
attach(efigi)


# III Examine spatial distribution of full and trimmed sample
# Divide T morphological measure into classical galaxy classes

par(mgp=c(2,0.5,0))
plot(RA, Dec, pch=20, cex=0.2, asp=1, xlim=c(360,0), ylim=c((-20),80),
     xlab='R.A.', ylab='Dec')
trim <- which(RA > 150 & RA < 250 & Dec > 17)  
RA.trim <- RA[trim]  ;  Dec.trim <- Dec[trim]  
T.trim <- T[trim]  ;  DCom.trim <- DCom[trim]
length(RA.trim)
plot(RA.trim, Dec.trim, pch=20, cex=0.3, asp=1, xlim=c(250,150), 
     ylim=c(15,75), xlab='R.A.', ylab='Dec')   # not shown

scatterplot3d(RA.trim, DCom[trim], Dec.trim, pch=20, cex.symbols=0.3, 
     xlab='RA', ylab='Dist (Mpc)', zlab='Dec')

galclass.trim <- cut(T.trim, breaks=c( (-7), (-4), (-1), 9, 10, 11), 
     labels=c('Ell', 'Len', 'Sp', 'Irr', 'Dw'))
summary(galclass.trim)


# IV Construct alpha-shape (non-convex hull) polygon around 2D point distribution 
# Sequence vertices in approximate clockwise order; correct mistake in polygon

ash <-  ashape(cbind(RA.trim, Dec.trim), alpha=5.)	
RA.extreme <- RA.trim[ash$alpha.extremes] 
Dec.extreme <- Dec.trim[ash$alpha.extremes] 
    angle <- function(ra,dec) { 			
	   tmp <- ra + 1i * dec
	   outang <- 90 - Arg(tmp) / pi * 180
	   outang %% 360
	 }
poly.ind <- rank(angle(RA.extreme - mean(RA.extreme), 
     Dec.extreme - mean(Dec.extreme)))
plot(RA.trim, Dec.trim, pch=20, cex=0.3)		# not shown
lines(RA.extreme[order(poly.ind)], Dec.extreme[order(poly.ind)])

RA.extreme <- RA.extreme[-35]  ;  Dec.extreme <- Dec.extreme[-35]
poly.ind <- rank(angle(RA.extreme - mean(RA.extreme), 
     Dec.extreme - mean(Dec.extreme)))
plot(RA.trim, Dec.trim, pch=20, cex=0.3)		# not shown
lines(RA.extreme[order(poly.ind)], Dec.extreme[order(poly.ind)])


# V Setup up spatstat infrastructure, construct convenient subsets, 
# and show smoothed spatial maps

trim.window <- owin(poly=list( x=RA.extreme[rev(order(poly.ind))], 
     y=Dec.extreme[rev(order(poly.ind))] ))
centroid.owin(trim.window)  ;  area.owin(trim.window)
gal.ppp <- ppp(RA.trim, Dec.trim, marks=galclass.trim, window=trim.window)
unitname(gal.ppp) <- 'deg'
str(gal.ppp)
plot(gal.ppp, cols=1, cex=0.5, main='')

Ell.gals <- split(gal.ppp)$Ell  ;  Len.gals <- split(gal.ppp)$Len
Sp.gals <- split(gal.ppp)$Sp  ;  Irr.gals <- split(gal.ppp)$Irr 
Dw.gals <- split(gal.ppp)$Dw
Ell.Sp.gals <- gal.ppp[gal.ppp$marks == 'Sp' | gal.ppp$marks == 'Ell',]
Irr.Sp.gals <- gal.ppp[gal.ppp$marks == 'Sp' | gal.ppp$marks == 'Irr',]

par(mfrow=c(3,2), mgp=c(2,0.5,0), mar=c(0,0,2,1)+0.1)
plot(density(Ell.gals, bw.ppl(Ell.gals)), col=gray(15:50/50), 
     main='Elliptical', ribsep=0.05)
plot(density(Len.gals, bw.ppl(Len.gals)), col=gray(15:50/50), 
     main='Lenticular', ribsep=0.05)
plot(density(Sp.gals, bw.ppl(Sp.gals)), col=gray(15:50/50), 
     main='Spiral', ribsep=0.05)
plot(density(Irr.gals, bw.ppl(Irr.gals)), col=gray(15:50/50), 
     main='Irregular', ribsep=0.05)
plot(density(Dw.gals, bw.ppl(Dw.gals)), col=gray(15:50/50), 
     main='Dwarf', ribsep=0.05)
plot(density(Irr.gals, bw.ppl(Irr.gals), se=TRUE)$SE, col=gray(15:50/50), 
     main='Irregular standard error', ribsep=0.05)


# VI Spatial visualizations: Smoothed, dominant class, and relative risk maps

bandw <- bw.ppl(gal.ppp)
par(mfrow=c(1,3), mar=c(1,0,0,1), mgp=c(0,0.5,0))
Irr.gals.smooth <- density(Irr.gals, sigma=bandw*3)
plot(Irr.gals.smooth, col=gray(15:50/50), ribsep=0.05,main='')
points(Irr.gals, pch=20, cex=0.5)
contour(Irr.gals.smooth, nlevels=5, xlab='RA', ylab='Dec', 
     zlab='Irregular galaxies', main='')
points(Irr.gals, pch=20, cex=0.5)
Irr.gals.adap <- adaptive.density(Irr.gals, f=0.2, nrep=100)
plot(Irr.gals.adap, col=gray(15:50/50), ribsep=0.05,main='')
points(Irr.gals, pch=20, cex=0.5)

class.prob <- relrisk(gal.ppp, sigma=0.5)
class.dominant <- im.apply(class.prob, which.max)
classes <- levels(marks(gal.ppp))
class.dominant <- eval.im(factor(class.dominant, levels=1:5, 
     labels=classes))
plot(class.dominant, main='', ribsep=0.05, col=gray(0:4/4))

plot(relrisk(Ell.Sp.gals, relative=TRUE, control='Sp'), main='') 
plot(relrisk(Ell.Sp.gals, relative=TRUE, control='Sp', se=TRUE)$SE, 
     main='')
plot(relrisk(Irr.Sp.gals, relative=TRUE, control='Sp'), main='') 
plot(relrisk(Irr.Sp.gals, relative=TRUE, control='Sp', se=TRUE)$SE, 
     main='')


# VII Nonparametric measures of spatial differences and correlations by class

nncorr(gal.ppp) 		# 1.00
nncorr(Ell.Sp.gals) 		# 1.01
nncorr(Irr.Sp.gals) 		# 1.01

segregation.test(gal.ppp, nsim=50)	# P<0.02
segregation.test(Ell.Sp.gals, nsim=50)	# P>0.05
segregation.test(Irr.Sp.gals, nsim=50)	# P<0.02

dixon(as.data.frame(gal.ppp))$tablaZ

plot(alltypes(gal.ppp, Kcross))		
plot(alltypes(gal.ppp, Lcross))		

plot(alltypes(gal.ppp, markconnect))


# VIII Hypothesis tests for Complete Spatial Randomness

summary(Irr.gals)
clarkevans.test(Irr.gals, alternative='clustered') 				# Clark-Evans test
hopskel.test(Irr.gals, alternative='clustered')				# Hopkins-Skellam test
mad.test(Irr.gals, interpolate=TRUE)  					# median absolute deviation test
dclf.test(Irr.gals, interpolate=TRUE) 						# Diggle-Cressie-Loosmore-Ford test
dg.test(Irr.gals, nsim=20, use.theory=TRUE, interp=TRUE) # Dao-Genton test
plot(dg.envelope(Irr.gals, Lest, nsim=50), main='')


# IX Clustering analysis for lenticular galaxies

par(mfrow=c(3,2), mgp=c(2,0.5,0), mar=c(3,4,2,1)+0.1)
plot(envelope(Len.gals, fun=Gest), legend=FALSE, main='G function')
plot(envelope(Len.gals, fun=Fest), legend=FALSE, main='F function')
plot(envelope(Len.gals, fun=Kest), legend=FALSE, main='K function')
plot(envelope(Len.gals, fun=Jest), legend=FALSE, main='J function')
plot(envelope(Len.gals, fun=pcf), legend=FALSE, 
     main='Pair correlation function')

pcf.CSRsim <- envelope(Len.gals, fun=pcf, log='xy')
str(pcf.CSRsim)
plot(pcf.CSRsim$r, pcf.CSRsim$obs, log='xy', type='l', lwd=2, 
     main='Log pair correlation function', xlab='log(r) deg', ylab='log(pcf)')
lines(pcf.CSRsim$r, pcf.CSRsim$theo, type='l', lty=2)
xx <- c(0, pcf.CSRsim$r, rev(pcf.CSRsim$r), 0)
yy <- c(c(0, pcf.CSRsim$lo), rev(c(0, pcf.CSRsim$hi)))
polygon(xx, yy, col='#11111120')
pcf.boot <- lohboot(Len.gals, fun='pcf', delta=0.25)
xx <- c(0, pcf.boot$r, rev(pcf.boot$r), 0)
yy <- c(c(0, pcf.boot$lo), rev(c(0, pcf.boot$hi)))
polygon(xx, yy, col='#11111150')
lines(c(0.5, 5), c(40, 6.34), lty=3, lwd=2)

plot(envelope(Len.gals, fun=Ginhom), legend=FALSE, main='G function')
plot(envelope(Len.gals, fun=Finhom), legend=FALSE, main='F function')
plot(envelope(Len.gals, fun=Kinhom), legend=FALSE, main='K function')
plot(envelope(Len.gals, fun=Jinhom), legend=FALSE, main='J function')
plot(envelope(Len.gals, fun=pcfinhom), legend=FALSE, 
     main='Pair correlation function')


# X Dirichlet-Voronoi tessellation with close galaxy pairs/groups

Vor.gal <- dirichlet(gal.ppp) 
plot(Vor.gal, main='')  
gal.area <- lapply(tiles(Vor.gal), area.owin) 
points(gal.ppp, pch=20, cex=0.2 )

dev.new() ; hist(unlist(gal.area), breaks=100, main='')  # not shown
gal.group <- cut(as.numeric(gal.area), breaks=c(0,1,100))
dev.set(dev.prev())
points(gal.ppp[gal.group == '(0,1]'], pch=20, cex=0.7)


# XI Spatially averaged maps of morphological class

gal2.ppp <- ppp(RA.trim, Dec.trim, marks=T.trim, window=trim.window)
bw.smoothppp(gal2.ppp)
Tmorph.sm <- Smooth(gal2.ppp, bw.smoothppp, diggle=TRUE)
Tmorph.sd <- sqrt(markvar(gal2.ppp, bw.smoothppp))
Tmorph.fit <- Smooth(gal2.ppp, bw.smoothppp, at='points')
Tmorph.resid <- abs(marks(gal2.ppp) - Tmorph.fit)
Tmorph.nn <- nnmark(gal2.ppp)

par(mfrow=c(2,2), mgp=c(2,0.5,0), mar=c(0,0,1,1)+0.1)
plot(density(gal2.ppp, sigma=1), col=gray(15:50/50), main='', ribsep=0.05)
plot(Tmorph.sm, col=gray(15:50/50), main='', ribsep=0.05)
plot(Tmorph.sd, col=gray(15:50/50), main='', ribsep=0.05)
points(gal2.ppp$x[which(Tmorph.resid > 8)], 
     gal2.ppp$y[which(Tmorph.resid > 8)], pch=4)
plot(Tmorph.nn, col=gray(15:50/50), main='', ribsep=0.05)


# XII Diagnostics for dependence of location and mark variables

gal2_pos.ppp <- gal2.ppp
gal2_pos.ppp$marks <- gal2.ppp$marks + 7
Tmorph.emark <- Emark(gal2_pos.ppp)
Tmorph.emark$iso <- Tmorph.emark$iso - 7
Tmorph.emark$theo <- Tmorph.emark$theo - 7
Tmorph.emark$trans <- Tmorph.emark$trans - 7
Tmorph.vmark <- Vmark(gal2_pos.ppp)
par(mfrow=c(1,2))
plot(Tmorph.emark, main='')
plot(Tmorph.vmark, main='')

nnmean(gal2.ppp)		# 0.95
nnvario(gal2.ppp)		# 0.86
nncorr(gal2.ppp)		# 0.145
corr.temp <- nncorr(gal2.ppp)[['correlation']]
sqrt(2 * (1 - corr.temp^2) / (npoints(gal2.ppp) - 4))	# approx 0.035