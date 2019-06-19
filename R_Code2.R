R_Code2

dev.copy2pdf(file='name.pdf')


# Density estimation and local regression 
# Eric Feigelson
# Penn State Summer School  June 2019


setwd('/Users/ericfeigelson/Desktop/Rdir')

# Regenerate bivariate dataset from first tutorial

set.seed(1)
x <- sample(seq(0.01, 3, length.out=500))
y <- 0.5*x + 0.3^(x^2) + rnorm(500, mean=0, sd=(0.05*(1+x^2)))
xy <- cbind(x, y)
plot(xy, pch=20)


# I. Kernel density estimator with two visualizations

par(mfrow=c(1,2))
library(KernSmooth)
dpik(x) ; dpik(y) 
smxy <- bkde2D(xy, bandwidth=c(dpik(x),dpik(y)))
image(smxy$x1, smxy$x2, smxy$fhat, col=topo.colors(30), xlab='Xvar', ylab='Yvar', cex.lab=1.2)
contour(smxy$x1, smxy$x2, smxy$fhat, add=T)
persp(smxy$x1,smxy$x2, smxy$fhat, theta=100, phi=40, shade=0.6, col='green1', xlab='Xvar', ylab='Yvar', zlab='Density')

# II. Two spline fits

# A cubic smoothing spline fit
# This function is based on code in the GAMFIT Fortran program by T. Hastie and R. Tibshirani (http://lib.stat.cmu.edu/general/), which makes use of spline code by Finbarr O'Sullivan.  [Chair of Statistics, Univ College Cork IE]

cubsplxy <- smooth.spline(log10(xy))			
plot(log10(xy), pch=20, cex=0.5, xlab='Xvar', ylab='Yvar')		# Plot points
plot(log10(xy), pch=20, cex=0.5, ylim=c(-0.7, 0.4), xlab='log(Xvar)', ylab='log(Yvar)', main='Two spline fits')		
lines(cubsplxy, lwd=2, col='darkgreen')		# Plot the spline fit
knotx <- cubsplxy$fit$knot*cubsplxy$fit$range + cubsplxy$fit$min 	# Find and plot the spline knots
rug(knotx,col='darkgreen')

# COnstrained B-Splines Nonparametric Regression Quantiles
# Bartels, R. and Conn A. (1980) Linearly Constrained Discrete L_1 Problems, ACM Transaction on Mathematical Software 6, 594Ã¢â‚¬â€œ608.
# Ng, P. (1996) An Algorithm for Quantile Smoothing Splines, Computational Statistics & Data Analysis 22, 99Ã¢â‚¬â€œ118.
# He, X. and Ng, P. (1999) COBS: Qualitatively Constrained Smoothing via Linear Programming; Computational Statistics 14, 315Ã¢â‚¬â€œ337.
# Ng, P. and Maechler, M. (2007) A Fast and Efficient Implementation of Qualitatively Constrained Quantile Smoothing Splines, Statistical Modelling 7(4), 315-328.


install.packages('cobs') ; library(cobs)  	
help(cobs)
cobsxy50 <- cobs(log10(x), log10(y), ic='BIC', tau=0.5)	#  Median regression fit
lines(sort(cobsxy50$x),cobsxy50$fitted[order(cobsxy50$x)], lwd=2, col=2)
cobsxy25 <- cobs(log10(x), log10(y), ic='BIC', tau=0.25)
lines(sort(cobsxy25$x),cobsxy25$fitted[order(cobsxy25$x)], lwd=1, col=2)
cobsxy75 <- cobs(log10(x), log10(y), ic='BIC', tau=0.75)
lines(sort(cobsxy75$x),cobsxy75$fitted[order(cobsxy75$x)], lwd=1, col=2)
rug(cobsxy50$knots, lwd=2, col=2)


# III. Five well-established bivariate semi-parametric local regression estimators

# LOESS,  W. S. Cleveland, `Visualizing Data', Hobart Press 1993

par(mfrow=c(2,2))
sortx <- x[order(x)] ; sorty <- y[order(x)]
local_fit <- loess(sorty ~ sortx, span=0.5, data.frame(x=x,y=y))	
summary(local_fit)
plot(x,y,pch=20, cex=0.5)
lines(sortx, predict(local_fit), lwd=2, col=2)
# Save evenly-spaced LOESS fit to a file 
x_seq <- seq(0.0, 3.0, by=0.03) 
loc_dat <- predict(local_fit, newdata=x_seq)
write(rbind(x_seq,loc_dat), sep=' ', ncol=2, file='localfit.txt')

# Nonlinear quantile regression  
# R. Koenker`Quantile Regression', Cambridge Univ Press 2005

install.packages('quantreg') ; library(quantreg)		
install.packages('MatrixModels') ; library(MatrixModels)
plot(sortx, sorty, pch=20, cex=0.5, xlab='Xvar', ylab='Yvar')
fit_rqss.25 <- rqss(sorty ~ qss(sortx), tau=0.25)
lines(sortx[-1],fit_rqss.25$coef[1]+fit_rqss.25$coef[-1], col='red')
fit_rqss.50 <- rqss(sorty ~ qss(sortx), tau=0.50)
lines(sortx[-1],fit_rqss.50$coef[1]+fit_rqss.50$coef[-1], lwd=2)
fit_rqss.75 <- rqss(sorty ~ qss(sortx), tau=0.75)
lines(sortx[-1],fit_rqss.75$coef[1]+fit_rqss.75$coef[-1], col='red')

# Nonparametric regression with bootstrap errors
# Hayfield, T. & Racine, J. S. Nonparametric Econometrics: The np package, 
# J. Statist. Software, 27(5), 2008   http://www.jstatsoft.org/v27/i05/

install.packages("np") ; library(np)	
help(npplot)
bw.NW <- npregbw(x, y, regtype='lc', bwtype='fixed') 
npplot(bws=bw.NW, ylim=c(0.0,2.5), plot.errors.method="bootstrap", 
    plot.errors.bar='I', plot.errors.type='quantiles') 
points(x, y, pch=20, cex=0.5)

# Locfit: local kernel regression methods including heteroscedastic weighting
# (unequal error bars), censoring (upper limits), and outliers.
# Loader, C. (1999). Local Regression and Likelihood. Springer, New York.
# >200 downloads/day

install.packages('locfit')  ;  library(locfit)
locfit_model <- locfit(y~lp(x, nn=0.7))
plot(locfit_model, band='local', ylim=c(0,2.5), col=2)  ;  points(xy, pch=20, cex=0.5)
locfit_model <- locfit(y~lp(x, nn=0.3))
plot(locfit_model, band='local', ylim=c(0,2.5), col=3)  ;  points(xy, pch=20, cex=0.5)


# Gaussian process regression (more commonly known as `kriging')
# Response variable errors and independent variable covariance assumed to be normal
# Maximum likelihood & Bayesian estimation
# Rasmussen & Williams, Gaussian Processes for Machine Learning, 2006

install.packages('kernlab') ; library(kernlab)
gpreg <- gausspr(x, y, variance.model=T, cross=10, kerne='polydot', kpar=list(5))
gpreg
xtest <- seq(from=min(x), to=max(x), length.out=200)
plot(x, y, pch=20, cex=0.5)
lines(xtest, predict(gpreg, xtest), col='red3', lwd=3)

# Sophisticated Gaussian processes regression codes are given in CRAN packages
# mlegp (Maximum Likelihood Estimates of Gaussian Processes) and gptk
# (Gaussian Processes tool-kit). See also the tutorial at
# http://www.r-bloggers.com/gaussian-process-regression-with-r/

###### Local regression comparison plot

npplot(bws=bw.NW, ylim=c(0.5,1.5), plot.errors.method="bootstrap",
plot.errors.bar='I', plot.errors.type='quantiles') 	# Nonparametric regression w/ bootstrap errors
points(x, y, pch=20, cex=0.5)
lines(xtest, predict(gpreg, xtest), col='red3', lwd=3)  #  Gaussian Processes regression
lines(sortx[-1],fit_rqss.25$coef[1]+fit_rqss.25$coef[-1], col='blue3')  #  Linear quantile regression
lines(sortx[-1],fit_rqss.50$coef[1]+fit_rqss.50$coef[-1], lwd=2, col='blue3')
lines(sortx[-1],fit_rqss.75$coef[1]+fit_rqss.75$coef[-1], col='blue3')
lines(sortx, predict(local_fit), lwd=2, col='darkgreen') # LOESS
locfit_values <- predict(locfit_model, seq(0,3,length.out=100))
lines(seq(0,3,length.out=100), locfit_values, lwd=2, col="chocolate")  # locfit
legend('topleft', lty=1, lwd=2, c("NP reg w/ bootstrap",'Gauss Proc reg', 'Quantile reg', 'LOESS', 'locfit'), col=c('black', 'red3', 'blue3', 'darkgreen', 'chocolate'))


