chptdat = read.table("http://www.stat.psu.edu/~mharan/MCMCtut/COUP551_rates.dat",skip=1) 
Note: This data set is just a convenient subset of the actual data set (see reference below.) 

We can begin with a simple time series plot as exploratory analysis. 
   Y=chptdat[,2] # store data in Y
   ts.plot(Y,main="Time series plot of change point data")
The plot suggests that the change point may be around 10. 