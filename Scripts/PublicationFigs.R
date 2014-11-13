#### plot of pi and segregating sites
source('~/Documents/Academics/CoopLab/Projects/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
source ( "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Scripts/run.ms.functions.R")

real.fs <- c ( 1 / 20000 , 0.001 , 0.01 , 0.05 , 0.1 )

T_f <- numeric ( length ( real.fs ) )
for ( i in 1 : length ( real.fs ) ) {
	
	T_f [ i ] <- log ( (2*N -1 ) * ( 1 - real.fs [ i ] ) / real.fs [ i ] ) / s
	
}