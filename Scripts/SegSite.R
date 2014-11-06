source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
source ( "~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R")



real.fs <- c ( 0.001 , 0.01 , 0.05 , 0.1 )
my.runs <- lapply ( real.fs , function ( x ) SweepFromStandingSim ( N = 10000 , s = 0.05 , f = x , reps = 1000 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T ) )
many.sim.freqs <- list ( my.runs , real.fs )
save ( many.sim.freqs  ,  file = "~/Documents/Academics/StandingSweeps/Sims/freqs.traj.s05.Rdata"  )
for ( i in 1 : length ( real.fs ) ) {
	run.ms.f ( runs = my.runs [[ i ]] [[ 1 ]] , n.sam = 12 , f = real.fs [ i ] , s = 0.05 , N = 10000 , path = "~/Documents/Academics/StandingSweeps/" )
}



plot ( c ( 0 , 200 ) , c ( 0 , 1 ) , type = "n" , xlab = "4NR" , ylab = expression ( S [ R ] / S [ 0 ] ) , cex.lab = 1.5 )

setwd("~/Documents/Academics/StandingSweeps/Sims")

s = 0.05
N = 10000
S.over.f <- list()
for ( i in 1 : length ( real.fs ) ) {
	
		f.lab <- strsplit ( as.character ( real.fs [ i ] ) , "\\." ) [[ 1 ]] [ 2 ]
		s.lab <- strsplit ( as.character ( s ) , "\\." ) [[ 1 ]][ 2 ]
		
		mut.density <- get.mut.density ( file = paste ( "mssel_f2" , f.lab , s.lab , N , ".out" , sep = "" ) )
		S.over.f [[ i ]] <- mut.density
		lines ( S.over.f [[ i ]]$x  , S.over.f [[ i ]]$y/(1000*20) , col = i )
}



my.logistic <- function ( x ) 1 / (2 * N  ) * exp(s * x ) / ( 1 + 1 / (2 * N  ) * ( exp(s * x )  - 1 )  )
Tsf <- numeric ( length ( real.fs ) )
T_f <- numeric ( length ( real.fs ) )



R <- 1:200
r <- R / ( 4*N )
n = 12
my.time <- list ()

for ( i in 1 : length ( real.fs ) ) {
	my.time [[ i ]] <- numeric ( length ( r ) )
	T_f [ i ] <- log ( (2*N -1 ) * ( 1 - real.fs [ i ] ) / real.fs [ i ] ) / s	
	Tsf [ i ]  <- integrate ( my.logistic , 0 , T_f [ i ] )$value
	P_NR <- exp ( - r * Tsf [ i ] )
	times <- 1 / ( 1 : ( n - 1 ) )
	tot.times <- c ( 0 , cumsum ( times ) )
	tot.times.mat <- matrix ( 0 , nrow = length ( tot.times ) , ncol = length ( tot.times ) )
	idx <- -1*(row ( tot.times.mat )-col ( tot.times.mat ))+n
	idx [ idx > n ] <- 1
	for ( k in 1 : ncol ( tot.times.mat ) ) {	
		tot.times.mat [ , k ] <- tot.times [ idx [ , k ]  ]
	}
	for ( j in 1 : length ( r ) ) {
		ewens <- EwensDist ( n = n , N = 10000 , r = 10^-8 , distance = r[ j ]/10^-8 , f = real.fs [ i ] )
		binom <- dbinom ( 1 : n , n , P_NR [ j ] )
		my.time [[ i ]] [ j ] <- ( sum ( binom * ewens * tot.times.mat ) + dbinom ( 0 , n , P_NR [ j ] ) * tail ( tot.times , 1 ) ) / tail ( tot.times , 1 )
		
	}
	
}




for ( i in 1 : length ( real.fs ) ) {
		
	lines ( R , my.time [[ i ]] , lty = 2 , lwd = 2 , col = i )

}





