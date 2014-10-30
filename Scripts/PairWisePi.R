#source('~/Documents/Academics/CoopLab/Projects/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
#source ( "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Scripts/run.ms.functions.R")
source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
source ( "~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R")

load ( "11000freq.trajectories.Rdata" )


N = 10000
s = 0.01
real.fs <- c ( 0.001 , 0.005 , 0.01 , seq ( 0.02 , 0.2 , length = 10 ) )
Ts <- list ()
sweep.starts <- list ()

for ( i in 1 : 13 ) {
	this.f <- my.runs [[ i ]] [[ 1 ]]
	Ts [[ i ]] <- numeric ( nrow ( this.f ) )
	sweep.starts [[ i ]] <- apply ( this.f , 1 , function ( x ) which ( x == real.fs [ i ] )  )
	
	
	for ( j in 1 : nrow ( this.f ) ) {
		Ts [[ i ]] [ j ] <- sum ( 1 - this.f [ j , 1 : sweep.starts [[ i ]] [ j ] ] )
	}
	
}

lapply ( Ts , mean )

T_f <- numeric ( length ( real.fs ) )
for ( i in 1 : length ( real.fs ) ) {
	
	T_f [ i ] <- log ( (2*N -1 ) * ( 1 - real.fs [ i ] ) / real.fs [ i ] ) / s
	
	
}

unlist ( lapply ( sweep.starts , mean ) ) - T_f



my.logistic <- function ( x ) 1 / (2 * N  ) * exp(s * x ) / ( 1 + 1 / (2 * N  ) * ( exp(s * x )  - 1 )  )


Tsf <- numeric ( length ( real.fs ) )
for ( i in 1 : length ( real.fs ) ) {
	
	Tsf [ i ]  <- integrate ( my.logistic , 0 , T_f [ i ] )$value
		
}

true.Tsf <- unlist ( lapply ( Ts , mean ) )
Tsf - unlist ( lapply ( Ts , mean ) )



plot ( c ( 0 , 200 ) , c ( 0 , 1 ) , type = "n" , xlab = "4NR" , ylab = expression ( pi [ R ] / pi [ 0 ] ) , cex.lab = 1.5 )
pi.over.f<-list()



setwd("~/Documents/Academics/StandingSweeps/Sims")

s = 0.05
N = 10000
pi.over.f <- list()
for ( i in 1 : length ( real.fs ) ) {
	
		f.lab <- strsplit ( as.character ( real.fs [ i ] ) , "\\." ) [[ 1 ]] [ 2 ]
		s.lab <- strsplit ( as.character ( s ) , "\\." ) [[ 1 ]][ 2 ]
		
		mut.density <- get.mut.density ( file = paste ( "mssel_f2" , f.lab , s.lab , N , ".out" , sep = "" ) )
		pi.over.f [[ i ]] <- mut.density
		lines ( pi.over.f [[ i ]]$x  , pi.over.f [[ i ]]$y/(1000*20) , col = i )
}






plot ( c ( 0 , 200 ) , c ( 0 , 1 ) , type = "n" , xlab = "4NR" , ylab = expression ( pi [ R ] / pi [ 0 ] ) , cex.lab = 1.5 )
for ( i in 1 : length ( real.fs ) ) {
		#run.ms.f(f.index)
		# mut.density <- get.mut.density ( file = paste ( "mssel_f2" , f.index , ".out" , sep = "" ) )
		# pi.over.f [[ f.index ]] <- mut.density
		lines ( pi.over.f [[ i ]]$x  , pi.over.f [[ i ]]$y/(1000*20) , col = i )
}


R <- 1:200
r <- R / ( 4*N )

for ( i in 1 : length ( real.fs ) ) {
		
	lines ( R , 1 - exp ( - 2 * r * Tsf [ i ]  ) / ( 1 + R * real.fs [ i ] * ( 1 - real.fs [ i ] ) ) , lty = 2 , lwd = 2 , col = i )

}





#source('~/Documents/Academics/CoopLab/Projects/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
#source ( "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Scripts/run.ms.functions.R")
source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
source ( "~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R")


real.fs <- c ( 0.001 , 0.01 , 0.05 , 0.1 )
my.runs <- lapply ( real.fs , function ( x ) SweepFromStandingSim ( N = 10000 , s = 0.05 , f = x , reps = 1000 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T ) )
many.sim.freqs <- list ( my.runs , real.fs )
save ( many.sim.freqs  ,  file = "~/Documents/Academics/StandingSweeps/Sims/freqs.traj.s05.Rdata"  )
for ( i in 1 : length ( real.fs ) ) {
	run.ms.f ( runs = my.runs [[ i ]] [[ 1 ]] , n.sam = 12 , f = real.fs [ i ] , s = 0.05 , N = 10000 , path = "~/Documents/Academics/StandingSweeps/" )
}






