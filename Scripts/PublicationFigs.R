#### plot of pi and segregating sites
setwd ( "~/Documents/Academics/StandingSweeps/Sims" )
source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R')
source ( "~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R")
options ( scipen = 400 )



#################
## extra functions ##
#################


# MyLogistic <- function ( x , N = 10000 , s  ) 1 / ( 2 * N  ) * exp(s * x ) / ( 1 + 1 / (2 * N  ) * ( exp(s * x )  - 1 )  )

MyLogistic <- function ( x , N = 10000 , s  ) 1 / ( 2*  N ) * exp(s * x ) / ( 1 + 1 / ( 2 * N  ) * ( exp ( s * x )  - 1 )  )


#################

#################



##### trajectories vs approx
f <- 0.01
s <- 0.01
# my.freq.trajs <- SweepFromStandingSim ( N = 10000 , s = s , f = f , reps = 10 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T )
# save ( my.freq.trajs , file = "Sims/10trajectoriesForTrajFigure.Robj")
load ( file = "Sims/10trajectoriesForTrajFigure.Robj")
freqs <- t ( my.freq.trajs[[1]] )
freqs <- apply ( freqs , 2 , rev )
#freqs <- freqs [ seq ( nrow ( freqs ) - 2000 , nrow ( freqs ) ) , ]
mid.points <- apply ( freqs , 2 , function ( x ) which.min ( abs ( x - 0.5 ) ) )
diff <- max ( mid.points ) - mid.points 

new.freqs <- list()
for ( i in 1 : ncol ( freqs ) ) {	
	
	new.freqs [[ i ]] <- c ( rep ( 1 , diff [ i ] ) , freqs [ , i ]  )
	
}

my.lengths <- sapply ( new.freqs , length )
new.diff <- max ( my.lengths ) - my.lengths 
for ( i in 1 : length ( new.freqs ) ) {
	
	new.freqs [[ i ]] <- c ( new.freqs [[ i ]] , rep ( 0 , new.diff [ i ] ) )
	
}
my.freqs <- do.call ( cbind , lapply ( new.freqs , rev ) )

i <- 1
det.freqs <- 1 / ( 20000 )
while ( det.freqs [ i ] < ( 1 - f ) ) {	
	det.freqs [ i + 1 ] <- MyLogistic ( i , N = 10000 , s = s)
	i <- i + 1 
}
det.freqs <- c ( rep ( f , nrow ( my.freqs) - length ( det.freqs ) ) , rev ( 1 - det.freqs ) )
A <- length ( det.freqs ) - which.min ( abs ( det.freqs - 0.5 ) )
B <- length ( det.freqs ) - which.min ( abs ( my.freqs[,1] - 0.5 ) )
det.freqs <- c ( rep ( f , (A-B) ) , det.freqs [ - seq ( length ( det.freqs ) - (A-B) , length ( det.freqs ) ) ] )




pdf ( "Paper/Paper_Figures/TrajectoryFigure.pdf" , width = 7 , height = 4 )
matplot ( my.freqs [ 1000:3617 ,  ] , type = "l" , lwd = 0.7 , col = "grey" , lty = 1 , bty = "n" , ylab = "Frequency" , xlab = "Generations" , xaxt = "n" )
lines ( det.freqs [ 1000:3617 ] , lwd = 1.7 )
my.at <- seq ( 2617 , 0 , -500)
axis ( 1 , at = my.at , labels = seq ( 0 , 2617 , 500) )
dev.off()

new.freqs <- matrix ( 0 , nrow = max ( my.freq.trajs [[ 2 ]] ) - my.freq.trajs [[ 2 ]] [ which ( my.freq.trajs [[ 1 ]] [ , ncol ( my.freq.trajs [[ 1 ]] ) ] != 0 ) ] + ncol ( my.freq.trajs [[ 1 ]] ) , ncol = nrow ( my.freq.trajs [[ 1 ]] ) )

for ( i in 1 : nrow ( my.freq.trajs [[ 1 ]] ) ) {
	temp <- c ( rep ( 1 , max ( my.freq.trajs [[ 2 ]] ) -  my.freq.trajs [[ 2 ]] [ i ] ) , my.freq.trajs [[ 1 ]] [ i , ]  )
	if ( length ( temp ) > nrow ( new.freqs ) ) {
		new.freqs [ , i ] <- rev ( temp [ 1 : nrow ( new.freqs ) ] )
	 } else { 
	 	new.freqs [ , i ] <- rev ( c ( temp , rep ( 0 , nrow ( new.freqs ) - length ( temp ) ) ) )
	 }
}

matplot ( new.freqs , type = "l" , lwd = 0.7 , col = "grey" , lty = 1 , bty = "n" , ylab = "Frequency" , xlab = "Generations" , xaxt = "n" )
lines ( det.freqs , lwd = 2 )








##################
#### Pairwise pi ####
##################

N = 10000
s = 0.05
real.fs <- c ( 1 / 20000 , 0.001 , 0.02 , 0.05 , 0.1 )

T_f <- numeric ( length ( real.fs ) )

for ( i in 1 : length ( real.fs ) ) {
	T_f [ i ] <- log ( (2*N -1 ) * ( 1 - real.fs [ i ] ) / real.fs [ i ] ) / s
}

my.logistic <- function ( x ) 1 / (2 * N  ) * exp(s * x ) / ( 1 + 1 / (2 * N  ) * ( exp(s * x )  - 1 )  )

Tsf <- numeric ( length ( real.fs ) )
for ( i in 1 : length ( real.fs ) ) {
	Tsf [ i ]  <- integrate ( my.logistic , 0 , T_f [ i ] )$value
}


pi.over.f <- list()
for ( i in 1 : length ( real.fs ) ) {
		f.lab <- strsplit ( as.character ( real.fs [ i ] ) , "\\." ) [[ 1 ]] [ 2 ]
		s.lab <- strsplit ( as.character ( s ) , "\\." ) [[ 1 ]][ 2 ]
		mut.density <- get.mut.density ( file = paste ( "mssel_f2" , f.lab , s.lab , N , ".out" , sep = "" ) )
		pi.over.f [[ i ]] <- mut.density
}

S.over.f <- list()
for ( i in 1 : length ( real.fs ) ) {
		f.lab <- strsplit ( as.character ( real.fs [ i ] ) , "\\." ) [[ 1 ]] [ 2 ]
		s.lab <- strsplit ( as.character ( s ) , "\\." ) [[ 1 ]][ 2 ]
		mut.density <- get.mut.density ( file = paste ( "mssel_f20" , f.lab , s.lab , N , ".out" , sep = "" ) )
		S.over.f [[ i ]] <- mut.density
}

##################
##### Seg Site #####
##################


R <- 1:200
r <- R / ( 4*N )
n = 20
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
std.time <- numeric()
P_NR <- exp ( - r * log ( 2 * N * s ) / s )
for ( j in 1 : length ( r ) ) {
	binom <- dbinom ( 0 : n , n , P_NR [ j ] )
	std.time [ j ] <-  sum ( c ( tail ( tot.times , 1 ) , rev ( tot.times ) ) * binom ) /  tail ( tot.times , 1 )
}







pdf ("../Paper/Paper_Figures/pi_and_S_density.pdf" , width = 7.874 , height = 5  )
par ( mfrow = c ( 1 , 2 ) )
plot ( c ( 0 , 200 ) , c ( 0 , 1 ) , type = "n" , xlab = "4NR" , ylab = "" , cex.lab = 1.5 , bty = "n" )
mtext ( expression ( pi [ R ] / pi [ 0 ] ) , side = 2 , line = 2 , cex = 1.5 )
for ( i in 1 : length ( real.fs ) ) {
		#run.ms.f(f.index)
		# mut.density <- get.mut.density ( file = paste ( "mssel_f2" , f.index , ".out" , sep = "" ) )
		# pi.over.f [[ f.index ]] <- mut.density
		lines ( pi.over.f [[ i ]]$x  , pi.over.f [[ i ]]$y/(1000*20) , col = i )
}
R <- 1:200
r <- R / ( 4*N )
for ( i in 2 : length ( real.fs ) ) {
		
	lines ( R , 1 - exp ( - 2 * r * Tsf [ i ]  ) / ( 1 + R * real.fs [ i ] * ( 1 - real.fs [ i ] ) ) , lty = 2 , lwd = 2 , col = i )
	
}
lines ( R , 1 - exp ( - 2 * r * log ( 2 * N * s ) / s ) , lty = 3 , lwd = 2 , col = 1 )
legend ( x = 125 , y = 0.35 , legend = c ( "1/2N" , real.fs [ 2 : length ( real.fs ) ] ) , col = 1 : 8 , lty = 1 , lwd = 1.5 , bty = "n")
legend ( x =20 , y = 0.2 , legend = c ( "Standard" , "Ours" , "Simulation") , lty = c ( 3 , 2 , 1 ) , lwd = 1.5 , bty = "n")

plot ( c ( 0 , 200 ) , c ( 0 , 1 ) , type = "n" , xlab = "4NR" , ylab = "" , cex.lab = 1.5 , bty = "n" )
mtext ( expression ( S [ R ] / S [ 0 ] ) , side = 2 , line = 2 , cex = 1.5 )
for ( i in 1 : length ( real.fs ) ) {
		# f.lab <- strsplit ( as.character ( real.fs [ i ] ) , "\\." ) [[ 1 ]] [ 2 ]
		# s.lab <- strsplit ( as.character ( s ) , "\\." ) [[ 1 ]][ 2 ]
		# mut.density <- get.mut.density ( file = paste ( "mssel_f20" , f.lab , s.lab , N , ".out" , sep = "" ) )
		# S.over.f [[ i ]] <- mut.density
		lines ( S.over.f [[ i ]]$x  , S.over.f [[ i ]]$y/(1000*20*sum ( 1/(1:19))) , col = i )
}
for ( i in 2 : length ( real.fs ) ) {
	lines ( R , my.time [[ i ]] , lty = 2 , lwd = 2 , col = i )
}
dev.off()







##### num haps
# EwensHaps40Sim <- StructuredCoalescentSweep ( N = 10000 , s = 0.01 , dominance = FALSE , f = 0.025 , reps = 1000 , n.tips = 40 , r = 10^-8 , sim.distance = 0.02 , interval.width = 1000 , no.sweep = TRUE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 2 )
# save ( EwensHaps40Sim , file = "Sims/EwensHaps40Sim.Robj" )
# load ( file = "Sims/EwensHaps40Sim.Robj" )


EwensHaps10Sim <- StructuredCoalescentSweep ( N = 10000 , s = 0.01 , dominance = FALSE , f = 0.025 , reps = 1000 , n.tips = 10 , r = 10^-8 , sim.distance = 0.02 , interval.width = 1000 , no.sweep = TRUE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 2 )
save ( EwensHaps10Sim , file = "Sims/EwensHaps10Sim.Robj" )

# # 
# ##### pi and seg sites
# real.fs <- c ( 1 / 20000 , 0.001 , 0.01 , 0.05 , 0.1 )

# T_f <- numeric ( length ( real.fs ) )
# for ( i in 1 : length ( real.fs ) ) {
	
	# T_f [ i ] <- log ( (2*N -1 ) * ( 1 - real.fs [ i ] ) / real.fs [ i ] ) / s
	
# }







##############################
###### Frequency Spectrum ######
##############################
source (  "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Scripts/freq_spectrum_standing_sweep_coal.R")
library ( RColorBrewer)


if ( FALSE ) {


### without sweep portion
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f02.specs.range.nosweep <- list ()
my.rs <- seq ( 0.0000001 , 0.006 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f02.specs.range.nosweep [[ i ]] <- expected.freq.times.standing ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.02 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f02.specs.range.nosweep , file = "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef02nosweepN10000n10.Robj")



## f = 0.05 
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f05.specs.range.nosweep <- list ()
my.rs <- seq ( 0.0000001 , 0.006 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f05.specs.range.nosweep [[ i ]] <- expected.freq.times.standing ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.05 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f05.specs.range.nosweep , file = "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef05nosweepN10000n10.Robj")



## f = 0.07
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f07.specs.range.nosweep <- list ()
my.rs <- seq ( 0.0000001 , 0.006 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f07.specs.range.nosweep [[ i ]] <- expected.freq.times.standing ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.07 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f07.specs.range.nosweep , file = "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef07nosweepN10000n10.Robj")

}

load ( "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef02nosweepN10000n10.Robj" )
load ( "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef05nosweepN10000n10.Robj" )
load ( "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef07nosweepN10000n10.Robj" )

load ( "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/sim.freq.spec.list.stoch.freq.nosweep.condloss.Rdata" )
# sims.f005.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[1]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f01.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[1]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f02.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[2]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f03.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[3]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f05.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[4]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f07.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[5]] , function ( x ) x [[5]]) ) [ , 1:9 ]



neutral <- 1 / ( 1:9 ) / sum ( 1 / ( 1:9 ) )
sim.rs <- c ( seq ( 0 , 0.0009 , length.out = 10 ) , seq ( 0.001 , 0.002 , length.out = 5 ) , c ( 0.0025 , 0.003 , 0.004))

# my.specs.f02.nosweep.relative <- t ( apply ( my.specs.f02.nosweep , 1 , function ( x )  x / neutral ) )
sims.f02.freq.spec.relative <- t ( apply ( sims.f02.freq.spec , 1 , function ( x ) x / neutral ) )
# my.specs.f05.nosweep.relative <- t ( apply ( my.specs.f05.nosweep , 1 , function ( x )  x / neutral ) )
sims.f05.freq.spec.relative <- t ( apply ( sims.f05.freq.spec , 1 , function ( x ) x / neutral ) )
# my.specs.f07.nosweep.relative <- t ( apply ( my.specs.f07.nosweep , 1 , function ( x )  x / neutral ) )
sims.f07.freq.spec.relative <- t ( apply ( sims.f07.freq.spec , 1 , function ( x ) x / neutral ) )



# adj
#my.specs.f01.nosweep.adj <- do.call ( rbind , lapply ( f01.specs.range.nosweep , function ( x ) x [[ 2 ]] ) )
my.specs.f02.nosweep.adj <- do.call ( rbind , lapply ( f02.specs.range.nosweep , function ( x ) x [[ 2 ]] ) )
#my.specs.f03.nosweep.adj <- do.call ( rbind , lapply ( f03.specs.range.nosweep , function ( x ) x [[ 2 ]] ) )
my.specs.f05.nosweep.adj <- do.call ( rbind , lapply ( f05.specs.range.nosweep , function ( x ) x [[ 2 ]] ) )
my.specs.f07.nosweep.adj <- do.call ( rbind , lapply ( f07.specs.range.nosweep , function ( x ) x [[ 2 ]] ) )


# sims.f02.freq.spec.relative <- t ( apply ( sims.f02.freq.spec , 1 , function ( x ) x / neutral ) )
# sims.f05.freq.spec.relative <- t ( apply ( sims.f05.freq.spec , 1 , function ( x ) x / neutral ) )
# sims.f07.freq.spec.relative <- t ( apply ( sims.f07.freq.spec , 1 , function ( x ) x / neutral ) )



my.specs.f02.nosweep.relative.adj <- t ( apply ( my.specs.f02.nosweep.adj , 1 , function ( x )  x / neutral ) )
my.specs.f05.nosweep.relative.adj <- t ( apply ( my.specs.f05.nosweep.adj , 1 , function ( x )  x / neutral ) )
my.specs.f07.nosweep.relative.adj <- t ( apply ( my.specs.f07.nosweep.adj , 1 , function ( x )  x / neutral ) )



## with sweep

load("/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/sim.freq.spec.list.stoch.freq.with.sweep.condloss.Rdata")

stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f02.specs.range.w.sweep <- list ()
my.rs <- seq ( 0.0000001 , 0.006 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f02.specs.range.w.sweep [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.02 , s = 0.05 ,  my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f02.specs.range.w.sweep , file = "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef02wsweepN10000n10.Robj" )



## f = 0.05 
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f05.specs.range.w.sweep <- list ()
my.rs <- seq ( 0.0000001 , 0.006 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f05.specs.range.w.sweep [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.05 , s = 0.05 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f05.specs.range.w.sweep , file = "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef05wsweepN10000n10.Robj")


## f = 0.07
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f07.specs.range.w.sweep <- list ()
my.rs <- seq ( 0.0000001 , 0.006 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f07.specs.range.w.sweep [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.07 , s = 0.05 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f07.specs.range.w.sweep , file = "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef07wsweepN10000n10.Robj" )


my.specs.f02.w.sweep.adj <- do.call ( rbind , f02.specs.range.w.sweep )
my.specs.f05.w.sweep.adj <- do.call ( rbind , f05.specs.range.w.sweep )
my.specs.f07.w.sweep.adj <- do.call ( rbind , f07.specs.range.w.sweep )

my.specs.f02.w.sweep.relative.adj <- t ( apply ( my.specs.f02.w.sweep.adj , 1 , function ( x )  x / neutral ) )
my.specs.f05.w.sweep.relative.adj <- t ( apply ( my.specs.f05.w.sweep.adj , 1 , function ( x )  x / neutral ) )
my.specs.f07.w.sweep.relative.adj <- t ( apply ( my.specs.f07.w.sweep.adj , 1 , function ( x )  x / neutral ) )






load ( "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/sim.freq.spec.list.stoch.freq.with.sweep.condloss.Rdata" )
#sims.f01.freq.spec.w.sweep <- do.call ( rbind , lapply ( sim.freq.spec.list[[1]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f02.freq.spec.w.sweep <- do.call ( rbind , lapply ( sim.freq.spec.list[[1]] , function ( x ) x [[5]]) ) [ , 1:9 ]
#sims.f03.freq.spec.w.sweep <- do.call ( rbind , lapply ( sim.freq.spec.list[[3]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f05.freq.spec.w.sweep <- do.call ( rbind , lapply ( sim.freq.spec.list[[2]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f07.freq.spec.w.sweep <- do.call ( rbind , lapply ( sim.freq.spec.list[[3]] , function ( x ) x [[5]]) ) [ , 1:9 ]


neutral <- 1 / ( 1:9 ) / sum ( 1 / ( 1:9 ) )
sim.rs <- c ( seq ( 0 , 0.0009 , length.out = 10 ) , seq ( 0.001 , 0.002 , length.out = 5 ) , c ( 0.0025 , 0.003 , 0.004))

sims.f02.freq.spec.w.sweep.relative <- t ( apply ( sims.f02.freq.spec.w.sweep , 1 , function ( x ) x / neutral ) )
sims.f05.freq.spec.w.sweep.relative <- t ( apply ( sims.f05.freq.spec.w.sweep , 1 , function ( x ) x / neutral ) )
sims.f07.freq.spec.w.sweep.relative <- t ( apply ( sims.f07.freq.spec.w.sweep , 1 , function ( x ) x / neutral ) )




# six panel theory and sims
pdf ( "Paper/Paper_Figures/freq_spec_nosweep_logfold_sixpanel_020507.pdf" , height = 10 , width = 16)
par ( mfrow = c ( 2, 3 ) )
matplot ( my.rs/2 , log ( my.specs.f02.nosweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -0.8 , 2 ))
matplot ( sim.rs , log ( sims.f02.freq.spec.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.02" , side = 3 )


matplot ( my.rs/2 , log ( my.specs.f05.nosweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n" , ylim = c ( -0.8 , 2 ))
matplot ( sim.rs , log ( sims.f05.freq.spec.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.05" , side = 3 )


matplot ( my.rs/2 , log ( my.specs.f07.nosweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n" , ylim = c ( -0.8 , 2 ))
matplot ( sim.rs , log ( sims.f07.freq.spec.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
legend ( "topright" , legend =  seq ( 1, 9 ) , col = brewer.pal ( 9 , "Set1" )  , pch = 20 , bty = "n" , cex = 1.5 )
mtext ( "f = 0.07" , side = 3 )



matplot ( my.rs/2 , log ( my.specs.f02.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -1.3 , 2.4 ))
matplot ( sim.rs , log ( sims.f02.freq.spec.w.sweep.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
#mtext ( "f = 0.02" , side = 3 )


matplot ( my.rs/2 , log ( my.specs.f05.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -1.3 , 2.4 ))
matplot ( sim.rs , log ( sims.f05.freq.spec.w.sweep.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
#mtext ( "f = 0.02" , side = 3 )

matplot ( my.rs/2 , log ( my.specs.f07.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -1.3 , 2.4 ))
matplot ( sim.rs , log ( sims.f07.freq.spec.w.sweep.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
#mtext ( "f = 0.02" , side = 3 )
dev.off()






pdf ( "Paper/Paper_Figures/freq_spec_nosweep_logfold_sixpanel_020507_theory_only.pdf" , height = 10 , width = 16)
par ( mfrow = c ( 2, 3 ) )
matplot ( my.rs/2 , log ( my.specs.f02.nosweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -0.8 , 2 ))
#matplot ( sim.rs , log ( sims.f02.freq.spec.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.02" , side = 3 )


matplot ( my.rs/2 , log ( my.specs.f05.nosweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n" , ylim = c ( -0.8 , 2 ))
#matplot ( sim.rs , log ( sims.f05.freq.spec.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.05" , side = 3 )


matplot ( my.rs/2 , log ( my.specs.f07.nosweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n" , ylim = c ( -0.8 , 2 ))
legend ( "topright" , legend =  seq ( 1, 9 ) , col = brewer.pal ( 9 , "Set1" )  , pch = 20 , bty = "n" , cex = 1.5 )
#matplot ( sim.rs , log ( sims.f07.freq.spec.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.07" , side = 3 )



matplot ( my.rs/2 , log ( my.specs.f02.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -1.3 , 2.4 ))
#matplot ( sim.rs , log ( sims.f02.freq.spec.w.sweep.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
#mtext ( "f = 0.02" , side = 3 )


matplot ( my.rs/2 , log ( my.specs.f05.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -1.3 , 2.4 ))
#matplot ( sim.rs , log ( sims.f05.freq.spec.w.sweep.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
#mtext ( "f = 0.02" , side = 3 )

matplot ( my.rs/2 , log ( my.specs.f07.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -1.3 , 2.4 ))
#matplot ( sim.rs , log ( sims.f07.freq.spec.w.sweep.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
#mtext ( "f = 0.02" , side = 3 )
dev.off()




pdf ( "Paper/Paper_Figures/freq_spec_nosweep_logfold_sixpanel_020507_sims_only.pdf" , height = 10 , width = 16)
par ( mfrow = c ( 2, 3 ) )
#matplot ( my.rs/2 , log ( my.specs.f02.nosweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -0.8 , 2 ))
matplot ( sim.rs , log ( sims.f02.freq.spec.relative , 2 ) , type = "b" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , ylim = c ( -0.8 , 2 ) , lty = 1 , lwd = 0.7 )
mtext ( "f = 0.02" , side = 3 )


#matplot ( my.rs/2 , log ( my.specs.f05.nosweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n" , ylim = c ( -0.8 , 2 ))
matplot ( sim.rs , log ( sims.f05.freq.spec.relative , 2 ) , type = "b" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , ylim = c ( -0.8 , 2 ) , lty = 1 , lwd = 0.7 )
mtext ( "f = 0.05" , side = 3 )


#matplot ( my.rs/2 , log ( my.specs.f07.nosweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n" , ylim = c ( -0.8 , 2 ))
matplot ( sim.rs , log ( sims.f07.freq.spec.relative , 2 ) , type = "b" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , ylim = c ( -0.8 , 2 ) , lty = 1 , lwd = 0.7 )
legend ( "topright" , legend =  seq ( 1, 9 ) , col = brewer.pal ( 9 , "Set1" )  , pch = 20 , bty = "n" , cex = 1.5 )
mtext ( "f = 0.07" , side = 3 )



#matplot ( my.rs/2 , log ( my.specs.f02.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -2 , 2 ))
matplot ( sim.rs , log ( sims.f02.freq.spec.w.sweep.relative , 2 ) , type = "b" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , ylim = c ( -2 , 2 ) , lty = 1 , lwd = 0.7 )
#mtext ( "f = 0.02" , side = 3 )


#matplot ( my.rs/2 , log ( my.specs.f05.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -2 , 2 ))
matplot ( sim.rs , log ( sims.f05.freq.spec.w.sweep.relative , 2 ) , type = "b" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , ylim = c ( -2 , 2 ) , lty = 1 , lwd = 0.7 )
#mtext ( "f = 0.02" , side = 3 )

#matplot ( my.rs/2 , log ( my.specs.f07.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -2 , 2 ))
matplot ( sim.rs , log ( sims.f07.freq.spec.w.sweep.relative , 2 ) , type = "b" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , ylim = c ( -2 , 2 ) , lty = 1 , lwd = 0.7 )
#mtext ( "f = 0.02" , side = 3 )
dev.off()







##### haplotype spectrum
directory="~/Documents/Academics/StandingSweeps"
# directory="~/Dropbox/Linked_selection_models/Soft_sweeps_coal/LinkedSelection/"

source (  paste(directory,"/Scripts/HapFreqSpecs.R",sep=""))

load ( "Sims/HapSims/one.side.hard.n100.denovo.s01.Robj" )
load ( "Sims/HapSims/one.side.standing.n100.f05.s01.Robj" )
load ( "Sims/HapSims/one.side.soft.n100.k3.s01.Robj" )
soft.sweep <- blah
load ( "Sims/HapSims/neutral.n100.Robj" )
## neutral <- Reduce ( "+" , blah ) / length ( blah)
soft.haps <- HapFreqs ( soft.sweep )
standing.haps <- HapFreqs ( standing.sweep [[ 2 ]] )
hard.haps <- HapFreqs ( hard.sweep [[ 2 ]] )
neutral.haps <- HapFreqs ( neutral [[ 2 ]] )

library("colorRamps")
coop.cols<-matlab.like(201)
 ramp1<-coop.cols[102:201]
 ramp2<- coop.cols[1:100] #rev(coop.cols[1:100]) 
 bland<-coop.cols[101]
 
 
upper.breaks <- c(2^(seq(log2(1.05),log2(2),length=100)),10^6)  #seq ( 1.05 , 2 ,length.out = 101)
lower.breaks<- c(0,2^seq(log2(1/2),log2(1/1.05),length=100))   # seq ( 0 , 0.95 ,length.out = 101)  
 
pdf ( "Paper/Paper_Figures/HapFreqRatiosCondExist.pdf" , height = 10 , width = 7.874 )

layout(matrix( c(1,1,2,3,4,5,6,7) , nrow=4,ncol=2,byrow=TRUE), heights=c(0.5,3,3,3))  

my.z<-c(lower.breaks,1,upper.breaks )

par(mar=c(2,10,1,10))
image(x= c(.00001,(my.z[-1])), z=cbind(my.z,my.z),col=c(ramp2,bland,bland,ramp1),breaks=c(lower.breaks,1,upper.breaks) ,xlim=(c(.5,2)),log="x",axes=FALSE)
axis(side=1, at=c(0.5,.75,1,1.25,1.5),labels=c("<0.5","0.75","1","0.25",">0.5"))

#par ( mfrow = c ( 3,2))
par(mar=c(3,3.5,1.5,1.2))

image ( t ( apply ( standing.haps [[ 2 ]] [ 1:50 , ] / hard.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = c ( bland , bland ) , breaks = seq ( 1/1.05 , 1.05 ,length.out = 3) ,xaxt = "n" , main = expression ( h[i]^stand/h[i]^hard ) , yaxt = "n" , ylab = "" )
axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
axis ( 2 , c ( 1 , seq ( 20 , 100 , length.out = 5 ) )/100 , labels = c ( seq ( 50 , 10 , length.out = 5 ) , 1 ) )
#mtext ( "Window Size (cM)" , side = 1 , line = 2.3 , cex = 0.8)
mtext ( expression ( h[i]) , side = 2 , line = 2 , cex = 0.8 )
image ( t ( apply ( standing.haps [[ 2 ]] [ 1:50 , ] / hard.haps [[ 2 ]] [ 1:50 , ] , 2 , rev) ) ,  col = ramp1 , breaks = upper.breaks , add = T )
image ( t ( apply ( standing.haps [[ 2 ]] [ 1:50 , ] / hard.haps [[ 2 ]] [ 1:50 , ] , 2 , rev) ) ,  col = ramp2 , breaks = lower.breaks, add = T )


image ( t ( apply ( standing.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = seq ( 1/1.05, 1.05 ,length.out=2), col = bland , xaxt = "n" , main = expression ( h[i]^stand/h[i]^neut ) ,yaxt = "n", ylab = "")
axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
axis ( 2 , c ( 1 , seq ( 20 , 100 , length.out = 5 ) )/100 , labels = c ( seq ( 50 , 10 , length.out = 5 ) , 1 ) )
#mtext ( "Window Size (cM)" , side = 1 , line = 2.3 , cex = 0.8)
#mtext ( expression ( h[i]) , side = 2 , line = 2 , cex = 0.8 )
image ( t ( apply ( standing.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = upper.breaks, col = ramp1,add = T )
image ( t ( apply ( standing.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = lower.breaks, col = ramp2, add =T )


image ( t ( apply ( soft.haps [[ 2 ]] [ 1 : 50 , ] / hard.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = c ( bland , bland ) , breaks = seq ( 1/1.05 , 1.05 ,length.out = 3) , xaxt = "n", main = expression ( h[i]^soft/h[i]^hard ) ,yaxt = "n" , ylab = "")
axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
axis ( 2 , c ( 1 , seq ( 20 , 100 , length.out = 5 ) )/100 , labels = c ( seq ( 50 , 10 , length.out = 5 ) , 1 ) )
#mtext ( "Window Size (cM)" , side = 1 , line = 2.3 , cex = 0.8)
mtext ( expression ( h[i]) , side = 2 , line = 2 , cex = 0.8 )
image ( t ( apply ( soft.haps [[ 2 ]] [ 1 : 50 , ] / hard.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = ramp1 , breaks = upper.breaks , add = T )
image ( t ( apply ( soft.haps [[ 2 ]] [ 1 : 50 , ] / hard.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = ramp2 , breaks =lower.breaks , add = T )

image ( t ( apply ( soft.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = seq ( 1/1.05, 1.05 ,length.out=2), col = bland , xaxt = "n", main = expression ( h[i]^soft/h[i]^neut ) ,yaxt = "n" , ylab = "")
axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
axis ( 2 , c ( 1 , seq ( 20 , 100 , length.out = 5 ) )/100 , labels = c ( seq ( 50 , 10 , length.out = 5 ) , 1 ) )
#mtext ( "Window Size (cM)" , side = 1 , line = 2.3 , cex = 0.8)
#mtext ( expression ( h[i]) , side = 2 , line = 2 , cex = 0.8 )
image ( t ( apply ( soft.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = upper.breaks, col = ramp1,add = T )
image ( t ( apply ( soft.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = lower.breaks, col = ramp2, add =T )


image ( t ( apply ( standing.haps [[ 2 ]] [ 1 : 50 , ] / soft.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = c ( bland , bland ) , breaks = seq ( 1/1.05 , 1.05 ,length.out = 3) , xaxt = "n", main = expression ( h[i]^stand/h[i]^soft ) , yaxt = "n" , ylab = "")
axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
axis ( 2 , c ( 1 , seq ( 20 , 100 , length.out = 5 ) )/100 , labels = c ( seq ( 50 , 10 , length.out = 5 ) , 1 ) )
mtext ( "Window Size (cM)" , side = 1 , line = 2 , cex = 0.8)
mtext ( expression ( h[i]) , side = 2 , line = 2 , cex = 0.8 )
image ( t ( apply ( standing.haps [[ 2 ]] [ 1 : 50 , ] / soft.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = ramp1 , breaks = upper.breaks , add = T )
image ( t ( apply ( standing.haps [[ 2 ]] [ 1 : 50 , ] / soft.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = ramp2 , breaks = lower.breaks, add = T )

image ( t ( apply ( hard.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = seq ( 1/1.05, 1.05 ,length.out=2), col = bland , xaxt = "n" , main = expression ( h[i]^hard/h[i]^neut ) , yaxt = "n" , ylab ="")
axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
axis ( 2 , c ( 1 , seq ( 20 , 100 , length.out = 5 ) )/100 , labels = c ( seq ( 50 , 10 , length.out = 5 ) , 1 ) )
mtext ( "Window Size (cM)" , side = 1 , line = 2 , cex = 0.8)
#mtext ( expression ( h[i]) , side = 2 , line = 2 , cex = 0.8 )
image ( t ( apply ( hard.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = upper.breaks, col = ramp1,add = T )
image ( t ( apply ( hard.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = lower.breaks, col = ramp2, add =T )
dev.off()




# pdf ( "Paper/Paper_Figures/HapFreqRatiosCondExist.pdf" , height = 10 , width = 7.874 )
# par ( mfrow = c ( 3,2))
# image ( t ( apply ( standing.haps [[ 2 ]] [ 1:50 , ] / hard.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = c ( "black" , "black" ) , breaks = seq ( 0.95 , 1.05 ,length.out = 3) ,xaxt = "n" , main = expression ( h[i]^stand/h[i]^hard ) , yaxt = "n" , ylab = "" )
# axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
# axis ( 2 , c ( 1 , seq ( 20 , 100 , length.out = 5 ) )/100 , labels = c ( seq ( 50 , 10 , length.out = 5 ) , 1 ) )
# mtext ( "Window Size (cM)" , side = 1 , line = 2.3 , cex = 0.8)
# mtext ( expression ( h[i]) , side = 2 , line = 2 , cex = 0.8 )
# image ( t ( apply ( standing.haps [[ 2 ]] [ 1:50 , ] / hard.haps [[ 2 ]] [ 1:50 , ] , 2 , rev) ) ,  col = heat.colors ( 100 ) , breaks = seq ( 1.05 , 2 ,length.out = 101) , add = T )
# image ( t ( apply ( standing.haps [[ 2 ]] [ 1:50 , ] / hard.haps [[ 2 ]] [ 1:50 , ] , 2 , rev) ) ,  col = cm.colors ( 100 ) , breaks = seq ( 0 , 0.95 ,length.out = 101) , add = T )


# image ( t ( apply ( standing.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = seq ( 0.95, 1.05 ,length.out=2), col = "black" , xaxt = "n" , main = expression ( h[i]^stand/h[i]^neut ) ,yaxt = "n", ylab = expression ( h[i]))
# axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
# axis ( 2 , c ( 1 , seq ( 20 , 100 , length.out = 5 ) )/100 , labels = c ( seq ( 50 , 10 , length.out = 5 ) , 1 ) )
# mtext ( "Window Size (cM)" , side = 1 , line = 2.3 , cex = 0.8)
# mtext ( expression ( h[i]) , side = 2 , line = 2 , cex = 0.8 )
# image ( t ( apply ( standing.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = seq ( 1.05 , 2 ,length.out=101), col = heat.colors ( 100 ),add = T )
# image ( t ( apply ( standing.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = seq ( 0 , 0.95 ,length.out=101), col = cm.colors ( 100 ), add =T )


# image ( t ( apply ( soft.haps [[ 2 ]] [ 1 : 50 , ] / hard.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = c ( "black" , "black" ) , breaks = seq ( 0.95 , 1.05 ,length.out = 3) , xaxt = "n", main = expression ( h[i]^soft/h[i]^hard ) ,yaxt = "n" , ylab = expression ( h[i]))
# axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
# axis ( 2 , c ( 1 , seq ( 20 , 100 , length.out = 5 ) )/100 , labels = c ( seq ( 50 , 10 , length.out = 5 ) , 1 ) )
# mtext ( "Window Size (cM)" , side = 1 , line = 2.3 , cex = 0.8)
# mtext ( expression ( h[i]) , side = 2 , line = 2 , cex = 0.8 )
# image ( t ( apply ( soft.haps [[ 2 ]] [ 1 : 50 , ] / hard.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = heat.colors ( 100 ) , breaks = seq ( 1.05 , 2 ,length.out = 101) , add = T )
# image ( t ( apply ( soft.haps [[ 2 ]] [ 1 : 50 , ] / hard.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = cm.colors ( 100 ) , breaks = seq ( 0 , 0.95 ,length.out = 101) , add = T )

# image ( t ( apply ( soft.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = seq ( 0.95, 1.05 ,length.out=2), col = "black" , xaxt = "n", main = expression ( h[i]^soft/h[i]^neut ) ,yaxt = "n" , ylab = expression ( h[i]))
# axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
# axis ( 2 , c ( 1 , seq ( 20 , 100 , length.out = 5 ) )/100 , labels = c ( seq ( 50 , 10 , length.out = 5 ) , 1 ) )
# mtext ( "Window Size (cM)" , side = 1 , line = 2.3 , cex = 0.8)
# mtext ( expression ( h[i]) , side = 2 , line = 2 , cex = 0.8 )
# image ( t ( apply ( soft.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = seq ( 1.05 , 2 ,length.out=101), col = heat.colors ( 100 ),add = T )
# image ( t ( apply ( soft.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = seq ( 0 , 0.95 ,length.out=101), col = cm.colors ( 100 ), add =T )



# image ( t ( apply ( standing.haps [[ 2 ]] [ 1 : 50 , ] / soft.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = c ( "black" , "black" ) , breaks = seq ( 0.95 , 1.05 ,length.out = 3) , xaxt = "n", main = expression ( h[i]^stand/h[i]^soft ) , yaxt = "n" , ylab = expression ( h[i]))
# axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
# axis ( 2 , c ( 1 , seq ( 20 , 100 , length.out = 5 ) )/100 , labels = c ( seq ( 50 , 10 , length.out = 5 ) , 1 ) )
# mtext ( "Window Size (cM)" , side = 1 , line = 2.3 , cex = 0.8)
# mtext ( expression ( h[i]) , side = 2 , line = 2 , cex = 0.8 )
# image ( t ( apply ( standing.haps [[ 2 ]] [ 1 : 50 , ] / soft.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = heat.colors ( 100 ) , breaks = seq ( 1.05 , 2 ,length.out = 101) , add = T )
# image ( t ( apply ( standing.haps [[ 2 ]] [ 1 : 50 , ] / soft.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) ,  col = cm.colors ( 100 ) , breaks = seq ( 0 , 0.95 ,length.out = 101) , add = T )

# image ( t ( apply ( hard.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = seq ( 0.95, 1.05 ,length.out=2), col = "black" , xaxt = "n" , main = expression ( h[i]^hard/h[i]^neut ) , yaxt = "n" , ylab = expression ( h[i]))
# axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
# axis ( 2 , c ( 1 , seq ( 20 , 100 , length.out = 5 ) )/100 , labels = c ( seq ( 50 , 10 , length.out = 5 ) , 1 ) )
# mtext ( "Window Size (cM)" , side = 1 , line = 2.3 , cex = 0.8)
# mtext ( expression ( h[i]) , side = 2 , line = 2 , cex = 0.8 )
# image ( t ( apply ( hard.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = seq ( 1.05 , 2 ,length.out=101), col = heat.colors ( 100 ),add = T )
# image ( t ( apply ( hard.haps [[ 2 ]] [ 1 : 50 , ] / neutral.haps [[ 2 ]] [ 1 : 50 , ] , 2 , rev) ) , breaks = seq ( 0 , 0.95 ,length.out=101), col = cm.colors ( 100 ), add =T )
# dev.off()







###################
#### Supplement ####
###################


### effective s
N <- 10000
my.r <- seq ( 10^-8 , 0.005 , by = 10^-6)
my.fs <- c ( 0.005 , 0.01 , 0.02 , 0.03 , 0.04 , 0.05 )
s <- 0.01

pi.reductions <- list()
pi.reductions.approx <- list()
eff.s <- sapply ( my.fs , function ( x ) EffectiveS ( N = N , s = s , x ) )
for ( i in 1 : length ( my.fs ) ) {
	
	pi.reductions [[ i ]] <- 1 - exp ( - my.r * log ( 1/ my.fs [ i ] ) / s ) / ( 1 + 4 * N * my.r * my.fs [ i ] * ( 1 - my.fs [ i ] ) )
	pi.reductions.approx [[ i ]] <- 1 - exp ( - my.r * log ( 2 * N * eff.s [ i ] ) / eff.s [ i ] )

}



plot ( my.r , pi.reductions [[ 1 ]] , type = "l" , lwd = 2 )
lines ( my.r , pi.reductions.approx [[ 1 ]] , col = "red" , lty = 5 )

lines ( my.r , pi.reductions [[ 2 ]] , type = "l" , lwd = 2 )
lines ( my.r , pi.reductions.approx [[ 2 ]] , type = "l" , col = "red" , lty = 2 )


lines ( my.r , pi.reductions [[ 3 ]] , type = "l" , lwd = 2 )
lines ( my.r , pi.reductions.approx [[ 3 ]] , type = "l" , col = "red" , lty = 2 )



##### coal events before reach standing phase
temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.05 , dominance = FALSE , f = 0.05 , reps = 10000 , n.tips = 10 , r = 10^-8 , sim.distance = 0.02 , interval.width = 1000 , no.sweep = F , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = FALSE , display.rep.count = TRUE ,   standing.haps = FALSE , time.factor = 1 )
probs <- colSums ( temp$coal.times-temp$sweep.start<0 )/10000





