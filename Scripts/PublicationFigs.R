#### plot of pi and segregating sites
setwd ( "~/Documents/Academics/StandingSweeps/Sims" )
source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R')
source ( "~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R")
options ( scipen = 400 )


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






#################
## extra functions ##
#################


# MyLogistic <- function ( x , N = 10000 , s  ) 1 / ( 2 * N  ) * exp(s * x ) / ( 1 + 1 / (2 * N  ) * ( exp(s * x )  - 1 )  )

MyLogistic <- function ( x , N = 10000 , s  ) 1 / ( 5 * N * s ) * exp(s * x ) / ( 1 + 1 / ( 5 * N * s  ) * ( exp ( s * x )  - 1 )  )


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

i <- 1
det.freqs <- 1 / ( 20000 * 0.01 )
while ( det.freqs [ i ] < ( 1 - f ) ) {	
	det.freqs [ i + 1 ] <- MyLogistic ( i , N = 10000 , s = s)
	i <- i + 1 
}
det.freqs <- c ( rep ( f , nrow ( freqs) - length ( det.freqs ) ) , rev ( 1 - det.freqs ) )
matplot ( freqs , type = "l" , lwd = 0.7 , col = "grey" , lty = 1 , bty = "n" , ylab = "Frequency" , xlab = "Generations" , xaxt = "n" )
lines ( det.freqs , lwd = 2 )

my.at <- length ( det.freqs ) - seq ( 0 , 3500 , 500)
axis ( 1 , at = my.at , labels = seq ( 0 , 3500 , 500) )


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




pdf ( "Figures/freq.spec.nosweep.logfold.threepanel.020507.pdf" , height = 5 , width = 16)
par ( mfrow = c ( 1, 3 ) )
matplot ( my.rs/2 , log ( my.specs.f02.nosweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -0.8 , 2 ))
matplot ( sim.rs , log ( sims.f02.freq.spec.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.02" , side = 3 )


matplot ( my.rs/2 , log ( my.specs.f05.nosweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n" , ylim = c ( -0.8 , 2 ))
matplot ( sim.rs , log ( sims.f05.freq.spec.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.05" , side = 3 )


matplot ( my.rs/2 , log ( my.specs.f07.nosweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n" , ylim = c ( -0.8 , 2 ))
matplot ( sim.rs , log ( sims.f07.freq.spec.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.07" , side = 3 )
dev.off()





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
save ( f02.specs.range.w.sweep , file = "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef02wsweepN10000n10.Robj")



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
pdf ( "Paper/Paper_Figures/freq.spec.nosweep.logfold.sixpanel.020507.pdf" , height = 10 , width = 16)
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



matplot ( my.rs/2 , log ( my.specs.f02.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -2 , 2 ))
matplot ( sim.rs , log ( sims.f02.freq.spec.w.sweep.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
#mtext ( "f = 0.02" , side = 3 )


matplot ( my.rs/2 , log ( my.specs.f05.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -2 , 2 ))
matplot ( sim.rs , log ( sims.f05.freq.spec.w.sweep.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
#mtext ( "f = 0.02" , side = 3 )

matplot ( my.rs/2 , log ( my.specs.f07.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -2 , 2 ))
matplot ( sim.rs , log ( sims.f07.freq.spec.w.sweep.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
#mtext ( "f = 0.02" , side = 3 )
dev.off()






pdf ( "Paper/Paper_Figures/freq.spec.nosweep.logfold.sixpanel.020507.theory.only.pdf" , height = 10 , width = 16)
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



matplot ( my.rs/2 , log ( my.specs.f02.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -4 , 2 ))
#matplot ( sim.rs , log ( sims.f02.freq.spec.w.sweep.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
#mtext ( "f = 0.02" , side = 3 )


matplot ( my.rs/2 , log ( my.specs.f05.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -4 , 2 ))
#matplot ( sim.rs , log ( sims.f05.freq.spec.w.sweep.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
#mtext ( "f = 0.02" , side = 3 )

matplot ( my.rs/2 , log ( my.specs.f07.w.sweep.relative.adj , 2 )  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = expression ( paste ( log[2] , "(Deviation from Neutral)" , sep = " " ) ) , xlab = "Genetic Distance" , bty = "n"  , ylim = c ( -4 , 2 ))
#matplot ( sim.rs , log ( sims.f07.freq.spec.w.sweep.relative , 2 ) , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
#mtext ( "f = 0.02" , side = 3 )
dev.off()




pdf ( "Paper/Paper_Figures/freq.spec.nosweep.logfold.sixpanel.020507.sims.only.pdf" , height = 10 , width = 16)
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










##### coal events before reach standing phase

temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.05 , dominance = FALSE , f = 0.05 , reps = 10000 , n.tips = 10 , r = 10^-8 , sim.distance = 0.02 , interval.width = 1000 , no.sweep = F , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = FALSE , display.rep.count = TRUE ,   standing.haps = FALSE , time.factor = 1 )
probs <- colSums ( temp$coal.times-temp$sweep.start<0 )/10000





