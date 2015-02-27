#### plot of pi and segregating sites
setwd ( "~/Documents/Academics/StandingSweeps" )
source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R')
source ( "~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R")


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
source (  "Scripts/freq_spectrum_standing_sweep_coal.R")
library ( RColorBrewer)


### without sweep portion

#load ( "Sims/sim.freq.spec.list.stoch.freq.no.sweep.nocondloss.Rdata" )
load ( "Sims/sim.freq.spec.list.stoch.freq.no.sweep.condloss.Rdata" )



## f = 0.005
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f005.specs.range.nosweep <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f005.specs.range.nosweep [[ i ]] <- expected.freq.times.standing ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.005 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f005.specs.range.nosweep , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef005nosweepN10000n10.Robj")


## f = 0.01 
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f01.specs.range.nosweep <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f01.specs.range.nosweep [[ i ]] <- expected.freq.times.standing ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.01 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f01.specs.range.nosweep , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef01nosweepN10000n10.Robj")



## f = 0.03
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f03.specs.range.nosweep <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f03.specs.range.nosweep [[ i ]] <- expected.freq.times.standing ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.03 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f03.specs.range.nosweep , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef03nosweepN10000n10.Robj" )




## f = 0.05 
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f05.specs.range.nosweep <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f05.specs.range.nosweep [[ i ]] <- expected.freq.times.standing ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.05 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f05.specs.range.nosweep , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef05nosweepN10000n10.Robj")



## f = 0.07
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f07.specs.range.nosweep <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f07.specs.range.nosweep [[ i ]] <- expected.freq.times.standing ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.07 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f07.specs.range.nosweep , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef07nosweepN10000n10.Robj")


sims.f005.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[1]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f01.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[2]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f03.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[3]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f05.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[4]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f07.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[5]] , function ( x ) x [[5]]) ) [ , 1:9 ]


my.specs.f005.nosweep <- do.call ( rbind , f005.specs.range.nosweep )
my.specs.f01.nosweep <- do.call ( rbind , f01.specs.range.nosweep )
my.specs.f03.nosweep <- do.call ( rbind , f03.specs.range.nosweep )
my.specs.f05.nosweep <- do.call ( rbind , f05.specs.range.nosweep )
my.specs.f07.nosweep <- do.call ( rbind , f07.specs.range.nosweep )



sim.rs <- seq ( 0 , 0.002 , length.out = 5 )

matplot ( my.rs , my.specs.f01.nosweep  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = "Frequency" , xlab = "Genetic Distance" , bty = "n" )
matplot ( sim.rs , sims.f01.freq.spec , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.01" , side = 3 )











###### sweeps from a few different freqs
load("/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/sim.freq.spec.list.s01.Rdata")
load("/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/sim.freq.spec.list.s05.Rdata")

## denovo
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
denovo.specs.range.s01 <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	denovo.specs.range.s01 [[ i ]] <- freq.spec.de.novo  ( nsam = 10 , N = 10000 , r = my.rs [ i ] , s = 0.01 )
	message ( i )
}
save ( denovo.specs.range.s01 , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.range.denovo.s01N10000n10.Robj")



## f = 0.005
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f005.specs.range.s01 <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f005.specs.range.s01 [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.005 , s = 0.01 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f005.specs.range.s01 , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef005s01N10000n10.Robj")


## f = 0.01 
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f01.specs.range.s01 <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f01.specs.range.s01 [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.01 , s = 0.01 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f01.specs.range.s01 , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef01s01N10000n10.Robj")



## f = 0.03
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f03.specs.range.s01 <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f03.specs.range.s01 [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.03 , s = 0.01 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f03.specs.range.s01 , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef03s01N10000n10.Robj" )




## f = 0.05 
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f05.specs.range.s01 <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f05.specs.range.s01 [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.05 , s = 0.01 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f05.specs.range.s01 , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef05s01N10000n10.Robj")



## f = 0.07
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f07.specs.range.s01 <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f07.specs.range.s01 [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.07 , s = 0.01 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f07.specs.range.s01 , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef07s01N10000n10.Robj")


load ( "Paper/Paper_Figures/Data_and_Robjs/f.specs.range.denovo.s01N10000n10.Robj" )
load ( "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef005s01N10000n10.Robj" )
load ( "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef01s01N10000n10.Robj" )
load ( "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef03s01N10000n10.Robj" )
load ( "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef05s01N10000n10.Robj" )
load ( "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef07s01N10000n10.Robj" )

my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )

my.f.specs.denovo.s01 <- do.call ( rbind , denovo.specs.range.s01 )
my.specs.f005.s01 <- do.call ( rbind , f005.specs.range.s01 )
my.specs.f01.s01 <- do.call ( rbind , f01.specs.range.s01 )
my.specs.f03.s01 <- do.call ( rbind , f03.specs.range.s01 )
my.specs.f05.s01 <- do.call ( rbind , f05.specs.range.s01 )
my.specs.f07.s01 <- do.call ( rbind , f07.specs.range.s01 )


sim.rs <- c ( seq ( 0 , 0.001 , length.out = 11 ) , 0.0015 , 0.002 , 0.003 )

sims.denovo.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[1]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f005.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[2]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f01.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[3]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f03.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[4]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f05.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[5]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f07.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[6]] , function ( x ) x [[5]]) ) [ , 1:9 ]


library ( RColorBrewer )
par ( mfrow = c ( 2, 3 ) )
matplot ( my.rs , my.f.specs.denovo.s01  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = "Frequency" , xlab = "Genetic Distance" , bty = "n"  )
matplot ( sim.rs , sims.denovo.freq.spec , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 1/2N" , side = 3 )
legend ( "topright" , legend = 1:9 , lty = 1 , col = brewer.pal ( 9 , "Set1" )  , bty = "n" )

matplot ( my.rs , my.specs.f005.s01  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = "Frequency" , xlab = "Genetic Distance" , bty = "n" )
matplot ( sim.rs , sims.f005.freq.spec , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.005" , side = 3 )


matplot ( my.rs , my.specs.f01.s01  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = "Frequency" , xlab = "Genetic Distance" , bty = "n" )
matplot ( sim.rs , sims.f01.freq.spec , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.01" , side = 3 )


matplot ( my.rs , my.specs.f03.s01  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = "Frequency" , xlab = "Genetic Distance" , bty = "n" , ylim = c ( 0 , 0.8))
matplot ( sim.rs , sims.f03.freq.spec , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.03" , side = 3 )


matplot ( my.rs , my.specs.f05.s01  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = "Frequency" , xlab = "Genetic Distance" , bty = "n" , ylim = c ( 0 , 0.8 ) )
matplot ( sim.rs , sims.f05.freq.spec , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.05" , side = 3 )


matplot ( my.rs , my.specs.f07.s01  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = "Frequency" , xlab = "Genetic Distance" , bty = "n" )
matplot ( sim.rs , sims.f07.freq.spec , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.07" , side = 3 )







denovo.specs.range.s05 <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	denovo.specs.range.s05 [[ i ]] <- freq.spec.de.novo  ( nsam = 10 , N = 10000 , r = my.rs [ i ] , s = 0.05 )
	message ( i )
}
save ( denovo.specs.range.s05 , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.range.denovo.s05N10000n10.Robj")



## f = 0.005
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f005.specs.range.s05 <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f005.specs.range.s05 [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.005 , s = 0.05 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f005.specs.range.s05 , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef005s05N10000n10.Robj")


## f = 0.01 
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f01.specs.range.s05 <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f01.specs.range.s05 [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.01 , s = 0.05 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f01.specs.range.s05 , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef01s05N10000n10.Robj")



## f = 0.03
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f03.specs.range.s05 <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f03.specs.range.s05 [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.03 , s = 0.05 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f03.specs.range.s05 , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef03s05N10000n10.Robj" )




## f = 0.05 
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f05.specs.range.s05 <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f05.specs.range.s05 [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.05 , s = 0.05 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f05.specs.range.s05 , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef05s05N10000n10.Robj")



## f = 0.07
stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f07.specs.range.s05 <- list ()
my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
for ( i in seq_along ( my.rs ) ) {

	f07.specs.range.s05 [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.07 , s = 0.05 , my.StirlingNumbers = stirlings )
	message ( i )
}
save ( f07.specs.range.s05 , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef07s05N10000n10.Robj")

load ( "Paper/Paper_Figures/Data_and_Robjs/f.specs.range.denovo.s05N10000n10.Robj" )
f.denovo.specs.range.s05 <- f.specs.range
load ( "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef005s05N10000n10.Robj" )
load ( "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef01s05N10000n10.Robj" )
load ( "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef03s05N10000n10.Robj" )
load ( "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef05s05N10000n10.Robj" )
load ( "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef07s05N10000n10.Robj" )


my.specs.denovo.s05 <- do.call ( rbind , f.denovo.specs.range.s05 )
my.specs.f005.s05 <- do.call ( rbind , f005.specs.range.s05 )
my.specs.f01.s05 <- do.call ( rbind , f01.specs.range.s05 )
my.specs.f03.s05 <- do.call ( rbind , f03.specs.range.s05 )
my.specs.f05.s05 <- do.call ( rbind , f05.specs.range.s05 )
my.specs.f07.s05 <- do.call ( rbind , f07.specs.range.s05 )



load("/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/sim.freq.spec.list.s05.Rdata")
sims.denovo.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[1]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f005.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[2]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f01.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[3]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f03.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[4]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f05.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[5]] , function ( x ) x [[5]]) ) [ , 1:9 ]
sims.f07.freq.spec <- do.call ( rbind , lapply ( sim.freq.spec.list[[6]] , function ( x ) x [[5]]) ) [ , 1:9 ]



my.rs <- seq ( 0.0000001 , 0.003 , by = 0.00001 )
sim.rs <- c ( seq ( 0 , 0.001 , length.out = 11 ) , 0.0015 , 0.002 , 0.003 , 0.004 , 0.005 )

par ( mfrow = c ( 2, 3 ) )
matplot ( my.rs , my.specs.denovo.s05  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = "Frequency" , xlab = "Genetic Distance" , bty = "n"  )
matplot ( sim.rs , sims.denovo.freq.spec , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 1/2N" , side = 3 )
legend ( "topright" , legend = 1:9 , lty = 1 , col = brewer.pal ( 9 , "Set1" )  , bty = "n" )

matplot ( my.rs , my.specs.f005.s05  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = "Frequency" , xlab = "Genetic Distance" , bty = "n" )
matplot ( sim.rs , sims.f005.freq.spec , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.005" , side = 3 )


matplot ( my.rs , my.specs.f01.s05  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = "Frequency" , xlab = "Genetic Distance" , bty = "n" )
matplot ( sim.rs , sims.f01.freq.spec , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.01" , side = 3 )


matplot ( my.rs , my.specs.f03.s05  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = "Frequency" , xlab = "Genetic Distance" , bty = "n" )
matplot ( sim.rs , sims.f03.freq.spec , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.03" , side = 3 )


matplot ( my.rs , my.specs.f05.s05  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = "Frequency" , xlab = "Genetic Distance" , bty = "n" )
matplot ( sim.rs , sims.f05.freq.spec , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.05" , side = 3 )


matplot ( my.rs , my.specs.f07.s05  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ) , ylab = "Frequency" , xlab = "Genetic Distance" , bty = "n" )
matplot ( sim.rs , sims.f07.freq.spec , type = "p" , pch = 20 , col = brewer.pal ( 9 , "Set1" ) , add = T )
mtext ( "f = 0.07" , side = 3 )
