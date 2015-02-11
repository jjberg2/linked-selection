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

stirlings <- StirlingNumbers(50)
neutral.f.spec <- 1 / ( 1 : 49 ) / ( sum ( 1 / ( 1:49 ) ) )
f.specs.range <- list ()
my.rs <- c ( 0.000001 , 0.00001 , 0.00005 , seq ( 0.0001 , 0.001 , by = 0.0001 ) , 0.005 , 0.01 )
for ( i in seq_along ( my.rs ) ) {

	f.specs.range [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 50 , N = 10000 , r = my.rs [ i ] , f = 0.05 , s = 0.01 , my.StirlingNumbers = stirlings )
	save ( f.specs.range , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.range.Robj")

}
my.f.specs <- do.call ( rbind , f.specs.range )
matplot ( log ( t ( t ( my.f.specs ) / neutral.f.spec ) ) , type = "l" )
my.f.specs / neutral.f.spec




stirlings <- StirlingNumbers(20)
neutral.f.spec <- 1 / ( 1 : 19 ) / ( sum ( 1 / ( 1:19 ) ) )
f.specs.range <- list ()
my.rs <- c ( 0.000001 , 0.00001 , 0.00005 , seq ( 0.0001 , 0.001 , by = 0.0001 ) , 0.005 , 0.01 )
for ( i in seq_along ( my.rs ) ) {

	f.specs.range [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 20 , N = 10000 , r = my.rs [ i ] , f = 0.05 , s = 0.01 , my.StirlingNumbers = stirlings )
	save ( f.specs.range , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef05s01N10000n20.Robj")

}

my.f.specs <- do.call ( rbind , f.specs.range )
matplot ( log ( t ( t ( my.f.specs ) / neutral.f.spec ) ) , type = "l" )





stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f.specs.range <- list ()
my.rs <- c ( 0.000001 , 0.00001 , 0.00005 , seq ( 0.0001 , 0.001 , by = 0.0001 ) , 0.005 , 0.01 )
for ( i in seq_along ( my.rs ) ) {

	f.specs.range [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.05 , s = 0.01 , my.StirlingNumbers = stirlings )
	save ( f.specs.range , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef05s01N10000n10.Robj")

}
load ( file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef05s01N10000n10.Robj")
my.f.specs <- do.call ( rbind , f.specs.range )

matplot ( my.f.specs  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ))
matplot ( t ( t ( my.f.specs ) / neutral.f.spec ) , type = "l" )
matplot ( log ( t ( t ( my.f.specs ) / neutral.f.spec ) ) , type = "l" )






stirlings <- StirlingNumbers(10)
neutral.f.spec <- 1 / ( 1 : 9 ) / ( sum ( 1 / ( 1:9 ) ) )
f.specs.range <- list ()
my.rs <- c ( 0.000001 , 0.00001 , 0.00005 , seq ( 0.0001 , 0.001 , by = 0.0001 ) , 0.005 , 0.01 )
for ( i in seq_along ( my.rs ) ) {

	f.specs.range [[ i ]] <- expected.freq.times.standing.w.sweep ( nsam = 10 , N = 10000 , r = my.rs [ i ] , f = 0.0001 , s = 0.01 , my.StirlingNumbers = stirlings )
	save ( f.specs.range , file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef0001s01N10000n10.Robj")

}
load ( file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef0001s01N10000n10.Robj")
my.f.specs <- do.call ( rbind , f.specs.range )

matplot ( my.f.specs  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ))
matplot ( t ( t ( my.f.specs ) / neutral.f.spec ) , type = "l" )
matplot ( log ( t ( t ( my.f.specs ) / neutral.f.spec ) ) , type = "l" )



load ( file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef05s05N10000n10.Robj")
my.f.specs.stand <- do.call ( rbind , f.specs.range )
load ( file = "Paper/Paper_Figures/Data_and_Robjs/f.specs.rangef0001s05N10000n10.Robj")
my.f.specs.denovo <- do.call ( rbind , f.specs.range )

par ( mfrow = c ( 1 , 2 ) )

matplot ( my.f.specs.stand  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ))
matplot ( my.f.specs.denovo  , type = "l"  , lty = 1 , col = brewer.pal ( 9 , "Set1" ))
