setwd ( "~/Documents/Academics/StandingSweeps/" )
source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
source ( '~/Documents/Academics/StandingSweeps/Scripts/freq_spectrum_standing_sweep_coal.R' )
source('~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R', chdir = TRUE)
options ( scipen = 400 )


my.N <- 10000
my.rs <- c ( seq ( 0 , 0.0009 , length.out = 10 ) , seq ( 0.001 , 0.002 , length.out = 5 ) , c ( 0.0025 , 0.003 , 0.004))
s <- 0.05
my.fs <- c ( 0.01 , 0.02 , 0.03 , 0.05 , 0.07 )
sim.freq.spec.list <- list ()
i <- 1
for ( f in my.fs ) {

	my.runs <-  SweepFromStandingSim ( N = 10000 , s = s , f = f , reps = 1000 , no.sweep = TRUE , cond.on.loss = T , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T )
	sim.freq.spec.list [[ i ]] <- list ()
	j <- 1
	for ( r in my.rs ) {

		
		sim.freq.spec <- run.ms.f ( runs = my.runs [[ 1 ]] , f = f , s = s , n.sam = 10 , N = 10000 , num.sims = 30 , path = "" , get.site.density = FALSE , recom = 4*10000*r )
		sim.freq.spec.list [[ i ]] [[ j ]] <- list ()
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 1 ]] <- r
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 2 ]] <- f
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 3 ]] <- s
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 4 ]] <- sim.freq.spec
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 5 ]] <- rowSums ( sim.freq.spec.list [[ i ]] [[ j ]] [[ 4 ]] ) / sum ( sim.freq.spec.list [[ i ]] [[ j ]] [[ 4 ]])
		save ( sim.freq.spec.list , file = "Sims/sim.freq.spec.list.stoch.freq.nosweep.condloss.Rdata" )
		
		message ( r )
		message ( f )
		j <- j + 1
	}
	i <- i + 1
}





my.N <- 10000
sim.freq.spec.list <- list ()
i <- 1
for ( f in my.fs ) {

	my.runs <-  SweepFromStandingSim ( N = 10000 , s = s , f = f , reps = 1000 , no.sweep = FALSE , cond.on.loss = T , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T )
	sim.freq.spec.list [[ i ]] <- list ()
	j <- 1
	for ( r in my.rs ) {

		
		sim.freq.spec <- run.ms.f ( runs = my.runs [[ 1 ]] , f = f , s = s , n.sam = 10 , N = 10000 , num.sims = 10 , path = "" , get.site.density = FALSE , recom = 4*10000*r )
		sim.freq.spec.list [[ i ]] [[ j ]] <- list ()
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 1 ]] <- r
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 2 ]] <- f
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 3 ]] <- s
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 4 ]] <- sim.freq.spec
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 5 ]] <- rowSums ( sim.freq.spec.list [[ i ]] [[ j ]] [[ 4 ]] ) / sum ( sim.freq.spec.list [[ i ]] [[ j ]] [[ 4 ]])
		save ( sim.freq.spec.list , file = "Sims/sim.freq.spec.list.stoch.freq.with.sweep.condloss.Rdata" )
		
		message ( r )
		message ( f )
		j <- j + 1
	}
	i <- i + 1
}






if ( FALSE ) {


###### checking freq specs with no sweep

load ( "Sims/sim.freq.spec.list.stoch.freq.no.sweep.nocondloss.Rdata" )
load ( "Sims/sim.freq.spec.list.stoch.freq.no.sweep.condloss.Rdata" )
load ( "Sims/sim.freq.spec.list.stoch.freq.with.sweep.nocondloss.Rdata" )







## full stochastic sims
my.N <- 10000
my.rs <- c ( seq ( 0 , 0.001 , length.out = 11 ) , 0.0015 , 0.002 , 0.003 , 0.004 , 0.005 )
s <- 0.05
my.fs <- c ( 1/20000 , 0.005 , 0.01 , 0.03 , 0.05 , 0.07 )
sim.freq.spec.list <- list ()
i <- 1
j <- 1
for ( f in my.fs ) {

	my.runs <-  SweepFromStandingSim ( N = 10000 , s = s , f = f , reps = 1000 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T )
	sim.freq.spec.list [[ i ]] <- list ()
	for ( r in my.rs ) {

		
		sim.freq.spec <- run.ms.f ( runs = my.runs [[ 1 ]] , f = f , s = s , n.sam = 10 , N = 10000 , num.sims = 10 , path = "" , get.site.density = FALSE , recom = 4*10000*r )
		sim.freq.spec.list [[ i ]] [[ j ]] <- list ()
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 1 ]] <- r
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 2 ]] <- f
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 3 ]] <- s
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 4 ]] <- sim.freq.spec
		sim.freq.spec.list [[ i ]] [[ j ]] [[ 5 ]] <- rowSums ( sim.freq.spec.list [[ i ]] [[ j ]] [[ 4 ]] ) / sum ( sim.freq.spec.list [[ i ]] [[ j ]] [[ 4 ]])
		save ( sim.freq.spec.list , file = "Sims/sim.freq.spec.list.s05.Rdata" )
		
		message ( r )
		message ( f )
		j <- j + 1
	}
	i <- i + 1
}

load("/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/sim.freq.spec.list.Rdata")
i <- 1 
for ( f in my.fs ) {
	
	for ( r in my.rs ) {
		
		sim.freq.spec.list [[ 8 ]] [[ i ]] <- expected.freq.times.standing.w.sweep (nsam = 12 , N = 10000 , r = r , f = f , s = 0.05 )
		i <- i + 1
	}
	
}
save ( sim.freq.spec.list , file = "/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/sim.freq.spec.list.Rdata" )	

f <- 0.025

these.f <- unlist ( sim.freq.spec.list[[2]] ) == f

unlist ( nosweep.sim.freq.spec.list [[ 1 ]] ) [ these.f ]

exact.sim.freq.spec <- do.call ( rbind , sim.freq.spec.list [[ 6 ]] ) [ these.f , ]
flat.sim.freq.spec <- do.call ( rbind , sim.freq.spec.list [[ 7 ]] ) [ these.f , ]
approx.sim.freq.spec <- do.call ( rbind , sim.freq.spec.list [[ 8 ]] ) [ these.f , ]

for ( i in 1 : nrow ( exact.sim.freq.spec ) ) {
	blah <- rbind ( exact.sim.freq.spec [ i , -ncol ( exact.sim.freq.spec ) ] , approx.sim.freq.spec [ i , -ncol ( exact.sim.freq.spec ) ] )
	
	pdf ( paste ( "Figures/FreqSpec/my.spec.f" , f , ".r" , unlist ( nosweep.sim.freq.spec.list [[ 1 ]] ) [ these.f ] [ i ] , "s" , 0.05 , ".pdf" , sep = "" ) )
	barplot ( blah  , beside = T , col = c ( "red" , "blue" ) , ylim = c ( 0 , 0.4 ) )
	#mtext ( "Counts" , side = 1 , line = 2 , cex = 1.5 )
	#mtext ( "Frequency" , side = 2 , line = 2 , cex = 1.5 )
	axis ( 1 , 1 + seq ( 1 , 31 , 3 ) , labels = 1:11 )
	dev.off()
}
	
	
	
	
	
}












