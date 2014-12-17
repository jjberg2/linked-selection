setwd ( "~/Documents/Academics/StandingSweeps/" )
source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
source('~/Documents/Academics/StandingSweeps/Scripts/freq_spectrum_standing_sweep_coal.R', chdir = TRUE)
source('~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R', chdir = TRUE)
options ( scipen = 400)



my.N <- 10000
my.rs <- c ( 0 , 0.0000001 , 0.00001 , 0.00001 , 0.001 , 0.004 , 0.008 , 0.012 , 0.016 , 0.02 ) #c ( 0.0001 , 0.001 , 0.01 , 0.05 , 0.1 , 0.5 ) 
my.fs <- c ( 0.001 , 0.01 , 0.025 , 0.05 , 0.075 , 0.1 )
my.s <- my.fs


#nosweep.freq.spec.list <- list ()
#for ( i in 1:8 ) {nosweep.freq.spec.list [[ i ]] <- list()}
load("/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/nosweep.freq.spec.list.Rdata")
i <- 31
for ( f in my.fs  [ 4:6 ] ) {

	my.runs <-  SweepFromStandingSim ( N = 10000 , s = 0.05 , f = f , reps = 1000 , no.sweep = TRUE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T )

	for ( r in my.rs ) {

		
		sim.freq.spec <- run.ms.f ( runs = my.runs [[ 1 ]] , f = f , s = 0.05 , n.sam = 12 , N = my.N , path = "" , ext = "fr.spec", get.site.density = FALSE , recom = 4*my.N*r )
		my.test <- matrix ( f , nrow = 1000 , ncol = 1000 )
		sim.freq.spec.const <- run.ms.f ( runs = my.test , f = f , s = 0.05 , n.sam = 12 , N = my.N , path = "" , ext = "fr.spec", get.site.density = FALSE , recom = 4*my.N*r )
		if ( r != 0 ) {
			approx.freq.spec <- expected.freq.times.standing(nsam=12,N=my.N,r = r , f = f )
		} else {
			approx.freq.spec <- array ( 0 , dim = c ( 12 , 12 , 12 ) )
		}
		
		nosweep.freq.spec.list [[ 1 ]] [[ i ]] <- r
		nosweep.freq.spec.list [[ 2 ]] [[ i ]] <- f
		nosweep.freq.spec.list [[ 3 ]] [[ i ]] <- sim.freq.spec
		nosweep.freq.spec.list [[ 4 ]] [[ i ]] <- sim.freq.spec.const
		nosweep.freq.spec.list [[ 5 ]] [[ i ]] <- approx.freq.spec
		nosweep.freq.spec.list [[ 6 ]] [[ i ]] <- rowSums ( nosweep.freq.spec.list [[ 3 ]] [[ i ]] ) / sum ( nosweep.freq.spec.list [[ 3 ]] [[ i ]] )
		nosweep.freq.spec.list [[ 7 ]] [[ i ]] <- rowSums ( nosweep.freq.spec.list [[ 4 ]] [[ i ]] ) / sum ( nosweep.freq.spec.list [[ 4 ]] [[ i ]] )
		nosweep.freq.spec.list [[ 8 ]] [[ i ]] <- apply ( nosweep.freq.spec.list [[ 5 ]] [[ i ]]  , 3 , sum )
		save ( nosweep.freq.spec.list , file = "Sims/nosweep.freq.spec.list.Rdata" )
		
		message ( r )
		message ( f )
		message ( i )
		
		i <- i + 1
	}
}

if ( FALSE ) {
load("/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/nosweep.freq.spec.list.Rdata")


for ( i in 1 : length ( my.fs [ 1:3] ) ) {
	
	this.f <- my.fs [ i ]
	
	these.lists <- which ( unlist ( nosweep.freq.spec.list [[ 2 ]] ) == this.f )
	these.true.specs <- lapply ( these.lists , function ( x ) nosweep.freq.spec.list [[ 6 ]] [[ x ]] )
	these.const.specs <- lapply ( these.lists , function ( x ) nosweep.freq.spec.list [[ 7 ]] [[ x ]] )
	these.approx.specs <- lapply ( these.lists , function ( x ) nosweep.freq.spec.list [[ 8 ]] [[ x ]] )
	
	neut.spec <- (1 / ( 1 : 11 ) ) / sum ( 1 / ( 1 : 11 ) )
	neut.scale.true.spec <- lapply ( these.true.specs , function ( x ) x [ 1 : length ( neut.spec ) ] / neut.spec )
	neut.scale.const.spec <- lapply ( these.const.specs , function ( x ) x [ 1 : length ( neut.spec ) ] / neut.spec )
	neut.scale.approx.spec <- lapply ( these.approx.specs , function ( x ) x [ 1 : length ( neut.spec ) ] / neut.spec )
	
	neut.scale.true.spec <- do.call ( rbind , neut.scale.true.spec )
	neut.scale.const.spec <- do.call ( rbind , neut.scale.const.spec )
	neut.scale.approx.spec <- do.call ( rbind , neut.scale.approx.spec )
	
	pdf ( paste ( "Figures/freq.spec.nosweep.f" , strsplit ( as.character ( this.f ) , "\\." )[[1]][2] , ".pdf" , sep = "" ) , width = 6 , height = 22 )
	par ( mfrow = c ( 11, 1 ) )
	for ( i in 1 : ncol ( neut.scale.true.spec ) ) {
	
		plot ( my.rs , neut.scale.true.spec [ , i ] , type = "l" , bty = "n" , col = "black" , lwd = 2 , lty = 1 , ylim = c ( 0 , 6 ) , xlim = c ( 0 , 0.02 )  , xlab = "r" , ylab = "" , yaxt = "n" )
		lines ( my.rs , neut.scale.const.spec [ , i ] , col = "black" , lwd = 2 , lty = 2 )
		lines ( my.rs [ -1 ] , neut.scale.approx.spec [ -1 , i ] , col = "red" , lwd = 2 , lty = 1 )
		abline ( h = 1 , lwd = 1 , lty = 3 , col = "blue" )
		# if ( r != 0 ) {
		# lines ( my.rs , neut.scale.approx.spec [ , i ] , col = "red" , lwd = 2 , lty = 1 )	
		# } else {

		# }
		axis ( 2 , at = seq ( 0 , 6 , 1 ) )
	
	}
	dev.off()
	
}








my.N <- 10000
my.rs <- c ( 0.0001 , 0.001 , 0.01 , 0.05 , 0.1 , 0.5 ) 
my.fs <- c ( 1/(2*my.N) , 0.001 , 0.01 , 0.025 , 0.05 , 0.075 , 0.1 )
my.s <- my.fs

my.runs <-  SweepFromStandingSim ( N = 10000 , s = 0.05 , f = 0.05 , reps = 1000 , no.sweep = TRUE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T )

nosweep.freq.spec.list <- list ()
for ( i in 1:5 ) nosweep.freq.spec.list [[ i ]] <- list()

i <- 1
for ( r in my.rs ) {

	for ( f in my.fs ) {
		
		sim.freq.spec <- run.ms.f ( runs = my.runs [[ 1 ]] , f = f , s = 0.05 , n.sam = 12 , N = my.N , path = "" , ext = "fr.spec", get.site.density = FALSE , recom = 4*my.N*r )
		my.test <- matrix ( f , nrow = 1000 , ncol = 1000 )
		sim.freq.spec.const <- run.ms.f ( runs = my.test , f = f , s = 0.05 , n.sam = 12 , N = my.N , path = "" , ext = "fr.spec", get.site.density = FALSE , recom = 4*my.N*r )
		approx.freq.spec <- expected.freq.times.standing(nsam=12,N=my.N,r = r , f = f )
		
		nosweep.freq.spec.list [[ 1 ]] [[ i ]] <- r
		nosweep.freq.spec.list [[ 2 ]] [[ i ]] <- f
		nosweep.freq.spec.list [[ 3 ]] [[ i ]] <- sim.freq.spec
		nosweep.freq.spec.list [[ 4 ]] [[ i ]] <- sim.freq.spec.const
		nosweep.freq.spec.list [[ 5 ]] [[ i ]] <- approx.freq.spec
		
		save ( nosweep.freq.spec.list , file = "Sims/nosweep.freq.spec.list.Rdata" )
		
		message ( r )
		message ( f )
		message ( i )
		
		i <- i + 1
	}
}




} # end of if ( FALSE )









