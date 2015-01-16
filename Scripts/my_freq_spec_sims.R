setwd ( "~/Documents/Academics/StandingSweeps/" )
source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
source('~/Documents/Academics/StandingSweeps/Scripts/freq_spectrum_standing_sweep_coal.R', chdir = TRUE)
source('~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R', chdir = TRUE)
options ( scipen = 400 )



my.N <- 10000
my.rs <- c ( 0 , 0.0000001 , 0.00001 , 0.00001 , 0.001 , 0.004 , 0.008 , 0.012 , 0.016 , 0.02 ) #c ( 0.0001 , 0.001 , 0.01 , 0.05 , 0.1 , 0.5 ) 
#my.fs <- c ( 0.001 , 0.01 , 0.025 , 0.05 , 0.075 , 0.1 )
my.fs <- c ( 0.025 , 0.05 , 0.1 )
s <- 0.05

if ( FALSE ) {
#nosweep.freq.spec.list <- list ()
#for ( i in 1:8 ) {nosweep.freq.spec.list [[ i ]] <- list()}
load("/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/nosweep.freq.spec.list.Rdata")
i <- 41
for ( f in my.fs  [ 5:6 ] ) {

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


load("/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/nosweep.freq.spec.list.Rdata")

f <- 0.025

these.f <- unlist ( nosweep.freq.spec.list[[2]] ) == f

unlist ( nosweep.freq.spec.list [[ 1 ]] ) [ these.f ]

exact.sim.freq.spec <- do.call ( rbind , nosweep.freq.spec.list [[ 6 ]] ) [ these.f , ]
flat.sim.freq.spec <- do.call ( rbind , nosweep.freq.spec.list [[ 7 ]] ) [ these.f , ]
approx.sim.freq.spec <- do.call ( rbind , nosweep.freq.spec.list [[ 8 ]] ) [ these.f , ]

for ( i in 1 : nrow ( exact.sim.freq.spec ) ) {
	blah <- rbind ( exact.sim.freq.spec [ i , ] , approx.sim.freq.spec [ i , ] )
	
	pdf ( paste ( "Figures/FreqSpec/my.spec.f" , f , ".r" , unlist ( nosweep.freq.spec.list [[ 1 ]] ) [ these.f ] [ i ] , ".nosweep.pdf" , sep = "" ) )
	barplot ( blah  , beside = T , col = c ( "red" , "blue" ) , ylim = c ( 0 , 0.4 ) )
	dev.off()
}

} # end of if ( FALSE )



withsweep.freq.spec.list <- list ()
for ( j in 1:8 ) {withsweep.freq.spec.list [[ j ]] <- list()}; rm ( j )
#load("/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/nosweep.freq.spec.list.Rdata")
i <- 1
for ( f in my.fs ) {

	my.runs <-  SweepFromStandingSim ( N = 10000 , s = s , f = f , reps = 1000 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T )

	for ( r in my.rs ) {
		
		
		sim.freq.spec <- run.ms.f ( runs = my.runs [[ 1 ]] , f = f , s = 0.05 , n.sam = 12 , N = my.N , path = "" , ext = "fr.spec", get.site.density = FALSE , recom = 4*my.N*r )
		my.test <- matrix ( f , nrow = 1000 , ncol = 1000 )
		sim.freq.spec.const <- run.ms.f ( runs = my.test , f = f , s = 0.05 , n.sam = 12 , N = my.N , path = "" , ext = "fr.spec", get.site.density = FALSE , recom = 4*my.N*r )
		if ( r != 0 ) {
			approx.freq.spec <- expected.freq.times.standing(nsam=12,N=my.N,r = r , f = f )
		} else {
			approx.freq.spec <- array ( 0 , dim = c ( 12 , 12 , 12 ) )
		}
		
		withsweep.freq.spec.list [[ 1 ]] [[ i ]] <- r
		withsweep.freq.spec.list [[ 2 ]] [[ i ]] <- f
		withsweep.freq.spec.list [[ 3 ]] [[ i ]] <- sim.freq.spec
		withsweep.freq.spec.list [[ 4 ]] [[ i ]] <- sim.freq.spec.const
		withsweep.freq.spec.list [[ 5 ]] [[ i ]] <- approx.freq.spec
		withsweep.freq.spec.list [[ 6 ]] [[ i ]] <- rowSums ( withsweep.freq.spec.list [[ 3 ]] [[ i ]] ) / sum ( withsweep.freq.spec.list [[ 3 ]] [[ i ]] )
		withsweep.freq.spec.list [[ 7 ]] [[ i ]] <- rowSums ( withsweep.freq.spec.list [[ 4 ]] [[ i ]] ) / sum ( withsweep.freq.spec.list [[ 4 ]] [[ i ]] )
		withsweep.freq.spec.list [[ 8 ]] [[ i ]] <- apply ( withsweep.freq.spec.list [[ 5 ]] [[ i ]]  , 3 , sum )
		save ( withsweep.freq.spec.list , file = "Sims/withsweep.freq.spec.list.Rdata" )
		
		message ( r )
		message ( f )
		message ( i )
		
		i <- i + 1
	}
}




if ( FALSE ) {
	
load("/Users/JeremyBerg/Documents/Academics/StandingSweeps/Sims/withsweep.freq.spec.list.Rdata")

f <- 0.025

these.f <- unlist ( nosweep.freq.spec.list[[2]] ) == f

unlist ( nosweep.freq.spec.list [[ 1 ]] ) [ these.f ]

exact.sim.freq.spec <- do.call ( rbind , nosweep.freq.spec.list [[ 6 ]] ) [ these.f , ]
flat.sim.freq.spec <- do.call ( rbind , nosweep.freq.spec.list [[ 7 ]] ) [ these.f , ]
approx.sim.freq.spec <- do.call ( rbind , nosweep.freq.spec.list [[ 8 ]] ) [ these.f , ]

for ( i in 1 : nrow ( exact.sim.freq.spec ) ) {
	blah <- rbind ( exact.sim.freq.spec [ i , -ncol ( exact.sim.freq.spec ) ] , approx.sim.freq.spec [ i , -ncol ( exact.sim.freq.spec ) ] )
	
	pdf ( paste ( "Figures/FreqSpec/my.spec.f" , f , ".r" , unlist ( nosweep.freq.spec.list [[ 1 ]] ) [ these.f ] [ i ] , "s" , 0.05 , ".pdf" , sep = "" ) )
	barplot ( blah  , beside = T , col = c ( "red" , "blue" ) , ylim = c ( 0 , 0.4 ) )
	#mtext ( "Counts" , side = 1 , line = 2 , cex = 1.5 )
	#mtext ( "Frequency" , side = 2 , line = 2 , cex = 1.5 )
	axis ( 1 , 1 + seq ( 1 , 31 , 3 ) , labels = 1:11 )
	dev.off()
}
	
	
	
	
	
}












