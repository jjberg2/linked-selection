setwd ( "~/Documents/Academics/StandingSweeps/" )
source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
source('~/Documents/Academics/StandingSweeps/Scripts/freq_spectrum_standing_sweep_coal.R', chdir = TRUE)
source('~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R', chdir = TRUE)
options ( scipen = 400)



my.N <- 10000
my.rs <- c ( 0.0001 , 0.001 , 0.01 , 0.05 , 0.1 , 0.5 ) 
my.fs <- c ( 1/(2*my.N) , 0.001 , 0.01 , 0.025 , 0.05 , 0.075 , 0.1 )
my.s <- my.fs






my.runs <-  SweepFromStandingSim ( N = 10000 , s = 0.05 , f = 0.05 , reps = 1000 , no.sweep = TRUE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T )
i <- 1
nosweep.freq.spec.list <- list ()
for ( i in 1:5 ) nosweep.freq.spec.list [[ i ]] <- list()
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
		
		save ( nosweep.freq.spec.list , file = "Sims/nosweep.freq.spec.list.Robj" )
		
		message ( r )
		message ( f )
		message ( i )
		
		i <- i + 1
	}
}