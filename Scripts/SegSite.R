real.fs <- c ( 0.001 , 0.01 , 0.05 , 0.1 )
my.runs <- lapply ( real.fs , function ( x ) SweepFromStandingSim ( N = 10000 , s = 0.05 , f = x , reps = 1000 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T ) )
many.sim.freqs <- list ( my.runs , real.fs )
save ( many.sim.freqs  ,  file = "~/Documents/Academics/StandingSweeps/Sims/freqs.traj.s05.Rdata"  )
for ( i in 1 : length ( real.fs ) ) {
	run.ms.f ( runs = my.runs [[ i ]] [[ 1 ]] , n.sam = 12 , f = real.fs [ i ] , s = 0.05 , N = 10000 , path = "~/Documents/Academics/StandingSweeps/" )
}