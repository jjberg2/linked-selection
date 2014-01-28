SoftSweepSim <- function ( N_A , N_a , s1 , s2 = NULL , h , switch = 0 , mu , gens , mut.index = 1001 , stop.at.fix = TRUE , reuse.mutant.classes = FALSE , suppress.output = FALSE , suppress.plot = FALSE , sort.mutation.time = FALSE ) {
	options ( error = recover )
	
	##########################
	####### Initialize #######
	##########################	
	
	num.A.class <- length ( N_A )
	num.a.class <- length ( N_a )
	pA.vect <- matrix ( c ( N_A / sum ( N_A , N_a ), rep ( 0, length ( N_A ) ) ) , nrow = num.A.class , ncol = 1 )	
	pa.vect <- matrix ( c ( N_a / sum ( N_A , N_a ), rep ( 0, length ( N_a ) ) ) , nrow = num.a.class , ncol = 1 )	
	A.classes <- 1 : num.A.class
	a.classes <- mut.index : ( mut.index + num.a.class -1 )
	if ( switch > 0 ){
		s <- s1
	} else {
		s <- s2
	}
	allele.classes.all <- c ( A.classes , mut.index : ( mut.index + ( num.a.class - 1 ) ) )
	pA.presel <- matrix ( N_A / sum ( N_A , N_a ) )
	pa.presel <- matrix ( N_a / sum ( N_A , N_a ) )
	
	
	##########################
	####### Simulation #######
	##########################
	i <- 1
	while ( ( sum ( pa.vect [ , i ] ) < 1 | stop.at.fix == FALSE ) & i < gens ) {
		if ( mu == 0 & sum ( N_a ) == 0 ) {
			cat ( "Allele lost from population. \n")
			return()
		}
		###############
		## selection ##
		###############
		pA.vect <- cbind ( pA.vect , 0 )
		pa.vect <- cbind ( pa.vect , 0 )
		pA.temp <- rowSums ( pA.presel %*% t ( pA.presel ) ) + rowSums ( pA.presel %*% t ( pa.presel ) ) * ( 1 + s )
		pa.temp <- colSums ( pa.presel %*% t ( pa.presel ) ) * ( 1 + 2*s ) + colSums ( pA.presel %*% t ( pa.presel ) ) * ( 1 + s )
		mean.fitness <- sum( pA.presel %*% t ( pA.presel ) ) + 2*sum ( pA.presel %*% t ( pa.presel ) )*( 1 + s ) + sum( pa.presel %*% t ( pa.presel ) ) * ( 1 + 2*s )
		pA.postsel <- pA.temp / mean.fitness
		pa.postsel <- pa.temp / mean.fitness
		pA.vect [ , i + 1 ] <- pA.postsel
		pa.vect [ , i + 1 ] <- c ( pa.postsel , rep ( 0 , nrow ( pa.vect ) - length ( pa.postsel ) ) )
		
		###########
		## drift ##
		###########
		next.gen <- sample ( c ( A.classes , a.classes ) , size = sum ( N_A , N_a ) , prob = c ( pA.vect [ , i + 1 ] , pa.vect [ , i + 1 ] ) , replace = TRUE )
		next.gen <- sort ( next.gen )
		allele.classes.present <- rle ( next.gen )$values
		allele.classes.absent <- allele.classes.all [ allele.classes.all %in% allele.classes.present == FALSE ]
		allele.abundance <- c ( rle ( next.gen )$lengths , rep ( 0 , length ( allele.classes.absent ) ) )
		allele.classes.all <- c ( allele.classes.present , allele.classes.absent )
		allele.abundance <- allele.abundance [ order ( allele.classes.all ) ]
		allele.classes.all <- allele.classes.all [ order ( allele.classes.all ) ]
		A.classes <- allele.classes.all [ allele.classes.all < mut.index ] 
		a.classes <- allele.classes.all [ allele.classes.all >= mut.index ] 
		N_A <- allele.abundance [ allele.classes.all < mut.index ]
		N_a <- allele.abundance [ allele.classes.all >= mut.index ]
		pA.premut <- matrix ( N_A / sum ( N_A , N_a ) )
		pa.premut <- matrix ( N_a / sum ( N_A , N_a ) )		
		
		##############
		## mutation ##
		##############
		if ( mu != 0 ) {
			## roll N_total uniform random numbers; the number that are less than 	
			mut.probs <- runif ( n = sum ( N_A , N_a ) )
			new.mutations <- sum ( mut.probs < mu )
			if ( new.mutations != 0 ) {
				# let the mutations pick parents in proportion to their frequency
				mut.from <- sort ( sample ( c ( rep ( A.classes , times = N_A ) , rep ( a.classes , times = N_a ) ) , size = new.mutations , replace = FALSE ) )
				mutating.alleles <- rle ( mut.from )$values
				non.mutating.alleles <- allele.classes.all [ allele.classes.all %in% mutating.alleles == FALSE ]
				mutation.counts.all.alleles <- c ( rle ( mut.from )$lengths , rep ( 0 , length ( non.mutating.alleles ) ) )
				mutation.id.all.alleles <- c ( mutating.alleles , non.mutating.alleles )
				mutation.counts.all.alleles <- mutation.counts.all.alleles [ order ( mutation.id.all.alleles ) ]
				mutation.id.all.alleles <- mutation.id.all.alleles [ order ( mutation.id.all.alleles ) ]
				mut.from.A <- mutation.counts.all.alleles [ mutation.id.all.alleles < mut.index ]
				mut.from.a <- mutation.counts.all.alleles [ mutation.id.all.alleles >= mut.index ]
				if ( reuse.mutant.classes == TRUE ) {
					idx <- c( a.classes [ N_a == 0 ] , ( max ( mutation.id.all.alleles ) + 1 ) : ( max ( mutation.id.all.alleles ) + new.mutations ) )
					new.mut.idx <- idx [ 1 : new.mutations ]
					new.IDs <- new.mut.idx [ new.mut.idx %in% a.classes == FALSE ]
				} else if ( reuse.mutant.classes == FALSE ) {
					new.mut.idx <- new.IDs <- ( max ( mutation.id.all.alleles ) + 1 ) : ( max ( mutation.id.all.alleles ) + new.mutations )
				}
				a.classes <- c ( a.classes, new.IDs )
				N_a <- c ( N_a , rep ( 0 , length ( new.IDs ) ) )
				pa.vect <- rbind ( pa.vect , matrix ( 0 , nrow = length ( new.IDs ) , ncol = ncol ( pa.vect ) ) )
				N_a [ a.classes %in% new.mut.idx ] <- 1
				N_A <- N_A - mut.from.A
				N_a [ a.classes %in% new.IDs == FALSE ] <- N_a [ a.classes %in% new.IDs == FALSE ] - mut.from.a
				allele.classes.all <- c ( A.classes , a.classes )
			}
		}
		
		#################
		## bookkeeping ##
		#################		
		if ( i == switch & is.null ( s2 ) == FALSE ) {
			s <-  s2 
		}
		i <- i + 1
		pA.presel <- matrix ( N_A / sum ( N_A , N_a ) )
		pa.presel <- matrix ( N_a / sum ( N_A , N_a ) )
		if ( i %% 20 == 0 & suppress.output == FALSE ){
			cat (  i , "p =" , sum ( pa.vect [ , i ] ) , "\n" )
		}
	}
	if ( switch > 100 ) {
		pre.switch.freqs <- pa.vect [ , 101:switch ]
	} else {
		pre.switch.freqs <- NULL
		mean.coal.prob.sim.mut.sel.drift <- NULL
		coal.prob.sim.mut.sel.drift <- NULL
	}
	
	
	##########################
	##### End Point Data #####
	##########################
	
	pA.order <- order ( pA.vect [ , ncol ( pA.vect ) ] , decreasing = FALSE ) 
	pa.order <- order ( pa.vect [ , ncol ( pa.vect ) ] , decreasing = FALSE )
	pa.lost <-  which ( pa.vect [ , ncol ( pa.vect ) ] == 0 )
	pa.fixed <- which ( pa.vect [ , ncol ( pa.vect ) ] != 0 )
	all.freqs <- rbind ( pa.vect [ pa.order , ] , pA.vect [ pA.order , ] )
	all.freqs.unsorted <- rbind ( pa.vect , pA.vect )



	##########################
	######## Plotting ########
	##########################


	if ( suppress.plot == FALSE ) {
		if ( reuse.mutant.classes == TRUE ) {	
			
			#recover()
			cum.freqs <- rbind ( 0 , apply ( all.freqs , 2 , cumsum ))
			green.mat <- matrix ( c ( rep ( 0 , length ( pa.order ) + 1 ) , seq ( from = .3 , to = 1 , length.out = length ( pa.order ) + 1 ) , rep ( 0 , length ( pa.order ) +1 ) ) , nrow = length ( pa.order ) + 1 )
			blue.mat <- matrix ( c ( rep ( 0 , 2 * length ( pA.order ) ) , seq ( from = .7 , to = 1 , length.out = length ( pA.order ) ) ) , nrow = length ( pA.order ) )
			col.vect <- rgb ( rbind ( green.mat , blue.mat) , alpha = .65 , maxColorValue = 1 )
		} else if ( sort.mutation.time == FALSE ) {
			
			cum.freqs <- rbind ( 0 , apply ( all.freqs , 2 , cumsum ))
			
			if ( length ( pa.lost ) != 0 ) {
				red.mat <- matrix ( c ( seq ( from = .3 , to = 1 , length.out = length ( pa.lost ) ) [ 1 : length ( pa.lost ) ] , rep ( 0 , 2 * length ( pa.lost ) ) ) , nrow = length ( pa.lost ) )
			} else {
				red.mat <- matrix ( nrow = 0 , ncol = 3 )
			}
			green.mat <- matrix ( c ( rep ( 0 , length ( pa.fixed ) ) , rev ( seq ( from = 1 , to = .4 , length.out = length ( pa.fixed ) ) ) [ 1 : length ( pa.fixed ) ] , rep ( 0 , length ( pa.fixed ) ) ) , nrow = length ( pa.fixed ) )
			blue.mat <- matrix ( c ( rep ( 0 , 2 * length ( pA.order ) ) , seq ( from = .7 , to = 1 , length.out = length ( pA.order ) ) ) , nrow = length ( pA.order ) )
			col.vect <- rgb ( rbind ( red.mat , green.mat , blue.mat) , alpha = .65 , maxColorValue = 1 )
		} else if ( sort.mutation.time == TRUE ) {
			cum.freqs <-  apply ( all.freqs.unsorted , 2 , cumsum )
			if ( length ( pa.lost ) != 0 ) {
				red.mat <- matrix ( c ( seq ( from = .3 , to = 1 , length.out = length ( pa.lost ) ) [ 1 : length ( pa.lost ) ] , rep ( 0 , 2 * length ( pa.lost ) ) ) , nrow = length ( pa.lost ) )
			} else {
				red.mat <- matrix ( nrow = 0 , ncol = 3 )
			}
			green.mat <- matrix ( c ( rep ( 0 , length ( pa.fixed ) ) , rev ( seq ( from = 1 , to = .4 , length.out = length ( pa.fixed ) ) ) [ 1 : length ( pa.fixed ) ] , rep ( 0 , length ( pa.fixed ) ) ) , nrow = length ( pa.fixed ) )
			blue.mat <- matrix ( c ( rep ( 0 , 2 * length ( pA.order ) ) , seq ( from = .7 , to = 1 , length.out = length ( pA.order ) ) ) , nrow = length ( pA.order ) )
							
			col.vect <- character ( nrow ( cum.freqs ) )
			col.vect [ pa.lost ] <- rgb ( red.mat , alpha = .65 , maxColorValue = 1 )
			col.vect [ pa.fixed ] <- rgb ( green.mat , alpha = .65 , maxColorValue = 1 )
			col.vect [ ( nrow ( pa.vect ) + 1 ) : nrow ( cum.freqs ) ] <- rgb ( blue.mat , alpha = .65 , maxColorValue = 1 )
		}
	#	recover()
		no.lines <- sum ( cum.freqs [ , ncol ( cum.freqs ) ] < 0.01 )
		line.vector <- c ( rep ( 0 , no.lines ) , rep ( .35 , nrow ( cum.freqs ) - no.lines ) )
		matplot ( t ( cum.freqs ) , type = "l" , lty = "solid" , lwd = line.vector , col = "black" , xlab = "Generations" , ylab = "Cumulative Frequency" , bty = "n" )
		for ( i in 1 : ( nrow ( cum.freqs ) - 1 ) ) {
			#i = i + 1
			X.ax <- which ( cum.freqs [ i , ] != cum.freqs [ i + 1 , ] )
			Y.ax1 <- cum.freqs [ i , X.ax ]
			Y.ax2 <- cum.freqs [ i + 1 , X.ax ]
			polygon ( x = c ( X.ax , rev ( X.ax ) ) , y = c ( Y.ax1 , rev ( Y.ax2 ) ) , lty = 0 , col = col.vect [ i ] )					
		}
		if ( switch > 0 ) {
			segments (x0 = switch , y0 = 0 , y1 = 1.0 , lty = 2 )
		}
	}	
	
	#mut.sel.max <- max(cum.freqs [ length (pa.lost ), ])
	#matplot ( t ( cum.freqs [1: length ( c ( pa.lost , pa.fixed ) ) , ] ) , type = "l" , lty = 1 , lwd = 1 , col = "black" , ylim = c ( 0 , mut.sel.max * 1.1 ) )
	
	
	
	if ( switch > 100 ) {
		coal.prob.sim.mut.sel.drift <- numeric ( ncol ( pre.switch.freqs ) )
		for ( i in 1 : ncol ( pre.switch.freqs ) ) {
				coal.prob.sim.mut.sel.drift [ i ] <- sum ( ( pre.switch.freqs [ , i ] / sum ( pre.switch.freqs [ , i ] ) )^2 )
		}
	mean.coal.prob.sim.mut.sel.drift <- mean ( coal.prob.sim.mut.sel.drift )
	}
	coal.prob.sim.fix <- sum( (pa.vect [ , ncol(pa.vect) ])^2 )
	coal.prob.theory.fix <- 1 / ( 1 + 2 * sum ( N_A , N_a ) * mu )
	if ( s1 < 0 ) {
		coal.prob.theory.mut.sel.determin <- 1 / ( 1 + 2 * sum ( N_A , N_a ) * mu / (1 + abs ( s1 ) ) )	
	} else {
		coal.prob.theory.mut.sel.determin = NULL
	}
	if ( suppress.output == FALSE ) {
		cat ( "Coalescence probability from simulation for lineages sampled at fixation =" , coal.prob.sim.fix , "\n" )
		cat ( "Coalescence probability from theory for lineages sampled at fixation =" , coal.prob.theory.fix , "\n" )
		if ( switch > 100 ) {
			cat ( "Coalescence probability from simulation for lineages sampled at mut-sel-drift balance =" , mean.coal.prob.sim.mut.sel.drift , "\n" )
		}
		cat ( "Coalescence probability from theory for lineages sampled at deterministic mut-sel balance =" , coal.prob.theory.mut.sel.determin , "\n" )
	}
	return ( list ( pA = pA.vect , pa = pa.vect , coal.prob.sim.fix = coal.prob.sim.fix , coal.prob.theory.fix = coal.prob.theory.fix , pre.switch.freqs = pre.switch.freqs , mean.coal.prob.sim.mut.sel.drift = mean.coal.prob.sim.mut.sel.drift , coal.prob.sim.mut.sel.drift = coal.prob.sim.mut.sel.drift) )
	
}

SweepReplicates <- function ( N_A , N_a , s1 , s2 = NULL , h , switch = 0 , shut.mut.off = FALSE , stop.at.fix = TRUE , reuse.mutant.classes = TRUE , mu , gens , mut.index = 1001 , reps , suppress.output = TRUE , suppress.plot = TRUE ){

	recover()
	#prob.at.fix <- numeric ( reps )
	#prob.at.balance <- numeric ( reps * ( switch - 100 ) )
	for ( r in 1 : reps ) {
		
		temp <- SoftSweepSim ( N_A = N_A , N_a = N_a , s1 = s1 , s2 = s2 , switch = switch , stop.at.fix = stop.at.fix , reuse.mutant.classes = reuse.mutant.classes , mu = mu , gens = gens , suppress.output = suppress.output , suppress.plot = suppress.plot )
		prob.at.fix [ r ] <- temp[[1]]
		prob.at.balance [ ( ( switch - 100 ) * ( r - 1 ) + 1 ) : ( ( r * ( switch - 100 ) - 1 ) + 1 ) ] <- temp[[2]]
		cat ( "rep =" , r , "\n" )
		if ( r %% 20 == 0 ) {	
			
			mean.prob.at.balance.only.present <- mean ( prob.at.balance [ 1 : r * ( switch - 100 ) ] , na.rm = TRUE )
			prob.at.balance.zeros <- prob.at.balance
			prob.at.balance.zeros [ is.nan ( prob.at.balance.zeros ) ] <- 0
			mean.prob.at.balance.zeros <- mean ( prob.at.balance.zeros [ 1 : r * ( switch - 100 ) ] )
			prob.hartl.theory <- 1 / ( 1 + 2 * sum ( N_a , N_A ) * mu * 1 / (1 + abs(s1) ) )
			prob.fix.herm.penn <- 1 / ( 1 + 2 * sum ( N_a , N_A ) * mu )
			prob.at.fix.sim <- mean ( prob.at.fix [1 : r ] )
			cat ( "sim prob at balance =" , mean.prob.at.balance.only.present , "\n")
			cat ( "hartl prob at balance =" , prob.hartl.theory , "\n")
			cat ( "sim prob at balance w/ zeros =" , mean.prob.at.balance.zeros , "\n")
			cat ( "herm penn prob at fix =" , prob.fix.herm.penn , "\n")
			cat ( "sim prob at fix =" , prob.at.fix.sim , "\n")
		}
	}	
	
	return ( list ( prob.at.fix = prob.at.fix , prob.at.balance = prob.at.balance ) )
	
	
}

DiversityPlot <- function ( N , mu.b , mu.n , l , r , s , threshold ){
	
	
	diversity.hard <- diversity.soft <- coal.prob.soft <- coal.prob.hard <- numeric ()
	no.mut.prob <- 1 / ( 1 + N * mu.b ) 
	no.rec.prob <- exp ( -0.5 * ( log ( N * s ) / s ) * ( l - 1 ) * r )			
		
	coal.prob.soft <- no.rec.prob * no.mut.prob
	coal.prob.hard <- no.rec.prob
		
	diversity.soft <- ( ( 2 * mu.n ) / ( 1 / ( 2 * N ) ) ) * ( 1 - coal.prob.soft [ l ] )
	diversity.hard <- ( ( 2 * mu.n ) / ( 1 / ( 2 * N ) ) ) * ( 1 - coal.prob.hard [ l ] )
	
	plot ( x = c ( - rev ( l ) , l  ) , y = c ( rev ( diversity.soft ) , diversity.soft ) , type = "l" , ylim = c ( 0 , ( ( 2 * mu.n ) / ( 1 / ( 2 * N ) ) ) * 1.1 ) , xlab = "Position" , ylab = "Diversity" )
	lines ( x = c ( - rev ( l ) , l  ) , y = c ( rev ( diversity.hard ) , diversity.hard ) , type = "l" , col = "red")

	legend ( "bottomright" , legend = c ( "Soft Sweep" , "Hard Sweep" ) , lty = 1 , col = c ( "black" , "red" ) , bty = "n" )

	recover()
	
	neutral.diversity <- 2*mu.n/(1/(2*N))
	min ( which ( diversity.hard > 0.99*neutral.diversity ) )
	
}

#DiversityPlot ( N = 20000 , mu.b = 10^-5 , mu.n = 10^-8 , l <- 1:9000000 , r = 10^-8 , s = 0.001 )




first.run <- SoftSweepSim ( N_A = c ( 10000 ) , N_a = c ( 0 ) , s1 = -0.01 , s2 = 0.01 , switch = 1000 , stop.at.fix = TRUE , reuse.mutant.classes = TRUE , mu = 10^-5 , gens = 50000 , suppress.output = FALSE , sort.mutation.time = FALSE )

temp <- SweepReplicates ( N_A = c ( 200000 ) , N_a = c ( 0 ) , s1 = -0.001 , s2 = 0.01 , switch = 0 , shut.mut.off = FALSE , stop.at.fix = TRUE , reuse.mutant.classes = TRUE , mu = 10^-5 , gens = 100000 , suppress.output = FALSE , suppress.plot = TRUE , reps = 500 )

#temp2 <- SweepReplicates ( N_A = c ( 20000 ) , N_a = c ( 0 ) , s1 = -0.1 , s2 = 0.1 , switch = 1000 , shut.mut.off = FALSE , stop.at.fix = TRUE , reuse.mutant.classes = TRUE , mu = 10^-5 , gens = 100000 , suppress.output = FALSE , suppress.plot = TRUE , reps = 500 )

#temp2 <- SweepReplicates ( N_A = c ( 20000 ) , N_a = c ( 0 ) , s1 = -0.01 , s2 = 0.1 , switch = 1000 , shut.mut.off = FALSE , stop.at.fix = TRUE , reuse.mutant.classes = TRUE , mu = 10^-5 , gens = 100000 , suppress.output = TRUE , suppress.plot = TRUE , reps = 500 )


#temp3 <- SweepReplicates ( N_A = c ( 20000 ) , N_a = c ( 0 ) , s1 = -0.01 , s2 = 0.1 , switch = 1000 , shut.mut.off = FALSE , stop.at.fix = TRUE , reuse.mutant.classes = TRUE , mu = 10^-3 , gens = 100000 , suppress.output = TRUE , suppress.plot = TRUE , reps = 500 )

#temp4 <- SweepReplicates ( N_A = c ( 1000000 ) , N_a = c ( 0 ) , s1 = -0.2 , s2 = 0.1 , switch = 400 , shut.mut.off = FALSE , stop.at.fix = TRUE , reuse.mutant.classes = TRUE , mu = 10^-5 , gens = 100000 , suppress.output = TRUE , suppress.plot = TRUE , reps = 500 )