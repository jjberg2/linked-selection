##install.packages("randtoolbox")
##install.packages("ape")
setwd ( "~/Documents/Academics/StandingSweeps/" )
library("randtoolbox")
library("ape")
turn.on.recovers=FALSE

StructuredCoalescentSweep <- function ( N , s , h , dominance = FALSE , f , reps , n.tips , r , sim.distance , interval.width , no.sweep = FALSE , constant.freq = FALSE, cond.on.loss = TRUE , cond.on.fix = TRUE , make.plot = FALSE , build.seq = TRUE , standing.haps = TRUE , display.rep.count = TRUE , time.factor = 1 ) {
	#options ( error = recover )
	
	
	if ( constant.freq == FALSE ) {
	
		temp <- SweepFromStandingSim ( N = N , s = s , h = h , dominance = dominance , f = f , time.factor = time.factor , reps = reps , no.sweep = no.sweep, cond.on.loss=cond.on.loss , cond.on.fix = cond.on.fix , display.rep.count )
		frequencies <- temp [[ 1 ]]
		if ( no.sweep == FALSE ) {	
			sweep.start <- temp [[ 2 ]]
			new.freqs <- temp [[ 1 ]]
		} else if ( no.sweep == TRUE ){
			sweep.start <- temp [[ 2 ]]
			new.freqs <- frequencies [ , 1 : ncol ( frequencies ) ]
			fixation.time <- rep ( 0 , reps )
		}
	} else if ( constant.freq == TRUE ) {
		#recover()
		new.freqs <- matrix ( f , nrow = reps , ncol = 4*N*f *10 )
		fixation.time <- 0
	}
	num.lineages <- rep ( n.tips , reps )
	coal.times <- matrix ( 0 , nrow = reps , ncol = n.tips - 1 )	
	num.gens.simulated <- ncol ( new.freqs )
	
	i = 1
	## Coalscense
	while ( any ( num.lineages > 1 ) ) {
		no.mrca <- num.lineages != 1
		coal.probs <- rep ( 0 , reps )
		coal.probs [ no.mrca ] <- choose ( num.lineages [ no.mrca ] , 2 ) / ( 2 * N * new.freqs [ no.mrca , i ] )
		r.nums <- runif ( reps )
		if ( any ( r.nums < coal.probs ) ) {
			coals <- r.nums < coal.probs
			num.lineages [ coals ] <- num.lineages [ coals ] - 1
			coal.rows <- which ( coals )
			if ( length ( coal.rows ) > 1 & ncol ( coal.times ) > 1 ) {
				coal.cols <- apply ( coal.times[coal.rows,] , 1 , which.min )
			} else if ( length ( coal.rows ) == 1 & ncol ( coal.times ) > 1 ) {
				coal.cols <- which.min ( coal.times [ coal.rows , ] )
			} else if ( ncol ( coal.times ) == 1) {
				coal.cols <- rep ( 1 , length ( coal.rows ) )
			}
			
			coal.times [ (coal.cols-1) * reps + coal.rows ] <- i
		}
		i <- i + 1
	}
	mean.coalescence.times <- colMeans ( coal.times )
	sd.coalescence.times <- apply ( coal.times , 2 , sd )
	se.coalescence.times <- sd.coalescence.times / sqrt ( reps )
	trees <- BuildTrees ( coal.times = coal.times )
	

	for ( i in 1 : reps ) { 
		trees [[ i ]] [[ "freqs" ]] <- new.freqs[i,new.freqs[i,] != 0 ]
		if ( !no.sweep ) trees [[ i ]] [[ "sweep.start"]] <- sweep.start [ i ]
	}

	if ( build.seq == TRUE ) {
	#recover()
		temp <- RecombinationEvents ( trees = trees , coal.times = coal.times , r = r , sim.distance = sim.distance , n.tips = n.tips )
		trees <- temp [[ 1 ]]
		T.total <- temp [[ 2 ]]
	
		#recover()
		trees <- BuildOnOffHaps ( trees = trees , freqs = new.freqs , sim.distance = sim.distance , r = r , n.tips = n.tips , f = f  )
	
		hap.dist <- HapCountDistribution ( input = trees , r = r , sim.distance = sim.distance , interval.width = interval.width , f = f , N = N , make.plot )

		if ( standing.haps ) {
			standing.hap.dist <- StandingHapCountDist ( input = trees , r = r , sim.distance = sim.distance , interval.width = interval.width , f = f , N = N , make.plot )
		} else {
			standing.hap.dist <- NULL
		}
		
		
	} else {
		T.total <- NULL
		hap.dist <- NULL
		standing.hap.dist <- NULL
	}
	return ( list ( coal.times = coal.times , new.freqs = new.freqs , sweep.start = sweep.start , mean.coalescence.times = mean.coalescence.times , sd.coalescence.times = sd.coalescence.times , trees = trees , hap.dist = hap.dist , standing.hap.dist = standing.hap.dist , T.total = T.total , sim.distance.bp = sim.distance/r ) )

}

s_of_u <- function ( u ) {
	
	exp ( - integrate ( 
				f = function ( x ) {4 * N * s * x * ( 1 - x ) * ( x + h * ( 1 - 2 * x ) ) / ( x * ( 1 - x ) )} ,
				lower = 0 ,
				upper = u 
				)$value 
		)
	
}


SweepFromStandingSim <- function ( N , s , h , dominance = FALSE , f , reps , no.sweep, cond.on.loss , cond.on.fix , display.rep.count , time.factor = 1  ) {
	#options ( error = recover )
	delta.T <- 1 / ( time.factor * 2 * N )
	sweep.freq.matrix <- list ( rep ( f , reps ) )
	neutral.freq.matrix <- list ( rep ( f , reps ) )
	not.all.sweeps.fixed <- TRUE
	if ( f == 1 / ( 2 * N ) ) {
		not.all.neutral.fixed <- FALSE
		neutral.freq.matrix [[ 2 ]] <- rep ( 0 , reps )
	} else {
		not.all.neutral.fixed <- TRUE
	}
	#recover()
	i = 1
	while ( not.all.sweeps.fixed  | not.all.neutral.fixed ) {
		if ( not.all.sweeps.fixed ) {
			update <- rep ( 0 , reps )
			sweep.not.fixed <- sweep.freq.matrix [[ i ]] %% 1 != 0
			sweep.fixed <- sweep.freq.matrix [[ i ]] %% 1 == 0
			current.seg.freqs <- sweep.freq.matrix [[ i ]] [ sweep.not.fixed ]
			#recover()
			mu.S <- ifelse ( rep ( dominance , reps ) , 
						ifelse ( rep ( cond.on.fix , reps ) ,
									2 * N * s * current.seg.freqs * ( 1 - current.seg.freqs )  * ( current.seg.freqs + h * ( 1 - 2 * current.seg.freqs ) ) + sapply ( current.seg.freqs , s_of_u ) * ( 2 * current.seg.freqs * ( 1 - current.seg.freqs ) ) / sapply ( current.seg.freqs , function ( x ) integrate ( function ( y ) sapply ( y , s_of_u ) , lower = 0 , upper = x )$value ) ,
									2 * N * s * current.seg.freqs * ( 1 - current.seg.freqs )  * ( current.seg.freqs + h * ( 1 - 2 * current.seg.freqs ) )
									)
			,
						ifelse ( rep ( cond.on.fix , reps ) ,
									2 * N * s * current.seg.freqs * ( 1 - current.seg.freqs ) / tanh ( 2 * N * s * current.seg.freqs ) ,
									2 * N * s * current.seg.freqs * ( 1 - current.seg.freqs )
									)
					)
			sel <- mu.S * delta.T
			update [ sweep.not.fixed ] <- rnorm ( sum ( sweep.not.fixed ) , sel , sd = sqrt ( current.seg.freqs * ( 1 - current.seg.freqs ) * delta.T ) )
		#	sweep.drift.mag <- sqrt ( current.seg.freqs * ( 1 - current.seg.freqs ) * delta.T)
		#	plus.minus <- sample ( c ( 0 , 1 ) , sum ( sweep.not.fixed ) , replace = TRUE )
		#	drift.sweep <- ifelse ( plus.minus == 1 , sweep.drift.mag , -1 * sweep.drift.mag )
		#	update [ sweep.not.fixed ] <- sel + drift.sweep			
			sweep.freq.matrix [[ i + 1 ]] <- sweep.freq.matrix [[ i ]] + update
			sweep.fixed.one <- sweep.freq.matrix [[ i + 1 ]] > ( 1 - ( 1 / ( 2 * N ) ) )
			sweep.freq.matrix [[ i + 1 ]] [ sweep.fixed.one ] <- 1
			sweep.fixed.zero <- sweep.freq.matrix [[ i + 1 ]] < 1 / ( 2 * N )
			sweep.freq.matrix [[ i + 1 ]] [ sweep.fixed.zero ] <- 1 / ( 2 * N )
			not.all.sweeps.fixed <- any ( sweep.freq.matrix [[ i + 1 ]] %% 1 != 0 )

		
		}
		if ( not.all.neutral.fixed & f != 1 / ( 2 * N ) ) {	
			update <- rep ( 0 , reps )
			neutral.not.fixed <- neutral.freq.matrix [[ i ]] %% 1 != 0
			neutral.fixed <- neutral.freq.matrix [[ i ]] %% 1 == 0
			#neutral.drift.mag <- sqrt ( neutral.freq.matrix [[ i ]] [ neutral.not.fixed ] * ( 1 - neutral.freq.matrix [[ i ]] [ neutral.not.fixed ] ) * delta.T )
			#plus.minus <- sample ( c ( 0 , 1 ) , sum ( neutral.not.fixed ) , replace = TRUE )	
			#drift.neutral <- ifelse ( plus.minus == 1 , neutral.drift.mag , -1 * neutral.drift.mag )
			
			cond.mean <- ifelse ( rep ( cond.on.loss , reps ) ,
			 				- neutral.freq.matrix [[ i ]] [ neutral.not.fixed ] * delta.T ,
			 				0)
			
			drift.neutral <- rnorm ( sum ( neutral.not.fixed ) , cond.mean , sd = sqrt ( neutral.freq.matrix [[ i ]] [ neutral.not.fixed ] * ( 1 - neutral.freq.matrix [[ i ]] [ neutral.not.fixed ] ) * delta.T ) )
			update [ neutral.not.fixed ] <- drift.neutral
			neutral.freq.matrix [[ i + 1 ]] <- neutral.freq.matrix [[ i ]] + update
			neutral.fixed.one <- neutral.freq.matrix [[ i + 1 ]] > ( 1 - ( 1 / ( 2 * N ) ) )
			neutral.freq.matrix [[ i + 1 ]] [ neutral.fixed.one ] <- 1
			neutral.fixed.zero <- neutral.freq.matrix [[ i + 1 ]] < 1 / ( 2 * N )
			neutral.freq.matrix [[ i + 1 ]] [ neutral.fixed.zero ] <- 0	
			not.all.neutral.fixed <- any ( neutral.freq.matrix [[ i ]] %% 1 != 0 )
		
		}
		if ( i %% 500 == 0 & display.rep.count) {
				
				if ( i <= length ( neutral.freq.matrix) ) {
					lineages.remaining <- sum ( neutral.freq.matrix [[ i + 1 ]] %% 1 != 0 )
					my.freq <- max ( neutral.freq.matrix [[ i + 1 ]] [ neutral.freq.matrix [[ i + 1 ]] < 1 ] )
					cat ( "p = " , my.freq , ",  " , sep = "" )
					cat ( lineages.remaining , "not fixed \n")
				} else {
					lineages.remaining <- sum ( sweep.freq.matrix [[ i + 1 ]] %% 1 != 0 )
					my.freq <- max ( sweep.freq.matrix [[ i + 1 ]] [ sweep.freq.matrix [[ i + 1 ]] < 1 ] )
					cat ( "p = " , my.freq , ",  " , sep = "" )
					cat ( lineages.remaining , "sweeps not fixed \n")
				}
				
		}		
		if ( i == time.factor * 16 * N ){
			break
		}
		i = i + 1
	}
	sweep.freq.matrix <- matrix ( unlist ( sweep.freq.matrix ) , nrow = reps )
	sweep.keep <- seq ( 1 , ncol ( sweep.freq.matrix ) , by = time.factor )
	if ( ncol ( sweep.freq.matrix ) %in% sweep.keep ) {
 		sweep.freq.matrix <- sweep.freq.matrix [ , sweep.keep ]
 	} else {
	 	sweep.freq.matrix <- cbind ( sweep.freq.matrix [ , sweep.keep ] , 1 ) 		
 	}
	sweep.start <- apply ( sweep.freq.matrix , 1 , function ( x ) which.max ( x ) / time.factor )
	neutral.freq.matrix <- matrix ( unlist ( neutral.freq.matrix ) , nrow = reps )
	neutral.keep <- seq ( 1 , ncol ( neutral.freq.matrix ) , by = time.factor )
	if ( ncol ( neutral.freq.matrix ) %in% neutral.keep ) {
		neutral.freq.matrix <- neutral.freq.matrix [ , neutral.keep ]
	} else {
		neutral.freq.matrix <- cbind ( neutral.freq.matrix [ , neutral.keep ] , 0 )
	}
	if ( no.sweep == FALSE ) {
		freq.traj.list <- mapply ( 	function ( X , Y ) {
											#recover()
											fixation <- which.max ( Y )
											mutation <- sum ( X > 0 )
											freq <- c ( rev ( Y [ 2 : fixation ] ) , X [ 1 : mutation ] )
											return ( freq )
										} ,
										X = split ( neutral.freq.matrix , 1 : nrow ( neutral.freq.matrix ) ) , 
										Y = split ( sweep.freq.matrix , 1 : nrow ( sweep.freq.matrix ) )
							)
		freq.trajectories <- matrix ( 0 , ncol = max ( unlist ( lapply ( freq.traj.list , length ) ) ) , nrow = reps )
		for ( i in seq_len ( nrow ( freq.trajectories ) ) ) {
			freq.trajectories [ i , 1 : length ( freq.traj.list [[ i ]] ) ] <- freq.traj.list [[ i ]]
		}
		return ( list ( freq.trajectories , sweep.start ) )		
	} else {
		
		return ( list ( neutral.freq.matrix , sweep.start = 0 ) )
	}
}


BuildTrees <- function ( coal.times ){
	#recover()
	#library ( ape )
	if ( is.matrix ( coal.times ) == FALSE ) {
		n.trees <- 1
		n.tips <- length ( coal.times ) + 1
		coal.times <- matrix ( coal.times , nrow = 1 )
	} else {
		n.trees <- nrow ( coal.times )
		n.tips <- ncol ( coal.times ) + 1
	}
	trees <- list ( )
	for ( j in 1 : n.trees ) {
		edge <- matrix ( 0 , nrow = 2 * n.tips - 2 , ncol = 2 )
		edge.length <- numeric ( 2 * n.tips - 2 )
		edge [ 1 : n.tips , 2 ] <- 1 : n.tips
		nodes <- ( 2 * n.tips - 1 ) : ( n.tips + 1 )
		node.depth <- numeric ( 2 * n.tips - 1 )
		Nnode <- n.tips - 1
		tip.label <- character ( n.tips )
		for ( l in 1 : length ( tip.label ) ){
			tip.label [ l ] <- paste ( "t" , l , sep = "")
		}
		k = 1
		for ( i in nodes ) {
			extant.lineages <- edge [ edge [ , 2] != 0 & edge [ , 1 ] == 0 , 2 ]
			coalescing.lineages <- sort ( sample ( extant.lineages , 2 , replace = FALSE ) )
			coal.index <- which ( edge [ , 2 ] %in% coalescing.lineages )
			edge [ coal.index , 1 ] <- i
			if ( i != tail ( nodes , 1 ) ) {
				edge [ i - 1 , 2 ] <- i
			}
			node.depth [ i ] <- coal.times [ j , k ] 
			edge.length [ coal.index ] <- coal.times [ j , k ] - node.depth [ coalescing.lineages ]
			k = k + 1
		}
		a.tree <- list ( edge = edge , edge.length = edge.length , tip.label = tip.label , Nnode = Nnode )
		class ( a.tree ) <- "phylo"
		my.tree <- list ( tree = a.tree , node.depth = node.depth )
		
		trees [[ j ]] <- my.tree
	}
	
	
	return ( trees )
	
}

RecombinationEvents <- function ( trees , coal.times , r , sim.distance , n.tips ) {
if(turn.on.recovers)	recover()
	if ( n.tips > 2 ) {
		internodes <- matrix ( nrow = nrow ( coal.times ) , ncol = n.tips - 1 )
		internodes [ , 1 ] <- coal.times [ , 1 ]
		for ( i in 2 : ( n.tips - 1 ) ) {
			internodes [ , i ] <- coal.times [ , i ] - coal.times [ , i - 1 ]
		}
	} else if ( n.tips == 2 ) {
		internodes <- coal.times
	}	
	
	T.total <- numeric ( length ( trees ) )
	#recover()
	cat ( "Laying down recombination events. \n \n")
	pb <- txtProgressBar ( min = 0 , max = length ( trees ) , style = 3 )
	for ( j in 1 : length ( trees ) ) {
		T.total [ j ] <- sum ( ( n.tips : 2 ) * internodes [ j , ] )
		sim.distance.bp <- sim.distance/r
		rec.right.temp <- data.frame ( sequence.location = 0 , branch = 0 , rec.depth = 0 )
		rec.left.temp <- data.frame ( sequence.location = 0 , branch = 0 , rec.depth = 0 )
		edges <- 1 : tail ( trees [[ j ]] [[ 1 ]] [[ 1 ]] [ , 2 ] , 1 )
		if ( ncol ( coal.times ) > 1 ) {
			edge.lengths <- c ( trees [[ j ]] [[ 1 ]] [[ 2 ]] [  1 : ( ( length ( edges ) + 1 ) / 2 ) ] , 0 , trees [[ j ]] [[ 1 ]] [[ 2 ]] [ ( ( ( length ( edges ) + 1 ) / 2 ) + 1 ) : ( length ( edges ) - 1 ) ] )
		} else {
			edge.lengths <- trees [[ j ]] [[ 1 ]] [[ 2 ]]
		}
		i = 1
		while ( rec.right.temp [ i , 1 ] < sim.distance.bp ) {
			rec.right.temp [ i + 1 , 1 ] <- rec.right.temp$sequence.location [ i ] + round ( rexp ( 1 , r * T.total [ j ] ) )
			rec.right.temp [ i + 1 , 2 ] <- sample ( edges , 1 , prob = edge.lengths )
			rec.right.temp [ i + 1 , 3 ] <- trees [[ j ]] [[2]] [ rec.right.temp [ i + 1 , 2 ] ] + sample ( seq ( 1 , edge.lengths [ rec.right.temp [ i + 1 , 2 ] ] - 1) , 1 )
			i = i + 1
		}
		
		i = 1
		while ( rec.left.temp [ i , 1 ] < sim.distance.bp ) {
			rec.left.temp [ i + 1 , 1 ] <- rec.left.temp$sequence.location [ i ] + round ( rexp ( 1 , r * T.total [ j ] ) )
			rec.left.temp [ i + 1 , 2 ] <- sample ( edges , 1 , prob = edge.lengths )
			rec.left.temp [ i + 1 , 3 ] <- trees [[ j ]] [[2]] [ rec.left.temp [ i + 1 , 2 ] ] + sample ( seq ( 1 , edge.lengths [ rec.left.temp [ i + 1 , 2 ] ] - 1 ) , 1 )
			i = i + 1
		}
		#recover()
		trees [[ j ]] [[ "T.total" ]] <- T.total [ j ]
		trees [[ j ]] [[ "rec.events" ]] <- recombination <-  list ( rec.right = rec.right.temp [ -c ( 1 , nrow ( rec.right.temp ) ), ] , rec.left = rec.left.temp [ -c ( 1 , nrow ( rec.left.temp ) ) , ] )
		setTxtProgressBar ( pb, j )
	}
	close ( pb )	
	
	return ( list ( trees, T.total ) )

}

BuildOnOffHaps <- function ( trees , freqs , r , sim.distance , n.tips , f , fixation.time ) {
	
	sim.distance.bp <- sim.distance / r
	#recover()
	cat ( "Building Haplotypes. \n \n")
	pb <- txtProgressBar ( min = 0 , max = length ( trees ) , style = 3 )
	for ( j in 1 : length ( trees ) ) {
		rec.right <- trees [[ j ]]$rec.events$rec.right
		rec.left <- trees[[ j ]]$rec.events$rec.left
		
	
		## build right side haplotype ##
		event.order <- order ( rec.right [ , 3 ] , decreasing = TRUE )
		right.sequence.temp <- matrix ( 0 , nrow = n.tips , ncol = nrow ( rec.right ) + 1 )
		sub.trees <- prop.part ( trees [[ j ]]$tree )
		to.remove <- numeric ( )
		h = 1
		l = 2
		if ( nrow ( rec.right ) != 0 ) {
			for ( i in event.order ) {
				this.event <- data.frame ( rec.right [ i , ] , hap.ID = h )
				if ( this.event$rec.depth == 0 ) {
					break
				} else {
					my.freq <- trees [[ j ]] [[ 3 ]] [ this.event$rec.depth ]
				}
				rec.roll <- runif ( 1 )
				if ( rec.roll < ( 1 - my.freq ) ) {
					if ( this.event$branch > n.tips ) {
						tips <- unlist ( sub.trees [ this.event$branch - n.tips ] )
						right.sequence.temp [ tips , ( i + 1 ) : ncol ( right.sequence.temp ) ] <- h
					} else {
						tip <- this.event$branch
						right.sequence.temp [ tip , ( i + 1 )  : ncol ( right.sequence.temp ) ] <- h
					}
				l = l + 1
				h = h + 1	
				}
			}
		
			
			for ( i in 2 : ncol ( right.sequence.temp ) ) {
				if ( length ( unique ( right.sequence.temp [ , i ] ) ) == length ( unique ( right.sequence.temp [ , i - 1 ] ) ) ) {
					to.remove [ length ( to.remove ) + 1 ] <- i
				}
			}
			
		} 
		if ( sum ( right.sequence.temp ) == 0 ) {
			right.sequence <- right.sequence.temp
			rec.right.off.background <- rec.right
		} else if ( length ( to.remove ) != 0 ){			
			right.sequence <- right.sequence.temp [ , -to.remove ]
			right.sequence <- MakeHapsPretty ( right.sequence )
			rec.right.off.background <- rec.right [ - ( to.remove - 1 ) , ]
		} else {
			right.sequence <- right.sequence.temp
			right.sequence <- MakeHapsPretty ( right.sequence )			
			rec.right.off.background <- rec.right
		}
		


		## build left side haplotype ## 
		event.order <- order ( rec.left [ , 3 ] , decreasing = TRUE )
		left.sequence.temp <- matrix ( 0 , nrow = n.tips , ncol = nrow ( rec.left ) + 1 )
		sub.trees <- prop.part( trees [[ j ]] [[ 1 ]] )
		to.remove <- numeric ( )
		l = 2
		if ( nrow ( rec.left ) != 0 ) {
			for ( i in event.order ) {
				this.event <- data.frame ( rec.left [ i , ] , hap.ID = h )
				if ( this.event$rec.depth == 0 ) {
					break
				} else {
					my.freq <- trees [[ j ]] [[ 3 ]] [ this.event$rec.depth ]
				}
				rec.roll <- runif ( 1 )
				if ( rec.roll < ( 1 - my.freq ) ) {
					if ( this.event$branch > n.tips ) {
						tips <- unlist ( sub.trees [ this.event$branch - n.tips ] )
						left.sequence.temp [ tips , ( i + 1 ) : ncol ( left.sequence.temp ) ] <- h
					} else {
						tip <- this.event$branch
						left.sequence.temp [ tip , ( i + 1 )  : ncol ( left.sequence.temp ) ] <- h
					}
				l = l + 1
				h = h + 1	
				}
			}
			#recover()
			for ( i in 2 : ncol ( left.sequence.temp ) ) {
				if ( length ( unique ( left.sequence.temp [ , i ] ) ) == length ( unique ( left.sequence.temp [ , i - 1 ] ) ) ) {
					to.remove [ length ( to.remove ) + 1 ] <- i
				}
			}
		} 
		if ( sum ( left.sequence.temp ) == 0 ) {
			left.sequence <- left.sequence.temp
			rec.left.off.background <- rec.left
		} else if ( length ( to.remove ) != 0 ) {
			left.sequence <- left.sequence.temp [ , -to.remove ]
			left.sequence <- MakeHapsPretty ( left.sequence )
			rec.left.off.background <- rec.left [ - ( to.remove - 1 ) , ]
		} else {
			left.sequence <- left.sequence.temp
			left.sequence <- MakeHapsPretty ( left.sequence )
			rec.left.off.background <- rec.left
		}

		setTxtProgressBar ( pb, j )
		trees [[ j ]] [[ "sequence.structure" ]] <- list ( right.seq = right.sequence , left.seq = left.sequence )
		trees [[ j ]] [[ "rec.events.off.background" ]] <- list ( rec.right.off.background = rec.right.off.background , rec.left.off.background = rec.left.off.background )
		trees [[ j ]] [[ "sim.distance.bp" ]] <- sim.distance.bp
	}
	close ( pb )
	return ( trees )
}


HapCountDistribution <- function ( input , r = 10^-8 , sim.distance , interval.width = 1000 , f , N , make.plot ) {
	#recover()
	sim.distance.bp <- sim.distance / r 
	intervals <- seq ( 0 , sim.distance.bp , interval.width )
	n.tips <- length ( input [[ 1 ]]$tree$tip.label )
	reps <- length ( input )
	# number of rows in "sequence" matrix = number of samples
	if ( turn.on.recovers ) {
		recover()
	}
	
	no.sing.haps.right <- no.sing.haps.left <- matrix ( nrow = length ( input ) , ncol = length ( intervals ) )
	n.haps.right <- n.haps.left <- matrix ( nrow = length ( input ) , ncol = length ( intervals ) )
	sing.haps.left <- sing.haps.right <- matrix ( nrow = length ( input ) , ncol = length ( intervals ) )
	#recover()
	cat ( "Counting up haplotypes. \n \n")
	pb <- txtProgressBar ( min = 0 , max = length ( intervals ) , style = 3 )
	for ( i in 1 : length ( intervals ) ) {
		
		k <- intervals [ i ]

		if ( k == 0 ) {
			
			# there is only one haplotype at the selected sight		
			n.haps.right [ , i ] <- n.haps.left [ , i ] <- 1
			no.sing.haps.right [ , i ] <- no.sing.haps.left [ , i ] <- 1
			
		} else {
			
			# now we loop through the simulated data to work out the number of haplotypes at various intervals away from the selected sight
			#recover ( )
			for ( j in 1 : length ( input ) ) {
				
				my.seqs <- input [[ j ]] $ sequence.structure
				my.rec.events <- input [[ j ]] $ rec.events.off.background
				
				# right side
				if ( sum ( my.rec.events$rec.right.off.background$sequence.location < k ) != 0 ) {
					last.rec.event <- sum ( my.rec.events$rec.right.off.background$sequence.location < k )
					n.haps.right [ j , i ] <-  length ( unique ( my.seqs$right.seq [ , last.rec.event + 1 ] ) )
					no.sing.haps.right [ j , i ] <- sum ( table ( my.seqs$right.seq [ , last.rec.event + 1 ] ) > 1 )
					sing.haps.right [ j , i ] <- sum ( table ( my.seqs$right.seq [ , last.rec.event + 1 ] ) == 1 )
				} else {
					n.haps.right [ j , i ] <- 1
					no.sing.haps.right [ j , i ] <- 1
					sing.haps.right [ j , i ] <- 0
				}
				
				# left.side
				if ( sum ( my.rec.events$rec.left.off.background$sequence.location < k ) != 0 ) {
					last.rec.event <- sum ( my.rec.events$rec.left.off.background$sequence.location < k )
					n.haps.left [ j , i ] <-  length ( unique ( my.seqs$left.seq [ , last.rec.event + 1 ] ) )
					no.sing.haps.left [ j , i ] <- sum ( table ( my.seqs$left.seq [ , last.rec.event + 1 ] ) > 1 )
					sing.haps.left [ j , i ] <- sum ( table ( my.seqs$left.seq [ , last.rec.event + 1 ] ) == 1 )
				} else {
					n.haps.left [ j , i ] <- 1
					no.sing.haps.left [ j , i ] <- 1
					sing.haps.left [ j , i ] <- 0
				}
			
			}	
		}	
	
		setTxtProgressBar(pb, i)
		
	}
	close(pb)

	#recover()
	n.haps <- rbind ( n.haps.right , n.haps.left )
	no.sing.haps <- rbind ( no.sing.haps.right , no.sing.haps.left )
	
	hap.counts.by.interval <- apply ( n.haps , 2 , function ( x ) table ( factor ( x , 1 : n.tips ) ) )
	hap.count.freqs.by.interval <- apply ( hap.counts.by.interval , 2 , function ( x ) x / nrow ( n.haps ) )
	
	no.sing.hap.counts.by.interval <- apply ( no.sing.haps , 2 , function ( x ) table ( factor ( x , 0 : n.tips ) ) )
	no.sing.hap.count.freqs.by.interval <- apply ( no.sing.hap.counts.by.interval , 2 , function ( x ) x / nrow ( no.sing.haps ) )
	
	if ( make.plot ) {
		MakeHapPlots ( hap.count.freqs.by.interval , N , f , sim.distance , r = 10^-8 , interval.width = 1000 )
	}
	
	return ( list ( hap.count.freqs.by.interval = hap.count.freqs.by.interval , no.sing.hap.count.freqs.by.interval = no.sing.hap.count.freqs.by.interval , n.haps = n.haps , no.sing.haps = no.sing.haps ) )
	
}


StandingHapCountDist <- function ( input , r = 10^-8 , sim.distance , interval.width = 1000 , f , N , make.plot ) {
	#recover()
	sim.distance.bp <- sim.distance / r 
	intervals <- seq ( 0 , sim.distance.bp , interval.width )
	n.tips <- length ( input [[ 1 ]]$tree$tip.label )
	reps <- length ( input )
	# number of rows in "sequence" matrix = number of samples
	if ( turn.on.recovers ) {
		recover()
	}
	n.haps.right <- n.haps.left <- matrix ( nrow = length ( input ) , ncol = length ( intervals ) )
	#recover()
	cat ( "Counting up haplotypes. \n \n")
	pb <- txtProgressBar ( min = 0 , max = length ( intervals ) , style = 3 )
	for ( i in 1 : length ( intervals ) ) {
		
		k <- intervals [ i ]

		if ( k == 0 ) {
			
			# there is only one haplotype at the selected sight		
			n.haps.right [ , i ] <- n.haps.left [ , i ] <- 1
			
		} else {
			
			# now we loop through the simulated data to work out the number of haplotypes at various intervals away from the selected sight
			#recover ( )
			for ( j in 1 : length ( input ) ) {
				
				my.seqs <- input [[ j ]] $ sequence.structure
				my.rec.events <- input [[ j ]] $ rec.events
				my.rec.events.off <- input [[ j ]] $ rec.events.off.background
				
				
				# right side
				sweep.recs <- my.rec.events$rec.right$rec.depth < input [[ j ]]$sweep.start 
				site.side <- my.rec.events$rec.right$sequence.location < k
				sweep.killed.branches <- my.rec.events$rec.right$branch [ site.side & sweep.recs ]
				sweep.killed.branches <- unique ( unlist ( sapply ( unique ( sweep.killed.branches ) , function ( x ) GetTips ( x , n.tips , input [[ j ]]$tree$edge) ) ) )
				# if ( any ( sweep.killed.branches > n.tips ) ) {
					# internal.recs <- sweep.killed.branches [ sweep.killed.branches > n.tips ]
					# for ( i in internal.recs ) {
						# temp <- extract.clade ( input [[ j ]]$tree , i )
						# my.tips <- as.numeric ( unlist ( lapply ( strsplit ( temp$tip.label , "t" ) , function ( x ) x [ 2 ] ) ) )
						# sweep.killed.branches <- c ( sweep.killed.branches , my.tips )
					# }
					# sweep.killed.branches <- unique ( sweep.killed.branches [ sweep.killed.branches <= n.tips ] )
				# }
				if ( sum ( my.rec.events.off$rec.right.off.background$sequence.location < k ) != 0 ) {
					
					last.rec.event <- sum ( my.rec.events.off$rec.right.off.background$sequence.location < k )
					my.tab <- table ( my.seqs$right.seq [ unlist ( ifelse ( is.null ( sweep.killed.branches ) , list(seq_len(n.tips)) ,  list(-sweep.killed.branches) ) ) , last.rec.event + 1 ] )
					n.haps.right [ j , i ] <- sum ( my.tab > 1 )
				} else {
					n.haps.right [ j , i ] <- 1
				}
				
				# left.side
				sweep.recs <- my.rec.events$rec.left$rec.depth < input [[ j ]]$sweep.start 
				site.side <- my.rec.events$rec.left$sequence.location < k
				sweep.killed.branches <- my.rec.events$rec.left$branch [ site.side & sweep.recs ]
				sweep.killed.branches <- unique ( unlist ( sapply ( unique ( sweep.killed.branches ) , function ( x ) GetTips ( x , n.tips , input [[ j ]]$tree$edge) ) ) )
				# if ( any ( sweep.killed.branches > n.tips ) ) {
					# internal.recs <- sweep.killed.branches [ sweep.killed.branches > n.tips ]
					# for ( i in internal.recs ) {
						# temp <- GetTips ( i , n.tips , input[[j]]$tree$edge )
						# my.tips <- as.numeric ( unlist ( lapply ( strsplit ( temp$tip.label , "t" ) , function ( x ) x [ 2 ] ) ) )
						# sweep.killed.branches <- c ( sweep.killed.branches , my.tips )
					# }
					# sweep.killed.branches <- unique ( sweep.killed.branches [ sweep.killed.branches <= n.tips ] )
				# }
				if ( sum ( my.rec.events.off$rec.left.off.background$sequence.location < k ) != 0 ) {
					last.rec.event <- sum ( my.rec.events.off$rec.left.off.background$sequence.location < k )
					my.tab <- table ( my.seqs$left.seq [ unlist ( ifelse ( is.null ( sweep.killed.branches ) , list(seq_len(n.tips)) ,  list(-sweep.killed.branches) ) ) , last.rec.event + 1 ] )
					n.haps.left [ j , i ] <-  sum ( my.tab > 1 )
				} else {
					n.haps.left [ j , i ] <- 1
				}
			
			}	
		}	
	
		setTxtProgressBar(pb, i)
		
	}
	close(pb)

	#recover()
	n.haps <- rbind ( n.haps.right , n.haps.left )
	
	hap.counts.by.interval <- apply ( n.haps , 2 , function ( x ) table ( factor ( x , 0 : n.tips ) ) )
	hap.count.freqs.by.interval <- apply ( hap.counts.by.interval , 2 , function ( x ) x / nrow ( n.haps ) )
	
	if ( make.plot ) {
		MakeHapPlots ( hap.count.freqs.by.interval , N , f , sim.distance , r = 10^-8 , interval.width = 1000 )
	}
	
	return ( list ( hap.count.freqs.by.interval = hap.count.freqs.by.interval , n.haps = n.haps ) )
	
}


SingHapCountDistribution <- function ( input , r = 10^-8 , sim.distance , interval.width = 1000 , f , N , make.plot ) {
	#recover()
	sim.distance.bp <- sim.distance / r 
	intervals <- seq ( 0 , sim.distance.bp , interval.width )
	n.tips <- length ( input [[ 1 ]]$tree$tip.label )
	reps <- length ( input )
	# number of rows in "sequence" matrix = number of samples
	if ( turn.on.recovers ) {
		recover()
	}
	
	sing.haps.right <- sing.haps.left <- matrix ( nrow = length ( input ) , ncol = length ( intervals ) )
	n.haps.right <- n.haps.left <- matrix ( nrow = length ( input ) , ncol = length ( intervals ) )
	#recover()
	cat ( "Counting up haplotypes. \n \n")
	pb <- txtProgressBar ( min = 0 , max = length ( intervals ) , style = 3 )
	for ( i in 1 : length ( intervals ) ) {
		
		k <- intervals [ i ]

		if ( k == 0 ) {
			
			# there is only one haplotype at the selected sight		
			n.haps.right [ , i ] <- n.haps.left [ , i ] <- 1
			no.sing.haps.right [ , i ] <- no.sing.haps.left [ , i ] <- 1
			
		} else {
			
			# now we loop through the simulated data to work out the number of haplotypes at various intervals away from the selected sight
			#recover ( )
			for ( j in 1 : length ( input ) ) {
				
				my.seqs <- input [[ j ]] $ sequence.structure
				my.rec.events <- input [[ j ]] $ rec.events.off.background
				
				# right side
				if ( sum ( my.rec.events$rec.right.off.background$sequence.location < k ) != 0 ) {
					last.rec.event <- sum ( my.rec.events$rec.right.off.background$sequence.location < k )
					
					n.haps.right [ j , i ] <-  length ( unique ( my.seqs$right.seq [ , last.rec.event + 1 ] ) )
					no.sing.haps.right [ j , i ] <- sum ( table ( my.seqs$right.seq [ , last.rec.event + 1 ] ) > 1 )
				} else {
					n.haps.right [ j , i ] <- 1
					no.sing.haps.right [ j , i ] <- 1
				}
				
				# left.side
				if ( sum ( my.rec.events$rec.left.off.background$sequence.location < k ) != 0 ) {
					last.rec.event <- sum ( my.rec.events$rec.left.off.background$sequence.location < k )
					n.haps.left [ j , i ] <-  length ( unique ( my.seqs$left.seq [ , last.rec.event + 1 ] ) )
					no.sing.haps.left [ j , i ] <- sum ( table ( my.seqs$left.seq [ , last.rec.event + 1 ] ) > 1 )
				} else {
					n.haps.left [ j , i ] <- 1
					no.sing.haps.left [ j , i ] <- 1
				}
			
			}	
		}	
	
		setTxtProgressBar(pb, i)
		
	}
	close(pb)

	#recover()
	n.haps <- rbind ( n.haps.right , n.haps.left )
	no.sing.haps <- rbind ( no.sing.haps.right , no.sing.haps.left )
	
	hap.counts.by.interval <- apply ( n.haps , 2 , function ( x ) table ( factor ( x , 1 : n.tips ) ) )
	hap.count.freqs.by.interval <- apply ( hap.counts.by.interval , 2 , function ( x ) x / nrow ( n.haps ) )
	
	no.sing.hap.counts.by.interval <- apply ( no.sing.haps , 2 , function ( x ) table ( factor ( x , 0 : n.tips ) ) )
	no.sing.hap.count.freqs.by.interval <- apply ( no.sing.hap.counts.by.interval , 2 , function ( x ) x / nrow ( no.sing.haps ) )
	
	if ( make.plot ) {
		MakeHapPlots ( hap.count.freqs.by.interval , N , f , sim.distance , r = 10^-8 , interval.width = 1000 )
	}
	
	return ( list ( hap.count.freqs.by.interval = hap.count.freqs.by.interval , no.sing.hap.count.freqs.by.interval = no.sing.hap.count.freqs.by.interval , n.haps = n.haps , no.sing.haps = no.sing.haps ) )
	
}





MakeHapPlots <- function ( hap.count.freqs.by.interval , N , f , sim.distance , r = 10^-8 , interval.width = 1000,plot.cumulative=TRUE,do.legend=FALSE) {
	library ( wesanderson)
	#recover()
	sim.distance.bp <- sim.distance / r 
	intervals <- seq ( 0 , sim.distance.bp , interval.width )
	n.tips <- max ( as.numeric ( rownames(hap.count.freqs.by.interval ) ) )
	if(plot.cumulative){	cum.probs <- rbind ( 0 , apply ( hap.count.freqs.by.interval , 2 , cumsum ) )}
	if(!plot.cumulative){ cum.probs <- rbind ( 0 ,hap.count.freqs.by.interval)}

	ewens.dist.matrix <- matrix ( nrow = n.tips , ncol = length ( intervals ) )

	#stirling.numbers <- StirlingNumbers ( n = n.tips ) [ n.tips , ]
	 for ( i in 1 : length ( intervals ) ) {
		
		 if ( i == 1 & intervals [ 1 ] == 0 ) {
			
			 ewens.dist.matrix [ , i ] <- c ( 1 , rep ( 0 , n.tips - 1 ) )
			
		 } else { 
		
			 ewens.dist.matrix [ , i ] <- EwensDist ( n = n.tips , N = N , r = r , distance = intervals [ i ] , f = f  ) [ n.tips , ]
			
		}
		
	 }
	#recover()
#recover()
	if(plot.cumulative){ ewens.cum.probs <-  apply ( ewens.dist.matrix , 2 , cumsum )}
	if(!plot.cumulative){ewens.cum.probs <-ewens.dist.matrix; }
#	recover()
	matplot ( 
		t ( ewens.dist.matrix ) , 
		type = "n" , 
		lty = 1 , 
		lwd = 0.7 , 
		col = "black" , 
		ylab = "Probability" , 
		xlab = "kb" , 
		#main = paste ( n.tips , "Lineages in a Sweep from f =" , f , "at s =" , s , "," , reps , "Reps" ) , 
		bty = "n" ,
		ylim = c ( 0 , 1 ) , 
		xaxt = "n"
	)
	axis ( 1 , seq ( 1 , ncol ( ewens.dist.matrix ) , by = 10e5/interval.width ) , 4*N*r * 1000*seq ( 0 , max ( intervals , 1 )/1000  , by = 1000 )  )  # seq ( 0 , tail ( intervals , 1 ) / 1000 , by = 1000 )  )
	#recover()
	#col.vect <- rainbow ( n.tips , s = 0.8  , v = 1 , start = 1/40 , end = 4/6  )
	#col.vect <- brewer.pal ( 10 , "Set3" )
	col.vect <- wes_palette ( "FantasticFox" , 10 , type = "continuous")
	if(do.legend) legend("topright", legend=1:n.tips, lty=1,col= col.vect,lwd=2, bty = "n")
	for ( i in  ( nrow ( cum.probs ) - 1 ):1 ) {
			#i = i + 1
			X.ax <- 1:ncol ( cum.probs ) #which ( cum.probs [ i , ] != cum.probs [ i + 1 , ] )
	if(plot.cumulative){	Y.ax1 <- cum.probs [ i , X.ax ]}
	if(!plot.cumulative){ Y.ax1 <- rep(0,ncol ( cum.probs ) )	}
			if(!plot.cumulative){ 
				lines(X.ax,cum.probs [ i + 1 , X.ax ], col = col.vect [ i ],lwd=2 )
				lines (ewens.cum.probs[i,], col = col.vect [ i ],lwd=2,lty=2 ) 
				}
			Y.ax2 <- cum.probs [ i + 1 , X.ax ]
			if(plot.cumulative) polygon ( x = c ( X.ax , rev ( X.ax ) ) , y = c ( Y.ax1 , rev ( Y.ax2 ) ) , lty = 0 , col = col.vect [ i ] )					
	}

	
	if(plot.cumulative){ 
		ewens.cum.probs <- ewens.cum.probs [ - nrow ( ewens.cum.probs ) , ]
		apply ( ewens.cum.probs , 1 , function ( x ) lines ( x , lty = 1 , lwd = 0.8 ) )
	}


}

GetTips <- function ( branch , n.tips , edges ) {
	#recover()
	if ( branch <= n.tips ) {
		return ( branch )
	}
	subtend <- edges [ edges [  , 1 ] == branch , 2 ]
	if ( all ( subtend <= n.tips ) ) {
		return ( subtend )
	} else {
		sapply ( subtend , function ( x ) GetTips ( x , n.tips , edges ) )
	}
		
}

StirlingNumbers <- function ( n ) {
	
	nums <- matrix ( NA , nrow = n + 1 , ncol = n + 1 )
	nums [ 1 ,  ] <- nums [  , 1 ] <- 0
	nums [ 1 , 1 ] <- 1
	for ( i in 1 : n ) {
		for ( k in 1 : n ) {	
			nums [ i + 1 , k + 1 ] <- ( i -1 ) * nums [ i  , k + 1 ] + nums [ i , k  ]
		}	
	}
	return ( nums [ -1 , - 1 ] )
}

EwensDist <- function ( n , N , r , distance , f ) {
	#recover()	
	param <- 4 * N * r * distance * f * ( 1- f )
	denom  <- cumprod ( param + 0 : ( n - 1 ) )
	stirling.numbers <- StirlingNumbers ( n )
	ewens.dist <- t ( param^(1:n) * t ( stirling.numbers / denom ) )
	return ( ewens.dist ) 

}


EwensCondDist <- function ( n , N , r , distance , f ) {
	#recover()	
	param <- 4 * N * r * distance * f * ( 1- f )
	my.prod  <- cumprod ( param + 0 : ( n - 1 ) )
	denom <- my.prod - param*factorial ( 0:( n - 1 ) )
	stirling.numbers <- StirlingNumbers ( n )
	ewens.dist <- t ( param^(1:n) * t ( stirling.numbers / denom ) )
	ewens.dist [ , 1 ] <- ewens.dist [ 1 , ] <- 0
	return ( ewens.dist ) 

}



MakeHapsPretty <- function ( seqs ) {
	if ( !is.numeric ( nrow ( seqs ) ) | !is.numeric ( ncol ( seqs ) ) ) recover()
	new.seqs <- matrix ( 0 , nrow = nrow ( seqs ) , ncol = ncol ( seqs ) )
	for ( i in 2 : ncol ( seqs ) ) {	
		j <- i - 1
		new.ids <- unique ( seqs [ seqs [ , i ] %in% seqs [ , i - 1 ] == FALSE , i ])
		for ( x in new.ids ){
			last.hap <- unique ( seqs [ seqs [ , i ] == x , i - 1 ] )
			if ( sum ( seqs [ , i ] == x ) != sum ( seqs [ , i - 1 ] == last.hap ) ) {
				new.hap <- x
				break 
			}	
		}
		new.seqs [ seqs [ , i ] == new.hap , i : ncol ( new.seqs ) ] <- j
	}
	
	return ( new.seqs )
	
}


EffectiveS <- function ( N , s , f ) {
	
	exp ( - lambertW_base ( ( 4*N*s*f^2 - 4*N*s*f - log ( 1/f ) ) / ( 2*N * s ) ) ) / (2 * N) 
	
}


EffectiveS2 <- function ( N , s , f ) {
	
	s*log ( 2*N ) / ( log ( 1 / f ) - 4 * ( f - 1 ) *f * N * s )
	
}


GetTermInExp <- function ( sims ) {

	my.exp <- numeric ()
	for ( i in 1 : nrow ( sims$coal.times ) ) {
		
		my.exp [ i ] <- 2*sum ( 1-sims$new.freqs [ 1 , 1 : sims$coal.times [ i , ] ] )
		
	}
	mean ( my.exp )
		
}



if(FALSE){
fs <- c ( 1/20000  , 0.01 , 0.05 , 0.1 )
ss <- c ( 0.001 , 0.01 , 0.05 )
fands <- expand.grid ( fs , ss )
colnames ( fands ) <- c ( "f" , "s")
temp <- apply ( fands , 1 , function ( x ) StructuredCoalescentSweep ( N = 10000 , s = x[2] , f = x[1] , reps = 200 , n.tips = 12 , r = 10^-8 , sim.distance = 0.01 , interval.width = 1000 , no.sweep = FALSE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = FALSE ,  time.factor = 1 ) )

#function to get haplotype distribution plots from function output

temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.05 , dominance = FALSE , f = 0.05 , reps = 1000 , n.tips = 10 , r = 10^-8 , sim.distance = 0.02 , interval.width = 1000 , no.sweep = T , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = FALSE , display.rep.count = TRUE ,   standing.haps = FALSE , time.factor = 1 )

sweep.coals <- rowSums ( temp$coal.times < temp$sweep.start )
when <- mapply ( function ( x , y ) y [ seq ( 0 , x ) ] , x = sweep.coals , y = split ( temp$coal.times , 1 : 1000 ) )

MakeHapPlots ( temp$hap.dist$hap.count.freqs.by.interval , N = 10000 , f = 0.025 , sim.distance = 0.02 )





MakeHapPlots ( temp$hap.dist$hap.count.freqs.by.interval , N = 10000, f = 1/20000, sim.distance = 0.01)

SequenceIBDPlots <- function ( trees ) {
	
	#recover()
	seq.structure <- trees$sequence.structure
	seq.resort <- do.call(what = order, as.data.frame(seq.structure))
	seq.structure <- lapply ( seq.structure , function ( x ) x [ seq.resort , ] )
	rec.points <- trees$rec.events.off.background
	scaled.rec.points <- list ()
	scaled.rec.points$right <- c ( 0 , rec.points$rec.right.off.background$sequence.location / trees$sim.distance , 1 )
	scaled.rec.points$left <- - c ( 0 , rec.points$rec.left.off.background$sequence.location / trees$sim.distance , 1 )
	#my.cols <- rainbow ( max ( unlist ( seq.structure ) ) + 1 , alpha = 0.7 )
	my.cols.right <- brewer.pal ( max ( unlist ( seq.structure ) ) + 1 , "Paired" )
	my.cols.left <- brewer.pal ( max ( unlist ( seq.structure ) ) + 1 , "Set3" )
	plot ( NA , bty = "n" , xlim = c ( -1 , 1 ) , ylim = c ( 0 , 12 ) , xaxt = "n" , yaxt = "n" , ylab = "" , xlab = "" )

	### right side
	for ( row in seq_len ( nrow ( seq.structure$right.seq ) ) ) {
		my.recs <- unique ( seq.structure$right.seq[row,] )
		recode.my.recs <-  c ( unique ( seq.structure$right.seq[row,] ) , max ( unlist ( seq.structure$right.seq ) ) + 1 ) + 1
		for ( i in seq_along ( my.recs ) ) {
			polygon ( x = c ( scaled.rec.points$right [ recode.my.recs [ i ] ] , scaled.rec.points$right [ recode.my.recs [ i ] ] , scaled.rec.points$right [ recode.my.recs [ i + 1 ] ] , scaled.rec.points$right [ recode.my.recs [ i + 1 ] ] ) , y = c ( row , row - 1 , row - 1 , row  ) , col = my.cols.right [ my.recs [ i ] + 1 ] , lty = 0 )
		}
	}


	### left side
	for ( row in seq_len ( nrow ( seq.structure$left.seq ) ) ) {
		my.recs <- unique ( seq.structure$left.seq[row,] )
		recode.my.recs <-  c ( unique ( seq.structure$left.seq[row,] ) , max ( unlist ( seq.structure$left.seq ) ) + 1 ) + 1
		for ( i in seq_along ( my.recs ) ) {
			polygon ( x = c ( scaled.rec.points$left [ recode.my.recs [ i ] ] , scaled.rec.points$left [ recode.my.recs [ i ] ] , scaled.rec.points$left [ recode.my.recs [ i + 1 ] ] , scaled.rec.points$left [ recode.my.recs [ i + 1 ] ] ) , y = c ( row , row - 1 , row - 1 , row  ) , col = my.cols.left [ my.recs [ i ] + 1 ] , lty = 0 )
		}
	}
	abline ( v = 0 )
}

}

setwd ( "~/Documents/Academics/StandingSweeps/" )
my.s <- 1:20 / 10000
my.sims <- list ()
term.in.exp <- list()
for ( i in 1:length ( my.s ) ) {
	term.in.exp [[ i ]] <- numeric(20)
	for ( j in 1 : 20 ) {	
		my.sims <- StructuredCoalescentSweep ( N = 10000 , s = my.s [ i ] , dominance = FALSE , f = 1/20000 , reps = 1000 , n.tips = 2 , r = 10^-8 , sim.distance = 0.02 , interval.width = 1000 , no.sweep = F , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = FALSE , display.rep.count = TRUE ,   standing.haps = FALSE , time.factor = 1 )
		term.in.exp [[ i ]] [ j ] <- GetTermInExp ( my.sims )
		rm ( my.sims )
		gc ( )
		message ( i )
		message ( j )
	}
}
save ( term.in.exp , "Sims/term.in.exp.Robj")