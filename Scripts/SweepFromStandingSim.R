##install.packages("randtoolbox")
##install.packages("ape")
library("randtoolbox")
library("ape")
turn.on.recovers=FALSE

StructuredCoalescentSweep <- function ( N , s , f , reps , n.tips , r , sim.distance , interval.width , no.sweep = FALSE , constant.freq = FALSE ) {
	options ( error = recover )
	#recover()
	
	
	
	
	if ( constant.freq == FALSE ) {
	
		temp <- SweepFromStandingSim ( N = N , s = s , f = f , time.factor = time.factor , reps = reps , no.sweep = no.sweep )
		frequencies <- temp [[ 1 ]]
	
		if ( no.sweep == FALSE ) {	
			sweep.start <- rep ( temp[[2]] , nrow ( frequencies ) )
			sweep.start.forward <- ncol ( frequencies ) - sweep.start[1] 
			# if ( nrow ( frequencies ) > 1 ) {
			fixation.time <- apply ( frequencies [ , sweep.start.forward : ncol ( frequencies ) ] , 1 , which.max ) + sweep.start.forward - 1
			zeros <- apply ( frequencies [ , 1 : sweep.start.forward ] %% 1 == 0 , 1 , which )
			entry <- numeric()
			for ( i in 1 : length ( zeros ) ){
				if ( length ( zeros [[ i ]] ) != 0 ) {
					entry [ i ] <- tail ( zeros [[ i ]] , 1 )
				} else {
					entry [ i ] <- 1
				}
			}
			transit.time <- fixation.time - entry
			new.freqs <- matrix ( 0 , nrow = reps , ncol = max ( transit.time + 1 ) )
			for ( i in 1 : nrow ( frequencies ) ) {
				new.freqs [ i , 1 : ( transit.time [ i ] + 1 ) ] <- frequencies [ i , fixation.time [ i ] : entry [ i ] ]
			}
	
		} else if ( no.sweep == TRUE ){
			#recover()
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
	}

	
	#recover()
	temp <- RecombinationEvents ( trees = trees , coal.times = coal.times , r = r , sim.distance = sim.distance , n.tips = n.tips )
	trees <- temp [[ 1 ]]
	T.total <- temp [[ 2 ]]

	#recover()
	trees <- BuildOnOffHaps ( trees = trees , freqs = new.freqs , sim.distance = sim.distance , r = r , n.tips = n.tips , f = f , fixation.time = fixation.time )

	hap.dist <- HapCountDistribution ( input = trees , r = r , sim.distance = sim.distance , interval.width = interval.width , f = f , s = s , reps = reps , N = N , n.tips = n.tips )
	
	return ( list ( coal.times = coal.times , new.freqs = new.freqs , mean.coalescence.times = mean.coalescence.times , sd.coalescence.times = sd.coalescence.times , trees = trees , hap.dist = hap.dist , fixation.time = fixation.time , T.total = T.total ) )
}

SweepFromStandingSim <- function ( N , s , f , time.factor ,  reps , no.sweep ) {
	
	delta.T <- 1 / ( 2 * N )
	sweep.freq.matrix <- list ( rep ( f , reps ) )
	neutral.freq.matrix <- list ( rep ( f , reps ) )
	not.all.sweeps.fixed <- TRUE
	not.all.neutral.fixed <- TRUE
	#recover()
	i = 1
	while ( not.all.sweeps.fixed  | not.all.neutral.fixed ) {
		if ( not.all.sweeps.fixed ) {
			update <- rep ( 0 , reps )
			sweep.not.fixed <- sweep.freq.matrix [[ i ]] %% 1 != 0
			sweep.fixed <- sweep.freq.matrix [[ i ]] %% 1 == 0
			mu.S <- 2 * N * s * sweep.freq.matrix [[ i ]] [ sweep.not.fixed ] * ( 1 - sweep.freq.matrix [[ i ]] [ sweep.not.fixed ] )
			sel <- mu.S * delta.T
			sweep.drift.mag <- sqrt ( sweep.freq.matrix [[ i ]] [ sweep.not.fixed ] * ( 1 - sweep.freq.matrix [[ i ]] [ sweep.not.fixed ] ) * delta.T)
			plus.minus <- sample ( c ( 0 , 1 ) , sum ( sweep.not.fixed ) , replace = TRUE )
			drift.sweep <- ifelse ( plus.minus == 1 , sweep.drift.mag , -1 * sweep.drift.mag )
			update [ sweep.not.fixed ] <- sel + drift.sweep			
			sweep.freq.matrix [[ i + 1 ]] <- sweep.freq.matrix [[ i ]] + update
			sweep.fixed.one <- sweep.freq.matrix [[ i + 1 ]] > ( 1 - ( 1 / ( 2 * N ) ) )
			sweep.freq.matrix [[ i + 1 ]] [ sweep.fixed.one ] <- 1
			sweep.fixed.zero <- sweep.freq.matrix [[ i + 1 ]] < 1 / ( 2 * N )
			sweep.freq.matrix [[ i + 1 ]] [ sweep.fixed.zero ] <- 0
			not.all.sweeps.fixed <- any ( sweep.freq.matrix [[ i + 1 ]] %% 1 != 0 )

		
		}
		if ( not.all.neutral.fixed ) {	
			update <- rep ( 0 , reps )
			neutral.not.fixed <- neutral.freq.matrix [[ i ]] %% 1 != 0
			neutral.fixed <- neutral.freq.matrix [[ i ]] %% 1 == 0
			#neutral.drift.mag <- sqrt ( neutral.freq.matrix [[ i ]] [ neutral.not.fixed ] * ( 1 - neutral.freq.matrix [[ i ]] [ neutral.not.fixed ] ) * delta.T )
			#plus.minus <- sample ( c ( 0 , 1 ) , sum ( neutral.not.fixed ) , replace = TRUE )	
			#drift.neutral <- ifelse ( plus.minus == 1 , neutral.drift.mag , -1 * neutral.drift.mag )
			drift.neutral <- rnorm ( sum ( neutral.not.fixed ) , - neutral.freq.matrix [[ i ]] [ neutral.not.fixed ] * 1 / ( 2 * N ) , sd = sqrt ( neutral.freq.matrix [[ i ]] [ neutral.not.fixed ] * ( 1 - neutral.freq.matrix [[ i ]] [ neutral.not.fixed ] ) * 1 / ( 2 * N ) ) )
			update [ neutral.not.fixed ] <- drift.neutral
			neutral.freq.matrix [[ i + 1 ]] <- neutral.freq.matrix [[ i ]] + update
			neutral.fixed.one <- neutral.freq.matrix [[ i + 1 ]] > ( 1 - ( 1 / ( 2 * N ) ) )
			neutral.freq.matrix [[ i + 1 ]] [ neutral.fixed.one ] <- 1
			neutral.fixed.zero <- neutral.freq.matrix [[ i + 1 ]] < 1 / ( 2 * N )
			neutral.freq.matrix [[ i + 1 ]] [ neutral.fixed.zero ] <- 0	
			not.all.neutral.fixed <- any ( neutral.freq.matrix [[ i ]] %% 1 != 0 )
		
		}
		if ( i %% 5000 == 0 ) {
				lineages.remaining <- sum ( neutral.freq.matrix [[ i + 1 ]] %% 1 != 0 )
				my.freq <- max ( neutral.freq.matrix [[ i + 1 ]] [ neutral.freq.matrix [[ i + 1 ]] < 1 ] )
				cat ( "p = " , my.freq , ",  " , sep = "" )
				cat ( lineages.remaining , "not fixed \n")
		}		
		if ( i == 16 * N ){
			break
		}
		i = i + 1
	}
	sweep.freq.matrix <- matrix ( unlist ( sweep.freq.matrix ) , nrow = reps )
	neutral.freq.matrix <- matrix ( unlist ( neutral.freq.matrix ) , nrow = reps )
if(turn.on.recovers)	recover()
	# if ( constant.freq == FALSE ) {
		# # if ( reps == 1 ) {
			# # freq.trajectories <- c ( neutral.freq.matrix [ length ( neutral.freq.matrix ) : 2 ] , sweep.freq.matrix [ 1 : length ( sweep.freq.matrix ) ] )
			# # #plot ( freq.trajectories , type = "l" )
			# # freq.trajectories <- matrix ( freq.trajectories , nrow = 1 )
			# # if ( freq.trajectories [ , ncol ( freq.trajectories ) ] == 1 ) {
				# # conditional.freq.trajectories <- freq.trajectories
				# # #generations <- seq ( 1 , ncol ( conditional.freq.trajectories ) , time.factor )
				# # #conditional.freq.trajectories <- conditional.freq.trajectories [ , generations ]
			# # } else {
				# # cat ( "Allele lost from population.\n")
				# # return ( )
			# # }
		# # } else {
		
	if ( no.sweep == FALSE ) {
		freq.trajectories <- cbind ( neutral.freq.matrix [ , ncol ( neutral.freq.matrix ) : 2 ] , sweep.freq.matrix [ , 1 : ncol ( sweep.freq.matrix ) ] )
	} else {
		freq.trajectories <- neutral.freq.matrix [ , 1 : ncol ( neutral.freq.matrix ) ]
		return ( list ( freq.trajectories , 0 ) )
	}
		# }
	# } else {
		# freq.trajectories <- sweep.freq.matrix [ , 1 : ncol ( sweep.freq.matrix ) ]
	# }
	#recover()
	keep.these <- freq.trajectories [ , ncol ( freq.trajectories ) ] == 1		
	conditional.freq.trajectories <- freq.trajectories [ keep.these , ]
	sweep.start <- ncol ( sweep.freq.matrix ) #/ time.factor
	return ( list ( conditional.freq.trajectories , sweep.start ) )	
}


BuildTrees <- function ( coal.times ){
	#recover()
	library ( ape )
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




# BuildTrees <- function ( coal.times , n.tips ){
	# #recover()
	# library(ape)
	# trees <- list ( )
	# for ( j in 1 : nrow ( coal.times ) ) {
		# edge <- matrix ( 0 , nrow = 2 * n.tips - 2 , ncol = 2 )
		# edge.length <- numeric ( 2 * n.tips - 2 )
		# edge [ 1 : n.tips , 2 ] <- 1 : n.tips
		# nodes <- ( 2 * n.tips - 1 ) : ( n.tips + 1 )
		# node.depth <- numeric ( 2 * n.tips - 1 )
		# Nnode <- n.tips - 1
		# tip.label <- character ( n.tips )
		# for ( l in 1 : length ( tip.label ) ){
			# tip.label [ l ] <- paste ( "t" , l , sep = "")
		# }
		# k = 1
		# for ( i in nodes ) {
			# extant.lineages <- edge [ edge [ , 2] != 0 & edge [ , 1 ] == 0 , 2 ]
			# coalescing.lineages <- sort ( sample ( extant.lineages , 2 , replace = FALSE ) )
			# coal.index <- which ( edge [ , 2 ] %in% coalescing.lineages )
			# edge [ coal.index , 1 ] <- i
			# if ( i != tail ( nodes , 1 ) ) {
				# edge [ i - 1 , 2 ] <- i
			# }
			# node.depth [ i ] <- coal.times [ j , k ] 
			# edge.length [ coal.index ] <- coal.times [ j , k ] - node.depth [ coalescing.lineages ]
			# k = k + 1
		# }
		# a.tree <- list ( edge = edge , edge.length = edge.length , tip.label = tip.label , Nnode = Nnode )
		# class ( a.tree ) <- "phylo"
		# my.tree <- list ( tree = a.tree , node.depth = node.depth )
		
		# trees [[ j ]] <- my.tree
	# }
	
	# ## note; the frequency path gets added on as trees [[ j ]] [[ 3 ]] in the top level function call immediately after returning from this function; bad writing that I should fix some day.
	# return ( trees )
	
# }

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

	}
		
	return ( list ( trees, T.total ) )

}

BuildOnOffHaps <- function ( trees , freqs , r , sim.distance , n.tips , f , fixation.time ) {
	
	sim.distance.bp <- sim.distance / r
	#recover()
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
		if ( length ( to.remove ) != 0 ){			
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
					my.freq <- rev ( trees [[ j ]] [[ 3 ]] ) [ this.event$rec.depth ]
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
		if ( length ( to.remove ) != 0 ) {
			left.sequence <- left.sequence.temp [ , -to.remove ]
			left.sequence <- MakeHapsPretty ( left.sequence )
			rec.left.off.background <- rec.left [ - ( to.remove - 1 ) , ]
		} else {
			left.sequence <- left.sequence.temp
			left.sequence <- MakeHapsPretty ( left.sequence )
			rec.left.off.background <- rec.left
		}

		
		trees [[ j ]] [[ "sequence.structure" ]] <- list ( right.seq = right.sequence , left.seq = left.sequence )
		trees [[ j ]] [[ "rec.events.off.background" ]] <- list ( rec.right.off.background = rec.right.off.background , rec.left.off.background = rec.left.off.background )
	}
	
	return ( trees )
}


HapCountDistribution <- function ( input , r , sim.distance , interval.width , f , s , reps , N , n.tips ) {
	
	sim.distance.bp <- sim.distance / r 
	intervals <- seq ( 0 , sim.distance.bp , interval.width )
	
	# number of rows in "sequence" matrix = number of samples

	
	n.haps.right <- n.haps.left <- matrix ( nrow = length ( input ) , ncol = length ( intervals ) )
	#recover()
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
				my.rec.events <- input [[ j ]] $ rec.events.off.background
				
				# right side
				if ( sum ( my.rec.events$rec.right.off.background$sequence.location < k ) != 0 ) {
					last.rec.event <- max ( which ( my.rec.events$rec.right.off.background$sequence.location < k ) )
					n.haps.right [ j , i ] <-  length ( unique ( my.seqs$right.seq [ , last.rec.event + 1 ] ) )
				} else {
					n.haps.right [ j , i ] <- 1
				}
				
				# left.side
				if ( sum ( my.rec.events$rec.left.off.background$sequence.location < k ) != 0 ) {
					last.rec.event <- max ( which ( my.rec.events$rec.left.off.background$sequence.location < k ) )
					n.haps.left [ j , i ] <-  length ( unique ( my.seqs$left.seq [ , last.rec.event + 1 ] ) )
				} else {
					n.haps.left [ j , i ] <- 1
				}
			
			}	
		}	
	}

	#recover()
	n.haps <- rbind ( n.haps.right , n.haps.left )
	
	hap.counts.by.interval <- apply ( n.haps , 2 , function ( x ) table ( factor ( x , 1 : n.tips ) ) )
	hap.count.freqs.by.interval <- apply ( hap.counts.by.interval , 2 , function ( x ) x / nrow ( n.haps ) )
	
	cum.probs <- rbind ( 0 , apply ( hap.count.freqs.by.interval , 2 , cumsum ) )
#	par ( mfrow = c ( 2 , 1 ) )
	#matplot ( t ( cum.probs ) , type = "l" , lty = 1 , lwd = 0.7 , col = "black" , ylab = "Cumulative Probability" , xlab = "kb" , main = paste ( n.tips , "Lineages in a Sweep from f =" , f , "at s =" , s , "," , reps , "Reps" ) , bty = "n")
	#recover()
	
	
	ewens.dist.matrix <- matrix ( nrow = n.tips , ncol = length ( intervals ) )
	
	stirling.numbers <- StirlingNumbers ( n = n.tips ) [ n.tips , ]
	for ( i in 1 : length ( intervals ) ) {
		
		if ( i == 1 & intervals [ 1 ] == 0 ) {
			
			ewens.dist.matrix [ , i ] <- c ( 1 , rep ( 0 , n.tips - 1 ) )
			
		} else { 
		
			ewens.dist.matrix [ , i ] <- EwensDist ( n = n.tips , N = N , r = r , distance = intervals [ i ] , f = f , stirling.numbers = stirling.numbers )
			
		}
		
	}
	#recover()
	ewens.cum.probs <-  apply ( ewens.dist.matrix , 2 , cumsum )
	matplot ( 
		t ( ewens.cum.probs ) , 
		type = "n" , 
		lty = 1 , 
		lwd = 0.7 , 
		col = "black" , 
		ylab = "Cumulative Probability" , 
		xlab = "kb" , 
		#main = paste ( n.tips , "Lineages in a Sweep from f =" , f , "at s =" , s , "," , reps , "Reps" ) , 
		bty = "n"
	)
	
	
	#recover()
	col.vect <- rainbow ( n.tips , s = 0.8  , v = 1 , start = 1/40 , end = 4/6  )
	for ( i in 1 : ( nrow ( cum.probs ) - 1 ) ) {
			#i = i + 1
			X.ax <- 1:ncol ( cum.probs ) #which ( cum.probs [ i , ] != cum.probs [ i + 1 , ] )
			Y.ax1 <- cum.probs [ i , X.ax ]
			Y.ax2 <- cum.probs [ i + 1 , X.ax ]
			polygon ( x = c ( X.ax , rev ( X.ax ) ) , y = c ( Y.ax1 , rev ( Y.ax2 ) ) , lty = 0 , col = col.vect [ i ] )					
	}
	ewens.cum.probs <- ewens.cum.probs [ - nrow ( ewens.cum.probs ) , ]
	apply ( ewens.cum.probs , 1 , function ( x ) lines ( x , lty = 1 , lwd = 0.8 ) )
	
	#recover()
	
	expected.num.haps <- colSums ( apply ( hap.counts.by.interval , 2 , function ( x ) x * 1 : n.tips ) / (2 * reps ) , 2 )
		
	#plot ( expected.num.haps , type = "l" , lty = 1 , lwd = 1.5 , xlab = "kb" , ylab = "Expected Number of Haplotypes" , bty = "n")
	
	return ( list ( hap.count.freq.by.interval = hap.count.freqs.by.interval , n.haps = n.haps ) )
	
}

StirlingNumbers <- function ( n ) {
	
	library ( randtoolbox )
	second.kind <- lapply ( 1 : n , stirling )
	second.kind.matrix <- matrix ( nrow = n , ncol = n )
	for ( i in 1 : n ) {
		if ( i < n ) {
			second.kind.matrix [ i , ] <- c ( second.kind [[ i ]] [ -1 ], rep ( 0 , n - length ( second.kind [[ i ]] ) + 1 ) )
		} else if ( i == n ) {
			second.kind.matrix [ i , ] <- second.kind [[ i ]] [ -1 ]
		}
	}
	#recover()
	first.kind.matrix <- abs ( solve ( second.kind.matrix ) )
	first.kind.matrix [ first.kind.matrix < 0.99 ] <- 0
		
	return ( first.kind.matrix )
	
}

EwensDist <- function ( n , N , r , distance , f , stirling.numbers ) {
	#recover()	
	param <- 4 * N * r * distance * f * ( 1- f )
	denom  <- prod ( param + 0 : ( n - 1 ) )
	#stirling.numbers <- StirlingNumbers ( n ) [ n , ]
	ewens.dist <- param^(1:n) * stirling.numbers / denom
	return ( ewens.dist ) 

}


MakeHapsPretty <- function ( seqs ) {
	
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


if(FALSE){

temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.5 , f = 0.02 , reps = 400 , n.tips = 12 , r = 10^-8 , sim.distance = 0.02 , interval.width = 1000 , no.sweep = TRUE , constant.freq = FALSE )






# Let's think about inference
my.seqs <- temp[["trees"]][[1]][["sequence.structure"]]$right.seq

InferenceFunction <- function ( seqs ) {
	
if(turn.on.recovers)	recover()
	hap.partitions <- apply ( seqs , 2 , function ( x ) table ( factor ( x , levels = 0 : ( nrow ( seqs ) - 1 ) ) ) )
	tree <- BuildTrees ( 1 : ( nrow ( seqs ) - 1 ) )
	
	
}


InferenceFunction ( seqs = my.seqs )


i = 1
par(mfrow=c(2,1))
plot ( temp$trees[[i]]$freqs , type = "l" , xlim = c ( length ( temp$trees[[i]][[3]] ) - max ( temp$trees[[i]][[2]] ) , length ( temp$trees[[i]][[3]] ) ) )
plot ( temp$trees[[i]][[1]] , x.lim = c ( 0 , max ( temp$trees[[i]][[2]] ) ) )
temp$trees[[i]][[5]]; i = i + 1

}