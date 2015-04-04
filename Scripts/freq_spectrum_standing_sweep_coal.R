setwd ( "~/Documents/Academics/StandingSweeps/" )
source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
source('~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R', chdir = TRUE)


####freq. spectrum w. no rec. during sweeps.

# expected.freq.times.standing<-function(nsam,N,r,distance,f,my.StirlingNumbers=NULL){
	# #recover()
	# #my.StirlingNumbers<-StirlingNumbers(n)
	# ESF.prob.k <- EwensDist( n=nsam , N =N, r=r , distance=1 , f=f)
	# ESF.cond.prob.k<-EwensCondDist( n=nsam , N =N, r=r , distance=1 , f=f) # ,stirling.numbers=my.StirlingNumbers)    ### is of form [n,k]
	# if ( all ( my.StirlingNumbers == 0 ) ) {
		# my.StirlingNumbers<-StirlingNumbers(nsam)    ##Usigned Stirling numbers of 1st kind. ma
	# }
	# expected.t.l<-rep(NA,nsam-1)
	# p_l_given_k <- array ( 0 , dim = c ( nsam , nsam , nsam ) )
	# freq.specs <- matrix ( 0 , nrow = nsam , ncol = nsam )
	# for ( i in 2 : nsam ) {
		# freq.specs [ 1 : ( i - 1 ) , i ] <- ( 1 / ( 1 : ( i - 1 ) ) ) / ( sum ( 1 / ( 1 : ( i - 1  )  ) ) )
	# }
	# freq.specs <- t ( freq.specs )
	# terms.in.sum <- array ( 0 , dim = c ( nsam , nsam , nsam ) )
	# for(l in 1:(nsam-1)){
	# #	recover()
	# #	terms.in.sum<-rep(0,nsam)
		# terms.given.j <- matrix ( 0 , ncol = nsam , nrow = nsam )
		# for(k in 2 : nsam ) {
			# ##runs from 2 otherwise there are no polymorphism

			# for(j in 1:(k-1)){
				# stirling.bit <- my.StirlingNumbers[l,j] * my.StirlingNumbers[nsam-l,k-j]  / my.StirlingNumbers[nsam,k]
				# p_l_given_k [ k , j , l ] <- stirling.bit*choose(nsam,l)/choose(k,j)
				# terms.given.j [ j , k ] <- p_l_given_k [ k , j , l ]*freq.specs[ k , j ]
				# terms.in.sum [ k , j, l  ] <- ESF.cond.prob.k [ nsam , k ] * p_l_given_k [ k , j , l ] * freq.specs [ k , j ]
		
			# }
		# }
	# }
	# anc.freq.spec <- colSums ( terms.in.sum , dims = 2 )
	# new.freq.spec <- 1 / ( 1 : ( nsam - 1 ) ) /  sum  ( 1 / ( 1 : ( nsam - 1 ) ) )
	# times <- sapply ( 1 : ( nsam - 1 ) , function ( x ) sum ( 1 / ( 1:x ) ) )
	# new <- ESF.prob.k [ nsam , 1 ] * 4 * f * sum ( 1 / ( 1 : 1 : ( nsam - 1 ) ) )
	# old <- ESF.prob.k [ nsam , 2:(nsam) ] * times
	# p.new.seg <- new / ( new + sum ( old ) )
	# my.freq.spec <- ESF.prob.k [ nsam , 1 ] * new.freq.spec + ( 1 -ESF.prob.k [ nsam , 1 ] ) * anc.freq.spec [ - length ( anc.freq.spec ) ]
	# adj.freq.spec <- p.new.seg * new.freq.spec + ( 1 - p.new.seg ) * anc.freq.spec [ - length ( anc.freq.spec ) ]
	# return ( list ( my.freq.spec , adj.freq.spec ) )
# }




### doing new mutations the correct way

expected.freq.times.standing<-function(nsam,N,r,distance,f,my.StirlingNumbers=NULL){
	#recover()
	#my.StirlingNumbers<-StirlingNumbers(n)
	ESF.prob.k <- EwensDist( n=nsam , N =N, r=r , distance=1 , f=f)
	ESF.cond.prob.k<-EwensCondDist( n=nsam , N =N, r=r , distance=1 , f=f) # ,stirling.numbers=my.StirlingNumbers)    ### is of form [n,k]
	if ( all ( my.StirlingNumbers == 0 ) ) {
		my.StirlingNumbers<-StirlingNumbers(nsam)    ##Usigned Stirling numbers of 1st kind. ma
	}
	expected.t.l<-rep(NA,nsam-1)
	p_l_given_k <- array ( 0 , dim = c ( nsam , nsam , nsam ) )
	freq.specs.anc <- matrix ( 0 , nrow = nsam , ncol = nsam )
	for ( i in 2 : nsam ) {
		freq.specs.anc [ 1 : ( i - 1 ) , i ] <- ( 1 / ( 1 : ( i - 1 ) ) )
	}
	freq.specs.anc <- t ( freq.specs.anc )	
	terms.in.sum <- array ( 0 , dim = c ( nsam , nsam , nsam ) )
	for(l in 1:(nsam-1)){
	#	recover()
	#	terms.in.sum<-rep(0,nsam)
		terms.given.j <- matrix ( 0 , ncol = nsam , nrow = nsam )
		for(k in 1 : nsam ) {
			##runs from 2 otherwise there are no polymorphism

			for(j in 1:(k-1)){
				stirling.bit <- my.StirlingNumbers[l,j] * my.StirlingNumbers[nsam-l,k-j]  / my.StirlingNumbers[nsam,k]
				
				p_l_given_k [ k , j , l ] <- ifelse ( 
																	k > 1 ,
																	stirling.bit*choose(nsam,l)/choose(k,j) ,
																	0
																)
				#terms.given.j [ j , k ] <- p_l_given_k [ k , j , l ]*freq.specs.anc[ k , j ]
				new.freq <- ifelse ( i == nsam & k == 1 , f / l , 0 )
				terms.in.sum [ k , j, l ] <- ESF.prob.k [ nsam , k ] * ( p_l_given_k [ k , j , l ] * freq.specs.anc [ k , j ] + new.freq )
		
			}
		}
	}
	temp <- colSums ( terms.in.sum , dims = 2 )
	freq.spec <- temp [ 1 : ( length ( temp ) - 1 ) ] / sum ( temp )
	# new.freq.spec <- 1 / ( 1 : ( nsam - 1 ) ) /  sum  ( 1 / ( 1 : ( nsam - 1 ) ) )
	# times <- sapply ( 1 : ( nsam - 1 ) , function ( x ) sum ( 1 / ( 1:x ) ) )
	# new <- ESF.prob.k [ nsam , 1 ] * 4 * f * sum ( 1 / ( 1 : 1 : ( nsam - 1 ) ) )
	# old <- ESF.prob.k [ nsam , 2:(nsam) ] * times
	# p.new.seg <- new / ( new + sum ( old ) )
	# my.freq.spec <- ESF.prob.k [ nsam , 1 ] * new.freq.spec + ( 1 -ESF.prob.k [ nsam , 1 ] ) * anc.freq.spec [ - length ( anc.freq.spec ) ]
	# adj.freq.spec <- p.new.seg * new.freq.spec + ( 1 - p.new.seg ) * anc.freq.spec [ - length ( anc.freq.spec ) ]
	return ( list ( 0 , freq.spec ) )
}










expected.freq.times.standing.anc.only<-function(nsam,N,r,distance,f,my.StirlingNumbers=NULL){
	#recover()
	#my.StirlingNumbers<-StirlingNumbers(n)
	ESF.prob.k <- EwensDist( n=nsam , N =N, r=r , distance=1 , f=f)
	ESF.cond.prob.k<-EwensCondDist( n=nsam , N =N, r=r , distance=1 , f=f) # ,stirling.numbers=my.StirlingNumbers)    ### is of form [n,k]
	if ( all ( my.StirlingNumbers == 0 ) ) {
		my.StirlingNumbers<-StirlingNumbers(nsam)    ##Usigned Stirling numbers of 1st kind. ma
	}
	expected.t.l<-rep(NA,nsam-1)
	p_l_given_k <- array ( 0 , dim = c ( nsam , nsam , nsam ) )
	freq.specs <- matrix ( 0 , nrow = nsam , ncol = nsam )
	for ( i in 2 : nsam ) {
		freq.specs [ 1 : ( i - 1 ) , i ] <- ( 1 / ( 1 : ( i - 1 ) ) ) / ( sum ( 1 / ( 1 : ( i - 1  )  ) ) )
	}
	freq.specs <- t ( freq.specs )
	terms.in.sum <- array ( 0 , dim = c ( nsam , nsam , nsam ) )
	for(l in 1:(nsam-1)){
	#	recover()
	#	terms.in.sum<-rep(0,nsam)
		terms.given.j <- matrix ( 0 , ncol = nsam , nrow = nsam )
		for(k in 2 : nsam ) {
			##runs from 2 otherwise there are no polymorphism

			for(j in 1:(k-1)){
				stirling.bit <- my.StirlingNumbers[l,j] * my.StirlingNumbers[nsam-l,k-j]  / my.StirlingNumbers[nsam,k]
				p_l_given_k [ k , j , l ] <- stirling.bit*choose(nsam,l)/choose(k,j)
				terms.given.j [ j , k ] <- p_l_given_k [ k , j , l ]*freq.specs[ k , j ]
				terms.in.sum [ k , j, l  ] <- ESF.cond.prob.k [ nsam , k ] * p_l_given_k [ k , j , l ] * freq.specs [ k , j ]
		
			}
		}
	}
	anc.freq.spec <- colSums ( terms.in.sum , dims = 2 )
	# new.freq.spec <- 1 / ( 1 : ( nsam - 1 ) ) /  sum  ( 1 / ( 1 : ( nsam - 1 ) ) )
	# my.freq.spec <- ESF.prob.k [ nsam , 1 ] * new.freq.spec + ( 1 -ESF.prob.k [ nsam , 1 ] ) * anc.freq.spec [ - length ( anc.freq.spec ) ]
	return ( anc.freq.spec )
}





##### a component function to make the code cleaner
plbarjkn <- function ( l , j , k , n , my.StirlingNumbers ) {
	#recover()
	####################
	#### special cases ####
	####################
	
	
	if ( l == 0 & j == 0 ) 	return ( 1 )
	if ( j == k & n == l )	return ( 1 )
	if ( k == n & j == l )	return ( 1 )
	if ( j == k & n != l )	return ( 0 )
	if ( l > 0 & j == 0 ) 	return ( 0 )
	if ( l == 0 & j > 0 ) 	return ( 0 )
	if ( l > j + n - k )	return ( 0 )
	if ( k == 1 ) {
		if ( n == l & j == 1 ) return ( 1 )
		if ( l != n ) return ( 0 )
		else stop ()
	}


	####################
	#### standard case ####
	####################
	stirling.bit <- my.StirlingNumbers[ l , j ] * my.StirlingNumbers[n - l , k - j ]  / my.StirlingNumbers[ n , k ]
	
	stirling.bit * choose ( n , l ) / choose ( k , j )
	
}


####freq. spectrum with rec. during sweeps.

# expected.freq.times.standing.w.sweep <- function ( nsam , N , r , f , s , my.StirlingNumbers = 0 ) {
	# #recover()
	# if ( all ( my.StirlingNumbers == 0 ) ) {
		# my.StirlingNumbers<-StirlingNumbers(nsam)    ##Usigned Stirling numbers of 1st kind. ma
	# }
	# ESF.prob.k <- EwensDist( n=nsam , N =N, r=r , distance=1 , f=f)
	# ESF.prob.k <- rbind ( c ( 1 , rep ( 0 , nsam ) ) , cbind ( rep ( 0 , nsam ) , ESF.prob.k ) )
	# ESF.condprob.k<-EwensCondDist( n=nsam , N =N, r=r , distance=1 , f=f)
	# ESF.condprob.k <- rbind ( c ( rep ( 0 , nsam + 1 ) ) , cbind ( rep ( 0 , nsam ) , ESF.condprob.k ) )
	# T_f <- log ( (2*N -1 ) * ( 1 - f ) / f ) / s
	# my.logistic <- function ( x ) 1 / (2 * N  ) * exp(s * x ) / ( 1 + 1 / (2 * N  ) * ( exp(s * x )  - 1 ) )
	# T_sf <- integrate ( my.logistic , 0 , T_f )$value
	# P_NR <- exp ( - r * T_sf )
	# expected.t.l<-rep(NA,nsam-1)
	# p_l_given_k <- array ( 0 , dim = c ( nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 ) )
	# H <- array ( 0 , dim = c ( nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 ) )
	# cond.freq.specs <- array ( 0 , dim = c ( nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 ) )
	# freq.specs <- matrix ( 0 , nrow = nsam , ncol = nsam )
	# for ( i in 2 : nsam ) {
		# freq.specs [ 1 : ( i - 1 ) , i ] <- ( 1 / ( 1 : ( i - 1 ) ) ) / ( sum ( 1 / ( 1 : ( i - 1  )  ) ) )
	# }
	# freq.specs <- t ( freq.specs )
	# my.freq.specs <- array ( 0 , dim = c ( nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 ) )
	# z <- 1
	# my.list <- list ()
	# for ( l in 1 : ( nsam - 1 ) ) {
	# #	recover()
	# #	terms.in.sum<-rep(0,nsam)
		# for ( i in 0 : nsam ) {
			
			# for ( k in min ( 1 , i ) : i ) {
				# terms.given.j <- matrix ( 0 , ncol = nsam , nrow = nsam )
				# for ( j in 1 : min ( k + nsam - i , l ) ) {
					# if ( max ( 0 , ( j - k ) ) > min ( l , nsam - i , j ) ) next
					# g.sec <- seq ( max ( 0 , ( j - k ) ) ,  min ( l , nsam - i , j ) , 1 )
					# for ( g in g.sec ) {
						# # if the number of lineages sitting under the beneficial mutatation is smaller than the number we need to get to l, the prob is zero
						# if ( i < ( l - g ) ) next
						
						
						# p_l_given_k [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] <- plbarjkn ( l - g , j - g , k , i , my.StirlingNumbers )
						

						# H [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] <- choose ( nsam - i , g ) * choose ( k , j - g ) / choose ( k + nsam - i , j )
						
						# cond.freq.specs [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] <- freq.specs [ k + nsam - i , j ] * H [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] * p_l_given_k [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ]
						
						# if ( g > 0 | i != nsam ) {
							# my.freq.specs [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] <- dbinom ( i , nsam , P_NR ) * ESF.condprob.k [ i + 1 , k + 1 ] * cond.freq.specs [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ]
						
						# } else{
							# my.list [[ z ]] <- c ( k + 1 , j + 1 , g + 1 , i + 1 , l + 1 )
							# my.freq.specs [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] <- dbinom ( i , nsam , P_NR ) * ESF.condprob.k [ i  + 1 , k + 1 ] * cond.freq.specs [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ]
							# #my.freq.specs [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] <- dbinom ( i , nsam , P_NR ) * ESF.condprob.k [ i  + 1 , k + 1 ] * cond.freq.specs [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ]
							# z <- z + 1
						# }
					# }
				# }
			# }
		# }
	# #	expected.t.l[l]<-sum(terms.in.sum)
	# }
	# singleton.reweighted.times <- 1/(1:(nsam-1)) + c ( nsam * T_f / ( 2 * N * f ) , rep ( 0 , nsam - 2 ) )
	# singleton.reweighted.freq.spec <- singleton.reweighted.times / sum ( singleton.reweighted.times )
	# this.freq.spec <- numeric ( nsam - 1 )
	# for ( l in 2 : ( dim ( my.freq.specs ) [ 5 ] - 1 ) ) {
		# this.freq.spec [ l - 1 ] <- sum ( my.freq.specs [ , , , , l ] ) #+ dbinom ( nsam , nsam , P_NR ) * ESF.prob.k [ nsam +1 , 2 ] * singleton.reweighted.freq.spec [  l - 1 ]
	# }
	# my.branch.sums <- cumsum ( singleton.reweighted.times )
	# num <- ESF.prob.k [ nsam + 1 , 2 ] * dbinom ( 10 , 10 , P_NR ) * f * my.branch.sums [ nsam - 1 ]
	# other <- numeric()
	# P_R <- 1 - P_NR
	# t <- 1
	# for ( j in 2 : ( nsam - 1 ) ) {
		# for ( k in 0 : j ) {
			# for ( i in 0 : ( j - k ) ) {
				# if ( i < 0 | (k + i) <= 1 ) { next }
				
				# other [ t ] <- dbinom ( i , nsam , P_R ) * ESF.prob.k [ nsam - i + 1  , k + 1 ] * my.branch.sums [ j ]
				# t <- t + 1
			# }
		# }
	# }
	# p_new <- num / ( num + sum ( other ) )
	# final.freq.spec <- p_new * singleton.reweighted.freq.spec + ( 1 - p_new ) * this.freq.spec
	
	# return ( final.freq.spec )
	
# }




### new mutations during sweep correct way

expected.freq.times.standing.w.sweep <- function ( nsam , N , r , f , s , my.StirlingNumbers = 0 ) {
	#recover()
	if ( all ( my.StirlingNumbers == 0 ) ) {
		my.StirlingNumbers<-StirlingNumbers(nsam)    ##Usigned Stirling numbers of 1st kind. ma
	}
	ESF.prob.k <- EwensDist( n=nsam , N =N, r=r , distance=1 , f=f)
	ESF.prob.k <- rbind ( c ( 1 , rep ( 0 , nsam ) ) , cbind ( rep ( 0 , nsam ) , ESF.prob.k ) )
	ESF.condprob.k<-EwensCondDist( n=nsam , N =N, r=r , distance=1 , f=f)
	ESF.condprob.k <- rbind ( c ( rep ( 0 , nsam + 1 ) ) , cbind ( rep ( 0 , nsam ) , ESF.condprob.k ) )
	T_f <- log ( (2*N -1 ) * ( 1 - f ) / f ) / s
	my.logistic <- function ( x ) 1 / (2 * N  ) * exp(s * x ) / ( 1 + 1 / (2 * N  ) * ( exp(s * x )  - 1 ) )
	T_sf <- integrate ( my.logistic , 0 , T_f )$value
	P_NR <- exp ( - r * T_sf )
	expected.t.l<-rep(NA,nsam-1)
	p_l_given_k <- array ( 0 , dim = c ( nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 ) )
	H <- array ( 0 , dim = c ( nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 ) )
	cond.freq.specs <- array ( 0 , dim = c ( nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 ) )
	freq.specs <- matrix ( 0 , nrow = nsam , ncol = nsam )
	for ( i in 2 : nsam ) {
		freq.specs [ 1 : ( i - 1 ) , i ] <- ( 1 / ( 1 : ( i - 1 ) ) ) 
	}
	freq.specs <- t ( freq.specs )
	my.freq.specs <- array ( 0 , dim = c ( nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 , nsam + 1 ) )
	z <- 1
	my.list <- list ()
	for ( l in 1 : ( nsam - 1 ) ) {
	#	recover()
	#	terms.in.sum<-rep(0,nsam)
		for ( i in 0 : nsam ) {
			
			for ( k in min ( 1 , i ) : i ) {
				terms.given.j <- matrix ( 0 , ncol = nsam , nrow = nsam )
				for ( j in 1 : min ( k + nsam - i , l ) ) {
					if ( max ( 0 , ( j - k ) ) > min ( l , nsam - i , j ) ) next
					g.sec <- seq ( max ( 0 , ( j - k ) ) ,  min ( l , nsam - i , j ) , 1 )
					for ( g in g.sec ) {
						# if the number of lineages sitting under the beneficial mutatation is smaller than the number we need to get to l, the prob is zero
						if ( i < ( l - g ) ) next
						
						
						p_l_given_k [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] <- plbarjkn ( l - g , j - g , k , i , my.StirlingNumbers )
						

						H [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] <- choose ( nsam - i , g ) * choose ( k , j - g ) / choose ( k + nsam - i , j )
						if ( i == nsam & k == 1 & l == 1) {
							new.freqs <- f + 1/(2*N)*nsam * T_f/2
						} else if ( i < nsam & k == 1 & l == 1 ) {
							new.freqs <- 1/(2*N)*nsam * T_f/2
						} else if ( i == nsam & k == 1 & l > 1 ) {
							new.freqs <- f / l
						} else {
							new.freqs <- 0
						}
						
						cond.freq.specs [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] <- freq.specs [ k + nsam - i , j ] * H [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] * p_l_given_k [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] + new.freqs
						
						if ( g > 0 | i != nsam ) {
							my.freq.specs [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] <- dbinom ( i , nsam , P_NR ) * ESF.prob.k [ i + 1 , k + 1 ] * cond.freq.specs [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ]
						
						} else{
							my.freq.specs [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ] <- dbinom ( i , nsam , P_NR ) * ESF.prob.k [ i  + 1 , k + 1 ] * cond.freq.specs [ k + 1 , j + 1 , g + 1 , i + 1 , l + 1 ]
						}
					}
				}
			}
		}
	#	expected.t.l[l]<-sum(terms.in.sum)
	}
	# singleton.reweighted.times <- 1/(1:(nsam-1)) + c ( nsam * T_f / ( 2 * N * f ) , rep ( 0 , nsam - 2 ) )
	# singleton.reweighted.freq.spec <- singleton.reweighted.times / sum ( singleton.reweighted.times )
	this.freq.spec <- numeric ( nsam - 1 )
	for ( l in 2 : ( dim ( my.freq.specs ) [ 5 ] - 1 ) ) {
		this.freq.spec [ l - 1 ] <- sum ( my.freq.specs [ , , , , l ] ) #+ dbinom ( nsam , nsam , P_NR ) * ESF.prob.k [ nsam +1 , 2 ] * singleton.reweighted.freq.spec [  l - 1 ]
	}
	final.freq.spec <- this.freq.spec / sum ( this.freq.spec )
	# my.branch.sums <- cumsum ( singleton.reweighted.times )
	# num <- ESF.prob.k [ nsam + 1 , 2 ] * dbinom ( 10 , 10 , P_NR ) * f * my.branch.sums [ nsam - 1 ]
	# other <- numeric()
	# P_R <- 1 - P_NR
	# t <- 1
	# for ( j in 2 : ( nsam - 1 ) ) {
		# for ( k in 0 : j ) {
			# for ( i in 0 : ( j - k ) ) {
				# if ( i < 0 | (k + i) <= 1 ) { next }
				
				# other [ t ] <- dbinom ( i , nsam , P_R ) * ESF.prob.k [ nsam - i + 1  , k + 1 ] * my.branch.sums [ j ]
				# t <- t + 1
			# }
		# }
	# }
	# p_new <- num / ( num + sum ( other ) )
	# final.freq.spec <- p_new * singleton.reweighted.freq.spec + ( 1 - p_new ) * this.freq.spec
	
	return ( final.freq.spec )
	
}









freq.spec.de.novo <- function ( nsam , N , r , s ) {
	#recover()
	p <- 1-exp ( - 2 * r / s *log ( 2*N*s) )
	spec.by.n.recombs <- matrix ( 0 , nrow = nsam -1 , ncol = nsam )
	for ( i in 1 : nrow ( spec.by.n.recombs ) ) {
		
		for ( j in 1 : ncol ( spec.by.n.recombs ) ) {
			
			p_i_over_j_plus1 <- (1 / i ) / sum ( 1 / ( 1 : ( j + 1 ) ) )
			p_core_not_hit <- ifelse ( i > j , 0 , choose ( j , i ) / choose ( min ( j + 1 , nsam ) , i ) )
			
			if ( j < nsam & i >= nsam -j ) {
				
				p_core_hit <- choose ( j , i-(nsam -j )+1 ) / choose ( j + 1 , i-(nsam -j )+1 )
				p_inj1_over_j_plus1 <- ( 1 / ( i - ( nsam - j ) + 1 ) ) / sum ( 1 / ( 1 : ( j + 1 ) ) )
			
			} else {
				
				p_core_hit <- 0
				p_inj1_over_j_plus1 <- 0
				
			}
			
			spec.by.n.recombs [ i , j ] <- p_i_over_j_plus1 * p_core_not_hit + p_inj1_over_j_plus1 * p_core_hit
			
		}
		
	}
	
	norm.spec.by.n.recombs <- t ( t ( spec.by.n.recombs ) / colSums ( spec.by.n.recombs ) )
	rec.probs <- dbinom ( 1:10 , 10 , p )
	freq.spec <- colSums ( t ( norm.spec.by.n.recombs ) * rec.probs )
	freq.spec[1] <- freq.spec[1] + dbinom ( 0 , 10 , p )
	
	return ( freq.spec )
	
}





if ( FALSE ) {
blah <- expected.freq.times.standing.w.sweep ( nsam = 50 , N = 10000 , r = 0.000001 , f = 0.05 , s = 0.05 )
blah2 <- freq.spec.de.novo ( 10 , 10000 , 0.0001 , 0.01 )

#### so, how many sweeps actually fit our model, anyway?


N <- 100000

spec.cond.on.fix <- function ( x , N , s ) ( 1 - exp ( -4*N*s*x) ) / (x * ( 1 - exp (-4*N*s) ) )
neutral.spec <- function ( x ) 1/x

SSVProb <- function ( N , s ) {
	#recover()
	blah1 <- list ()
	blah2 <- list ()
	normalizing.constant <- list ()
	my.prob1 <- list ()
	my.prob2 <- list ()
	upper <- c ( seq ( max ( 0.05 , 5 / (4*N*s) + 0.01 ) , 0.95 , by = 0.01 ) , ( 2 * N - 1 ) / ( 2 * N ) )
	for ( i in 1 : length ( upper ) ) {
		blah1 [[ i ]] <- integrate ( function ( y ) spec.cond.on.fix ( x = y , N = N , s = s ), max ( 5 / ( 4 * N * s ) , 500 / (2*N) ) , upper [ i ] )$value
		blah2 [[ i ]] <- integrate ( function ( y ) spec.cond.on.fix ( x = y , N = N , s = s ), 1 / ( 2 * N ) , max ( 5 / ( 4 * N * s ) , 500 / ( 2 * N ) ) )$value
		normalizing.constant [[ i ]] <- integrate ( function ( y ) spec.cond.on.fix ( x = y , N = N , s = s ), 1 / ( 2 * N ) , upper  [ i ] )$value
		my.prob1 [[ i ]] <- ( blah1 [[ i ]] ) / ( normalizing.constant [[ i ]] )
		my.prob2 [[ i ]] <- ( blah2 [[ i ]] ) / ( normalizing.constant [[ i ]] )	
	}
	my.soft.probs <- cbind ( upper , unlist ( my.prob1 ) , unlist ( my.prob2 ) )
	
}
my.soft.probs <- list ()
my.s <- c ( 0.001 , seq ( 0.01 , 0.1 , by = 0.02 ) )
for ( i in 1 : length ( my.s ) ) {

	my.soft.probs [[ i ]] <- SSVProb ( N = 100000 , s = my.s [ i ] )

}
plot ( NA , xlim = c ( 0 , 0.15 ) , ylim = c ( 0 , 1 ) , bty = "n" , xlab = "Upper Bound on f" , ylab = "Probability our Model Applies" )
for ( i in 1 : length ( my.soft.probs ) ) {
	
	lines ( x = my.soft.probs [[ i ]] [ , 1 ] , my.soft.probs [[ i ]] [ , 2 ] , lwd = 2 , col = i )
	
}
legend ( "bottomright" , legend = my.s , col = 1:length(my.s) , bty = "n" , lwd = 2 )



my.freqs <- seq ( 1 / ( 2 * N ) , ( 2 * N -1 ) / ( 2 * N ) , by = 0.0001 )
cond.fix.spec <- spec.cond.on.fix ( my.freqs , N , s )
neutral.spec <- neutral.spec ( my.freqs )

plot ( x = my.freqs , y = cond.fix.spec / sum ( cond.fix.spec ) , type = "l" , lwd = 2 , xlim = c ( 0 , 1 ) )
lines ( x = my.freqs , y = neutral.spec / sum ( neutral.spec ) , type = "l" , lwd = 2 , col = "red" )








my.freqs.specs <- run.ms.f ( runs = my.test , f = 0.05 , s = 0.05 , n.sam = 12 , N = 10000 , path = "" , ext = "fr.spec", get.site.density = FALSE , recom = 100 )





##setup for BATCH put inside if FALSE after
recoms<- seq(1,200,length=20)
f.freq.spec<-list()
for(f.index in 1:length(real.fs)){
	print(f.index)
	collect.freq.spec<-numeric()
	for(recom in recoms){
		cat("recom",recom,"\n")
		my.freqs.specs<- run.ms.f ( runs = my.runs [[ f.index ]] [[ 1 ]] ,f.index=1, n.sam = 10 , N = 10000 , path = path,get.site.density = FALSE , recom = recom)
		collect.freq.spec<-cbind(collect.freq.spec,rowMeans(my.freqs.specs))
		save(file=paste(path,"Sims/Freq_spec.Robj",sep=""),collect.freq.spec,f.freq.spec)
	}
	f.freq.spec[[f.index]]<-collect.freq.spec
	save(file=paste(path,"Sims/Freq_spec.Robj",sep=""),collect.freq.spec,f.freq.spec)
}


norm.freqs<-apply(collect.freq.spec,2,function(x){x/sum(x)})

N<-10000

exp.times<-matrix(NA,nrow=200,ncol=10)

recs<-seq(1,200,length=200)/(4*N)

for(my.rec.index in 1:length(recs)){
	exp.times[my.rec.index,]<- expected.freq.times.standing(n=10,N=N,r=recs[my.rec.index],distance=1,f=0.01)
}

}
