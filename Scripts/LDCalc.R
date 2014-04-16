library ( expm )
library ( combinat )
library ( plyr )
library ( RColorBrewer )

f <- 0.10
N <- 100000
r_x <- 10^-4
r_y <- 10^-4
s <- 0.01
#### selected portion for ld across site
SelTransMat <- function ( r_x , r_y , s = f , N ) {
	p_x <- 1 - exp ( - r_x / s * log ( ( N - 1 ) / f - ( N - 1 ) ) )
	q_x <- 1 - p_x
	p_y <- 1 - exp ( - r_y / s * log ( ( N - 1 ) / f - ( N - 1 ) ) )
	q_y <- 1 - p_y
	
	trans.mat <- matrix ( 	c ( 
								q_x^2*q_y^2 , 0 , 2*p_x*q_x*q_y^2 , 0 , 0 , p_x^2*q_y^2 , ## 6
								2*q_x^2*p_y*q_y , 0 , 2*p_x*q_x*p_y*q_y , 0 , 0 , 2*p_x*q_x*p_y*q_y , ## 12
								0 , 2*p_x^2*p_y*q_y , q_x^2*p_y^2 , 0 , 2*p_x*q_x*p_y^2 , 0 , ## 18
								0 , p_x^2*p_y^2
								) ,
							c ( 0 , q_x^2*q_y^2 , p_x*q_x*q_y^2 , 0 , p_x*q_x*q_y^2 , p_x^2*q_y^2 , ## 6
								q_x^2*q_y*p_y , 0 , p_x*q_x*p_y*q_y , q_x^2*q_y*p_y , 0 , 3*q_x*p_x*q_y*p_y , ## 12
								0 , 2*p_x^2*q_y*p_y , q_x^2*p_y^2 , 0 , 2*p_x*q_x*p_y^2 , 0 , ## 18
								0 , p_x^2*p_y^2
								)
						)	
	
}



#### netural portion for ld across site.
MakeMatrices <- function ( states ) {
	#recover()
	x <- y <- numeric ( 8 )
	x [ states[[1]] ] <- y [ states[[2]] ] <- 1
	S <- as.data.frame ( cbind ( x [ 1:4 ] , y [ 1:4 ] ) )
	W <- as.data.frame ( cbind ( x [ 5:8 ] , y [ 5:8 ] ) )
	S <- S [ order ( S$V1 , S$V2 , decreasing = T ) , ]
	W <- W [ order ( W$V1 , W$V2 , decreasing = T ) , ]
	return ( list ( S = S , W = W ) )
}

one.locus.states <- c ( combn ( seq ( 1,8,1) , 2 , simplify = F ) , list ( 1,2,5,6 ) )
#one.locus.states <- split ( one.locus.states , rownames ( one.locus.states ) ); names ( one.locus.states ) <-NULL
two.locus.states <- expand.grid ( one.locus.states , one.locus.states )
state.matrices <- apply ( t ( two.locus.states ) , 2 , MakeMatrices )

i = 1 
while ( i < ( length ( state.matrices ) + 1 ) ) {
	j = i + 1
	current <- state.matrices [[ i ]]
	while ( j < length ( state.matrices ) + 1 ) {
		test <- state.matrices [[ j ]]
		if ( all ( c ( current$S == test$S ) & c (current$W == test$W) ) ) {
			state.matrices [[ j ]] <- NULL
		} else {
			j = j + 1
		}
	}
	#if ( i == 44 ) stop()
	if ( sum ( unlist ( current ) ) < 4 ) {
		state.matrices [[ i ]] <- NULL
	}
	i = i + 1
	#head(state.matrices [ - c (1:i)])
}
GetNumerator <- function ( state.one , state.two , f, N , r_x , r_y ) {
	#recover()
	test.state <- rbind ( state.one$S [ 1 , ] , state.one$W )
	test.state <- test.state [ order ( test.state$V1 , test.state$V2 , decreasing = T ) , ]
	test.state <- test.state [ 1 : 4 , ]
	if (
				# if only one lineage is left on selected background, mutate
				sum ( rowSums ( state.one$S ) > 0 ) == 1 &
				sum ( rowSums ( state.two$S ) > 0 ) == 0 &
				all ( test.state == state.two$W )
			) {
		return ( 1 )
	} else if (
				sum ( rowSums ( state.one$S ) > 0 ) == 1
			) {
		return ( 0 )
			}
	else if ( 
				# if all lineages reside on wild type background, do nothing w.p. one
				all ( state.one$W == state.two$W ) &
				sum ( state.one$S ) == 0 & 
				sum ( state.two$S ) == 0
			) {
		return ( 1 )
	} else if ( 
				# otherwise, cannot transition to identical state; probability zero
				all ( c ( state.one$S == state.two$S ) , c ( state.one$W == state.two$W ) ) 
			) {
		return ( 0 )
	} else if ( 
				## recombination at x locus off of the selected background into wt background
				all ( colSums ( state.one$S ) - c ( 1 , 0 ) == ( colSums ( state.two$S ) ) ) &
				all ( colSums ( state.one$W ) + c ( 1 , 0 ) == ( colSums ( state.two$W ) ) ) &
				sum ( rowSums ( state.one$W ) > 0 ) + 1 == sum ( rowSums ( state.two$W ) > 0 ) &
				(
					sum ( rowSums ( state.one$S ) > 0 ) == sum ( rowSums ( state.two$S ) > 0 ) |
					sum ( rowSums ( state.one$S ) > 0 ) - 1 == sum ( rowSums ( state.two$S ) > 0 )
				)
				## when there are exchangeable alleles at the x locus
			) {
				if ( 
					(
						all ( rowSums ( state.one$S ) == c ( 2 , 2 , 0 , 0 ) ) & 
						all ( colSums ( state.one$S ) == c ( 2 , 2 ) )
					) |
					(	all ( rowSums ( state.one$S ) == c ( 1 , 1 , 1 , 0 ) ) & 
						all ( colSums ( state.one$S ) == c ( 2 , 1 ) )
					) |
					(	all ( rowSums ( state.one$S ) == c ( 1 , 1 , 1 , 1 ) ) & 
						all ( colSums ( state.one$S ) == c ( 2 , 2 ) )
					)
				) {
					return ( 2 * r_x * ( 1 - f ) )
				} else {
					# and when there are not
					return ( r_x * ( 1 - f ) )
				}
	} else if ( 
				## recombination at y locus off of the selected background into wt background
				all ( colSums ( state.one$S ) - c ( 0 , 1 ) == ( colSums ( state.two$S ) ) ) &
				all ( colSums ( state.one$W ) + c ( 0 , 1 ) == ( colSums ( state.two$W ) ) ) &
				sum ( rowSums ( state.one$W ) > 0 ) + 1 == sum ( rowSums ( state.two$W ) > 0 ) &
				(
					sum ( rowSums ( state.one$S ) > 0 ) == sum ( rowSums ( state.two$S ) > 0 ) |
					sum ( rowSums ( state.one$S ) > 0 ) - 1 == sum ( rowSums ( state.two$S ) > 0 )
				)
				## when there are exchangeable alleles at the y locus
			) {
				if ( 
					(
						all ( rowSums ( state.one$S ) == c ( 2 , 2 , 0 , 0 ) ) & 
						all ( colSums ( state.one$S ) == c ( 2 , 2 ) )
					) |
					(	all ( rowSums ( state.one$S ) == c ( 1 , 1 , 1 , 0 ) ) & 
						all ( colSums ( state.one$S ) == c ( 1 , 2 ) )
					) |
					(	all ( rowSums ( state.one$S ) == c ( 1 , 1 , 1 , 1 ) ) & 
						all ( colSums ( state.one$S ) == c ( 2 , 2 ) )
					)
				) {
					return ( 2 * r_y * ( 1 - f ) )
				} else {
					# and when there are not
					return ( r_y * ( 1 - f ) )
				}
	} else if ( 
				## coalescence when 4 possible events give the same transition
				sum ( rowSums ( state.one$S ) > 0 ) == 4 &
				sum ( rowSums ( state.two$S ) > 0 ) == 3 &
				all ( colSums ( state.one$S ) == c ( 2 , 2 ) ) &
				all ( colSums ( state.two$S ) == c ( 2 , 2 ) ) &
				all ( state.one$W == state.two$W )
			) {
		return ( 4 / ( 2 * N * f ) )
	} else if (
				## coalescence when two possible events give the same transition
				sum ( rowSums ( state.one$S ) > 0 ) == 3 &
				sum ( rowSums ( state.two$S ) > 0 ) == 2 &
				( (
					all ( colSums ( state.one$S ) == c ( 2 , 1 ) ) &
					all ( colSums ( state.two$S ) == c ( 2 , 1 ) ) 
				) |
				(
					all ( colSums ( state.one$S ) == c ( 1 , 2 ) ) &
					all ( colSums ( state.two$S ) == c ( 1 , 2 ) )
				) ) &
				all ( state.one$W == state.two$W )
			) {
		return ( 2 / ( 2 * N * f ) )
	} else if (
				# all other coalescent events on the selected background occur w.p. 1/2Nf
				sum ( rowSums ( state.one$S ) > 0 ) - 1 == sum ( rowSums ( state.two$S ) > 0 ) &
				all ( state.one$W == state.two$W )
			) {
		return ( 1 / ( 2 * N * f ) )		
	} else if (
				# coalescence on the wild type background when alleles split evenly across backgrounds at both loci
				all ( rowSums ( state.one$W ) == c ( 1 , 1 , 0 , 0 ) ) &
				all ( colSums ( state.one$W ) == c ( 1 , 1 ) ) &
				all ( rowSums ( state.two$W ) == c ( 2 , 0 , 0 , 0 ) ) &
				all ( colSums ( state.two$W ) == c ( 1 , 1 ) ) &
				all ( state.one$S == state.two$S )
			){
		return ( 1 / ( 2 * N * ( 1 - f ) ) )
	} else {
		return (0)
	}
}

GetDenominator <- function ( state.one , state.two , f, N , r_x , r_y ) {
	#recover()
	test.state <- rbind ( state.one$S [ 1 , ] , state.one$W )
	test.state <- test.state [ order ( test.state$V1 , test.state$V2 , decreasing = T ) , ]
	test.state <- test.state [ 1 : 4 , ]
	if (
				sum ( rowSums ( state.one$S ) > 0 ) == 1 &
				sum ( rowSums ( state.two$S ) > 0 ) == 0 &
				all ( test.state == state.two$W )
			) {
		return ( 1 )
	} else if ( 
				# if all lineages reside on wild type background, do nothing w.p. one
				all ( state.one$W == state.two$W ) &
				sum ( state.one$S ) == 0 & 
				sum ( state.two$S ) == 0
			) {
		return ( 1 )
	} else {
		rec_x <- colSums ( state.one$S ) [ 1 ] * r_x * ( 1 - f )
		rec_y <- colSums ( state.one$S ) [ 2 ] * r_y * ( 1 - f )
		coal_s <- choose ( sum ( rowSums ( state.one$S ) > 0 ) , 2 ) / ( 2 * N * f )
		coal_w <- choose ( sum ( rowSums ( state.one$W ) > 0 ) , 2 ) / ( 2 * N * ( 1 - f ) )
		return ( rec_x + rec_y + coal_s + coal_w )
	}
}
LDCalc <- function ( f , N , r_x , r_y , state.matrices ) {
	
	#recover()
	
	
	
	
	## Neutral phase transition matrix
	transition.mat <- matrix ( NA , nrow = 20 , ncol = 20 )
	for ( i in 1 : 20 ) {
		for ( j in 1 : 20 ) {
			num <- GetNumerator (state.matrices [[ i ]] , state.matrices [[ j ]] , f , N , r_x , r_y )
			denom <- GetDenominator (state.matrices [[ i ]] , state.matrices [[ j ]] , f , N , r_x , r_y )
			transition.mat [ i , j ] <- ifelse ( num == 0 & denom == 0 , 0 , num / denom )
		}
	}
	
	# default neutral values
	tt_ijij_w <- ( 36 + 14 * N * ( r_x + r_y ) + ( N * ( r_x + r_y ) )^2 ) / ( 18 + 13 * N * ( r_x + r_y ) )
	tt_ijik_w <- ( 24 + 13 * N * ( r_x + r_y ) + ( N * ( r_x + r_y ) )^2 ) / ( 18 + 13 * N * ( r_x + r_y ) )
	tt_ijkl_w <- ( 22 + 13 * N * ( r_x + r_y ) + ( N * ( r_x + r_y ) )^2 ) / ( 18 + 13 * N * ( r_x + r_y ) )
	
	
	phi_a <- ( c ( 1 , rep ( 0 , 19 ) ) %*% ( transition.mat %^%10 ) ) [ 18 : 20 ]
	phi_b <- ( c ( 0 , 1 , rep ( 0 , 18 ) ) %*% ( transition.mat %^%10 ) ) [ 18 : 20 ]
	phi_c <- ( c ( 0 , 0 , 0 , 1 , rep ( 0 , 16 ) ) %*% ( transition.mat %^%10 ) ) [ 18 : 20 ]
	# reweighted selected values
	tt_ijij_s <- sum ( phi_a * c ( tt_ijij_w , tt_ijik_w , tt_ijkl_w ) )
	tt_ijik_s <- sum ( phi_b * c ( tt_ijij_w , tt_ijik_w , tt_ijkl_w ) )
	tt_ijkl_s <- sum ( phi_c * c ( tt_ijij_w , tt_ijik_w , tt_ijkl_w ) )
	
	( tt_ijij_s - 2 * tt_ijik_s + tt_ijkl_s ) / ( tt_ijkl_s )
	
}


f <- 0.10
N <- 100000
r_x <- 10^-4
r_y <- 10^-4
LDCalc ( f = 0.05 , N = 100000 , r_x = 10^-4 , r_y = 10^-4 , state.matrices )

my.freqs <- c ( 0.001 , seq ( 0.01 , 0.1 , 0.01 ) )
my.dist <- seq ( 0 , 10^-3 , length.out = 100 ) 
my.ld <- list ()
for ( j in 1 : 11 ) {
	my.ld [[ j ]] <- numeric ( 100 )
	for ( i in 1 : 100 ) {
		my.ld [[ j ]] [ i ] <-  LDCalc ( f = my.freqs [ j ] , N = 100000 , r_x = my.dist [ i ] , r_y = my.dist [ i ] , state.matrices )
	}
}


tt_ijij_w <- ( 36 + 14 * N * ( r_x + r_y ) + ( N * ( r_x + r_y ) )^2 ) / ( 18 + 13 * N * ( r_x + r_y ) )
tt_ijik_w <- ( 24 + 13 * N * ( r_x + r_y ) + ( N * ( r_x + r_y ) )^2 ) / ( 18 + 13 * N * ( r_x + r_y ) )
tt_ijkl_w <- ( 22 + 13 * N * ( r_x + r_y ) + ( N * ( r_x + r_y ) )^2 ) / ( 18 + 13 * N * ( r_x + r_y ) )
neutral.ld <- ( tt_ijij_w - 2 * tt_ijik_w + tt_ijkl_w ) / ( 1 + tt_ijkl_w - N^2)





save ( my.ld , file = "LDvals.Robj")
cols <- rainbow ( 10 )
plot ( my.ld [[ 1 ]] , type = "l"  )
matplot ( do.call ( cbind, my.ld ) , type = "l" , xlab = "4Nr" , xaxt = "n" , col = brewer.pal ( 11 , "Spectral" ) , lty = 1 , bty = "n" , ylab = "sigma^2_d" )
axis ( 1 , at = c ( 0 , 49 , 99 ) , labels = c ( "0" , "100" , "200" ) )


on <- function ( ) recover.on <<- T
off <- function ( ) recover.on <<- F







