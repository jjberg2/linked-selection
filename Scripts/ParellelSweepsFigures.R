

MakeParellelSweepFigure <- function ( N , B , f , r.vect , T_bot , s ) {
	
	#recover()
	T_sel <- log ( N*B * (1-f)/f) / s
	
	# probability of recombining out in ancestral pop before coalescing, conditional on suriviving that far back
	P_A <- 4 * N * r.vect * f * ( 1 - f ) / (1 + 4 * N * r.vect * f * ( 1 - f ) )
	
	# probability of either lineage recombining off the selected background during the sweep
	T_int <- integrate ( function ( t ) {1 - exp ( s*t ) / (N*B - 1 + exp ( s * t ) )} , 0 , T_sel )$value
	P_S <- 1 - exp ( - 2*r.vect*T_int )

	# probability of both lineages escaping the bottleneck without recombining if sampled from different populations
	P_BD <- 1 - exp ( -2 * r.vect * T_bot * ( 1- f ) )
	
	# probability that one of two lineages sampled from the same population recombine during the bottleneck without first coalescing
	P_BS <- r.vect * ( 1 - f ) * sapply ( r.vect , function ( r ) integrate ( function ( t ) { exp ( -2 * t * r * ( 1 - f ) ) * exp ( -t/(2*N*f*B) ) } , 0 , T_bot )$value )

	# probability that both lineages escape bottleneck without coalescing or recombining
	P_EBS <- exp ( - 2 * T_bot * r.vect * ( 1 - f ) - T_bot / ( 2 * N * B * f ) )
	
	T_same <- P_S + ( 1 - P_S ) *P_BS + ( 1 - P_S ) * P_EBS * P_A
	T_diff <- P_S + ( 1 - P_S ) *P_BD + ( 1 - P_S ) * (1 - P_BD ) * P_A
	F_st <- ( T_diff - T_same ) / ( T_same + T_diff )
	plot ( r.vect , F_st , type = "l")
	return ( F_st )

	
}

temp <- MakeParellelSweepFigure ( N = 2000000 , B = 1 , f = 0.01 , r.vect = seq ( 0 , 0.02 , 10^-7) , T_bot = 1 , s = 0.1 )

temp <- MakeParellelSweepFigure ( N = 2000000 , B = .1 , f = 0.01 , r.vect = seq ( 0 , 0.02 , 10^-7) , T_bot = 100 , s = 0.1 )

temp <- MakeParellelSweepFigure ( N = 2000000 , B = .01 , f = 0.01 , r.vect = seq ( 0 , 0.02 , 10^-7) , T_bot = 100 , s = 0.1 )


MakeParellelSweepFigure2 <- function ( N , B , f , r.vect , T_0 , s ) {
	
	recover()
	T_sel <- log ( N*B * (1-f)/f) / s
	
	# probability of recombining out in ancestral pop before coalescing, conditional on suriviving that far back
	P_A <- 4 * N * r.vect * f * ( 1 - f ) / (1 + 4 * N * r.vect * f * ( 1 - f ) )
	
	# probability of either lineage recombining off the selected background during the sweep
	T_int <- integrate ( function ( t ) {1 - exp ( s*t ) / (N*B - 1 + exp ( s * t ) )} , 0 , T_sel )$value
	P_S <- 1 - exp ( - 2*r.vect*T_int )

	P_C = exp ( -T_0 / (2*N*B))/(2*N*B)

	T_same <- P_S*(P_A + P_C - P_A*P_C)
	T_diff <- (2-P_C)*(P_S+P_A-P_S*P_A)
	F_st <- ( T_diff - T_same ) / ( T_same + T_diff )
	plot ( r.vect , F_st , type = "l")
	return ( F_st )

	
}

temp <- MakeParellelSweepFigure2 ( N = 2000000 , B = 1 , f = 0.01 , r.vect = seq ( 0 , 0.02 , 10^-7) , T_0 = 1 , s = 0.1 )

