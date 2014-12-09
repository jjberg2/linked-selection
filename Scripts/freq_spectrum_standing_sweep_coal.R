setwd ( "~/Documents/Academics/StandingSweeps/" )
source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
source('~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R', chdir = TRUE)
####frequency spectrum

#real.fs <- c ( 0.001 ,  0.01 , 0.05 , 0.1 )

#my.runs <- lapply ( real.fs , function ( x ) SweepFromStandingSim ( N = 10000 , s = 0.01 , f = x , reps = 10 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T ) )


#save ( my.runs , file = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Sims/11000freq.trajectories.Rdata")


##load( file = "~/../Shared/11000freq.trajectories.Rdata")  #graham's work machine
#load ( file = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Sims/11000freq.trajectories.Rdata" )    ##Jeremy's machine

#this.comp <-Sys.info()["nodename"]
#setwd ( "~/Documents/Academics/StandingSweeps/" )  ##Jeremy's machine
##path = "~/Dropbox/Linked_selection_models/Soft_sweeps_coal/LinkedSelection/" #graham's work machine






####freq. spectrum w. no rec. during sweeps.

expected.freq.times.standing<-function(nsam,N,r,distance,f){
	#recover()
	#my.StirlingNumbers<-StirlingNumbers(n)
	ESF.prob.k<-EwensCondDist( n=nsam , N =N, r=r , distance=1 , f=f) # ,stirling.numbers=my.StirlingNumbers)    ### is of form [n,k]
	my.StirlingNumbers<-StirlingNumbers(nsam)    ##Usigned Stirling numbers of 1st kind. ma
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
				terms.in.sum [ k , j, l  ] <- ESF.prob.k [ nsam , k ] * p_l_given_k [ k , j , l ] * freq.specs [ k , j ]
		
			}
		}
	}

	return(terms.in.sum)
}
if ( FALSE ) {

pdf ( file = "Figures/Spec_ProblemsCond.pdf" , width = 15 , height = 40 )
par ( mfrow = c ( 5 , 3 ) )

my.N <- 10000
my.rs <- rev ( c ( 0.0001 , 0.001 , 0.01 , 0.1 , 0.5 ) )
my.r <- 0.001
blah <- list ()
specs <- list ()
norm.spec <- list ()


my.runs <-  SweepFromStandingSim ( N = 10000 , s = 0.05 , f = 0.05 , reps = 1000 , no.sweep = TRUE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T )
my.test <- matrix ( 0.05 , nrow = 1000 , ncol = 1000 )




	my.freqs.specs <- run.ms.f ( runs = my.runs [[ 1 ]] , f = 0.05 , s = 0.05 , n.sam = 12 , N = my.N , path = "" , ext = "fr.spec", get.site.density = FALSE , recom = my.r / ( 4 * my.N ) )
	save ( my.freqs.specs , file = "Sims/freq.spec.05.05.10000.nosweep.n12.Robj" )

expected.freq.times.standing(nsam=12,N=my.N,r = my.r , f = 0.05 )



for ( i in 1 : length ( my.rs ) ) {
	blah [[ i ]] <- expected.freq.times.standing(nsam=10,N=10000,r = my.rs [ i ] , f = 0.01 )
	specs [[ i ]] <- colSums ( blah [[ i ]] , dims = 2)

	# my.spec <- ifelse ( specs [[ i ]] == 0 , NA , specs [[ i ]] )
	# matplot ( t ( my.spec ) [ 1:9 , 2 : 10 ] , type = "o" , lty = 1 , col = 1:9 , pch = 20 , ylim = c ( 0 , 0.4 ) , xlim = c ( 1 , 9 ) )

	# norm.spec [[ i ]] <- specs [[ i ]] / rowSums ( specs [[ i ]] )
	# norm.spec [[ i ]] <- ifelse ( norm.spec [[ i ]] == 0 , NA , norm.spec [[ i ]] )
	# matplot ( t ( norm.spec [[ i ]]) [ 1:9 , 2 : 10 ] , type = "o" , lty = 1 , col = 1:9 , pch = 20 )

	# norm.spec [[ i ]] <- specs [[ i ]] / sum ( specs [[ i ]] )
	# norm.spec [[ i ]] <- ifelse ( norm.spec [[ i ]] == 0 , NA , norm.spec [[ i ]] )
	# matplot ( t ( norm.spec [[ i ]] ) [ 1:9 , 2 : 10 ] , type = "o" , lty = 1 , col = 1:9 , pch = 20 )

	# if ( i == 1 ) {
		# legend ( "topright" , legend = sapply ( c ( 2 : 10 ) , function ( x ) paste ( "k = " , x , sep = "" ) ) , col = 1:9 , bty = "n" , lty = 1 )
	# }
}
dev.off()




pdf ( )


##### some component function to make the code cleaner
plbarjkn <- function ( l , j , k , n ) {
	#recover()
	####################
	#### special cases ####
	####################
	
	
	if ( l == 0 & j == 0 ) 	return ( 1 )
	if ( j == k & n == l )	return ( 1 )
	if ( j == k & n != l )	return ( 0 )
	if ( l > 0 & j == 0 ) 	return ( 0 )
	if ( l == 0 & j > 0 ) 	return ( 0 )
	if ( l > j + n - k )		return ( 0 )
	if ( k == 1 ) {
		if ( n == l & j == 1 ) return ( 1 )
		if ( l != n ) return ( 0 )
		else stop ()
	}


	####################
	#### standard case ####
	####################
	stirling.bit <- my.StirlingNumbers[ l , j ] * my.StirlingNumbers[n - l , k - j ]  / my.StirlingNumbers[ n , k ]
	
	stirling.bit * choose ( nsam , l ) / choose ( k , j )
	
}





####freq. spectrum with rec. during sweeps.

expected.freq.times.standing.w.sweep<-function(nsam,N,r,distance,f){
	recover()
	#my.StirlingNumbers<-StirlingNumbers(n)
	ESF.prob.k<-EwensDist( n=nsam , N =N, r=r , distance=1 , f=f) # ,stirling.numbers=my.StirlingNumbers)    ### is of form [n,k]
	ESF.condprob.k<-EwensCondDist( n=nsam , N =N, r=r , distance=1 , f=f) # ,stirling.numbers=my.StirlingNumbers)    ### is of form [n,k]
	my.StirlingNumbers<-StirlingNumbers(nsam)    ##Usigned Stirling numbers of 1st kind. ma
	expected.t.l<-rep(NA,nsam-1)
	p_l_given_k <- array ( 0 , dim = c ( nsam , nsam , nsam , nsam + 1 , nsam ) )
	freq.specs <- matrix ( 0 , nrow = nsam , ncol = nsam )
	for ( i in 2 : nsam ) {
		freq.specs [ 1 : ( i - 1 ) , i ] <- ( 1 / ( 1 : ( i - 1 ) ) ) / ( sum ( 1 / ( 1 : ( i - 1  )  ) ) )
	}
	freq.specs <- t ( freq.specs )
	terms.in.sum <- array ( 0 , dim = c ( nsam , nsam , nsam , nsam ) )
	H <- array ( 0 , dim = c ( nsam , nsam , nsam , nsam , nsam ) )
	jg.tracker <- array ( NA , dim = c ( nsam , nsam , nsam , nsam , nsam ) )
	kg.tracker <- array ( NA , dim = c ( nsam , nsam , nsam , nsam , nsam ) )
	for ( l in 1 : ( nsam - 1 ) ) {
	#	recover()
	#	terms.in.sum<-rep(0,nsam)
		for ( i in 0 : nsam ) {
			for ( k in 1 : i ) {
				terms.given.j <- matrix ( 0 , ncol = nsam , nrow = nsam )
				for ( j in 1 : min ( k + nsam - i  - 1 , l ) ) {
					if ( max ( 0 , ( j - k ) ) > min ( l , nsam - i , j ) ) next
					g.sec <- seq ( max ( 0 , ( j - k ) ) ,  min ( l , nsam - i , j ) , 1 )
					for ( g in g.sec ) {
						# if the number of lineages sitting under the beneficial mutatation is smaller than the number we need to get to l, the prob is zero
						if ( i < ( l - g ) ) next
						
						p_l_given_k [ k , j , g , i + 1 , l ] <- plbarjkn ( l - g , j - g , k , i )
						

						H [ k , j , g , i + 1 , l ] <- choose ( nsam - i , g ) * choose ( k , j - g ) / choose ( k + nsam - i , j )
						
						freq.specs [ k + nsam - i , j ] * H [ k , j , g , i + 1 , l ] * p_l_given_k [ k , j , g , i , l ]
						
						if ( FALSE ) {
						terms.given.j [ j , k ] <-  [ k , j , l ]*freq.specs[ k + nsam - i  , j ] * H [ g , j , k , nsam - i ]
					
						terms.in.sum [ k , j, l , i  ] <- ESF.prob.k [ nsam , k ] * p_l_given_k [ k , j , l ] * freq.specs [ k , j ]
						stopifnot( H <= 1 , H >= 0 )
						}

						jg.tracker [ k , j , g , i ,  l ]	 <- j - g
						kg.tracker [ k , j , g , i , l ]	 <- k - ( j - g )
					}
				}
			}
		}
	#	expected.t.l[l]<-sum(terms.in.sum)
	}

	return(terms.in.sum)
}

expected.freq.times.standing.w.sweep ( nsam = 12 , N = 10000 , r = 0.0001 , f = 0.05 )












my.freqs.specs<- run.ms.f ( runs = my.test , f = 0.05 , s = 0.05 , n.sam = 12 , N = 10000 , path = "" , ext = "fr.spec", get.site.density = FALSE , recom = 100 )





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
