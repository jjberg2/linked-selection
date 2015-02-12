setwd("~/Documents/Academics/StandingSweeps")
source ( "Scripts/SweepFromStandingSim.R")


msHapSims <- function ( runs , n.sam = 2  , f , s , N , path , get.site.density = TRUE , num.sims , len.bp , r.bp , mu.bp , ext = "out" , hap.count.interval , both.sides = FALSE ) {
	#recover()
	options ( "scipen" = 100 , "digits" = 4 )
	f.lab <- strsplit ( as.character ( f ) , "\\." ) [[ 1 ]] [ 2 ]
	s.lab <- strsplit ( as.character ( s ) , "\\." ) [[ 1 ]][ 2 ]
	
	my.file <- paste ( path , "/mssel_f" , n.sam , "." ,  f.lab , "."  , s.lab , "." , N  , "." , ext ,  sep = "" )
	system ( paste ( "rm " , my.file ) )
	#for ( run in 1:5 ) {
	#	load ( paste ( "run_cond_lost_" , run , ".Robj" , sep = "" ) )
	hap.counts <- list ()
	hap.counts.no.sing <- list ()
	window.hap.counts <- list ()
	for ( i in 1:nrow ( runs ) ) {
		
		my.freqs <- c ( runs [ i , runs [ i , ] > 0 ] , 0 )
		my.times <- 0 : ( length ( my.freqs ) - 1 ) / ( 4*N  )

		#recover()
		header.material <- c ( "1" , "1" , paste ( "n:" , length ( my.times ) ) )
		traj.file <- paste ( path , "/my.standing" , "." , f.lab , "." , s.lab , "." , N , ".traj" , sep = "" )
		write ( file = traj.file , header.material )
		write.table ( file = traj.file , cbind ( my.times , my.freqs ) , append = TRUE , sep = "\t" , quot = FALSE , col.nam = FALSE , row.name = FALSE )
		cat( i ," " )
		if ( both.sides == FALSE ) {
			system ( paste ( "Scripts/msseldir/mssel " , n.sam , " " , num.sims , " 0 " , n.sam ,  " " , traj.file ,  " 0 -t " , 2 * N * len.bp * mu.bp , " -r " , 2 * N * len.bp * r.bp , " " , len.bp , " > " , path, "/myseqdata" , sep = "" ) ) 
			seqs <- GetSeqs ( n.sam , num.sims , path )
			hap.counts [ ( i - 1 ) * num.sims + 1:num.sims ] <- lapply ( seqs , CountHaps , len.bp , hap.count.interval )
			hap.counts.no.sing [ ( i - 1 ) * num.sims + 1:num.sims ] <- lapply ( seqs , CountHapsNoSing , len.bp , hap.count.interval )
			marginal.hap.freqs <- Reduce ( "+" , hap.counts ) / length ( hap.counts )
		
		} else {
		
			system ( paste ( "Scripts/msseldir/mssel " , n.sam , " " , num.sims , " 0 " , n.sam ,  " " , traj.file , " " , len.bp/2 , " -r " , 2 * N * len.bp * r.bp , " " , len.bp ," -t " , 2 * N * len.bp * mu.bp  , " " , " > " , path, "/myseqdata" , sep = "" ) ) 
			seqs <- GetSeqs ( n.sam , num.sims , path )
			window.hap.counts [ ( i - 1 ) * num.sims + 1:num.sims ] <- lapply ( seqs  , TwoSideCountHaps , len.bp , hap.count.interval )
			marginal.hap.freqs <- Reduce ( "+" , window.hap.counts ) / length ( window.hap.counts )
		
		}
	}
	
	return ( list ( marginal.hap.freqs , hap.counts ) )
}






GetSeqs <- function ( n , num.sims , path ) {
	#recover()
	a<-system(paste("grep segsites ", path , "/myseqdata",sep=""),intern=TRUE)
	seg.sites<-sapply(a,function(b){as.numeric(strsplit(b,":")[[1]][2])})
	polymorph<- seg.sites>0
	seq.lines<-c(0,cumsum(polymorph*n)[-length(polymorph)])	
	my.seqs <- lapply ( 0 : ( num.sims - 1 ) , function ( iter ) {		
		#recover()
		positions <- read.table ( paste ( path , "/myseqdata" ,sep = "" ) , skip = 5 + 4 * iter + seq.lines [ iter + 1 ] , nrow = 1 )
		seqs.raw <- scan ( paste ( path , "/myseqdata" , sep = "" ) , skip = 6 + 4 * iter + seq.lines [ iter + 1 ] , nline = n , what = character ( ) , quiet = TRUE )
		seqs <- sapply ( seqs.raw , function ( seq ) { as.numeric ( strsplit ( seq , "" ) [[ 1 ]] ) } )
		colnames ( seqs ) <- NULL
		seqs <- t ( seqs )
		list ( positions [ - 1 ] , seqs )
	})
	#freq.specs <- rowSums ( freq.specs )
	return ( my.seqs )
}


## needs work
CountHaps <- function ( these.seqs , len.bp , hap.count.interval ) {
	#recover ()
	positions <- as.numeric ( these.seqs [[ 1 ]] )
	my.part <- list ( seq ( 1 : nrow ( these.seqs [[ 2 ]] ) ) )
	hap.freqs <- matrix ( 0 , nrow = nrow ( these.seqs [[ 2 ]] ) , ncol = len.bp / hap.count.interval + 1 )
	pos.cuts <- seq ( 0 , 1 , by = hap.count.interval / len.bp )
	hap.freqs [ 1 , pos.cuts < positions [ 1 ] ] <- nrow ( these.seqs [[ 2 ]] )
	for ( i in 1 : length ( positions ) ) {
		#if ( i == 18 ) break
		this.site <- which ( these.seqs [[ 2 ]] [ , i ] == 1 )
		if ( any ( unlist ( lapply ( my.part ,  function ( x ) all ( this.site %in% x ) ) ) ) ) { 
			if ( i < length ( positions ) ) {
				hap.freqs [ 1 : length ( my.part ) , pos.cuts >= positions [ i ] & pos.cuts < positions [ i + 1 ] ] <- sort ( unlist ( lapply ( my.part , length ) ) , d = T )
			}
			if ( i == length ( positions ) ) {
				hap.freqs [ 1 : length ( my.part ) , pos.cuts >= positions [ i ] ] <- sort ( unlist ( lapply ( my.part , length ) ) , d = T )
			}
			next
		} else {
			new.part <- list ()
			for ( j in 1 : length ( my.part ) ) {
				x <- my.part [[ j ]]
				if ( length ( x [ x %in% this.site ] ) > 0 ) {
					new.part [[ length ( new.part ) + 1 ]] <- x [ x %in% this.site ]
					new.part [[ length ( new.part ) + 1 ]] <- x [ !(x %in% this.site) ]
				} else {
					new.part [[ length ( new.part ) + 1 ]] <- x
				}
			}
			my.part <- new.part
		}
		if ( i < length ( positions ) ) {
			hap.freqs [ 1 : length ( my.part ) , pos.cuts >= positions [ i ] & pos.cuts < positions [ i + 1 ] ] <- sort ( unlist ( lapply ( my.part , length ) ) , d = T )
		} 
		
		if ( i == length ( positions ) ) {
			hap.freqs [ 1 : length ( my.part ) , pos.cuts >= positions [ i ] ] <- sort ( unlist ( lapply ( my.part , length ) ) , d = T )
		}
		
		if ( length ( my.part ) == nrow ( hap.freqs ) ) {
			hap.freqs [ , pos.cuts >= positions [ i ] ] <- 1
			break
		}
	}
	return ( hap.freqs )
}



CountHapsNoSing <- function ( these.seqs , len.bp , hap.count.interval ) {
	#recover ()
	positions <- as.numeric ( these.seqs [[ 1 ]] )
	my.part <- list ( seq ( 1 : nrow ( these.seqs [[ 2 ]] ) ) )
	hap.freqs <- matrix ( 0 , nrow = nrow ( these.seqs [[ 2 ]] ) , ncol = len.bp / hap.count.interval + 1 )
	pos.cuts <- seq ( 0 , 1 , by = hap.count.interval / len.bp )
	hap.freqs [ 1 , pos.cuts < positions [ 1 ] ] <- nrow ( these.seqs [[ 2 ]] )
	for ( i in 1 : length ( positions ) ) {
		#if ( i == 18 ) break
		this.site <- which ( these.seqs [[ 2 ]] [ , i ] == 1 )
		if ( i == 94 ) break
		if ( any ( unlist ( lapply ( my.part ,  function ( x ) all ( x %in% this.site ) & all ( this.site %in% x ) ) ) ) |  length ( this.site ) == 1 ) { 
			hap.freqs [ 1 : length ( my.part ) , pos.cuts >= positions [ i ] & pos.cuts < positions [ i + 1 ] ] <- sort ( unlist ( lapply ( my.part , length ) ) , d = T )
			next
		} else {
			new.part <- list ()
			for ( j in 1 : length ( my.part ) ) {
				x <- my.part [[ j ]]
				if ( length ( x [ x %in% this.site ] ) > 0 ) {
					new.part [[ length ( new.part ) + 1 ]] <- x [ x %in% this.site ]
					new.part [[ length ( new.part ) + 1 ]] <- x [ !(x %in% this.site) ]
				} else {
					new.part [[ length ( new.part ) + 1 ]] <- x
				}
			}
			my.part <- new.part
		}
		#these.hap.freqs <- numeric ( length ( unlist ( my.part ) ) )
		
		if ( i == length ( positions ) ) {
			hap.freqs [ 1 : length ( my.part ) , pos.cuts >= positions [ i ] ] <- sort ( unlist ( lapply ( my.part , length ) ) , d = T )
		} else {
			hap.freqs [ 1 : length ( my.part ) , pos.cuts >= positions [ i ] & pos.cuts < positions [ i + 1 ] ] <- sort ( unlist ( lapply ( my.part , length ) ) , d = T )
		}
		if ( length ( my.part ) == nrow ( hap.freqs ) ) {
			hap.freqs [ , pos.cuts >= positions [ i ] ] <- 1
			break
		}
	}
	return ( hap.freqs )
}






TwoSideCountHaps <- function ( these.seqs , len.bp , hap.count.interval ) {
	#recover ()
	positions <- as.numeric ( these.seqs [[ 1 ]] )
	c.positions <- positions - 1/2
	my.part <- list ( seq ( 1 : nrow ( these.seqs [[ 2 ]] ) ) )
	hap.freqs <- matrix ( 0 , nrow = nrow ( these.seqs [[ 2 ]] ) , ncol = len.bp / hap.count.interval + 1 )
	pos.cuts <- seq ( 0 , 0.5 , by = hap.count.interval / (2*len.bp) )
	pb <- txtProgressBar ( min = 0 , max = ncol ( hap.freqs ) , style = 3 )
	for ( i in seq_along ( pos.cuts ) ) {
		my.window <- c.positions > -pos.cuts [ i ] & c.positions < pos.cuts [ i ]
		
		if ( all ( my.window == FALSE ) ) {
			
			hap.freqs [ 1 , i ]	<- nrow ( these.seqs [[ 2 ]] )
		
		} else if ( sum ( my.window ) == 1 ){
			
			freqs <- as.numeric ( sort ( table ( these.seqs [[ 2 ]] [ , my.window ] ) , decreasing = T ) )
			hap.freqs [ seq_along ( freqs ) , i ] <- freqs
			
		} else {
			#break
			blah <- apply ( these.seqs [[ 2 ]] [ , my.window ] , 1 , paste , collapse = " " )
			freqs <- sort ( rle ( sort ( as.numeric ( factor ( blah , unique ( blah ) , ordered = T ) ) ) )$lengths , decreasing = T )
			hap.freqs [ seq_along ( freqs ) , i ] <- freqs
	
		}
		setTxtProgressBar ( pb, i )
	}
	
	return ( hap.freqs )
}




hard.runs <- SweepFromStandingSim ( N = 10000 , s = 0.01 , f = 1/20000 , reps = 1000 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE  , display.rep.count = TRUE , time.factor = 1  )
hard.sweep <- msHapSims ( hard.runs [[ 1 ]] , n.sam = 100 , f = 1/20000 , s = 0.01 , N = 10000 , path = "Sims/HapSims" , num.sims = 1 , len.bp = 2000000 , r.bp = 10^-8 , mu.bp = 10^-8 , ext = "hapSims" , hap.count.interval = 5000 , both.side = T )
save ( hard.sweep , file = "Sims/HapSims/both.sides.n100.f05.Robj" )

standing.runs <- SweepFromStandingSim ( N = 10000 , s = 0.01 , f = 0.05 , reps = 1000 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE  , display.rep.count = TRUE , time.factor = 1  )
standing.sweep <- msHapSims ( standing.runs [[ 1 ]] , n.sam = 100 , f = 0.05 , s = 0.01 , N = 10000 , path = "Sims/HapSims" , num.sims = 1 , len.bp = 2000000 , r.bp = 10^-8 , mu.bp = 10^-8 , hap.count.interval = 5000 , both.sides = TRUE )
save ( standing.sweep , file = "Sims/HapSims/both.sides.n100.f05.Robj" )





