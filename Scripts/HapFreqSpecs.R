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
		message ( i )
		if ( both.sides == FALSE ) {
			system ( paste ( "Scripts/msseldir/mssel " , n.sam , " " , num.sims , " 0 " , n.sam ,  " " , traj.file ,  " 0 -t " , 2 * N * len.bp * mu.bp , " -r " , 2 * N * len.bp * r.bp , " " , len.bp , " > " , path, "/myseqdata" , sep = "" ) ) 
			seqs <- GetSeqs ( n.sam , num.sims , path )
			hap.counts [ ( i - 1 ) * num.sims + 1:num.sims ] <- lapply ( seqs , CountHaps , len.bp , hap.count.interval )
			#hap.counts.no.sing [ ( i - 1 ) * num.sims + 1:num.sims ] <- lapply ( seqs , CountHapsNoSing , len.bp , hap.count.interval )
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


GetSeqsMultMut <- function ( n , num.sims , path , muts ) {
	#recover()
	a<-system(paste("grep segsites ", path , "/myseqdata.mult.mut",sep=""),intern=TRUE)
	seg.sites<-sapply(a,function(b){as.numeric(strsplit(b,":")[[1]][2])})
	polymorph<- seg.sites>0
	seq.lines<-c(0,cumsum(polymorph*n)[-length(polymorph)])	
	my.seqs <- list ()
	for ( iter in 0 : ( num.sims - 1 ) ) {
		positions <- read.table ( paste ( path , "/myseqdata.mult.mut" ,sep = "" ) , skip = 5 + 5 * iter + seq.lines [ iter + 1 ] , nrow = 1 )
		seqs.raw <- scan ( paste ( path , "/myseqdata.mult.mut" , sep = "" ) , skip = 6 + 5 * iter + seq.lines [ iter + 1 ] , nline = ( n + 2 ) , what = character ( ) , quiet = TRUE )
		n.muts <- strsplit ( seqs.raw [ length ( seqs.raw ) ] , ":" )[[1]][2]
		if ( n.muts != muts & n.muts != ( muts + 1 ) ) {
			next
		}
		seqs <- sapply ( seqs.raw [ - length ( seqs.raw ) ] , function ( seq ) { as.numeric ( strsplit ( seq , "" ) [[ 1 ]] ) } )
		colnames ( seqs ) <- NULL
		seqs <- t ( seqs )
		
		my.seqs [[ length ( my.seqs ) + 1 ]] <- list ( positions [ - 1 ] , seqs , n.muts )
	}
	#freq.specs <- rowSums ( freq.specs )
	return ( my.seqs )
}






CountHaps <- function ( these.seqs , len.bp , hap.count.interval ) {
	#recover ()
	positions <- as.numeric ( these.seqs [[ 1 ]] )
	hap.freqs <- matrix ( 0 , nrow = nrow ( these.seqs [[ 2 ]] ) , ncol = len.bp / hap.count.interval + 1 )
	pos.cuts <- seq ( 0 , 1 , by = hap.count.interval / len.bp )
	if ( positions [ 1 ] == 0 ) {
		my.part <- sapply ( unique ( these.seqs [[ 2 ]] [ , 1 ] ) , function ( x ) which ( these.seqs [[ 2 ]] [ , 1 ] == x ) )	
		hap.freqs [  seq_along ( unique ( these.seqs [[ 2 ]] [ , 1 ] ) ) , 1 ] <- sort ( sapply ( my.part , length ) , d = T )	
	} else {
		my.part <- list ( seq ( 1 : nrow ( these.seqs [[ 2 ]] ) ) )
		hap.freqs [ 1 , pos.cuts < positions [ 1 ] ] <- nrow ( these.seqs [[ 2 ]] )
	}
	pb <- txtProgressBar ( min = 0 , max = length ( positions ) , style = 3 )
	for ( i in 1 : length ( positions ) ) {
		setTxtProgressBar ( pb, i )
		this.site <- which ( these.seqs [[ 2 ]] [ , i ] == 1 )	
		#lapply ( my.part , function ( x ) all ( sort ( x ) == sort ( this.site ) ) )
		#lapply ( my.part , function ( x ) length ( x ) == length ( this.site ) )
		
		if ( any ( unlist ( lapply ( my.part , function ( x ) all ( sort ( x ) == sort ( this.site ) ) ) ) ) ) { 
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
				
				if ( all ( my.part [[ j ]] %in% this.site ) | length ( x [ x %in% this.site ] ) == 0 ) {
					new.part [[ length ( new.part ) + 1 ]] <- x
				} else {
					new.part [[ length ( new.part ) + 1 ]] <- x [ x %in% this.site ]
					new.part [[ length ( new.part ) + 1 ]] <- x [ !(x %in% this.site) ]
				}
			}
			my.part <- new.part [ lapply ( my.part , length ) != 0 ]
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
	cat ( "\n" )
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


HapFreqs <- function ( hap.counts ) {
	
	sum.counts <- Reduce ( "+" , hap.counts )
	expect.freqs <- sum.counts / length ( hap.counts )
	exist.counts <- Reduce ( "+" , lapply ( hap.counts , function ( x ) x > 0 ) )
	expect.freqs.cond.exist <- sum.counts / exist.counts
	prob.exist <- exist.counts / length ( hap.counts )
	
	return ( list ( expect.freqs , expect.freqs.cond.exist , prob.exist ) )
	
}


### multiple mutations
load ( "Sims/HapSims/one.side.soft.n100.k3.s01.manywindows.Robj")
while ( length ( blah ) < 5000 ) {
	message ( length ( blah ) )
	system ( paste ( "Scripts/msms/bin/msms -ms 100 5 -N 10000 -t 200 -r 200 500000 -SAA 400 -SAa 200 -Smu 0.4 -Sp 0 -oOC -SF 0 -Smark  > Sims/myseqdata.mult.mut" , sep = "" ) )
	seqs <- GetSeqsMultMut ( 100 , 5 , "Sims" ,  3 )
	if ( length ( seqs ) == 0 ) next
	keep <- which ( unlist ( lapply ( seqs , function ( x ) length ( unique ( x[[2]][,1] ) ) == 3 ) ) )
	seqs <- seqs [ keep ]
	if ( length ( seqs ) == 0 ) next
	tmp <- lapply ( seqs , CountHaps , 500000 , 100 )
	blah [ length ( blah ) + 1 : length ( tmp ) ] <- tmp
	if ( length ( blah ) %% 1000 == 0 )
	save ( blah , file = "Sims/HapSims/one.side.soft.n100.k3.s01.manywindows.Robj" )
}



## one side
hard.runs <- SweepFromStandingSim ( N = 10000 , s = 0.01 , f = 1/20000 , reps = 2000 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE  , display.rep.count = TRUE , time.factor = 1  )
hard.sweep.n100.N10000.s01 <- msHapSims ( hard.runs [[ 1 ]] , n.sam = 100 , f = 1/20000 , s = 0.01 , N = 10000 , path = "Sims/HapSims" , num.sims = 5 , len.bp = 500000 , r.bp = 10^-8 , mu.bp = 10^-8 , ext = "hapSims" , hap.count.interval = 100 , both.side = F )
save ( hard.sweep.n100.N10000.s01 , file = "Sims/HapSims/one.side.hard.n100.denovo.s01.manywindows.Robj" )


hard.runs <- SweepFromStandingSim ( N = 10000 , s = 0.006 , f = 1/20000 , reps = 2000 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE  , display.rep.count = TRUE , time.factor = 1  )
hard.sweep.n100.N10000.s006 <- msHapSims ( hard.runs [[ 1 ]] , n.sam = 100 , f = 1/20000 , s = 0.006 , N = 10000 , path = "Sims/HapSims" , num.sims = 5 , len.bp = 500000 , r.bp = 10^-8 , mu.bp = 10^-8 , ext = "hapSims" , hap.count.interval = 100 , both.side = F )
save ( hard.sweep.n100.N10000.s006 , file = "Sims/HapSims/one.side.hard.n100.denovo.s006.manywindows.Robj" )



standing.runs <- SweepFromStandingSim ( N = 10000 , s = 0.01 , f = 0.05 , reps = 2000 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE  , display.rep.count = TRUE , time.factor = 1  )
standing.sweep.f05.n100.N10000.s01 <- msHapSims ( standing.runs [[ 1 ]] , n.sam = 100 , f = 0.05 , s = 0.01 , N = 10000 , path = "Sims/HapSims" , num.sims = 5 , len.bp = 500000 , r.bp = 10^-8 , mu.bp = 10^-8 , hap.count.interval = 100 , both.sides = F )
save ( standing.sweep.f05.n100.N10000.s01 , file = "Sims/HapSims/one.side.standing.n100.f05.s01.manywindows.Robj" )




if ( FALSE ) {
	
	
	
blah <- list ()
while ( length ( blah ) < 1000 ) {
	system ( paste ( "Sims/msms/bin/msms -ms 100 10 -N 10000 -t 200 -r 200 500000 -SAA 2000 -SAa 1000 -Smu 0.4 -Sp 0 -oOC -SF 0 -Smark  > Sims/myseqdata.mult.mut" , sep = "" ) )
	seqs <- GetSeqsMultMut ( 100 , 10 , "Sims" ,  3 )
	if ( length ( seqs ) == 0 ) next
	keep <- which ( unlist ( lapply ( seqs , function ( x ) length ( unique ( x[[2]][,1] ) ) == 3 ) ) )
	seqs <- seqs [ keep ]
	if ( length ( seqs ) == 0 ) next
	tmp <- lapply ( seqs , CountHaps , 500000 , 500 )
	blah [ length ( blah ) + 1 : length ( tmp ) ] <- tmp
	message ( length ( blah ) )
	save ( blah , file = "Sims/HapSims/one.side.soft.n100.k3.s05.Robj" )
}	
	
	

### neutral sims
blah <- list ()
for ( i in 1 : 100 ) {
	message ( i )
	system ( paste ( "Scripts/msdir/ms 100 10 -t 200 -r 200 500000 > Sims/myseqdata" , sep = "" ) )
	seqs <- GetSeqs ( 100 , 10 , "Sims" )
	blah [ ( i - 1 )*10 + ( 1 : 10 ) ] <- lapply ( seqs , CountHaps , 500000 , 5000 )
	save ( blah , file = "Sims/HapSims/neutral.n100.Robj" )
}
neutral <- list ()
neutral [[ 1 ]] <- Reduce ( "+" , blah ) / length ( blah )
neutral [[ 2 ]] <- blah
save ( neutral , file = "Sims/HapSims/neutral.n100.Robj" )





load ( "Sims/HapSims/one.side.hard.n100.denovo.s01.Robj" )
load ( "Sims/HapSims/one.side.standing.n100.f05.s01.Robj" )
load ( "Sims/HapSims/one.side.soft.n100.k3.s01.Robj" )
soft.sweep <- blah
load ( "Sims/HapSims/neutral.n100.Robj" )
## neutral <- Reduce ( "+" , blah ) / length ( blah)


soft.haps <- HapFreqs ( soft.sweep )
standing.haps <- HapFreqs ( standing.sweep [[ 2 ]] )
hard.haps <- HapFreqs ( hard.sweep [[ 2 ]] )
neutral.haps <- HapFreqs ( neutral [[ 2 ]] )


pdf ( "Figures/HapFreqsExpectAllThree.pdf" , width = 10 , height = 6 )
matplot ( t ( hard.haps [[ 1 ]] ) [ , 1:9 ] , type = "l" , lwd = 2 , lty = 1 , ylim = c ( 0 , 40 ) , col = brewer.pal( 9 , "Set1" ) , )
matplot ( t ( standing.haps [[ 1 ]] ) [ , 1:9 ] , type = "l" , lwd = 2 , lty = 2 , ylim = c ( 0 , 80 ) , col = brewer.pal( 9 , "Set1" ) , add = T )
matplot ( t ( soft.haps [[ 1 ]] ) [ , 1:9 ] , type = "l" , lwd = 2 , lty = 3 , ylim = c ( 0 , 80 ) , col = brewer.pal( 9 , "Set1" ) , add = T )
dev.off()


pdf ( "Figures/HapFreqsExpectCondExistAllThree.pdf" , width = 10 , height = 6 )
matplot ( t ( hard.haps [[ 2 ]] ) [ , 1:9 ] , type = "l" , lwd = 2 , lty = 1 , ylim = c ( 0 , 20 ) , col = brewer.pal( 9 , "Set1" ) , )
matplot ( t ( standing.haps [[ 2 ]] ) [ , 1:9 ] , type = "l" , lwd = 2 , lty = 2 , ylim = c ( 0 , 20 ) , col = brewer.pal( 9 , "Set1" ) , add = T )
matplot ( t ( soft.haps [[ 2 ]] ) [ , 1:9 ] , type = "l" , lwd = 2 , lty = 3 , ylim = c ( 0 , 20 ) , col = brewer.pal( 9 , "Set1" ) , add = T )
dev.off()


pdf ( "Figures/HapFreqsExistProbAllThree.pdf" , width = 10 , height = 6 )
matplot ( t ( hard.haps [[ 3 ]] ) [ , 1:9 ] , type = "l" , lwd = 2 , lty = 1 , ylim = c ( 0 , 1 ) , col = brewer.pal( 9 , "Set1" ) , xlim = c ( 0 , 20 ) )
matplot ( t ( standing.haps [[ 3 ]] ) [ , 1:9 ] , type = "l" , lwd = 2 , lty = 2 , ylim = c ( 0 , 1 ) , col = brewer.pal( 9 , "Set1" ) , add = T )
matplot ( t ( soft.haps [[ 3 ]] ) [ , 1:9 ] , type = "l" , lwd = 2 , lty = 3 , ylim = c ( 0 , 1 ) , col = brewer.pal( 9 , "Set1" ) , add = T )
dev.off()


pdf ( "Figures/HapFreqRatiosCondExist.pdf" , height = 10 , width = 8 )
par ( mfrow = c ( 3,2))
image ( t ( apply ( standing.haps [[ 2 ]] / hard.haps [[ 2 ]] , 2 , rev) ) ,  col = c ( "black" , "black" ) , breaks = seq ( 0.95 , 1.05 ,length.out = 3) ,xaxt = "n" , main = expression ( h[i]^stand/h[i]^hard ) , yaxt = "n" , ylab = expression ( h[i]))
axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
axis ( 2 , c ( seq ( 0 , 80 , length.out = 5 ) , 99 )/100 , labels = c ( seq ( 100 , 20 , length.out = 5 ) , 1 ) )
image ( t ( apply ( standing.haps [[ 2 ]] / hard.haps [[ 2 ]] , 2 , rev) ) ,  col = heat.colors ( 100 ) , breaks = seq ( 1.05 , 2 ,length.out = 101) , add = T )
image ( t ( apply ( standing.haps [[ 2 ]] / hard.haps [[ 2 ]] , 2 , rev) ) ,  col = cm.colors ( 100 ) , breaks = seq ( 0 , 0.95 ,length.out = 101) , add = T )


image ( t ( apply ( standing.haps [[ 2 ]] / neutral.haps [[ 2 ]] , 2 , rev) ) , breaks = seq ( 0.95, 1.05 ,length.out=2), col = "black" , xaxt = "n" , main = expression ( h[i]^stand/h[i]^neut ) ,yaxt = "n", ylab = expression ( h[i]))
axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
axis ( 2 , c ( seq ( 0 , 80 , length.out = 5 ) , 99 )/100 , labels = c ( seq ( 100 , 20 , length.out = 5 ) , 1 ) )
image ( t ( apply ( standing.haps [[ 2 ]] / neutral.haps [[ 2 ]] , 2 , rev) ) , breaks = seq ( 1.05 , 2 ,length.out=101), col = heat.colors ( 100 ),add = T )
image ( t ( apply ( standing.haps [[ 2 ]] / neutral.haps [[ 2 ]] , 2 , rev) ) , breaks = seq ( 0 , 0.95 ,length.out=101), col = cm.colors ( 100 ), add =T )


image ( t ( apply ( soft.haps [[ 2 ]] / hard.haps [[ 2 ]] , 2 , rev) ) ,  col = c ( "black" , "black" ) , breaks = seq ( 0.95 , 1.05 ,length.out = 3) , xaxt = "n", main = expression ( h[i]^soft/h[i]^hard ) ,yaxt = "n" , ylab = expression ( h[i]))
axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
axis ( 2 , c ( seq ( 0 , 80 , length.out = 5 ) , 99 )/100 , labels = c ( seq ( 100 , 20 , length.out = 5 ) , 1 ) )
image ( t ( apply ( soft.haps [[ 2 ]] / hard.haps [[ 2 ]] , 2 , rev) ) ,  col = heat.colors ( 100 ) , breaks = seq ( 1.05 , 2 ,length.out = 101) , add = T )
image ( t ( apply ( soft.haps [[ 2 ]] / hard.haps [[ 2 ]] , 2 , rev) ) ,  col = cm.colors ( 100 ) , breaks = seq ( 0 , 0.95 ,length.out = 101) , add = T )

image ( t ( apply ( soft.haps [[ 2 ]] / neutral.haps [[ 2 ]] , 2 , rev) ) , breaks = seq ( 0.95, 1.05 ,length.out=2), col = "black" , xaxt = "n", main = expression ( h[i]^soft/h[i]^neut ) ,yaxt = "n" , ylab = expression ( h[i]))
axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
axis ( 2 , c ( seq ( 0 , 80 , length.out = 5 ) , 99 )/100 , labels = c ( seq ( 100 , 20 , length.out = 5 ) , 1 ) )
image ( t ( apply ( soft.haps [[ 2 ]] / neutral.haps [[ 2 ]] , 2 , rev) ) , breaks = seq ( 1.05 , 2 ,length.out=101), col = heat.colors ( 100 ),add = T )
image ( t ( apply ( soft.haps [[ 2 ]] / neutral.haps [[ 2 ]] , 2 , rev) ) , breaks = seq ( 0 , 0.95 ,length.out=101), col = cm.colors ( 100 ), add =T )

#plot ( c ( 0,1), type = "n",bty = "n" ,xaxt ="n",yaxt ="n",xlab = "",ylab ="")

image ( t ( apply ( standing.haps [[ 2 ]] / soft.haps [[ 2 ]] , 2 , rev) ) ,  col = c ( "black" , "black" ) , breaks = seq ( 0.95 , 1.05 ,length.out = 3) , xaxt = "n", main = expression ( h[i]^stand/h[i]^soft ) , yaxt = "n" , ylab = expression ( h[i]))
axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
axis ( 2 , c ( seq ( 0 , 80 , length.out = 5 ) , 99 )/100 , labels = c ( seq ( 100 , 20 , length.out = 5 ) , 1 ) )
image ( t ( apply ( standing.haps [[ 2 ]] / soft.haps [[ 2 ]] , 2 , rev) ) ,  col = heat.colors ( 100 ) , breaks = seq ( 1.05 , 2 ,length.out = 101) , add = T )
image ( t ( apply ( standing.haps [[ 2 ]] / soft.haps [[ 2 ]] , 2 , rev) ) ,  col = cm.colors ( 100 ) , breaks = seq ( 0 , 0.95 ,length.out = 101) , add = T )

image ( t ( apply ( hard.haps [[ 2 ]] / neutral.haps [[ 2 ]] , 2 , rev) ) , breaks = seq ( 0.95, 1.05 ,length.out=2), col = "black" , xaxt = "n" , main = expression ( h[i]^hard/h[i]^neut ) , yaxt = "n" , ylab = expression ( h[i]))
axis ( 1 , seq ( 0 , 1, length.out = 5 ) , seq ( 0 , 0.005 , length.out = 5 ))
axis ( 2 , c ( seq ( 0 , 80 , length.out = 5 ) , 99 )/100 , labels = c ( seq ( 100 , 20 , length.out = 5 ) , 1 ) )
image ( t ( apply ( hard.haps [[ 2 ]] / neutral.haps [[ 2 ]] , 2 , rev) ) , breaks = seq ( 1.05 , 2 ,length.out=101), col = heat.colors ( 100 ),add = T )
image ( t ( apply ( hard.haps [[ 2 ]] / neutral.haps [[ 2 ]] , 2 , rev) ) , breaks = seq ( 0 , 0.95 ,length.out=101), col = cm.colors ( 100 ), add =T )
dev.off()







image ( t ( apply ( standing.haps [[ 1 ]] / hard.haps [[ 1 ]] , 2 , rev) ) , col = heat.colors ( 100 ) , breaks = seq ( 0,5,length.out=101) )
my.range <- range ( standing.haps [[ 2 ]] / hard.haps [[ 2 ]] ,na.rm = T )
image ( t ( apply ( standing.haps [[ 2 ]] / hard.haps [[ 2 ]] , 2 , rev) ) ,  col = heat.colors ( 1000 ), breaks = seq ( 0.24 , 3.5 ,length.out=1001))
image ( t ( apply ( standing.haps [[ 3 ]] / hard.haps [[ 3 ]] , 2 , rev) ) , breaks = seq ( 0,5,length.out=101) , col = heat.colors ( 100 ))








image ( t ( apply ( standing.haps [[ 1 ]] / neutral.haps [[ 1 ]] , 2 , rev) ) , breaks = seq ( 0,5,length.out=101) ,  col = heat.colors ( 100 ) )
my.range <- range ( standing.haps [[ 2 ]] / neutral.haps [[ 2 ]] ,na.rm = T )
image ( t ( apply ( standing.haps [[ 2 ]] / neutral.haps [[ 2 ]] , 2 , rev) ) , breaks = seq ( my.range [1 ] , 4 ,length.out=101), col = heat.colors ( 100 ))
image ( t ( apply ( standing.haps [[ 3 ]] / neutral.haps [[ 3 ]] , 2 , rev) ) , breaks = seq ( 0,1,length.out=101), col = heat.colors ( 100 ))






image ( t ( apply ( hard.haps [[ 1 ]] / neutral.haps [[ 1 ]] , 2 , rev) ) , breaks = seq ( 0,5,length.out=101) ,  col = heat.colors ( 100 ) )
my.range <- range ( hard.haps [[ 2 ]] / neutral.haps [[ 2 ]] ,na.rm = T )
image ( t ( apply ( hard.haps [[ 2 ]] / neutral.haps [[ 2 ]] , 2 , rev) ) , breaks = seq ( my.range [1 ] , 4 ,length.out=101), col = heat.colors ( 100 ))
image ( t ( apply ( hard.haps [[ 3 ]] / neutral.haps [[ 3 ]] , 2 , rev) ) , breaks = seq ( 0,1,length.out=101), col = heat.colors ( 100 ))





image ( t ( apply ( soft.haps [[ 1 ]] / hard.haps [[ 1 ]] , 2 , rev) ) , breaks = seq ( 0,4,length.out=13) )
my.range <- range ( soft.haps [[ 2 ]] / hard.haps [[ 2 ]] ,na.rm = T )
image ( t ( apply ( soft.haps [[ 2 ]] / hard.haps [[ 2 ]] , 2 , rev) ) , breaks = seq ( 0.24 , 3.5 ,length.out=1001),col=heat.colors(1000))
image ( t ( apply ( soft.haps [[ 3 ]] / hard.haps [[ 3 ]] , 2 , rev) ) , breaks = seq ( 0,4,length.out=13))









pdf ( "Figures/HapFreqsExpect.pdf" , width = 10 , height = 6 )
par ( mfrow = c ( 2 , 2 ) )
matplot ( t ( hard.haps [[ 1 ]] ) [ , 1:50 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 20 ) )
matplot (  t ( standing.haps [[ 1 ]] ) [ , 1 : 50 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 20 ) )
matplot ( t ( soft.haps [[ 1 ]] ) [ , 1:50 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 20 ) )
matplot ( t ( neutral.haps [[ 1 ]] ) [ , 1:50 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 20 ) )
dev.off()

pdf ( "Figures/HapFreqsExpectCondExist.pdf" , width = 10 , height = 6 )
par ( mfrow = c ( 2 , 2 ) )
matplot ( t ( hard.haps [[ 2 ]] ) [ , 1:50 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 8 ) )
matplot (  t ( standing.haps [[ 2 ]] ) [ , 1 : 50 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 8 ) )
matplot ( t ( soft.haps [[ 2 ]] ) [ , 1:50 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 8 ) )
matplot ( t ( neutral.haps [[ 2 ]] ) [ , 1:50 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 8 ) )
dev.off()

pdf ( "Figures/HapFreqsExistProbs.pdf" , width = 10 , height = 6 )
par ( mfrow = c ( 2 , 2 ) )
matplot ( t ( hard.haps [[ 3 ]] ) [ , 1:50 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 8 ) )
matplot (  t ( standing.haps [[ 3 ]] ) [ , 1 : 50 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 8 ) )
matplot ( t ( soft.haps [[ 3 ]] ) [ , 1:50 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 8 ) )
matplot ( t ( neutral.haps [[ 3 ]] ) [ , 1:50 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 8 ) )
dev.off()




pdf ( "Figures/HapFreqRatiosUncond.pdf" , width = 10 , height = 6 )
par ( mfrow = c ( 2 , 3 ) )
matplot (  t ( standing.haps [[ 1 ]] / hard.haps [[ 1 ]] ) [ , 1 : 80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( standing.haps [[ 1 ]] / soft.haps [[ 1 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( standing.haps [[ 1 ]] / soft.haps [[ 1 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( standing.haps [[ 1 ]] / neutral.haps [[ 1 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( hard.haps [[ 1 ]] / neutral.haps [[ 1 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( soft.haps [[ 1 ]] / neutral.haps [[ 1 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
dev.off()

pdf ( "Figures/HapFreqRatiosCondExistence.pdf" , width = 10 , height = 6 )
par ( mfrow = c ( 2 , 3 ) )
matplot (  t ( standing.haps [[ 2 ]] / hard.haps [[ 2 ]] ) [ , 1 : 80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( standing.haps [[ 2 ]] / soft.haps [[ 2 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( hard.haps [[ 2 ]] / soft.haps [[ 2 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( standing.haps [[ 2 ]] / neutral.haps [[ 2 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( hard.haps [[ 2 ]] / neutral.haps [[ 2 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( soft.haps [[ 2 ]] / neutral.haps [[ 2 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
dev.off()

pdf ( "Figures/HapExistenceRatios.pdf" , width = 10 , height = 6 )
par ( mfrow = c ( 2 , 3 ) )
matplot (  t ( standing.haps [[ 3 ]] / hard.haps [[ 3 ]] ) [ , 1 : 80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( standing.haps [[ 3 ]] / soft.haps [[ 3 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( hard.haps [[ 3 ]] / soft.haps [[ 3 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( standing.haps [[ 3 ]] / neutral.haps [[ 3 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( hard.haps [[ 3 ]] / neutral.haps [[ 3 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
matplot ( t ( soft.haps [[ 3 ]] / neutral.haps [[ 3 ]] ) [ , 1:80 ] , type = "l" , lwd = 1 , lty = 1 , ylim = c ( 0 , 5 ) )
dev.off()






par ( mfrow = c ( 1 , 3 ) )
matplot (  t ( hard.sweep [[ 1 ]] ) [ , 1 : 20 ] , type = "l" , lwd = 1 , lty = 1  )
matplot (  t ( standing.sweep [[ 1 ]] ) [ , 1 : 20 ] , type = "l" , lwd = 1 , lty = 1 )
matplot (  t ( neutral ) [ , 1 : 20 ] , type = "l" , lwd = 1 , lty = 1 )




## both sides
hard.runs <- SweepFromStandingSim ( N = 10000 , s = 0.01 , f = 1/20000 , reps = 1000 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE  , display.rep.count = TRUE , time.factor = 1  )
hard.sweep <- msHapSims ( hard.runs [[ 1 ]] , n.sam = 100 , f = 1/20000 , s = 0.01 , N = 10000 , path = "Sims/HapSims" , num.sims = 10 , len.bp = 2000000 , r.bp = 10^-8 , mu.bp = 10^-8 , ext = "hapSims" , hap.count.interval = 5000 , both.side = T )
save ( hard.sweep , file = "Sims/HapSims/both.sides.hard.n100.denovo.s01.Robj" )

standing.runs <- SweepFromStandingSim ( N = 10000 , s = 0.01 , f = 0.05 , reps = 1000 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE  , display.rep.count = TRUE , time.factor = 1  )
standing.sweep <- msHapSims ( standing.runs [[ 1 ]] , n.sam = 100 , f = 0.05 , s = 0.01 , N = 10000 , path = "Sims/HapSims" , num.sims = 10 , len.bp = 2000000 , r.bp = 10^-8 , mu.bp = 10^-8 , hap.count.interval = 5000 , both.sides = TRUE )
save ( standing.sweep , file = "Sims/HapSims/both.sides.standing.n100.f05.s01.Robj" )


}


