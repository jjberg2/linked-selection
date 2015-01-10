source ( "Scripts/SweepFromStandingSim.R")


msHapSims <- function ( runs , n.sam = 2  , f , s , N , path , get.site.density = TRUE , recom , num.sims , len.bp , r.bp , mu.bp , hap.count.interval ) {
	recover()
	options ( "scipen" = 100 , "digits" = 4 )
	f.lab <- strsplit ( as.character ( f ) , "\\." ) [[ 1 ]] [ 2 ]
	s.lab <- strsplit ( as.character ( s ) , "\\." ) [[ 1 ]][ 2 ]
	counter <- 1
	
	my.file <- paste ( path , "/mssel_f" , n.sam , "." ,  f.lab , "."  , s.lab , "." , N  , ".out" , sep = "" )
	num.sims<-20
	system ( paste ( "rm " , my.file ) )
	#for ( run in 1:5 ) {
	#	load ( paste ( "run_cond_lost_" , run , ".Robj" , sep = "" ) )
	
	for ( i in 1:nrow ( runs ) ) {

		my.freqs <- runs [ i , runs [ i , ] > 0 ]
		my.times <- 0 : length ( my.freqs )
		my.freqs <- c ( my.freqs , 0 )
		
		my.times <- my.times / ( 4*N  )

		#recover()
		header.material <- c ( "1" , "1" , paste ( "n:" , length ( my.times ) ) )
		traj.file <- paste ( path , "/my.standing" , "." , f.lab , "." , s.lab , "." , N , ".traj" , sep = "" )
		write ( file = traj.file , header.material )
		write.table ( file = traj.file , cbind ( my.times , my.freqs ) , append = TRUE , sep = "\t" , quot = FALSE , col.nam = FALSE , row.name = FALSE )
		cat( i ," " )
		# if ( get.site.density ) { 
			# system ( paste ( path , "Scripts/msseldir/mssel " , n.sam , " 20 0 " , n.sam , " " , path , "Sims/my.standing" , "." , f.lab , "." , s.lab , "." , N, ".traj 0 -t 200. -r 200. 20000 | grep pos | cut -f 2 -d : >> " , my.file , sep = "" ) )
		# }	else	{   ##setup for the mo. to do freq. spectrum
		system ( paste ( "Scripts/msseldir/mssel " , n.sam , " " , num.sims , " 0 " , n.sam ,  " " , traj.file ,  " 0 -t " , 2 * N * len.bp * mu.bp , " -r " , 2 * N * len.bp * r.bp , " " , len.bp , " > " , path, "/myseqdata" , sep = "" ) ) 
			
			seqs <- GetSeqs ( n.sam , num.sims , path )
			lapply ( seqs , CountHaps )
			
			
			my.specs[,(1+(counter-1)*num.sims):(counter*num.sims)]<-spec
			counter<-counter+1
			
			#recover()
		# }
	}
	if (! get.site.density ) return(my.specs)
}






GetSeqs <- function ( n , num.sims , path ) {
	#recover()
	a<-system(paste("grep segsites ", path , "/myseqdata",sep=""),intern=TRUE)
	seg.sites<-sapply(a,function(b){as.numeric(strsplit(b,":")[[1]][2])})
	polymorph<- seg.sites>0
	seq.lines<-c(0,cumsum(polymorph*n)[-length(polymorph)])	
	my.seqs <- lapply ( 0 : ( num.sims - 1 ) , function ( iter ) {		
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
CountHaps <- function ( these.seqs ) {
	recover ()
	positions <- as.numeric ( these.seqs [[ 1 ]] )
	my.part <- list ( seq ( 1 : nrow ( these.seqs [[ 2 ]] ) ) )
	hap.freqs <- matrix ( 0 , nrow = nrow ( these.seqs [[ 2 ]] ) , ncol = len.bp / hap.count.interval + 1 )
	pos.cuts <- seq ( 0 , 1 , by = hap.count.interval / len.bp )
	hap.freqs [ 1 , pos.cuts < positions [ 1 ] ] <- nrow ( these.seqs [[ 2 ]] )
	for ( i in 1 : length ( positions ) ) {
		#if ( i == 18 ) break
		this.site <- which ( these.seqs [[ 2 ]] [ , i ] == 1 )
		if ( any ( unlist ( lapply ( my.part ,  function ( x ) all ( x %in% this.site ) ) ) ) ) { 
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
		these.hap.freqs <- numeric ( length ( unlist ( my.part ) ) )
		
		hap.freqs [ 1 : length ( my.part ) , pos.cuts >= positions [ i ] & pos.cuts < positions [ i + 1 ] ] <- sort ( unlist ( lapply ( my.part , length ) ) , d = T )
		if ( length ( my.part ) == nrow ( hap.freqs ) ) break
	}
	
	
	
	
}


SplitPartition <- function ( )





get.freq.spec<-function(n,num.sims, path){
	recover()
	a<-system(paste("grep segsites ", path , "/myseqdata",sep=""),intern=TRUE)
	seg.sites<-sapply(a,function(b){as.numeric(strsplit(b,":")[[1]][2])})
	polymorph<- seg.sites>0
	seq.lines<-c(0,cumsum(polymorph*n)[-length(polymorph)])	
	freq.specs<-sapply(0:(num.sims-1),function(iter){		
		#recover()
		
		#####
		# if(!polymorph[1+iter]) {freq.spec<-rep(0,n);return(freq.spec)}
		
		positions<-read.table(paste(path, "/myseqdata",sep=""),skip=5+4*iter+seq.lines[iter+1],nrow=1)
#		print(positions[1])
#		if(length(positions)==1){freq.spec<-rep(0,n);return(freq.spec)}		
		seqs.raw<-scan(paste(path, "/myseqdata",sep=""),skip=6+4*iter+seq.lines[iter+1],nline=n,what=character(),quiet=TRUE)
		seqs<-sapply(seqs.raw,function(seq){as.numeric(strsplit(seq,"")[[1]])})
		colnames(seqs)<-NULL
		seqs<-t(seqs)
		these.pos<-positions[-1]>0.5   ###why the -1 here? oh because positions has label
		if(sum(these.pos)==0){freq.spec<-rep(0,n);return(freq.spec)}
		seqs<-seqs[,these.pos] ##throw out first 1/2 of seq.
		if(sum(these.pos)==1){freq.spec<-(1:n==sum(seqs)); return(freq.spec)}
		mut.freq<-colSums(seqs)
		freq.spec<- sapply(1:n,function(i){sum(mut.freq==i)})
		return(freq.spec)
	})
	#freq.specs <- rowSums ( freq.specs )
	return(freq.specs)
}





my.runs <- SweepFromStandingSim ( N = 10000 , s = 0.05 , f = 1/20000 , reps = 100 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE  , display.rep.count = TRUE , time.factor = 1  )



msHapSims ( my.runs [[ 1 ]] , n.sam = 20 , f = 1/20000 , s = 0.05 , N = 10000 , path = "Sims/HapSims" , len.bp = 2000000 , r.bp = 10^-8 , mu.bp = 10^-8 , hap.count.interval = 1000 )







