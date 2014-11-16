

#args <- commandArgs(trailingOnly=T)
run.ms.f <- function ( runs , n.sam = 2  , f , s , N , path , ext = "", get.site.density = TRUE , recom = FALSE ) {
	#recover()
	options ( "scipen" = 100 , "digits" = 4 )
	f.lab <- strsplit ( as.character ( f ) , "\\." ) [[ 1 ]] [ 2 ]
	s.lab <- strsplit ( as.character ( s ) , "\\." ) [[ 1 ]][ 2 ]
	counter <- 1
	
	my.file <- paste ( path , "Sims/mssel_f" , n.sam ,  f.lab  , s.lab , N  , ".out" , sep = "" )
	num.sims<-20
	system ( paste ( "rm " , my.file ) )
	#for ( run in 1:5 ) {
	#	load ( paste ( "run_cond_lost_" , run , ".Robj" , sep = "" ) )

	if (! get.site.density ) {
		my.specs <- matrix ( NA , nrow = n.sam , ncol = num.sims * nrow ( runs ) )
	}
	
	for ( i in 1: nrow ( runs ) ) {

		my.freqs <- runs [ i , runs [ i , ] > 0 ]
		my.times <- 0 : length ( my.freqs )
		my.freqs <- c ( my.freqs , 0 )
		
		my.times <- my.times / ( 4*N  )

		#recover()
		header.material <- c ( "1" , "1" , paste ( "n:" , length ( my.times ) ) )
		write ( file = paste ( path , "Sims/my.standing" , "." , f.lab , "." , s.lab , "." , N, "." ,ext , ".traj" , sep = "" ) , header.material )
		write.table ( file = paste ( path , "Sims/my.standing" , "." , f.lab , "." , s.lab , "." , N , "." ,ext, ".traj" , sep = "" ) , cbind ( my.times , my.freqs ) , append = TRUE , sep = "\t" , quot = FALSE , col.nam = FALSE , row.name = FALSE )
		cat( i ," " )
		if ( get.site.density ) { 
			system ( paste ( path , "Scripts/msseldir/mssel " , n.sam , " 20 0 " , n.sam , " " , path , "Sims/my.standing" , "." , f.lab , "." , s.lab , "." , N, ".traj 0 -t 200. -r 200. 20000 | grep pos | cut -f 2 -d : >> " , my.file , sep = "" ) )
		}	else	{   ##setup for the mo. to do freq. spectrum
			system ( paste ( "Scripts/msseldir/mssel " , n.sam , " " , 20 , " 0 " , n.sam , " Sims/my.standing" , "." , f.lab , "." , s.lab , "." , N, "." ,ext, ".traj 0 -t 200. -r " , recom , " 2 >",path, "Sims/myseqdata" , sep = "" ) ) 
			
			spec <- get.freq.spec ( n.sam , num.sims = num.sims, path=path )
			my.specs[,(1+(counter-1)*num.sims):(counter*num.sims)]<-spec
			counter<-counter+1
			
			#recover()
		}
	}
	if (! get.site.density ) return(my.specs)
}


get.mut.density<-function(file){
	myres = scan(file)
	##myres*theta, note range of smoothing
	myres<-myres*200
	
	mydens = density(myres,from=3,to=190,bw=2.0,na.rm=T)
	mydens$y<-mydens$y*length(myres)
	return(mydens)
}

get.freq.spec<-function(n,num.sims, path){
	#recover()
	a<-system(paste("grep segsites ","Sims/myseqdata",sep=""),intern=TRUE)
	seg.sites<-sapply(a,function(b){as.numeric(strsplit(b,":")[[1]][2])})
	polymorph<- seg.sites>0
	seq.lines<-c(0,cumsum(polymorph*n)[-length(polymorph)])	
	freq.specs<-sapply(0:(num.sims-1),function(iter){		
		if(!polymorph[1+iter]) {freq.spec<-rep(0,n);return(freq.spec)}
		positions<-read.table(paste(path, "Sims/myseqdata",sep=""),skip=5+4*iter+seq.lines[iter+1],nrow=1)
#		print(positions[1])
#		if(length(positions)==1){freq.spec<-rep(0,n);return(freq.spec)}		
		seqs.raw<-scan(paste(path, "Sims/myseqdata",sep=""),skip=6+4*iter+seq.lines[iter+1],nline=n,what=character(),quiet=TRUE)
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
return(freq.specs)
}


