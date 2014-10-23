

#args <- commandArgs(trailingOnly=T)
run.ms.f <- function ( runs , n.sam = 2  ,f.index, N , path , get.site.density = TRUE , recom = FALSE ) {
	#recover()
	my.file <- paste ( path , "Sims/mssel_f" , n.sam , f.index , ".out" , sep = "" )
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

#		recover()
		header.material <- c ( "1" , "1" , paste ( "n:" , length ( my.times ) ) )
		write ( file = paste ( path , "Sims/my.standing" , f.index , ".traj" , sep = "" ) , header.material )
		write.table ( file = paste ( path , "Sims/my.standing" , f.index , ".traj" , sep = "" ) , cbind ( my.times , my.freqs ) , append = TRUE , sep = "\t" , quot = FALSE , col.nam = FALSE , row.name = FALSE )
		cat( i ," " )
		if ( get.site.density ) { 
			system ( paste ( path , "Scripts/msseldir/mssel " , n.sam , " 20 0 " , n.sam , " " , path , "Sims/my.standing" , f.index , ".traj 0 -t 200. -r 200. 20000 | grep pos | cut -f 2 -d : >> " , my.file , sep = "" ) )
		}	else	{   ##setup for the mo. to do freq. spectrum
			system ( paste ( "Sims/msseldir/mssel " , n.sam , " " , 20 , " 0 " , n.sam , " my.standing" , f.index , ".traj 0 -t 200. -r " , recom , " 2 >",path, "Sims/myseqdata" , sep = "" ) ) 
			
			spec <- get.freq.spec ( n.sam , num.sims = num.sims, path=path )
			my.specs[,(1+(counter-1)*num.sims):(counter*num.sims)]<-spec
			counter<-counter+1
			
			#recover()
		}
	}
	if (! get.site.density ) return(my.specs)
}
