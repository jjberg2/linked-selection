<<<<<<< HEAD
=======
source('~/Documents/Academics/CoopLab/Projects/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
#source('~/Documents/Academics/StandingSweep/Scripts/SweepFromStandingSim.R', chdir = TRUE)
source('~/Documents/Academics/StandingSweep/Scripts/run.ms.f.R', chdir = TRUE)
>>>>>>> FETCH_HEAD


run.ms.f <- function ( runs , n.sam = 2  ,f.index, N , path , get.site.density = TRUE , recom = FALSE ) {
	#recover()
	my.file <- paste ( path , "Sims/mssel_f" , n.sam , f.index , ".out" , sep = "" )
	num.sims<-20
	system ( paste ( "rm " , my.file ) )
	#for ( run in 1:5 ) {
	#	load ( paste ( "run_cond_lost_" , run , ".Robj" , sep = "" ) )
	if (! get.site.density ){ my.specs<-matrix(NA,nrow=n.sam,ncol=num.sims*nrow ( runs )); counter<-1}
	
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
			system ( paste (path, "Scripts/msseldir/mssel " , n.sam , " " , 20 , " 0 " , n.sam  , " " , path , "Sims/my.standing"  , f.index , ".traj 0 -t 200. -r " , recom , " 2 > ",path, "Sims/myseqdata" , sep = "" ) ) 
			
			spec <- get.freq.spec ( n.sam , num.sims = num.sims, path=path )
			my.specs[,(1+(counter-1)*num.sims):(counter*num.sims)]<-spec
			counter<-counter+1
			
			#recover()
		}
	}
	if (! get.site.density ) return(my.specs)
}



# example of running one set of frequencies
#run.ms.f ( runs = my.runs [[ 1 ]] [[ 1 ]] , n.sam = 2 , N = 10000 , path = "~/Documents/Academics/StandingSweep" )

#load ( file = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Sims/11000freq.trajectories.Rdata" )

#for ( f.index in 1 : length ( real.fs ) ) {
#	run.ms.f ( runs = my.runs [[ f.index ]] [[ 1 ]] , n.sam = 12 , N = 10000 , path = "~/Documents/Academics/StandingSweep/" )
#}
get.freq.spec<-function(n,num.sims, path){
	a<-system(paste("grep segsites ", path,"Sims/myseqdata",sep=""),intern=TRUE)
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



if(FALSE){
run.ms.f ( runs = my.runs [[ 1 ]] [[ 1 ]] ,f.index=1, n.sam = 2 , N = 10000 , path = path )


for ( f.index in 1 : length ( real.fs ) ) {
	run.ms.f ( runs = my.runs [[ f.index ]] [[ 1 ]] ,f.index=f.index ,n.sam = 20 , N = 10000 , path = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/" )
}


my.freqs.specs<- run.ms.f ( runs = my.runs [[ 1 ]] [[ 1 ]] ,f.index=1, n.sam = 10 , N = 10000 , path = path,get.site.density = FALSE , recom = 100)

pdf(file="pi_density.pdf")
	plot(c(0,200),c(0,1),type="n",xlab="4NR",ylab=expression(pi[R]/pi[0]),cex.lab=1.5)
	pi.over.f<-list()
	for ( f.index in 1 : length ( real.fs ) ) {
	#run.ms.f(f.index)
		mut.density <- get.mut.density ( file = paste ( path,  "Sims/mssel_f2" , f.index , ".out" , sep = "" ) )
		pi.over.f [[ f.index ]] <- mut.density
		lines ( mut.density$x , mut.density$y/(1000*20) , col = f.index )
	}
	fs<-1:5/100
	rho=seq(0,200,length=1000)
	for(i in 1:5){
	
	f=fs[i]
	s=0.01
	
	lines(rho,1-1/(1+rho*f*(1-f)),col=i,lwd=2,lty=2)
	# lines(r,1-1/(1+4*N*r*f*(1-f))*exp(-2*r/s*log(1/f)),col=i,lwd=2)
	# lines(r,1-exp(-2*r*(1/s*log(1/f)+2*N*f*(1-f))),col=i,lty=2)
	
	}
	legend("bottomright",legend= paste("f=",fs),col=1:5,lty=1,lwd=2,cex=1.5)
dev.off()
 
 

N=1e4
 r = 10^-8 ; interval.width = 1000; sim.distance = 0.05 
	sim.distance.bp <- sim.distance / r 
	intervals <- seq ( 0 , sim.distance.bp , interval.width )
 
 pdf(file="segsite_density.pdf")
plot(c(0,200),c(0,2.5),type="n",xlab="4NR",ylab="Num. seg. sites",cex.lab=1.5)
pi.over.f<-list()
 fs<-1:5/100
for(f.index in 1:5){
#run.ms.f(f.index)
	mut.density<-get.mut.density(file=paste("mssel_f10",f.index,".out",sep=""))
	pi.over.f[[f.index]]<-mut.density
	lines(mut.density$x,mut.density$y/(5*100*20),col=f.index,lty=1,lwd=2)

	f<-fs[f.index]
	ESF.prob.k<-sapply(intervals,function(distance){EwensDist( n=n , N =N, r=r , distance=distance , f=f )[n,]})
 	
	num.seg.sites.cond.i<-c(0,sapply(1:9,function(i){ sum(1/(1:i)) }))
	num.seg.sites.r<-colSums(apply(ESF.prob.k,2,function(x){x*num.seg.sites.cond.i}))

	lines(4*N*intervals*r,num.seg.sites.r,col=f.index,lty=2,lwd=2) 	
 	
 }
 legend("bottomright",legend= paste("f=",fs),col=1:5,lty=1,lwd=2,cex=1.5)
dev.off()

}
 