source('~/Documents/Academics/CoopLab/Projects/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
#source('~/Documents/Academics/StandingSweep/Scripts/SweepFromStandingSim.R', chdir = TRUE)
source('~/Documents/Academics/StandingSweep/Scripts/run.ms.f.R', chdir = TRUE)


real.fs <- c ( 0.001 , 0.005 , 0.01 , seq ( 0.02 , 0.2 , length = 10 ) )

my.runs <- lapply ( real.fs , function ( x ) SweepFromStandingSim ( N = 10000 , s = 0.01 , f = x , reps = 10 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T ) )


#save ( my.runs , file = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Sims/11000freq.trajectories.Rdata")

#my.runs <- lapply ( real.fs , function ( x ) SweepFromStandingSim ( N = 10000 , s = 0.01 , f = x , reps = 1000 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T ) )
#save ( my.runs , file = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Sims/11000freq.trajectories.Rdata")

#load( file = "~/../Shared/11000freq.trajectories.Rdata")  #graham's work machine
load ( file = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Sims/11000freq.trajectories.Rdata" )    ##Jeremy's machine

#Sys.info()["nodename"]
path = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/"  ##Jeremy's machine
path = "~/Dropbox/Linked_selection_models/Soft_sweeps_coal/LinkedSelection/" #graham's work machine




# example of running one set of frequencies
run.ms.f ( runs = my.runs [[ 1 ]] [[ 1 ]] , n.sam = 2 , N = 10000 , path = "~/Documents/Academics/StandingSweep" )

#load ( file = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Sims/11000freq.trajectories.Rdata" )

for ( f.index in 1 : length ( real.fs ) ) {
	run.ms.f ( runs = my.runs [[ f.index ]] [[ 1 ]] , n.sam = 12 , N = 10000 , path = "~/Documents/Academics/StandingSweep/" )
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
 