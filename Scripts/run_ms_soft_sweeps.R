source ("SweepFromStandingSim.R")
run.ms.f<-function(f.index,n.sam=2, get.site.density=TRUE,recom=FALSE){

	my.file<-paste("mssel_f",n.sam,f.index,".out",sep="")

	system(paste("rm ",my.file))
	for(run in 1:5){
		load(paste("run_cond_lost_",run,".Robj",sep=""))
		for(i in 1:100){
			my.freqs<-my.runs[[f.index]]$trees[[i]]$freqs
			my.times<-0:length(my.freqs)
			my.freqs<-c(my.freqs,0)
			
			my.times<-my.times / (4*10e3)
	
	#		recover()
			header.material<-c("1","1",paste("n:",length(my.times) ))
			write(file=paste("my.standing",f.index,".traj",sep=""),header.material)
			write.table(file=paste("my.standing",f.index,".traj",sep=""),cbind(my.times,my.freqs), append=TRUE, sep="\t",quot=FALSE,col.nam=FALSE,row.name=FALSE)
			cat(i," ")
			if(get.site.density){ 
				system(paste("msseldir/mssel ",n.sam," 20 0 ",n.sam," my.standing",f.index,".traj 0 -t 200. -r 200. 20000 | grep pos | cut -f 2 -d : >> ",my.file,sep=""))
				}else{   ##setup for the mo. to do freq. spectrum
				system(paste("msseldir/mssel ",n.sam," ",20," 0 ",n.sam," my.standing",f.index,".traj 0 -t 200. -r ",recom, " 2 > myseqdata",sep="")) 
				spec<-get.freq.spec(n.sam,num.sims=20)
				recover()
			}
		}
	}
}

real.fs <- c ( 0.001 , seq ( 0.01 , 0.1 , length = 10 ) )
my.runs <- lapply ( real.fs , function ( x ) SweepFromStandingSim ( N = 10000 , s = 0.01 , f = x , reps = 100 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 4 ) )


get.mut.density<-function(file){
	myres = scan(file)
	##myres*theta, note range of smoothing
	myres<-myres*200
	
	mydens = density(myres,from=3,to=190,bw=2.0,na.rm=T)
	mydens$y<-mydens$y*length(myres)
	return(mydens)
}



pdf(file="pi_density.pdf")
	plot(c(0,200),c(0,1),type="n",xlab="4NR",ylab=expression(pi[R]/pi[0]),cex.lab=1.5)
	pi.over.f<-list()
	for(f.index in 1:5){
	#run.ms.f(f.index)
	mut.density<-get.mut.density(file=paste("mssel_f",f.index,".out",sep=""))
	pi.over.f[[f.index]]<-mut.density
	lines(mut.density$x,mut.density$y/(5*100*20),col=f.index)
	
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
 



 
 for(f.index in 1:5){
 run.ms.f(f.index=f.index,n.sam=10)
 }
 

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
 