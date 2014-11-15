source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R', chdir = TRUE)
####frequency spectrum

real.fs <- c ( 0.001 , 0.005 , 0.01 , seq ( 0.02 , 0.2 , length = 10 ) )

#my.runs <- lapply ( real.fs , function ( x ) SweepFromStandingSim ( N = 10000 , s = 0.01 , f = x , reps = 10 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T ) )


#save ( my.runs , file = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Sims/11000freq.trajectories.Rdata")


load( file = "~/../Shared/11000freq.trajectories.Rdata")  #graham's work machine
#load ( file = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Sims/11000freq.trajectories.Rdata" )    ##Jeremy's machine

this.comp <-Sys.info()["nodename"]
setwd ( "~/Documents/Academics/StandingSweeps/" )  ##Jeremy's machine
##path = "~/Dropbox/Linked_selection_models/Soft_sweeps_coal/LinkedSelection/" #graham's work machine




expected.freq.times.standing(n=10,N=10000,r = 0.0025 , f = 0.05 )

####freq. spectrum w. no rec. during sweeps.

expected.freq.times.standing<-function(n,N,r,distance,f){
	recover()
	#my.StirlingNumbers<-StirlingNumbers(n) 
	ESF.prob.k<-EwensDist( n=n , N =N, r=r , distance=1 , f=f) # ,stirling.numbers=my.StirlingNumbers)    ### is of form [n,k]

	my.StirlingNumbers<-StirlingNumbers(n)    ##Usigned Stirling numbers of 1st kind. ma


	expected.t.l<-rep(NA,n-1)
	p_l_given_k <- array ( NA , dim = c ( n , n  ,n ) )
	for(l in 1:(n-1)){
		freq.specs<-sapply(1:n,function(k){
			freq.spec<-rep(NA,n)
			freq.spec[1:k]<-(1/(1:k))
			return(freq.spec)
		}
		)	
		freq.specs<-t(freq.specs)
	#	recover()
		terms.in.sum<-rep(0,n)
		for(k in 2:n) {
			##runs from 2 otherwise there are no polymorphism
			for(j in 1:(k-1)){
				if ( l > k ) next
				
				stirling.bit<-	my.StirlingNumbers[l,j] * my.StirlingNumbers[n-l,k-j]  / my.StirlingNumbers[n,k]
				p_l_given_k [ l , k , j ] <- stirling.bit*choose(n,l)/choose(k,l)
				if(!is.finite(p_l_given_k [ l , k , j ])){ stop ("is infinite") }  ##cat("problem",l,k,j," "); ##is this right?
				terms.in.sum[k]<-terms.in.sum[k]+ESF.prob.k[n,k] *p_l_given_k [ l , k , j ]*freq.specs[k,j]
				stopifnot(is.finite(terms.in.sum[k])) 
			}		
		}
		expected.t.l[l]<-sum(terms.in.sum)
	}
	
	return(expected.t.l)
}




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
	

