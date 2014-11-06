####frequency spectrum

real.fs <- c ( 0.001 , 0.005 , 0.01 , seq ( 0.02 , 0.2 , length = 10 ) )

#my.runs <- lapply ( real.fs , function ( x ) SweepFromStandingSim ( N = 10000 , s = 0.01 , f = x , reps = 10 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T ) )


#save ( my.runs , file = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Sims/11000freq.trajectories.Rdata")


load( file = "~/../Shared/11000freq.trajectories.Rdata")  #graham's work machine
#load ( file = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/Sims/11000freq.trajectories.Rdata" )    ##Jeremy's machine

#Sys.info()["nodename"]
#path = "~/Documents/Academics/CoopLab/Projects/StandingSweeps/"  ##Jeremy's machine
path = "~/Dropbox/Linked_selection_models/Soft_sweeps_coal/LinkedSelection/" #graham's work machine




expected.freq.times.standing(n=10,N=N,r=,distance,f)

####freq. spectrum w. no rec. during sweeps.

expected.freq.times.standing<-function(n,N,r,distance,f){

#my.StirlingNumbers<-StirlingNumbers(n) 
ESF.prob.k<-EwensDist( n=n , N =N, r=r , distance=distance , f=f) # ,stirling.numbers=my.StirlingNumbers)    ### is of form [n,k]

my.StirlingNumbers<-StirlingNumbers(n)    ##Usigned Stirling numbers of 1st kind. ma


expected.t.l<-rep(NA,n)
for(l in 1:(n-1)){

	freq.specs<-sapply(1:n,function(k){
		freq.spec<-rep(NA,n)
		freq.spec[1:k]<-(1/(1:k))
		return(freq.spec)
	})
	
	freq.specs<-t(freq.specs)
#	recover()
	terms.in.sum<-rep(0,n)
	for(k in 2:n){   ##runs from 2 otherwise there are no polymorphism
			for(j in 1:(k-1)){
			stirling.bit<-	my.StirlingNumbers[l,j] * my.StirlingNumbers[n-l,k-j]  / my.StirlingNumbers[n,k]
			p_l_given_k<-stirling.bit*choose(n,l)/choose(k,l)
			if(!is.finite(p_l_given_k)){ p_l_given_k<-0 }  ##cat("problem",l,k,j," "); ##is this right?
			terms.in.sum[k]<-terms.in.sum[k]+ESF.prob.k[n,k] *p_l_given_k*freq.specs[k,j]
			stopifnot(is.finite(terms.in.sum[k])) 
		}		
	}

	expected.t.l[l]<-sum(terms.in.sum)
}

return(expected.t.l)
}


#####get frequency spectrum from ms.
#modified from 
### function get.freq.spec() in freq_spec.R
##get.freq.spec<-function(n,num.sims)

# run.ms.f<-function(f.index,n.sam=2, get.site.density=TRUE,recom=FALSE){
	# num.sims<-20
	# my.file<-paste("mssel_f",n.sam,f.index,".out",sep="")
	# counter<-1
	
	# my.specs<-matrix(NA,nrow=n.sam,ncol=10000)
	# system(paste("rm ",my.file))
	# for(run in 1:5){
		# load(paste("run_cond_lost_",run,".Robj",sep=""))
		# for(i in 1:100){
			# my.freqs<-my.runs[[f.index]]$trees[[i]]$freqs
			# my.times<-0:length(my.freqs)
			# my.freqs<-c(my.freqs,0)
			
			# my.times<-my.times / (4*10e3)
	
	# #		recover()
			# header.material<-c("1","1",paste("n:",length(my.times) ))
			# write(file=paste("my.standing",f.index,".traj",sep=""),header.material)
			# write.table(file=paste("my.standing",f.index,".traj",sep=""),cbind(my.times,my.freqs), append=TRUE, sep="\t",quot=FALSE,col.nam=FALSE,row.name=FALSE)
			# cat(i," ")
			# if(get.site.density){ 
				# system(paste("msseldir/mssel ",n.sam," ",num.sims," 0 " ,n.sam," my.standing",f.index,".traj 0 -t 200. -r 200. 20000 | grep pos | cut -f 2 -d : >> ",my.file,sep=""))
				# }else{   ##setup for the mo. to do freq. spectrum
				# system(paste("msseldir/mssel ",n.sam," ",num.sims," 0 " ,n.sam," my.standing",f.index,".traj 0 -t 200. -r ",recom, " 2 > myseqdata",sep="")) 
				# spec<-get.freq.spec(n.sam,num.sims)
# #				recover()
				
				# my.specs[,(1+(counter-1)*num.sims):(counter*num.sims)]<-spec
				# counter<-counter+1
			# }
		# }
	# }
	# return(my.specs)
# }



# get.freq.spec<-function(n,num.sims){
	# a<-system("grep segsites myseqdata",intern=TRUE)
	# seg.sites<-sapply(a,function(b){as.numeric(strsplit(b,":")[[1]][2])})
	# polymorph<- seg.sites>0
	# seq.lines<-c(0,cumsum(polymorph*n)[-length(polymorph)])	
	# freq.specs<-sapply(0:(num.sims-1),function(iter){		
		# if(!polymorph[1+iter]) {freq.spec<-rep(0,n);return(freq.spec)}
		# positions<-read.table("myseqdata",skip=5+4*iter+seq.lines[iter+1],nrow=1)
# #		print(positions[1])
# #		if(length(positions)==1){freq.spec<-rep(0,n);return(freq.spec)}		
		# seqs.raw<-scan("myseqdata",skip=6+4*iter+seq.lines[iter+1],nline=n,what=character(),quiet=TRUE)
		# seqs<-sapply(seqs.raw,function(seq){as.numeric(strsplit(seq,"")[[1]])})
		# colnames(seqs)<-NULL
		# seqs<-t(seqs)
		# these.pos<-positions[-1]>0.5   ###why the -1 here? oh because positions has label
		# if(sum(these.pos)==0){freq.spec<-rep(0,n);return(freq.spec)}
		# seqs<-seqs[,these.pos] ##throw out first 1/2 of seq.
		# if(sum(these.pos)==1){freq.spec<-(1:n==sum(seqs)); return(freq.spec)}
		# mut.freq<-colSums(seqs)
		# freq.spec<- sapply(1:n,function(i){sum(mut.freq==i)})
		# return(freq.spec)
	# })
# return(freq.specs)
# }


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
	

my.freq.specs<-list(20)
recs<-seq(0,200,length=20)

for(my.rec.index in 1:length(recs)){
	my.freq.specs[[my.rec.index]] <-run.ms.f(f.index=f.index,n.sam=10, get.site.density=FALSE,recom=recs[my.rec.index])
	save(file=paste("ms_standing_freq_spectrum_",f.index,".Robj",sep=""),my.freq.specs)
}

## reads ms output from myseqdata
###  run.ms.f<-function(f.index,n.sam=2, get.site.density=TRUE)


load(file=paste("~/Dropbox/Linked_selection_models/Soft_sweeps_coal/LinkedSelection/Scripts/ms_standing_freq_spectrum_",f.index,".Robj",sep=""))
