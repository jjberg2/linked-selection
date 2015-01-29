

##extract multiple runs results
worked<-numeric()
	hap.counts<-lapply(1:10,function(i){0})
	coal.times<-lapply(1:10,function(i){numeric()})
for(run in 1:15){
	 load(paste("run_cond_lost_",run,".Robj",sep=""))
	 worked<-rbind(worked,sapply(lapply(my.runs,names),length)>0)
	 for(f.index in 1:10){
	 	if(length(names(my.runs[[f.index]]))){ ##did run run
	 		##if hap.counts is currently blank, then use 0 ifelse(length( hap.counts[[f.index]]),
	 		hap.counts[[f.index]]<- hap.counts[[f.index]]+ 
	 							my.runs[[f.index]]$hap.dist$hap.count.freqs.by.interval
	 		
	 		coal.times[[f.index]]<-rbind(coal.times[[f.index]],my.runs[[f.index]]$coal.times)
	 		}
	 }
	 print(colSums(wosrked))
	} 
worked<-colSums(worked)	 
	 
	 save(file="Coal_Sims_w_traj.Robj",hap.counts,worked,coal.times)
	 
	 pdf(file="Ewens_vs_Jeremy_many_runs_cond_on_loss.pdf"); for(i in 1:10){ MakeHapPlots ( hap.counts[[i]]/worked[i] , N = 10000, f = fs[i], sim.distance = 0.02,plot.cumulative=FALSE);title=paste("f=",fs[i])}
	 
	 
	 
	my.mean<- numeric()
	coeff.var<-numeric()
for(i in 1:10){
coal.diffs<-apply(cbind(0,coal.times[[i]]),1,diff)
my.std.dev<-apply(coal.diffs,1,sd)/(2*10000)
my.mean<-rbind(my.mean,apply(coal.diffs,1,mean)/(2*10000))
coeff.var<-rbind(coeff.var,my.std.dev/my.mean[i,])

}

	my.cols<- rainbow(10)
	
fs<-(1:10)/100
	pdf(file="mean_coal_times_derived.pdf")
	plot(range(fs),c(0,.1),type="n",xlab="f",ylab="Expected time while k lineages")
	sapply(1:9,function(i){points(fs,my.mean[,i],bg=my.cols[i],type="b",pch=21)})
expect.times<-sapply(fs,function(f){j<-10:2;(f*2/(j*(j-1)))})
apply(expect.times,1,lines,x=fs)
legend("topleft",legend=paste("k=",10:2),pch=21,pt.bg=my.cols) 
dev.off()

	
	  pdf(file="coeff_var_coal_times.pdf")
	plot(range(fs),c(0.8,3),type="n",xlab="f",ylab="Coeff. of var. of time while k lineages")
	sapply(1:9,function(i){points(fs,coeff.var[,i],bg=my.cols[i],type="b",pch=21)})
abline(h=1)
legend("topright",legend=paste("k=",10:2),pch=21,pt.bg=my.cols) 
dev.off()


pdf(file="n_two_coal_time.pdf")
g<-function(f,x){x/(f*(x+(1-x)*f))*(2-f+2*(1-f)*log(1-f)/f) }

x<-seq(0,.1,length=100); 
plot(x,sapply(x,function(x){2*integrate(g,0,1,x=x)$value}),ylim=c(0,.1),lwd=2,type="l",xlab="f",ylab="Expected pairwise coal. time")
k<-10:2
prob.coal.when.i<-c(1,cumprod((1-1/(choose(k[-length(k)],2)))))*1/(choose(k,2))
for(i in 1:10){
	coal.cum.sum<-cumsum(my.mean[i,])
	points(fs[i],sum(coal.cum.sum*prob.coal.when.i),bg="red",pch=21)
	}
abline(0,1,lty=2)
dev.off()

for(i in 1:10){
		 pdf(file=paste("Ewens_vs_Jeremy_many_runs_cond_on_loss_f_",i,".pdf",sep="")); 
	 MakeHapPlots ( hap.counts[[i]]/worked[i] , N = 10000, f = fs[i], sim.distance = 0.02,plot.cumulative=FALSE);title=paste("f=",fs[i]); dev.off()}
	 }



###########old stuff for single run

fs<-(1:10)/100

plot(range(fs),c(0,0.15),type="n",xlab="f",ylab="expected time till there's k lineages")
for(i in 1:10){points(rep(fs[i],9),freq.coal[[i]]$mean.coalescence.times/(2*10000))}

expect.times<-sapply(fs,function(f){j<-10:2;cumsum(f*2/(j*(j-1)))})
apply(expect.times,1,lines,x=fs)

plot(range(fs),c(0,0.2),type="n",xlab="f")


plot(c(0,9),c(0,3),type="n")
for(i in 1:10){

coal.diffs<-apply(cbind(0,freq.coal[[i]]$coal.times),1,diff)
my.std.dev<-apply(coal.diffs,1,sd)/(2*10000)
my.mean<-apply(coal.diffs,1,mean)/(2*10000)

coeff.var<-my.std.dev/my.mean
print(coeff.var)
lines(coeff.var)
}


g<-function(f,x){x/(f*(x+(1-x)*f))*(2-f+2*(1-f)*log(1-f)/f) }

x<-seq(0,.1,length=100); plot(x,sapply(x,function(x){2*integrate(g,0,1,x=x)$value}))

prob.coal.when.i<-c(1,cumprod((1-1/(choose(k[-length(k)],2)))))*1/(choose(k,2))
for(i in 1:10){points(fs[i],sum(freq.coal[[i]]$mean.coalescence.times/(2*10000)*prob.coal.when.i),col="red")}
abline(0,1)


pdf(file="Ewens_vs_Jeremy.pdf"); for(i in 1:10){ MakeHapPlots ( freq.coal[[i]]$hap.dist$hap.count.freq.by.interval , N = 10000, f = fs[i], sim.distance = 0.02,plot.cumulative=FALSE);title=paste("f=",fs[i])}

