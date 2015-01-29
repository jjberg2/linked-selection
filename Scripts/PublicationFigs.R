#### plot of pi and segregating sites
setwd ( "~/Documents/Academics/StandingSweeps" )
source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R')
source ( "~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R")


##### num haps
# EwensHaps40Sim <- StructuredCoalescentSweep ( N = 10000 , s = 0.01 , dominance = FALSE , f = 0.025 , reps = 1000 , n.tips = 40 , r = 10^-8 , sim.distance = 0.02 , interval.width = 1000 , no.sweep = TRUE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 2 )
# save ( EwensHaps40Sim , file = "Sims/EwensHaps40Sim.Robj" )
# load ( file = "Sims/EwensHaps40Sim.Robj" )


EwensHaps10Sim <- StructuredCoalescentSweep ( N = 10000 , s = 0.01 , dominance = FALSE , f = 0.025 , reps = 1000 , n.tips = 10 , r = 10^-8 , sim.distance = 0.02 , interval.width = 1000 , no.sweep = TRUE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 2 )
save ( EwensHaps10Sim , file = "Sims/EwensHaps10Sim.Robj" )

# # 
# ##### pi and seg sites
# real.fs <- c ( 1 / 20000 , 0.001 , 0.01 , 0.05 , 0.1 )

# T_f <- numeric ( length ( real.fs ) )
# for ( i in 1 : length ( real.fs ) ) {
	
	# T_f [ i ] <- log ( (2*N -1 ) * ( 1 - real.fs [ i ] ) / real.fs [ i ] ) / s
	
# }


##made by Coal_Sims_w_traj.Robj
show(load("~/Dropbox/Linked_selection_models/Soft_sweeps_coal/LinkedSelection/Scripts/Coal_Sims_w_traj.Robj"))

fs<-(1:10)/100



layout(t(1:2))
par(mar=c(3,3,1,0))
i=1; MakeHapPlots ( hap.counts[[i]]/worked[i] , N = 10000, f = fs[i], sim.distance = 0.02,plot.cumulative=FALSE,do.legend=TRUE)
mtext(side=2,line=2,text="Probability")
mtext(side=1,line=2,text="4Nr")
mtext(side=3,line=0,text="f=1%")

 i=5; MakeHapPlots ( hap.counts[[i]]/worked[i] , N = 10000, f = fs[i], sim.distance = 0.02,plot.cumulative=FALSE)
#mtext(side=2,line=2,text="Probability")
mtext(side=1,line=2,text="4Nr")
mtext(side=3,line=0,text="f=5%")
dev.copy2pdf(file="~/Dropbox/Linked_selection_models/Soft_sweeps_coal/LinkedSelection/Paper/Paper_Figures/Prob_hap_distribution.pdf")

	my.cols<- rainbow(10)

	
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


