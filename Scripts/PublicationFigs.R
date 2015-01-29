setwd ( "~/Documents/Academics/StandingSweeps" )
source('~/Documents/Academics/StandingSweeps/Scripts/SweepFromStandingSim.R')
source ( "~/Documents/Academics/StandingSweeps/Scripts/run.ms.functions.R")


#################
## extra functions ##
#################


# MyLogistic <- function ( x , N = 10000 , s  ) 1 / ( 2 * N  ) * exp(s * x ) / ( 1 + 1 / (2 * N  ) * ( exp(s * x )  - 1 )  )

MyLogistic <- function ( x , N = 10000 , s  ) 1 / ( 5 * N * s ) * exp(s * x ) / ( 1 + 1 / ( 5 * N * s  ) * ( exp ( s * x )  - 1 )  )


#################
#################



##### trajectories vs approx
f <- 0.01
s <- 0.01
# my.freq.trajs <- SweepFromStandingSim ( N = 10000 , s = s , f = f , reps = 10 , no.sweep = FALSE , cond.on.loss = TRUE , cond.on.fix = TRUE , time.factor = 1 , display.rep.count = T )
# save ( my.freq.trajs , file = "Sims/10trajectoriesForTrajFigure.Robj")
load ( file = "Sims/10trajectoriesForTrajFigure.Robj")
freqs <- t ( my.freq.trajs[[1]] )
freqs <- apply ( freqs , 2 , rev )

i <- 1
det.freqs <- 1 / ( 20000 * 0.01 )
while ( det.freqs [ i ] < ( 1 - f ) ) {	
	det.freqs [ i + 1 ] <- MyLogistic ( i , N = 10000 , s = s)
	i <- i + 1 
}
det.freqs <- c ( rep ( f , nrow ( freqs) - length ( det.freqs ) ) , rev ( 1 - det.freqs ) )
matplot ( freqs , type = "l" , lwd = 0.7 , col = "grey" , lty = 1 , bty = "n" , ylab = "Frequency" , xlab = "Generations" , xaxt = "n" )
lines ( det.freqs , lwd = 2 )

my.at <- length ( det.freqs ) - seq ( 0 , 3500 , 500)
axis ( 1 , at = my.at , labels = seq ( 0 , 3500 , 500) )


new.freqs <- matrix ( 0 , nrow = max ( my.freq.trajs [[ 2 ]] ) - my.freq.trajs [[ 2 ]] [ which ( my.freq.trajs [[ 1 ]] [ , ncol ( my.freq.trajs [[ 1 ]] ) ] != 0 ) ] + ncol ( my.freq.trajs [[ 1 ]] ) , ncol = nrow ( my.freq.trajs [[ 1 ]] ) )

for ( i in 1 : nrow ( my.freq.trajs [[ 1 ]] ) ) {
	temp <- c ( rep ( 1 , max ( my.freq.trajs [[ 2 ]] ) -  my.freq.trajs [[ 2 ]] [ i ] ) , my.freq.trajs [[ 1 ]] [ i , ]  )
	if ( length ( temp ) > nrow ( new.freqs ) ) {
		new.freqs [ , i ] <- rev ( temp [ 1 : nrow ( new.freqs ) ] )
	 } else { 
	 	new.freqs [ , i ] <- rev ( c ( temp , rep ( 0 , nrow ( new.freqs ) - length ( temp ) ) ) )
	 }
}

matplot ( new.freqs , type = "l" , lwd = 0.7 , col = "grey" , lty = 1 , bty = "n" , ylab = "Frequency" , xlab = "Generations" , xaxt = "n" )
lines ( det.freqs , lwd = 2 )


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


