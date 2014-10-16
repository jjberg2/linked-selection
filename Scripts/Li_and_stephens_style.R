source ("SweepFromStandingSim.R")

N=10000
rBP=1e-8
f=.01

real.fs<-seq(0.01,.1,length=9)
all.like.surfs<-list()
for(real.f.index in 1:length(real.fs)){
like.surfs<-numeric()
print(real.f.index)
temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.5 , f = real.fs[real.f.index] , reps = 20 , n.tips = 10 , r = 10^-8 , sim.distance = 0.02 , interval.width = 1000 , no.sweep = TRUE , constant.freq = FALSE , cond.on.loss = TRUE )

fs<-seq(.005,.2,length=40)
for(tree.index in 1:20){
tree<-temp$trees[[tree.index]]
perm.like<-CalcLike(tree,fs,30)
like.surfs<- cbind(like.surfs,rowMeans(exp(perm.like)))

}
all.like.surfs[[real.f.index]]<-like.surfs
}

pdf(file="../../../temp.pdf")
sapply(1:9,function(i){
	plot(range(fs),c(0,1),type="n"); 
	apply(all.like.surfs[[i]],2,lines,x=fs)
	abline(v=mean(fs[apply(all.like.surfs[[i]],2,which.max)]))
	abline(v=real.fs[i],col="red")
	})
dev.off()


pdf(file="../../../MLEs_and_CI.pdf")
 plot(c(1,6),c(0,.15),type="n",xlab="true f",ylab="estimated f",cex.lab=1.5,axes=FALSE);
 axis(2,cex.ax=1.5)
 box()
 axis(1,at=c(1:5)+.5,labels=c(1:5)/100)
sapply(1:5,function(i){
lines(c(i,i+1),rep(real.fs[i],2),co="grey")
sapply(1:10,function(j){
	log.like<-log(all.like.surfs[[i]][,j]);log.like<-log.like-max(log.like)
	points(x=i+.1*j,fs[which.max(log.like)])
lines(x=rep(i+.1*j,2) , y=fs[range(which(log.like-max(log.like)>-2))],col=i)
})
})
dev.off()



CalcLike<-function(tree,fs,perms){

	rec.loc<-tree$rec.events.off.background$rec.right.off.background
	sim.distance.bp <- tree$sim.distance.bp
	num.haps<-nrow(tree$sequence.structure$right.seq)
	perm.like<-matrix(NA,nrow=length(fs),ncol=perms)
	
	for(perm in 1:perms){
		haps<-tree$sequence.structure$right.seq
		perm.haps<-haps[sample(1:num.haps),]
		for(f.index in 1:length(fs)){
			f<-fs[f.index]
			rec.param<-  2*rBP*N*f*(1-f)
	
			hap.probs<-sapply(1:(num.haps-1),function(k){
				new.hap<-perm.haps[k+1,]
				old.haps<-perm.haps[1:k,,drop=FALSE]
				hap.prob<-CalcProbNewHap(k,new.hap,old.haps,rec.loc,rec.param,sim.distance.bp)
				return(hap.prob)
			})
	
			perm.like[f.index,perm]<-sum(hap.probs)
	#cat(f,sum(hap.probs),"\n")
		}
	}
	perm.like<-perm.like - max(perm.like)
	return(perm.like)
}





CalcProbNewHap<-function(k,new.hap,old.haps,rec.loc,rec.param,sim.distance.bp){

	recom.off<-apply(old.haps,1,function(old.hap){max(which(new.hap == old.hap))})
	
	last.recom<-max(recom.off)
	these.rec<-which(recom.off==last.recom)
	
	choose.this<- length(these.rec)
	if(last.recom<=nrow(rec.loc)){ 
		my.like<-(-rec.param*rec.loc$seq[last.recom]/k)+log(rec.param*choose.this) 
	}else{
		my.like<-(-rec.param*sim.distance.bp/k)+log(choose.this)
		}
	
	return(my.like)
}

