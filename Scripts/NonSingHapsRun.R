
## s = 0.01
temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.01 , f = 0.01 , reps = 200 , n.tips = 20 , r = 10^-8 , sim.distance = 0.05 , interval.width = 5000 , no.sweep = FALSE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 1 )
MakeHapPlots ( temp$standing.hap.dist$hap.count.freqs.by.interval , N = 10000, f = 0.01, sim.distance = 0.05 , interval.width = 5000)

temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.01 , f = 0.02 , reps = 200 , n.tips = 20 , r = 10^-8 , sim.distance = 0.05 , interval.width = 5000 , no.sweep = FALSE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 1 )
MakeHapPlots ( temp$standing.hap.dist$hap.count.freqs.by.interval , N = 10000, f = 0.02, sim.distance = 0.05 , interval.width = 5000)

temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.01 , f = 0.03 , reps = 200 , n.tips = 20 , r = 10^-8 , sim.distance = 0.05 , interval.width = 5000 , no.sweep = FALSE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 1 )
MakeHapPlots ( temp$standing.hap.dist$hap.count.freqs.by.interval , N = 10000, f = 0.03, sim.distance = 0.05 , interval.width = 5000 )

temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.01 , f = 0.04 , reps = 200 , n.tips = 20 , r = 10^-8 , sim.distance = 0.05 , interval.width = 5000 , no.sweep = FALSE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 1 )
MakeHapPlots ( temp$standing.hap.dist$hap.count.freqs.by.interval , N = 10000, f = 0.04, sim.distance = 0.05 , interval.width = 5000 )

temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.01 , f = 0.05 , reps = 200 , n.tips = 20 , r = 10^-8 , sim.distance = 0.05 , interval.width = 5000 , no.sweep = FALSE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 1 )
MakeHapPlots ( temp$standing.hap.dist$hap.count.freqs.by.interval , N = 10000, f = 0.05, sim.distance = 0.05 , interval.width = 1000 )



## s = 0.03
temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.03 , f = 0.01 , reps = 200 , n.tips = 20 , r = 10^-8 , sim.distance = 0.05 , interval.width = 10000 , no.sweep = FALSE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 1 )
MakeHapPlots ( temp$standing.hap.dist$hap.count.freqs.by.interval , N = 10000, f = 0.01, sim.distance = 0.05 , interval.width = 10000)

temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.03 , f = 0.02 , reps = 200 , n.tips = 20 , r = 10^-8 , sim.distance = 0.05 , interval.width = 10000 , no.sweep = FALSE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 1 )
MakeHapPlots ( temp$standing.hap.dist$hap.count.freqs.by.interval , N = 10000, f = 0.02, sim.distance = 0.05 , interval.width = 10000)

temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.03 , f = 0.03 , reps = 200 , n.tips = 20 , r = 10^-8 , sim.distance = 0.05 , interval.width = 10000 , no.sweep = FALSE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 1 )
MakeHapPlots ( temp$standing.hap.dist$hap.count.freqs.by.interval , N = 10000, f = 0.03, sim.distance = 0.05 , interval.width = 10000)

temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.03 , f = 0.04 , reps = 200 , n.tips = 20 , r = 10^-8 , sim.distance = 0.05 , interval.width = 10000 , no.sweep = FALSE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 1 )
MakeHapPlots ( temp$standing.hap.dist$hap.count.freqs.by.interval , N = 10000, f = 0.04, sim.distance = 0.05 , interval.width = 10000)














## f = 1/2N
temp <- StructuredCoalescentSweep ( N = 10000 , s = 0.01 , f = 1/20000 , reps = 200 , n.tips = 40 , r = 10^-8 , sim.distance = 0.1 , interval.width = 10000 , no.sweep = FALSE , constant.freq = FALSE , cond.on.loss = TRUE , build.seq = TRUE , display.rep.count = TRUE ,  time.factor = 1 )
MakeHapPlots ( temp$standing.hap.dist$hap.count.freqs.by.interval , N = 10000, f = 1/20000, sim.distance = 0.1 , interval.width = 10000)