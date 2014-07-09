

args <- commandArgs(trailingOnly=T)
run.ms.f<-function( n.sam , f , s , no.sweep , reps , N , n.loc ){
        #recover()
        source ( "Scripts/SweepFromStandingSim.R")
        #setwd("/Users/JeremyBerg/Documents/Academics/StandingSweep")
        freqs <- SweepFromStandingSim ( N , s , f , reps , no.sweep , cond.on.loss = T , cond.on.fix = T , display.rep.count = T ) [[ 1 ]]

        my.trajectories <- apply ( freqs , 1 , function ( x ) cbind ( 0 : ( length ( x [ x != 0 ] ) ) / ( 4*N ) , c ( rev ( x [ x != 0 ]  ) , 0 ) ) )
        #recover()
        ld.stats <- list ()
        for ( i in 1 : length ( my.trajectories ) ) {
            header.material <- c ( "1" , "1" , paste ( "n:" , nrow ( my.trajectories [[ i ]] ) ) )
            write ( header.material , file = paste ( "Working/freq.traj." , i , sep = "" ) )
            write.table ( my.trajectories [[  i ]] , file = paste ( "Working/freq.traj." , i , sep = "" ) , quote = F , row.names = F , col.names = F , append = T , sep = "\t" )
            system(paste("Scripts/msseldir/mssel ",n.sam," 1 0 ",n.sam," Working/freq.traj.",i," 20 -t 10000.0 -r 2000.0 " , n.loc , " > Output/myseqdata" , i ,sep="") )
        }
}

#setwd("/Users/JeremyBerg/Documents/Academics/StandingSweep")
run.ms.f ( n.sam = as.numeric(args[1]) , f = as.numeric(args[2]) , s = as.numeric(args[3]) , no.sweep = args[4] , reps = as.numeric(args[5]) , N = as.numeric(args[6]) , n.loc = as.numeric(args[7]) )
for ( i in 1 : args[5] ) {
    system ( paste ( "Rscript Scripts/LDSimCalcFunc.R Output/myseqdata" , i , args[7] , sep = "" ) )
    # ld.stats [[ i ]] <- LDSimCalc ( paste ( "Output/myseqdata" , i , sep = "" ) )
    # ld.stats [[ i ]] [[ 4 ]] <- cut ( ld.stats [[ i ]] [[ 3 ]] , 0:n.loc/n.loc , include.lowest = T )
    #cat ( i )
}
