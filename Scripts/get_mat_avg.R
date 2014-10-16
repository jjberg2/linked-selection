# R script to load LD matrix objects and take the average
args <- commandArgs(trailingOnly=T)
my.names <- list.files ( pattern = args[1] )
my.mats <- array ( data = NA , dim = c ( args [ 2 ] , args [ 2 ] , length ( my.names ) ) )
for ( i in seq_along ( my.names ) ) {
   my.mats [ , , i ] <- get ( load ( my.names [ i ] ) ) 
}
mat.avg <- rowMeans ( my.mats , dims = 2 , na.rm = T )
save ( mat.avg , file =  paste ( args [ 1 ] , "_mean.Robj", sep = "" ) )
