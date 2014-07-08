args <- commandArgs ( trailingOnly = T )
LDSimCalc <- function ( file.name , n.loc ) {
	
	recover()
	my.con <- file ( file.name )
	open ( my.con )
	my.lines <- readLines ( my.con )
	close ( my.con )
	pos.line.id <- grep ( "pos" , my.lines )
	pos.line <- my.lines [ pos.line.id ]
	positions <- as.numeric ( strsplit ( pos.line , " ")[[ 1 ]] [ -1 ] )
	geno.lines <- my.lines [ ( pos.line.id + 1 ) : length ( my.lines ) ]
	split.geno.lines <-  lapply ( strsplit ( geno.lines , "" ) , as.numeric )
	geno.mat <- do.call ( rbind , split.geno.lines )
	#D.stat.sq <- cov ( geno.mat )^2
	#freqs <- colMeans ( geno.mat )
	#rt.vars <- sqrt ( freqs * ( 1 - freqs ) )
	#var.prod.mat <- rt.vars %o% rt.vars
	pos.bins <- cut ( positions , 0:n.loc/n.loc , include.lowest = T )
	r.sq <- matrix ( nrow = n.loc, ncol = n.loc )
	sig.sq.d <- matrix ( nrow = n.loc, ncol = n.loc )
	for ( k in 1:length ( levels ( pos.bins ) ) ) {
		row.bin <- levels ( pos.bins )[ k ]
		for ( j in 1: length ( levels ( pos.bins ) ) ) {
			col.bin <- levels ( pos.bins ) [ j ] 
			D.stat.sq <- cov ( geno.mat [ , pos.bins == col.bin ] , geno.mat [ , pos.bins == row.bin ] )^2
			freqs1 <- colMeans ( geno.mat [ , pos.bins == col.bin ] )
			freqs2 <- colMeans ( geno.mat [ , pos.bins == row.bin ] )
			rt.vars1 <- sqrt ( freqs1 * ( 1 - freqs1 ) )
			rt.vars2 <- sqrt ( freqs2 * ( 1 - freqs2 ) )
			var.prod.mat <- rt.vars1 %o% rt.vars2
			if ( k == j ) {
				temp <- as.matrix ( D.stat.sq / var.prod.mat )
				temp_num <-  as.matrix ( D.stat.sq )
				temp_denom <- as.matrix ( var.prod.mat )
				r.sq [ k, j ] <- mean ( temp [ col ( temp ) != row ( temp ) ] )
				sig.sq.d [ k , j ] <- mean ( temp_num [ col ( temp_num ) != row ( temp_num ) ] ) / mean ( temp_denom [ col ( temp_num) != row ( temp_num ) ] )
				
			} else {
				temp <- as.matrix ( D.stat.sq / var.prod.mat )
				r.sq [ k, j ] <- mean ( temp )
				sig.sq.d [ k, j ] <- mean ( as.matrix ( D.stat.sq ) ) / mean ( as.matrix ( var.prod.mat ) )
			}		
		}
	}
	save ( r.sq , file = paste ( "Output/rsq" , tail(strsplit(file.name,"a")[[1]],1) , ".Robj" , sep = "" ) )
	save ( sig.sq.d , file = paste ( "Output/sig_sq" , tail(strsplit(file.name,"a")[[1]],1) , ".Robj" , sep = "" ) )
}
LDSimCalc ( args [ 1 ] , as.numeric ( args [ 2 ] ) )