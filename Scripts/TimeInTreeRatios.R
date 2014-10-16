TimeInTreeRatio <- function ( f , s , n , Ne ) {
	( s * 2 * Ne * f * sum ( 1 / seq ( 1 , n - 1 ) ) ) / ( n*log ( (Ne - 1)/f - Ne + 1  ) )
}

FixRatioForf <- function ( f , y , n , Ne ) {
	y*n*log ( (Ne - 1)/f - Ne + 1  ) / ( 2 * Ne * f * sum ( 1 / seq ( 1 , n - 1 ) ) )
}
TimeInTreeFrac <- function ( f , s , n , Ne ) {
	( 2 * Ne * f * sum ( 1 / seq ( 1 , n - 1 ) ) ) / ( n*log ( (Ne - 1)/f - Ne + 1  )/s + 2 * Ne * f * sum ( 1 / seq ( 1 , n - 1 ) ) )
}

FixRatioForfFrac <- function ( f , y , n , Ne ) {
	y*n*log ( (Ne - 1)/f - Ne + 1  ) / ( 2 * Ne * f * sum ( 1 / seq ( 1 , n - 1 ) ) * ( 1 - y ) )
}
image.scale <- function(z, zlim, col = heat.colors(12), breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...) {
	#recover()
	if(!missing(breaks)){
 		if(length(breaks) != (length(col)+1)) {
			stop("must have one more break than colour")
		}
 	}
	if(missing(breaks) & !missing(zlim)) {
		breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
	}
 	if ( missing(breaks) & missing(zlim) ) {
  		zlim <- range(z, na.rm=TRUE)
  		zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  		zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  		breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 	}
	poly <- vector(mode="list", length(col))
 	for(i in seq(poly)){
  		poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 	}
 	# xaxt <- ifelse(horiz, "s", "n")
 	# yaxt <- ifelse(horiz, "n", "s")
 	if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
 	if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
 	if(missing(xlim)) xlim=XLIM
 	if(missing(ylim)) ylim=YLIM
 	plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt="n", yaxt="n", xaxs="i", yaxs="i",ylab="", ...)  
 	for(i in seq(poly)) {
  		if(horiz) {
   			polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
  		}
  		if(!horiz){
  			polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  		}
 	}
}


setwd("Documents/Academics/StandingSweep/Figures/")
TtotRatioPlot <- function ( my.s , my.f , my.N , my.n ) {
	recover()
	#my.f <- seq ( 0.0001 , 0.060 , by = 0.0001)
	#my.s <- seq ( 0.0001 , 0.060 , by = 0.0001)
	my.s.f <- expand.grid ( f = my.f , s = my.s )
	## get ratio of time in tree for many different combinations of s and f
	# y <- numeric ( nrow ( my.s.f ) )
	# for ( i in 1:nrow ( my.s.f ) ) {
		
		# y [ i ] <- TimeInTreeRatio ( my.s.f [ i, 1 ] , my.s.f [ i , 2 ] , my.n , my.N )
		
	# }
	# ylog <- log ( y )
	# ylog.mat <- matrix ( ylog , nrow = length ( my.f ) , ncol = length ( my.s ) )
	
	
	# # Get selection coefficients for different time in tree ratios
	# s.1 <- numeric ( length ( my.f ) )
	# s.0.1 <- numeric ( length ( my.f ) )
	# s.0.01 <- numeric ( length ( my.f ) )
	# s.10 <- numeric ( length ( my.f ) )
	# s.100 <- numeric ( length ( my.f ) )
	# for ( i in 1 : length ( my.f ) ) {
		
		# s.1 [ i ] <- FixRatioForf ( my.f [ i ] , 1 , my.n , my.N )
		# s.0.1 [ i ] <- FixRatioForf ( my.f [ i ] , 0.1 , my.n , my.N )	
		# s.10 [ i ] <- FixRatioForf ( my.f [ i ] , 10 , my.n , my.N )	
		# s.0.01 [ i ] <- FixRatioForf ( my.f [ i ] , 0.01 , my.n , my.N )	
		# s.100 [ i ] <- FixRatioForf ( my.f [ i ] , 100 , my.n , my.N )	
	# }
	
	
	# quartz(height=5,width=7)
		# par(mar=c(0,0,0,0),oma=c(2,1,0,0.5))
		# layout(matrix(c(1,2,3,4),2,2),heights=c(1.2,6) , widths = c ( 5.57 , 1.43))
		# par(mar=c(2,2,1.5,0.5))
		# image.scale ( ylog.mat , col = rainbow ( 2000  , end = 0.9 ) )
		# axis ( 1 , log ( 1/rev ( c ( 0.01 , 0.1 , 1 , 10 , 100 , 1000 , 1000 , 10000 , 100000 , 1000000 ) ) ) , 1/rev ( c ( 0.01 , 0.1 , 1 , 10 , 100 , 1000 , 1000  , 10000 , 100000 , 1000000 ) ) )
		# mtext ( expression ( T[paste ( "tot" , "," , "stand" , sep = "" )]/T[paste ( "tot",",","sweep" , sep = "" )] ), 1 , 2.3  )
		
		
		# image ( y = my.s , x = my.f , z=ylog.mat , col = rainbow ( 2000 , end = 0.9 ), bty = "n" , xlab = "" , ylab = "", xaxt = "n" , yaxt = "n" , xlim = range ( c ( 0 , my.s ) ) , ylim = range ( c ( 0 , my.f ) ) )
		# axis ( 1 , seq ( 0 , 0.1 , 0.02 ) )
		# axis ( 2 , seq ( 0 , 0.1 , 0.02 ) )
		# mtext ( "Selection Coefficient" ,  1 , 2 )
		# mtext ( "Frequency" ,  2 , 2 )
		# lines ( s.1 , my.f , lty = 2 )
		# lines ( s.0.1 , my.f , lty = 3 )
		# lines ( s.10 , my.f , lty = 3 )
		# lines ( s.0.01 , my.f , lty = 4 )
		# lines ( s.100 , my.f , lty = 4 )
		
		# plot ( 1,1,type = "n" , bty = "n" , xaxt = "n" , yaxt = "n" )
		
		# plot ( 1 , 1,type = "n" , bty = "n" ,  xaxt = "n" , yaxt = "n" )
		# text ( x = 0.9 , y = 1.3 , sprintf ("%s %d" , "N = ",  my.N ) )
		# text ( x = 0.83 , y = 1.25 , sprintf ("%s %d" , "n = ",  my.n ) )
		# legend ( 0.6 , 1.2 , legend = c ( 1 , 10 , 100 ) , lty = c ( 2 , 3 , 4 ) , bty = "n")
	# #dev.off()
	
	
	
	my.s.f <- expand.grid ( f = my.f , s = my.s )
	# frac of time spent in standing
	x <- numeric ( nrow ( my.s.f ) )
	for ( i in 1:nrow ( my.s.f ) ) {
		
		x [ i ] <- TimeInTreeFrac ( my.s.f [ i, 1 ] , my.s.f [ i , 2 ] , my.n , my.N )
		
	}
	#ylog <- log ( y )
	x.mat <- matrix ( x , nrow = length ( my.f ) , ncol = length ( my.s ) )
	s.1 <- numeric ( length ( my.f ) )
	s.0.1 <- numeric ( length ( my.f ) )
	s.0.01 <- numeric ( length ( my.f ) )
	s.0.9 <- numeric ( length ( my.f ) )
	s.0.99 <- numeric ( length ( my.f ) )
	for ( i in 1 : length ( my.f ) ) {
		
		s.1 [ i ] <- FixRatioForfFrac ( my.f [ i ] , 0.5 , my.n , my.N )
		s.0.1 [ i ] <- FixRatioForfFrac ( my.f [ i ] , 0.1 , my.n , my.N )	
		s.0.9 [ i ] <- FixRatioForfFrac ( my.f [ i ] , 0.9 , my.n , my.N )	
		s.0.01 [ i ] <- FixRatioForfFrac ( my.f [ i ] , 0.01 , my.n , my.N )	
		s.0.99 [ i ] <- FixRatioForfFrac ( my.f [ i ] , 0.99 , my.n , my.N )	
	}
	quartz(height=5,width=7)
		par(mar=c(0,0,0,0),oma=c(2,1,0,0.5))
		layout(matrix(c(1,2,3,4),2,2),heights=c(1.2,6) , widths = c ( 5.57 , 1.43))
		par(mar=c(2,2,1.5,0.5))
		image.scale ( x.mat , col = rainbow ( 2000  , end = 0.9 ) )
		axis ( 1 , seq ( 0 , 1 , 0.1 )  )
		mtext ( expression ( T[paste ( "tot" , "," , "stand" , sep = "" )]/T[paste ( "tot",",","sweep+stand" , sep = "" )] ), 1 , 2.3  )
		image ( y = my.s , x = my.f , z=x.mat , col = rainbow ( 2000 , end = 0.9 ), bty = "n" , xlab = "" , ylab = "", xaxt = "n" , yaxt = "n" , xlim = range ( c ( 0 , my.s ) ) , ylim = range ( c ( 0 , my.f ) ) )
		axis ( 1 , seq ( 0 , 0.1 , 0.02 ) )
		axis ( 2 , seq ( 0 , 0.1 , 0.02 ) )
		mtext ( "Selection Coefficient" ,  1 , 2 )
		mtext ( "Frequency" ,  2 , 2 )
		lines ( s.1 , my.f , lty = 2 )
		lines ( s.0.1 , my.f , lty = 3 )
		lines ( s.0.9 , my.f , lty = 3 )
		lines ( s.0.01 , my.f , lty = 4 )
		lines ( s.0.99 , my.f , lty = 4 )
		plot ( 1,1,type = "n" , bty = "n" , xaxt = "n" , yaxt = "n" )
		plot ( 1 , 1,type = "n" , bty = "n" ,  xaxt = "n" , yaxt = "n" )
		text ( x = 0.9 , y = 1.3 , sprintf ("%s %d" , "N = ",  my.N ) )
		text ( x = 0.83 , y = 1.25 , sprintf ("%s %d" , "n = ",  my.n ) )
		legend ( 0.6 , 1.2 , legend = c ( 1 , 10 , 100 ) , lty = c ( 2 , 3 , 4 ) , bty = "n")

	
	
	
	#setEPS()
	#pdf ( file = paste ( save.loc , ".pdf",sep = "") ,height=5,width=7)
	

}

TtotRatioPlot ( seq ( 0.0001 , 0.1 , by = 0.0001) , seq ( 0.0001 , 0.1 , by = 0.0001) , 20000 , 20 )


quartz(height=5,width=7)
par(mar=c(0,0,0,0),oma=c(2,1,0,0.5))
layout(matrix(c(1,2,3,4),2,2),heights=c(1.2,6) , widths = c ( 5.57 , 1.43))
par(mar=c(2,2,1.5,0.5))
image.scale ( ylog.mat , col = rainbow ( 2000 ) )
axis ( 1 , log ( 1/rev ( c ( 0.25 , 1 , 10 , 100 , 1000 , 1000 , 10000 , 100000 ) ) ) , 1/rev ( c ( 0.25 , 1 , 10 , 100 , 1000 , 1000  , 10000 , 100000 ) ) )
mtext ( expression ( T[paste ( "tot" , "," , "stand" , sep = "" )]/T[paste ( "tot",",","sweep" , sep = "" )] ), 1 , 2.3  )


image ( y = my.s , x = my.f , z=ylog.mat , col = rainbow ( 2000 ), bty = "n" , xlab = "" , ylab = "", xaxt = "n" , yaxt = "n" , xlim = range ( c ( 0 , my.s ) ) , ylim = range ( c ( 0 , my.f ) ) )
axis ( 1 , seq ( 0 , 0.1 , 0.02 ) )
axis ( 2 , seq ( 0 , 0.1 , 0.02 ) )
mtext ( "Selection Coefficient" ,  1 , 2 )
mtext ( "Frequency" ,  2 , 2 )
lines ( s.1 , my.f , lty = 2 )
lines ( s.0.25 , my.f , lty = 3 )
lines ( s.4 , my.f , lty = 3 )
lines ( s.0.1 , my.f , lty = 4 )
lines ( s.10 , my.f , lty = 4 )

plot ( 1,1,type = "n" , bty = "n" , xaxt = "n" , yaxt = "n" )

plot ( 1 , 1,type = "n" , bty = "n" ,  xaxt = "n" , yaxt = "n" )
text ( x = 0.9 , y = 1.3 , "N = 20000")
text ( x = 0.83 , y = 1.25 , "n = 20" )
legend ( 0.6 , 1.2 , legend = c ( 1 , expression ( 1/4 ) , expression ( 1/10 ) ) , lty = c ( 2 , 3 , 4 ) , bty = "n")








