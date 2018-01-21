##This function was taken from the dtw package
##Ref: Toni Giorgino (2009). Computing and Visualizing Dynamic Time Warping Alignments in R: The dtw Package. Journal of Statistical Software, 31(7), 1-24, doi:10.18637/jss.v031.i07.

##
## Document in dtwWindowingFunctions.Rd
##




## no warping window: no restrictions

`noWindow` <-
function(iw,jw,...) {
  return(TRUE);
}



## A band around the diagonal. The band includes the diagonal +-
## window.size.

`sakoeChibaWindow` <-
function(iw,jw,window.size,...) {
  return(abs(jw-iw)<=window.size);
}


## A band around the diagonal. The band includes the segment
## connecting [1,1] to [query.size,reference.size] window.size,
## measured along the second index (columns)

`slantedBandWindow` <-
function(iw,jw,query.size,reference.size,window.size,...) {
  diagj<-(iw*reference.size/query.size);
  return(abs(jw-diagj)<=window.size);
}



## "Itakura" parallelogram: see documentation.
##

`itakuraWindow` <-
function(iw,jw,query.size,reference.size,...) {
	## Shorthands
  	n<-query.size;
  	m<-reference.size;

	ok<- 	(jw <  2*iw) &
 		(iw <= 2*jw) &
		(iw >= n-1-2*(m-jw)) &
		(jw >  m-1-2*(n-iw)) ;

	return(ok);
}


## Plot a sample of the windowing function

dtwWindow.plot <- function(fun,query.size=200,reference.size=220,...) {
	n<-query.size;
  	m<-reference.size;

	mm<-matrix(0,n,m);
	mm[fun(row(mm),col(mm),query.size=n,reference.size=m,...)]<-1;

	image(  z=mm,
		x=1:n,
		y=1:m,
		xlab=sprintf("Query: samples 1..%d",n),
		ylab=sprintf("Reference: samples 1..%d",m)
	 );
}

