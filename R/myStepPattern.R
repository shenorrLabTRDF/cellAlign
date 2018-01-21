##This function was taken from the dtw package
##Ref: Toni Giorgino (2009). Computing and Visualizing Dynamic Time Warping Alignments in R: The dtw Package. Journal of Statistical Software, 31(7), 1-24, doi:10.18637/jss.v031.i07.


## For pre-defined step patterns see below.


#############################
## Methods for accessing and creating step.patterns
## TODO: validate norm

stepPattern <- function(v,norm=NA) {
  obj <- NULL;
  if(is.vector(v)) {
    obj <- matrix(v,ncol=4,byrow=TRUE);
  } else if(is.matrix(v)) {
    obj <- v;
  } else {
    stop("stepPattern constructor only supports vector or matrix");
  }
  class(obj)<-"stepPattern";
  attr(obj,"npat") <- max(obj[,1]);
  attr(obj,"norm") <- norm;
  return(obj);
}


is.stepPattern <- function(x) {
  return(inherits(x,"stepPattern"));
}



## Transpose - exchange role of query and reference
t.stepPattern <- function(x) {

  # exchange dx <-> dy
  tsp <- x[,c(1,3,2,4)];
  tsp <- stepPattern(tsp);

  # fix normalization, if available
  on <- attr(x,"norm");
  if(! is.na(on) ) {
    if(on == "N") {
      attr(tsp,"norm") <- "M";
    } else if(on == "M") {
      attr(tsp,"norm") <- "N";
    }
  }

  return(tsp);
}


## plot the step pattern

plot.stepPattern <- function(x,...) {
  pats <- unique(x[,1]);                #list of patterns
  xr <- max(x[,2]);
  yr <- max(x[,3]);

  #for weight labels
  fudge <- c(-.5,1.2);
  alpha <- .5;                          # 1 start, 0 end

  ## dummy plot to fix the plot limits
  plot(-x[,2],-x[,3],type="n",
       xlab="Query index",ylab="Reference index",
       asp=1,lab=c(xr+1,yr+1,1),
       ax=FALSE,
       ...);

  for(i in pats) {
    ss <- x[,1]==i;
    lines(-x[ss,2],-x[ss,3],type="o", ...);

    if(sum(ss)==1) {
      next;
    }


    xh <- alpha*head(x[ss,2],-1) + (1-alpha)*x[ss,2][-1];
    yh <- alpha*head(x[ss,3],-1) + (1-alpha)*x[ss,3][-1];

    text(-xh,-yh,
         labels=round(x[ss,4][-1],2),
         adj=fudge,
         ...);
  }

  axis(1,at=c(-xr:0), ...)
  axis(2,at=c(-yr:0), ...)

  endpts <- x[,4]==-1;
  points(-x[endpts,2],-x[endpts,3],pch=16, ...);
}




## pretty-print the matrix meaning,
## so it will not be as write-only as now

print.stepPattern <-function(x,...) {

  step.pattern<-x;                      # for clarity
  np<-max(step.pattern[,1]);            #no. of patterns

  head<-"g[i,j] = min(\n";
  body<-"";

  ## cycle over available step patterns
  for(p in 1:np) {
    steps<-.extractpattern(step.pattern,p);
    ns<-dim(steps)[1];

    ## restore row order
    steps<-matrix(steps[ns:1,],ncol=3); # enforce a matrix

    ## cycle over steps s in the current pattern p
    for(s in 1:ns) {
      di<-steps[s,1];                   # delta in query
      dj<-steps[s,2];                   # delta in templ
      cc<-steps[s,3];                   # step cost multiplier

      ## make pretty-printable negative increments
      dis<-ifelse(di==0,"",-di);        # 4 -> -4; 0 -> .
      djs<-ifelse(dj==0,"",-dj);        #  0 maps to empty string

      ## cell origin, as coordinate pair
      dijs<-sprintf("i%2s,j%2s",dis,djs);

      if(cc==-1) {                      # g
        gs<-sprintf("   g[%s]",dijs);
        body<-paste(body,gs);
      } else {
        ## prettyprint step cost multiplier in ccs:  1 -> .; 2 -> 2 *
        ccs<-ifelse(cc==1,"    ",sprintf("%2.2g *",cc));
        ds<-sprintf("+%s d[%s]",ccs,dijs);
        body<-paste(body,ds);
      }

    }

    body<-paste(body,",\n",s="");
  }

  tail<-")\n\n";

  norm <- attr(x,"norm");
  ntxt <- sprintf("Normalization hint: %s\n",norm);

  rv<-paste(head,body,tail,ntxt);

  cat("Step pattern recursion:\n");
  cat(rv);

}



## TODO: sanity check on the step pattern definition

.checkpattern <- function(sp) {
  ## must have 4 x n elements
  ## all integers
  ## first column in ascending order from 1, no missing steps
  ## 2nd, 3rd row non-negative
  ## 4th: first  for each step is -1
}


# Auxiliary function to easily map pattern -> delta

.mkDirDeltas <- function(dir) {
  m1 <- dir[ dir[,4]==-1, ,drop=FALSE ];
  m1 <- m1[,-4];
  m1 <- m1[,-1];
  return(m1);
}


## Extract rows belonging to pattern no. sn
## with first element stripped
## in reverse order

.extractpattern <- function(sp,sn) {
	sbs<-sp[,1]==sn;	# pick only rows beginning by sn
	spl<-sp[sbs,-1,drop=FALSE];
                                # of those: take only column Di, Dj, cost
                                # (drop first - pattern no. column)

	nr<-nrow(spl);	# how many are left
        spl<-spl[nr:1,,drop=FALSE];	# invert row order

	return(spl);
}





##################################################
##################################################


## Utility inner functions to manipulate
## step patterns. Could be implemented as
## a grammar, a'la ggplot2

.Pnew <- function(p,subt,smoo) {
  sp <- list();
  sp$i <- 0;
  sp$j <- 0;
  sp$p <- p;
  sp$subt <- subt;
  sp$smoo <- smoo;
  return(sp);
}

.Pstep <- function(sp,di,dj) {
  sp$i <- c(sp$i,di);
  sp$j <- c(sp$j,dj);
  return(sp);
}

.Pend <- function(sp,subt,smoo) {
  sp$si <- cumsum(sp$i);
  sp$sj <- cumsum(sp$j);
  sp$ni <- max(sp$si)-sp$si;
  sp$nj <- max(sp$sj)-sp$sj;

  w <- NULL;

  # smallest of i,j jumps
  if(sp$subt=="a") {
    w <- pmin(sp$i,sp$j);
  } else if(sp$subt=="b") {
    # largest of Di, Dj
    w <- pmax(sp$i,sp$j);
  } else if(sp$subt=="c") {
    # Di exactly
    w <- sp$i;
  } else if(sp$subt=="d") {
    # Di+Dj
    w <- sp$i+sp$j;
  } else {
    stop("Unsupported subtype");
  }

                                        # drop first element in w
  w <- w[-1];

  if(sp$smoo)
    w <- rep(mean(w),length(w));

                                        # prepend -1
  w <- c(-1,w);
  sp$w <- w;

  return(sp);
}

.PtoMx <- function(sp) {
  nr <- length(sp$i);
  mx <- matrix(nrow=nr,ncol=4)
  mx[,1] <- sp$p;
  mx[,2] <- sp$ni;
  mx[,3] <- sp$nj;
  mx[,4] <- sp$w;
  return(mx);
}



rabinerJuangStepPattern <- function(type,slope.weighting="d",smoothed=FALSE) {

  sw <- slope.weighting;
  sm <- smoothed;

  ## Actually build the step
  r <- switch(type,
              .RJtypeI(sw,sm),
              .RJtypeII(sw,sm),
              .RJtypeIII(sw,sm),
              .RJtypeIV(sw,sm),
              .RJtypeV(sw,sm),
              .RJtypeVI(sw,sm),
              .RJtypeVII(sw,sm)
              );

  norm <- NA;
  if(sw=="c") {
    norm <- "N";
  } else if(sw=="d") {
    norm <- "N+M";
  }

  # brain-damaged legacy
  rv <- as.vector(t(r));
  rs <- stepPattern(rv);
  attr(rs,"norm") <- norm;
  attr(rs,"call") <- match.call();
  return(rs);
}



.RJtypeI <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,0,1)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  return(rbind(m1,m2,m3));
}


.RJtypeII <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,0,1)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  return(rbind(m1,m2,m3));
}



.RJtypeIII <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,2,1)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,1,2)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  return(rbind(m1,m2,m3));
}



.RJtypeIV <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,2)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  t <- .Pnew(4,s,m)
  t <- .Pstep(t,1,2)
  t <- .Pend(t);
  m4 <- .PtoMx(t)

  return(rbind(m1,m2,m3,m4));
}



.RJtypeV <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  t <- .Pnew(4,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,0,1)
  t <- .Pend(t);
  m4 <- .PtoMx(t)

  t <- .Pnew(5,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,0,1)
  t <- .Pstep(t,0,1)
  t <- .Pend(t);
  m5 <- .PtoMx(t)

  return(rbind(m1,m2,m3,m4,m5));
}



.RJtypeVI <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,0,1)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  return(rbind(m1,m2,m3));
}




.RJtypeVII <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,2)
  t <- .Pstep(t,1,0)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,1,3)
  t <- .Pstep(t,1,0)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  t <- .Pnew(4,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m4 <- .PtoMx(t)

  t <- .Pnew(5,s,m)
  t <- .Pstep(t,1,2)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m5 <- .PtoMx(t)

  t <- .Pnew(6,s,m)
  t <- .Pstep(t,1,3)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m6 <- .PtoMx(t);

  t <- .Pnew(7,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m7 <- .PtoMx(t)

  t <- .Pnew(8,s,m)
  t <- .Pstep(t,1,2)
  t <- .Pend(t);
  m8 <- .PtoMx(t)

  t <- .Pnew(9,s,m)
  t <- .Pstep(t,1,3)
  t <- .Pend(t);
  m9 <- .PtoMx(t)

  return(rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9));
}














##################################################
##################################################


##
## Various step patterns, defined as internal variables
##
## First column: enumerates step patterns.
## Second   	 step in query index
## Third	 step in reference index
## Fourth	 weight if positive, or -1 if starting point
##
## For \cite{} see dtw.bib in the package
##



## Widely-known variants

## White-Neely symmetric (default)
## aka Quasi-symmetric \cite{White1976}
## normalization: no (N+M?)
symmetric1 <- stepPattern(c(
                            1,1,1,-1,
                            1,0,0,1,
                            2,0,1,-1,
                            2,0,0,1,
                            3,1,0,-1,
                            3,0,0,1
                            ));


## Normal symmetric
## normalization: N+M
symmetric2 <- stepPattern(c(
                            1,1,1,-1,
                            1,0,0,2,
                            2,0,1,-1,
                            2,0,0,1,
                            3,1,0,-1,
                            3,0,0,1
                            ),"N+M");


## classic asymmetric pattern: max slope 2, min slope 0
## normalization: N
asymmetric <-  stepPattern(c(
                             1,1,0,-1,
                             1,0,0,1,
                             2,1,1,-1,
                             2,0,0,1,
                             3,1,2,-1,
                             3,0,0,1
                           ),"N");



## normalization: max[N,M]
## note: local distance matrix is 1-d
## \cite{Velichko}
.symmetricVelichkoZagoruyko <- stepPattern(c(
		1, 0, 1, -1,
		2, 1, 1, -1,
		2, 0, 0, -1.001,
		3, 1, 0, -1 ));



## Itakura slope-limited asymmetric \cite{Itakura1975}
## Max slope: 2; min slope: 1/2
## normalization: N
.asymmetricItakura <-  stepPattern(c(
                        1, 1, 2, -1,
			1, 0, 0, 1,
			2, 1, 1, -1,
			2, 0, 0, 1,
			3, 2, 1, -1,
			3, 1, 0, 1,
			3, 0, 0, 1,
			4, 2, 2, -1,
			4, 1, 0, 1,
			4, 0, 0, 1
                       ));







#############################
## Slope-limited versions
##
## Taken from Table I, page 47 of "Dynamic programming algorithm
## optimization for spoken word recognition," Acoustics, Speech, and
## Signal Processing, vol.26, no.1, pp. 43-49, Feb 1978 URL:
## http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055
##
## Mostly unchecked



## Row P=0
symmetricP0 <- symmetric2;

## normalization: N ?
asymmetricP0 <- stepPattern(c(
                                  1,0,1,-1,
                                  1,0,0,0,
                                  2,1,1,-1,
                                  2,0,0,1,
                                  3,1,0,-1,
                                  3,0,0,1
                                ),"N");


## alternative implementation
.asymmetricP0b <- stepPattern(c(
                                  1,0,1,-1,
                                  2,1,1,-1,
                                  2,0,0,1,
                                  3,1,0,-1,
                                  3,0,0,1
                                ),"N");



## Row P=1/2
symmetricP05 <-  stepPattern(c(
                        1  ,  1, 3 , -1,
                        1  ,  0, 2 ,  2,
                        1  ,  0, 1 ,  1,
                        1  ,  0, 0 ,  1,
                        2  ,  1, 2 , -1,
                        2  ,  0, 1 ,  2,
                        2  ,  0, 0 ,  1,
                        3  ,  1, 1 , -1,
                        3  ,  0, 0 ,  2,
                        4  ,  2, 1 , -1,
                        4  ,  1, 0 ,  2,
                        4  ,  0, 0 ,  1,
                        5  ,  3, 1 , -1,
                        5  ,  2, 0 ,  2,
                        5  ,  1, 0 ,  1,
                        5  ,  0, 0 ,  1
                               ),"N+M");

asymmetricP05 <-  stepPattern(c(
                        1  , 1 , 3 , -1,
                        1  , 0 , 2 ,1/3,
                        1  , 0 , 1 ,1/3,
                        1  , 0 , 0 ,1/3,
                        2  , 1 , 2 , -1,
                        2  , 0 , 1 , .5,
                        2  , 0 , 0 , .5,
                        3  , 1 , 1 , -1,
                        3  , 0 , 0 , 1 ,
                        4  , 2 , 1 , -1,
                        4  , 1 , 0 , 1 ,
                        4  , 0 , 0 , 1 ,
                        5  , 3 , 1 , -1,
                        5  , 2 , 0 , 1 ,
                        5  , 1 , 0 , 1 ,
                        5  , 0 , 0 , 1
                               ),"N");



## Row P=1
## Implementation of Sakoe's P=1, Symmetric algorithm

symmetricP1 <- stepPattern(c(
                              1,1,2,-1,	# First branch: g(i-1,j-2)+
                              1,0,1,2,	#            + 2d(i  ,j-1)
                              1,0,0,1,	#            +  d(i  ,j)
                              2,1,1,-1,	# Second branch: g(i-1,j-1)+
                              2,0,0,2,	#              +2d(i,  j)
                              3,2,1,-1,	# Third branch: g(i-2,j-1)+
                              3,1,0,2,	#            + 2d(i-1,j)
                              3,0,0,1	#            +  d(  i,j)
                        ),"N+M");

asymmetricP1 <- stepPattern(c(
                              1, 1 , 2 , -1 ,
                              1, 0 , 1 , .5 ,
                              1, 0 , 0 , .5 ,
                              2, 1 , 1 , -1 ,
                              2, 0 , 0 ,  1 ,
                              3, 2 , 1 , -1 ,
                              3, 1 , 0 ,  1 ,
                              3, 0 , 0 ,  1
                              ),"N");


## Row P=2
symmetricP2 <- stepPattern(c(
	1, 2, 3, -1,
	1, 1, 2, 2,
	1, 0, 1, 2,
	1, 0, 0, 1,
	2, 1, 1, -1,
	2, 0, 0, 2,
	3, 3, 2, -1,
	3, 2, 1, 2,
	3, 1, 0, 2,
	3, 0, 0, 1
),"N+M");

asymmetricP2 <- stepPattern(c(
	1, 2 , 3  , -1,
	1, 1 , 2  ,2/3,
	1, 0 , 1  ,2/3,
	1, 0 , 0  ,2/3,
	2, 1 , 1  ,-1 ,
	2, 0 , 0  ,1  ,
	3, 3 , 2  ,-1 ,
	3, 2 , 1  ,1  ,
	3, 1 , 0  ,1  ,
	3, 0 , 0  ,1
),"N");






################################
## Taken from Table III, page 49.
## Four varieties of DP-algorithm compared

## 1st row:  asymmetric

## 2nd row:  symmetricVelichkoZagoruyko

## 3rd row:  symmetric1

## 4th row:  asymmetricItakura




#############################
## Classified according to Rabiner
##
## Taken from chapter 2, Myers' thesis [4]. Letter is
## the weighting function:
##
##      rule       norm   unbiased
##   a  min step   ~N     NO
##   b  max step   ~N     NO
##   c  x step     N      YES
##   d  x+y step   N+M    YES
##
## Mostly unchecked

# R-Myers     R-Juang
# type I      type II
# type II     type III
# type III    type IV
# type IV     type VII


typeIa <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0,  1,
                         1, 0, 0,  0,
                         2, 1, 1, -1,
                         2, 0, 0,  1,
                         3, 1, 2, -1,
                         3, 0, 1,  1,
                         3, 0, 0,  0
 ));

typeIb <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0,  1,
                         1, 0, 0,  1,
                         2, 1, 1, -1,
                         2, 0, 0,  1,
                         3, 1, 2, -1,
                         3, 0, 1,  1,
                         3, 0, 0,  1
 ));

typeIc <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0,  1,
                         1, 0, 0,  1,
                         2, 1, 1, -1,
                         2, 0, 0,  1,
                         3, 1, 2, -1,
                         3, 0, 1,  1,
                         3, 0, 0,  0
 ),"N");

typeId <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0,  2,
                         1, 0, 0,  1,
                         2, 1, 1, -1,
                         2, 0, 0,  2,
                         3, 1, 2, -1,
                         3, 0, 1,  2,
                         3, 0, 0,  1
 ),"N+M");

## ----------
## smoothed variants of above

typeIas <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0, .5,
                         1, 0, 0, .5,
                         2, 1, 1, -1,
                         2, 0, 0,  1,
                         3, 1, 2, -1,
                         3, 0, 1, .5,
                         3, 0, 0, .5
 ));


typeIbs <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0,  1,
                         1, 0, 0,  1,
                         2, 1, 1, -1,
                         2, 0, 0,  1,
                         3, 1, 2, -1,
                         3, 0, 1,  1,
                         3, 0, 0,  1
 ));


typeIcs <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0,  1,
                         1, 0, 0,  1,
                         2, 1, 1, -1,
                         2, 0, 0,  1,
                         3, 1, 2, -1,
                         3, 0, 1, .5,
                         3, 0, 0, .5
 ),"N");


typeIds <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0, 1.5,
                         1, 0, 0, 1.5,
                         2, 1, 1, -1,
                         2, 0, 0,  2,
                         3, 1, 2, -1,
                         3, 0, 1, 1.5,
                         3, 0, 0, 1.5
 ),"N+M");






## ----------

typeIIa <- stepPattern(c(
                        1,  1,  1, -1,
                        1,  0,  0, 1,
                        2,  1,  2, -1,
                        2,  0,  0, 1,
                        3,  2,  1, -1,
                        3,  0,  0, 1
                        ));

typeIIb <- stepPattern(c(
                        1,  1,  1, -1,
                        1,  0,  0, 1,
                        2,  1,  2, -1,
                        2,  0,  0, 2,
                        3,  2,  1, -1,
                        3,  0,  0, 2
                        ));

typeIIc <- stepPattern(c(
                        1,  1,  1, -1,
                        1,  0,  0, 1,
                        2,  1,  2, -1,
                        2,  0,  0, 1,
                        3,  2,  1, -1,
                        3,  0,  0, 2
                        ),"N");

typeIId <- stepPattern(c(
                        1,  1,  1, -1,
                        1,  0,  0, 2,
                        2,  1,  2, -1,
                        2,  0,  0, 3,
                        3,  2,  1, -1,
                        3,  0,  0, 3
                        ),"N+M");

## ----------

## Rabiner [3] discusses why this is not equivalent to Itakura's

typeIIIc <-  stepPattern(c(
                        1, 1, 2, -1,
			1, 0, 0, 1,
			2, 1, 1, -1,
			2, 0, 0, 1,
			3, 2, 1, -1,
			3, 1, 0, 1,
			3, 0, 0, 1,
			4, 2, 2, -1,
			4, 1, 0, 1,
			4, 0, 0, 1
                       ),"N");



## ----------

## numbers follow as production rules in fig 2.16

typeIVc <-  stepPattern(c(
                          1,  1,  1,  -1,
                          1,  0,  0,   1,
                          2,  1,  2,  -1,
                          2,  0,  0,   1,
                          3,  1,  3,  -1,
                          3,  0,  0,   1,
                          4,  2,  1,  -1,
                          4,  1,  0,   1,
                          4,  0,  0,   1,
                          5,  2,  2,  -1,
                          5,  1,  0,   1,
                          5,  0,  0,   1,
                          6,  2,  3,  -1,
                          6,  1,  0,   1,
                          6,  0,  0,   1,
                          7,  3,  1,  -1,
                          7,  2,  0,   1,
                          7,  1,  0,   1,
                          7,  0,  0,   1,
                          8,  3,  2,  -1,
                          8,  2,  0,   1,
                          8,  1,  0,   1,
                          8,  0,  0,   1,
                          9,  3,  3,  -1,
                          9,  2,  0,   1,
                          9,  1,  0,   1,
                          9,  0,  0,   1
 ),"N");






#############################
##
## Mori's asymmetric step-constrained pattern. Normalized in the
## reference length.
##
## Mori, A.; Uchida, S.; Kurazume, R.; Taniguchi, R.; Hasegawa, T. &
## Sakoe, H. Early Recognition and Prediction of Gestures Proc. 18th
## International Conference on Pattern Recognition ICPR 2006, 2006, 3,
## 560-563
##

mori2006 <-  stepPattern(c(
                           1, 2, 1, -1,
                           1, 1, 0,  2,
                           1, 0, 0,  1,
                           2, 1, 1, -1,
                           2, 0, 0,  3,
                           3, 1, 2, -1,
                           3, 0, 1,  3,
                           3, 0, 0,  3
 ),"M");


## Completely unflexible: fixed slope 1. Only makes sense with
## open.begin and open.end
rigid <- stepPattern(c(1,1,1,-1,
                       1,0,0,1  ),"N")
