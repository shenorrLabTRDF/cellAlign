###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino,gmail.com>               #
#       Istituto di Ingegneria Biomedica (ISIB-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id: globalCostMatrix.R 343 2013-12-11 11:04:41Z tonig $
#                                                             #
###############################################################


########################################
## Compute the cost matrix from a local distance matrix

## Wrapper to the native function

`globalCostMatrix` <-
function(lm,
         step.matrix=symmetric1,
         window.function=noWindow,
         native=TRUE,
         seed=NULL,
         ...) {


  ## sanity check - be extra cautions w/ binary
  if (!is.stepPattern(step.matrix))
    stop("step.matrix is no stepMatrix object");



  # i = 1 .. n in query sequence, on first index, ie rows
  # j = 1 .. m on reference sequence, on second index, ie columns
  #   Note:  reference is usually drawn vertically, up-wise

  n <- nrow(lm);
  m <- ncol(lm);


  # number of individual steps (counting all patterns)
  nsteps<-dim(step.matrix)[1];


  # clear the cost and step matrix
  # these will be the outputs of the binary
  # for  cm use  seed if given
  if(!is.null(seed)) {
    cm <- seed;
  } else {
    cm <- matrix(NA,nrow=n,ncol=m);
    cm[1,1] <- lm[1,1];                 # Questionable.
  }

  sm <- matrix(NA,nrow=n,ncol=m);


  if(is.loaded("computeCM") && native){
    ## precompute windowing
    wm <- matrix(FALSE,nrow=n,ncol=m);
    wm[window.function(row(wm),col(wm),
                       query.size=n, reference.size=m,
                       ...)]<-TRUE;

    if(FALSE) {
      ## this call could be optimized. Copies are killing perf.
      out<-.C(computeCM,
            NAOK=TRUE,
            PACKAGE="dtw",
            ## IN
            as.integer(dim(cm)),               # int *s
            as.logical(wm),                    # int *wm
            as.double(lm),                     # double *lm
            as.integer(nsteps),                # int *nstepsp
            as.double(step.matrix),            # double *dir
            ## IN+OUT
            costMatrix=as.double(cm),                 # double *cm
            ## OUT
            directionMatrix=as.integer(sm));               # int *sm

      ## Hopefully avoids a copy
      dim(out$costMatrix) <- c(n,m);
      dim(out$directionMatrix) <- c(n,m);
      warning("You should not be here");
    } else {
      storage.mode(wm) <- "logical";
      storage.mode(lm) <- "double";
      storage.mode(cm) <- "double";
      storage.mode(step.matrix) <- "double";
      out <- .Call("computeCM_Call",
                   wm,lm,cm,step.matrix);
    }

  } else {

    ####################
    ## INTERPRETED PURE-R IMPLEMENTATION
    warning("Native dtw implementation not available: using (slow) interpreted fallback");
                                        # now walk through the matrix, column-wise and row-wise,
                                        # and recursively compute the accumulated distance. Unreachable
                                        # elements are handled via NAs (removed)
    dir <- step.matrix;
    npats <- attr(dir,"npat");
    for (j in 1:m) {
      for (i in 1:n) {
        ## It is ok to window on the arrival point (?)
        if(!window.function(i,j, query.size=n, reference.size=m, ...)) { next; }

        ## Skip if already initialized
        if(!is.na(cm[i,j])) { next; }

        clist<-numeric(npats)+NA;
        for (s in 1:nsteps) {
          ## current pattern
          p<-dir[s,1];
          ## ii,jj is the cell from where potentially we could
          ## have come from.
          ii<-i-dir[s,2];                 # previous step in inp
          jj<-j-dir[s,3];                 # previous step in tpl
          if(ii>=1 && jj>=1) {            # element exists?
            cc<-dir[s,4];                 # step penalty
            if(cc == -1) {                #  -1? cumulative cost:
              clist[p]<-cm[ii,jj];	#  there must be exactly 1 per pattern
            } else {			#  a cost for
              clist[p]<-clist[p]+cc*lm[ii,jj];
            }
          }
        }


        ## no NAs in clist at this point BUT clist can be empty
        ## store in cost matrix
        minc<-which.min(clist);           # pick the least cost
        if(length(minc)>0) {          	# false if clist has all NAs
          cm[i,j]<-clist[minc];
          sm[i,j]<-minc;			# remember the pattern picked
        }
      }
    }
    out <- list(costMatrix=cm,directionMatrix=sm);
  }

  ## END PURE-R IMPLEMENTATION
  ####################

  ## At this point out$cmo and out$smo should be set
  out$stepPattern <- step.matrix;
  return(out);
}

