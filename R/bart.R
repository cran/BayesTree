bart <- function(x.train,y.train,x.test=matrix(0.0,0,0),
	sigest=NA,sigdf=3, sigquant=.90, k=2.0, ntree=200,
	ndpost=1000,nskip=100,
        printevery=100,keepevery=1,keeptrainfits=TRUE,
	numcut=100) 
{

        #check input arguments:
        if((!is.matrix(x.train)) || (typeof(x.train)!="double")) stop("argument x.train must be a double matrix")
        if((!is.matrix(x.test)) || (typeof(x.test)!="double")) stop("argument x.test must be a double matrix")
        if((!is.vector(y.train)) || (typeof(y.train)!="double")) stop("argument y.train must be a double vector")
        if(nrow(x.train) != length(y.train)) stop("number of rows in x.train must equal length of y.train")
        if((nrow(x.test) >0) && (ncol(x.test)!=ncol(x.train))) stop("input x.test must have the same number of columns as x.train")
        if((!is.na(sigest)) && (typeof(sigest)!="double")) stop("input sigest must be double")
        if((!is.na(sigest)) && (sigest<0.0)) stop("input sigest must be positive")
        if((mode(sigdf)!="numeric") || (sigdf<0)) stop("input sigdf must be a positive number")
        if((mode(printevery)!="numeric") || (printevery<0)) stop("input printevery must be a positive number")
        if((mode(keepevery)!="numeric") || (keepevery<0)) stop("input keepevery must be a positive number")
        if((mode(sigquant)!="numeric") || (sigquant<0)) stop("input sigquant must be a positive number")
        if((mode(ntree)!="numeric") || (ntree<0)) stop("input ntree must be a positive number")
        if((mode(ndpost)!="numeric") || (ndpost<0)) stop("input ndpost must be a positive number")
        if((mode(nskip)!="numeric") || (nskip<0)) stop("input nskip must be a positive number")
        if((mode(k)!="numeric") || (k<0)) stop("input k must be a positive number")
        if((mode(numcut)!="numeric") || (numcut<0)) stop("input numcut must be a positive number")
        if(typeof(keeptrainfits)  != "logical") stop("input keeptrainfits must a logical variable")

	# note: training data is split into x and y to avoid mixing up cols
	# of data.  Will run a total of ndpost+nskip steps of algorithm.

	traindata <- as.matrix(cbind(x.train,y.train))
	# scale the data
	if (1) {
   	   rg <- apply(traindata,2,range)
	   for (i in 1:ncol(traindata)) 
	      traindata[,i] <- -.5 + (traindata[,i]-rg[1,i])/(rg[2,i]-rg[1,i])
	   if (nrow(x.test))
	      for (i in 1:ncol(x.test)) 
	         x.test[,i] <- -.5 + (x.test[,i]-rg[1,i])/(rg[2,i]-rg[1,i])
	}
	# if sigest=NA, fit a lm to training data to get the value of sigest...
	# sigest is on the scale of the transformed y, so we do the lm after
	# the scaling above...
	if (is.na(sigest)){
		templm <- lm(traindata[,ncol(traindata)]~traindata[,1:(ncol(traindata)-1)])
		sigest <- summary(templm)$sigma
	} else {
		# we're assuming sigest is given on y scale, but BART code needs
		# it on the transformed scale.
		rgy <- range(y.train)
		sigest <- sigest/(rgy[2]-rgy[1])
	}

        ncskip = floor(nskip/keepevery)
        ncpost = floor(ndpost/keepevery)
        nctot = ncskip + ncpost
        totnd = keepevery*nctot

        cres = .C('mbart',as.integer(nrow(traindata)), as.integer(ncol(x.train)), as.integer(nrow(x.test)),
                   as.double(traindata[,1:(ncol(traindata)-1)]), as.double(traindata[,ncol(traindata)]),
                   as.double(x.test),
                   as.double(sigest),   as.integer(sigdf), as.double(sigquant),
                   as.double(k),
		   as.integer(ntree),      as.integer(totnd),
                   as.integer(printevery), as.integer(keepevery),  as.integer(keeptrainfits),
                   as.integer(numcut),
                   sdraw=double(nctot), 
                   trdraw=double(nrow(x.train)*nctot),
                   tedraw=double(nrow(x.test)*nctot),
                   vcdraw=integer(ncol(x.train)*nctot))
	
	# now read in the results...
	rgy <- range(y.train)
        sigma = cres$sdraw*(rgy[2]-rgy[1])
	first.sigma <- sigma[1:ncskip] # we often want the sigma draws 
	sigma <- sigma[ncskip+(1:ncpost)]

	# put sigest on the original y scale for output purposes
	sigest <- sigest*(rgy[2]-rgy[1]) 

        yhat.train <- yhat.test <- yhat.train.mean <- yhat.test.mean <- NULL
        yhat.train.quantiles <- yhat.test.quantiles <- NULL
        varcount = NULL

	if (keeptrainfits) {
                yhat.train = matrix(cres$trdraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
		yhat.train <- (rgy[2]-rgy[1])*(yhat.train+.5) + rgy[1]
		yhat.train.mean <- apply(yhat.train,2,mean)
		yhat.train.quantiles <- apply(yhat.train,2,quantile,c(.05,.5,.95))
	}
	if (nrow(x.test)){
                yhat.test = matrix(cres$tedraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
		yhat.test <- (rgy[2]-rgy[1])*(yhat.test+.5) + rgy[1]
		yhat.test.mean <- apply(yhat.test,2,mean)
		yhat.test.quantiles <- apply(yhat.test,2,quantile,c(.05,.5,.95))
	}

        varcount = matrix(cres$vcdraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
	
	return(invisible(list(
		call=match.call(),
		first.sigma=first.sigma,
		sigma=sigma,
		sigest=sigest,
		yhat.train=yhat.train,
		yhat.train.mean=yhat.train.mean,
		yhat.train.quantiles=yhat.train.quantiles,
		yhat.test=yhat.test,
		yhat.test.mean=yhat.test.mean,
		yhat.test.quantiles=yhat.test.quantiles,
                varcount=varcount
	)))

}
