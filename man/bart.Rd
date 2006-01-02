\name{bart}
\alias{bart}
\title{Bayesian Additive Regression Trees}
\description{ BART is a Bayesian ``sum-of-trees'' model:
\eqn{y = f(x) + \epsilon}{y = f(x) + e},
where f is the sum of many tree models and \eqn{\epsilon~N(0,\sigma^2)}{e ~ N(0,sigma^2)}.
Each tree is constrained by a prior to be a weak leaner.
Fitting and inference are accomplished via an iterative
backfitting MCMC algorithm. This model is motivated by ensemble
methods in general, and boosting algorithms in particular. Like
boosting, each weak learner (i.e., each weak tree) contributes a
small amount to the overall model, and the training of a weak
learner is conditional on the estimates for the other weak
learners. The differences from the boosting algorithm are just as
striking as the similarities: BART is defined by a statistical
model: a prior and a likelihood, while boosting is defined by an
algorithm. MCMC is used both to fit the model and to qualify
predictive inference. 
}
\keyword{nonparametric}
\keyword{tree}
\usage{
bart(x.train,y.train,x.test=matrix(0.0,0,0),
        sigest=NA,sigdf=3, sigquant=.90, k=2.0, ntree=200,
        ndpost=1000,nskip=100,
        printevery=100,keepevery=1,keeptrainfits=TRUE,
        numcut=100)
}
\arguments{
\item{x.train}{explanatory variables for training (in sample) data - must be a matrix
with (as usual) rows corresponding to observations and columns to variables.
bart will generate draws of f(x) for each x which is a row of x.train.
Right now, factors are not supported so you have to use dummies.
Note that for bart, if there are more than two levels you put in all the dummies.}
\item{y.train}{dependent variable for training (in sample) data - must be a numeric vector with length
equal to the number of observations.}
\item{x.test}{explanatory variables for test (out of sample) data - must be a matrix with the same
number of columns as x.train, bart will generate draws of f(x) for each x which is a row of x.test.}
\item{sigest}{The prior for the variance of the error is inverted chi-squared (the standard conditionally
conjugate prior).  The prior is specified by choosing the degrees of freedom, a rough estimate of the 
corresponding standard deviation and a quantile to put this rough estimate at.
If sigest=NA then the rough estimate will be obtained from the usual least squares estimator.
Otherwise the supplied value will be used.}
\item{sigdf}{degrees of freedom for error variance prior.}
\item{sigquant}{the quantile of the prior that the rough estimate (see sigest) is placed at.
The closer the quantile is to 1, the more aggresive the fit will be as you are putting more prior weight
on error standard deviations (sigma) less than the rough estimate.}
\item{k}{the number of prior standard deviations E(Y|x) is away from +/-.5.  The bigger this value is,
the more conservative the fitting will be.}
\item{ntree}{the number of trees in the sum.}
\item{ndpost}{the number of posterior draws after burn in, ndpost/keepevery will actually be returned.}
\item{nskip}{number of mcmc iterations to be treated as burn in.}
\item{printevery}{as the mcmc runs, a message is printed every printevery draws.}
\item{keepevery}{every keepevery draw is kept to be returned to the user, a "draw" will consist
of values of the error standard deviation and f(x) at x = rows from the train(optionally) and test data.}
\item{keeptrainfits}{if true the draws of f(x) for x=rows of x.train are returned.}
\item{numcut}{the number of equally spaced values between the min and max of each explanatory variable
used as cut-points in the tree decision rules.}
}
\value{
a list containing components:
\item{yhat.train}{a matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
Each row is the current draws of f(x) for each x which is a row of x.train, burn-in is dropped.}
\item{yhat.test}{same as yhat.train but now the x's are the rows of the test data.}
\item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
\item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
\item{sigma}{post burn in draws of sigma, length = ndpost/keepevery.}
\item{first.sigma}{burn-in draws of sigma.}
\item{yhat.train.quantile}{matrix with 3 rows and nrow(x) columns, ith column gives
5\%, 50\%, and 95\% quantile for f(x) draws where x is ith row of train x.}
\item{yhat.test.quantile}{same as yhat.train.quantile except for test x.}
\item{varcount}{a matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
Each row is for a draw. For each variable (corresponding to the columns), the total count of the number of times
that variable is used in a tree decision rule (over all trees) is given.}
}
\author{
Hugh Chipman and Robert McCulloch.
See \url{http://gsbwww.uchicago.edu/fac/robert.mcculloch/research/code/BART/BART_Code.html}.
}
\examples{
data(cheese)
##weekly data.
##first three columns are dummies indicating which of three New York retailers the data is from.
##the fourth column is the log of price.
##the fifth column is measure of how often the item is advertised through an in-store display.
##the sixth column is the log of weekly sales volume, which we treat as the dependent variable.

##fit bart
##note that you use all the dummies in the x for bart!
x=as.matrix(cheese[,1:5]) 
y=cheese[,6] 
set.seed(99) 
bartFit = bart(x,y)

##fit linear model
##drop the first dummy for linear regression
lmFit = lm(y~.,cheese[,2:6]) 

##compare fits
cat("Squared correlation between y and fits (R^2) from linear model:", cor(y,lmFit$fitted)^2,"\n",
"Squared correlation from bart model:", cor(y,bartFit$yhat.train.mean)^2,"\n")

##I got .75 for the linear model and .89 for bart, so bart has better in-sample fit.
##Of course, out-of-sample is always another matter.
}
