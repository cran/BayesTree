\name{bart}
\title{Bayesian Additive Regression Trees}
\alias{bart}
\alias{plot.bart}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model:
\eqn{y = f(x) + \epsilon}{y = f(x) + e},
where \eqn{f} is the sum of many tree models
and \eqn{\epsilon \sim N(0,\sigma^2)}{e ~ N(0,sigma\^2)}.

The goal is to have very flexible inference for the uknown
function \eqn{f}.

In the spirit of \dQuote{ensemble models}, 
each tree is constrained by a prior to be a weak learner
so that it contributes a
small amount to the overall fit.
}
\usage{
bart(
   x.train, y.train, x.test=matrix(0.0,0,0),
   sigest=NA, sigdf=3, sigquant=.90,
   k=2.0,
   power=2.0, base=.95,
   ntree=200,
   ndpost=1000, nskip=100,
   printevery=100, keepevery=1, keeptrainfits=TRUE,
   usequants=FALSE, numcut=100, printcutoffs=0,
   verbose=TRUE)
\method{plot}{bart}(
   x,
   plquants=c(.05,.95), cols =c('blue','black'),
   \dots)
}
\arguments{
   \item{x.train}{
   Explanatory variables for training (in sample) data.\cr
   Must be a matrix (typeof double)
   with (as usual) rows corresponding to observations and columns to variables.\cr
   Note that for a categorical variable you need to use dummies and if there
   are more than two categories, you need to put all the dummies in
   (unlike linear regression).\cr
   \code{bart} will generate draws of \eqn{f(x)} for each \eqn{x}
   which is a row of x.train.
   }

   \item{y.train}{
   Dependent variable for training (in sample) data.\cr
   Must be a vector (typeof double) with length
   equal to the number of observations
   (equal to the number of rows of x.train).
   }

   \item{x.test}{
   Explanatory variables for test (out of sample) data.\cr
   Must be a (typeof double) matrix with the same
   number of columns as x.train.\cr
   \code{bart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of x.test.
   }

   \item{sigest}{
   The prior for the error variance  (\eqn{\sigma^2}{sigma\^2}) is inverted chi-squared
   (the standard conditionally conjugate prior).
   The prior is specified by choosing the degrees of freedom, a rough estimate of the
   corresponding standard deviation and a quantile to put this rough estimate at.
   If sigest=NA then the rough estimate will be the usual least squares estimator.
   Otherwise the supplied value will be used.
   }

   \item{sigdf}{
   Degrees of freedom for error variance prior.
   }

   \item{sigquant}{
   The quantile of the prior that the rough estimate (see sigest) is placed at.
   The closer the quantile is to 1,
   the more aggresive the fit will be as you are putting more prior weight
   on error standard deviations (\eqn{\sigma}{sigma}) less than the rough estimate.
   }

   \item{k}{
   The number of prior standard deviations \eqn{E(Y|x)} is away from +/-.5.
   The bigger this value is, the more conservative the fitting will be.
   The response (y.train) is internally scaled to range from -.5 to .5.
   }

   \item{power}{
   Power parameter for tree prior.
   }

   \item{base}{
   Base parameter for tree prior.
   }

   \item{ntree}{
   The number of trees in the sum.
   }

   \item{ndpost}{
   The number of posterior draws after burn in, ndpost/keepevery will actually be returned.
   }

   \item{nskip}{
   Number of MCMC iterations to be treated as burn in.
   }

   \item{printevery}{
   As the MCMC runs, a message is printed every printevery draws.
   }

   \item{keepevery}{
   Every keepevery draw is kept to be returned to the user.\cr
   A \dQuote{draw} will consist of values of the error standard deviation (\eqn{\sigma}{sigma})
   and \eqn{f^*(x)}{f*(x)}
   at \eqn{x} = rows from the train(optionally) and test data, where \eqn{f^*}{f*} denotes
   the current draw of \eqn{f}.
   }

   \item{keeptrainfits}{
   If true the draws of \eqn{f(x)} for \eqn{x} = rows of x.train are returned.
   }

   \item{usequants}{
   Decision rules in the tree are of the form
   \eqn{x \le c}{x <= c}  vs. \eqn{x > c}
   for each variable corresponding to a column of x.train.
   usequants determines how the set of possible c is determined.
   If usequants is true, then the c are
   a subset of the values (xs[i]+xs[i+1])/2 where xs is
   unique sorted values obtained from
   the corresponding column of x.train.
   If usequants is false, the cutoffs are equally spaced across the
   range of values taken on by the corresponding
   column of x.train.
   }
   
   \item{numcut}{
   The number of possible values of c (see usequants).
   If a single number if given, this is used for all variables.
   Otherwise a vector with length equal to ncol(x.train) is required,
   where the \eqn{i^{th}}{i^th} element gives the number of c used for
   the \eqn{i^{th}}{i^th} variable in x.train.
   If usequants is false, numcut equally spaced cutoffs
   are used covering the range of values in the corresponding
   column of x.train.  If usequants is true, then  min(numcut, the number of unique values in the
   corresponding columns of x.train - 1) c values are used.
   }

   \item{printcutoffs}{
   The number of cutoff rules c to printed to screen before the MCMC is run.
   Give a single integer, the same value will be used for all variables.
   If 0, nothing is printed.
   }

   \item{verbose}{
   Logical, if FALSE supress printing.
   }

   \item{x}{
   Value returned by \code{bart} which contains the information to be plotted.
   }

   \item{plquants}{
   In the plots, beliefs about \eqn{f(x)} are indicated by plotting the
   posterior median and a lower and upper quantile.
   plquants is a double vector of length two giving the lower and upper quantiles.
   }

   \item{cols}{
   Vector of two colors.  First color is used to plot the median of \eqn{f(x)} and
   the second color is used to plot the lower and upper quantiles.
   }

   \item{\dots}{
   Additional arguments passed on to plot.
   }
}
\details{
   BART is an Bayesian MCMC method for the model \eqn{y = f(x) + \epsilon}{y = f(x) + e},
   \eqn{\epsilon \sim N(0,\sigma^2)}{e ~ N(0,sigma\^2)}.\cr
   At each MCMC interation, we produce a draw from the joint posterior
   \eqn{(f,\sigma) | (x,y)}{(f,sigma) \| (x,y)}.

   Thus, unlike a lot of other modelling methods in R, we do not produce a single model object
   from which fits and summaries may be extracted.  The output consists of values
   \eqn{f^*(x)}{f*(x)} and \eqn{\sigma^*}{sigma*} where * denotes a particular draw.
   The \eqn{x} is either a row from the training data (x.train) or the test data (x.test).
}
\value{
   The \code{plot} method sets mfrow to c(1,2) and makes two plots.
   
   The first plot is the sequence of kept draws of \eqn{\sigma}{sigma}
   including the burn-in draws.  Initially these draws will decline as BART finds fit
   and then level off when the MCMC has burnt in.\cr
   The second plot has \eqn{y} on the horizontal axis and posterior intervals for
   the corresponding \eqn{f(x)} on the vertical axis. 
   

   \code{bart} returns a list assigned class \sQuote{bart}.
   The list has components:

   \item{yhat.train}{
   A matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }
   
   \item{yhat.test}{Same as yhat.train but now the x's are the rows of the test data.}
   \item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
   \item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
   \item{sigma}{post burn in draws of sigma, length = ndpost/keepevery.}
   \item{first.sigma}{burn-in draws of sigma.}
   \item{varcount}{a matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}

   \item{sigest}{
   The rough error standard deviation (\eqn{\sigma}{sigma}) used in the prior.
   }
 
   \item{y}{
   The input dependent vector of values for the dependent variable.\cr
   This is used in plot.bart.
   }
}
\references{
Chipman, H., George, E., and McCulloch, R. (2006)
   BART: Bayesian Additive Regression Trees.

Chipman, H., George, E., and McCulloch R. (2006)
   Bayesian Ensemble Learning.

both of the above at:
\url{http://faculty.chicagogsb.edu/robert.mcculloch/research/rob-mcculloch-cv.html}

Friedman, J.H. (1991)
   Multivariate adaptive regression splines.
        \emph{The Annals of Statistics}, \bold{19}, 1--67.
}
\note{
\emph{There was a bug in BayesTree\_0.1-0} (now fixed of course).\cr
If the number of test observations was less than the number of trees
(200 is the default), the yhat.test and yhat.test.mean components were suspect.
}
\author{
Hugh Chipman: \email{hugh.chipman@acadiau.ca}\cr
Robert McCulloch: \email{robert.mcculloch@chicagogsb.edu}.
}
\seealso{
\code{\link{pdbart}}
}
\examples{
##simulate data (example from Friedman MARS paper)
f = function(x){
10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0  #y = f(x) + sigma*z , z~N(0,1)
n = 100      #number of observations
set.seed(99)
x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(x)
y=Ey+sigma*rnorm(n)
lmFit = lm(y~.,data.frame(x,y)) #compare lm fit to BART later
##run BART
set.seed(99)
bartFit = bart(x,y)
plot(bartFit) # plot bart fit
##compare BART fit to linear matter and truth = Ey
fitmat = cbind(y,Ey,lmFit$fitted,bartFit$yhat.train.mean)
colnames(fitmat) = c('y','Ey','lm','bart')
print(cor(fitmat))
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
