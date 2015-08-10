#' LC50 lethal concentrations in the presence of additional stressors.
#'
#' The common 'critical effect concentration' used in toxicology bioassays is the Lethal concentration of 50% of the testing populaiton (LC50).
#' This package provides the facilities for estimating LC50s for a
#' toxin from sampled survival data at a single time-point in the presence of
#' additional stressors and non-ignorable control mortality. The result of which is the simultaneous calculations of multiple LC50s within a single experiemental framework where toxicants are tested within a framework of a factorial experiment.
#' Examples of facotrial arragement of multiple stressors include temperature, salinity, 
#'
#'
#' @name LC50-package
#' @docType package
#' @author S. Wotherspoon, A. Proctor.
NULL


##' Simulated toxicity data
##'
##' A simulated dataset showing individual survival following exposure
##' to a known toxin in the presence of the additional stressors
##' temperature and salinity.
##'
##' Samples of an aqueous solution of a toxin, of varying
##' concentrations, are prepared for each of three salinities and are
##' held at one of three constant temperatures.  Approximately 10
##' individuals were added to each sample and the number of survivors
##' recorded at the end of four days exposure.
##'
##' @format A data frame with 225 rows and 8 variables:
##' \describe{
##'   \item{vial}{the replicate}
##'   \item{temperature}{temperature of solution}
##'   \item{salinity}{salinity of solution}
##'   \item{group}{additinal stressor treatment groups}
##'   \item{conc}{concentration of toxin in solution}
##'   \item{total}{number of individuals in vial}
##'   \item{alive}{number of survivors after 4 days}
##'   \item{dead}{number of dead individuals after 4 days}
##' }
##'
"toxicity"


##' Estimate LC50 from survival data in the presence of additional
##' stressors and non-ignorable control mortality
##'
##' DESCRIBE the model.
##'
##' \code{lc50.fit} is the workhorse function: it is not normally
##' called directly but can be more efficient where the response
##' vector, design matrix and family have already been calculated.
##'
##'
##' @title Estimate LC50 for a toxin
##' @param formula a formula relating log LC50 to covariates
##' describing the aditional stressors.
##' @param concentration the name of variable that is the
##' concentration of the toxin.
##' @param group a factor distinguishing treatment groups for the additional stressors.
##' @param data data frame containing variables in the model.
##' @param start Starting values used to initialize the model.  If
##' \code{start=NULL} these parameters are determined by
##' \code{\link{lc50.initialize}}.
##' @param link the link function for survival fractions
##' @param common.background should a common background survival be
##' estimated for each treatment group.
##' @param X a design matrix
##' @param Y a two column matrix of responses
##' @param conc a vector of toxin concentrations
##' @param alpha vector of starting rate parameters
##' @param beta vector of starting coefficients
##' @param gamma vector of background survival parameters
##'
##' @return \code{lc50} returns an object of class inheriting from
##' "lc50". See later in this section.
##'
##' The function \code{\link{summary}} (i.e.,
##' \code{link{summary.lc50}}) can be used to obtain or print a
##' summary of the results and the function \code{\link{anova}} (i.e.,
##' \code{\link{anova.lc50}}) to produce an analysis of deviance table
##' for the tests of additional stressor effects.
##'
##' An LC50 model has several sets of coefficients, the generic
##' accessor function \code{\link{coef}} returns only the beta
##' coeffients.
##'
##' An object of class "lc50" is a list containing at least the
##' following components:
##'
##' \item{\code{logLik}}{the maximized log likelihood.}
##' \item{\code{aic}}{Akaike's information criteria.}
##' \item{\code{alpha}}{a vector of rate coefficients.}
##' \item{\code{alpha.cov}}{covariance of the rate coefficients.}
##' \item{\code{beta}}{a named vector of lc50 model coefficients.}
##' \item{\code{beta.cov}}{covariance of the lc50 model coefficients.}
##' \item{\code{gamma}}{a vector of background survival coefficients.}
##' \item{\code{gamma.cov}}{covariance of the background survival coefficients.}
##' \item{\code{coefficients}}{a named vector of lc50 model coefficients.}
##' \item{\code{cov.scaled}}{covariance of the lc50 model coefficients.}
##' \item{\code{loglc50}}{a named vector of log lc50s for the treatment groups.}
##' \item{\code{loglc50.cov}}{covariance of the lc50s for the treatment groups.}
##' \item{\code{concentration}}{a vector of taxin concentrations.}
##' \item{\code{group}}{a factor distinguishing treatment groups.}
##' \item{\code{x}}{a design matrix relating log lc50 to factors describing the additional stressors.}
##' \item{\code{y}}{a two column matrix of responses, giving the survivals and mortalities.}
##' \item{\code{fitted.values}}{the fitted probability of survival.}
##' \item{\code{deviance}}{the deviance.}
##' \item{\code{df.residual}}{the residual degrees of freedom.}
##' \item{\code{null.deviance}}{the deviance of the null model, which fits a single mortality rate to all data.}
##' \item{\code{df.null}}{the degrees of freedom for the null model.}
##' \item{\code{optim}}{the result of the call to \code{optim}.}
##' \item{\code{xlevels}}{a record of the levels of the factors used in fitting.}
##' \item{\code{contrasts}}{the contrasts used.}
##' \item{\code{call}}{the matched call.}
##' \item{\code{link}}{the link function.}
##' \item{\code{terms}}{the terms object used.}
##' \item{\code{model}}{the model frame.}
##'
##' @export
lc50 <- function(formula,concentration,group,data,start=NULL,link=c("probit","logit"),common.background=FALSE) {

  ## Record call and link function
  cl <- match.call()
  link <- match.arg(link)

  ## Create the model frame and terms
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "concentration", "group", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf,"terms")

  ## Create the model matrix and response
  X <- model.matrix(mt,mf)
  Y <- model.response(mf)

  ## Extract concentrations
  conc <- mf[,"(concentration)"]

  ## Determine the treatment groups
  group <- as.factor(mf[,"(group)"])

  ## Fit separate models to each group to generate initial parameter estimates
  if(is.null(start)) start <- lc50.initialize(Y,conc,group)
  alpha <- start$alpha
  gamma <- start$gamma
  beta <- qr.solve(X,start$loglc50[group])
  r <- lc50.fit(X,Y,conc,group,alpha,beta,gamma,link,common.background)
  r <- c(r,
         list(
           xlevels=.getXlevels(mt, mf),
           contrasts=attr(X,"contrasts"),
           call=cl,
           link=link,
           common.background=common.background,
           terms=mt,
           model=mf))
  class(r) <- "lc50"
  r
}


##' @rdname lc50
##' @importFrom MASS ginv
lc50.fit <- function(X,Y,conc,group,alpha,beta,gamma,link,common.background) {

  ## Decompose response
  y <- Y[,1]
  N <- rowSums(Y)

  ng <- nlevels(group)
  alpha.k <- seq_len(ng)
  beta.k <- ng+seq_len(ncol(X))
  gamma.k <- ng+ncol(X)+if(common.background) 1 else seq_len(ng)

  ## Index of first row of X for each group
  k <- match(levels(group),group)
  Xg <- X[k,,drop=FALSE]

  ## Select inverse link function
  ilink <- switch(link,probit=pnorm,logit=plogis)

  fitted.pq <- function(alpha,beta,gamma) {
    p <- ilink(alpha[group]*(log(conc)-X%*%beta))
    q <- ilink(if(!common.background) gamma[group] else gamma)
    ifelse(conc>0,p*q,q)
  }

  ## Negative log likelihood
  nlogL <- function(pars) {
    alpha <- pars[alpha.k]
    beta <- pars[beta.k]
    gamma <- rep(pars[gamma.k],length.out=ng)
    pq <- fitted.pq(alpha,beta,gamma)
    nll <- -sum(dbinom(y,N,pq,log=TRUE))
    if(!is.finite(nll)) nll <- .Machine$double.xmax
    nll
  }
  ## Minimize negative log likelihood
  mn <- optim(c(alpha,beta,gamma),nlogL,method="BFGS",hessian=TRUE)

  ## Basic parameters
  alpha <- mn$par[alpha.k]
  names(alpha) <- levels(group)
  beta <- mn$par[beta.k]
  names(beta) <- colnames(X)
  gamma <- mn$par[gamma.k]
  names(gamma) <- if(!common.background) levels(group) else "Common"

  ## Covariance of the beta (is subset of inverse hessian)
  V <- ginv(mn$hessian)
  alpha.cov <- V[alpha.k,alpha.k,drop=FALSE]
  colnames(alpha.cov) <- rownames(alpha.cov) <- levels(group)
  beta.cov <- V[beta.k,beta.k,drop=FALSE]
  colnames(beta.cov) <- rownames(beta.cov) <- colnames(X)
  gamma.cov <- V[gamma.k,gamma.k,drop=FALSE]
  colnames(gamma.cov) <- rownames(gamma.cov) <- if(!common.background) levels(group) else "Common"

  ## Compute lc50 and covariance by group
  loglc50 <- as.numeric(Xg%*%beta)
  names(loglc50) <- levels(group)
  loglc50.cov <- Xg%*%beta.cov%*%t(Xg)
  colnames(loglc50.cov) <- rownames(loglc50.cov) <- levels(group)

  ## Compute the deviance
  fitted <- fitted.pq(alpha,beta,gamma)
  deviance <- -2*sum(dbinom(y,N,fitted,log=T)-dbinom(y,N,y/N,log=T))
  df.residual <- nrow(X)-(ncol(X)+ng)
  null.deviance <- -2*sum(dbinom(y,N,sum(y)/sum(N),log=T)-dbinom(y,N,y/N,log=T))
  df.null <- nrow(X)-(1+ng)
  aic <- 2*(length(mn$par)+mn$value)

  r <- list(logLik=-mn$value,
            aic=aic,
            alpha=alpha,
            alpha.cov=alpha.cov,
            beta=beta,
            beta.cov=beta.cov,
            gamma=gamma,
            gamma.cov=gamma.cov,
            coefficients=beta,
            cov.scaled=beta.cov,
            loglc50=loglc50,
            loglc50.cov=loglc50.cov,
            concentration=conc,
            group=group,
            x=X,
            y=Y,
            fitted.values=fitted,
            deviance=deviance,
            df.residual=df.residual,
            null.deviance=null.deviance,
            df.null=df.null,
            optim=mn)
  class(r) <- "lc50"
  r
}


##' Estimate starting parameters for and LC50 model fit
##'
##' This is the default method for computing the starting values used
##' to initialize an \code{\link{lc50}} model.
##'
##' @title Starting parameters for an LC50 model fit
##' @param Y a two column matrix of the number of survivals and mortalities in each sample.
##' @param conc a vector of tixin concentrations
##' @param group a factor delineating treatment groups
##' @param link  the link function for survival fractions
##' @return Return a list of with components
##' \item{\code{alpha}}{the rate parameter for each treatment group}
##' \item{\code{gamma}}{the probit of the control surival for each treatment group}
##' \item{\code{loglc50}}{the log lc50 for each treatment group}
##' @export
lc50.initialize <- function(Y,conc,group,link=c("probit","logit")) {
  link <- match.arg(link)

  init <- function(Y,conc) {
    X <- cbind(1,ifelse(conc>0,1,0),ifelse(conc>0,log(conc),0))
    glm.fit(X,Y,family=binomial(link=link))$coefficient
  }

  cfs <- lapply(levels(group),function(g) init(Y[group==g,],conc[group==g]))
  list(alpha=sapply(cfs,function(cs) cs[3]),
       gamma=sapply(cfs,function(cs) cs[1]),
       loglc50=sapply(cfs,function(cs) -sum(cs[1:2])/cs[3]))
}


##' @export
print.lc50 <- function(x,digits = max(3L, getOption("digits") - 3L),...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits=digits),print.gap=2L,quote=FALSE)
  cat("log LC50:\n")
  print.default(format(x$loglc50, digits=digits),print.gap=2L,quote=FALSE)
  cat("Probit Control Survival:\n")
  print.default(format(x$gamma, digits=digits),print.gap=2L,quote=FALSE)
  invisible(x)
}


##' Summary method for class "\code{lc50}".
##'
##'
##' @title Summmarizing LC50 model fits
##' @param object an object of class \code{lc50}, obtained as the
##' result of a call to \code{\link{lc50}}
##' @param x an object of class \code{summary.lc50}, usually, a result
##' of a call to \code{summary.lc50}.
##' @param digits the number of significant digits to use when printing.
##' @param signif.stars logical. If \code{TRUE}, 'significance stars'
##' are printed for each coefficient.
##' @param ... additional parameters are ignored.
##' @return Returns an object of class \code{summary.lc50}, with components
##' \item{\code{coefficients}}{a table of coefficients.}
##' \item{\code{lc50}}{a table of LC50 for each treatment group.}
##' \item{\code{csurv}}{a table of control survival for each treatment group.}
##' @export
summary.lc50 <- function(object,...) {

  keep <- match(c("call","deviance","aic","contrasts","df.residual","null.deviance","df.null"),names(object),0L)
  cf <- object$coefficients
  cf.se <- sqrt(diag(object$cov.scaled))
  zvalue <- abs(cf)/cf.se
  pvalue <- pvalue <- 2*pnorm(-abs(zvalue))
  coef.table <- cbind(cf, cf.se, zvalue, pvalue)
  dimnames(coef.table) <- list(names(cf), c("Estimate","Std. Error","z value","Pr(>|z|)"))

  loglc50 <- object$loglc50
  loglc50.se <- sqrt(diag(object$loglc50.cov))
  lc50.table <- cbind(loglc50, loglc50.se, exp(loglc50), exp(loglc50-1.96*loglc50.se), exp(loglc50+1.96*loglc50.se))
  dimnames(lc50.table) <- list(names(loglc50), c("Estimate","Std. Error", "LC50", "Lwr 95%", "Upr 95%"))

  ilink <- switch(object$link,probit=pnorm,logit=plogis)
  gamma <- object$gamma
  gamma.se <- sqrt(diag(object$gamma.cov))
  csurv.table <- cbind(gamma,gamma.se,ilink(gamma),ilink(gamma-1.96*gamma.se),ilink(gamma-1.96*gamma.se))
  dimnames(csurv.table) <- list(names(gamma), c("Estimate","Std. Error", "C Surv", "Lwr 95%", "Upr 95%"))

  r <- c(list(coefficients=coef.table,
              lc50=lc50.table,
              csurv=csurv.table),
         object[keep])
  class(r) <- c("summary.lc50")
  r
}


##' @rdname summary.lc50
##' @export
print.summary.lc50 <- function(x,digits=max(3L,getOption("digits")-3L),
                               signif.stars=getOption("show.signif.stars"),...) {

  cat("\nCall:\n",paste(deparse(x$call),sep="\n",collapse="\n"),"\n\n",sep="")
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients,digits=digits,signif.stars=signif.stars,na.print="NA",...)
  cat("\n",
      apply(cbind(paste(format(c("Null", "Residual"), justify = "right"), "deviance:"),
                  format(unlist(x[c("null.deviance", "deviance")]), digits = max(5L, digits + 1L)),
                  " on",
                  format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
            1L, paste, collapse = " "), sep = "")
  cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)),"\n")
  cat("\nLC50:\n")
  printCoefmat(x$lc50,digits=digits,cs.ind=1:2,tst.ind=NULL,has.Pvalue=FALSE,na.print="NA",...)
  cat("\nControl Survival:\n")
  printCoefmat(x$csurv,digits=digits,cs.ind=1:2,tst.ind=NULL,has.Pvalue=FALSE,na.print="NA",...)
  cat("\n")
  invisible(x)
}






##' Compute an analysis of deviance table for an LC50 model fit.
##'
##' Specifying a single object gives a sequential analysis of deviance
##' table for that fit. That is, the reductions in the residual
##' deviance as each term of the formula is added in turn are given in
##' as the rows of a table, plus the residual deviances themselves.
##'
##' If more than one object is specified, the table has a row for the
##' residual degrees of freedom and deviance for each model. For all
##' but the first model, the change in degrees of freedom and deviance
##' is also given. (This only makes statistical sense if the models
##' are nested.) It is conventional to list the models from smallest
##' to largest, but this is up to the user.
##'
##' When \code{test} is "LRT" or "Chissq" the table will contain test
##' statistics (and P values) comparing the reduction in deviance for
##' the row to the residuals.
##'
##' @title Analysis of Deviance for lc50 model fits
##' @param object an object of class \code{lc50}, usually obtained as the
##' results from a call to \code{\link{lc50}}
##' @param ... additional objects of class \code{lc50}.
##' @param test a character string, partially matching one of
##' "\code{LRT}", "\code{Chisq}", or "\code{Cp}". See
##' \code{link{stat.anova}}.
##' @return An object of class \code{anova} inheriting from class
##' \code{data.frame}.
##' @export
anova.lc50 <- function(object,...,test = NULL)  {
  ## Handle multiple fits
  dotargs <- list(...)
  named <- if (is.null(names(dotargs))) rep_len(FALSE, length(dotargs)) else (names(dotargs) != "")
  if (any(named))
    warning("the following arguments to 'anova.lc50' are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
  dotargs <- dotargs[!named]
  is.lc50 <- vapply(dotargs, function(x) inherits(x, "lc50"), NA)
  dotargs <- dotargs[is.lc50]
  if (length(dotargs))
    return(anova.lc50list(c(list(object),dotargs),test = test))

  ## Passed single fit object
  varlist <- attr(object$terms, "variables")
  x <- model.matrix(object)
  varseq <- attr(x, "assign")
  nvars <- max(0, varseq)
  resdev <- resdf <- NULL
  if(nvars > 0) {
    for (i in seq_len(nvars)) {
      fit <- eval(call("lc50.fit",X=x[,varseq<i,drop=FALSE],
                       Y=object$y,conc=object$concentration,group=object$group,
                       alpha=object$alpha,beta=object$beta[varseq<i],gamma=object$gamma,
                       link=object$link,common.background=object$common.background))
      resdev <- c(resdev, fit$deviance)
      resdf <- c(resdf, fit$df.residual)
    }
  }
  resdf <- c(resdf, object$df.residual)
  resdev <- c(resdev, object$deviance)
  table <- data.frame(c(NA, -diff(resdf)), c(NA, pmax(0, -diff(resdev))), resdf, resdev)
  tl <- attr(object$terms, "term.labels")
  if (length(tl) == 0L) table <- table[1, , drop = FALSE]
  dimnames(table) <- list(c("NULL", tl), c("Df", "Deviance", "Resid. Df", "Resid. Dev"))
  title <- paste0("Analysis of Deviance Table",
                  "\n\nResponse: ", as.character(varlist[-1L])[1L],
                  "\n\nTerms added sequentially (first to last)\n\n")
  df.dispersion <- object$df.residual
  if (!is.null(test))
    table <- stat.anova(table=table,test=test,scale=1,df.scale=df.dispersion,n=NROW(x))
  structure(table, heading = title, class = c("anova", "data.frame"))
}

##' @rdname anova.lc50
##' @export
anova.lc50list <- function (object, ..., test = NULL) {
  responses <- as.character(lapply(object, function(x) deparse(formula(x)[[2L]])))
  sameresp <- responses == responses[1L]
  if (!all(sameresp)) {
    object <- object[sameresp]
    warning(gettextf("models with response %s removed because response differs from model 1",
                     sQuote(deparse(responses[!sameresp]))), domain = NA)
  }
  ns <- sapply(object, function(x) length(x$residuals))
  if (any(ns != ns[1L]))
    stop("models were not all fitted to the same size of dataset")
  nmodels <- length(object)
  if (nmodels == 1)
    return(anova.lc50(object[[1L]], test = test))
  resdf <- as.numeric(lapply(object, function(x) x$df.residual))
  resdev <- as.numeric(lapply(object, function(x) x$deviance))
  table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), c(NA, -diff(resdev)))
  variables <- lapply(object, function(x) paste(deparse(formula(x)), collapse = "\n"))
  dimnames(table) <- list(1L:nmodels, c("Resid. Df", "Resid. Dev", "Df", "Deviance"))
  title <- "Analysis of Deviance Table\n"
  topnote <- paste("Model ", format(1L:nmodels), ": ", variables, sep = "", collapse = "\n")
  if (!is.null(test)) {
    bigmodel <- object[[order(resdf)[1L]]]
    df.dispersion <- min(resdf)
    table <- stat.anova(table=table,test=test,scale=1,
                        df.scale=df.dispersion,n=length(bigmodel$residuals))
  }
  structure(table, heading = c(title, topnote), class = c("anova", "data.frame"))
}




## This is ripped off coef.glm
##' @export
coef.lc50 <- function(object,...) {
  object$coefficients
}

## This is ripped off from vcov.glm
##' @export
vcov.lc50 <- function(object,...) {
  object$cov.scaled
}

## This is ripped off from model.matrix.lm
##' @export
model.matrix.lc50 <- function(object,...) {
  object$x
}


## This is ripped off from model.frame.lm
##' @export
model.frame.lc50 <- function(formula, ...)  {
  dots <- list(...)
  nargs <- dots[match(c("data","concentration","group"),names(dots),0)]
  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    m <- match(c("formula","data","concentration","group"),names(fcall),0L)
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    fcall[[1L]] <- quote(stats::model.frame)
    fcall$xlev <- formula$xlevels
    fcall$formula <- terms(formula)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms)
    if (is.null(env))
      env <- parent.frame()
    eval(fcall, env)
  }
  else formula$model
}



## This is ripped off from simulate.lm
##' @importFrom stats simulate
##' @export
simulate.lc50 <- function(object, nsim=1, seed=NULL, ...) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  sim <- function(.) {
    y <- rbinom(n,N,ftd)
    y <- cbind(y,N-y)
    colnames(y) <- colnames(object$y)
  }

  ftd <- object$fitted.values
  n <- length(ftd)
  N <- rowSums(object$y)
  val <- vector("list", nsim)
  for (i in seq_len(nsim)) {
    y <- rbinom(n,N,ftd)
    y <- cbind(y,N-y)
    colnames(y) <- colnames(object$y)
    val[[i]] <- y
  }
  class(val) <- "data.frame"
  names(val) <- paste("sim", seq_len(nsim), sep = "_")
  row.names(val) <- rownames(object$y)
  attr(val, "seed") <- RNGstate
  val
}


##' Predicted survival from a fitted LC50 object
##'
##' If \code{newdata} is omitted the predictions are based on the data
##' used for the fit.  For \code{type="adjusted"}, it is assumed there
##' is no background mortality and predictions are made based purely
##' on the mortality due to the toxin.
##'
##' @title Predict method for LC50 fits
##' @param object object an object of class \code{lc50}, obtained as the
##' result of a call to \code{\link{lc50}}
##' @param newdata optionally, a data frame in which to look for
##' variables with which to predict. If omitted, the fitted linear
##' predictors are used.
##' @param type the type of prediction required. If
##' \code{type="response"} predictions include both the background
##' mortality and the mortality due to the toxin, but if
##' \code{type="adjusted"} predictions only reflect mortality due to
##' the toxin.
##' @param ... further arguments passed to or from other methods.
##' @return a vector of predicted survival fractions.
##' @export
predict.lc50 <- function (object, newdata, type = c("response", "adjusted"),...) {

  type <- match.arg(type)

  ## Get model frame and design matrix
  if (missing(newdata) || is.null(newdata)) {
    mf <- object$model
    X <- object$x
  } else {
    mt <- delete.response(object$terms)
    mf <- substitute(model.frame(mt,newdata,xlev=object$xlevels,concentration=concentration,group=group),
                     as.list(object$call[c("concentration","group")]))
    mf <- eval(mf)
    X <- model.matrix(mt,mf,contrasts.arg=object$contrasts)
  }

  ## Extract concentrations
  conc <- mf[,"(concentration)"]

  ## Determine the treatment groups
  group <- factor(mf[,"(group)"],levels(object$group))

  ## Extract coefficients
  alpha <- object$alpha
  beta <- object$beta
  gamma <- object$gamma

  ## Select inverse link function
  ilink <- switch(object$link,probit=pnorm,logit=plogis)

  ## Predict survival fraction
  p <- ilink(alpha[group]*(log(conc)-X%*%beta))
  q <- if(type=="adjusted") 1 else ilink(if(!object$common.background) gamma[group] else gamma)
  ifelse(conc>0,p*q,q)
}

