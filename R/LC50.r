#' LC50 lethal concentrations in the presence of additional stressors.
#'
#' Provides facilities for estimating lethal concentrations for a
#' toxin from destructively sampled survival data in the presence of
#' additional stressors and non-ignorable control mortality.
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
##' @param quasi should a quasibinomial model be fitted.
##' @param common.background should a common background survival be
##' estimated for each treatment group.
##' @param optim.control control parameters for \code{optim}
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
##' \item{\code{cov.unscaled}}{covariance of the lc50 model coefficients.}
##' \item{\code{loglc50}}{a named vector of log lc50s for the treatment groups.}
##' \item{\code{loglc50.cov}}{covariance of the lc50s for the treatment groups.}
##' \item{\code{concentration}}{a vector of taxin concentrations.}
##' \item{\code{group}}{a factor distinguishing treatment groups.}
##' \item{\code{x}}{a design matrix relating log lc50 to factors describing the additional stressors.}
##' \item{\code{y}}{a two column matrix of responses, giving the survivals and mortalities.}
##' \item{\code{fitted.values}}{the fitted probability of survival.}
##' \item{\code{residuals}}{the deviance residuals for the fit.}
##' \item{\code{deviance}}{the deviance for the fit.}
##' \item{\code{df.residual}}{the residual degrees of freedom.}
##' \item{\code{dispersion}}{the dispersion.}
##' \item{\code{null.deviance}}{the deviance of the null model, which fits a single mortality rate to all data.}
##' \item{\code{df.null}}{the degrees of freedom for the null model.}
##' \item{\code{optim}}{the result of the call to \code{optim}.}
##' \item{\code{xlevels}}{a record of the levels of the factors used in fitting.}
##' \item{\code{contrasts}}{the contrasts used.}
##' \item{\code{call}}{the matched call.}
##' \item{\code{link}}{the link function.}
##' \item{\code{quasi}}{is the dispersion estimated.}
##' \item{\code{terms}}{the terms object used.}
##' \item{\code{model}}{the model frame.}
##'
##' @export
lc50 <- function(formula,concentration,group,data,start=NULL,
                 link=c("probit","logit"),quasi=FALSE,
                 common.background=FALSE,optim.control=list()) {

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
  gamma <- if(common.background) mean(start$gamma) else start$gamma
  beta <- qr.solve(X,start$loglc50[group])
  r <- lc50.fit(X,Y,conc,group,alpha,beta,gamma,link,quasi,common.background,optim.control)
  r <- c(r,
         list(
           xlevels=.getXlevels(mt, mf),
           contrasts=attr(X,"contrasts"),
           call=cl,
           link=link,
           quasi=quasi,
           common.background=common.background,
           terms=mt,
           model=mf))
  class(r) <- "lc50"
  r
}


##' @rdname lc50
##' @importFrom MASS ginv
lc50.fit <- function(X,Y,conc,group,alpha,beta,gamma,link,
                     quasi=FALSE,common.background=FALSE,optim.control=list()) {

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
  mn <- optim(c(alpha,beta,gamma),nlogL,method="BFGS",hessian=TRUE,control=optim.control)

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
  residuals <- sign(y/N-fitted)*sqrt(abs(2*(dbinom(y,N,fitted,log=T)-dbinom(y,N,y/N,log=T))))
  deviance <- -2*sum(dbinom(y,N,fitted,log=T)-dbinom(y,N,y/N,log=T))
  df.residual <- nrow(X)-(length(alpha.k)+length(beta.k)+length(gamma.k))
  null.deviance <- -2*sum(dbinom(y,N,sum(y)/sum(N),log=T)-dbinom(y,N,y/N,log=T))
  df.null <- nrow(X)-(length(alpha.k)+1+length(gamma.k))
  dispersion <- if(quasi) sum((y-N*fitted)^2/(N*fitted*(1-fitted)))/df.residual else 1
  aic <- if(quasi) NA else 2*(length(mn$par)+mn$value)

  r <- list(logLik=-mn$value,
            aic=aic,
            alpha=alpha,
            alpha.cov=alpha.cov,
            beta=beta,
            beta.cov=beta.cov,
            gamma=gamma,
            gamma.cov=gamma.cov,
            coefficients=beta,
            cov.unscaled=beta.cov,
            loglc50=loglc50,
            loglc50.cov=loglc50.cov,
            concentration=conc,
            group=group,
            x=X,
            y=Y,
            fitted.values=fitted,
            residuals=residuals,
            deviance=deviance,
            dispersion=dispersion,
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
  cat("Baseline Survival:\n")
  print.default(format(x$gamma, digits=digits),print.gap=2L,quote=FALSE)
  cat("Rate:\n")
  print.default(format(x$alpha, digits=digits),print.gap=2L,quote=FALSE)
  invisible(x)
}


##' Summary method for class "\code{lc50}".
##'
##'
##' @title Summmarizing LC50 model fits
##' @param object an object of class \code{lc50}, obtained as the
##' result of a call to \code{\link{lc50}}
##' @param background if \code{TRUE} a summary table for the background survival is calculated
##' @param rate if \code{TRUE} a summary table for the rate parameters is calculated
##' @param ... additional parameters are ignored.
##' @param x an object of class \code{summary.lc50}, usually, a result
##' of a call to \code{summary.lc50}.
##' @param digits the number of significant digits to use when printing.
##' @param signif.stars logical. If \code{TRUE}, 'significance stars'
##' are printed for each coefficient.
##' @return Returns an object of class \code{summary.lc50}, with components
##' \item{\code{coefficients}}{a table of coefficients.}
##' \item{\code{lc50}}{a table of LC50 for each treatment group.}
##' \item{\code{bsurv}}{optionally, a table of background survival for each treatment group.}
##' \item{\code{rate}}{optionally, a table of rates for each treatment group.}
##' @export
summary.lc50 <- function(object,background=TRUE,rate=FALSE,...) {

  keep <- match(c("call","deviance","aic","dispersion","contrasts",
                  "df.residual","null.deviance","df.null"),names(object),0L)
  cf <- object$coefficients
  cf.se <- sqrt(diag(object$dispersion*object$cov.unscaled))
  zvalue <- abs(cf)/cf.se
  pvalue <- pvalue <- 2*pnorm(-abs(zvalue))
  coef.table <- cbind(cf, cf.se, zvalue, pvalue)
  dimnames(coef.table) <- list(names(cf), c("Estimate","Std. Error","z value","Pr(>|z|)"))

  loglc50 <- object$loglc50
  loglc50.se <- sqrt(diag(object$dispersion*object$loglc50.cov))
  lc50.table <- cbind(loglc50, loglc50.se, exp(loglc50), exp(loglc50-1.96*loglc50.se), exp(loglc50+1.96*loglc50.se))
  dimnames(lc50.table) <- list(names(loglc50), c("Estimate","Std. Error", "LC50", "Lwr 95%", "Upr 95%"))

  r <- c(list(coefficients=coef.table,
              lc50=lc50.table),
         object[keep])

  if(background) {
    ilink <- switch(object$link,probit=pnorm,logit=plogis)
    gamma <- object$gamma
    gamma.se <- sqrt(diag(object$dispersion*object$gamma.cov))
    bsurv.table <- cbind(gamma,gamma.se,ilink(gamma),ilink(gamma-1.96*gamma.se),ilink(gamma+1.96*gamma.se))
    dimnames(bsurv.table) <- list(names(gamma), c("Estimate","Std. Error", "Survival", "Lwr 95%", "Upr 95%"))
    r$bsurv <- bsurv.table
  }

  if(rate) {
    alpha <- object$alpha
    alpha.se <- sqrt(diag(object$dispersion*object$alpha.cov))
    rate.table <- cbind(alpha,alpha.se,alpha-1.96*alpha.se,alpha+1.96*alpha.se)
    dimnames(rate.table) <- list(names(alpha), c("Estimate","Std. Error", "Lwr 95%", "Upr 95%"))
    r$rate <- rate.table
  }

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
  cat("\n(Dispersion parameter taken to be ", format(x$dispersion), ")\n\n",
      apply(cbind(paste(format(c("Null", "Residual"), justify = "right"), "deviance:"),
                  format(unlist(x[c("null.deviance", "deviance")]), digits = max(5L, digits + 1L)),
                  " on",
                  format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
            1L, paste, collapse = " "), sep = "")
  cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)),"\n")
  cat("\nLC50:\n")
  printCoefmat(x$lc50,digits=digits,cs.ind=1:2,tst.ind=NULL,has.Pvalue=FALSE,na.print="NA",...)
  if(!is.null(x$bsurv)) {
    cat("\nBackground Survival:\n")
    printCoefmat(x$bsurv,digits=digits,cs.ind=1:2,tst.ind=NULL,has.Pvalue=FALSE,na.print="NA",...)
    cat("\n")
  }
  if(!is.null(x$rate)) {
    cat("\nRate:\n")
    printCoefmat(x$rate,digits=digits,cs.ind=1:2,tst.ind=NULL,has.Pvalue=FALSE,na.print="NA",...)
    cat("\n")
  }
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
##' When \code{test} is "LRT" or "Chisq" the table will contain test
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
  df.dispersion <- if(object$quasi) object$df.residual else Inf
  if (!is.null(test))
    table <- stat.anova(table=table,test=test,scale=object$dispersion,df.scale=df.dispersion,n=NROW(x))
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
    df.dispersion <- if(bigmodel$quasi) bigmodel$df.residual else Inf
    table <- stat.anova(table=table,test=test,scale=bigmodel$dispersion,
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
  object$dispersion*object$cov.unscaled
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
  if(object$quasi) stop("Cannot simulate from a quasi model")
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



##' Bayesian estimates of LC50 from survival data in the presence of additional
##' stressors and non-ignorable control mortality
##'
##' This function is an analog of \code{\link{lc50}} that produces an
##' object of class \code{jags} which can be used to draw samples from
##' the posterior using \code{update} and \code{coda.samples} from
##' \pkg{rjags}.
##'
##' The model assumes half Normal priors for \code{alpha} and Normal
##' priors for \code{beta} and \code{gamma}.  For \code{alpha} and
##' \code{gamma}, a single prior mean and precision is assumed for all
##' groups, for \code{beta} individual prior means and precisions can be specified
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
##' @param n.adapt parameter passed to \code{jags.model}
##' @param n.chains parameter passed to \code{jags.model}
##' @param alpha.mu prior mean for alpha
##' @param alpha.tau prior precision for alpha
##' @param beta.mu either a single prior mean for all beta parameters,
##' or a vector of prior means, one for each parameter.
##' @param beta.tau either a single prior precision for all beta
##' parameters, or a vector of prior precisions, one for each
##' parameter.
##' @param gamma.mu prior mean for gamma
##' @param gamma.tau prior precision for gamma
##' @return Returns an object inheriting from class \code{jags} which
##' can be used to generate dependent samples from the posterior
##' distribution of the parameters
##' @export
lc50JAGS <- function(formula,concentration,group,data,start=NULL,link=c("probit","logit"),
                     common.background=FALSE,n.adapt=500,n.chains=4,
                     alpha.mu=0,alpha.tau=0.0001,
                     beta.mu=0,beta.tau=0.0001,
                     gamma.mu=0,gamma.tau=0.0001) {

  if(!requireNamespace("rjags",quietly=TRUE)) {
    stop("lc50Jags requires the rjags package to be installed", call. = FALSE)
  }

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

  ## Check lengths of prior hyper parameters for beta
  if(!(length(beta.mu)==1 || length(beta.mu)==ncol(X)))
    stop("Vector of prior means for beta is wrong length")
  beta.mu <- rep(beta.mu,length.out=ncol(X))

  if(!(length(beta.tau)==1 || length(beta.tau)==ncol(X)))
    stop("Vector of prior precisions for beta is wrong length")
  beta.tau <- rep(beta.tau,length.out=ncol(X))

  ## Fit separate models to each group to generate initial parameter estimates
  if(is.null(start)) start <- lc50.initialize(Y,conc,group)
  start$beta <- qr.solve(X,start$loglc50[group])
  if(common.background) start$gamma <- mean(start$gamma)

  ## Index of first row of X for each group
  k <- match(levels(group),group)
  Xg <- X[k,,drop=FALSE]


  ## Define BUGS model
  if(!common.background) {

    bugs.model <- paste("
model {
  ## Likelihood
  for(i in 1:N) {
    alive[i] ~ dbin(pi[i],total[i])
    pi[i] <- q[group[i]]*ifelse(zero[i]>0,1,p[i])
    ",link,"(p[i]) <- alpha[group[i]]*(log(conc[i]+zero[i]) - loglc50[group[i]])
  }

  loglc50 <- X %*% beta
  for(i in 1:Ngroup) {
    log(lc50[i]) <- loglc50[i]
    ",link,"(q[i]) <- gamma[i]
  }

  ## Priors
  for(i in 1:Ngroup) {
    ## alpha's must be negative
    alpha[i] ~ dnorm(alpha.mu,alpha.tau)T(,0)
  }
  for(i in 1:Ncoef) {
    beta[i] ~ dnorm(beta.mu[i],beta.tau[i])
  }
  for(i in 1:Ngroup) {
    gamma[i] ~ dnorm(gamma.mu,gamma.tau)
  }
}",sep="")

  } else {

    bugs.model <- paste("
model {
  ## Likelihood
  for(i in 1:N) {
    alive[i] ~ dbin(pi[i],total[i])
    pi[i] <- q*ifelse(zero[i]>0,1,p[i])
    ",link,"(p[i]) <- alpha[group[i]]*(log(conc[i]+zero[i]) - loglc50[group[i]])
  }

  loglc50 <- X %*% beta
  for(i in 1:Ngroup) {
    log(lc50[i]) <- loglc50[i]
  }
  ",link,"(q) <- gamma

  ## Priors
  for(i in 1:Ngroup) {
    ## alpha's must be negative
    alpha[i] ~ dnorm(alpha.mu,alpha.tau)T(,0)
  }
  for(i in 1:Ncoef) {
    beta[i] ~ dnorm(beta.mu[i],beta.tau[i])
  }
  gamma ~ dnorm(gamma.mu,gamma.tau)
}",sep="")

  }

  model <- rjags::jags.model(textConnection(bugs.model),
                      data = list(
                        "alive" = Y[,1],
                        "total" = rowSums(Y),
                        "conc" = conc,
                        "group" = group,
                        "zero" = as.numeric(conc==0),
                        "X" = Xg,
                        "N" = nrow(Y),
                        "Ncoef" = ncol(X),
                        "Ngroup" = nlevels(group),
                        "alpha.mu" = alpha.mu,
                        "alpha.tau" = alpha.tau,
                        "beta.mu" = beta.mu,
                        "beta.tau" = beta.tau,
                        "gamma.mu" = gamma.mu,
                        "gamma.tau" = gamma.tau),
                      inits=start[c("alpha","beta","gamma")],
                      n.chains = n.chains,
                      n.adapt = n.adapt)

  model
}
