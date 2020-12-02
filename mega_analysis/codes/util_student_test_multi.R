## ##############################################################
## Apply Welch's version of Student test to each row of a multivariate
## table
##
## Developed by Jacques van Helden <Jacques.van-Helden@univ-amu.fr>
## First version: September 2003
## Last modification: December 2003
##
## source(file.path(dir.util, "util_student_test_multi.R"))
source(file.path(dir.util, "multitesting_corrections.R"))

t.test.multi <- function (x, ## a matrix or data.frame
                          cl, ##  vector of classes
                          P.threshold=NA, ## select objects below a given threshold of P.value
                          E.threshold=NA, ## select objects below a given threshold of E.value
                          FDR.threshold=NA, ## select objects below a given threshold of FDR
			  robust.est=F, ## use robust estimators for central tendency and dispersion
			  verbosity=1, ## Level of verbosity
                          volcano.plot = T, ## Draw a volcano plot
                          alternative = "two.sided", ## Supported: "two.sided", "less", "greater"
                          ... ## additional parameters are passed to the function volcano.plot()
                          ) {
  ## Report starting time
  if (verbosity >= 1) { print (paste(date(), "Multiple t-test started", sep= " - ")) }
  
  ## Dimensions of the data set
  n <- nrow(x)
  p <- ncol(x)

  ## check the dimensions of the class vector
  if (length(cl) != p) {
    stop (paste ('The number of columns of the dataset (',
                 p,
                 ') should be the same as the length of the class vector (',
                 length(cl),
                 ')'))
  }

  classes <- unique(cl)
  if (robust.est) {
    means.1 <- apply(x[,cl==classes[1]],1,median,na.rm=T)
    means.2 <- apply(x[,cl==classes[2]],1,median,na.rm=T)
    iqr.1 <- apply (x[,cl==classes[1]], 1, IQR, na.rm=T)
    iqr.2 <- apply (x[,cl==classes[2]], 1, IQR, na.rm=T)
    sd.est.1 <- iqr.1/(qnorm(0.75) - qnorm(0.25))
    sd.est.2 <- iqr.2/(qnorm(0.75) - qnorm(0.25))
    var.est.1 <- sd.est.1^2
    var.est.2 <- sd.est.2^2
  } else {
    means.1 <- apply(x[,cl==classes[1]],1,mean,na.rm=T)
    means.2 <- apply(x[,cl==classes[2]],1,mean,na.rm=T)
    var.est.1 <- apply(x[,cl==classes[1]],1,var,na.rm=T)
    var.est.2 <- apply(x[,cl==classes[2]],1,var,na.rm=T)
    sd.est.1 <- sqrt(var.est.1)
    sd.est.2 <- sqrt(var.est.2)
  }


  ## 2012-05-02: restored means.diff = means.1 - means.2 for the sake
  ## of consistency with R funciton t.test()
  means.diff <- means.1 - means.2
  ##  means.diff <- means.2 - means.1
  ##  hist(means.diff,breaks=50)
  
  n.1 <- sum(cl == classes[1])
  n.2 <- sum(cl == classes[2])

  ## Calculate observed t value
  st.err.diff <- sqrt(var.est.1/n.1 + var.est.2/n.2)
  t.obs <- means.diff/st.err.diff
  ## Calculate degrees of freedom with Welch's formula  
  df.welch <- (var.est.1/n.1 + var.est.2/n.2)^2 / ((var.est.1/n.1)^2/(n.1-1) + (var.est.2/n.2)^2/(n.2-1))

  ## Calculate P-value and E-value
#  P.value.normal.approx <- 2*pnorm(abs(t.obs),lower.tail=F)
  if (alternative == "greater") {
    P.value <- pt(t.obs,df.welch,lower.tail=FALSE)
  } else if (alternative == "less") {
    P.value <- pt(t.obs,df.welch,lower.tail=TRUE)
  } else if (alternative == "two.sided") {
    P.value <- 2*pt(abs(t.obs),df.welch,lower.tail=F)
  } else {
    stop('Invalid alternative option for t.test.multi(). Supported: "two.sided", "less", or "greater"')
  }
  E.value <- P.value*nrow(x)
  sig <- -log(E.value)/log(10)

  multi.corr <- multitest.corrections(P.value, plots=FALSE)

  ## Collect all statistics in a data frame
  result <- data.frame(
                       means.1,
                       means.2,
                       means.diff,
                       var.est.1,
                       var.est.2,
                       sd.est.1,
                       sd.est.2,
                       st.err.diff,
                       t.obs,
                       df.welch,
#                       P.value.normal.approx,
                       P.value,
                       E.value,
		       sig)
  result$fwer <- multi.corr$multitest.table$fwer
  result$q.value <- multi.corr$multitest.table$qval.Storey
  result$fdr <- multi.corr$multitest.table$fdr
  result$rank <- multi.corr$multitest.table$rank

  if (robust.est) {
    names(result)[1:7] <- c(
                            paste("median.", classes[1], sep=""),
                            paste("median.", classes[2], sep=""),
                            "medians.diff",
                            paste("var.est.", classes[1], sep=""),
                            paste("var.est.", classes[2], sep=""),
                            paste("sd.est.", classes[1], sep=""),
                            paste("sd.est.", classes[2], sep="")
                            )
  } else {
    names(result)[1:7] <- c(
                            paste("mean.", classes[1], sep=""),
                            paste("mean.", classes[2], sep=""),
                            "means.diff",
                            paste("var.est.", classes[1], sep=""),
                            paste("var.est.", classes[2], sep=""),
                            paste("sd.est.", classes[1], sep=""),
                            paste("sd.est.", classes[2], sep="")
                            )
  }

  
  ## Filtering on P-value and E-value thresholds
  if (!is.na(P.threshold)) {
    result <- result[result$P.value < P.threshold,]
  }
  if (!is.na(E.threshold)) {
    result <- result[result$E.value < E.threshold,]
  }
  if (!is.na(FDR.threshold)) {
    result <- result[result$fdr < FDR.threshold,]
  }

  ## Report done time
  if (verbosity >= 1) { print (paste(date(), "Multiple t-test done", sep= " - ")) }

  ## Draw a volcano plot
  # if (volcano.plot) {
  #   x11(width=7, height=5)
  #   volcano.plot(result, ...)
  # }

  ##  plot(P.value.normal.approx,P.value,panel.first=grid(col='#0000ff'),log="xy")
  return (result)
}



## ##############################################################
## Volcano plot.
## Abcissa represents the mean difference,
## ordinate the significance (sig=-log10(E-value))
# volcano.plot <- function(multi.t, ## Must be the data frame obtained with t.test.multi()
#                          plot.col = c('grid'='#CCCCCC', 'lines'='blue', 'points'='#888888', 'positive'='#00BB00'), ## Colors for the plot
#                          x.max = max(abs(multi.t$means.diff)),
#                          y.max = NULL, ## max(c(1, multi.t$sig)),
#                          y.min = NULL, ## min(c(0, multi.t$sig)),
# #                         Y.score = "sig", ## Score to plot on the Y axis. Supported: sig (default), pval eval
#                          ... ## additional parameters are passed to the plot function
#                          ) {
# 
# 
# #  if (Y.score == "pval") {
# #    y.values <- as.vector(multi.t$P.value)
# #    ylab <- "nominal p-value"
# #    plot.log <- "y"
#     
#   ## } else if (Y.score == "eval") {
#   ##   y.values <- multi.t$E.value
#   ##   y.values[y.values <= 0] <- 1e-320
#   ##   ylab <- "e-value"
#   ##   plot.log <- "y"
# 
# #  } else {
#     y.values <- multi.t$sig
#     ylab <- "sig = -log10(E-value)"
#     plot.log <- ""
# #  }
# 
#   if (is.null(y.min)) {
#     y.min <- min(y.values)
#   }
# 
#   if (is.null(y.max)) {
#     y.max <- max(y.values)
#   }
# 
#     
#     
#   # 
#   # plot(multi.t$means.diff,
#   #      y.values,
#   #      xlab="Difference between the means",
#   #      ylab=ylab,
#   #      xlim=c(-x.max, x.max),
#   #      ylim=c(y.min, y.max),
#   #      col=plot.col['points'],
#   #      panel.first=grid(lty='solid', col=plot.col['grid']),
#   #      log=plot.log,
#   #      ...)
#   # 
#   # abline(v=0,col=plot.col['lines'], lwd=1)
#   # abline(h=0,col=plot.col['lines'],lwd=2)
#   # positive <- multi.t$sig >= 0
#   # lines(multi.t[positive, c('means.diff','sig')],col=plot.col['positive'], type='p',pch='+')
#   }
#   
#   
# 
# 
# ## ##############################################################
# ## Complement to the volcano plots: compare
# ## - standard error and significance
# ## - mean difference, standard error and significance
# stder.plot <- function(multi.t, ## Must be the data frame obtained with t.test.multi()
#                          plot.col = c('grid'='#CCCCCC', 'lines'='blue', 'points'='#888888', 'positive'='#00BB00'), ## Colors for the plot
#                          ... ## additional parameters are passed to the plot function
#                          ) {
#   x11(width=10, height=5)
#   par(mfrow=c(2,1))
#   
#   ## Draw standard error as a function of Welch significance
#   plot(multi.t[,c('st.err.diff','sig')],
#        xlab='Standard error of mean difference',
#        ylab='Welch significance (sig=-log10(E-value)',
#        col=plot.col['points'],
#        panel.first=grid(lty='solid',col=plot.col['grid']))
#   abline(h=0,col=plot.col['lines'],lwd=2)
#   points(multi.t[multi.t$sig > 0,c('st.err.diff','sig')],col=plot.col['positive'])
# 
#   ## Plot standard error versus mean differences, with circles
#   ## proportional to significance
#   plot(multi.t[,c('means.diff', 'st.err.diff')],
#        cex=multi.t$sig-min(multi.t$sig),
#        ylab='Standard error on the difference',
#        xlab='Difference between the means',
#        col=plot.col['points'],
#        panel.first=grid(lty='solid',col=plot.col['grid']))
#   abline(v=0, col=plot.col['lines'])
#   points(multi.t[multi.t$sig > 0,c('means.diff', 'st.err.diff')],
#          cex=multi.t[multi.t$sig > 0,'sig']-min(multi.t$sig),
#          col=plot.col['positive'])
# 
#   ## Restore default graphical parameters
#   par(mfrow=c(1,1))
# 
# }
# 
# 
# 
