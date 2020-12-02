################################################################
## Apply various multi-testing corrections on a list of P-values, and
## plot illutrative graphs.
##
## Developed by Jacques van Helden <jvhelden@univ-amu.fr>
## First version: Feb 10, 2012
## Last modification: Feb 10, 2012
##
##   source(file.path(dir.util, "multitesting_corrections.R"))

## For the sake of comparison / validation, we also generate the
## graphs with Storey's library for the q-value (downloaed from
## http://genomics.princeton.edu/storeylab/qvalue/), which includes an
## estimation of the proportion of true and alternative hypotheses.
library(qvalue)

multitest.corrections <- function(pval, ## list of p-values
                                  names=NULL, ## Names associated to the pval list (will be exported to result table if the attribute file.prefix is specified)
                                  plots=FALSE, ## If true, generate illutrsative plots
                                  plot.pch=c(pval=2, eval=1, fwer=20, qval.0=3,qval.Storey=4), ## point type (character) associated to each statistics for the plot
                                  plot.col=c(pval='#000000', eval='#BBBBBB', fwer='#444444', qval.0='#666666',qval.Storey='#888888'), ## color for the plot
                                  plot.elements=c("pval", "eval", "fwer", "qval.Storey", "qval.0"), ## Elements to draw on the plots
                                  lambda=0.5, ## lambda parameter for estimating pi0 = m0/m1
                                  threshold=0.05, ## Threshold of significance for the plots
                                  main='Multitesting corrections', ## PRefix for the main title of the plots
                                  file.prefix=NULL ## Prefix for storing the plots and results
) {
  
  m <- length(pval)
  
  ## Create a table for storing the results
  multitest.table <- data.frame(pval=pval)
  if (!is.null(multitest.table)) {
    rownames(multitest.table) <- names
  }    
  multitest.table$rank <- rank(pval)
  
  ## E-value
  multitest.table$eval <- pval * m
  
  ## FDR
  multitest.table$fdr <- p.adjust(pval, "fdr")
  
  ## Family-wise error rate (fwer)
  multitest.table$fwer <- pbinom(q=0, size=m, prob=multitest.table$pval, lower.tail=FALSE)
  #  multitest.table$fwer <- 1 - ( 1 - multitest.table$pval)^m
  
  ## Compute q-value with Storey's library
  qobj <- qvalue(multitest.table$pval, lambda=lambda)
  multitest.table$qval.Storey <-  qobj$qvalues
  
  ## q-value_0 (q-value corresponding to the case where all the features
  ## are trully null).
  multitest.table$qval.0 <- multitest.table$pval * m/multitest.table$rank
  head(multitest.table, n=10)
  
  ## Export the multitest correction table
  if (!is.null(file.prefix)) {
    export.object(x=multitest.table,file.prefix=paste(file.prefix, '_multitest_scores', sep=''), export.formats='table')
  }
  
  ## Plot reproducing Fig 1 from Storey and Tibshirani (2003),
  ## illustrating the estimation of the parameter PI0 = m0/m1
  if (plots) {
    x11(width=7,height=6)
    distrib <- hist(pval, breaks=seq(from=0,to=1, by=0.05), plot=TRUE,
                    main=paste(main, "p-value distribution", sep=' ; '))
  } else {
    distrib <- hist(pval, breaks=seq(from=0,to=1, by=0.05), plot=FALSE)
  }
  m.per.class <- sum(distrib$counts)/length(distrib$counts)
  
  ## Estimate numbers of null and alternative hypotheses
  m0.est <- min(m, sum(pval >= lambda) / (1-lambda)) ## m0 cannot be larger than m
  m1.est <- m - m0.est
  pi0 <- m0.est/m
  
  #  print(c(lambda, m0.est))
  m0.per.class <- m0.est / length(distrib$count)
  
  if (plots) {
    abline(h=m.per.class, col='black', lty='dashed', lwd=2)
    abline(h=m0.per.class, col='black', lty='dotted', lwd=2)
    arrows(lambda, m.per.class*1.5, lambda, m.per.class*1.1, angle = 30, length=0.05, code = 2, col='black', lwd = 2)
    legend('topright', bty='o', bg='white', border=NA,
           legend=(c(paste("lambda =", lambda),
                     paste("m0.est =", m0.est),
                     paste("pi0 =", round(pi0, digits=3)))
           ))
    if (!is.null(file.prefix)) {
      export.plot(file.prefix=paste(file.prefix, "_pval_distrib"), export.formats=export.formats.plots, width=7,height=6);
    }
  }
  
  ################################################################
  ## Draw plots
  
  if (plots) {
    
    ## Draw plots of the function qvalue::qplot
    x11(width=8, height=6)
    qplot(qobj)
    if (!is.null(file.prefix)) {
      export.plot(file.prefix=paste(file.prefix, "_qvalue_plots"), export.formats=export.formats.plots, width=8,height=6);
    }
    
    ## Plot the multi-testing corrected statistics as a function of P-value
    pval.min = min(pval)
    X11(width=7, height=9)
    plot(multitest.table$pval,
         multitest.table$pval,
         main=paste(main, "significance scores", sep=' ; '),
         xlab=("P-value"),
         ylab=("Multi-testing corrected statistics"),
         log='xy',
         panel.first=c(
           abline(h=10^c(-10:3),col="#CCCCCC"),
           abline(v=10^c(-10:0),col="#CCCCCC"),
           abline(h=1,col="#666666")
         ),
         col=plot.col["pval"], pch=plot.pch["pval"],
         xlim=c(pval.min, 1),
         ylim=c(pval.min, m))
    lines(multitest.table$pval, multitest.table$qval.Storey, col=plot.col["qval.Storey"], pch=plot.pch["qval.Storey"], type='p')
    lines(multitest.table$pval, multitest.table$qval.0, col=plot.col["qval.0"], pch=plot.pch["qval.0"], type='p')
    lines(multitest.table$pval, multitest.table$eval, col=plot.col["eval"], pch=plot.pch["eval"], type='p')
    lines(multitest.table$pval, multitest.table$fwer, col=plot.col["fwer"], pch=plot.pch["fwer"], type='p')
    legend('topleft', bg='#FFFFFF', bty="o",
           legend=plot.elements,
           col=plot.col[plot.elements],
           pch=plot.pch[plot.elements],
    )
    abline(h=threshold, col='black', lty='dashed')
    
    ## Save plot file
    if (!is.null(file.prefix)) {
      export.plot(file.prefix=paste(file.prefix, "_pvalue_corrections"), export.formats=export.formats.plots, width=7,height=9);
    }
    
    ## Plot the number of significant features as a function of the control criterion
    # X11(width=7, height=7)
    # plot(multitest.table$rank,
    #      multitest.table$pval,
    #      main=paste(main, "significant features", sep=' ; '),
    #      xlab=("Number of significant features"),
    #      ylab=("Multi-testing corrected statistics"),
    #      log='xy',
    #      panel.first=c(
    #        grid(equilogs=F, lty='solid', col='#CCCCCC'),
    #        abline(h=10^c(-11:3), col='#CCCCCC'),
    #        abline(h=1,col="#666666")
    #      ),
    #      col=plot.col["pval"], pch=plot.pch["pval"],
    #      xlim=c(1,nrow(multitest.table)),
    #      ylim=c(pval.min, m))
    # abline(h=threshold, col='black', lty='dashed')
    # lines(multitest.table$rank, multitest.table$qval.Storey, col=plot.col["qval.Storey"], pch=plot.pch["qval.Storey"], type='p')
    # lines(multitest.table$rank, multitest.table$fwer, col=plot.col["fwer"], pch=plot.pch["fwer"], type='p')
    # lines(multitest.table$rank, multitest.table$eval, col=plot.col["eval"], pch=plot.pch["eval"], type='p')
    # legend('topleft', bg='#FFFFFF', bty="o",
    #        legend=plot.elements,
    #        col=plot.col[plot.elements],
    #        pch=plot.pch[plot.elements],
    # )
    # 
    ## Save result file
    # if (!is.null(file.prefix)) {
    #   export.plot(file.prefix=paste(file.prefix, "_nb_signif_features"), export.formats=export.formats.plots, width=7,height=7);
    # }
    
    
    
  }
  
  
  ## Prepare the result
  result <- list()
  result$qobj <- qobj
  result$multitest.table <- multitest.table
  result$m <- m
  result$lambda <- lambda
  result$m0.est <- m0.est
  result$m1.est <- m1.est
  result$pi0 <- pi0
  result$nb.signif <- c(
    pval=sum(multitest.table$pval <= threshold), 
    bonferoni=sum(multitest.table$eval <= threshold),
    eval.le.1=sum(multitest.table$eval <= 1),
    fwer=sum(multitest.table$fwer <= threshold),
    qval.0=sum(multitest.table$qval.0 <= threshold),
    qval.Storey=sum(multitest.table$qval.Storey <= threshold),
    fdr=sum(multitest.table$fdr <= threshold)
  )
  
  return(result)
}