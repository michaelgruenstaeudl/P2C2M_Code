#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl","Noah Reid")
#email = "mi.gruenstaeudl@gmail.com"
#version = "2014.10.29.1200"


statshelpers.qntls = function (ind, alpha, tail) {
  # Descr:  generates a distribution of quantiles;
  #         quantile levels are hardcoded in variable "qntlLevels"
  # Deps:   (none)
  #         ind = a data frame
  #         alpha
  #         tail

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> statshelpers.qntls",
        fg="red"), sep="")
  }

  #################################
  # Set quantile levels via alpha #
  #################################
  if (tail=="1l" | tail=="1r") {
    qntlLevels = c(alpha, 1-alpha)
  }
  # Adjusting alpha value for two-tailed test
  if (tail=="2") {
    qntlLevels = c((alpha/2), 1-(alpha/2))
  }

  qntls = t(apply(ind, MARGIN=2, quantile, qntlLevels, na.rm=T))
  #qntls = rbind(qntls, quantile(acrossGenes, qntlLevels, na.rm=T))

  return(qntls)
}


statshelpers.cv = function(diff) {
  # Descr:    Calculates the coefficient of variance
  # Deps:     (none)
  # I/p:      inD
  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> statshelpers.cv", fg="red"), sep="")
  }

  # Structure of "diff" at this point:
  #        [gene1]   [gene2]  [gene3]
  # [1,]  -2.72395  -9.48870 45.79565
  # [2,]  14.97560  38.06380 88.60285
  # [3,]   8.76280  21.09895 50.51950

#  inD = as.data.frame(diff)
#  inD = ifelse(is.nan(diff), NA, inD)
  Stdv = apply(diff, MARGIN=1, sd, na.rm=T)
  Mean = rowMeans(diff)
  cv = Stdv/Mean

  # TFL is CRITICAL, because there are occasional "Inf" in the matrix for 
  # descriptive statistics "NDC"
  is.na(cv) <- do.call(cbind, lapply(cv, is.infinite))
  return(cv)
}


statshelpers.diffrnce = function (post_dist, post_pred_dist) {
  ####################################
  # Function "statshelpers.diffrnce" #
  ####################################
  # Descr:    calculates the difference between post_distirical and post_pred_distulated
  # Deps:     -
  # I/p:      post_dist
  #           post_pred_dist

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> statshelpers.diffrnce", 
        fg="red"), sep="")
  }

  # Note: The difference is "post_distirical-post_pred_distulated". Since "post_distirical" is worse 
  # than "post_pred_distulated" whenever it does not conform to the coalescent model (i.e.
  # has larger values), significant differences will be more positive than 
  # non-significant differences.
  diff = post_dist - post_pred_dist
  # TFL converts from type "list" to type "double"; is important, because
  # is.infinite and is.nan can only work on type "double"
  diff = as.matrix(diff)

  # Removing diff. values that are infinite ("Inf")
  diff = ifelse(is.infinite(diff), NA, diff)
  # Removing diff. values that are "NaN" ("not a number", i.e. 0/0)
  diff = ifelse(is.nan(diff), NA, diff)

  return(diff)
}


statshelpers.sigsgn = function (qntls, tail) {
  ##################################
  # Function "statshelpers.sigsgn" #
  ##################################
  # Descr:    applies significance signs
  # Deps:     -
  # I/p:      qntls = a data frame of two columns

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> statshelpers.sigsgn", 
        fg="red"), sep="")
  }

  sigSgns = rep(0, length(qntls[,1]))

  # Schemes for one-tailed test
  if (tail=="1l") {
    # Read: "Whichever elements of quants are smaller than zero,
    #        are considered significant and, hence, receive a star 
    #        in the first column."
    sigSgns[which(qntls[,1] > 0)] = "*"
    sigSgns[which(qntls[,1] < 0)] = "n.s."
    # rare cases
    sigSgns[which(qntls[,1] == 0)] = "n.s."
  }
  if (tail=="1r") {
    sigSgns[which(qntls[,2] < 0)] = "*"
    sigSgns[which(qntls[,2] > 0)] = "n.s."
    # rare cases
    sigSgns[which(qntls[,2] == 0)] = "n.s."
  }

  # Scheme for two-tailed test
  if (tail=="2") {
    sigSgns[which(qntls[,1] > 0 | qntls[,2] < 0)] = "*"
    sigSgns[which(qntls[,1] < 0 & qntls[,2] > 0)] = "n.s."
    # rare cases
    sigSgns[which(qntls[,1] == 0 & qntls[,2] > 0)] = "n.s."
    sigSgns[which(qntls[,1] < 0 & qntls[,2] == 0)] = "n.s."
  }

  return(sigSgns)
}
