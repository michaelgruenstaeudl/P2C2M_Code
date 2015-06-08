#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl","Noah Reid")
#email = "mi.gruenstaeudl@gmail.com"
#version = "2015.02.20.1700"

stats.main = function (path, xml.file, loci, resultData, prmFile) {
  # Descr:    coordinates executing of modules
  # Deps:     (various)
  # I/p:      path = absolute path to Working directory
  #           xml.file = name of infile
  #           loci = names of loci
  #           resultData
  #           prmFile

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> stats.main", fg="red"), "\n",
        sep="")
  }

##########################
# 1. Summarizing results #
##########################
  loghelpers.prntmngr("Summarizing results", uprFlg=T)
                                               
  # Setting outdata
  outData = list()
  # Setting alpha values
  alphaValues = c(0.1, 0.05, 0.01)
  for (val in alphaValues) {
    # Coordinating the calculation of the statistics
    results = stats.coord(resultData, loci, val)
    valStr = paste("alpha", as.character(val), sep="")
    outData[[valStr]] = results
  }

#####################
# 2. Writing legend #
#####################
  legend = "Differences between the posterior and the posterior predictive distributions. Each cell contains the following information in said order: mean, standard deviation, significance level. Codes in square brackets indicate the number of tails. Alpha values are automatically adjusted for the number of tails."
  outData$legend = legend
return(outData)
}


stats.coord = function(ind, loci, alpha) {
  # Descr:  generate significance tables
  # Deps:   stats.outmaker
  # I/p:    ind
  #         loci
  #         alpha
  # Note:   CV = coefficient of variance
  #         avgs = arithmetic means
  #         qntls = quantiles

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  slctn = get("P2C2M_flg_dscrStat", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> stats.coord", fg="red"), "\n",
        sep="")
  }

##############################
# 1. Setting number of tails #
##############################

  # T-tests for all descriptive statistics are two-tailed, because there is no
  # a priori reason in which direction they should differ.
  tailL = list()
  tailL$LCWT = "2"
  tailL$COAL = "2"
  tailL$NDC = "2"
  tailL$GSI = "2"
  ## T-tests for NDC should be left one-tailed, because trees not compliant 
  ## with the coalescent model have a higher number of deep coalescences.
  #NDCtail = "1l"
  ## T-tests for GSI should be right one-tailed, because trees not compliant 
  ## with the coalescent model have lower values.
  #GSItail = "1r"

#############################
# 2. Inferring significance #
#############################
  perGene = acrGenes = list()
  for (s in slctn) {
    perGene[[s]] = stats.perGene(ind[[s]]$dif, alpha, tailL[[s]])
    acrGenes[[s]] = stats.acrGenes(ind[[s]]$dif, alpha, tailL[[s]])
  }

#######################
# 3. Combining output #
#######################
  outList = list()
  # perGene output
  outList$perGene = sapply(perGene, cbind)
  rownames(outList$perGene) = c(loci)
  # acrossGene output
  outList$acrGenes = sapply(acrGenes, cbind)
  rownames(outList$acrGenes) = c("Sum", "Mean", "Median", "Mode", "CV")
  # Naming rows of output
  names = c()
  for (stat in slctn) {
    names = c(names, paste(stat, "[", tailL[[stat]], "]", sep=""))
  }
  colnames(outList$perGene) = colnames(outList$acrGenes) = names
  
  return(outList)
}


stats.perGene = function (diff, alpha, tail) {
  # Descr:    calculates statistics per gene
  # Deps:     statshelpers.qntls
  # I/p:      diff
  #           alpha
  #           tail

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> stats.perGene", fg="red"), 
        sep="")
  }

  # Structure of "diff" at this point:
  #        [gene1]   [gene2]  [gene3]
  # [1,]  -2.72395  -9.48870 45.79565
  # [2,]  14.97560  38.06380 88.60285
  # [3,]   8.76280  21.09895 50.51950
  qntls = statshelpers.qntls(diff, alpha, tail)
  sigSgns = statshelpers.sigsgn(qntls, tail)
  # Calculation of means
  # Note: "MARGIN=2" means "across a column"
  subCol1 = c(apply(diff, MARGIN=2, mean))
  # Calculation of stdv
  subCol2 = c(apply(diff, MARGIN=2, sd))
  # Assignment of signif. values
  subCol3 = c(sigSgns)
  # "\u00B1" is unicode sign of "plusminus"
  # Output format:   mean(+-1 SD)[:space:]signif.level
  #                  Mean and SD are rounded to two decimal places.
  outData = paste(round(subCol1, 2), 
               paste("(", "\u00B1", round(subCol2, 2), ")", sep=""), 
               subCol3)

  return(outData)
}


stats.acrGenes = function (diff, alpha, tail) {
  # Descr:    calculates statistics across genes
  # Deps:     statshelpers.qntls
  #           statshelpers.cv
  # I/p:      diff
  #           alpha
  #           tail

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> stats.acrGenes", fg="red"),
        sep="")
  }

  # Perform calculations per row (i.e per MCMC generation)
  acrossG_Sum = rowSums(diff)
  acrossG_Mean = rowMeans(diff)
  acrossG_Median = rowMedians(diff)
  acrossG_Mode = rowModes(diff)
  acrossG_CV = statshelpers.cv(diff)

  # "acrossG" is a matrix with four columns hereafter
  acrossG = cbind(acrossG_Sum, 
                  acrossG_Mean, 
                  acrossG_Median,
                  acrossG_Mode,
                  acrossG_CV)

  qntls = statshelpers.qntls(acrossG, alpha, tail)
  sigSgns = statshelpers.sigsgn(qntls, tail)

  subCol1 = c(apply(acrossG, MARGIN=2, mean))
  subCol2 = c(apply(acrossG, MARGIN=2, sd))
  subCol3 = c(sigSgns)
  outData = paste(round(subCol1, 2), 
               paste("(", "\u00B1", round(subCol2, 2), ")", sep=""), 
               subCol3)

  return(outData)
}
