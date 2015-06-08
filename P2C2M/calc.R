#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl","Noah Reid")
#email = "mi.gruenstaeudl@gmail.com"
#version = "2014.10.05.1700"


calc.lcwt = function(gTree, sTree, assoc, ploidy) {
  # Descr:    calculates the probability of a whole gene tree
  #           by calculating the probability of subtrees defined 
  #           by their MRCA (i.e. their node); done for all 
  #           nodes in a tree, values are summed up thereafter
  # Deps:     calc.parse
  #           calchelpers.brprob
  # I/p:      sTree
  #           gTree
  #           assoc
  #           ploidy
  # Note:     gtp = "gene tree probability" (a form of coalescent likelihood)

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n",xtermStyle::style("DEBUG> calc.lcwt",fg="red"),sep="")
  }

  handle = calc.parse(sTree, assoc)
  nodes = handle$nodes
  dmvD = handle$dmvD
  
  # DEBUGLINES:
  #cat("\nnodes\n"); print(nodes)
  #cat("\ndmvD\n"); print(dmvD)

  lnP = c()
  for(node in nodes) {
    lnP = c(lnP, log(calchelpers.brprob(sTree, gTree, assoc, 
                                        ploidy, dmvD, node)))
  }

  return(sum(lnP))
}


calc.ndc = function(gTree, sTree, assoc) {
  # Descr:  returns the number of deep coalescences for an entire tree
  # Deps:   calc.parse
  #         calchelpers.gtreeparse
  # I/p:    sTree
  #         gTree
  #         assoc
  # Note:   ndc = "number of deep coalescences"

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n",xtermStyle::style("DEBUG> calc.ndc",fg="red"),sep="")
  }

  handle = calc.parse(sTree, assoc)
  nodes = handle$nodes
  dmvD = handle$dmvD
  ndc = c()
  for (node in nodes) {
    tempData = calchelpers.gtreeparse(sTree, gTree, assoc, dmvD, node)
    ndc = c(ndc, tempData$lNd-1)
  }

  return(sum(ndc))
}


calc.gsi = function(gTree, assoc, singleAllele) {
  # Descr:  returns the Genealogical Sorting Index for a given tree
  # Deps:   -
  # I/p:    gTree
  #         assoc
  #         singleAllele
  # Note:   gsi = "genealogical sorting index"

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calc.gsi", fg="red"), sep="")
  }

  # matching the correct associations
  speciesAssoc = assoc[,1][match(gTree$tip.label, assoc[,2])]

## 1. Generate list of all species (i.e. unique elements of "speciesAssoc")
  species = unique(speciesAssoc)
  # Remove single allele species
  species = suppressWarnings(species[species != singleAllele])

## 2. Enclosing data in a list
  data = list()
  data$species = species
  data$gTree = gTree
  data$speciesAssoc = speciesAssoc

### 3. DEBUG mode
#    logdata = list(list(GSI=data))
#    names(logdata) = get("P2C2M_flg_repID", envir=P2C2M_globalVars)
#    loghelpers.dbg(logdata, "DescrStatsRawInput", "INPUT OF 'GSI'")

## 4. Calculating descriptive statistic
  # Calculating the gsi for every "group" (i.e. species other than 
  # single allele species)
  outL = c()
  for (sp in data$species) {
    outL = c(outL, genealogicalSorting::gsi(data$gTree, sp, data$speciesAssoc))
  }

  # Ensuring that the gsi values together add up to a number between 0 and 1
  outD = lapply(outL, function(x){x/length(outL)})
  outD = sum(unlist(outD))

  return(outD)
}


calc.parse = function(sTree, assoc) {
  # Descr:  parses species tree nodes for metric calculation
  # Deps:   calchelpers.dmvparse
  # I/p:    sTree
  #         assoc

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n",xtermStyle::style("DEBUG> calc.parse",fg="red"),sep="")
  }

  # DEBUGLINES:
  #cat("\nassoc\n"); print(assoc)
  #cat("\nsTree$tip.label\n"); print(sTree$tip.label)
  
  spNames = sTree$tip.label
  n_sp = length(spNames)
  #nBr = (2*n_sp)-1
  tiplist = list()
  for(i in 1:n_sp) {
    tiplist[[spNames[i]]] = assoc[which(assoc[,1]==spNames[i]),2]
  }
  dmvD = calchelpers.dmvparse(sTree, n_sp)                              # returns demographic info of the species tree
  n_tips_per_sp = lapply(tiplist, length)                               # calculate number of tips per species

  if(any(n_tips_per_sp==1)) {                                           # evaluate number of tips per species tree
    tmp = dmvD[-which(n_tips_per_sp==1),]                               # Remove dmv values for terminals (i.e., where n_tips_per_sp==1)
  } else {
    tmp = dmvD
  }
  nodes = tmp[,"node"]
  names(nodes) = NULL

  outD = list()
  outD$nodes = nodes
  outD$dmvD = dmvD

  return(outD)
}


calc.coal = function (gTree, sTree, stNames, assoc) {
  # Descr:  calculating the Ranalla&Yang index
  # Deps:   (various)     
  # I/p:    gTree
  #         sTree
  #         stNames
  #         assoc

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n",xtermStyle::style("DEBUG> calc.coal",fg="red"),sep="")
  }

## 1. Loading of gene trees and gene tree attributes
  # LEGACY: gTree.treeshape = apTreeshape::as.treeshape(gTree)
  # LEGACY: gTree.treeshape = apTreeshape::as.treeshape.phylo(gTree)
  # LEGACY: gTreeTaxa = gTree.treeshape$names
  gTreeTaxa = gTree$tip.label
  gTree.string = ape::write.tree(gTree)

## 2. Loading of species trees and species tree attributes
  sTreeTaxa = stNames
  sTree.string = sTree

## 3. Generating "species structure matrix"
  sp = unique(assoc[,1])
  taxa = assoc[,2]
  spStrMtrx = matrix(0, nrow=length(sp), ncol=length(taxa))
  rownames(spStrMtrx) = sp
  colnames(spStrMtrx) = taxa
  for (s in sp) {
    spStrMtrx[s, assoc[,2][which(assoc[,1]==s)]] = 1
  }

## 4. Enclosing data in a list
  data = list()
  data$gTree.string = gTree.string
  data$sTree.string = sTree.string
  data$gTreeTaxa = gTreeTaxa
  data$sTreeTaxa = sTreeTaxa
  data$spStrMtrx = spStrMtrx

  # DEBUGLINES:
  #cat("\ndata\n"); print(data)

### 5. DEBUG mode
#  logdata = list(list(COAL=data))
#  names(logdata) = get("P2C2M_flg_repID", envir=P2C2M_globalVars)
#  loghelpers.dbg(logdata, "DescrStatsRawInput", "INPUT OF 'COAL'")

## 6. Calculating descriptive statistic
  ray = phybase::loglikeSP(data$gTree.string, data$sTree.string,
                           data$gTreeTaxa, data$sTreeTaxa, 
                           data$spStrMtrx, strict=F)

  return(ray)
}
