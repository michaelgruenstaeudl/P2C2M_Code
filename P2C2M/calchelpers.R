#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl","Noah Reid")
#email = "mi.gruenstaeudl@gmail.com"
#version = "2014.07.11.1700"

calchelpers.descend = function (sTree, gTree, assoc, cNode) {
##################################
# Function "calchelpers.descend" #
##################################
# Descr:    returns a tree containing descendants of all gene lineages 
#           that pass through a node; here this function is used to 
#           obtain the gene tree that is part of a species tree
# Deps:     calchelpers.nodetips
# I/p:      sTree
#           gTree
#           assoc
#           cNode = coalescence node

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.descend", fg="red"),
        sep="")
  }

    spNames = sTree$tip.label
    if (!is.na(cNode) && cNode != (length(spNames)+1))                  # Unless cNode is 'NA' AND unless the coalescence node is the root, do ...
        {
        subtreeTaxa = spNames[calchelpers.nodetips(sTree, cNode)]
        gTreeDscr = c()
        for(taxon in subtreeTaxa)
            {
            nodesRaw = assoc[which(assoc[,1] == taxon),2]
            if (is.character(nodesRaw[1]) == T)
                {nodes = nodesRaw}
            if (is.character(nodesRaw[1]) == F)
                {nodes = as.numeric(nodesRaw)}
            gTreeDscr = c(gTreeDscr, nodes)
            }
        taxaIndex = match(gTreeDscr, gTree$tip.label)
        # TFL removes all terminal branches of the gene tree 
        # that are not in the subtree; gTree must be in DNAbin-format
        subtree = ape::drop.tip(gTree, gTree$tip.label[-taxaIndex])
        }
    else {subtree = gTree}

return(subtree)
}


calchelpers.nodetips = function(tree, node) {
  # Descr:    returns the tip numbers for a given node
  # Deps:     calchelpers.nodetips
  # I/p:      tree = a phylog. tree with numbered internal nodes and numbered tips
  #           node = a node in that tree, specified as an integer

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.nodetips", fg="red"),
        sep="")
  }

    # DEBUGLINES:
    #cat("\ntree$edge\n"); print(tree$edge)
    #cat("\nnode\n"); print(node)

    n_tips = length(tree$tip.label)
    if (node <= n_tips) {node}
    #if (node <= n_tips) {                                              # Potential improvements than line above
    #    outD = node
    #    outD
    #}
    else 
        {
        outD = numeric()
        k = tree$edge[which(tree$edge[,1] == node),2]                   # CURRENT PROBLEM UNDER LCWT: When diff. allele numbers, node nicht in tree$edge[,1]
        for (j in k)
            {
            if (j <= n_tips) {outD = c(outD, j)}
            else {outD = c(outD, calchelpers.nodetips(tree, j))}
            }
        outD        
        }
}


calchelpers.dmvparse = function(sTree, nSp) {
  # Descr:    parsing demographic data (dmv and dmt) of a branch
  # Deps:     -
  # I/p:      sTree
  #           nSp

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.dmvparse", fg="red"),
        sep="")
  }
  
    # DEBUGLINES:
    #cat("\nsTree\n"); print(sTree)
    #cat("\nsTree$edge\n"); print(sTree$edge)
    #cat("\nsTree$edge.length\n"); print(sTree$edge.length)
    #cat("\nsTree$dmv\n"); print(sTree$dmv)

    dmvD = cbind(sTree$edge[,2], sTree$edge.length)                     # generating rows of node number - branch length pairs
    dmvD = rbind(dmvD, c((nSp+1), Inf))                                 # adding another branch of length Inf
    dmvD = cbind(dmvD, sTree$dmv)                                       # adding dmv values
    dmvD = dmvD[order(dmvD[,1]),]                                       # order the matrix by the first column
    stBt = ape::branching.times(sTree)                                  # calc branching times of the species tree (via ape-function)

    # TFL may not be necessary, as stBt may already be sorted
    stBt = stBt[order(as.numeric(names(stBt)))]                         # sort the branching times by their names (which are numbers)

    pre = structure(rep(0, nSp), .Names=c(1:nSp))                       # add x zeros to the beginning of branching times list, where x = number of species in sTree
    stBt = c(pre, stBt)

    dmvD = cbind(dmvD, stBt)                                            # add the column "stBt" to the matrix "dmvD"
    colnames(dmvD) = c("node", "length", "dmv", "sbt")
    rownames(dmvD) = c(1:length(stBt))

return(dmvD)
}

 
calchelpers.brprob = function(sTree, gTree, assoc, ploidy, dmvD, node) {
# Descr:    organising likelihood calculations for particular branches
#           Calculate prob of gene tree part, which is a subset 
#           of the species tree.
#           More general: Calculate the prob of subtrees defined 
#           by their MRCA (i.e. their node). 
# Deps:     calchelpers.gtreeparse
#           calchelpers.probcalc
# I/p:      sTree
#           gTree
#           assoc
#           ploidy
#           dmvD = demographic data
#           node
# Note:     branching times = distance from each node to the tips, 
#            under the assumption that the tree is ultrametric

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.brprob", fg="red"),
        sep="")
  }

    ind = calchelpers.gtreeparse(sTree, gTree, assoc, dmvD, node)

    gBt = ind$gBt
    fBt = ind$fBt
    lBt = ind$lBt
    fNd = ind$fNd
    lNd = ind$lNds

    if(fBt!=0) {prob=1}
    else
        {
        # Retain those gTree branching times (i.e. column 1) whose node names 
        # (i.e. column 2) fall within "fBt" and "lBt"
        gBt = gBt[gBt[,1]>=fBt & gBt[,1]<=lBt,]
        # Add lNd row at bottom of gBt-matrix
        gBt = rbind(gBt, c(lBt, lNd))
        # Perform actual prob calculation
        prob = calchelpers.probcalc(gBt, dmvD, ploidy, node, fBt, lBt, fNd, lNd)
        }

return(prob)
}


calchelpers.gtreeparse = function(sTree, gTree, assoc, dmvD, node) {
# Descr:    parsing branches for likelihood calculations
# Deps:     calchelpers.descend
# I/p:      sTree
#           gTree
#           assoc
#           dmvD
#           node
# Note:     branching times = distance from each node to the tips, 
#           under the assumption that the tree is ultrametric
#           gBt = gene tree branching times (plural!)
#           fBt = first branching time
#           lBt = last branching time
#           fNd = first node
#           lNd = last node

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.gtreeparse", fg="red"),
        sep="")
  }

    # returns the subtree that starts at node <node>
    subtree = calchelpers.descend(sTree, gTree, assoc, node)
    # infers the branching times for a subtree (here the gene tree)
    gBt = sort(ape::branching.times(subtree))
    # generate a matrix with branching times in first column and node 
    # names in second column
    gBt = c(0, gBt)
    gBt = cbind(gBt, length(gBt):1)
    # get the branching times for a node in the species tree 
    fBt = dmvD[node,"sbt"]
    # get the branching times for a node in the species tree + the branch 
    # lengths for that node
    lBt = (dmvD[node,"sbt"] + dmvD[node,"length"])
    # get those node IDs (column 2), whose gene tree branching times (column 1) 
    # are equal the largest of the species tree branching times
    fBtMax = max(gBt[gBt[,1]<=fBt,1])
    fNd = gBt[gBt[,1]==fBtMax,2]
    # get those node IDs (column 2), whose  gene tree branching times (column 1) 
    # that equal the largest of the species tree branching times + branch lengths
    lBtMax = max(gBt[gBt[,1]<=lBt,1])
    lNd = gBt[gBt[,1]==lBtMax,2]

    outd = list("gBt"=gBt, "fBt"=fBt, "lBt"=lBt, "fNd"=fNd, "lNd"=lNd)

return(outd)
}



calchelpers.probcalc = function(gBt, dmvD, ploidy, node, fBt, lBt, fNd, lNd) {
# Descr:    calculating likelihoods for particular branches
# Deps:     -
# I/p:      gBt
#           dmvD
#           ploidy
#           node
#           fBt
#           lBt
#           fNd
#           lNd
# Note:     The input variable "gBt" contains two columns: 
#           branching time differences and node ID
#           branching times = distance from each node to the tips, 
#           under the assumption that the tree is ultrametric

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.probcalc", fg="red"),
        sep="")
  }

# 1. Calculating differences in branching times
    value = length(gBt[,1])
    # waitTms = branching time differences between two branches
    waitTms = gBt[2:value,1] - gBt[1:value-1,1]
    # append node IDs again (so that "waitTms" mimicks "gBt"), 
    # with the exception of last node ID in list
    waitTms = cbind(waitTms, gBt[1:value-1,2])
# 2. Calculating lists of values
    dmv = (2*dmvD[node,"dmv"])*ploidy
    # lambda is a list of values, each calculated according to the following formula
    lambda = (waitTms[,2]*(waitTms[,2]-1))/dmv
    # exponent is a list of values; Euler-constant to the power of (-lambda*waitTms[,1])
    exponent = exp(-lambda*waitTms[,1])
# 3. Calculating products for each list
    # Calculate product of exponent list, given that lambda is never 0
    exponent = prod(exponent[lambda!=0])
    # Calculate product of lambda list, however disregarding last list element
    lambda = prod(lambda[1:(length(lambda)-1)])
# 4. Calculating overall probability
    prob = lambda * exponent

return(prob)
}

