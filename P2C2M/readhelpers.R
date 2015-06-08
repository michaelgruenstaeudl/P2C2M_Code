#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl","Noah Reid")
#email = "mi.gruenstaeudl@gmail.com"
#version = "2015.04.02.1500"


readhelpers.extract_br_and_meta = function(tree, locus, sTreeFlg) {
  # Descr:  extracting tree and branch rate vectors, 
  #         and order vectors appropriately
  # Deps:   -
  # I/p:    tree
  #         locus
  #         sTreeFlg

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> readhelpers.extract_br_and_meta", fg="red"),
        sep="")
  }
  
  # DEBUGLINE:  print(tree)

###########################
# 1. Generate node labels #
###########################
  n_nodes = length(unlist(strsplit(tree, "\\[")))                       # Count the number of nodes by counting number of square brackets (b/c each taxon is followed by metadata)
  if (sTreeFlg==F) {
    n_internNodes = (n_nodes+1)/2                                       # Count the number of terminals by dividing the number of nodes by two
    labels_internNodes = (n_internNodes+2):((2*n_internNodes)-1)        # Generate list of terminal labels
  }
  if (sTreeFlg==T) {
    n_internNodes = n_nodes/2
    labels_internNodes = (n_internNodes+1):((2*n_internNodes)-1)
  }
  last_internNode = labels_internNodes[length(labels_internNodes)]

##############################################
# 2. Add internal node labels to tree string #
##############################################
  for (label in labels_internNodes) {
    if (sTreeFlg==F) {
      repl_str = paste(")", label, ":[", sep="")
      tree = sub("):\\[", repl_str, tree)
    }

    if (sTreeFlg==T) {
      repl_str = paste(")", label, "[", sep="")
      tree = sub(")\\[", repl_str, tree)
    }
  }
  tree_wo_MetaD = gsub("\\[[^]]*\\]", "", tree)                         # Remove any metadata from tree string
  pure_treestring = ape::read.tree(text=tree_wo_MetaD)                  # Convert tree from string to phylo object
  # DEBUGLINES:
  #cat("\npure_treestring\n"); print(pure_treestring)

#####################################################
# 3. Extract all lines with "dmv/dmt" | "rate" info #
#####################################################
  wntdList = strsplit(tree, "\\[|\\]")                                  # Split tree string at every square bracket
  if (sTreeFlg==F) {wntdLines = grep("rate", unlist(wntdList))}         # Extract those list elements that contain keyword "dmv/dmt" | "rate"
  if (sTreeFlg==T) {wntdLines = grep("dmv", unlist(wntdList))}
  metaD_raw = unlist(wntdList)[wntdLines]

  metaD_n_cols = length(unlist(strsplit(metaD_raw[1], ",")))            # Count seperate pieces of info (i.e. columns) related to "dmv/dmt" | "rate" values
  if (sTreeFlg==F) {metaD_unformat = gsub("&|rate=|\\{|\\}", "", metaD_raw)}  # Parse "dmv/dmt" | "rate" values as raw numbers
  if (sTreeFlg==T) {metaD_unformat = gsub("&|dm.=|\\{|\\}", "", metaD_raw)}

  metaD_unformat = gsub(paste(locus, ".", sep=""), "", metaD_unformat)  # Necessary for starBeast.v.1.8. and higher: Remove locus name plus trailing period from metaD_unformat element

############################################################
# 4. Extract node names and ultrametric branch length info #
############################################################
  brLens = gsub("\\[[^]]*\\]", "\\[\\]", tree)                          # Remove any metaD_unformat (i.e., anything inside square brackets) from tree string
  brLens = unlist(strsplit(brLens, ",|)"))                              # Split by comma or closed parenthesis
  brLens = gsub("\\(|\\)|;|\\[|\\]", "", brLens)                        # Remove all parentheses or square brackets
  if (sTreeFlg==T) {
    #cat("\nNote: 'NA' being added.\n")
    brLens[length(brLens)] = paste(brLens[length(brLens)], NA, sep=":") # In species trees, the root branch receives a "NA" (why?)
  }
  # DEBUGLINES:
  #cat("\nmetaD_n_cols\n"); print(metaD_n_cols)
  #cat("\nmetaD_unformat\n"); print(metaD_unformat)
  #cat("\nbrLens\n"); print(brLens)

##########################
# 5. Metadata extraction #
##########################
  metaD_formatted = readhelpers.formatmeta(metaD_unformat, metaD_n_cols, brLens, sTreeFlg)

#####################################
# 6. Generate correspondence matrix #
#####################################

  if (sTreeFlg==F) {
    correspTable_internNodes = cbind(pure_treestring$node.label[-1],    # Remove root node in gene trees
                                     (n_internNodes+2):last_internNode) # Generate correspondence table (to understand relation btw. tree terminals and extracted brLen info)
  }
  if (sTreeFlg==T) {
    correspTable_internNodes = cbind(pure_treestring$node.label,
                                     (n_internNodes+1):last_internNode)
  }
  correspTable_terminNodes = cbind(pure_treestring$tip.label, 
                                   1:n_internNodes)
  correspTable = rbind(correspTable_internNodes, correspTable_terminNodes) 
  rownames(correspTable) = correspTable[,2]                             # Name the rows by the second column

  # DEBUGLINES:
  #cat("\ncorrespTable\n"); print(correspTable)
  
  ord_correspTable = correspTable[as.character(pure_treestring$edge[,2]),]  # Orders the items in "correspTable" to the order specified in "pure_treestring"
  if (sTreeFlg==T) {
    ord_correspTable = ord_correspTable[which(ord_correspTable[,1]!=""),]  # Remove metaD_formatted entries that do not contain a branch name (i.e., that represent nodes)
  }
  # DEBUGLINES:
  #cat("\nord_correspTable\n"); print(ord_correspTable)

################################
# 7. Extract relevant metadata #
################################
  outD = metaD_formatted[ord_correspTable[,1],]                         # Obtain these entries in "metaD_formatted" that are labelled as given in "ord_correspTable[,1]"
  if (sTreeFlg==T) {
    outD = rbind(outD, metaD_formatted[length(metaD_formatted[,2]),])   # Add the last entry of the matrix "metaD_formatted"
  }
  rownames(outD) = NULL                                                 # remove rownames
  # DEBUGLINES:
  #cat("\noutD\n"); print(outD)

  return(outD)
}


readhelpers.formatmeta = function(metaD_unformat, metaD_n_cols, brLens, sTreeFlg) {
  # Descr:  formatting metadata
  # Deps:   -
  # I/p:    metaD_unformat
  #         metaD_n_cols
  #         brLens
  #         sTreeFlg

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> readhelpers.formatmeta", fg="red"), 
        sep="")
  }
  
  # DEBUGLINES:
  #cat("\nmetaD_unformat\n"); print(metaD_unformat)
  #cat("\nmetaD_n_cols\n"); print(metaD_n_cols)

################################
# 1. Setting up empty matrices #
################################
  if(metaD_n_cols==1) {
    metaD_formatted = array(dim=c(length(metaD_unformat), metaD_n_cols+2))
  }
  if(sTreeFlg && metaD_n_cols==3) {
    metaD_formatted = array(dim=c(length(metaD_unformat), metaD_n_cols+3))
  }

###################################################################
# 2. Save node name, branch length and demographic info to matrix #
###################################################################
  for(i in 1:length(metaD_unformat)) {
    brInfo = unlist(strsplit(brLens[i], ":"))
    metaInfo = unlist(strsplit(metaD_unformat[i], ","))

    if(metaD_n_cols==1) {
      metaD_formatted[i,] = c(brInfo, metaInfo)
    }
    if(sTreeFlg && metaD_n_cols==3) {
      dmvValue = mean(as.numeric(c(metaInfo[2], metaInfo[3])))          # Calculation of dmvValue as mean between dmv_start and dmv_end
      metaD_formatted[i,] = c(brInfo, metaInfo, dmvValue)               # Fill the first two columns of the matrix "metaD_formatted" with node name and branch length info,and fill the remaining three columns with the metadata info
    }
  }

  # DEBUGLINES:
  #cat("\nmetaD_formatted\n"); print(metaD_formatted)

##########################
# 3. Formatting metadata #
##########################
  metaD_formatted = metaD_formatted[which(metaD_formatted[,1]!=""),]    # Remove all those metaD_formatted entries that do not contain a branch name (i.e. that represent nodes)

  if(sTreeFlg && metaD_n_cols==1) {
      colnames(metaD_formatted) = c("br", "length", "dmv")              # Name the columns of matrix "metaD_formatted" depending on the number of metadata infos
  }
  if(!sTreeFlg && metaD_n_cols==1) {
      colnames(metaD_formatted) = c("br", "length", "rate")
  }
  if(sTreeFlg && metaD_n_cols==3) {
      colnames(metaD_formatted) = c("br", "length", "dmt",
                                    "dmv_start", "dmv_end", "dmv")
  }

  rownames(metaD_formatted) = metaD_formatted[,1]                       # name the rows of of matrix "metaD_formatted" by the node name

  return(metaD_formatted)
}


readhelpers.parsetrees = function(inFn) {
  # Descr:  reading in trees
  # Deps: -
  # I/p:  inFn

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n",xtermStyle::style("DEBUG> readhelpers.parsetrees",fg="red"),sep="")
  }

  outD = list()
  tree = ape::read.nexus(inFn)
  outD[["full"]] = tree                                                 # Pass on the full tree data

  lines = scan(file=inFn, what="", sep="\n", quiet=T)                   # Load inFn line by line
  outD[["treestrings"]] = lines[grep("tree STATE", lines)]              # Extract all those lines that contain keyword "tree STATE"

  return(outD)
}



readhelpers.makeCFtable = function(ind) {
  # Descr:  generating a constituent-frequency table
  # Deps: -
  # I/p:  ind
  #
  # The command will perform the following conversion:
  # FROM:
  #         Var1            Freq
  #    [1,] "alpinus"        "4"
  # TO:
  #        spec              V2 
  #   Var1 "alpinus"         "1"
  #   Var1 "alpinus"         "2"
  #   Var1 "alpinus"         "3"
  #   Var1 "alpinus"         "4"

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> readhelpers.makeCFtable", fg="red"), 
        sep="")
  }

  ind = as.matrix(as.data.frame(table(ind[,1])))
  aList = c()
  for (i in 1:length(ind[,1])) {
    aList = c(aList, rep(ind[i,1], times=ind[i,2]))
  }
  outd = as.matrix(as.data.frame(cbind(aList, 1:length(aList))))

  return(outd)
}
