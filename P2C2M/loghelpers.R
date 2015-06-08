#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl","Noah Reid")
#email = "mi.gruenstaeudl@gmail.com"
#version = "2015.02.24.1300"


#loghelpers.prntmngr = function(inD, prmFHndl, tofileFlg, color) {
#  ###############################
#  # Function "loghelpers.prntmngr" #
#  ###############################
#  # Descr:  printing log; with canvas of pound signs
#  # Deps:   -
#  # I/p:    inD
#  #         prmFHndl
#  #         tofileFlg
#  #         color
#
#  metadataBool = get("P2C2M_flg_metadataBool", envir=P2C2M_globalVars)
#  verboseBool = get("P2C2M_flg_vrbBool", envir=P2C2M_globalVars)
#
#  if (verboseBool) {
#    # Logging to screen
#    cat("\n", xtermStyle::style(toupper(inD), fg=color), "\n", sep="")
#
#    # Writing to parameter file
#    if (metadataBool) {
#      if (tofileFlg) {
#        #tmpD = paste(rep("#", length(unlist(strsplit(inD, split="")))), collapse="")
#        #writeLines(paste("\n", "##", tmpD, "##", sep=""), prmFHndl)
#        writeLines(paste("\n## ", toupper(inD), "\n", sep=""), prmFHndl)
#        #writeLines(paste("##", tmpD, "##", sep=""), prmFHndl)
#      }
#    }
#  }
#}


loghelpers.prntmngr = function(inText, uprFlg=F, nwlFlg=T) {
  ##############################
  # Function "loghelpers.prnt" #
  ##############################
  # Descr:  print text to screen
  # Deps:   -
  # I/p:    inText

  #metadataBool = get("P2C2M_flg_metadataBool", envir=P2C2M_globalVars)
  verboseBool = get("P2C2M_flg_vrbBool", envir=P2C2M_globalVars)

  # Logging to screen
  if (verboseBool) {
      if (uprFlg & nwlFlg) {cat("\n", toupper(inText), "\n", sep="")}
      if (!uprFlg & nwlFlg) {cat("\n", inText, "\n", sep="")}
      if (!uprFlg & !nwlFlg) {cat(inText)}
  }

  ## Logging to parameter file
  #if (metadataBool) {
  #  if (tofileFlg) {
  #    writeLines(paste("\n## ", toupper(inText), "\n", sep=""),
  #               prmFHndl)
  #  }
  #}

}


#loghelpers.wrt.treedata = function(inD, prmFHndl, colFlg) {
#  ######################################
#  # Function "loghelpers.wrt.treedata" #
#  ######################################
#  # Descr:  write all kinds of tree data
#  # Deps:   -
#  # I/p:    inD
#  #         prmFHndl
#  #         colFlg
#
#  metadataBool = get("P2C2M_flg_metadataBool", envir=P2C2M_globalVars)
#
#  if (metadataBool) {
#    write.table(inD, prmFHndl, append=T, row.names=F, col.names=colFlg, 
#                quote=F, sep=",")
#  }
#}


#loghelpers.wrt2 = function(inD, prmFHndl) {
#  ##############################
#  # Function "loghelpers.wrt2" #
#  ##############################
#  # Descr:  write table under formatting scheme 2
#  # Deps:   -
#  # I/p:    inD
#  #         prmFHndl
#  #         colFlg
#
#  metadataBool = get("P2C2M_flg_metadataBool", envir=P2C2M_globalVars)
#
#  if (metadataBool) {
#    sink(prmFHndl, append=T)
#    cat("\n")
#    print(inD)
#    sink()
#  }
#}


#loghelpers.dbg = function(inList, varName, title) {
#  #############################
#  # Function "loghelpers.dbg" #
#  #############################
#  # Descr:  write debug data to file
#  # Deps:   -
#  # I/p:    inList
#  #         title
#
#  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
#
#  if (debugBool) {
#    cat("\n", xtermStyle::style(paste(
#        "DEBUG> *** SAVING", title, "TO FILE ***"), 
#        fg="red"), "\n", sep="")
#
#    # Readying the data to be saved
#    outVarName = paste("p2c2m.dbg.", varName, sep="")
#    assign(outVarName, inList)
#
#    # Specifying the outFileName
#    xmlFile = get("P2C2M_flg_xmlFile", envir=P2C2M_globalVars)
#    outFileName = paste(rmext(xmlFile), ".", outVarName, ".rda", sep="")
#
#    # Saving initial R object
#    if (!file.exists(outFileName)) {
#      save(outVarName, file=outFileName)
#    }
#
#    # Appending to existing R object
#    if (file.exists(outFileName)) {
#      old.objects = load(outFileName)
#      new.objects = c(old.objects, outVarName)
#      save(new.objects, file = outFileName)
#    }
#  #  if (file.exists(outFileName)) {
#  #    old.objects = load(outFileName)
#  #    save(list=c(old.objects, outVarName), file = outFileName)
#  #  }
#  }
#}
