rmext = function (instring, delimiter="\\.") {
  # function for returning instring without extension (i.e. ".txt" or ".trees")
  # get all instances of "."
  alist = stringr::str_locate_all(instring, delimiter)
  # get position of last element
  pos2 = tail(alist[[1]], 1)[1]-1
  return(substr(instring, 1, pos2))
}

str2lst = function (instring) {
  # function for returning a list of letters from a string
  unlist(strsplit(instring, ""))
}

is.integer0 = function (instring) {
  # function to check if return is integer(0)
  is.integer(instring) && length(instring) == 0L
}

is.character0 = function (instring) {
  # function to check if return is character(0)
  is.character(instring) && length(instring) == 0
}

frmtMntse = function (inNum, mntse) {
  # function to control the significands of a number
return(as.numeric(format(round(inNum, mntse), nsmall=mntse)))
}

rmOutlrs = function(ind) {
  # function to remove outlier values
  outd = ind[!ind %in% boxplot.stats(ind)$out]
  return(outd)
}

rtnlc = function (m, x) {
  # function for returning listcolumn x in a matrix of lists 
  # ("ReTurN List Column")
  # m = inMatrix
  aFun = function(r, c, m, x) {as.numeric(unlist(strsplit(m[r, c],",")))[x]}
  vect_aFun = Vectorize(aFun, vectorize.args = c('r','c'))
  return(outer(1:nrow(m), 1:ncol(m) , FUN = vect_aFun, m, x))
}

Mode = function(x) {
  # function for calculating the mode
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

rowMedians = function(x) {
  # function for calculating the median of all row values
  # analogous to rowSums and rowMeans
  apply(x, MARGIN=2, median, na.rm=T)
}

rowModes = function(x) {
  # function for calculating the mode of all row values
  # analogous to rowMeans and rowMedians
  apply(x, MARGIN=2, Mode)
}
