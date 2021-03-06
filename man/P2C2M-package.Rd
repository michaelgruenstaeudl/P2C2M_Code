\name{P2C2M-package}
\alias{P2C2M-package}
\alias{P2C2M}
\docType{package}

\title{
  Posterior Predictive Checks of Coalescent Models
}

\description{
  \pkg{P2C2M} provides functions to read default output from BEAST 
  (Drummond and Rambaut 2007) and *BEAST (Heled and Drummond 2010) and 
  conduct posterior predictive checks of coalescent models (Reid et al. 
  2014) with the help of data simulation and summary statistics under 
  various settings.
}

\note{

  \bold{Installation Instructions}

  To use \pkg{P2C2M}, the default version of Python must be set to 
  Python 2.7. Users of unix-like operating systems can insure that 
  this requirement is fulfilled by setting the following alias:

  \code{echo 'alias python=python2.7' >> ~/.bashrc}

  Mandatory and optional dependencies of \pkg{P2C2M} can be installed 
  automatically via two installation scripts that are co-supplied with 
  the package. These scripts were designed for unix-like operating 
  systems and are located in folder /exec. To use these installation 
  scripts, a correct configuration of python2-setuptools is required. 
  Users of unix-like operating systems can insure a correct 
  configuration by setting the following alias:

  \code{echo 'alias python-config=python2-config' >> ~/.bashrc}

  To execute the R installer, please run the following commands in R:

  \code{source('/path_to_P2C2M/exec/P2C2M.installRlibs.R')}; 
  \code{p2c2m.install()}
  
  To execute the Python installer, please run the following command in 
  a terminal:

  \code{python /path_to_P2C2M/exec/P2C2M.installPylibs.py}


  \emph{Special Note for MacOS}

  Users of the MacOS operating system need to install the dependencies 
  manually. Prior to their installation, please confirm that file 
  '/usr/bin/python2-config' exists in your file system and that it 
  points to the Python 2.7 executable. Please refer to \url{http://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html} 
  on how to install R packages manually. For the manual installation of 
  Python libraries, please refer to 
  \url{http://docs.python.org/2/using/mac.html}


  \bold{Study Design Requirements}

  In the user-supplied data set, every species should be represented by 
  at least two alleles. Species that are represented by only a single 
  allele, by contrast, must be specified via option 
  \code{"single.allele"} and thereby are not included in the 
  calculation of the summary statistic 'GSI'; misspecifications 
  causes \pkg{P2C2M} to print the error message \emph{'Error: given 
  group represents one or fewer taxa. Cannot compute index.'}).


  \bold{Input File Requirements}

  In order to execute \pkg{P2C2M}, a user must provide a directory with 
  three different types of input files: (a) a file that contains 
  species trees, (b) a file that contains gene trees for each gene 
  under study, and (c) an XML-formatted file generated by BEAUTi, the 
  input generator of BEAST (Drummond and Rambaut 2007). A species tree 
  file contains a draw of s generations from the posterior distribution 
  of species trees. Each gene tree file contains an equally large draw 
  from the respective posterior distribution of ultrametric 
  genealogies. Please note that the generations recorded in the species 
  tree file must match those in the gene tree files exactly. The input 
  file generated by BEAUTi is formatted in XML markup language and 
  represents the starting point for a species tree inference in *BEAST. 
  Here, it provides information on allele and species names, the 
  association between alleles and species, and ploidy levels to 
  \pkg{P2C2M}.


  \bold{File Name Requirements}

  The following requirements for input file names are in place: The 
  species tree file must be named 'species.trees'. Each gene tree file 
  must be named 'g.trees', where the letter g is substituted with the 
  actual name of the gene. The name of the xml-formatted input file is 
  not constrained and at the discretion of the user. Please be aware 
  that \pkg{P2C2M} uses the name of the xml-formatted input file name 
  to label all subsequent output of the package.
}

\author{
  Michael Gruenstaeudl, Noah Reid

  Maintainer: Michael Gruenstaeudl \email{mi.gruenstaeudl@gmail.com}
}

\references{
  Drummond, A.J. and Rambaut, A. (2007) BEAST: Bayesian evolutionary analysis by sampling trees. \emph{BMC Evolutionary Biology}, \bold{7}, 214.

  Gruenstaeudl, M., Reid, N.M., Wheeler, G.R. and Carstens, B.C., submitted. Posterior Predictive Checks of Coalescent Models: P2C2M, an R package.

  Heled, J. and Drummond, A.J. (2010) Bayesian inference of species trees from multilocus data. \emph{Molecular Biology And Evolution}, \bold{27}, 570--580.

  Reid, N.M., Brown, J.M., Satler, J.D., Pelletier, T.A., McVay, J.D., Hird, S.M. and Carstens, B.C. (2014) Poor fit to the multi-species coalescent model is widely detectable in empirical data. \emph{Systematic Biology}, \bold{63}, 322--333.
}

