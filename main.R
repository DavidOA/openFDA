gCurrentDir = dirname(parent.frame(2)$ofile)
setwd(gCurrentDir)

source('importLibraries.R')
source('fetchData.R')
source('runStats.R')
source('plotReports.R')
