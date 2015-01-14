library(RUnit)
library(AmpCNVrstudio)

source('./R/FolderAnalysis_functions.R')
source('./inst/unitTests/test_testing.R')
test.bootList()
#test.examples()











#library('RUnit')
#source('sample.R')
#test.suite <- defineTestSuite("AmpCNVrstudio",                               dirs = file.path("tests"),                               testFileRegexp = '^\\d+\\.R') 
#test.result <- runTestSuite(test.suite)
#printTextProtocol(test.result)
