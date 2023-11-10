# this requires NeEstimator V2.1 software
# https://www.molecularfisherieslaboratory.com.au/neestimator-software/

# files -------------------------------------------------------------------
# upload files from this project to linux 
# server to run parallel neEstimator analysis

   # pherm_ne_prep.rds 
   # syoung_ne_prep.rds
   # ecological_data.csv

   # functions_neEstimate.R


# libraries ---------------------------------------------------------------

library(doParallel)
library(dartR)
library(tidyverse)
library(lubridate)


# source ------------------------------------------------------------------

source('functions_neEstimate.R')

# load --------------------------------------------------------------------
prepNePh <- readRDS("pherm_ne_prep.rds")        # data from 01_curate.R 
prepNeSy <- readRDS("syoung_ne_prep.rds") 
rain_caps <- readRDS("rain_captures.rds")


# ind and loci ------------------------------------------------------------
sapply(prepNePh, nInd)
sapply(prepNePh, nLoc)

sapply(prepNeSy, nInd)
sapply(prepNeSy, nLoc)

# parallel ----------------------------------------------------------------


#Setup backend to use many processors
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makeCluster(21) 
registerDoParallel(cluster)

x <- foreach(i = 1:20) %dopar% {
  sqrt(i)
}
x

neEstimate <- foreach(i = 1:length(prepNePh)) %dopar% {
  library(dartR)
  gl.LDNe(prepNePh[[i]],
          neest.path = '/your_linux_root_directory/neEstimator/Zip_Folder_64_bit_191125/', 
          save2tmp = T, singleton.rm = T)
}


neEstimateSy <- foreach(i = 1:length(prepNeSy)) %dopar% {
  library(dartR)
  gl.LDNe(prepNeSy[[i]],
          neest.path = '/your_linux_root_directory/neEstimator/Zip_Folder_64_bit_191125/', 
          save2tmp = T, singleton.rm = T)
}

#Stop cluster
stopCluster(cluster)



# save ne estimates -------------------------------------------------------
saveRDS(neEstimatePh, 'pherm_neEstimates.rds')
saveRDS(neEstimateSy, 'syoung_neEstimates.rds')

# summarise ---------------------------------------------------------------
syne <- em.ne.summarise(neEstimateSy, prepNeSy, rain_caps, spp = 'sy') %>% 
  mutate(ogCaptures = captures,
         captures = ifelse(is.na(meanCaps), ogCaptures, meanCaps),
         captures.log = log(captures))


phne <- em.ne.summarise(neEstimatePh, prepNePh, rain_caps, spp = 'ph') %>% 
  mutate(ogCaptures = captures,
         captures = ifelse(is.na(meanCaps), ogCaptures, meanCaps),
         captures.log = log(captures))


## save summary ----------------------------------------------------------
write.csv(syne, '/home/stringer2/syoung_ne_estimates_summary.csv', row.names = F)
write.csv(phne, '/home/stringer2/pherm_ne_estimates_summary.csv', row.names = F)
