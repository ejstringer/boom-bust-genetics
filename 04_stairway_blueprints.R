

# source ------------------------------------------------------------------

source('./code/functions_stairway.R')

# load --------------------------------------------------------------------

# run 01_curate.R for ph and sy objects and functions.

# define samples ----------------------------------------------------------
table(ph@other$ind.metrics$trip)
table(sy@other$ind.metrics$trip)

phTrips <- list(phBoom07 = '2007-09-01',
                phBust06 = '2006-05-01',
                phBoom10 = '2010-11-01',
                phBust09 = c('2009-08-01', '2009-09-01'))
syTrips <- list(syRain07 = '2007-09-01',
                syDry06 = c('2006-05-01','2006-06-01','2006-08-01'),
                syRain10 = c('2010-06-01', '2010-11-01'),
                syDry09 = c('2009-04-01', '2009-08-01', '2009-09-01'))

phsub <- lapply(phTrips, function(x) gl.keep.pop(ph, x, 'trip'))
sysub <- lapply(syTrips, function(x) gl.keep.pop(sy, x, 'trip'))

# filter ------------------------------------------------------------------
nedata <- lapply(c(phsub, sysub), em.filtering.stairway) 
sapply(nedata, nInd)

# downsample --------------------------------------------------------------

down <- lapply(nedata, em.downsample, ndown = 32, seed = 444)
down2 <- lapply(nedata, em.downsample, ndown = 32, seed = 555)

names(down) <- paste0('down_', names(nedata))
names(down2) <- paste0('down2_', names(nedata))

sapply(down, nInd)
sapply(down2, nInd)

# selection ---------------------------------------------------------------
nedataUse <- c(nedata, down[1:2], down2[1:2]) # only down sampling ph06 and ph07
names(nedataUse)

# blueprints --------------------------------------------------------------
dir.create('./data_processed/blueprints/') # create director if needed
for (i in 1:length(nedataUse)) {
  print(i)
  ne <- nedataUse[[i]]
  folder <- names(nedataUse)[i]
  filename <- paste0(folder, '.blueprint')

  loci <- nLoc(ne)*69
  mutation_rate <- ifelse(grepl('sy', folder), 
                          1.394808e-08, # sy: https://doi.org/10.1093/molbev/msz191
                          0.39e-8 # ph https://doi.org/10.1038/s41467-019-12023-w
                          )
                          
                        
  gl2stairwayplot2(ne, 
                   stairway_plot_dir = 'stairway_plot_es',   # this is relevant when running the blueprints using the bash script (putty_stairway_code)
                   outpath = './data_processed/blueprints/', # directory for blueprints to be saved
                   verbose = T, 
                   mu = mutation_rate,
                   gentime = 1,         # years
                   simfolder = folder,
                   outfile= filename, 
                   run=FALSE, 
                   nreps = 200,
                   L=loci)
  
}

# stairway blue prints run using bash on linux server 
#     (next step: file 05_putty_stairway_code)

