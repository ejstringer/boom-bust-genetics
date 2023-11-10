
# source ------------------------------------------------------------------
source("./code/libraries.R")        # load libraries
source("./code/functions_create.R") # load functions
source('./code/functions_neEstimate.R') # prep for linux analysis
# load --------------------------------------------------------------------

ph <- readRDS("./data/pherm_genotypes.rds")        # DArT genotype data
sy <- readRDS("./data/ALL_syoung_genotypes.rds")   # DArT genotype data
rain_caps <- read.csv("./data/ecological_data.csv")  # rainfall and captures
grids <- read.csv("./data/em_gridcoortable.csv")   # grid coordinates

rain_caps$trip <- ymd(rain_caps$trip)

# grid distance  ----------------------------------------------------------
pairwise_grids <- em.distance.grids(grids) # distance between grids


# filtering ---------------------------------------------------------------
pherm <- em.filtering(ph, maf.cutoff = 0.005) # sandy inland mouse
syoung <- em.filtering(sy, maf.cutoff = 0.01) # lesser hairy-footed dunnart


# heterozygosity ----------------------------------------------------------
pherm_heterozygosity <- em.ind.ho(pherm, rain_caps)     # individual Ho
syoung_heterozygosity <- em.ind.ho(syoung, rain_caps)   #   calculations 

# Fst ---------------------------------------------------------------------
pherm_fst <- em.fst(pherm,             # if there is a saved fst file it will 
                    pairwise_grids,    # load the file instead of recalculating
                    rain_caps,         # fst again. 
                    min.sample = 4,    # Minimum individuals per grid per trip
                    filecode = 'pherm')


# Ne by trip --------------------------------------------------------------

prepNePh <- em.genlight.ne.prep(pherm, spp = 'ph') 
prepNeSy <- em.genlight.ne.prep(syoung, spp = 'sy')

saveRDS(prepNePh, './data_processed/pherm_ne_prep.rds')   # for neEstimator analysis
saveRDS(prepNePh, './data_processed/syoung_ne_prep.rds')  # for neEstimator analysis
saveRDS(rain_caps, './data_processed/rain_captures.rds') # for neEstimator analysis
 

