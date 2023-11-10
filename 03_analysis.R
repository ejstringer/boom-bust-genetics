


# source ------------------------------------------------------------------

source('./code/functions_models.R')


# load data ---------------------------------------------------------------

# get data by running 01_curate.R and 02_neEstimator_linux.R then 
# downloading the data created by 02_neEstimator_linux.R and putting
# them into the data_processed folder

phne <- read.csv('./data_processed/pherm_ne_estimates_summary.csv')
syne <- read.csv('./data_processed/syoung_ne_estimates_summary.csv') %>% 
  filter(!is.infinite(lower), !is.infinite(upper))

pherm_fst <- em.fst(pherm, pairwise_grids, rain_caps,        
                    min.sample = 4, filecode = 'pherm')


# data summary ------------------------------------------------------------

lapply(split(rain_caps, as.character(rain_caps$sinceEvent)), 
       function(x) c(min(ymd(x$trip)), max(ymd(x$trip))))   # period time frames


pherm %>% nInd()    # sampled individuals
syoung %>% nInd()   # sampled individuals


sapply(pherm_heterozygosity, nrow)    # mice samples per period
sapply(syoung_heterozygosity, nrow)   # dunnart samples per period

sapply(pherm_fst, nrow)               # grid pairwise comparisons per period

nlevels(pherm_heterozygosity$`period 1`$gridId) # number of grids

# mean sd ho
em.mean.CI(pherm_heterozygosity,'ho')
em.mean.CI(syoung_heterozygosity,'ho')
em.mean.CI(pherm_fst, 'fst')

em.mean.CI(list(phne), 'ne')
em.mean.CI(list(syne), 'ne')

psych::harmonic.mean(phne$ne)
psych::harmonic.mean(syne$ne)

# models -----------------------------------------------------------------

genetic_change <- c(pherm_heterozygosity, 
                    syoung_heterozygosity,
                    pherm_fst)


drift_models <- em.run_models(genetic_change)      # linear mixed effects models

sapply(drift_models,
       function(x) MuMIn::r.squaredGLMM(x)[2]) %>% 
  round(., 3)                                      # condition R squared


# model summaries ---------------------------------------------------------

summary(drift_models$Ph_ho_period1_lm)
summary(drift_models$Ph_ho_period2_quad)
summary(drift_models$Ph_ho_period3_lm)
summary(drift_models$Sy_ho_period1_quad)
summary(drift_models$Sy_ho_period2_lm)
summary(drift_models$Sy_ho_period3_lm)
summary(drift_models$Ph_fst_period1_lm)
summary(drift_models$Ph_fst_period2_lm)
summary(drift_models$Ph_fst_period3_lm)

# ne and capture ----------------------------------------------------------


ne_cap_models <- lapply(list(Ph_ne_caps_lm = phne, Sy_ne_caps_lm = syne), 
                        function(x) lm(ne.log ~ captures.log, data = x)) 

lapply(ne_cap_models, summary)


