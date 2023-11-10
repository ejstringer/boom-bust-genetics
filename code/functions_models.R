
# model functions 
em.run_models <- function(gdata){
  
  models <- list()
  modelnames <- NA
  j <- 1
  for (i in 1:length(gdata)) {
    g <- gdata[[i]]
    
    
    if("ho" %in% names(g)){
      models[[j]] <- lmer(ho.log ~ yearsSince + (1|gridId), data = g)
      modelnames[j] <- paste(unique(g$species), "ho", sub(" ", "", unique(g$sinceEvent)),
                             "lm", sep = '_')
      
      j = j+1
      
      models[[j]] <- lmer(ho.log ~ yearsSince + yearsSince2 + (1|gridId), 
                          data = g)
      modelnames[j] <- paste(unique(g$species), "ho", sub(" ", "", unique(g$sinceEvent)),
                             "quad", sep = '_')
      
      j = j+1
      
    }# end if
    
    if('fst' %in% names(g)){
      models[[j]] <- lmer(fst.log ~ yearsSince +  log(metres) +
                            (1|pairs), data = g)
      modelnames[j] <- paste('Ph', "fst", sub(" ", "", unique(g$sinceEvent)),
                             "lm", sep = '_')
      j = j+1
    }#end if 
    
    
  }# end i
  
  names(models) <- modelnames
  return(models)
}

em.estimates_ci <- function(model)   {
  
  R2 <- MuMIn::r.squaredGLMM(model)
  s <- summary(model)
  co <- data.frame(s$coefficients)
  
  t95 <-  qt(0.975, df = co$df)
  cil <- co$Estimate - co$Std..Error * t95
  ciu <- co$Estimate + co$Std..Error * t95
  
  co$lower <- cil
  co$upper <- ciu
  co$marginalR2[nrow(co)] <- R2[1]
  co$conditionalR2[nrow(co)] <- R2[2]
  co$response <- as.character(formula(s))[2]
  return(co) 
}


