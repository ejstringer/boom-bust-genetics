# functions
em.genlight.ne.prep <- function(gl, spp = 'sy'){
  
  glsub <- em.remove.perc.relatives(gl, threshold = 0.25,
                                    removePercent = 0.75)
  pop(glsub) <- glsub@other$ind.metrics$trip
  table(glsub@pop)
  
  popsy <- ymd(pop(glsub))
  popsyPlain <- popsy
  if(spp == 'sy'){
    datechanges <- data.frame(defined = c('2006 Feb to May',
                                          '2006 Jun to Nov',
                                          '2008 Feb to Oct',
                                          '2011 Jun to Nov',
                                          '2009 Aug to Sep',
                                          '2013 Apr to Feb 2014',
                                          '2017 Sep to Jun 2018'),
                              trip = ymd(c('2006-02-01',
                                           '2006-08-01',
                                           '2008-06-01',
                                           '2011-09-01',
                                           '2009-09-01',
                                           '2013-09-01',
                                           '2018-04-01')))
    
    popsyPlain[popsy > ymd('2006-01-01') & popsy < ymd('2006-06-01')] <- ymd('2006-02-01')
    popsyPlain[popsy > ymd('2006-05-01') & popsy < ymd('2006-12-01')] <- ymd('2006-08-01')
    popsyPlain[year(popsy) == '2008'] <- ymd('2008-06-01')
    popsyPlain[year(popsy) == '2011'] <- ymd('2011-09-01')
    popsyPlain[popsy > '2017-06-01'] <- ymd('2018-04-01')
    popsyPlain[popsy > ymd('2009-07-01') & popsy < ymd('2009-10-01')] <- ymd('2009-09-01')
    popsyPlain[popsy > ymd('2013-03-01') & popsy < ymd('2014-03-01')] <- ymd('2013-09-01')
  }
  if(spp == 'ph'){
    datechanges <- data.frame(defined = c('2006 Jun to Nov',
                                          '2008 Jul to Apr 2009',
                                          '2012 April to Nov',
                                          '2014 Feb to Nov'),
                              trip = ymd(c('2006-09-01',
                                           '2008-09-01',
                                           '2012-04-01', 
                                           '2014-06-01')))
    
    popsyPlain[popsy > ymd('2006-05-01') & popsy < ymd('2006-12-01')] <- ymd('2006-09-01')
    popsyPlain[popsy > ymd('2008-06-01') & popsy < ymd('2009-05-01')] <- ymd('2008-09-01')
    popsyPlain[year(popsy) == '2012'] <- ymd('2012-04-01')
    popsyPlain[year(popsy) == '2014'] <- ymd('2014-06-01')
    
  }
  pop(glsub) <- popsyPlain
  
  # filtering 
  ntrip <- table(glsub@pop)
  ntrip20 <-ntrip[ntrip>19 & ymd(names(ntrip)) > ymd('2006-01-01')]
  
  glsub20 <- gl.keep.pop(glsub, 
                         pop.list = names(ntrip20))
  glsub20 <- gl.recalc.metrics(glsub20)
  glsub20 <- gl.filter.monomorphs(glsub20)
  
  
  glpop <- seppop(glsub20)
  
  minAf <- 4/(sapply(glpop, nInd)*2)
  popMaf <- list() 
  for(i in 1:length(glpop)){
    glx <- glpop[[i]]
    
    glxx <- gl.filter.maf(glx, threshold = minAf[i])
    glxx <- gl.filter.callrate(glxx, threshold = 0.99)
    
    popMaf[[i]] <- glxx
    
    
  }
  names(popMaf) <- names(minAf)
  return(popMaf)
  
}

em.ne.summarise <- function(neEstimates, glprep, rain_caps, spp = 'sy'){
  
  if(spp == 'sy'){
    datechanges <- data.frame(defined = c('2006 Feb to May',
                                          '2006 Jun to Nov',
                                          '2008 Feb to Oct',
                                          '2011 Jun to Nov',
                                          '2009 Aug to Sep',
                                          '2013 Apr to Feb 2014',
                                          '2017 Sep to Jun 2018'),
                              trip = ymd(c('2006-02-01',
                                           '2006-08-01',
                                           '2008-06-01',
                                           '2011-09-01',
                                           '2009-09-01',
                                           '2013-09-01',
                                           '2018-04-01')),
                              from = ymd(c('2006-02-01',
                                           '2006-06-01',
                                           '2008-02-01',
                                           '2011-06-01',
                                           '2009-08-01',
                                           '2013-04-01',
                                           '2017-09-01')),
                              to = ymd(c('2006-05-01',
                                         '2006-11-01',
                                         '2008-10-01',
                                         '2011-11-01',
                                         '2009-09-01',
                                         '2014-02-01',
                                         '2018-06-01')))
    
    datechanges$meanCaps <- NA
    for(i in 1:nrow(datechanges)){
      caps <- filter(rain_caps, 
                     trip >= datechanges$from[i] & trip <= datechanges$to[i]) 
      datechanges$meanCaps[i] <- mean(caps$capturesSy, na.rm = T)
      
    }
  }
  if(spp == 'ph'){
    datechanges <- data.frame(defined = c('2006 Jun to Nov',
                                          '2008 Jul to Apr 2009',
                                          '2012 April to Nov',
                                          '2014 Feb to Nov'),
                              trip = ymd(c('2006-09-01',
                                           '2008-09-01',
                                           '2012-04-01', 
                                           '2014-06-01')),
                              from = ymd(c('2006-06-01',
                                           '2008-07-01',
                                           '2012-04-01',
                                           '2014-02-01')),
                              to = ymd(c('2006-11-01',
                                         '2009-04-01',
                                         '2012-11-01',
                                         '2014-11-01')))
    datechanges$meanCaps <- NA
    for(i in 1:nrow(datechanges)){
      caps <- filter(rain_caps, trip >= datechanges$from[i] & trip <= datechanges$to[i]) 
      datechanges$meanCaps[i] <- mean(caps$captures, na.rm = T)
    }
  }
  
  
  
  ne <- as.numeric(sapply(neEstimates, function(x) x[[1]][6,2]))
  lower <- as.numeric(sapply(neEstimates, function(x) x[[1]][7,2]))
  upper <- as.numeric(sapply(neEstimates, function(x) x[[1]][8,2]))
  trip <- ymd(sapply(neEstimates, function(x) names(x)))
  
  n <- sapply(glprep, nInd)
  loci <- sapply(glprep, nLoc)
  
  df <- data.frame(trip, ne,lower, upper, n, loci) %>% 
    left_join(rain_caps, by = 'trip') %>% 
    left_join(datechanges) %>% 
    rowwise() %>% 
    mutate(defined = ifelse(is.na(defined), as.character(trip), defined),
           trip = ymd(trip),
           captures = ifelse(spp == 'sy', capturesSy, captures),
           ciInfinite = ifelse(is.infinite(lower) | is.infinite(upper),
                               'inf', 'no'),
           ne.log = log(ne),
           captures.log = log(captures),
           species = spp) %>% 
    dplyr::select(-capturesSy, -year, -monthsSince)
  return(df)
}

em.remove.perc.relatives <- function(pherm, threshold = 0.25, removePercent = 0.75){
  #pherm <- readRDS('./output/pherm_filtered.rds')
  
  
  pherm@other$ind.metrics$id <- paste0('x', pherm@other$ind.metrics$id)
  pherm@ind.names <- paste0('x', pherm@ind.names)
  pherm@other$ind.metrics$id == pherm@ind.names
  
  pop(pherm) <- pherm@other$ind.metrics$trip
  phpop <- seppop(pherm)
  
  
  # find related inidividuals
  # setup 
  removeList <- list()
  triprelated <- NA
  
  
  r <- 1 #remove list
  
  for(y in 1:length(phpop)){
    # choose trip
    tripx <- phpop[[y]]
    related <- FALSE
    # relativeness
    if(nInd(tripx)>1){ 
      ibdecent <- gl.grm(tripx,plotheatmap = F) 
      related <- sum(as.dist(ibdecent) >= threshold) != 0
      
    }
    # if related add to remove list
    if(related){
      df <- data.frame(ibdecent)
      for(k in 1:nrow(df)) df[k,k] <- NA
      df[(df) < threshold] <- NA
      df[1:6, 1:5]
      
      
      indRel <- apply(df, MARGIN = 1, function(x) sum(x, na.rm = T))
      
      sibs <- names(indRel[indRel>0])
      
      removeInd <- NA
      x = 1
      df2 <- df[sibs, sibs]
      nKeep <- length(sibs)- floor(length(sibs)*removePercent)
      while(sum(df2 >=threshold, na.rm = T) > nKeep){
        
        nsibs <- rowSums(!is.na(as.matrix(df2))) #remove based on most siblings
        
        topsib <- which.max(nsibs)
        
        print(x)
        removeInd[x] <- names(topsib)
        x = x+1
        df2 <- df2[-topsib, -topsib]
        
      }
      
      removeList[[r]] <- removeInd
      triprelated[r] <- names(phpop)[y]
      
      r = r+1
    }
    
  }
  
  names(removeList) <- triprelated
  removeALL <- do.call('c', removeList)
  
  
  
  df2<- data.frame(indn = sapply(phpop, nInd)) %>%
    mutate(trip = rownames(.)) %>%
    left_join(data.frame(rel = sapply(removeList, length),
                         trip = triprelated)) %>%
    mutate(prop = round(rel/indn,2),
           remove75 = floor(rel*removePercent),
           reprop = round(remove75/rel,2))
  propR <- df2$remove75
  names(propR) <- df2$trip
  print(propR)
  
  
  pherm2 <- gl.drop.ind(pherm, ind.list = removeALL)
  
  pherm2@other$ind.metrics$id <- sub('x','', pherm2@other$ind.metrics$id)
  pherm2@ind.names<- sub('x','', pherm2@ind.names)
  cat("number of relatives removed:", length(removeALL), "\n", 
      nInd(pherm2), "individuals remaining", '\n')
  return(pherm2)
  
}
