# create functions -------------------------------------------------------------

em.rain_caps_period <- function(session, boom_dates){
  
  meanTripSummary <- session %>% 
    group_by(trip) %>%
    summarise(rain = mean(rain),
              captures = mean(phAbundance, na.rm = T),
              capturesSy = mean(syAbundance, na.rm = T)) %>%
    ungroup() %>% 
    mutate(year = lubridate::year(lubridate::ymd(trip)),
           trip = ymd(trip)) %>% 
    bind_rows(data.frame(trip = max(.$trip)+months(1:12)))

  ## period data
  rain <- meanTripSummary$rain
  
  trip <- meanTripSummary$trip 
  period <- c(1, which(trip %in% boom_dates), length(trip))
  periodLength <- length(period)
  
  periodName <- data.frame(row1 = period[-periodLength], 
                           row2 = c(period[2:(periodLength-1)]-1, 
                                    period[periodLength]),
                           periodName = paste("period", 0:(periodLength-2))) %>% 
    mutate(months = row2-row1 + 1) %>% 
    split(seq(nrow(.)))
  
  periodSeq <- reduce(lapply(periodName,
                             function(x) rep(x$periodName,x$months)), c)
  periodSince <- reduce(lapply(periodName, function(x) 0:(x$months-1)), c)
  
  rain_caps <- meanTripSummary %>% 
    mutate(sinceEvent = periodSeq,
           monthsSince = periodSince,
           yearsSince = monthsSince/12)
  
  return(rain_caps)
}

em.distance.grids <- function (dfcoor, between = "gridId", proj = "+proj=longlat", 
                            transform = "+init=epsg:32754") {
  cord.dec = sp::SpatialPoints(cbind(dfcoor$lon, dfcoor$lat), 
                               proj4string = sp::CRS(proj))
  cord.UTM <- sp::spTransform(cord.dec, sp::CRS(transform))
  mat <- as.matrix(dist(cord.UTM@coords))
  matNames <- dfcoor[, between]
  rownames(mat) <- matNames
  colnames(mat) <- matNames
  distmat <- as.dist(mat)
  dfdist <- otuSummary::matrixConvert(distmat, colname = c("pop1", 
                                                           "pop2", "metres"))
  dfdist_pairs <- dfdist %>%
    rowwise() %>% 
    mutate(pairs = paste(sort(c(pop1, pop2)), collapse = "-"),
           km = metres/1000) %>% 
    dplyr::select(pairs, metres, km)
  
  return(dfdist_pairs)
}

em.filtering <- function (glx, maf.cutoff) 
{
  
  glO <- glx
  glx <- gl.filter.callrate(glx, method = 'ind', threshold = 0.5)
  glx <- gl.filter.callrate(glx, method = "loc", threshold = 0.95, verbose = 0)
  glx <- gl.filter.callrate(glx, method = "ind", threshold = 0.95, verbose = 0)
  glx <- gl.filter.monomorphs(glx, verbose = 0)
  glx <- gl.filter.maf(glx, threshold = maf.cutoff, verbose = 0)
  glx <- gl.filter.secondaries(glx, verbose = 0, method = 'best')
  glx <- gl.filter.rdepth(glx, verbose =0)
  glx <- gl.filter.reproducibility(glx, threshold = 0.95)
  glx <- em.remove.relatives(glx, threshold = 0.99)
  
  cat("\n")
  cat("\n")
  cat("pre filtering loc:", nLoc(glO), "and ind:", 
      nInd(glO), "\n")
  cat("post filtering loc:", nLoc(glx), "and ind:", 
      nInd(glx), "\n")
  return(glx)
}

em.ind.ho <- function (gl, boomData) 
{
  
  pop(gl) <- gl@ind.names
  ind_as_pop <- seppop(gl)
  
  mean.Ho <- function(x)  mean(gl.Ho(x), na.rm = TRUE)
  heterozygosity_Inds <- sapply(ind_as_pop, mean.Ho)
  
  df <- data.frame(id = names(heterozygosity_Inds), 
                   ho = heterozygosity_Inds, 
                   row.names = NULL)
  
  heterozygosity_metadata <- left_join(df, 
                                   gl@other$ind.metrics)
  heterozygosity_metadata$ho.log <- log(heterozygosity_metadata$ho)
  

  
  final_het <- heterozygosity_metadata %>% 
    dplyr::select(species, id, pop, gridId, trip, lat, lon, sex, ho, ho.log) %>%
    mutate(trip = ymd(trip),
           species =   gsub('\\b(\\pL)\\pL{4,}|.','\\1', 
                            species,perl = TRUE)) %>% 
    left_join(boomData) %>% 
  filter(sinceEvent != "period 0") %>%
    mutate(yearsSince2 = yearsSince^2) %>% 
    as.data.frame() %>% 
    split.data.frame(., .$sinceEvent)
  # save results as csv
  save_het <- do.call('rbind', final_het)
  het_filename <- paste0(sub(' ', '_', gl@other$ind.metrics$species[1]), 
                         "_ind_Ho.csv")
  
  write.csv(save_het, paste0("./data_processed/", het_filename), row.names = F)
  
  return(final_het)
}


em.fst <- function(gl, pairwiseDistance, boomData, min.sample = 4, filecode = 'pherm'){
  
  fst_filename <- paste0(filecode, "_fst_min", min.sample, ".csv")
  
  if(fst_filename %in% list.files("./data_processed/")){
    print("loading from data_processed folder...")
    save_fst <- read.csv(paste0("./data_processed/", fst_filename))
  }else{
  # min sample 
  gl@other$ind.metrics$gridtrip <- paste(gl@other$ind.metrics$gridId,
                                         gl@other$ind.metrics$trip)
  
  meta <- gl@other$ind.metrics %>% 
    mutate(poptrip = paste(gridId, trip),
           n = 1) 
  
  grid_trips_include <- data.frame(n = tapply(meta$n, meta$poptrip, sum)) %>%  
    filter(ifelse(n >= min.sample, TRUE, FALSE)) %>% 
    mutate(pop = rownames(.)) %>%
    separate(pop, into = c("gridId", "trip"), sep = " ", remove = FALSE) %>% 
    filter(trip %in% trip[duplicated(trip)]) %>% 
    dplyr::select(pop)
  
  
  fst_subset <- gl.keep.pop(gl, 
                            pop.list = grid_trips_include$pop,
                            as.pop = "gridtrip")
  
  pop(fst_subset) <- fst_subset@other$ind.metrics$trip
  pop_sep <- seppop(fst_subset)
  
  em.assign.pop <- function(x) {
    pop(x) <- x@other$ind.metrics$gridtrip
    return(x)
  }
  
  pre_fst <- lapply(pop_sep, em.assign.pop)
  
  # fst 
  fst <- pbapply::pblapply(pre_fst,
                           function(x) gl.fst.pop(x, nboots = 1))
  
  
  fst.dataframe <- function(x) {
    otuSummary::matrixConvert(as.dist(x),
                              colname = c("pop1", "pop2", "fst")) %>% 
      separate(pop1, into = c("pop1", "trip"), sep = " ") %>%
      separate(pop2, into = c("pop2", "trip2"), sep = " ") %>% 
      mutate(trip = ymd(trip)) %>% 
      rowwise() %>% 
      mutate(pairs = paste(sort(c(pop1, pop2)), collapse = "-"),
             fst.log = log(fst + abs(min(fst))*1.1)) %>% 
      dplyr::select(trip, pairs, pop1, pop2, fst)
    
  }
  
  clean_fst <- lapply(fst, fst.dataframe)
  
  all_fst <- reduce(clean_fst, bind_rows)
  
  save_fst <- all_fst %>% 
    left_join(boomData) %>% 
    left_join(pairwiseDistance) 
  
  write.csv(save_fst, paste0("./data_processed/", fst_filename), row.names = F)
  
  }
  
  load_fst <- save_fst %>% 
    filter(sinceEvent != "period 0") %>%
    mutate(yearsSince2 = yearsSince^2,
           trip = ymd(trip),
           fst.log = log(fst + abs(min(fst))*1.1),
           species = gl@other$ind.metrics$species[1],
           species =   gsub('\\b(\\pL)\\pL{4,}|.','\\1', 
                            species,perl = TRUE)) %>% 
    as.data.frame() %>% 
    split.data.frame(., .$sinceEvent)
  
  return(load_fst)
}


em.remove.relatives <- function(pherm, threshold = 0.25){
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
      
      while(sum(df2 >=threshold, na.rm = T) >0){
        
        nsibs <- rowSums(!is.na(as.matrix(df2)))
        topsib <-  which.max(nsibs)
        
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
  pherm2 <- pherm
  nRemove <- 0
  if(length(removeList)>0){
    
   names(removeList) <- triprelated
   removeALL <- do.call('c', removeList)
   pherm2 <- gl.drop.ind(pherm, ind.list = removeALL)
   nRemove <- length(removeALL)
  }
  
  pherm2@other$ind.metrics$id <- sub('x','', pherm2@other$ind.metrics$id)
  pherm2@ind.names<- sub('x','', pherm2@ind.names)
  
  cat("number of relatives removed:", nRemove, "\n", 
      nInd(pherm2), "individuals remaining", '\n')
  return(pherm2)
  
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

em.mean.CI<- function(x, stat) {
  y <- do.call('rbind', x)[,stat] %>% mean
  s <- do.call('rbind', x)[,stat] %>% sd
  n <- do.call('rbind', x)[,stat] %>% length
  
  se <- s/sqrt(n)
  tstar <- qt(p = 0.975, df =n-1)
  c(mean = y, lower = y-se*tstar, upper = y+se*tstar)
}
