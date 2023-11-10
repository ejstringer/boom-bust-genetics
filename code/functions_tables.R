
em.model_table <- function(listmodels = estimates, decPts = 3){
  flex <- list()
  
  for(i in 1:length(listmodels)){
    m <- listmodels[[i]]
    m.name <- paste("period", parse_number(names(listmodels)[i]))
    row.names(m) <- paste0('b', 1:nrow(m)) #gsub("\\(|\\)", "", rownames(m))
    
    est <- paste0(round(m$Estimate, decPts),
                  "  (", paste(round(m$lower,decPts), "-", round(m$upper,decPts)), ")")
    r2 <- round(m$conditionalR2[nrow(m)], decPts)
    response <- m$response[nrow(m)]
    significance <- paste(m$Pr...t.. < 0.05, collapse = "-")
    
    m.species <- sub("\\_.*", "", names(listmodels)[i])
    m.type <- str_sub(names(listmodels)[i], - 2, - 1)
    
    str_extract( names(listmodels)[i], "(?<=_)[^_]")
    
    sub('\\_.*', '\1',         # Extract last word
        names(listmodels)[i])
    
    m.sum <- c(m.type,m.species, response, m.name, est, r2, significance)
    names(m.sum) <- c('type', 'species', "response", "period",
                      rownames(m), "R2", "significant")
    flex[[i]] <- data.frame(t(data.frame(m = m.sum)))
  }
  # flextable
  df <- reduce(flex, bind_rows) %>% 
    arrange(desc(response), species,
            desc(type), period) %>% 
    dplyr::select(-type) %>% 
    relocate(R2, significant, .after = last_col())  
  
  df$significant <- sub("TRUE", "FALSE", df$significant)
  
  return(df)
  
}


em.flextable_table <- function(df_ho, df_fst, df_ne){
  df_ho$species[duplicated(df_ho$species)] <- ''
  df_fst$species[duplicated(df_fst$species)] <- ''
  
  headfst <- "C) Genetic differentiation: log(Fst)"
  headne <- 'A) Effective population size estimates: log(Ne)'
  headho <- 'B) Genetic diversity: log(H)'
  neNob3 <-  names(df_ho)
  neNob3[neNob3=='b3'] <- ''
  
  
  df <- rbind(neNob3,
              df_ne,
              c(headho, rep('', ncol(df_ho)-1)),
              names(df_ho),
              df_ho,
              c(headfst, rep('', ncol(df_ho)-1)), 
              names(df_ho),
              df_fst
  ) %>% 
    dplyr::select(-response) %>% 
    rowwise() %>% 
    mutate(species = ifelse(species == 'Ph', 'P. hermanns', species),
           species = ifelse(species == 'Sy', 'S.youngsoni', species),
           significant = ifelse(significant == 'significant', NA,
                                significant)) %>% 
    rename(`A) Effective population size estimates: log(Ne)` = species)  # head ho
  
  df
  
  
  
  boldtb <- strsplit(df$significant, "-") %>% 
    lapply(., function(x) x[1:max(sapply(., length))]) %>% 
    do.call("rbind", .)
  
  
  #df$response[duplicated(df$response)] <- ""
  ft <- flextable(df[,-grep("significant", names(df))])
  ft <- width(ft, width = 2.2)
  
  ft <- color(ft, part = "footer", color = "#666666")
  
  coladjust <- grep('b1', names(df))-1
  for(r in 1:nrow(boldtb)){
    for (v in 1:ncol(boldtb)) {
      ft <- bold(ft, i = r, j = v+coladjust, bold = boldtb[r,v])
    }
  }
  ft <- set_table_properties(ft, layout = "autofit")
  
  ft_composition <- ft %>% 
    border_remove()%>% 
    hline(i = 11,
          border = fp_border(color = "grey40",
                             width = 2)) %>% 
    hline(i = c(3,17,22), 
          border = fp_border(color = "grey40", width = 4)) %>% 
    hline(i = c(1, 5, 19), 
          border = fp_border(color = "grey40", width = 2)) %>% 
    merge_at(i = 4 , j = 1:5, part = "body") %>% 
    merge_at(j = 1:7, part = "header") %>% 
    merge_at(i = 18 , j = 1:5, part = "body") %>%
    italic(i = c(2,3,4, 6, 12,18,20), j = 1) %>% 
    italic(part = "header") %>% 
    hline(i = c(8,14), 
          border = fp_border(color = "grey40", 
                             width = 1.5, style = "dotted")) %>% 
    
    flextable::align(i = 3, j = 1, align = "right", part = "body")
  
  ft_composition
  
  return(ft_composition)
}



em.model_predict <- function(estimates){
  
  e <-  estimates
  if(sum(grepl("metres", rownames(e)))>0) e <- e[-grep("metres", rownames(e)),]
  which(rownames(e) == "(Intercept)") 
  x <- seq(0,6, 0.01)
  xdata <- data.frame(cbind(int = rep(1, length(x)), x))
  if("yearsSince2" %in% rownames(e)) xdata$x2 <- xdata$x^2
  
  evalue  <-e$Estimate  
  up <- e$upper
  low <- e$lower
  
  y <- apply(xdata, MARGIN = 1, function(i) sum(i*evalue))
  uppred <- apply(xdata, MARGIN = 1, function(i) sum(i*up))
  lowpred <- apply(xdata, MARGIN = 1, function(i) sum(i*low))
  
  signrow <- ifelse("yearsSince2" %in% rownames(e),
                    "yearsSince2", "yearsSince")
  
  df <- data.frame(x = x, y = y, lower = lowpred, upper = uppred,
                   sig = e[grep(signrow, rownames(e)),]$Pr...t.. < 0.05)
  
  return(df)
  
}

em.ne.model.extract <- function(m, sp = 'Ph', p = '2006-2017'){
  modelNE <- summary((m))
  
  R2 <- modelNE$r.squared
  s <- modelNE
  co <- data.frame(s$coefficients)
  
  t95 <-  qt(0.975, df = modelNE$df[2])
  cil <- co$Estimate - co$Std..Error * t95
  ciu <- co$Estimate + co$Std..Error * t95
  
  co$lower <- cil
  co$upper <- ciu
  co$marginalR2 <- R2[1]
  co$conditionalR2[nrow(co)] <- R2[1]
  co$response <- as.character(formula(s))[2]
  
  df_neSy <- em.model_table(list(Ph = co, q=co))[1,] %>% 
    mutate(b3 = '', period = p) %>% 
    relocate(b3, .after = b2)
  df_neSy$species <- sp
  
  return(df_neSy)
}