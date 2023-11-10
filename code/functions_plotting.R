
em.plot_level1 <- function(g_data, model_pred, model = 'Ph_ho_period1_lm'){
  
  filterby <- unlist(strsplit(model, '_'))

  
  dfpred <- model_pred[[model]]
  
  g_data2 <- g_data[sapply(sapply(g_data, names), 
                           function(x) filterby[2] %in% x)]
  
  g_data3 <- reduce(g_data2[sapply(g_data2,
                                   function(x) filterby[1] %in% x$species)],
                    bind_rows)
  
  
  d <- filter(g_data3, sinceEvent == paste('period',
                                           sub('period',"", filterby[3])))
  xName <- "yearsSince"
  

  dfpred2 <- filter(dfpred, x <= max(d$yearsSince))
  dfpred2$xVar <- dfpred2$x
  
  
  
  
  rGenetic <- ifelse(sum(grepl("fst", names(d))), "fst.log", "ho.log")
  d$response <- data.frame(d)[, rGenetic]
  d$xVar <- data.frame(d)[, xName]
  
  dmax <- max(d$xVar)
  gmin <- min(exp(g_data3[, rGenetic]))
  gmax <- max(exp(g_data3[, rGenetic]))
  
  dfpred2$y <- ifelse(dfpred2$y >= log(gmax), NA, dfpred2$y)
  dfpred2$lower <- ifelse(dfpred2$lower >= log(gmax), NA, dfpred2$lower)
  dfpred2$upper <- ifelse(dfpred2$upper >= log(gmax), log(gmax), dfpred2$upper)
  dfpred2$upper <- ifelse(is.na(dfpred2$lower), NA, dfpred2$upper)
  
  dfpred2$y <- ifelse(dfpred2$y <= log(gmin), NA, dfpred2$y)
  dfpred2$upper <- ifelse(dfpred2$upper <= log(gmin), NA, dfpred2$upper)
  dfpred2$lower <- ifelse(dfpred2$lower <= log(gmin), log(gmin), dfpred2$lower)
  dfpred2$lower <- ifelse(is.na(dfpred2$upper), NA, dfpred2$lower)
  
  l1 <- ifelse('fst.log' %in% names(d), 0, 4.5)
  r1 <- ifelse('fst.log' %in% names(d), 2, 2)
  p <- ggplot(dfpred2, aes(x = xVar, y = exp(y)))+
    geom_point(data = d, aes(x = xVar, y = exp(response)),
                size = 2, alpha = 0.2, width = 0.05)+
    # geom_ribbon(aes(ymin = exp(lower), ymax = exp(upper)), alpha = 0.1)+
    geom_line(data = filter(dfpred2, dfpred2$xVar <= dmax), 
              aes(x = xVar, y = exp(y)), 
              size=1, 
              linetype = ifelse(dfpred2$sig[1], "solid", "dashed"), 
              colour = "red", alpha = 0.9)+
    theme_bw()+
    scale_y_continuous(trans = "log", limits = c(gmin*0.99, gmax*1.01),
                       labels = scales::number_format(accuracy = ifelse('fst.log' %in% names(d), 0.001, 0.01)))+
    xlab(expression(italic("years since resource pulse (+3)"))) +
    ylab(expression(paste("", F[S*T], "")))+
    theme(axis.title.x = element_blank(),
          axis.ticks.length.y = unit(-1, "mm"),
          axis.title.y = element_blank(),
          plot.margin = margin(0.01, 0.1, 0.01, 0.01, "cm"),
          plot.background = element_rect(
            fill = "white"
          ),
          axis.text.y = element_text(margin = margin(l = l1, r = r1)),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  p
  
  
  return(p)
}


em.plot_level2 <- function(plots_level1, xticks = FALSE,
                           ylabFst = FALSE){
  
  if(length(plots_level1) != 3) stop('need three plots to combine')
    
  plot3 <- plots_level1
    
    p1fin <- plot3[[1]] +
      scale_x_continuous(limits = c(0, 3), breaks = c(0:4))
      
     
    p2fin <- plot3[[2]] +
      scale_x_continuous(limits = c(0, 5), breaks = c(0:6))+
      theme(axis.text.y = element_blank())
    
    p3fin <- plot3[[3]] +
      scale_x_continuous(limits = c(0, 2.5), breaks = c(0:4))+
      theme(axis.text.y = element_blank())
    
   

    ytitle <-grid::textGrob(expression(paste("Individual  ",
                                             H[` log-scale`])),
                            rot = 90,
                            gp = grid::gpar(fontsize = 12))   
    
    if(ylabFst){
      ytitle <- grid::textGrob(expression(paste(F[`ST  log-scale`])),
                               rot = 90,
                               gp = grid::gpar(fontsize = 12))
    }
    
      if(!xticks){
        p1fin <- p1fin +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        p2fin <- p2fin +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        p3fin <- p3fin +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        
  
    
      }

    
    gp <- gridExtra::grid.arrange(p1fin, p2fin, p3fin, nrow = 1,
                                  left = ytitle,
                                  widths = c(3, 5, 2.5))
    
 return(gp)
}


em.plot_labels <- function(x){ 
  ggplot() + 
    theme_void() +
    annotation_custom(grobTree(textGrob(paste("          Period", x))),
                      xmin = 0, xmax = 1,
                      ymin = 0, ymax = 1)
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

em.historic.Ne <- function(sy67, subsetspp. = 'ph', bw = TRUE){
  
  colcondition <- c(`dry period` = "white", 
                    `resource pulse` = "black")
  colcondition2 <- c(`dry period` = "black", 
                     `resource pulse` = "grey50")
  
  year1 <- c(`2006` = 'white',
             `2007` = 'grey70',
             `2009` = 'grey30',
             `2010` = 'black')
  
  year2 <- c(`2006` = 'grey50',
             `2007` = 'grey70',
             `2009` = 'grey30',
             `2010` = 'black')
  
  #myfillbooms <- c(`period 1` = '#FFC107',`period 2` = '#D81B60')
  
  myfillbooms <- c(`2006-2007` = 'black',`2009-2010` = 'grey50')
  
  myfillbooms2 <- c(`2006-2007` = 'black',
                    `2009-2010` = 'white')
  
  if (!bw) {
    colcondition <- colcondition2 <- c(`dry period` = "#E69F00", 
                                       `resource pulse` = "#009E73")
    
    year1 <- year2 <- c("#999999", "#E69F00", "#56B4E9",
                        "#009E73", "#F0E442", "#0072B2", 
                        "#D55E00", "#CC79A7")[1:4]
    
  }
  
  plotdata <- filter(sy67, year > 1, species == subsetspp.) %>% 
    mutate(nyear = str_sub(taken, 3,6),
           boom = ifelse(nyear %in% c('2006', '2009'), '(bust)', '(boom)'),
           taken2 = taken,
           yearboom = paste(nyear, boom))
  plotdata$nyear %>% table
  table(plotdata$nyear, plotdata$boom)
  
  plotdata$taken <- factor(plotdata$yearboom, levels = rev(sort(unique(plotdata$yearboom))))
  
  nehistoric<-ggplot(plotdata,
                     aes(year, Ne_median, group =taken, 
                         colour = taken)) +
    scale_x_log10(labels = scales::comma,
                  breaks = c(1,10,100,200,1000,3000,10000,100000), 
                  limits= c(1, 350000))+
    scale_y_log10(labels = scales::comma, limits = c(27,1500000),
                  breaks = c(100,1000,10000,100000, 1000000))+
    geom_vline(xintercept = c(200,3000), linetype = c(2,3), alpha = 0.8)+
    # geom_rect(aes(xmin=18000, xmax=25000, ymin=0, ymax=1000000),
    #           colour = NA, fill = 'grey90', alpha = 0.9)+
    geom_ribbon(aes(ymin=Ne_2.5., ymax=Ne_97.5., fill = condition, colour = NULL),
                alpha = 0.1, lwd = 0.1, show.legend = FALSE)+
    scale_color_manual('Sample year', 
                       #values = c('white', 'black', 'lightpink3', 'firebrick'),
                       values = (c('grey35', 'black', 'mistyrose', 'lightpink3')),
                       breaks=c('2006 (bust)', '2007 (boom)','2009 (bust)', '2010 (boom)'))+
    scale_fill_manual('Population phase', values = colcondition2, guide="none")+
    # scale_linetype_manual(values=c("dotdash", "solid"),labels = c('bust', 'boom'))+
    scale_size_manual('Sampling', values = c(0.75, 0.75), guide="none")+
    geom_line(aes(size = sample_point))+ 
    #  facet_grid(~species2)+
    theme_bw()+
    theme(legend.position = 'right',
          #axis.title.x = element_blank(),
          #axis.ticks.length.y = unit(-1, "mm"),
          #axis.title.y = element_blank(),
          axis.text = element_text(size = 8),
          plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm"),
          plot.background = element_rect(
            fill = "white"
          ),
          legend.key = element_rect(fill = "grey80"),
          # axis.text.y = element_text(margin = margin(l = 0, r = 0)),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.background = element_rect(#fill="grey75",
            size=0.5, linetype="solid",
            colour ="grey65")
    )+
    guides(colour = guide_legend(override.aes = list(size=rep(1.25, 4))))+
    ylab(expression("Historical N"[e]*" Estimates"[` log-scale`]))+
    xlab(expression(`Years before sample year`[` log-scale`])); nehistoric
  
  
  
  return(nehistoric)
}

em.ne_trip_fig <- function(necaps2, visual_data, yrain = 230){
  
  
  myfillbooms <- c(`period 1` = '#FFC107',`period 2` = '#D81B60',
                   `period 3` ='#1E88E5')
  myfillbooms <- virid(3)
  names(myfillbooms) <- paste('period', 1:3)
  start_dates <- visual_data$dates
  grob_periods <- grobTree(textGrob((start_dates$subset),
                                    x=c(0.14, 0.365,0.79)+0.005,  
                                    y= 0.9825,
                                    hjust=0, rot = 0,     #was fontsize = 12
                                    gp=gpar(col="black", fontsize=10, fontface="italic")))
  subset_data_visualise <- visual_data$subset
  
  adjustRain <- yrain
  
  neggSy<- ggplot(necaps2, aes(x = ymd(trip), y = ne)) +
    geom_area(data = subset_data_visualise, 
              aes(x = trip, y = rain*adjustRain,
                  group = subset, fill = subset), 
              alpha = 0.2) +
    geom_area(alpha = 0.9, fill = "grey90", col = "grey80")+
    # geom_bar(stat = 'identity', width = 40*1.5,
    #         fill = 'grey90', col = 'grey80', alpha = 0.9)+
    geom_bar(stat = "identity", alpha = 0.8, width = 40,
             col="grey20", fill = "black")+
    #geom_line()+
    scale_color_manual(values = c(inf = 'grey20',no = 'grey90'))+
    scale_fill_manual(values = c(myfillbooms, inf = 'grey20', no = 'black'))+
    theme_classic()+
    geom_errorbar(aes(ymin = lower, ymax = upper), 
                  colour = 'grey60', 
                  width = 0,
                  size = 0.6) +
    ylim(0, 10000)+
    #geom_point(colour = 'black', size = 3)+
    # scale_y_sqrt(breaks = c(20, 100, 1000,3000, 6000, 10000))+
    # scale_y_log10()+
    scale_x_date(date_breaks = "years" , date_labels = "%Y",
                 limits = ymd(c('2006-01-01', '2018-12-01')))+
    geom_vline(xintercept = start_dates$trip, 
               lty = 2, col = "black", alpha = 0.4, lwd = 0.5) +
    # ylim(0,10000)+
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm")) +
    xlab('Year')+
    #ylab(expression("N"[e]*" Estimates"[` sqrt-scale`]))
    ylab(expression("N"[e]*" Estimates"))
  return(neggSy)
  
}
