
# plotting

em.visualisation_subset_period <- function(rain_caps, syoung){
  
  max_y <- max(rain_caps$captures, na.rm = T)
  
  subset_dates <- rain_caps %>% 
    filter(trip >= min(ymd(syoung@other$ind.metrics$trip)),
           sinceEvent != 'period 0') %>%
    group_by(sinceEvent) %>% 
    summarise(trip = min(trip)) %>% 
    #filter(complete.cases(sinceEvent)) %>% 
    arrange(trip) %>% 
    rename(subset = sinceEvent)
  
  
  vis_dates <- c(subset_dates$trip, 
                 max(ymd(syoung@other$ind.metrics$trip))+months(2))
  
  vis_comb <- combn(1:length(vis_dates),2) 
  vis_comb_seq <-t(vis_comb[,vis_comb %>% diff == 1])
  
  subset_data_visualise <- data.frame(subset = subset_dates$subset,
                                      trip1 = ymd(vis_dates[vis_comb_seq[,1]]),
                                      trip2 = ymd(vis_dates[vis_comb_seq[,2]])) %>% 
    mutate(trip3 = ymd(trip1) - days(1),
           trip4 = ymd(trip2) - days(1)) %>% 
    pivot_longer(cols = trip1:trip4, values_to = "trip") %>% 
    arrange(subset, trip) %>% 
    mutate(rain = rep(c(0, rep(max_y, 2), 0), nrow(subset_dates))) %>% 
    dplyr::select(-name)
  
  return(list(subset = subset_data_visualise,
              dates = subset_dates))
  
}


em.conceptual_plot <- function(rain_caps,
                               start_dates,
                               subset_data_visualise,
                               imgph, imgsy,
                               subsetCol = c('#FFC107', '#D81B60','#1E88E5')){
myfillbooms <- subsetCol #rev(c('#1E88E5', '#D81B60','#FFC107'))


grob_periods <- grobTree(textGrob((start_dates$subset),
                                  x=c(0.14, 0.365,0.79)+0.005,  
                                  y= 0.97,
                                  hjust=0, rot = 0,     #was fontsize = 12
                                  gp=gpar(col="black", fontsize=10, fontface="italic")))



g1 <- ggplot(filter(rain_caps, trip > ymd("2006-01-01")), aes(ymd(trip), captures))+
  geom_area(data = subset_data_visualise, 
            aes(x = trip, y = rain*1.16, group = subset, fill = subset), 
            alpha = 0.2) +
  geom_area(alpha = 0.9, fill = "grey90", col = "grey80")+
  geom_bar(stat = "identity", alpha = 0.8, col="grey20", fill = "black")+
  scale_fill_manual(values = myfillbooms) +
  ylim(0,50)+
  theme_classic()+
  xlab("Year") +
  ylab("mean captures")+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.y = element_text(vjust = 0.3,
                               margin = margin(l = 15, r = -17)),
    plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm")) +
  scale_x_date(date_breaks = "years" , date_labels = "%Y") +
  geom_vline(xintercept = start_dates$trip, 
             lty = 2, col = "black", alpha = 0.4, lwd = 0.5) +
  annotation_custom(grob_periods)


g2 <- ggplot(filter(rain_caps, trip > ymd("2006-01-01")), aes(ymd(trip), capturesSy))+
  geom_area(data = subset_data_visualise, 
            aes(x = trip, y = rain/6, group = subset, fill = subset), 
            alpha = 0.2) +
  geom_area(alpha = 0.9, fill = "grey90", col = "grey80")+
  geom_bar(stat = "identity", alpha = 0.8, col="grey20", fill = "black")+
   scale_fill_manual(values = myfillbooms) +
  theme_classic()+
  xlab("Year") +
  ylab("mean captures")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.position = "none",
        axis.text.y = element_text(vjust = -0.05,
                                   margin = margin(l = 15, r = -12)),
        plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm")) +
  scale_x_date(date_breaks = "years" , date_labels = "%Y") +
  scale_y_continuous(breaks = c(0,2,4))+
  geom_vline(xintercept = start_dates$trip, 
             lty = 2, col = "black", alpha = 0.4, lwd = 0.5)


myfillboomsWhite <- rep("white", 3)
g3 <- ggplot(filter(rain_caps, trip > ymd("2006-01-01")), 
             aes(ymd(trip), rain))+
  geom_area(data = subset_data_visualise, 
            aes(x = trip, y = rain*7.5, fill = "white"), 
            alpha = 0.2) +
  geom_bar(stat = "identity", alpha = 0.8, col="lightblue", fill = "darkblue")+
  scale_fill_manual(values = myfillboomsWhite) +
  theme_classic()+
  xlab("Year") +
  ylab("Rainfall (mm)")+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.length.y = unit(1, "mm"),
    axis.text.y = element_text(vjust = 0.3,
                               margin = margin(l = 13, r = -20)),
    legend.position = "none",
    plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm")) +
  scale_x_date(date_breaks = "years" , date_labels = "%Y") 


speciesName <- grobTree(textGrob('P. hermannsburgensis',
                                  x=0,  
                                  y= 0.28-0.18,
                                  hjust=0, rot = 0,     #was fontsize = 12
                                  gp=gpar(col="black", fontsize=9, fontface="italic")))
speciesName2 <- grobTree(textGrob('S. youngsoni',
                                 x=0.15+0.02,  
                                 y= 0.34-0.18,
                                 hjust=0, rot = 0,     #was fontsize = 12
                                 gp=gpar(col="black", fontsize=9, fontface="italic")))
mammalName <- grobTree(textGrob('rodent',
                                x=0.0+0.25,  
                                y= 0.34 -0.14,
                                hjust=0, rot = 0,     #was fontsize = 12
                                gp=gpar(col="black", fontsize=10, fontface = 'bold')))
mammalName2 <- grobTree(textGrob('marsupial',
                                 x=0.15+0.05,  
                                 y= 0.32 - 0.09,
                                 hjust=0, rot = 0,     #was fontsize = 12
                                 gp=gpar(col="black", fontsize=10, fontface = 'bold')))



g <- rasterGrob(imgph, interpolate=TRUE) 
gsy <- rasterGrob(imgsy, interpolate=TRUE) 
pic <- ggplot() + theme_void() +
  annotation_custom(g) + 
  annotation_custom(speciesName)+
  annotation_custom(mammalName)


pic
pic2 <- ggplot() + theme_void() +
  annotation_custom(gsy)+
  annotation_custom(speciesName2)+
  annotation_custom(mammalName2)
pic2

lay <- rbind(c(rep(1, 6),NA, NA),
             c(rep(1, 6),NA, NA),
             c(rep(1, 6),NA, NA),
             c(rep(1, 6),NA, NA),
             c(rep(2, 6),4, 4),
             c(rep(2, 6),4, 4),
             c(rep(2, 6),NA, NA),
             c(rep(2, 6),5, 5),
             c(rep(3, 6),5, 5),
             c(rep(3, 6),5, 5))

lay <- rbind(c(rep(2, 6),5,5),
             c(rep(2, 6),5,5),
             c(rep(2, 6),NA, NA),
             c(rep(2, 6),3,3),
             c(rep(4, 6),3,3),
             c(rep(4, 6),3,3))

lay2 <- rbind(c(rep(1, 6),NA, NA),
             c(rep(1, 6),NA, NA),
             c(rep(1, 6),NA, NA),
             c(rep(1, 6),NA, NA),
             c(rep(2, 8)),
             c(rep(2, 8)),
             c(rep(2, 8)),
             c(rep(2, 8)),
             c(rep(2, 8)),
             c(rep(2, 8)))

fsize <- 16

left <- grid::textGrob(expression('Mean captures (per 100 trap nights)'),
                       gp = grid::gpar(fontsize = 11), rot = 90)
cap_plot <- grid.arrange(grobs = list(g1+ labs(tag = 'C')+ 
                                  theme(plot.tag = element_text(size = fsize),
                                        axis.title.y = element_blank()),
                                  pic2,
                                      g2+ labs(tag = 'C')+ 
                              theme(plot.tag = element_text(size = fsize, colour = 'white'),
                                    axis.title.y = element_blank()),
                              pic),
                         left = left, layout_matrix = lay)

left2 <- grid::textGrob('Rainfall (mm)',
                       gp = grid::gpar(fontsize = 12), rot = 90)
rain_plot <- grid.arrange(grobs = list(g3+labs(tag = "B")+
                                         theme(plot.tag = element_text(size = fsize),
                                               axis.title.y = element_blank())),
                          left = left2)
ppraw <-grid.arrange(grobs = list(rain_plot,
                                  # g1+labs(tag = "B")+
                                  #   theme(plot.tag = element_text(size = fsize)),
                                  # g2+labs(tag = "")+
                                  #   theme(plot.tag = element_text(size = fsize)), 
                                  # pic, pic2
                                  cap_plot
                                  ),
                     layout_matrix = lay2)

return(ppraw)
}





