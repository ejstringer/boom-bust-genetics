
# model table 

# source ------------------------------------------------------------------

source('./code/functions_tables.R')
source('./code/functions_plotting.R')
source('./code/functions_conceptual.R')

# load --------------------------------------------------------------------

# start by running 01_curate.R

# ne processed data
phne <- read.csv('./data_processed/pherm_ne_estimates_summary.csv')
syne <- read.csv('./data_processed/syoung_ne_estimates_summary.csv')

imgsy <- readPNG("./data/images/syoung_close.png")
imgph <- readPNG("./data/images/pherm_close.png")

syDistribution <- st_read("./data/syoung_distribution/data_0.shp")
phDistribution <- st_read("./data/phermann_distribution/data_0.shp")
ausOutline <- st_read("./data/australia/Aus.outline.shp")
grids <- read.csv("./data/em_gridcoortable.csv")

siteBoarder <- raster("./data/images/ANPP_simpson_2000-01-01.tif") 

visual_data <- em.visualisation_subset_period(rain_caps, syoung)

# ho and fst --------------------------------------------------------------
genetic_change <- c(pherm_heterozygosity, 
                    syoung_heterozygosity,
                    pherm_fst)




drift_models <- em.run_models(genetic_change)

estimates <- lapply(drift_models, em.estimates_ci)
estimates$Ph_fst_period1_lm


change_n <- sapply(genetic_change, nrow) # n 

df_ho <- em.model_table(estimates[grep('ho', names(estimates))])  %>% 
  mutate(n = c(rep(change_n[1:3], 2),
               rep(change_n[4:6], 2))) %>% 
  relocate(species, response, period, n)

df_fst <- em.model_table(estimates[grep('fst', names(estimates))]) %>% 
  mutate(n = change_n[7:9]) %>% 
  relocate(species, response, period, n)


df_ho$significant[12] <- sub('TRUE', 'FALSE', df_ho$significant[12]) 

# ne ----------------------------------------------------------------------
syne <- filter(syne,!is.infinite(lower),
               !is.infinite(upper))
ne_cap_models <- lapply(list(Ph_ne_caps_lm = phne, Sy_ne_caps_lm = syne), 
                        function(x) lm(ne.log ~ captures.log, data = x)) 

df_ne <- rbind(em.ne.model.extract(ne_cap_models$Ph_ne_caps_lm),
               em.ne.model.extract(ne_cap_models$Sy_ne_caps_lm,
                                   sp = 'Sy', p = '2006-2017')) %>% 
  mutate(n = c(length(phne$captures.log),
               length(syne$captures.log))) %>% 
  relocate(n, .after = period)

# table 1 ------------------------------------------------------------
fxtb <- em.flextable_table(df_ho, df_fst, df_ne)
fxtbNice <- fxtb %>% 
  compose(
    part = "body", i = c(5,19), j = 2,
    value = as_paragraph(as_b("period"), as_sub("1,2,3"))
  ) %>% 
  compose(
    part = "body", i = c(1,5, 19), j = 1,
    value = as_paragraph(as_b("species"))
  ) %>% 
  compose(
    part = "body", i = c(1,5, 19), j = 3,
    value = as_paragraph(as_b("n"))
  ) %>% 
  compose(
    part = "body", i = c(1,5, 19), j = 4,
    value = as_paragraph(as_b("β"), as_sub("0"))
  )%>% 
  compose(
    part = "body", i = c(1,5, 19), j = 5,
    value = as_paragraph(as_b("β"), as_sub("1"))
  )%>% 
  compose(
    part = "body", i = c(1,5,19), j = 6,
    value = as_paragraph(as_b("β"), as_sub("2"))
  )%>% 
  color(i = 1, j = 6, 'grey80', part = "body", source = j) %>% 
  compose(
    part = "body", i = c(1,5, 19), j = 7,
    value = as_paragraph(as_b("R"), as_sup("2"))
  )  %>% 
  compose(
    part = "body", i = c(1), j = 2,
    value = as_paragraph(as_b("period"), as_sub("YEARS"))
  )  %>% 
  fontsize(part = 'body', c(1,4, 5, 18,19), size = 12) %>% 
  fontsize(part = 'header', size = 12)
fxtbNice


fxtbNice %>% save_as_docx(path = "./figures/table1_models.docx")


# figure 1 ----------------------------------------------------------------

## conceptual -------------
myvisdates <- visual_data$dates
myvisdates$subset <- paste('Period', 1:3)
conceptual_figure <- em.conceptual_plot(rain_caps, 
                                        myvisdates, 
                                        visual_data$subset,
                                        imgph,
                                        imgsy, subsetCol = virid(3))

## aus map ------------------

### create sample site box ----

e <- as(extent(siteBoarder), 'SpatialPolygons')
ebox <- st_bbox(e)
Poly_Coord_df <- data.frame(lon = c(ebox[1], ebox[3]),
                            lat = c(ebox[2], ebox[4]))

pol = st_polygon(
  list(
    cbind(
      Poly_Coord_df$lon[c(1,2,2,1,1)], 
      Poly_Coord_df$lat[c(1,1,2,2,1)])
  )
)


polc <- st_sfc(pol, crs=st_crs(ausOutline))

nt <- st_crop(ausOutline[7,1], (st_bbox(polc)))

### legend ---------
ausmap_legend <- ggplot(data.frame(a = c(1:2), 
                                   fieldsite = 2:3,
                                   g = c(' field site'),
                                   distribution = c('a','b')),
                        aes(x = fieldsite, y = a))+
  geom_bar(stat = 'identity', aes(fill = distribution), alpha = 0.5)+
  geom_point(shape = 15,aes(x = fieldsite, y =a, colour = g), size = 2)+
  scale_colour_manual(values = 'red', name = NULL)+
  scale_fill_manual(values = c('grey', 'orange'),
                    labels = c(expression(italic('P. hermannsburgensis')),
                               expression(italic('S. youngsoni'))),
                    name = 'Distribution:')+
  theme_classic()+
  theme(legend.text.align = 0,
        legend.spacing = unit(0, 'cm'),
        legend.margin = margin(0,0,0.15,0, unit = 'cm'),
        legend.key.size = unit(0.4, 'cm'))
### aus map ---------
ausMap_con <-ggplot() + theme_void() +
  geom_sf(data = phDistribution, size = 0.1, alpha = 1,
          color = "grey", fill = "grey") + 
  geom_sf(data = syDistribution, size = 0.1, alpha = 0.5,
          color = "orange", fill = "orange") + 
  geom_sf(data = ausOutline[1:8,1], size = 0.5,
          color = "black", fill = "red", alpha =0.01) + 
  geom_sf(data = polc, size = 1, 
          color = "red", fill = "red", alpha = 0.9) + 
  coord_sf()+
  labs(tag = 'A')+
  theme(plot.tag = element_text(size = 16))



## big plot ----------
ausleg <-g_legend(ausmap_legend)

lay <- rbind(c(NA,NA,NA,4,4,4,4,3,3),
             c(1,1,1,4,4,4,4,3,3),
             c(1,1,1,1,4,4,4,4,4),
             c(1,1,1,1,4,4,4,4,4),
             c(1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1))

tiff("./figures/figure1_conceptual.tiff", units="cm", width=18,
     height=16, res=600)
grid.arrange(conceptual_figure,ausleg, ausMap_con,  layout_matrix = lay)
dev.off()


# figure 2 ----------------------------------------------------------------

necaps2 <- bind_rows(phne, syne)
necaps2$abline <- NA
necaps2$abline[necaps2$species == 'ph'] <- predict(ne_cap_models$Ph_ne_caps_lm,
                                                   necaps2[necaps2$species == 'ph', c('trip', 'captures.log')])
necaps2$abline[necaps2$species == 'sy'] <- predict(ne_cap_models$Sy_ne_caps_lm,
                                                   necaps2[necaps2$species == 'sy', c('trip', 'captures.log')])
necaps2$species2 <- paste0('z', necaps2$species,2)
necaps2$species2 <- ifelse(necaps2$species=='ph', 'P. hermannsburgensis','S. youngsoni')
necaps2$line <- 1
necaps2$ratio <- log(necaps2$ne/necaps2$captures)
necaps2$ratio <- (necaps2$ne.log/necaps2$captures.log); necaps2$ratio <- necaps2$ratio + abs(min(necaps2$ratio, na.rm = T))

necaps2$ratio[necaps2$species == 'sy'] <- mean(necaps2$ratio, na.rm = T)

nemodelgg<- ggplot(necaps2, aes(x = captures, y = (ne), group = species2))+
  theme_bw()+
  geom_point(alpha = 0.6, size = 2)+
  geom_line(data = necaps2, aes(x = captures, y = exp(abline), linetype= species2, #colour = species2,
  ), size = 1, colour = 'red'
  )+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text.align = 0,
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        # plot.margin = margin(0,0,0,0, "cm"),
        strip.text = element_text(face = "italic"),
        strip.background = element_blank())+
  scale_y_continuous(trans = "log",
                     labels = scales::number_format(accuracy= 1))+
  scale_x_continuous(trans = "log",
                     labels = scales::number_format(accuracy= 0.1))+
  ylab(expression("N"[e]*" Estimates "[` log-scale`]))+
  scale_linetype(guide="none") +
  scale_size(guide="none") +
  facet_wrap(~species2, ncol = 2) +
  xlab(expression(`Mean captures (per 100 trap nights) `[`log-scale`]))
nemodelgg

tiff("./figures/figure2_NeCap_model.tiff", units="cm", width=11,
     height=7, res=600)

nemodelgg

dev.off()


# figure 3 ----------------------------------------------------------------

## model predict ----------------------------------------------------------------
names(estimates)
best_models <- estimates[c('Ph_ho_period1_lm',
                           'Ph_ho_period2_quad',
                           'Ph_ho_period3_lm',
                           'Sy_ho_period1_quad',
                           'Sy_ho_period2_lm',
                           'Sy_ho_period3_lm',
                           'Ph_fst_period1_lm',
                           'Ph_fst_period2_lm',
                           'Ph_fst_period3_lm')]

length(best_models)
model_pred <- lapply(best_models, em.model_predict)
names(model_pred) <- names(best_models)

lapply(model_pred, function(x) x$sig[1])

## plots ------------------------------------------------------------------------
## level 1 ---------------------------------------------------------------------


g_data <- c(pherm_heterozygosity, 
            syoung_heterozygosity, 
            pherm_fst)
plots_level1 <- list()

for(i in 1:length(model_pred)){
  plots_level1[[i]] <-  em.plot_level1(g_data, model_pred, 
                                       model = names(model_pred)[i])
}

names(plots_level1) <- names(model_pred)
## level 2 ---------------------------------------------------------------------


gp1 <- em.plot_level2(plots_level1[1:3], xticks = F)
gp2 <- em.plot_level2(plots_level1[4:6], xticks = F)
gp3 <- em.plot_level2(plots_level1[7:9], xticks = T, ylabFst = T)

## level 3 ---------------------------------------------------------------------

periods <- lapply(1:3, em.plot_labels)


bottom <- grid::textGrob(expression(italic(paste("Years since boom (t)                 "))),
                         gp = grid::gpar(fontsize = 12))


speciesName <- grobTree(textGrob('P. hermannsburgensis',
                                 x=0.01,  
                                 y= 0.34-0.05,
                                 hjust=0, rot = 0,     #was fontsize = 12
                                 gp=gpar(col="black", fontsize=9, fontface="italic")))
speciesName2 <- grobTree(textGrob('S. youngsoni',
                                  x=0.2,  
                                  y= 0.32-0.05,
                                  hjust=0, rot = 0,     #was fontsize = 12
                                  gp=gpar(col="black", fontsize=9, fontface="italic")))

mammalName <- grobTree(textGrob('rodent',
                                x=0.05+0.25,  
                                y= 0.34,
                                hjust=0, rot = 0,     #was fontsize = 12
                                gp=gpar(col="black", fontsize=10, fontface = 'bold')))
mammalName2 <- grobTree(textGrob('marsupial',
                                 x=0.2+0.05,  
                                 y= 0.32,
                                 hjust=0, rot = 0,     #was fontsize = 12
                                 gp=gpar(col="black", fontsize=10, fontface = 'bold')))




pic <- ggplot() + theme_void() +
  annotation_custom(rasterGrob(imgph, interpolate=TRUE))+
  annotation_custom(speciesName)+
  annotation_custom(mammalName)
pic2 <- ggplot() + theme_void() +
  annotation_custom(rasterGrob(imgsy, interpolate=TRUE)) +
  annotation_custom(speciesName2)+
  annotation_custom(mammalName2)


w <- 12
lay <- rbind(c(7,7, 7, 7 ,8,8,8, 8,9,9,9,9,NA, NA, NA),
             c(rep(1, w), 4,4,4),
             c(rep(1, w), 4,4,4),
             c(rep(1, w), 4,4,4),
             c(rep(1, w), 4,4,4),
             c(rep(1, w), 4,4,4),
             c(rep(1, w), 4,4,4),
             c(rep(2, w), 5,5,5),
             c(rep(2, w), 5,5,5),
             c(rep(2, w), 5,5,5),
             c(rep(2, w), 5,5,5),
             c(rep(2, w), 5,5,5),
             c(rep(2, w), 5,5,5),
             c(rep(3, w), 6,6,6),
             c(rep(3, w), 6,6,6),
             c(rep(3, w), 6,6,6),
             c(rep(3, w), 6,6,6),
             c(rep(3, w), 6,6,6),
             c(rep(3, w), 6,6,6))



biggerPlot <-grid.arrange(grobs = list( gp1, gp2,gp3,
                                        pic, pic2, pic, 
                                        periods[[1]], 
                                        periods[[2]],
                                        periods[[3]]),
                          bottom = bottom,
                          layout_matrix = lay)

tiff("./figures/figure3_Ho_Fst.tiff", units="cm", width=18,
     height=20, res=600)
grid.arrange(grobs = list( gp1, gp2,gp3,
                           pic, pic2, pic, 
                           periods[[1]], 
                           periods[[2]],
                           periods[[3]]),
             bottom = bottom,
             layout_matrix = lay)

dev.off()


# figure 4 ----------------------------------------------------------------

sy07 <- read.delim('./data_processed/stairway_summaries/Ne against time.final.syRain07.summary', sep = '\t')
sy06 <- read.delim('./data_processed/stairway_summaries/Ne against time.final.syDry06.summary', sep = '\t')
sy09 <- read.delim('./data_processed/stairway_summaries/Ne against time.final.syDry09.summary', sep = '\t')
sy10 <- read.delim('./data_processed/stairway_summaries/Ne against time.final.syRain10.summary', sep = '\t')

down07 <-  read.delim('./data_processed/stairway_summaries/Ne against time.final.down_phBoom07.summary', sep = '\t')
down06 <-  read.delim('./data_processed/stairway_summaries/Ne against time.final.down_phBust06.summary', sep = '\t')
down207 <- read.delim('./data_processed/stairway_summaries/Ne against time.final.down2_phBoom07.summary', sep = '\t')
down206 <- read.delim('./data_processed/stairway_summaries/Ne against time.final.down2_phBust06.summary', sep = '\t')

ph10 <- read.delim('./data_processed/stairway_summaries/Ne against time.final.phBoom10.summary', sep = '\t')
ph09 <- read.delim('./data_processed/stairway_summaries/Ne against time.final.phbust09.summary', sep = '\t')
ph07 <- read.delim('./data_processed/stairway_summaries/Ne against time.final.phBoom07.summary', sep = '\t')
ph06 <- read.delim('./data_processed/stairway_summaries/Ne against time.final.phBust06.summary', sep = '\t')


sy07$taken <- 'sy2007_n32'
sy06$taken <- 'sy2006_n32'
sy09$taken <- 'sy2009_n48'
sy10$taken <- 'sy2010_n37'
ph09$taken <- 'ph2009_n32'
ph06$taken <- 'ph2006_n39'
ph07$taken <- 'ph2007_n183'
ph10$taken <- 'ph2010_n95'


down07$taken  <- 'down07_n32'
down06$taken  <- 'down06_n32'
down207$taken <- 'downtwo07_n32'
down206$taken <- 'downtwo06_n32'


sy67 <- rbind(sy06, sy07, sy10, sy09, ph09, ph06,ph07,ph10, down06,down07, down206, down207)

sy67$species <- 'ph' 
sy67$species[grep('sy',sy67$taken)] <- 'sy'
sy67$species[grep('down',sy67$taken)] <- 'ph down'


sy67$species2 <- ifelse(grepl('ph', sy67$species), 'ph', 'sy')

sy67$condition <- ifelse(grepl('2006', sy67$taken) |grepl('2009', sy67$taken),
                         'dry period', 'resource pulse')

sy67$sample_year <- str_sub(sy67$taken, start = 3, end = 6)
sy67$sample_point <- ifelse(sy67$sample_year < 2008, '2006-2007', '2009-2010' )


colref <- c(ph = 'black',`ph down` = 'grey90', sy = 'black')
colcondition <- c(`dry period` = "#E69F00", 
                  `resource pulse` = "#009E73")

hist.ph <- em.historic.Ne(sy67, subsetspp. = 'ph',bw = T);hist.ph
hist.ph.down <-em.historic.Ne(sy67, subsetspp. = 'ph down')
hist.sy <- em.historic.Ne(sy67, subsetspp. = 'sy', bw = T)



speciesName <- grobTree(textGrob('P. hermannsburgensis',
                                 x=0.04,  
                                 y= 0.12-0.03+0.1,
                                 hjust=0, rot = 0,     #was fontsize = 12
                                 gp=gpar(col="black", fontsize=8, fontface="italic")))
speciesName2 <- grobTree(textGrob('S. youngsoni',
                                  x=0.21,  
                                  y= 0.06-0.01+0.1,
                                  hjust=0, rot = 0,     #was fontsize = 12
                                  gp=gpar(col="black", fontsize=8, fontface="italic")))


mammalName <- grobTree(textGrob('rodent',
                                x=0.05+0.24,  
                                y= 0.180+0.1,
                                hjust=0, rot = 0,     #was fontsize = 12
                                gp=gpar(col="black", fontsize=9, fontface = 'bold')))
mammalName2 <- grobTree(textGrob('marsupial',
                                 x=0.2+0.04,  
                                 y= 0.14+0.1,
                                 hjust=0, rot = 0,     #was fontsize = 12
                                 gp=gpar(col="black", fontsize=9, fontface = 'bold')))




pic <- ggplot() + theme_void() +
  annotation_custom(rasterGrob(imgph, interpolate=TRUE))+
  annotation_custom(speciesName)+
  annotation_custom(mammalName)
pic2 <- ggplot() + theme_void() +
  annotation_custom(rasterGrob(imgsy, interpolate=TRUE)) +
  annotation_custom(speciesName2)+
  annotation_custom(mammalName2)
w <- 7
lay <- rbind(c(rep(1,w), 4,4),
             c(rep(1,w), 4,4),
             c(rep(1,w), 4,4),
             c(rep(1,w), 4,4),
             c(rep(1,w), NA,NA),
             c(rep(1,w), 5,5),
             c(rep(2,w), 5,5),
             c(rep(2,w), NA,NA),
             c(rep(2,w), 3,3),
             c(rep(2,w), 3,3),
             c(rep(2,w), 3,3),
             c(rep(2,w), 3,3),
             c(rep(2,w), NA,NA))


mylegend<-g_legend(hist.ph)


left <- grid::textGrob(expression("Historical N"[e]*" Estimates"[` log-scale`]),
                       gp = grid::gpar(fontsize = 12), rot = 90)

tiff("./figures/figure4_historicalNe.tiff", units="cm", width=18,
     height=12, res=600)

grid.arrange(hist.ph + theme(legend.position="none",
                             axis.text.x = element_blank(),
                             axis.title= element_blank(),
                             axis.ticks.x = element_blank(),) , 
             hist.sy + theme(legend.position="none",
                             axis.title.y = element_blank()), 
             pic2,  pic, 
             mylegend, left = left,
             layout_matrix = lay)

dev.off()

