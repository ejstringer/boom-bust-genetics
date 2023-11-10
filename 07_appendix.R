

# load --------------------------------------------------------------------

# run previous R files prior to plotting appendix figures and tables



# figure S1 ---------------------------------------------------------------

ausMap<-ggplot() + theme_void() +
  theme(plot.background = element_rect(fill = 'white', colour = 'black'))+
  geom_sf(data = ausOutline[1:8,1], size = 0.5,
          color = "black", fill = "red", alpha =0.01) + 
  geom_sf(data = polc, size = 1, 
          color = "red", fill = "red", alpha = 0.9) + 
  coord_sf()

studySite <-ggplot() + theme_void() +
  geom_sf(data = nt, size = 0.5, color = "black", fill = "grey90") + 
  geom_sf(data = polc, size = 1, color = "red", fill = "grey90", alpha = 0.1) + 
  geom_point(data = grids, aes(x = lon, y = lat, colour = 'grid   '),
             alpha = 0.5, size = 2.5)+
  annotation_scale(location = "bl", width_hint = 0.47, pad_x = unit(0.60, "cm"))+
  annotation_north_arrow(pad_x = unit(0.5, "cm"),
                         pad_y = unit(0.75, "cm"),
                         style = north_arrow_fancy_orienteering)+
  scale_colour_manual(values = 'black', name = NULL)+
  annotation_custom(grobTree(textGrob(expression(italic("Northern Territory")),
                                      x= 0.23, y= 0.3,
                                      rot = 90, 
                                      gp = gpar(fontsize=11)),
                             textGrob(expression(italic("Queensland")),
                                      x= 0.28, y= 0.55,
                                      rot = 90, 
                                      gp = gpar(fontsize=11))))+
  theme(plot.margin = margin(0.25,0.25, 0.25, 0.25, "cm"),
        legend.position = 'none',
        legend.background = element_rect(colour = 'black'));studySite

study_legend<- ggplot(data.frame(a = c(1:2), 
                                 fieldsite = 2:3,
                                 g = c('Field site'),
                                 distribution = c('a','b')),
                      aes(x = fieldsite, y = a))+
  geom_point(shape = 15,aes(x = fieldsite, y =a, colour = g), size = 2)+
  geom_point(shape = 21, aes(x = fieldsite, y =a, fill= 'Grid'),colour = 'black',
             size = 2.5)+
  scale_colour_manual(values = c('red'), name = NULL)+
  scale_fill_manual(values = c('grey50'), name = NULL)+
  theme_classic()+
  theme(legend.text.align = 0,
        legend.text = element_text(size = 8.5),
        legend.spacing = unit(0, 'cm'),
        legend.margin = margin(0,0,0,0, unit = 'cm'),
        legend.key.size = unit(0.2, 'cm'),
        legend.key.height = unit(0, 'cm'));study_legend

mylegend<-g_legend(study_legend)

lay <- rbind(c(2,2,2,2,2,2,1,1,1,1,1,1),
             c(2,2,2,2,2,2,1,1,1,1,1,1),
             c(2,2,2,2,2,2,1,1,1,1,1,1),
             c(3,3,3,3,2,2,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1,1,1,1))
# Small: approximately 9 cm x 6 cm
# Medium: approximately 11 cm x 11 cm
# Large: approximately 18 cm x 22 cm

tiff("./figures/figureS1_study_site.tiff", units="cm", width=14,
     height=14, res=600)
grid.arrange(studySite, ausMap,mylegend, layout_matrix = lay)
dev.off()

# figure S2 ---------------------------------------------------------------

netripSy <-  em.ne_trip_fig(filter(syne, !is.infinite(lower),
                                   !is.infinite(upper)), visual_data, 230)
netripSy



netripPh <- em.ne_trip_fig(phne, visual_data, 225)+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
netripPh 

grob_periods <- grobTree(textGrob(paste('Period', 1:3),
                                  x=c(0.14, 0.365,0.77)+0.005,  
                                  y= 0.9825-0.05,
                                  hjust=0, rot = 0,     #was fontsize = 12
                                  gp=gpar(col="black", fontsize=10, fontface="italic")))



speciesName <- grobTree(textGrob('P. hermannsburgensis',
                                 x=0.0,  
                                 y= 0.12+0.06,
                                 hjust=0, rot = 0,     #was fontsize = 12
                                 gp=gpar(col="black", fontsize=9, fontface="italic")))
speciesName2 <- grobTree(textGrob('S. youngsoni',
                                  x=0.2,  
                                  y= 0.06+0.06,
                                  hjust=0, rot = 0,     #was fontsize = 12
                                  gp=gpar(col="black", fontsize=9, fontface="italic")))


mammalName <- grobTree(textGrob('rodent',
                                x=0.05+0.25,  
                                y= 0.270,
                                hjust=0, rot = 0,     #was fontsize = 12
                                gp=gpar(col="black", fontsize=9, fontface = 'bold')))
mammalName2 <- grobTree(textGrob('marsupial',
                                 x=0.2+0.05,  
                                 y= 0.22,
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
lay <- rbind(c(rep(1, w), NA,NA),
             c(rep(1, w), 4,4),
             c(rep(1, w), 4,4),
             c(rep(1, w),4,4),
             c(rep(1, w), 4,4),
             c(rep(2, w), NA,NA),
             c(rep(2, w),3,3),
             c(rep(2, w),3,3),
             c(rep(2, w),3,3),
             c(rep(2, w), 3,3),
             c(rep(2, w), NA,NA))


tiff("./figures/figureS2_NeEstimates.tiff", units="cm", width=18,
     height=12, res=600)
grid.arrange(netripPh + 
             annotation_custom(grob_periods),
             netripSy,
             pic2, pic,
             layout_matrix = lay)


dev.off()


# figure S3 ---------------------------------------------------------------

fx <- function(x){
  x$ntrip <- as.numeric(factor(x$trip))
  return(x)
}
pherm_fst1 <- lapply(pherm_fst, fx) 
fst1 <- reduce(pherm_fst1, bind_rows)

best_fst_model <- best_models[grep('fst', names(best_models))]

slopem <- sapply(best_fst_model, function(x) x[grep('m', rownames(x)), 'Estimate'])
interceptm <- sapply(best_fst_model, function(x) x[grep('Int', rownames(x)), 'Estimate'])

mmodel <- data.frame(sinceEvent = paste('period', 1:3), 
                     intercept = interceptm,
                     slope = slopem)
ggdist <- ggplot(fst1, aes(km, exp(fst.log))) +
  geom_point(size = 2, alpha = 0.2)+
  facet_grid(~sinceEvent)+
  scale_fill_gradientn(colours = grey.colors(4))+
  theme_bw()+
  scale_x_log10()+
  scale_y_continuous(tran = 'log',
                     labels = scales::number_format(accuracy =  0.001))+
  geom_abline(aes(intercept = intercept, slope = slope), mmodel,
              colour = 'red', size = 1, lty = 'dashed')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'grey95'))+
  ylab(expression(paste("", F[`ST log-scale`], "")))+
  xlab(expression(paste("", km[` log-scale`], "")))

ggdist

tiff("./figures/figureS3_fst_km.tiff", units="cm", width=16,
     height=9, res=600)

ggdist

dev.off()


# figure S4 ---------------------------------------------------------------

downComp <- filter(sy67, 
                   year > 1,
                   taken %in% c('ph2007_n183',
                                "ph2006_n39") | grepl('down',taken)) %>% 
  mutate(groupyear = ifelse(grepl('06', taken), '2006', '2007'),
         groupn = ifelse(grepl('down', taken), 'down', 'all'),
         grouped = paste(groupn, groupyear))

downComp$grouped %>% table

colcondition <- c('black', 'red')
nehistoric<-ggplot(downComp,
                   aes(year, Ne_median, group =taken, 
                       colour = grouped)) +
  scale_x_log10(labels = scales::comma,
                breaks = c(1,10,100,1000,10000,100000), 
                limits= c(1, 350000))+
  scale_y_log10(labels = scales::comma, limits = c(27,1500000),
                breaks = c(100,1000,10000,100000, 1000000))+
  # geom_vline(xintercept = c(235,3500), linetype = 2, alpha = 0.4)+
  # geom_rect(aes(xmin=18000, xmax=25000, ymin=0, ymax=1000000),
  #           colour = NA, fill = 'grey90', alpha = 0.9)+
  geom_ribbon(aes(ymin=Ne_2.5., ymax=Ne_97.5., fill = grouped, colour = NULL),
              alpha = 0.1, lwd = 0.1,
              # fill = 'grey50', 
              show.legend = FALSE)+
  scale_color_manual('Samples', 
                     #values = c('white', 'black', 'lightpink3', 'firebrick')
                     #values = c('grey35', 'black', 'mistyrose', 'lightpink3')
                     values = c('grey35', 'black', 'mistyrose', 'indianred2')
  )+
  scale_fill_manual(values = c('grey35', 'black', 'indianred3', 'indianred4'))+
  geom_line(lwd = 0.75)+ 
  #  facet_grid(~species2)+
  theme_bw()+
  theme(legend.position = c(0.86, 0.25),
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
          colour ="grey65"))+
  guides(fill = 'none',
         lwd = 'none',
         colour = guide_legend(override.aes = list(size=rep(1.25, 4))))+
  ylab(expression("Historical N"[e]*" Estimates "[` log-scale`]))+
  xlab(expression(`Years before sample year`[` log-scale`])); nehistoric


tiff("./figures/figureS4_historicalNe_down.tiff", units="cm", width=16,
     height=10, res=600)

nehistoric

dev.off()



# table S1 ----------------------------------------------------------------

# source appendix sample table raw

metaph <- ph@other$ind.metrics
metasy <- sy@other$ind.metrics
names(metaph)== names(metasy)

metaBoth <- bind_rows(metaph, metasy)



head(metaBoth)

gridSum <- data.frame(grids = colSums(table(metaBoth$gridId, metaBoth$trip) > 0))
gridSum$trip <- ymd(rownames(gridSum))
meta <- metaBoth %>% 
  dplyr::select(trip, gridId, species) %>% 
  mutate(trip = ymd(trip), year = year(trip),
         monYear = as.yearmon(trip,  "%m/%Y"),
         month = month(trip)) %>% 
  filter(trip > ymd("2006-01-01")) %>% 
  group_by(trip, species) %>% 
  summarise(
    n = n()) %>%
  pivot_wider(names_from = species, values_from = n) %>%
  left_join(gridSum) %>%
  ungroup() %>% 
  # group_by(trip) %>% 
  # summarise(#grids = sum(grids),
  #           `Pseudomys hermannsburgensis` = sum(`Pseudomys hermannsburgensis`, na.rm = T),
  #           `Sminthopsis youngsoni` = sum(`Sminthopsis youngsoni`, na.rm = T)) %>% 
  mutate(trip = as.yearmon(trip,  "%m/%Y"),
         year = year(trip),
         month = month(trip),
         month = month.abb[month],
         year = factor(ifelse(duplicated(year), NA, year)),
         `Pseudomys hermannsburgensis`= ifelse(is.na(`Pseudomys hermannsburgensis`), 0,
                                               `Pseudomys hermannsburgensis`),
         `Sminthopsis youngsoni` = ifelse(is.na(`Sminthopsis youngsoni`), 0,
                                          `Sminthopsis youngsoni`)
  ) %>% 
  relocate(year, month, grids) %>% 
  dplyr::select(-trip)
meta %>% tail

flextable(meta) %>% 
  #  set_caption("Table S1. Number of tissue samples used in analysis for P.hermanns(Ph) and S.youngsoni (Sy") %>% 
  width(width = c(2,2,2,4.5,4.5), unit = "cm") %>% 
  bg(i = seq(1,nrow(meta),2), bg = "grey90") %>% 
  set_header_labels(year = 'Year', month = 'Month', grids = 'Grids sampled') %>% 
  fontsize(i = NULL, j = NULL, size = 9, part = "body") %>% 
  fontsize(i = NULL, j = NULL, size = 11, part = "header") %>% 
  height(i = NULL, height = 0.10, part = "body", unit = "cm") %>%
  # bold(j = 1:2, part = 'header') %>% 
  italic(j = 4:5, part = "header") -> fxtb 
fxtb
save_as_docx(fxtb, path = './figures/tableS1_samples.docx')


# table S2 ----------------------------------------------------------------

necaps2$speciesName <- ifelse(necaps2$species == 'ph',
                              'P.hermannsburgensis',
                              'S.youngsoni')
necaps2$speciesName <- ifelse(duplicated(necaps2$speciesName),
                              NA, necaps2$speciesName)
netable <- necaps2[,c('speciesName','defined', 
                      'trip', 'ne', 'lower',
                      'upper', 'n', 'loci')]

#(Hill 1981; Waples 2006; Waples & Do 2010)

nefxtb <- netable %>% 
  mutate(lower = ifelse(is.na(lower), Inf, lower),
         upper = ifelse(is.na(upper), Inf, upper),
         speciesName = ifelse(grepl('herm', speciesName), 
                              'P.hermanns', speciesName),
         defined = sub('April', 'Apr', defined)) %>% 
  flextable() %>% 
  autofit %>% 
  width(width = c(3,3.5,2,2,1.5,1.5,1,2), unit = "cm") %>% 
  fontsize(i = NULL, j = NULL, size = 9, part = "body") %>% 
  fontsize(i = NULL, j = NULL, size = 11, part = "header") %>% 
  height(i = NULL, height = 0.10, part = "body", unit = "cm") %>% 
  italic(j = 1) %>% 
  hline(i = 21,  border = fp_border(color = "grey40", width = 2)) %>% 
  set_header_labels(speciesName = "species",
                    tripDefined = 'trip definition')
nefxtb
save_as_docx(nefxtb, path = './figures/tableS2_neEstimator.docx')


# table S3 ----------------------------------------------------------------

ndf <-  data.frame(n = sapply(nedataUse, nInd)) %>% 
  cbind(data.frame(loci = sapply(nedataUse, nLoc))) %>% 
  mutate(blueprint = rownames(.),
         species = rep(c('P. hermanns', 'S. youngsoni',
                         'P. hermanns'), each = 4),
         trips = c(sapply(phTrips, paste, collapse = '; '),
                   sapply(syTrips, paste, collapse = '; '),
                   sapply(phTrips[1:2], paste),
                   sapply(phTrips[1:2], paste))) %>% 
  relocate(n, loci, .after = trips)


fx_stair <- flextable(ndf) %>% autofit() %>% 
  flextable::theme_alafoli() %>% italic(j = 2) %>% 
  hline(i = c(4,8,12),
        border = fp_border(color = "grey40", width = 2)) %>% 
  hline(part = 'header', border = fp_border(color = "grey40", width = 3)) %>% 
  bold(part = 'header') 

  save_as_docx(fx_stair, path = './figures/tableS3_stairway.docx')
