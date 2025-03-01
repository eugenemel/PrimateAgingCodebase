# Define species list for survival plot output
SELECTED_SPECIES <- c("Pygmy Marmoset",
            "Lesser Bushbaby",   
            "Red-ruffed Lemur", 
            "Ringtailed Lemur",  
            "Cottontop Tamarin", 
            "Common Marmoset",  
            "Vervet",  
            "Rhesus Indian-derived",
            "Orangutan",
            "Western Lowland Gorilla",
            "Chimpanzee")

# Create survival plot for selected species
pdf("plots/top10_survival_plot.pdf", width=9, height=8)
  # Filter data for selected species
  tmp <- ALL %>% 
  filter(species %in% SELECTED_SPECIES)
  
  tmp$species <- factor(tmp$species, levels=SELECTED_SPECIES)
  
  # Create survival object and fit
  SurvObj <- Surv(tmp$lifespan, tmp$event)
  fit <- survfit(SurvObj ~ species, data=tmp)
  
  # Clean up labels
  names(fit$strata) <- sub("species=", "", names(fit$strata))
  
  # Order species by median survival
  speciesOrder <- rownames(as.data.frame(quantile(fit, probs=c(0.5))) %>% 
               arrange(X50))
  tmp$species <- factor(tmp$species, levels=speciesOrder)
  
  # Plot survival curves
  p0 <- survminer::ggsurvplot(fit, 
               conf.int=FALSE, 
               risk.table=FALSE,
               break.time.by=5, 
               fontsize=1.5, 
               legend="bottom",
               legend.title="", 
               xlab="Age",
               palette="", 
               title="")
  print(p0)
dev.off()

# Create hazard rate plot
# Get top species by event count
TOPSPECIES <- ALL %>% 
  group_by(species) %>% 
  summarise(N=n(), nevents=sum(event)) %>%  
  arrange(desc(nevents))

pdf("plots/primate_hazards.pdf", width=7.5, height=5.5)
  theme_set(theme_bw())
  
  # Set up phylogenetic order
  TREE_ORDER <- data.frame(species=paste(primate_tree$tip.label))
  TREE_ORDER$tree_order <- 1:nrow(TREE_ORDER)
  
  # Filter data for species with >10 observations
  LT$species_tree_order <- factor(LT$species, level=paste(rev(TREE_ORDER$species)))
  tmp <- TOPSPECIES %>% dplyr::filter(N > 10)
  tmp <- LT %>% dplyr::filter(species %in% tmp$species)
  
  # Custom formatting function for y-axis
  scaleFUN <- function(x) sprintf("%.0f", x)
  
  # Create hazard plot
  p0 <- ggplot(tmp, aes(x=time, y=log_rate, color=sex, fill=species_tree_order)) + 
  geom_point(alpha=0.7) + 
  scale_y_continuous(labels=scaleFUN, limits=c(-5,0)) + 
  xlim(1,50) + 
  xlab("Age") + 
  ylab("Log Hazard") + 
  geom_line(aes(x=time, y=rlm_fitted)) + 
  facet_wrap(~species_tree_order, scales="free", strip.position="left") +
  theme_bw(base_size=8) +
  theme(strip.background=element_rect(fill="white"), legend.position="none") + 
  scale_color_brewer(palette="Set1")

  print(p0)
  remove(p0)
dev.off()

