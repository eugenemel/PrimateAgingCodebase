ESTIMATE=FALSE; #set to TRUE to re-estimate models

LTsubset =  LT
primate_tree_subset = keep.tip(primate_tree,tip = unique(LTsubset$species))
tree_order=data.frame(species=paste(primate_tree_subset$tip.label),
                      tree_order=1:(primate_tree_subset$Nnode+1))

LTsubset=  LTsubset %>% left_join(tree_order)
LTsubset$speciesId = as.factor(LTsubset$tree_order)
LTsubset = LTsubset %>% filter(rate>0)
LTsubset$sex[ LTsubset$sex == "Both"] <- "Female"; #only a single factor
LTsubset$sd = sqrt(LTsubset$enter*LTsubset$rate*(1-LTsubset$rate))

VCV=ape::vcv.phylo(phy=primate_tree_subset)

if (ESTIMATE) {
  
  model_simple <- brm(
    log_rate ~ 1 + time + sex + time*sex + (1 + time + sex + time*sex | gr(species)),
    data = LTsubset, 
    family = gaussian(),
    #family = student(),
    cores = 4,
    control = list(adapt_delta = 0.8)
  )
  
  
  model_phylo <- brm(
    log_rate ~ 1 + time + sex + time*sex + (1 + time + sex + time*sex | gr(phy=species,cov=A)),
    data = LTsubset, 
    family = gaussian(),
    #family = student(),
    data2 = list(A = VCV),
    cores = 4,
    control = list(adapt_delta = 0.8)
  )
  saveRDS(model_simple, file="outputs/brms_simple.Rds")
  saveRDS(model_phylo, file="outputs/brms_phylo.Rds")
} else {
  model_simple=readRDS(file="outputs/brms_simple.Rds");
  model_phylo=readRDS(file="outputs/brms_phylo.Rds");
} 

#extract simple (noVCV)
brms_coef=coef(model_simple)$'species' %>% as.data.frame()
brms_coef$species = rownames(brms_coef)
brms_alpha= brms_coef %>% mutate(term = "alpha0", model="LifeTable_brms" ) %>% 
          dplyr::select(species, term, 
                            estimate=Estimate.Intercept,
                            std.error=Est.Error.Intercept,
     
                                               model)
brms_beta= brms_coef %>% mutate(term = "beta0", model="LifeTable_brms" ) %>% 
  dplyr::select(species, term, 
                estimate=Estimate.time,
                std.error=Est.Error.time,
                model)

#extract simple (with VCV)
brms_coef=coef(model_phylo)$species %>% as.data.frame()
brms_coef$species = rownames(brms_coef)

brms_alpha_phylo= brms_coef %>% mutate(term = "alpha0", model="LifeTable_brms_phylo" ) %>% 
  dplyr::select(species, term, 
                estimate=Estimate.Intercept,
                std.error=Est.Error.Intercept,
                model)
brms_beta_phylo= brms_coef %>% mutate(term = "beta0", model="LifeTable_brms_phylo" ) %>% 
  dplyr::select(species, term, 
                estimate=Estimate.time,
                std.error=Est.Error.time,
                model)


tmp=predict(model_simple);
LTsubset$brms_simple_fitted=tmp[, 'Estimate' ]
tmp=predict(model_phylo);
LTsubset$brms_phylo_fitted=tmp[, 'Estimate' ]

pdf("plots/primate_hazards_brms.pdf",width=7.5, height=8)
TREE_ORDER=data.frame(species=paste(primate_tree_subset$tip.label));
TREE_ORDER$tree_order =1:nrow(TREE_ORDER);
LTsubset$species_tree_order=factor(LTsubset$species,level=paste(rev(TREE_ORDER$species)))
tmp=SPECIES_SUMMARY %>% dplyr::filter(N>10)
tmp = LTsubset %>% dplyr::filter(species %in% tmp$species )
scaleFUN <- function(x) sprintf("%.0f", x)

p0=ggplot(tmp, aes(x=time,y=log_rate,color=sex,fill=species_tree_order)) + geom_point() + 
  scale_y_continuous(labels=scaleFUN,limits = c(-5,0)) + 
  xlim(1,50) +
  xlab("Age") +
  ylab("Log Hazard") + 
  #geom_errorbar(aes(ymin=log(rate-sd), ymax=log(rate+sd),alpha=1/sd), width=.2, alpha=0.5, position=position_dodge(0.05)) + 
  geom_line(aes(x=time,y=brms_simple_fitted),linetype=3) +
  geom_line(aes(x=time,y=brms_phylo_fitted),linetype=1) + 
  facet_wrap(~species_tree_order, scales="free") +
  theme_bw(base_size = 8) +
  theme(strip.background =element_rect(fill="white"), legend.position = "none") +
  scale_color_brewer(palette="Set1")
print(p0);
remove(p0)
dev.off()

if(0) {
    pp_check(model_simple, type='error_scatter_avg')
    pp_check(model_phylo)
    pp_check(model_simple)
    pp_check(model_phylo,x='time',group='species',type='intervals')
    pp_check(model_phylo,x='time', type='error_scatter_avg_vs_x')
    pp_check(model_phylo,type='stat_2d')
    pp_check(model_simple,type='loo_pit')
    pp_check(model_phylo,type='loo_pit')
    pp_check(model_phylo, type = "stat_grouped", stat = "mean", group = "species")
    pp_check(model_simple, type = "ribbon_grouped", x="time", group = "species")
    pp_check(model_simple, type = "scatter_avg_grouped", x="time", group = "species")
    pp_check(model_simple, type = "ecdf_overlay_grouped", x="time", group = "species")
}

# sensitivity check to tree topology 
if(0) {
  for(i in 1:10) { 
     VCV=get_random_vcv(primate_tree_subset);
     #fit phylogenetic model
     model_phylo <- brm(
       log_rate ~ 1 + time + sex + time*sex + (1 + time + sex + time*sex | gr(phy=species,cov=A)),
       data = LTsubset, 
       family = gaussian(),
       data2 = list(A = VCV),
       cores = 4,
       control = list(adapt_delta = 0.8)
     );
     saveRDS(model_phylo, file=sprintf("outputs/testvcv_brms_phylo_%d.Rds",i));
  }
  
  #load all test vcv models, extract aging parameters, and plot correlation
  all_models=list();
  for(i in 1:10) { 
    model=readRDS(file=sprintf("outputs/testvcv_brms_phylo_%d.Rds",i));
    brms_coef=coef(model)$species %>% as.data.frame()
    brms_coef$species = rownames(brms_coef)
    brms_coef$ittr <- i;
    all_models = bind_rows(all_models, brms_coef);
  }
  tmp=all_models %>% reshape2::dcast(species ~ ittr,value.var = "Estimate.Intercept");
  pairs(tmp[,-1])
  tmp=all_models %>% reshape2::dcast(species ~ ittr,value.var = "Estimate.time");
  pairs(tmp[,-1])
}
     
  
