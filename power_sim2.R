rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
NumMCMCChains=4

##setup model data frame
primate_tree_subset = keep.tip(primate_tree,tip = unique(ALL$species))
dd = ALL %>% filter(species %in% primate_tree_subset$tip.label) %>%  filter(lifespan>=sexMature)
dd$lifespan = dd$lifespan - dd$sexMature #adjust t0 for sexual maturity

tree_oder=data.frame(species=paste(primate_tree_subset$tip.label), tree_order=1:(primate_tree_subset$Nnode+1))
dd=dd %>% left_join(tree_oder)

imalanced=dd %>% group_by(species) %>% 
  summarise(n=n(), nFemale=sum(sex=="Female"),  nMale=sum(sex=="Male")) %>% 
  filter(nFemale/nMale < 0.5 | n < 100)

## can't fit..assumes sex==0 
for(sp in unique(imalanced$species)) {
  print(sp);
  dd$sex[dd$species == sp] <- "Female"
}

# stan fits
stanFitSurvivalCurveLogLik =function(inputdf, Sigma, sample_n=100) {
  
  modeldf = inputdf %>% filter(species != "Rhesus Indian-derived")
  RH = inputdf %>%
    filter(species=="Rhesus Indian-derived") %>% 
    sample_n(n, replace=FALSE)
    
  modeldf = bind_rows(modeldf, RH);

  model_data = list(
    N=nrow(modeldf),
    K=length(unique(modeldf$species)),
    species = modeldf$tree_order,
    event=modeldf$event,
    time=modeldf$lifespan,
    sex=ifelse(modeldf$sex == "Female",0,1),
    VCV=Sigma
  );
  
  #starting parameters for alpha and beta
  initf = function() { list("alpha0"=-3, "beta0"=-1) }
  
  #run MCMC
  model <- stan_model(file="gompertzFitLogLikUnion.stan");
  fit <- rstan::sampling(model,
                         data = model_data,  
                         chains = NumMCMCChains, 
                         iter = 500,
                         init = initf
  );
  return(fit)
}



#FIT BAYESIAN MODEL

if (ESTIMATE) {
  #Fit STAN model, ( diag VCV matrix, ie no information sharing)
  VCV=vcv(primate_tree_subset) #BROWNIAN 
  nspp=nrow(VCV);
  VCV <- VCV/(det(VCV)^(1/nspp)); #rescale vcv
  
  for (n in c(10,20,30,50,100,200,500,1000)) {
    stanfitBrownian  = stanFitSurvivalCurveLogLik(dd, VCV, sample_n=n)
    saveRDS(stanfitBrownian,sprintf("outputs/gompertzFitLogLikUnion.brownian.fit.power_%d.rds",n))
  }
}

#load results 
allfits = data.frame()
  for( n in c(10,20,30,50,100,200,500,1000)) {
    fit=readRDS(sprintf("outputs/gompertzFitLogLikUnion.brownian.fit.power_%d.rds",n))
    fit = getStanFit(fit) %>% filter(species == "Rhesus Indian-derived")
    fit$n = n
    allfits = bind_rows(allfits, fit)
  }

#forest plot 
p0=ggplot(allfits, aes(x=n, y=beta)) +
  geom_point() +
  geom_errorbar(aes(ymin=beta - 1.96*betaStdDev, ymax=beta + 1.96*betaStdDev), width=0.2) +
  xlab("Sample Size") +
  ylab("Estimated Aging Rate (beta)") +
  ggtitle("Power Simulation: Aging Rate Estimate - Subsampling Rhesus Macaque") +
  geom_hline(yintercept=0.17, linetype="dashed", color = "red") +
  theme(base_size = 18) + theme_bw()


p1=ggplot(allfits, aes(x=n, y=logA)) +
  geom_point() +
  geom_errorbar(aes(ymin=logA - 1.96*alphaStdDev, ymax=logA + 1.96*alphaStdDev), width=0.2) +
  xlab("Sample Size") +
  ylab("Estimated Baseline Hazard ln(alpha)") +
  ggtitle("Power Simulation: Baseline Hazard Estimates- Subsampling Rhesus Macaque") +
  geom_hline(yintercept=-6.1, linetype="dashed", color = "red") +
  theme(base_size = 18) + theme_bw()

pdf("plots/power_simulation_rhesus.pdf",width=7.5,height=8)
ggarrange(p0,p1, ncol=1)
dev.off()

