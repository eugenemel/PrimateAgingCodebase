source("support_functions.R")

primate_tree_subset = keep.tip(primate_tree,tip = unique(ALL$species))
#Life-table construction only post Sex Maturity
dd = ALL %>% filter(species %in% primate_tree_subset$tip.label) %>% 
  filter(lifespan>=sexMature)
dd$lifespan = dd$lifespan - dd$sexMature
levels(dd$sex) = c("Female","Male","Both")

#adjust for truncation
tree_oder=data.frame(species=paste(primate_tree_subset$tip.label), 
                     tree_order=1:(primate_tree_subset$Nnode+1))
dd=dd %>% left_join(tree_oder)
imalanced=dd %>% group_by(species) %>% 
  summarise(n=n(), 
            nFemale=sum(sex=="Female"), 
            nMale=sum(sex=="Male")) %>% filter(nFemale/nMale < 0.5 | n < 100)

## can't fit..assumes sex==0 
for(sp in unique(imalanced$species)) {
  print(sp);
  dd$sex[dd$species == sp] <- "Both"
}

#
# Lifetable precedure
#

LT <- dd %>% group_by(species,sex,source) %>% do(lifetable(.))
LT =  LT %>% filter(rate >0)
write.table(LT, file="outputs/primate_life_tables.tsv",sep="\t")

#filer out species with fewer than 5 qx estimates
EXCLUDE =LT %>% group_by(species,sex) %>% summarise(np=n()) %>% filter(np<3)
LT = LT%>% filter( !(paste(species,sex) %in% paste(EXCLUDE$species,EXCLUDE$sex)) )
rm(EXCLUDE)


#
# wrapper for linear model fits against lifetable
#

LMFIT=LT %>% group_by(species,source) %>%  do(lmbroom(.))
LMFIT$model = "LifeTable_lm"

RLMFIT=LT %>% group_by(species,source) %>%  do(rlmbroom(.))
RLMFIT$model = "LifeTable_rlm"

LT=LT %>% group_by(species,source) %>%  do(rlmpredict(.))
LT=LT %>% group_by(species,source) %>%  do(lmpredict(.))

lmbeta=LMFIT %>% dplyr::filter(term=="time") %>% dplyr::select(species,beta=estimate)
lmalpha=LMFIT %>% dplyr::filter(term=="(Intercept)") %>% dplyr::select(species,alpha=estimate)
lmbeta=lmbeta %>% group_by(species) %>% summarise(beta=mean(beta,na.rm=TRUE))
lmalpha=lmalpha %>% group_by(species) %>% summarise(alpha=mean(alpha,na.rm=TRUE))

#
# FLEXSURV
#
#install.packages("flexsurv")
#species   term  estimate std.error model
# library("flexsurv")
# FLEXFITS=data.frame();
# for( cohortName in unique(dd$species)) {
#   print(cohortName)
#   Z = dd %>% dplyr::filter(species == cohortName) %>% arrange(lifespan)
#   Z$sex  = ifelse(Z$sex == "Female", 0, 1)
#   if(length(unique(Z$sex)) == 1) { Z$sex=0 };
#   gfit=flexsurvreg(Surv(lifespan,event) ~ sex, dist="gompertz", data=Z, 
#                    control=list(fnscale = 2500, method="L-BFGS-B")); 
#   gfit_b=gfit$res.t["shape",]
#   gfit_a=gfit$res.t["rate",]
#   df=data.frame(species=cohortName, term="alpha0", estimate=gfit_a[1], std.error=gfit_a[4], model="Direct_FlexSurv");
#   FLEXFITS = bind_rows(FLEXFITS,df)
#   df=data.frame(species=cohortName, term="beta0", estimate=gfit_b[1], std.error=gfit_b[4], model="Direct_FlexSurv");
#   FLEXFITS = bind_rows(FLEXFITS,df)
# }
#
#clean up

SPECIES_SUMMARY$adultMaxLifespan=SPECIES_SUMMARY$maxLifespan - SPECIES_SUMMARY$SexMature
SPECIES_SUMMARY$adultMedLifespan=SPECIES_SUMMARY$quantile.50 - SPECIES_SUMMARY$SexMature

lmfit=lm(adultMaxLifespan ~ SexMature,data=SPECIES_SUMMARY); summary(lmfit)
p0=ggplot(SPECIES_SUMMARY,aes(y=adultMaxLifespan, x=SexMature, label=species, color=family_desc)) + 
  geom_point(size=5,alpha=0.9) + theme_bw(base_size = 16) + 
  geom_abline(slope = lmfit$coefficients[2], intercept =  lmfit$coefficients[1]) + 
  #geom_label_repel(max.overlaps = 9) + scale_color_brewer(palette="Set1") + 
  geom_abline(slope = 1, intercept=0, linetype = "dashed") + 
  ylab("Adult Max Lifespan") + xlab("Age of Sexual Maturity")

lmfit=lm(adultMedLifespan ~ SexMature,data=SPECIES_SUMMARY); summary(lmfit)
p1=ggplot(SPECIES_SUMMARY,aes(y=adultMedLifespan, x=SexMature, label=species, color=family_desc)) + 
  geom_point(size=5,alpha=0.9) + theme_bw(base_size = 16) + 
  geom_abline(slope = lmfit$coefficients[2], intercept =  lmfit$coefficients[1]) + 
#  geom_label_repel(max.overlaps = 9) + scale_color_brewer(palette="Set1") + 
  geom_abline(slope = 1, intercept=0, linetype = "dashed") + 
  ylab("Adult Median Lifespan") + xlab("Age of Sexual Maturity")

lmfit=lm(adultMaxLifespan~adultMedLifespan,data=SPECIES_SUMMARY); summary(lmfit)
p2=ggplot(SPECIES_SUMMARY,aes(x=adultMedLifespan, y=adultMaxLifespan, label=species,color=family_desc)) + 
  geom_point(size=5,alpha=0.9) + theme_bw(base_size = 16) + 
  geom_abline(slope = lmfit$coefficients[2], intercept =  lmfit$coefficients[1]) + 
  geom_abline(slope = 1, intercept=0, linetype = "dashed") + 
  geom_label_repel(max.overlaps = 8) + scale_color_brewer(palette="Set1") + 
  ylab("Adult Max Lifespan") + xlab("Adult Median Lifespan")


pdf("plots/lifespan_agesex_mature.pdf",width=12,height=6); ggarrange(p0,p1); dev.off();

##cleanup
KEEP_OBJECTS=c(KEEP_OBJECTS, "lifetable", "LT",  "LMFIT", "RLMFIT")
remove(list=c(setdiff(ls(), KEEP_OBJECTS)))#


