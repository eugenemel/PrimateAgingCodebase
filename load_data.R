#ALL=read_excel("lifespan_data/
BW=read_excel("data/subjects_measurements_quarterly_50th_BodyWeight.xlsx")
PAD<- read.csv("data/PAD_subjects.tsv",sep="\t",stringsAsFactors = F) %>% 
   dplyr::select("subject","birthdate","deathdate","status","species","sex","housing","diet",
                 "experimental", "last_measurement_date","last_transfer_date","last_observation_date")

dim(PAD)
table(PAD$status)

PAD$source = "PAD";
PAD$birthdate = as.Date(PAD$birthdate, format="%m/%d/%Y")
PAD$deathdate = as.Date(PAD$deathdate, format="%m/%d/%Y")
PAD$last_measurement_date= as.Date(PAD$last_measurement_date, format="%m/%d/%Y")
PAD$last_transfer_date=as.Date(PAD$last_transfer_date, format="%m/%d/%Y")
PAD$last_observation_date=as.Date(PAD$last_observation_date, format="%m/%d/%Y")
PAD$experimental=gsub(" ","",PAD$experimental)
PAD$experimental[ PAD$experimental=="" ] <- "UNKNOWN"


PAD$event=0
PAD$status[ PAD$status %in% c("alive","Alive") ] <- "alive"
PAD$status[ PAD$status %in% c("lost","Lost","transferred", "Transferred")]<- "lost"
PAD$status[ PAD$status %in% c("died", "Died", "Dead","death","dead")] <- "died"
PAD$status[ !is.na(PAD$deathdate) ] <- "died"
PAD = PAD %>% filter(!status == "")
PAD$status = paste(PAD$status)
table(PAD$status,PAD$experimental)

#Animals that died
PAD$lifespan=0
tmpset = which(PAD$status == "died")
PAD$lifespan[tmpset] =as.numeric((PAD$deathdate[tmpset]-PAD$birthdate[tmpset])/365.25)
PAD$event[tmpset]=1

#Animals that died within short time period
source("strange_death.R")
tmpset = which(PAD$strange_death == "YES")
PAD$event[tmpset]=0

#Animals that are alive or lost, set event = 0
tmpset = which(PAD$status %in% c("alive","lost"))
PAD$lifespan[tmpset] =as.numeric(PAD$last_observation_date[tmpset]-PAD$birthdate[tmpset])/365.25
PAD$event[tmpset]=0

#Animals that are manually marked as experimental endpoint 
tmpset = which(PAD$experimental == "YES")
PAD$lifespan[tmpset] =as.numeric(PAD$last_observation_date[tmpset]-PAD$birthdate[tmpset])/365.25
PAD$event[tmpset]=0

LAST_OBS=BW %>% dplyr::select(subject,starts_with("y")) %>% reshape2::melt(id="subject") %>% 
  filter(!is.na(value)) %>% mutate(variable=as.numeric(gsub("^y","",variable))) %>% 
  arrange(subject,variable) %>% group_by(subject) %>% 
  summarise(last_measurment_bw=max(variable), max_bw=max(value))

ALL= PAD %>% left_join(LAST_OBS, by="subject")
dim(ALL)
ALL$sex=factor(paste(ALL$sex),level=c("Female","Male"))

### infint mortaltaity
#ALL = ALL %>% filter(lifespan >0) #filter out strange lifespans

ALL = ALL %>% dplyr::select(subject,birthdate, deathdate, lifespan, species, sex, housing, diet, event,max_bw, experimental,strange_death,status,"source")
table(ALL$species, completeLifespan=!is.na(ALL$deathdate) & !is.na(ALL$birthdate))
ALL$status[ALL$experimental == "YES" | ALL$strange_death == "YES"] <- "died_experimental"
table(ALL$species,ALL$status)

#Merge strains
#ALL=ALL %>% filter(species != "Rhesus Chinese-derived");
ALL$species[ ALL$species == "Blue-eyed Black Lemur" ] <- "Common Black Lemur";
ALL$species[ ALL$species == "Rhesus Chinese-derived" ] <- "Rhesus Indian-derived";
ALL$species[ ALL$species == "Golden-Headed Lion Tamarin" ] <- "Golden Lion Tamarin";

################################
### LOAD STUD BOOK DATA
############################### 
STUDBOOKS=read.csv("data/STUD_DATA_20220204.tsv",sep="\t")
STUDBOOKS$birthdate = as.Date(STUDBOOKS$birthdate)
STUDBOOKS$deathdate = as.Date(STUDBOOKS$deathdate)

ALL=bind_rows(ALL,STUDBOOKS);
#Remove suspect species
SUSPECT = c("Canine", "Human", "Feline", "Aye-aye")

ALL = ALL %>% filter(!(species %in% SUSPECT)) 

### remove lifespan less than half a year
ALL = ALL %>% filter(lifespan >= 1 & lifespan <= 100)
ALL$sex=factor(paste(ALL$sex),level=c("Female","Male"))
table(ALL$species,ALL$source)
dim(ALL)

#summary of species 
SPECIES_SUMMARY=ALL %>% group_by(species) %>% summarise(
  N=n(), 
  Ndeath=sum(event), 
  maxLifespan=max(lifespan,na.rm=TRUE),
  Nfemales=sum(ifelse(sex=="Female",1,0)),
  Nmales=sum(ifelse(sex=="Male",1,0)),
  medianBodyWt=median(max_bw,na.rm=TRUE)
)

SPECIES_SUMMARY= SPECIES_SUMMARY %>% 
  left_join(ALL %>% filter(sex=="Female") %>% 
   group_by(species) %>% 
   summarise(medianBodyWtFemale=median(max_bw,na.rm=TRUE)))

SPECIES_SUMMARY = SPECIES_SUMMARY %>% 
  left_join(ALL %>% filter(sex=="Male") %>% 
   group_by(species) %>% 
   summarise(medianBodyWtMale=median(max_bw,na.rm=TRUE)))


SPECIES_SUMMARY$maxLifespan [ is.infinite(SPECIES_SUMMARY$maxLifespan) ] = NA

kmfit=survfit(Surv(lifespan,event)~species,data=ALL)
medianKMSurv=as.data.frame(quantile(kmfit,probs = c(0.5,0.95))) %>%  dplyr::select(quantile.50,quantile.95)
medianKMSurv$species = gsub("species=","",rownames(medianKMSurv))

kmfitsex=survfit(Surv(lifespan,event)~species+sex,data=ALL)
medianSexSurv=as.data.frame(quantile(kmfitsex,probs = c(0.5,0.95))) %>%  dplyr::select(quantile.50,quantile.95)
medianSexSurv$species = gsub("species=","",rownames(medianSexSurv))
kmfitsexMale=medianSexSurv %>% dplyr::filter(grepl("Male",species))
kmfitsexFemale=medianSexSurv %>% dplyr::filter(grepl("Female",species))
kmfitsexMale$species = gsub(", sex=.*","",kmfitsexMale$species)
kmfitsexFemale$species = gsub(", sex=.*","",kmfitsexFemale$species)

SPECIES_SUMMARY= SPECIES_SUMMARY %>% left_join(medianKMSurv,by="species", suffix = c("","both")) %>% 
                                    left_join(kmfitsexMale,by="species",suffix = c("","Male")) %>% 
                                    left_join(kmfitsexFemale, by="species",suffix = c("","Female"))

## remove species with too few events
LOWCOUNTS=SPECIES_SUMMARY$species[ SPECIES_SUMMARY$Ndeath < 5]
ALL = ALL %>% filter(!species %in% LOWCOUNTS); 
SPECIES_SUMMARY  = SPECIES_SUMMARY %>% filter(!species %in% LOWCOUNTS)



#MAP species to evalutionary tree
common_sci_names <-read.csv("data/CommonNameMapping.tsv",sep="\t",stringsAsFactors = FALSE)
common_sci_names$tip.label = gsub(" ", "_", common_sci_names$sci_name)


primates <- read.nexus("data/10kTrees_consensusTree_version3.nex",force.multi = FALSE)
primates = primates[[1]]
primate_tree = keep.tip(primates,common_sci_names$tip.label)

#order table by tree labels
common_sci_names  = common_sci_names[match(primate_tree$tip.label,common_sci_names$tip.label), ]
primate_tree$tip.label = paste(common_sci_names$common_name)

#clean up
KEEP_OBJECTS=c("ALL", "SPECIES_SUMMARY","KEEP_OBJECTS","primate_tree","common_sci_names","CenterData")
remove(list=c(setdiff(ls(), KEEP_OBJECTS)))

setdiff(common_sci_names$common_name, unique(ALL$species))
setdiff(SPECIES_SUMMARY$species,common_sci_names$common_name)
SPECIES_SUMMARY=SPECIES_SUMMARY %>% left_join(common_sci_names,by=c("species" = "common_name"))

## sexual maturity
ALL$sexMature=0;
for (i in 1:nrow(SPECIES_SUMMARY)) {
  ALL$sexMature[ ALL$species == SPECIES_SUMMARY$species[i]] =  SPECIES_SUMMARY$SexMature[i];
}
### save final dataframe into tsv
write.csv(SPECIES_SUMMARY,"outputs/primates_table_ndeath.csv")
ALL %>% group_by(species) %>% arrange(desc(lifespan)) %>%  readr::write_csv("outputs/ALL_SUBJECTS_FINAL.tsv")
