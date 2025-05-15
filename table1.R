library(readxl)
library(dplyr)
d<-read_excel("../PrimateAgingManuscript/Primate\ Aging\ Data\ Tables\ Final.xlsx", sheet=8)
d$beta_s <- log(2)/d$mrdr
d$alpha_s <- exp(d$alpha)
d$alpha_s_up <- exp(d$alpha+d$alpha_sd)
d$alpha_s_low <- exp(d$alpha-d$alpha_sd)

d$lifepan_q99 <- qgomp(0.99,shape=d$beta_s,rate=d$alpha_s)
d$lifepan_q99l <- qgomp(0.99,shape=d$beta_s,rate=d$alpha_s_up)
d$lifepan_q99u <- qgomp(0.99,shape=d$beta_s,rate=d$alpha_s_low)


d_filtered <- d %>% filter(node %in% c("46n","30n","23n","14n","19n")) %>%
  select(node, alpha, alpha_sd, mrdr, mrdr_stddev, lifepan_q99, lifepan_q99l, lifepan_q99u) %>% arrange(mrdr)

xtable::xtable(d_filtered, digits = c(0, 0, 2, 2, 2, 2, 2, 2, 2), 
         label = "tab:table1", 
         caption = "Table 1. Lifespan estimates for the primate species in the study. The table includes the species name (node), the alpha and beta parameters of the Gompertz model, the mean relative risk of death (MRDR), and the estimated lifespan at the 99th percentile (lifepan_q99) with its lower and upper bounds.") %>%
  print(include.rownames = FALSE, floating = TRUE)
