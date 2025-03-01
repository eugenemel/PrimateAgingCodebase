## Censorship of experimental end points

# Clean up site names and extract year-week
PAD$site <- gsub("\\-.*", "", PAD$subject)
PAD$death_yearweak <- strftime(PAD$deathdate, format = "%Y:%V")
PAD$birth_yearweak <- strftime(PAD$birthdate, format = "%Y:%V")

# Identify deaths that occurred in patterns (same week or same day)
SAMEWEEK <- PAD %>% 
    filter(event == 1) %>% 
    group_by(site, species, death_yearweak) %>% 
    mutate(week_event = n()) %>% 
    ungroup() %>% 
    filter(week_event >= 3)

SAMEDAY <- PAD %>% 
    filter(event == 1) %>% 
    group_by(site, species, deathdate) %>% 
    mutate(day_event = n()) %>% 
    ungroup() %>% 
    filter(day_event >= 2)

# Display summary tables
table(SAMEDAY$day_event)
table(SAMEWEEK$week_event)

# Flag strange deaths in the main dataset
PAD$strange_death <- "NO"
PAD$strange_death[PAD$subject %in% c(SAMEWEEK$subject, SAMEDAY$subject)] <- "YES"

# Export results
write.csv(SAMEWEEK, file = "outputs/strange_death_week.csv", row.names = FALSE)
write.csv(SAMEDAY, file = "outputs/strange_death_day.csv", row.names = FALSE)
