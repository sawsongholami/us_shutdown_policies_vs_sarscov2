# Authors: Sawson Gholami, Grant Pointon
# Created: 06/12/2020
# Updated: 08/25/2020
# Description: A (too-crowded) script for importing, cleaning, and exploring data on states' responses to COVID-19


# Importing libraries -----------------------------------------------------

library(data.table)
library(janitor)
library(lubridate)
library(tidyverse)

library(cluster)
library(stats)

library(lme4)
library(lmerTest)

library(cowplot)
library(factoextra)
library(ggrepel)


# Custom functions --------------------------------------------------------

## A useful little function
`%ni%` <- Negate(`%in%`)


## This function finds the sum of all cases/deaths at date_target for each state.
## If the user wishes to find the sum of all cases/deaths at the shutdown or opening dates, which vary by state, then
##   then they should use date_special instead of date_target.
counts_at_date <- function(states_all, cases_or_deaths, date_target = NULL, date_special = FALSE, date_increment = 0, counts_long. = counts_long, counts_wide. = counts_wide) {
  
  result <- vector("list", length(states_all))
  
  for(idx in seq_along(states_all)){
    
    state_target <- states_all[[idx]]
    
    if (date_special == "shutdown") {
      date_target <- counts_wide. %>% filter(state == state_target) %>% select(shutdown_date_official) %>% .[[1]] %>% as_date()
    } else if (date_special == "opening") {
      date_target <- counts_wide. %>% filter(state == state_target) %>% select(shutdown_date_official) %>% .[[1]] %>% as_date()
    } else {
      date_target <- as_date(date_target)
    }
    
    date_target <- date_target + date_increment
    date_increment <- 0 ### This ensures that date_increment is added to date_target only on the first iteration of the for-loop.
    
    date_initial <- as_date(counts_long.$date[[1]])
    
    idx_col <- min(which(str_detect(colnames(counts_wide.), cases_or_deaths)))
    
    result[[idx]] <- counts_wide. %>%
      filter(state == state_target) %>%
      select(idx_col + abs(as.integer(date_target - date_initial)))
    
  }
  
  return(unlist(result))
  
}




# Importing data ----------------------------------------------------------

## Names and abbreviations for all states
state_names_raw <- read_csv("data/geo_crosswalk.csv")

## Population info for all states
state_pop_raw <- read_csv("data/state_populations.csv")

## Policy data
policies_raw <- read_csv("data/state_responses.csv")

## Daily counts of cases and deaths from NYT: https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv
counts_raw <- read_csv("data/us-states.csv")

## Daily estimates of basic reproduction numbers (R0) from rt.live: https://d14wlfuexuxgcm.cloudfront.net/covid/rt.csv
r0_raw <- read_csv("data/rt.csv")

## Daily mobility from Google
mobility_raw <- read_csv("data/Global_Mobility_Report.csv")




# Cleaning data -----------------------------------------------------------

state_names <- state_names_raw %>% 
  select(state = state_name, state_abbr, state_code) %>%
  distinct() %>% 
  arrange(state)


state_pop <- state_pop_raw %>% 
  group_by(state) %>% 
  summarise(state_pop = sum(tot_pop)) %>% 
  ungroup()


policies <- policies_raw %>%
  clean_names() %>% 
  filter(state != "Washington D.C." & state != "District of Columbia") %>% 
  left_join(state_names, by = "state") %>%
  left_join(state_pop, by = "state") %>%
  select(state,
         state_abbr,
         state_pop,
         governor_gender,
         shutdown_date_official = date_of_first_shutdown,
         shutdown_ordertype = type_of_order,
         shutdown_level = level_of_first_shutdown,
         shutdown_gathlim = gathering_limit_under_lockdown,
         shutdown_gathlim_ordertype = guideline_or_requirement,
         shutdown_details = details_of_first_shutdown,
         opening_date_official = date_of_start_of_first_opening,
         opening_date_initial = initial_expiration_date,
         opening_level = level_of_first_opening,
         opening_gathlim = gathering_limit_at_reopening,
         asof0605_gathlim = gathering_limit_6_5) %>%
  ### Replacing the missing initial opening dates with the official opening dates
  mutate(opening_date_initial = ifelse(is.na(opening_date_initial), opening_date_official, opening_date_initial)) %>% 
  ### Correcting the type of all columns with dates
  mutate(shutdown_date_official = mdy(shutdown_date_official),
         opening_date_official = mdy(opening_date_official),
         opening_date_initial = mdy(opening_date_initial)) %>%
  ### Converting all lists to doubles 
  mutate(opening_gathlim = as.double(unlist(opening_gathlim)),
         asof0605_gathlim = as.double(unlist(asof0605_gathlim))) %>%
  ### Converting all text to lowercase
  mutate(shutdown_level = str_to_lower(shutdown_level),
         shutdown_ordertype = str_to_lower(shutdown_ordertype),
         shutdown_gathlim_ordertype = str_to_lower(shutdown_gathlim_ordertype),
         opening_level = str_to_lower(opening_level)) %>% 
  ### Replacing the abbreviations in shutdown_gathlim_ordertype with more explicit values
  mutate(shutdown_gathlim_ordertype = case_when(shutdown_gathlim_ordertype == "g" ~ "guideline",
                                                shutdown_gathlim_ordertype == "r" ~ "requirement")) %>%
  ### Replacing unspecified gathering limits with Inf
  mutate(shutdown_gathlim = ifelse(is.na(shutdown_gathlim), Inf, shutdown_gathlim),
         opening_gathlim = ifelse(is.na(opening_gathlim), Inf, opening_gathlim),
         asof0605_gathlim = ifelse(is.na(asof0605_gathlim), Inf, asof0605_gathlim)) %>%
  ### Creating a new column to identify the states whose legislatures challenged their governor's executive order
  mutate(legislature_challenged_governor = case_when(state_abbr %in% c("NC", "KS", "WI") ~ TRUE,
                                                     TRUE ~ FALSE),
         .after = governor_gender) %>%
  arrange(state)


counts_long <- tibble(date = seq(counts_raw$date %>% min(), counts_raw$date %>% max(), by = "days")) %>%
  uncount(length(policies$state)) %>%
  group_by(date) %>% 
  mutate(state = policies$state) %>%
  left_join(counts_raw %>% select(-fips), by = c("date", "state")) %>% 
  mutate(cases = ifelse(is.na(cases), 0, cases),
         deaths = ifelse(is.na(deaths), 0, deaths)) %>% 
  ungroup()


counts_wide <- counts_raw %>%
  filter(state %in% policies$state) %>% 
  left_join(policies, by = "state") %>% 
  select(date, state, state_abbr, shutdown_date_official, opening_date_official, cases, deaths) %>%
  pivot_wider(names_from = date, values_from = c(cases, deaths), values_fill = 0) %>%
  clean_names() %>%
  arrange(state)


r0_long <- r0_raw %>%
  clean_names() %>% 
  mutate(date = as_date(date)) %>% 
  select(date, state_abbr = region, r0 = mean, median, lower_80, upper_80) %>% 
  filter(!is.na(r0)) %>% 
  filter(state_abbr %in% policies$state_abbr)



mobility_all <- mobility_raw %>%
  clean_names() %>%
  rename(state = sub_region_1) %>% 
  filter(country_region == "United States" & !is.na(state)) %>% 
  filter(state %in% policies$state) %>% 
  mutate(date = mdy(date)) %>% 
  select(state, date:13) %>% 
  group_by(state, date) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
  ungroup()


mobility_weekdays <- mobility_all %>%
  mutate(day_of_week = weekdays(date)) %>% 
  filter(day_of_week != "Saturday" & day_of_week != "Sunday") %>%
  rename(retail_rec_mobility_weekdays = retail_and_recreation_percent_change_from_baseline,
         workplaces_mobility_weekdays = workplaces_percent_change_from_baseline,
         residential_mobility_weekdays = residential_percent_change_from_baseline) %>% 
  group_by(state) %>% 
  mutate(retail_rec_mobility_rolling_avg = frollmean(retail_rec_mobility_weekdays, n = 5, align = "left", na.rm = TRUE),
         workplaces_mobility_rolling_avg = frollmean(workplaces_mobility_weekdays, n = 5, align = "left", na.rm = TRUE),
         residential_mobility_rolling_avg = frollmean(residential_mobility_weekdays, n = 5, align = "left", na.rm = TRUE)) %>%
  arrange(state, date) %>% 
  ungroup() %>% 
  select(-c(day_of_week, parks_percent_change_from_baseline,transit_stations_percent_change_from_baseline, grocery_and_pharmacy_percent_change_from_baseline))


mobility_shutdown <- mobility_weekdays %>%
  group_by(state) %>%
  mutate(in_mobility_shutdown = ifelse(retail_rec_mobility_rolling_avg <= -20, TRUE, FALSE)) %>%
  filter(in_mobility_shutdown == TRUE) %>% 
  mutate(shutdown_date_mobility = first(date)) %>% 
  ungroup() %>% 
  select(state, shutdown_date_mobility) %>% 
  distinct()


## This is the fully joined dataset.
responses_all <- counts_long %>%  ## cases and deaths
  left_join(policies, by = "state" ) %>%  ## policies
  left_join(r0_long, by = c("date", "state_abbr")) %>%  ### R0
  left_join(mobility_all, by = c("date", "state")) %>%  ### mobility, raw values
  left_join(mobility_weekdays, by = c("date", "state")) %>%  ### mobility, rolling averages
  left_join(mobility_shutdown, by = "state") %>%  ### mobility-defined shutdown dates
  mutate(cases_per_capita = cases / state_pop,
         deaths_per_capita = deaths / state_pop) %>% 
  mutate(currently_in_shutdown_official = case_when(date >= shutdown_date_official & date < opening_date_official ~ TRUE,
                                                    TRUE ~ FALSE),
         currently_in_shutdown_mobility = case_when(date >= shutdown_date_mobility & date < opening_date_official ~ TRUE,
                                                    TRUE ~ FALSE)) %>%
  ### Creating a new column to represent the number of days since shutdown/opening for each state
  group_by(state) %>% 
  mutate(shutdown_day_count_official = case_when(date >= shutdown_date_official & date < opening_date_official ~ as.double(date - shutdown_date_official),
                                                 TRUE ~ -1),
         shutdown_day_count_mobility = case_when(date >= shutdown_date_mobility & date < opening_date_official ~ as.double(date - shutdown_date_mobility),
                                                 TRUE ~ -1),
         opening_day_count_official = case_when(date >= opening_date_official ~ as.double(date - opening_date_official),
                                                TRUE ~ -1),
         opening_day_count_initial = case_when(date >= opening_date_initial ~ as.double(date - opening_date_initial),
                                               TRUE ~ -1)) %>%
  ### Creating a new column to represent an approximation of the second derivatives of retail/rec. mobility and r0 values
  arrange(state, date) %>%
  mutate(retail_rec_mobility_second_deriv = lead(retail_and_recreation_percent_change_from_baseline) -2*retail_and_recreation_percent_change_from_baseline + lag(retail_and_recreation_percent_change_from_baseline),
         retail_rec_second_deriv_rolling = frollmean(retail_rec_mobility_second_deriv, n = 5, align = "left", na.rm = TRUE),
         r0_second_deriv = lead(r0) - 2*r0 + lag(r0)) %>%
  ungroup() %>%
  filter(date <= as_date(max(mobility_all$date)))


## Similar to the above, this is the fully joined dataset with weekends excluded.
responses_weekdays <- responses_all %>% 
  mutate(day_of_week = weekdays(date)) %>% 
  filter(day_of_week != "Saturday" & day_of_week != "Sunday") %>%
  select(-day_of_week) %>%
  ### Creating a new column to represent an approximation of the second derivatives of retail/rec. mobility and r0 values
  group_by(state) %>% 
  arrange(state, date) %>%
  mutate(retail_rec_mobility_second_deriv = lead(retail_and_recreation_percent_change_from_baseline) -2*retail_and_recreation_percent_change_from_baseline + lag(retail_and_recreation_percent_change_from_baseline),
         r0_second_deriv = lead(r0) - 2*r0 + lag(r0)) %>%
  ungroup() %>%
  filter(date <= as_date(max(mobility_weekdays$date)))




# EDA: Cases and deaths ---------------------------------------------------

## WARNING: The following is SLOW.
counts_at_main_intervals <- counts_wide %>%
  ### Cases for interval up to shutdown
  mutate(cases_at_shutdown = counts_at_date(state, "cases", date_special = "shutdown"),
         cases_at_shutdown_plus14days = counts_at_date(state, "cases", date_special = "shutdown", date_increment = 14),
         cases_at_opening = counts_at_date(state, "cases", date_special = "opening"),
         cases_at_opening_plus14days = counts_at_date(state, "cases", date_special = "opening", date_increment = 14),
         cases_at_0605 = counts_at_date(state, "cases", date_target = "2020-06-05"),
         cases_at_0605_plus14days = counts_at_date(state, "cases", date_target = "2020-06-05", date_increment = 14),
         cases_at_present = counts_at_date(state, "cases", date_target = tail(counts_raw$date, n=1))) %>%
  ### Deaths for interval up to shutdown
  mutate(deaths_at_shutdown = counts_at_date(state, "deaths", date_special = "shutdown"),
         deaths_at_shutdown_plus14days = counts_at_date(state, "deaths", date_special = "shutdown", date_increment = 14),
         deaths_at_opening = counts_at_date(state, "deaths", date_special = "opening"),
         deaths_at_opening_plus14days = counts_at_date(state, "deaths", date_special = "opening", date_increment = 14),
         deaths_at_0605 = counts_at_date(state, "deaths", date_target = "2020-06-05"),
         deaths_at_0605_plus14days = counts_at_date(state, "deaths", date_target = "2020-06-05", date_increment = 14),
         deaths_at_present = counts_at_date(state, "deaths", date_target = tail(counts_raw$date, n=1))) %>% 
  ### Cases for interval between shutdown and opening
  mutate(cases_gained_between_shutdown_and_opening = cases_at_opening - cases_at_shutdown,
         cases_gained_between_opening_and_0605 = cases_at_0605 - cases_at_opening,
         cases_gained_between_0605_and_present = cases_at_present - cases_at_0605) %>% 
  ### Deaths for interval between shutdown and opening
  mutate(deaths_gained_between_shutdown_and_opening = deaths_at_opening - deaths_at_shutdown,
         deaths_gained_between_opening_and_0605 = deaths_at_0605 - deaths_at_opening,
         deaths_gained_between_0605_and_present = deaths_at_present - deaths_at_0605)


## Here we find the dates when the earliest/latest shutdowns/openings occured, and the corresponding states.
date_latest_overall <- max(responses_all$date)
### Official dates
date_shutdown_earliest_official <- min(responses_all$shutdown_date_official)
date_shutdown_latest_official <- max(responses_all$shutdown_date_official)
date_opening_earliest_official <- min(responses_all$opening_date_official)
date_opening_latest_official <- max(responses_all$opening_date_official)
### Official states
state_shutdown_earliest_official <- responses_all %>% filter(shutdown_date_official == date_shutdown_earliest_official) %>% select(state_abbr) %>% distinct() %>% .[[1]] %>% str_c(collapse = ", ")
state_shutdown_latest_official <- responses_all %>% filter(shutdown_date_official == date_shutdown_latest_official) %>% select(state_abbr) %>% distinct() %>% .[[1]] %>% str_c(collapse = ", ")
state_opening_earliest_official <- responses_all %>% filter(opening_date_official == date_opening_earliest_official) %>% select(state_abbr) %>% distinct() %>% .[[1]] %>% str_c(collapse = ", ")
state_opening_latest_official <- responses_all %>% filter(opening_date_official == date_opening_latest_official) %>% select(state_abbr) %>% distinct() %>% .[[1]] %>% str_c(collapse = ", ")
### Mobility/Initial dates
date_shutdown_earliest_mobility <- min(responses_all$shutdown_date_mobility)
date_shutdown_latest_mobility <- max(responses_all$shutdown_date_mobility)
date_opening_earliest_initial <- min(responses_all$opening_date_initial)
date_opening_latest_initial <-max(responses_all$opening_date_initial)
### Mobility/Initial states
state_shutdown_earliest_mobility <- responses_all %>% filter(shutdown_date_mobility == date_shutdown_earliest_mobility) %>% select(state_abbr) %>% distinct() %>% .[[1]] %>% str_c(collapse = ", ")
state_shutdown_latest_mobility <- responses_all %>% filter(shutdown_date_mobility == date_shutdown_latest_mobility) %>% select(state_abbr) %>% distinct() %>% .[[1]] %>% str_c(collapse = ", ")
state_opening_earliest_initial <- responses_all %>% filter(opening_date_initial == date_opening_earliest_initial) %>% select(state_abbr) %>% distinct() %>% .[[1]] %>% str_c(collapse = ", ")
state_opening_latest_initial <- responses_all %>% filter(opening_date_initial == date_opening_latest_initial) %>% select(state_abbr) %>% distinct() %>% .[[1]] %>% str_c(collapse = ", ")


## Plotting cases per capita vs. time.
ggplot() +
  geom_vline(xintercept = c(date_shutdown_earliest_official, date_shutdown_latest_official, date_opening_earliest_official, date_opening_latest_official),
             linetype = "solid",
             color = "grey30",
             alpha = 0.6) +
  annotate("text",
           label = paste0("earliest shutdown (", state_shutdown_earliest_official, " on ", date_shutdown_earliest_official, ")"),
           x = (date_shutdown_earliest_official + 2),
           y = 0.012,
           angle = 90,
           size = 4) +
  annotate("text",
           label = paste0("latest shutdown (", state_shutdown_latest_official, " on ", date_shutdown_latest_official, ")"),
           x = (date_shutdown_latest_official + 2),
           y = 0.012,
           angle = 90,
           size = 4) +
  annotate("text",
           label = paste0("earliest opening (", state_opening_earliest_official, " on ", date_opening_earliest_official, ")"),
           x = (date_opening_earliest_official + 2),
           y = 0.012,
           angle = 90,
           size = 4) +
  annotate("text",
           label = paste0("latest opening (", state_opening_latest_official, " on ", date_opening_latest_official, ")"),
           x = (date_opening_latest_official + 2),
           y = 0.012,
           angle = 90,
           size = 4) +
  ### Not in official shutdown
  geom_line(aes(date, cases_per_capita, group = state, linetype = currently_in_shutdown_official),
            responses_all %>% filter(currently_in_shutdown_official == FALSE,
                                     date >= as_date("2020-03-01"),
                                     state_abbr %ni% c("NJ", "NY")),
            alpha = 0.6) +
  ### In official shutdown
  geom_line(aes(date, cases_per_capita, group = state, linetype = currently_in_shutdown_official),
            responses_all %>% filter(currently_in_shutdown_official == TRUE,
                                     date >= as_date("2020-03-01"),
                                     state_abbr %ni% c("NJ", "NY")),
            alpha = 0.9) +
  labs(title = "Total cases per capita between March 01 and June 27, 2020",
       subtitle = "(excluding NY, NJ)",
       x = "date",
       y = "cases per capita",
       linetype = "In official shutdown") +
  scale_linetype_manual(values = c("dotted", "solid")) +
  guides(linetype = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(plot.title =  element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title = element_text(size = 15, face = "plain")) +
  ggsave("figures/counts_01.png", height = 8, width = 11, units = "in")




# EDA: R0 -----------------------------------------------------------------

## Subplot for official shutdown
plot_r0_official <- ggplot() +
  geom_hline(yintercept = 1, color = "grey30") +
  geom_line(aes(date, r0, group = state, linetype = "open"),
            responses_all %>% filter(currently_in_shutdown_official == FALSE),
            alpha = 0.6) +
  geom_line(aes(date, r0, group = state, linetype = "shutdown"),
            responses_all %>% filter(currently_in_shutdown_official == TRUE),
            alpha = 1) +
  labs(title = "Official shutdown",
       x = "date",
       y = expression(transmissibility~(R[0])),
       linetype = "In shutdown") +
  scale_linetype_manual(values = c("solid", "dotted"),
                        breaks = c("shutdown", "open")) + 
  theme_classic() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        legend.position = c(0.75, 0.85),
        legend.justification = c(1, 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        legend.background = element_blank()) +
  guides(linetype = guide_legend(title.position = "top"))


## Subplot for mobility-defined shutdown
plot_r0_mobility <- ggplot() +
  geom_hline(yintercept = 1, color = "grey30") +
  geom_line(aes(date, r0, group = state, linetype = "open"),
            responses_all %>% filter(currently_in_shutdown_mobility == FALSE),
            alpha = 0.6) +
  geom_line(aes(date, r0, group = state, linetype = "shutdown"),
            responses_all %>% filter(currently_in_shutdown_mobility == TRUE),
            alpha = 1) +
  labs(title = "Mobility-defined shutdown",
       x = "date",
       y = expression(transmissibility~(R[0])),
       linetype = "In shutdown") +
  scale_linetype_manual(values = c("solid", "dotted"),
                        breaks = c("shutdown", "open")) + 
  theme_classic() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        legend.position = c(0.75, 0.85),
        legend.justification = c(1, 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        legend.background = element_blank()) +
  guides(linetype = guide_legend(title.position = "top"))

## Combining the subplots
plot_grid(plot_r0_official, plot_r0_mobility, ncol = 1) +
  ggsave("figures/brn_01.png", height = 12, width = 11, units = "in")  




# EDA: Length of shutdown -------------------------------------------------

r0_across_shutdown <- responses_all %>% 
  group_by(state_abbr) %>% 
  mutate(cases_per_capita_at_opening_official = ifelse(date == opening_date_official, cases_per_capita, NA_real_),
         r0_at_shutdown_official = ifelse(date == shutdown_date_official, r0, NA_real_),
         r0_at_opening_official = ifelse(date == opening_date_official, r0, NA_real_),
         shutdown_length_official = max(shutdown_day_count_official)) %>%
  summarise(state_abbr = first(state_abbr),
            cases_per_capita_at_opening_official = max(cases_per_capita_at_opening_official, na.rm = TRUE),
            r0_at_shutdown_official = max(r0_at_shutdown_official, na.rm = TRUE),
            r0_at_opening_official = max(r0_at_opening_official, na.rm = TRUE),
            r0_difference_across_shutdown_official = r0_at_opening_official - r0_at_shutdown_official,
            shutdown_length_official = max(shutdown_length_official, na.rm = TRUE))



## Length of shutdown vs. Cases at opening
kendall_shutdown_cases <- cor.test(r0_across_shutdown$shutdown_length_official,
                                   r0_across_shutdown$cases_per_capita_at_opening_official,
                                   method = "kendall")

ggplot(r0_across_shutdown, aes(shutdown_length_official, cases_per_capita_at_opening_official)) +
  geom_text(aes(label = state_abbr),
            position = position_dodge(0.5),
            size = 4) +
  geom_smooth(formula = y ~ x,
              method = "lm",) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 14)) +
  labs(title = paste0("Official shutdown length vs. Cases at opening\n(Kendall's tau = ",
                      signif(kendall_shutdown_cases$estimate, 3),
                      ", p = ",
                      signif(kendall_shutdown_cases$p.value, 3),
                      ")"),
       x = "length of official shutdown period (days)",
       y = "cases per capita at official opening") +
  ggsave("figures/shutdown_length_cases.png", height = 8, width = 10, units = "in")



## Length of shutdown vs. Change in R0
kendall_shutdown_r0 <- cor.test(r0_across_shutdown$shutdown_length_official,
                                r0_across_shutdown$r0_difference_across_shutdown_official,
                                method = "kendall")

ggplot(r0_across_shutdown, aes(shutdown_length_official, r0_difference_across_shutdown_official)) +
  geom_text(aes(label = state_abbr),
            position = position_dodge(0.5),
            size = 4) +
  geom_smooth(formula = y ~ x,
              method = "lm",
              color = "green4") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 14)) +
  labs(title = paste0("Official shutdown length vs. Change in R0\n(Kendall's tau = ",
                      signif(kendall_shutdown_r0$estimate, 3),
                      ", p = ",
                      signif(kendall_shutdown_r0$p.value, 3),
                      ")"),
       x = "length of official shutdown period (days)",
       y = expression((R[0]~at~official~opening)-(R[0]~at~official~shutdown))) +
  ggsave("figures/shutdown_length_r0.png", height = 8, width = 10, units = "in")




# EDA: Mobility in retail and recreation ----------------------------------

mobility_weekdays_national_avg <- mobility_weekdays %>% 
  group_by(date) %>% 
  summarise(retail_rec_mobility_rolling_avg = mean(retail_rec_mobility_rolling_avg))


ggplot() +
  geom_vline(xintercept = c(date_shutdown_earliest_mobility, date_shutdown_earliest_official),
             alpha = 0.4) +
  annotate("text",
           label = paste0("earliest MOBILITY shutdown\n(", date_shutdown_earliest_mobility, ")"),
           x = date_shutdown_earliest_mobility - 3.5,
           y = -40,
           angle = 90,
           size = 2) +
  annotate("text",
           label = paste0("earliest OFFICIAL shutdown\n(", date_shutdown_earliest_official, ")"),
           x = date_shutdown_earliest_official + 2.5,
           y = 5,
           angle = 90,
           size = 2) +
  ### All states
  geom_line(aes(x = date, y = retail_rec_mobility_rolling_avg, group = state, linetype = "foo"),
            mobility_weekdays,
            linetype = "dotted",
            alpha = 0.4) +
  ### National average
  geom_line(aes(x = date, y = retail_rec_mobility_rolling_avg, linetype = "national\naverage"),
            mobility_weekdays_national_avg,
            size = 1.5) +
  theme_classic() +
  theme(plot.title =  element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title = element_text(size = 13, face = "plain"),
        legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.8, 0.97)) +
  labs(title = "Mobility in areas of retail and recreation",
       subtitle = "(weekdays only)",
       x = "date",
       y = "% change from baseline mobility\n(rolling averages with window of +4 days)") +
  ggsave("figures/mobility_retail_rec.png", height = 5, width = 10)




# EDA: Mobility second derivatives ----------------------------------------

## Plotting second derivative of mobility vs. time
ggplot() +
  geom_line(aes(date, retail_rec_second_deriv_rolling, group = state, color = state),
            responses_all %>% filter(date >= min(mobility_all$date))) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 14),
        legend.position = "none") +
  labs(title = "\"Acceleration\" of mobility in areas of recreation and retail",
       y = "approx. second derivative of mobility ") +
  ggsave("figures/mobility_second_deriv.png", height = 6, width = 10, units = "in")



## Plotting R0 vs. time, colored by second derivative of mobility
ggplot() +
  geom_line(aes(date, r0, group = state, color = retail_rec_mobility_second_deriv >= 0),
            responses_all %>% filter(!is.na(retail_rec_mobility_second_deriv))) + ### Filter here to see individual states.
  geom_hline(yintercept = 1, size = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        legend.position = c(0.65, 0.75),
        legend.text = element_text(size = 10)) +
  labs(y = expression(transmissibility~(R[0])),
       color = "\"acceleration\" of mobility (retail and recreation)") +
  guides(color = guide_legend(reverse = TRUE)) +
  scale_color_manual(values = c("TRUE" = "orangered1", "FALSE" = "black"), labels = c("TRUE" = "positive", "FALSE" = "negative")) +
  ggsave("figures/mobility_second_deriv_r0.png", height = 6, width = 11, units = "in")




# PCA: R0 -----------------------------------------------------------------

## Widening R0 data in order to perform PCA
r0_wide <- r0_long %>% 
  select(date, state_abbr, r0) %>% 
  group_by(date) %>% 
  arrange(state_abbr) %>% 
  ungroup() %>% 
  pivot_wider(names_from = state_abbr, values_from = r0) %>% 
  filter(date >= "2020-03-08", ### This is the first date when all states have a non-missing R0 value.
         date >= date_shutdown_earliest_mobility,
         date < date_opening_latest_official) %>% 
  arrange(date) %>% 
  column_to_rownames(var = "date")


## PCA
r0_pca <- prcomp(t(r0_wide), scale = TRUE)
r0_pca_as_tibble <- r0_pca$x %>%
  as_tibble() %>% 
  add_column(state_abbr = r0_pca$x %>% row.names(), .before = 1)


## PCA scree plot
r0_pca_variance <- r0_pca$sdev^2
r0_pca_variance_percent <- round(r0_pca_variance/sum(r0_pca_variance)*100, 1)


## K-means clustering and silhouette analysis

for (idx in seq(2, 7)) {
  
  r0_kmeans <- kmeans(t(r0_wide), idx)
  
  fviz_silhouette(silhouette(r0_kmeans$cluster, dist(t(r0_wide)))) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_blank(),
          panel.border = element_rect(size = 0.7, fill = NA)) +
    labs(title = paste0("Silhouette analysis with K = ", idx),
         subtitle = "K-means clustering of R0 across mobility-defined shutdown period",
         x = "",
         y = "") +
    ggsave(paste0("figures/kmeans_silhouette_0", idx, ".png"), height = 8, width = 10, units = "in")
  
}


r0_kmeans <- kmeans(t(r0_wide), 3)
r0_kmeans_as_tibble <- r0_kmeans$cluster %>%
  as_tibble() %>%
  add_column(state_abbr = r0_kmeans$cluster %>% names(), .before = 1) %>% 
  rename(cluster_kmeans = value)


## Plotting k-means clusters along PCs
r0_kmeans_for_plotting <- left_join(r0_pca_as_tibble, r0_kmeans_as_tibble)

ggplot(r0_kmeans_for_plotting, aes(PC1, PC2, color = factor(cluster_kmeans))) +
  geom_text(aes(label = state_abbr),
            hjust = 0,
            vjust = 0,
            size = 3.5,
            key_glyph = "point") +
  scale_color_manual(values = c("tomato", "purple2", "forest green", "goldenrod4", "forest green")) +
  labs(title = expression(`K-means`~clustering~of~R[0]~across~`mobility-defined`~shutdown~period),
       x = paste0("Principal Component 1 ", "(", r0_pca_variance_percent[1], "%)"),
       y = paste0("Principal Component 2 ", "(", r0_pca_variance_percent[2], "%)"),
       color = "cluster") +
  theme_classic() +
  theme(title = element_text(size = 13),
        axis.title = element_text(size = 13),
        legend.title = element_text(size = 11),
        panel.border = element_rect(size = 0.7, fill = NA),
        axis.line.x.bottom = element_line(size = 0),
        axis.line.y.left = element_line(size = 0)) +
  ggsave("figures/kmeans_clusters.png", height = 8, width = 11, units = "in")




# Modeling: R0 vs. time, with official shutdown ---------------------------

responses_for_modeling_official <- responses_all %>% 
  group_by(state) %>% 
  filter(date > shutdown_date_official & date < opening_date_official) %>% 
  ungroup() %>%
  arrange(date) %>% 
  rename(shutdown_day_count = shutdown_day_count_official)


## Mixed model linear, random intercept only
r0_lin_base_official <- lmer(r0 ~ shutdown_day_count +
                               (1 | state), 
                             data = responses_for_modeling_official,
                             REML = TRUE, 
                             lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))


## Mixed model linear
r0_lin_lmer_official <- lmer(r0 ~ shutdown_day_count +
                               (1 + shutdown_day_count | state), 
                             data = responses_for_modeling_official,
                             REML = TRUE, 
                             lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))


## Mixed model polynomial
r0_poly_lmer_official <- lmer(r0 ~ shutdown_day_count + I(shutdown_day_count^2) + 
                                (1 + shutdown_day_count | state), 
                              data = responses_for_modeling_official,
                              REML = TRUE, 
                              lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))


## Mixed model cubic
r0_cubic_lmer_official <- lmer(r0 ~ shutdown_day_count + I(shutdown_day_count^2)  + I(shutdown_day_count^3) +
                                 (1 + shutdown_day_count | state), 
                               data = responses_for_modeling_official,
                               REML = TRUE, 
                               lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=4e5)))


anova(r0_lin_base_official, r0_lin_lmer_official, r0_poly_lmer_official, r0_cubic_lmer_official)

## Cubic model shows best fit 
summary(r0_cubic_lmer_official)




# Modeling: R0 vs. time, with mobility-defined shutdown -------------------

responses_for_modeling_mobility <- responses_all %>% 
  group_by(state) %>% 
  filter(date > shutdown_date_mobility & date < opening_date_official) %>% 
  ungroup() %>%
  arrange(date) %>% 
  rename(shutdown_day_count = shutdown_day_count_mobility)


## Mixed model linear, random intercept only
r0_lin_base_mobility <- lmer(r0 ~ shutdown_day_count +
                               (1 | state), 
                             data = responses_for_modeling_mobility,
                             REML = TRUE, 
                             lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))


## Mixed model linear, random intercept and slope
r0_lin_lmer_mobility <- lmer(r0 ~ shutdown_day_count +
                               (1 + shutdown_day_count | state), 
                             data = responses_for_modeling_mobility,
                             REML = TRUE, 
                             lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))


# Mixed model polynomial
r0_poly_lmer_mobility <- lmer(r0 ~ shutdown_day_count + I(shutdown_day_count^2) + 
                                (1 + shutdown_day_count | state),
                              data = responses_for_modeling_mobility,
                              REML = TRUE, 
                              lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))


## Mixed model cubic
r0_cubic_lmer_mobility <- lmer(r0 ~ shutdown_day_count + I(shutdown_day_count^2)  + I(shutdown_day_count^3) +
                                 (1 + shutdown_day_count | state), 
                               data = responses_for_modeling_mobility,
                               REML = TRUE, 
                               lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=4e5)))


anova(r0_lin_base_mobility, r0_lin_lmer_mobility, r0_poly_lmer_mobility, r0_cubic_lmer_mobility)

## Cubic model shows best fit 
summary(r0_cubic_lmer_mobility)




# Modeling: Predicting future R0 values -----------------------------------

r0_future <- data.frame(state = sort(rep(unique(responses_all$state), 240)), 
                        shutdown_day_count = rep(c(1:240), length(unique(responses_all$state)))) 


## Predicting values for official shutdown
r0_future_official <- r0_future %>% 
  mutate(r0_predicted = predict(r0_cubic_lmer_official, newdata = r0_future, allow.new.levels = TRUE)) %>% 
  ### For graphing purposes
  left_join(responses_for_modeling_official %>% select(state, shutdown_date_official, opening_date_official), by = "state") %>% 
  distinct()


## Predicting values for mobility shutdown
r0_future_mobility <- r0_future %>%
  mutate(r0_predicted = predict(r0_cubic_lmer_mobility, newdata = r0_future, allow.new.levels = TRUE)) %>% 
  ### For graphing purposes
  left_join(responses_for_modeling_mobility %>% select(state, shutdown_date_mobility, opening_date_initial), by = "state") %>% 
  distinct()




# Modeling: Plotting cubic model with official shutdown -------------------

r0_future_official_national_average <- r0_future_official %>% 
  group_by(shutdown_day_count) %>% 
  summarise(r0_predicted = mean(r0_predicted))

day_when_official_r0_intercepts_1 <- r0_future_official_national_average %>% 
  mutate(is_intercept_day = r0_predicted <= 1.0 & lag(r0_predicted) > 1.0) %>% 
  filter(is_intercept_day) %>% 
  filter(shutdown_day_count == max(shutdown_day_count)) %>% 
  select(shutdown_day_count) %>% 
  pull()

ggplot() +
  geom_hline(yintercept = 1, color = "grey30", linetype = "dashed") +
  geom_vline(xintercept = day_when_official_r0_intercepts_1, color = "red") +
  annotate(geom = "text",
           label = paste0("day ", day_when_official_r0_intercepts_1),
           x = day_when_official_r0_intercepts_1 + 8,
           y = 4,
           color = "red") +
  ### All states
  geom_line(aes(x = shutdown_day_count, y = r0_predicted, group = state),
            r0_future_official,
            alpha = 0.15) +
  ### National average
  geom_line(aes(x = shutdown_day_count, y = r0_predicted, linetype = "national average"),
            r0_future_official_national_average,
            size = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        legend.position = c(0.1, 0.03),
        legend.justification = c(0, 0),
        legend.background = element_blank()) +
  scale_linetype_manual(values = c("solid")) +
  scale_x_continuous(limits = c(0, 205), breaks = seq(0, 200, 25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 4.5), breaks = seq(0, 4.5, 0.5), expand = c(0, 0)) +
  labs(title = expression(Cubic~model~predictions~`for`~R[0]~using~official~shutdown),
       x = "length of official shutdown period (days)",
       y = expression(predicted~R[0]),
       linetype = "") +
  guides(linetype = guide_legend(reverse = TRUE)) +
  ggsave("figures/lmer_official.png", height = 6, width = 10, units = "in")




# Modeling: Plotting cubic model with mobility-defined shutdown -----------

r0_future_mobility_national_average <- r0_future_mobility %>%
  group_by(shutdown_day_count) %>% 
  summarise(r0_predicted = mean(r0_predicted))

day_when_mobility_r0_intercepts_1 <- r0_future_mobility_national_average %>% 
  mutate(is_intercept_day = r0_predicted <= 1.0 & lag(r0_predicted) > 1.0) %>% 
  filter(is_intercept_day) %>% 
  filter(shutdown_day_count == max(shutdown_day_count)) %>% 
  select(shutdown_day_count) %>% 
  pull()

ggplot() +
  geom_hline(yintercept = 1, color = "grey30", linetype = "dashed") +
  geom_vline(xintercept = day_when_mobility_r0_intercepts_1, color = "blue") +
  annotate(geom = "text",
           label = paste0("day ", day_when_mobility_r0_intercepts_1),
           x = day_when_mobility_r0_intercepts_1 + 8,
           y = 2.15,
           color = "blue") +
  ### All states
  geom_line(aes(x = shutdown_day_count, y = r0_predicted, group = state),
            r0_future_mobility,
            alpha = 0.15) +
  ### National average
  geom_line(aes(x = shutdown_day_count, y = r0_predicted, linetype = "national average"),
            r0_future_mobility_national_average,
            size = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        legend.position = c(0.1, 0.03),
        legend.justification = c(0, 0),
        legend.background = element_blank()) +
  scale_linetype_manual(values = c("solid")) +
  scale_x_continuous(limits = c(0, 130), breaks = seq(0, 150, 25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5), expand = c(0, 0)) +
  labs(title = expression(Cubic~model~predictions~`for`~R[0]~using~`mobility-defined`~shutdown),
       x = "length of mobility-defined shutdown period (days)",
       y = expression(predicted~R[0]),
       linetype = "") +
  guides(linetype = guide_legend(reverse = TRUE)) +
  ggsave("figures/lmer_mobility.png", height = 6, width = 10, units = "in")



