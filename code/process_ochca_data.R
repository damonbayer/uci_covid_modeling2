library(tidyverse)
library(lubridate)
library(here)
library(data.table)
library(dtplyr)

line_list_path <- "data/from_OCHCA/8.2.21 release to UCI team.csv"
negative_line_list_path <- "data/from_OCHCA/All PCR tests updated 8.2.21.csv"

metadata_zip <- tibble(
  zip = c(
    90620L, 90621L, 90623L, 90630L, 90631L,
    90680L, 90720L, 90740L, 90742L, 92602L, 92603L, 92604L, 92606L,
    92610L, 92612L, 92614L, 92617L, 92618L, 92620L, 92624L, 92625L,
    92626L, 92627L, 92629L, 92630L, 92637L, 92646L, 92647L, 92648L,
    92649L, 92651L, 92653L, 92655L, 92656L, 92657L, 92660L, 92661L,
    92662L, 92663L, 92672L, 92673L, 92675L, 92677L, 92679L, 92683L,
    92688L, 92691L, 92692L, 92694L, 92701L, 92703L, 92704L, 92705L,
    92706L, 92707L, 92708L, 92780L, 92782L, 92801L, 92802L, 92804L,
    92805L, 92806L, 92807L, 92808L, 92821L, 92823L, 92831L, 92832L,
    92833L, 92835L, 92840L, 92841L, 92843L, 92844L, 92845L, 92861L,
    92865L, 92866L, 92867L, 92868L, 92869L, 92870L, 92886L, 92887L
  ),
  SeptIncid = c(
    34.35816727, 35.55884277, 19.28764305,
    13.54364178, 41.40847987, 31.72482885, 36.77991816, 10.53563151,
    0, 4.372349263, 9.908838684, 11.1719361, 30.2395906, 13.33570413,
    18.16728435, 18.18328754, 7.120478496, 48.88182818, 37.67603804,
    20.69536424, 4.007052412, 23.30718875, 22.76052674, 34.94331418,
    32.10435605, 0, 9.959437926, 13.97501965, 25.37679017, 26.18365524,
    16.74971735, 15.36308081, 35.98416697, 13.25286466, 10.26588646,
    18.67977124, 13.35470085, 0, 13.857453, 24.6634169, 27.29537002,
    43.1890818, 16.58846391, 16.86547484, 21.17062409, 22.83522104,
    21.01635072, 20.11774173, 13.67116296, 66.78044075, 51.95202078,
    57.30626511, 23.4867803, 27.42957457, 98.33254891, 9.820727091,
    23.38026705, 30.39249739, 50.75078946, 51.51139104, 40.73841283,
    67.47063252, 24.21112097, 27.64645711, 27.44647937, 19.69999719,
    55.35566012, 24.85089463, 22.22042663, 28.01012228, 52.06164098,
    50.84777102, 30.4460344, 35.38726943, 6.171061834, 24.49029572,
    8.64902266, 50.75111652, 57.10446758, 22.46433786, 47.23665564,
    24.20395869, 40.35900294, 22.54960914, 54.98350495
  )
) %>%
  left_join(as_tibble(fread("data/oc_zips.csv")))

metadata_city <- tibble(
  city = c("Aliso Viejo", "Anaheim", "Brea", "Buena Park", "Corona Del Mar",
           "Costa Mesa", "Cypress", "Dana Point", "Foothill Ranch", "Fountain Valley",
           "Fullerton", "Garden Grove", "Huntington Beach", "Irvine", "La Habra",
           "La Palma", "Ladera Ranch", "Laguna Beach", "Laguna Hills", "Laguna Niguel",
           "Laguna Woods", "Lake Forest", "Los Alamitos", "Midway City",
           "Mission Viejo", "Newport Beach", "Newport Coast", "Orange",
           "Placentia", "Rancho Santa Margarita", "San Clemente", "San Juan Capistrano",
           "Santa Ana", "Seal Beach", "Stanton", "Trabuco Canyon", "Tustin",
           "Villa Park", "Westminster", "Yorba Linda"),
  cases = c(3, 53, 2, 16, 0, 17, 2, 0, 0, 7, 11, 18, 17, 31, 7, 0, 2, 1,
            4, 1, 3, 5, 2, 0, 5, 4, 0, 24, 3, 4, 4, 4, 59, 5, 3, 1, 15, 2,
            10, 5),
  tests = c(55, 333, 29, 75, 11, 126, 34, 18, 10, 71, 117, 162, 213, 296,
            38, 9, 26, 12, 26, 51, 12, 61, 22, 3, 97, 55, 9, 192, 51, 45,
            37, 45, 325, 22, 16, 25, 103, 11, 79, 52)) %>%
  left_join(fread("data/oc_zips.csv") %>%
              group_by(city) %>%
              summarize(population = sum(population)) %>%
              as_tibble())

deaths_tbl <-
  fread(line_list_path, select = c("DtDeath", "Date.death.posted", "DeathDueCOVID", "Zip")) %>%
  filter(DeathDueCOVID == "Y") %>%
  rename(zip = Zip) %>%
  select(-DeathDueCOVID)

deaths_tbl_county <-
  deaths_tbl %>%
  select(-zip) %>%
  as_tibble() %>%
  pivot_longer(everything()) %>%
  count(name, value) %>%
  pivot_wider(names_from = name, values_from = n) %>%
  arrange(value) %>%
  mutate(DtDeath = replace_na(DtDeath, 0),
         Date.death.posted = replace_na(Date.death.posted, 0)) %>%
  select(date = value,
         deaths = DtDeath,
         deaths_calcat = Date.death.posted)


neg_line_list <-
  fread(negative_line_list_path,
        select = c("unique_num", "Specimen.Collected.Date", "TestResult", "Zip"),
        na.strings = "") %>%
  select(id = unique_num, date = Specimen.Collected.Date, zip = Zip, test_result = TestResult) %>%
  filter(!(is.na(id) | is.na(test_result))) %>%
  mutate(test_result = fct_collapse(str_to_lower(test_result),
                                    negative = "negative",
                                    positive = "positive",
                                    other = c("invalid", "inconclusive"))) %>%
  filter(date >= lubridate::ymd("2020-01-01")) %>%
  group_by(id) %>%
  arrange(date) %>%
  ungroup()

if(length(levels(as.data.table(neg_line_list)$test_result)) != 3) stop(str_c(c("New test result category not accounted for.",levels(neg_line_list$test_result)), collapse = "\n"))


# Create OC Data ----------------------------------------------------------
first_pos <-
  neg_line_list %>%
  as.data.table() %>%
  filter(test_result == "positive") %>%
  group_by(id) %>%
  summarise(first_pos_date = min(date))

neg_line_list_filtered <-
  left_join(neg_line_list, first_pos) %>%
  mutate(first_pos_date = replace_na(first_pos_date, lubridate::ymd("9999-12-31"))) %>%
  filter(date <= first_pos_date) %>%
  select(-first_pos_date) %>%
  distinct()

oc_data <-
  neg_line_list_filtered %>%
  count(date, test_result) %>%
  pivot_wider(names_from = test_result, values_from = n, values_fill = 0) %>%
  full_join(deaths_tbl_county) %>%
  as_tibble() %>%
  replace(is.na(.), 0) %>%
  mutate(cases = positive, tests = negative + positive + other) %>%
  select(date, cases, tests, deaths, deaths_calcat) %>%
  mutate(date = as_date(date))

# Create OC City Data -----------------------------------------------------
deaths_tbl_city <-
  deaths_tbl %>%
  select(zip, date = DtDeath) %>%
  count(date, zip, name = "deaths") %>%
  arrange(date) %>%
  as_tibble() %>%
  left_join(select(metadata_zip, zip, city)) %>%
  count(date, city, wt = deaths, name = "deaths") %>%
  drop_na()

neg_line_list_filtered_city <-
  neg_line_list_filtered %>%
  left_join(metadata_zip %>%
              mutate(zip = as.character(zip)) %>%
              select(zip, city)) %>%
  as_tibble() %>%
  drop_na() %>%
  count(date, test_result, city) %>%
  pivot_wider(names_from = test_result, values_from = n) %>%
  replace(is.na(.), 0)

oc_city_data <-
  full_join(neg_line_list_filtered_city, deaths_tbl_city) %>%
  replace(is.na(.), 0) %>%
  mutate(cases = positive, tests = negative + positive + other) %>%
  select(date, city, cases, tests, deaths) %>%
  arrange(city, date) %>%
  mutate(date = as_date(date))


# Create OC Zip Month Data ------------------------------------------------
oc_zip_month_data <-
  neg_line_list %>%
  mutate(zip = str_sub(zip, end = 5)) %>%
  filter(zip %in% metadata_zip$zip) %>%
  filter(date >= lubridate::ymd("2020-03-01")) %>%
  as_tibble() %>%
  drop_na() %>%
  mutate(year = lubridate::year(date),
         month = format.Date(date, "%m")) %>%
  count(year, month, test_result, zip) %>%
  pivot_wider(names_from = test_result, values_from = n) %>%
  replace_na(list(other = 0, positive = 0, negative = 0)) %>%
  mutate(cases = positive, tests = other + positive + negative) %>%
  select(zip, year, month, cases, tests) %>%
  pivot_longer(-c(zip, year, month)) %>%
  pivot_wider(names_from = c(name, year, month), values_from = value) %>%
  right_join(metadata_zip %>% mutate(zip = as.character(zip))) %>%
  select(zip, city, population, SeptIncid, starts_with("tests"), starts_with("cases")) %>%
  arrange(zip)


# Refresh map data --------------------------------------------------------
source(here::here("code", "process_map_data.R"))

prep_and_save_map_data(
  neg_line_list_file = negative_line_list_path,
  line_list_file = line_list_path,
  zip_code_file = here("data", "map_zipcodes.csv"),
  cases_per = 100000, # Number of cases per cases_per people in zip per time frame
  tests_per = 100000, # Number of tests per tests_per people in geog_level per time frame
  reporting_delay = 5, # Don't use the most recent 5 days in data
  path_to_save_folder = "~/research/uci_covid19_dashboard/data" # Path to folder to save map data in
)


# OC City Incid -----------------------------------------------------------
oc_city_incid <-
  metadata_city %>%
  mutate(pop_incid = population * cases / tests) %>%
  mutate(prop_pop = population / sum(population),
         prop_incid = pop_incid / sum(pop_incid)) %>%
  arrange(city) %>%
  select(-pop_incid)

# Create ECDF -------------------------------------------------------------
deaths_tbl <- as_tibble(deaths_tbl)
death_delay_ecdf <- ecdf(as.numeric(deaths_tbl$Date.death.posted - deaths_tbl$DtDeath))


# Plot Data ---------------------------------------------------------------
source('code/helper_functions.R')
dat <- oc_data %>%
  lump_oc_data(time_interval_in_days = 7) %>%
  mutate(prop_deaths_reported = death_delay_ecdf(as.numeric(max(oc_data$date) - end_date)))

# Plot data
dat %>%
  mutate(pos = cases / tests) %>%
  select(date = end_date, cases, tests, deaths, pos) %>%
  pivot_longer(-date) %>%
  ggplot(aes(date, value)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = scales::comma)


# Write Data --------------------------------------------------------------
write_csv(oc_data, "data/oc_data.csv")
write_csv(oc_city_data, "data/oc_city_data.csv")
write_csv(oc_city_data, "~/research/uci_covid19_dashboard/data/oc_city_data.csv")
write_csv(oc_zip_month_data, "~/research/uci_covid19_dashboard/data/oc_zip_month_data.csv")
write_csv(oc_city_incid, "data/oc_city_incidence.csv")
write_rds(death_delay_ecdf, file = "data/death_delay_ecdf.rds")
