library(tidyverse)
library(lubridate)
library(fs)
n_models_to_archive <- 4

# dput(unique(read_csv("data/oc_city_data.csv", col_types = cols(.default = col_skip(), city = col_character()))$city))

city_folder_dict <- tibble(city_name = c(
  "Aliso Viejo", "Anaheim", "Brea", "Buena Park", "Capistrano Beach",
  "Corona Del Mar", "Costa Mesa", "Cypress", "Dana Point", "Foothill Ranch",
  "Fountain Valley", "Fullerton", "Garden Grove", "Huntington Beach",
  "Irvine", "La Habra", "La Palma", "Ladera Ranch", "Laguna Beach",
  "Laguna Hills", "Laguna Niguel", "Laguna Woods", "Lake Forest",
  "Los Alamitos", "Midway City", "Mission Viejo", "Newport Beach",
  "Newport Coast", "Orange", "Placentia", "Rancho Santa Margarita",
  "San Clemente", "San Juan Capistrano", "Santa Ana", "Seal Beach",
  "Stanton", "Sunset Beach", "Trabuco Canyon", "Tustin", "Villa Park",
  "Westminster", "Yorba Linda"
)) %>%
  mutate(city_folder = str_replace_all(str_to_lower(city_name), "\\s", "-"))


all_models <-
  tibble(path = dir_ls("code/results", recurse = T)) %>%
  filter(path_ext(path) == "") %>%
  mutate(start_date = ymd(str_sub(path, start = 14, end = 23)),
         end_date = ymd(str_sub(path, start = 25, end = 34))) %>%
  mutate(city_folder = path_file(path)) %>%
  left_join(city_folder_dict) %>%
  select(-city_folder) %>%
  left_join(city_folder_dict) %>%
  mutate(location_name = city_name) %>%
  replace_na(list(location_name = "Orange County"))

template <- read_lines("code/analysis_template.Rmd")

# Rmd's -------------------------------------------------------------------
rmd_tbl <-
  all_models %>%
  mutate(rmd = pmap(
    list(path = path,
         start_date = start_date,
         end_date = end_date,
         city_name = city_name,
         location_name = location_name),
    function(path, start_date, end_date, city_name, location_name) {
      rmd <- template

      rmd[rmd == "title: \"LOCATION_HERE, CA COVID Situation Report START_DATE_HERE - END_DATE_HERE\""] <-
        rmd[rmd == "title: \"LOCATION_HERE, CA COVID Situation Report START_DATE_HERE - END_DATE_HERE\""] %>%
        str_replace("LOCATION_HERE", location_name) %>%
        str_replace("START_DATE_HERE", format(start_date, "%b %e, %Y")) %>%
        str_replace("END_DATE_HERE", format(end_date, "%b %e, %Y"))

      rmd[rmd == "results_folder <- RESULTS_FOLDER_HERE"] <-
        rmd[rmd == "results_folder <- RESULTS_FOLDER_HERE"] %>%
        str_replace("RESULTS_FOLDER_HERE", str_c('\"', path, '\"'))

      rmd[rmd == "## LOCATION_HERE, CA COVID-19 Situation Report, `r format(last_day + 5, \"%B %e, %Y\")`"] <-
        rmd[rmd == "## LOCATION_HERE, CA COVID-19 Situation Report, `r format(last_day + 5, \"%B %e, %Y\")`"] %>%
        str_replace("LOCATION_HERE", location_name)

      rmd[rmd == "location_name <- NA"] <- rmd[rmd == "location_name <- NA"] %>%
        str_replace("NA", str_c('\"', location_name, '\"'))

      if(!is.na(city_name)) {
        rmd[rmd == "city_name <- NA"] <-
          rmd[rmd == "city_name <- NA"] %>%
          str_replace("NA", str_c('\"', city_name, '\"'))
      }

      rmd
    })) %>%
  mutate(rmd_path = path("analysis",
                         str_c(start_date, "_", end_date,
                               str_replace_na(str_c("_", city_folder), "")),
         ext = "Rmd"))

# Write Rmd's
invisible(map2(rmd_tbl$rmd, rmd_tbl$rmd_path, ~write_lines(x = .x, file = .y)))

# YML ---------------------------------------------------------------------
site_yml_old <- read_lines("analysis/_site.yml")
archive_menu_begin <- which(site_yml_old == "    icon: fas fa-archive") + 1
archive_menu_end <- which(site_yml_old == "  - text: About")

all_models_for_yml <-
  all_models %>%
  filter(is.na(city_folder)) %>%
  mutate(folder = path_file(path)) %>%
  select(start_date, end_date, folder) %>%
  tail(n_models_to_archive)

archive_menu_new <-
  str_c(
    str_c('      - text: \"', format(rev(all_models_for_yml$start_date), "%b %e"), " - ", format(rev(all_models_for_yml$end_date), "%b %e"), '\"'),
    str_c("        href: ", rev(all_models_for_yml$folder), ".html"),
    sep = "\n") %>%
  str_split(pattern = "\\n") %>%
  unlist()

site_yml <- c(site_yml_old[1:archive_menu_begin], archive_menu_new, site_yml_old[archive_menu_end:length(site_yml_old)])

# Write new site yml
write_lines(site_yml, file = "analysis/_site.yml")

# Build Site --------------------------------------------------------------
rmd_to_build <- c(path("analysis", all_models_for_yml$folder, ext = "Rmd"),
  dir_ls("analysis")[dir_ls("analysis") %>% path_file() %>% str_sub(end = 21) == tail(all_models_for_yml$folder, 1)],
  dir_ls("analysis")[dir_ls("analysis") %>% path_file() %>% str_starts("\\d|_", negate = T)]) %>%
  unique() %>%
  unname()

workflowr::wflow_build(files = rmd_to_build)

# Replace Homepage --------------------------------------------------------
file_copy(path = path("docs", tail(all_models_for_yml, 1)$folder, ext = "html"),
          new_path = path("docs", "index", ext = "html"),
          overwrite = T)


# Replace Cities ----------------------------------------------------------
rmd_tbl %>%
  drop_na() %>%
  group_by(city_name) %>%
  filter(end_date == max(end_date)) %>%
  mutate(html_path = path("docs", path_ext_set(path_file(rmd_path), "html"))) %>%
  mutate(html_dest = path("docs", city_folder, ext = "html")) %>%
  map2(.x = .$html_path, .y = .$html_dest, .f = ~file_copy(path = .x, new_path = .y, overwrite = T))

# View Site ---------------------------------------------------------------
workflowr::wflow_view()
