library(tidyverse)
library(broom)
library(here)
library(lubridate)
library(viridis)
library(TTR)
library(ggrepel)
library(sf)
library(ggspatial)
library(zoo)
library(transformr)
library(animation)
library(gridExtra)

# prep-map-data-function ---------------------------------------------
## Reads in
## - line list data with all positive cases in OC
## - negative line list data with all tests done in OC
## - zip code file with zip codes, cities, and population data for OC
##   UCI zip code pop from: https://www.zipdatamaps.com/92697
## - shape file for area including and surrounding OC
## Returns list including file with new cases, new tests and new deaths by date for each city and zip

prep_map_data <- function(
    neg_line_list_file,
    line_list_file,
    zip_code_file,
    shp_file,
    ca_shp_file,
    return_covid_data = TRUE
  ){

  oc_zips <- read_csv(zip_code_file, col_types = cols(Zip = col_character())) %>%
    rename_all(str_to_lower)

  oc_cities <- oc_zips %>%
    group_by(city) %>%
    summarize(zip = list(zip), population = sum(population))

  shp <- st_read(shp_file, quiet = TRUE)
  ca_shp <- st_read(ca_shp_file, quiet = TRUE)


  if (return_covid_data) {
    neg_line_list <- read_csv(neg_line_list_file) %>%
      mutate(Specimen.Collected.Date = as.Date(
        Specimen.Collected.Date,
        format = "%m-%d-%Y"
      )) %>%
      select(
        id = IncidentID,
        posted_date = Specimen.Collected.Date,
        test_result = TestResult,
        zip = Zip
      ) %>%
      mutate(test_result = tolower(test_result)) %>%
      mutate(test_result = factor(case_when(
        test_result == "positive" ~ "positive",
        test_result == "negative" ~ "negative",
        test_result == "inconclusive" | test_result == "invalid" ~ "other"
      ))) %>%
      filter(!is.na(test_result)) %>%
      mutate(zip = str_sub(zip, end = 5)) %>%
      filter(posted_date >= ymd("2020-01-01")) %>%
      drop_na() %>%
      group_by(id) %>%
      arrange(posted_date) %>%
      ungroup()

    new_deaths_tbl <- read_csv(
      line_list_file,
      col_types = cols(
        .default = col_skip(),
        `DtDeath` = col_date("%Y-%m-%d"),
        `DeathDueCOVID` = col_character(),
        Zip = col_character()
      )) %>%
      drop_na() %>%
      select(posted_date = `DtDeath`, zip = Zip) %>%
      count(posted_date, zip, name = "new_deaths") %>%
      arrange(posted_date)

    first_pos <- neg_line_list %>%
      filter(test_result == "positive") %>%
      group_by(id) %>%
      summarise(first_pos = min(posted_date))

    neg_line_list_filtered <- left_join(neg_line_list, first_pos) %>%
      mutate(first_pos = replace_na(first_pos, lubridate::ymd("9999-12-31"))) %>%
      filter(posted_date <= first_pos) %>%
      select(-first_pos) %>%
      distinct()

    neg_line_list_filtered_city <- neg_line_list_filtered %>%
        right_join(oc_cities %>%  unnest(zip)) %>%
        count(posted_date, test_result, city) %>%
        pivot_wider(names_from = test_result, values_from = n)

    new_deaths_tbl_city <- new_deaths_tbl %>%
        right_join(oc_cities %>%  unnest(zip)) %>%
        drop_na() %>%
        count(posted_date, city, wt = new_deaths, name = "new_deaths")

    covid_city_data <- full_join(neg_line_list_filtered_city, new_deaths_tbl_city) %>%
        replace(is.na(.), 0) %>%
        mutate(new_cases = positive, new_tests = negative + positive + other) %>%
        select(posted_date, city, new_cases, new_tests, new_deaths) %>%
        arrange(city, posted_date) %>%
        mutate(month_date = as.yearmon(posted_date, "%m/%Y"))

    neg_line_list_filtered_zip <- neg_line_list_filtered %>%
        right_join(oc_zips) %>%
        count(posted_date, test_result, zip) %>%
        pivot_wider(names_from = test_result, values_from = n)

    new_deaths_tbl_zip <- new_deaths_tbl %>%
        right_join(oc_zips) %>%
        drop_na() %>%
        count(posted_date, zip, wt = new_deaths, name = "new_deaths")

    covid_zip_data <- full_join(neg_line_list_filtered_zip, new_deaths_tbl_zip) %>%
        replace(is.na(.), 0) %>%
        mutate(new_cases = positive, new_tests = negative + positive + other) %>%
        select(posted_date, zip, new_cases, new_tests, new_deaths) %>%
        arrange(zip, posted_date) %>%
        mutate(month_date = as.yearmon(posted_date, "%m/%Y"))
  } else {
    covid_city_data = NULL
    covid_zip_data = NULL
  }

  return(list(
    "covid_city_data" = covid_city_data,
    "covid_zip_data" = covid_zip_data,
    "oc_cities" = oc_cities,
    "oc_zips" = oc_zips,
    "shp" = shp,
    "ca_shp" = ca_shp,
    "min_date" = min(covid_zip_data$posted_date),
    "max_date" = max(covid_zip_data$posted_date)
  ))
}







# gen-city-map-labeled -----------------------------------------------------
gen_city_map_labeled <- function(
    shp_file = here("data/shape-files/socal-zip", "Zipcode_boundary_scag_2009.shp"),
    ca_shp_file = here("data/shape-files/ca-counties", "cnty19_1.shp"),
    road_shp_file = here("data/shape-files/ca-interstates", "tl_2015_06_prisecroads.shp"),
    zip_code_file = here("data", "map_zipcodes.csv")
  ) {

  # Load OC zip code file
  oc_zips <- read_csv(zip_code_file, col_types = cols(Zip = col_character())) %>%
    rename_all(str_to_lower)

  oc_cities <- oc_zips %>%
    group_by(city) %>%
    summarize(population = sum(population))

  # Load CA roads shp
  roads_shp <- st_read(road_shp_file, quiet = TRUE) %>%
    subset(RTTYP == "I")

  # Load usa shp file to get cover below southern CA
  ca_shp <- st_read(ca_shp_file, quiet = TRUE)

  ca_shp_oc_only <- ca_shp %>% filter(COUNTY_NAM == "Orange")

  # Load southern CA shp
  all_shp <- st_read(shp_file, quiet = TRUE)

  all_shp1 <- all_shp %>%
    group_by(NAME) %>%
    summarize()

  all_city_shp1 <- all_shp1 %>% # Plot for entire area of shape file
    left_join(oc_cities, by = c("NAME" = "city")) %>%
    mutate(ID = NAME)

  oc_city_shp1 <- all_shp1 %>% # Only want names for OC cities
    right_join(oc_cities, by = c("NAME" = "city")) %>%
    mutate(ID = NAME)

  names_to_drop <- c(
    "Buena Park",
    "Capistrano Beach",
    "Corona del Mar",
    "Cypress",
    "Fountain Valley",
    "La Palma",
    "Ladera Ranch",
    "Laguna Hills",
    "Laguna Niguel",
    "Laguna Woods",
    "Midway City",
    "Newport Beach",
    "Newport Coast",
    "Placentia",
    "Rancho Santa Margarita",
    "San Juan Capistrano",
    "Stanton",
    "Sunset Beach",
    "Trabuco Canyon",
    "Villa Park",
    "Westminster"
  )

  oc_city_shp_renamed <- oc_city_shp1 %>%  # Drop names that don't place well
    mutate(NAME = ifelse(NAME %in% names_to_drop, "", NAME))


  # Plot of cities with roads and UCI
  ggplot(all_city_shp1) +
    geom_sf(fill = "khaki1", color = "gray60") +
    geom_sf(data = ca_shp, fill = "khaki1", color = "gray60") +
    geom_sf(data = ca_shp_oc_only, fill = "khaki1", color = "black", size = 1.25) +
    geom_sf(data = oc_city_shp1, color = "black", fill = "khaki1") +
    geom_sf(data = roads_shp, color = "gray", fill = "white", size = 1) +
    coord_sf( # Outlines Orange County
      xlim = c(396639.2, 461568.2),
      ylim = c(3694363, 3759819),
      expand = FALSE
    ) +
    annotation_scale(
      location = "br",
      width_hint = 0.2,
      pad_x = unit(0.25, "cm"),
      pad_y = unit(0.25, "cm")
    ) +
    annotation_north_arrow(
      location = "br",
      which_north = "true",
      pad_x = unit(0.25, "cm"), pad_y = unit(0.5, "cm"),
      style = north_arrow_fancy_orienteering
    )  +
    theme_void() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA),
      panel.background = element_rect(fill = "darkslategray2")
    ) +
    geom_sf_text(
      data = oc_city_shp_renamed,
      mapping = aes(label = NAME),
      size = 4,
      color = "darkgreen"
    ) +
    annotate( # UCI label
      geom = "text",
      x = 422000,
      y = 3722000,
      label = "UC Irvine",
      color = "blue",
      size = 4,
      fontface = "bold"
    ) +
    annotate( # UCI point
      geom = "point",
      x = 421873.3807,
      y = 3723199.945,
      color = "blue",
      size = 4,
      pch = 17
    ) +
    annotate( # I-15 label
      geom = "text",
      x = 456000,
      y = 3740000,
      label = "I-15",
      color = "black",
      size = 4,
      fontface="bold"
    ) +
    annotate( # I-405 label
      geom = "text",
      x = 422000,
      y = 3727500,
      label = "I-405",
      color = "black",
      size = 4,
      fontface="bold"
    ) +
    annotate( # I-5 label
      geom = "text",
      x = 430000,
      y = 3730000,
      label = "I-5",
      color = "black",
      size = 4,
      fontface="bold"
    ) +
    annotate( # I-5 label
      geom = "text",
      x = 403250,
      y = 3753500,
      label = "I-5",
      color = "black",
      size = 4,
      fontface = "bold"
    ) +
    annotate( # I-605 label
      geom = "text",
      x = 400000,
      y = 3750000,
      label = "I-605",
      color = "black",
      size = 4,
      fontface="bold"
    )
}








# gen-map-function-gif --------------------------------------------------------
gen_map_gif <- function(
  plot_data,
  shp,
  ca_shp,
  legend_label,
  geog_level,
  month_seq,
  var_type
){

  ca_shp_oc_only <- ca_shp %>% filter(COUNTY_NAM == "Orange")

  map_labeled <- gen_city_map_labeled()

  prev_directory <- getwd()
  setwd("docs/")

  saveHTML({

    ani.options(verbose = FALSE)

    for (i in 1:length(month_seq)) {
      curr_time <- month_seq[i]
      plot_data_sub <- plot_data %>% filter(month_date == curr_time)

      if (geog_level == "city") {
        shp1 <- shp %>%
          group_by(NAME) %>%
          summarize()

        all_shp1 <- shp1 %>% # Plot for entire area of shape file
          full_join(plot_data_sub, by = c("NAME" = "city")) %>%
          mutate(ID = NAME)

        oc_shp1 <- shp1 %>% # Only want names for OC cities
          right_join(plot_data_sub, by = c("NAME" = "city")) %>%
          mutate(ID = NAME)
      } else if (geog_level == "zip") {
        shp1 <- shp
        all_shp1 <- shp %>% # Plot for entire area of shape file
          full_join(plot_data_sub, by = c("KEY_" = "zip")) %>%
          mutate(ID = KEY_)

        oc_shp1 <- shp %>% # Only want names for OC cities
          right_join(plot_data_sub, by = c("KEY_" = "zip")) %>%
          mutate(ID = KEY_)
      }

      curr_map <- ggplot(shp1) +
        geom_sf(fill = "khaki1", color = "gray60") +
        geom_sf(data = ca_shp, fill = "khaki1", color = "gray60") +
        geom_sf(data = ca_shp_oc_only, fill = "khaki1", color = "black", size = 1.25) +
        geom_sf(
          data = oc_shp1,
          mapping = aes(fill = plot_var, group = ID),
          color = "black"
        ) +
        coord_sf( # Outlines Orange County
          xlim = c(396639.2, 461568.2),
          ylim = c(3694363, 3759819),
          expand = FALSE
        ) +
        labs(fill = legend_label) +
        theme_void() +
        theme(
          panel.background = element_rect(fill = "darkslategray2"),
          #legend.position = "bl",
          legend.justification = c(0, 0),
          legend.position = c(0.01, 0.01),
          legend.background = element_rect(fill = "white"),
          # legend.title = element_text(size = 5, color = "black"),
          # legend.text = element_text(size = 5, color = "black"),
          legend.margin = margin(1, 1, 1, 1),
          panel.border = element_rect(colour = "black", fill = NA)
        ) +
        annotation_scale(
          location = "br",
          width_hint = 0.2,
          pad_x = unit(0.25, "cm"),
          pad_y = unit(0.25, "cm")
        ) +
        annotation_north_arrow(
          location = "br",
          which_north = "true",
          pad_x = unit(0.25, "cm"),
          pad_y = unit(0.5, "cm"),
          style = north_arrow_fancy_orienteering
        ) +
        annotate(
          geom = "text",
          x = 450000,
          y = 3755000,
          label = curr_time,
          color = "black",
          size = 7,
          fontface="bold"
        ) +
        scale_fill_viridis(drop = FALSE, discrete = TRUE, direction = -1)

      grid.arrange(
        curr_map,
        map_labeled,
        ncol = 2
      )

      ani.pause()
    }
  },
  ani.width = 800,
  ani.height = 400,
  img.name = paste0(var_type, "-map-gif.png"),
  htmlfile = paste0(var_type, "-map-gif.html"),
  verbose = FALSE
  )

  setwd(prev_directory)
}


# map-cases-function-gif ------------------------------------------------------
## Map number of COVID-19 cases to OC geog_level ("city" or "zip")
map_cases_gif <- function(
  map_data_list,
  geog_level,
  month_seq,
  cases_per = 100000 # Number of cases per cases_per people in zip per time frame
){

  legend_label <- paste0(
    "Reported cases per\n",
    prettyNum(cases_per, big.mark = ",", scientific = FALSE),
    " people"
  )

  if (is.null(map_data_list$oc_cities) | is.null(map_data_list$oc_zips)) {
    stop("Missing zip/city file in map_data_list")
  } else if (geog_level == "city") {
    oc_data <- map_data_list$oc_cities
    covid_data <- map_data_list$covid_city_data
  } else if (geog_level == "zip") {
    oc_data <- map_data_list$oc_zips
    covid_data <- map_data_list$covid_zip_data
  }

  shp <- map_data_list$shp
  ca_shp <- map_data_list$ca_shp

  oc_data$scaled_pop <- oc_data$population / cases_per


  if (geog_level == "city") {
    plot_data <- covid_data %>%
      group_by(city, month_date) %>%
      summarize(new_cases_in_frame = sum(new_cases)) %>%
      inner_join(oc_data, by = "city") %>%
      mutate(new_cases_scaled = new_cases_in_frame / scaled_pop) %>%
      mutate(plot_var = factor(
        case_when(
          new_cases_scaled <= 5 ~ "0-5",
          new_cases_scaled <= 30 ~ "5-30",
          new_cases_scaled <= 60 ~ "30-60",
          new_cases_scaled <= 120 ~ "60-120",
          new_cases_scaled <= 240 ~ "120-240",
          new_cases_scaled <= 480 ~ "240-480",
          new_cases_scaled > 480 ~ ">480"
        ),
        levels = c(
          "0-5",
          "5-30",
          "30-60",
          "60-120",
          "120-240",
          "240-480",
          ">480"
        ))) %>%
      filter(month_date %in% month_seq) %>%
      select(city, plot_var, month_date)

  } else if (geog_level == "zip") {
    plot_data <- covid_data %>%
      group_by(zip, month_date) %>%
      summarize(new_cases_in_frame = sum(new_cases)) %>%
      inner_join(oc_data, by = "zip") %>%
      mutate(new_cases_scaled = new_cases_in_frame / scaled_pop) %>%
      mutate(plot_var = factor(
        case_when(
          new_cases_scaled <= 5 ~ "0-5",
          new_cases_scaled <= 30 ~ "5-30",
          new_cases_scaled <= 60 ~ "30-60",
          new_cases_scaled <= 120 ~ "60-120",
          new_cases_scaled <= 240 ~ "120-240",
          new_cases_scaled <= 480 ~ "240-480",
          new_cases_scaled > 480 ~ ">480"
        ),
        levels = c(
          "0-5",
          "5-30",
          "30-60",
          "60-120",
          "120-240",
          "240-480",
          ">480"
        ))) %>%
      filter(month_date %in% month_seq) %>%
      select(zip, plot_var, month_date)
  }

  gen_map_gif(
    plot_data = plot_data,
    shp = shp,
    ca_shp = ca_shp,
    legend_label = legend_label,
    geog_level = geog_level,
    month_seq = month_seq,
    var_type = "cases"
  )
}


# map-tests-function-gif ------------------------------------------------------
## Map number of COVID-19 tests to geog_level ("city", or "zip")
map_tests_gif <- function(
  map_data_list,
  geog_level,
  month_seq,
  tests_per = 100000 # Number of cases per tests_per people in geog_level per time frame
){

  legend_label <- paste0(
    "Tests per\n",
    prettyNum(tests_per, big.mark = ",", scientific = FALSE),
    " people"
  )


  if (is.null(map_data_list$oc_cities) | is.null(map_data_list$oc_zips)) {
    stop("Missing zip/city file in map_data_list")
  } else if (geog_level == "city") {
    oc_data <- map_data_list$oc_cities
    covid_data <- map_data_list$covid_city_data
  } else if (geog_level == "zip") {
    oc_data <- map_data_list$oc_zips
    covid_data <- map_data_list$covid_zip_data
  }

  shp <- map_data_list$shp
  ca_shp <- map_data_list$ca_shp

  oc_data$scaled_pop <- oc_data$population / tests_per


  if(geog_level == "city"){
    plot_data <- covid_data %>%
      group_by(city, month_date) %>%
      summarize(new_tests_in_frame = sum(new_tests)) %>%
      inner_join(oc_data, by = "city") %>%
      mutate(new_tests_scaled = new_tests_in_frame / scaled_pop) %>%
      mutate(plot_var = factor(
        case_when(
          new_tests_scaled <= 100 ~ "0-100",
          new_tests_scaled <= 200 ~ "100-200",
          new_tests_scaled <= 600 ~ "200-600",
          new_tests_scaled <= 1200 ~ "600-1200",
          new_tests_scaled <= 2400 ~ "1200-2400",
          new_tests_scaled > 2400 ~ ">2400"
        ),
        levels = c(
          "0-100",
          "100-200",
          "200-600",
          "600-1200",
          "1200-2400",
          ">2400"
        )
      )) %>%
      filter(month_date %in% month_seq) %>%
      select(city, plot_var, month_date)

  } else if(geog_level == "zip"){
    plot_data <- covid_data %>%
      group_by(zip, month_date) %>%
      summarize(new_tests_in_frame = sum(new_tests)) %>%
      inner_join(oc_data, by = "zip") %>%
      mutate(new_tests_scaled = new_tests_in_frame / scaled_pop) %>%
      mutate(plot_var = factor(
        case_when(
          new_tests_scaled <= 100 ~ "0-100",
          new_tests_scaled <= 200 ~ "100-200",
          new_tests_scaled <= 600 ~ "200-600",
          new_tests_scaled <= 1200 ~ "600-1200",
          new_tests_scaled <= 2400 ~ "1200-2400",
          new_tests_scaled > 2400 ~ ">2400"
        ),
        levels = c(
          "0-100",
          "100-200",
          "200-600",
          "600-1200",
          "1200-2400",
          ">2400"
        )
      )) %>%
      filter(month_date %in% month_seq) %>%
      select(zip, plot_var, month_date)
  }

  gen_map_gif(
    plot_data = plot_data,
    shp = shp,
    ca_shp = ca_shp,
    legend_label = legend_label,
    geog_level = geog_level,
    month_seq = month_seq,
    var_type = "tests"
  )
}








# map-per-pos-function-gif ----------------------------------------------------
## Map percent of positive COVID-19 cases in geog_level to geog_level ("city", or "zip")
map_per_pos_gif <- function(
  map_data_list,
  geog_level,
  month_seq
){

  legend_label <- paste0("Percent of COVID-19\ntest positive")

  if (is.null(map_data_list$oc_cities) | is.null(map_data_list$oc_zips)) {
    stop("Missing zip/city file in map_data_list")
  } else if (geog_level == "city") {
    oc_data <- map_data_list$oc_cities
    covid_data <- map_data_list$covid_city_data
  } else if (geog_level == "zip") {
    oc_data <- map_data_list$oc_zips
    covid_data <- map_data_list$covid_zip_data
  }

  shp <- map_data_list$shp
  ca_shp <- map_data_list$ca_shp


  if (geog_level == "city") {
    plot_data <- covid_data %>%
      group_by(city, month_date) %>%
      summarize(per_pos = 100 * sum(new_cases) / sum(new_tests)) %>%
      inner_join(oc_data, by = "city") %>%
      mutate(plot_var = factor(
        case_when(
          per_pos <= 6 ~ "0-6",
          per_pos <= 12 ~ "6-12",
          per_pos <= 18 ~ "12-18",
          per_pos <= 24 ~ "18-24",
          per_pos <= 30 ~ "24-30",
          per_pos > 30 ~ ">30"
        ),
        levels = c(
          "0-6",
          "6-12",
          "12-18",
          "18-24",
          "24-30",
          ">30"
        )
      )) %>%
      filter(month_date %in% month_seq) %>%
      select(city, plot_var, month_date)

  } else if (geog_level == "zip") {
    plot_data <- covid_data %>%
      group_by(zip, month_date) %>%
      summarize(per_pos = 100 * sum(new_cases) / sum(new_tests)) %>%
      inner_join(oc_data, by = "zip") %>%
      mutate(plot_var = factor(
        case_when(
          per_pos <= 6 ~ "0-6",
          per_pos <= 12 ~ "6-12",
          per_pos <= 18 ~ "12-18",
          per_pos <= 24 ~ "18-24",
          per_pos <= 30 ~ "24-30",
          per_pos > 30 ~ ">30"
        ),
        levels = c(
          "0-6",
          "6-12",
          "12-18",
          "18-24",
          "24-30",
          ">30"
        )
      )) %>%
      filter(month_date %in% month_seq) %>%
      select(zip, plot_var, month_date)
  }

  gen_map_gif(
    plot_data = plot_data,
    shp = shp,
    ca_shp = ca_shp,
    legend_label = legend_label,
    geog_level = geog_level,
    month_seq = month_seq,
    var_type = "positivity"
  )
}




# refresh-all-3-map-gifs ------------------------------------------------------
refresh_all_3_map_gifs <- function(
    neg_line_list_file,
    line_list_file,
    zip_code_file = here("data", "map_zipcodes.csv"),
    shp_file = here("data/shape-files/socal-zip", "Zipcode_boundary_scag_2009.shp"),
    ca_shp_file = here("data/shape-files/ca-counties", "cnty19_1.shp"),
    first_date = as.Date("2020/03/01")
  ) {

  map_covid_data <- prep_map_data(
      neg_line_list_file = neg_line_list_file,
      line_list_file = line_list_file,
      zip_code_file = zip_code_file,
      shp_file = shp_file,
      ca_shp_file = ca_shp_file
  )


  # Do not use the most recent 5 days in the data
  max_date_to_use <- map_covid_data$max_date - days(5)

  map_covid_data$covid_city_data <- map_covid_data$covid_city_data %>%
    filter(posted_date <= max_date_to_use)

  map_covid_data$covid_zip_data <- map_covid_data$covid_zip_data %>%
    filter(posted_date <= max_date_to_use)

  date_seq <- seq(first_date, max_date_to_use, by = "months")


  map_cases_gif(
    map_data_list = map_covid_data,
    geog_level = "zip",
    month_seq = as.yearmon(date_seq, "%m/%Y")
  )

  map_tests_gif(
    map_data_list = map_covid_data,
    geog_level = "zip",
    month_seq = as.yearmon(date_seq, "%m/%Y")
  )

  map_per_pos_gif(
    map_data_list = map_covid_data,
    geog_level = "zip",
    month_seq = as.yearmon(date_seq, "%m/%Y")
  )

}

