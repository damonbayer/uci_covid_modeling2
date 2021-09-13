library(tidyverse)


# prep_and_save_map_data-function ---------------------------------------------
prep_and_save_map_data <- function(
  neg_line_list_file,
  line_list_file,
  zip_code_file,
  cases_per = 100000, # Number of cases per cases_per people in zip per time frame
  tests_per = 100000, # Number of tests per tests_per people in geog_level per time frame
  reporting_delay = 5, # Drop most recent 5 days due to reporting delay
  path_to_save_folder = "map-covid-data" # Path to folder to save map data in
){

  oc_zips <- read_csv(zip_code_file, col_types = cols(zip = col_character()))


  # Read in and summarize OCHCA data by zip
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
    filter(posted_date >= ymd("2020-03-01")) %>%
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
    mutate(month_date = zoo::as.yearmon(posted_date, "%m/%Y")) %>%
    filter(posted_date < max(posted_date) - days(reporting_delay))


  # Summarize case data
  oc_zips$cases_scaled_pop <- oc_zips$population / cases_per

  cases_plot_data <- covid_zip_data %>%
    group_by(zip, month_date) %>%
    summarize(new_cases_in_frame = sum(new_cases)) %>%
    inner_join(oc_zips, by = "zip") %>%
    mutate(new_cases_scaled = new_cases_in_frame / cases_scaled_pop) %>%
    mutate(plot_var_cont = new_cases_scaled) %>%
    mutate(plot_var_dis = factor(
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
    select(zip, plot_var_dis, month_date, plot_var_cont)

  cases_legend_label <- paste0(
    "Reported cases per\n",
    prettyNum(cases_per, big.mark = ",", scientific = FALSE),
    " people"
  )


  # Summarize tests data
  oc_zips$tests_scaled_pop <- oc_zips$population / tests_per

  tests_plot_data <- covid_zip_data %>%
    group_by(zip, month_date) %>%
    summarize(new_tests_in_frame = sum(new_tests)) %>%
    inner_join(oc_zips, by = "zip") %>%
    mutate(new_tests_scaled = new_tests_in_frame / tests_scaled_pop) %>%
    mutate(plot_var_cont = new_tests_scaled) %>%
    mutate(plot_var_dis = factor(
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
    select(zip, plot_var_dis, month_date, plot_var_cont)

  tests_legend_label <- paste0(
    "Tests per\n",
    prettyNum(tests_per, big.mark = ",", scientific = FALSE),
    " people"
  )


  # Summarize positivity data
  pos_plot_data <- covid_zip_data %>%
    group_by(zip, month_date) %>%
    summarize(per_pos = 100 * sum(new_cases) / sum(new_tests)) %>%
    mutate(plot_var_cont = per_pos) %>%
    inner_join(oc_zips, by = "zip") %>%
    mutate(plot_var_dis = factor(
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
    select(zip, plot_var_dis, month_date, plot_var_cont)

  pos_legend_label <- paste0("Percent of COVID-19\ntest positive")


  # Save aggregated OCHCA data to dashboard
  writeLines(paste0(max(covid_zip_data$posted_date)), paste0(path_to_save_folder, "/max_date_in_map_data.txt"))

  write.csv(cases_plot_data, file = paste0(path_to_save_folder, "/cases_map_data.csv"))
  writeLines(cases_legend_label, paste0(path_to_save_folder, "/cases_legend_label.txt"))

  write.csv(tests_plot_data, file = paste0(path_to_save_folder, "/tests_map_data.csv"))
  writeLines(tests_legend_label, paste0(path_to_save_folder, "/tests_legend_label.txt"))

  write.csv(pos_plot_data, file = paste0(path_to_save_folder, "/pos_map_data.csv"))
  writeLines(pos_legend_label, paste0(path_to_save_folder, "/pos_legend_label.txt"))
}
