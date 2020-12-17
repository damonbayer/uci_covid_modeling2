library(tidyverse)
library(lubridate)
library(fs)

all_models <-
  tibble(path = dir_ls("code/results")) %>%
  mutate(folder = path_file(path)) %>%
  separate(folder, c("start_date", "end_date"), sep = "_", remove = F) %>%
  mutate(across(ends_with("date"), ymd))

template = read_lines("code/analysis_template.Rmd")


# Rmd's -------------------------------------------------------------------
rmd_tbl <-
  all_models %>%
  mutate(rmd = pmap(
    list(path = path,
         start_date = start_date,
         end_date = end_date),
    function(path, start_date, end_date) {
      rmd <- template
      rmd[2] <- str_replace(rmd[2], "START_DATE - END_DATE", str_c(format(start_date, "%b %e, %Y"), "-", format(end_date, "%b %e, %Y"), sep = " "))
      rmd[11] <- str_replace(rmd[11], "path_here", path)
      rmd
  })) %>%
  mutate(rmd_path = path("analysis", str_c(start_date, "_", end_date), ext = "Rmd"))

# Write Rmd's
invisible(map2(rmd_tbl$rmd, rmd_tbl$rmd_path, ~write_lines(x = .x, path = .y)))

# YML ---------------------------------------------------------------------
site_yml_old <- read_lines("analysis/_site.yml")
archive_menu_begin <- which(site_yml_old == "  - text: Archived Reports") + 1
archive_menu_end <- which(site_yml_old == "  - text: About")

archive_menu_new <-
  str_c(
    str_c('      - text: \"', format(rev(all_models$start_date), "%b %e"), " - ", format(rev(all_models$end_date), "%b %e"), '\"'),
    str_c("        href: ", rev(all_models$folder), ".html"),
    sep = "\n") %>%
  str_split(pattern = "\\n") %>%
  unlist()

site_yml <-c(site_yml_old[1:archive_menu_begin], archive_menu_new, site_yml_old[archive_menu_end:length(site_yml_old)])

# Write new site yml
write_lines(site_yml, path = "analysis/_site.yml")

# Build Site --------------------------------------------------------------
# workflowr::wflow_build(files = "analysis/2020-03-29_2020-05-03.Rmd")
workflowr::wflow_build()

# Replace Homepage --------------------------------------------------------
file_copy(path = path("docs", tail(all_models, 1)$folder, ext = "html"),
          new_path = path("docs", "index", ext = "html"),
          overwrite = T)

# View Site ---------------------------------------------------------------
wflow_view()
