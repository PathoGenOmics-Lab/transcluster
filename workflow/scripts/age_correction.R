library(tidyverse)

correct.age <- function(df, agecol) {
    df %>%
    # Re-parse values
    mutate(
        !!agecol := case_when(
        # e.g. 50
        !is.na(as.numeric(!!agecol)) ~
            as.numeric(!!agecol),
        # e.g. 2 months / 1 month / 6 mos
        grepl("^[0-9]+ ([mM]onths?|mos?)$", !!agecol) ~
            parse_number(!!agecol) / 12,
        # e.g. 1 year / 2 Years
        grepl("^[0-9]+ [yY]ears?$", !!agecol) ~
            parse_number(!!agecol),
        # e.g. 2 years 3 months
        grepl("^[0-9]+ [yY]ears? [0-9]+ ([mM]onths?|mos?)$", !!agecol) ~
            as.numeric(
                gsub("^([0-9]+) [yY]ears? [0-9]+ ([mM]onths?|mos?)$", "\\1", !!agecol)
            ) +
            as.numeric(
                gsub("^[0-9]+ [yY]ears? ([0-9]+) ([mM]onths?|mos?)$", "\\1", !!agecol)
            ) / 12,
        ),
        # all the rest
        TRUE ~ NA
    ) %>%
    # Remove unrealistic age values
    mutate(!!agecol := ifelse(0 <= !!agecol & !!agecol <= 120, !!agecol, NA))
}
