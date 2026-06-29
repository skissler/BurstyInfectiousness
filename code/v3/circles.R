 d_min <- .05
d_max <- .20

  list.files("output", pattern = "^pest_periodic_.*\\.csv$", full.names = TRUE) %>%
    map_dfr(read_csv, show_col_types = FALSE) %>%
    group_by(pathogen) %>%
    mutate(
      p_min  = min(p_extinct),
      p_max  = max(p_extinct),
      diameter = sqrt(d_min^2 + (d_max^2 - d_min^2) * (p_extinct - p_min) / (p_max - p_min))
    ) %>%
    ungroup() %>%
    select(pathogen, c_amp, psi, p_extinct, diameter) %>%
    arrange(pathogen, c_amp, psi) %>%
    print(n = Inf)


getdiameter <- function(p_extinct, d_min=.05, d_max=.20, p_min=0, p_max=1){
  diameter = sqrt(d_min^2 + (d_max^2 - d_min^2) * (p_extinct - p_min) / (p_max - p_min))
  return(diameter)
}