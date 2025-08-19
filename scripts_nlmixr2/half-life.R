library(dplyr)
library(purrr)
library(broom)

df <- readr::read_rds("data/pkpd_3cmt_transit.rds")

# helper to pick terminal points after Tmax (simple heuristic)
pick_terminal <- function(d){
  tmax <- d |> slice_max(DV_CP, n = 1, with_ties = FALSE) |> pull(time)
  cand <- d |> filter(time > tmax, DV_CP > 0)
  # take last 4 points; refine if needed
  cand |> slice_tail(n = min(4, nrow(cand)))
}

est_half_life <- function(df){
  df |>
    group_by(id) |>
    arrange(time, .by_group = TRUE) |>
    group_modify(~{
      term <- pick_terminal(.x)
      if(nrow(term) < 3) return(tibble(t_half_h = NA_real_))
      fit <- lm(log(DV_CP) ~ time, data = term)
      lam_z <- -coef(fit)[["time"]]
      tibble(t_half_h = log(2) / lam_z)
    }) |>
    ungroup() |> 
    select(-id) |> 
    summarise(
      mean_half_life = median(t_half_h)
    )
}

est_half_life(df)
