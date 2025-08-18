# Solarized palette
solarized <- list(
  base03 = "#002b36", base02 = "#073642", base01 = "#586e75",
  base00 = "#657b83", base0  = "#839496", base1  = "#93a1a1",
  base2  = "#eee8d5", base3  = "#fdf6e3",
  yellow = "#b58900", orange = "#cb4b16", red   = "#dc322f",
  magenta= "#d33682", violet = "#6c71c4", blue  = "#268bd2",
  cyan   = "#2aa198", green  = "#859900"
)

# Discrete palette (ordered to read well and stay legible on both modes)
solarized_discrete <- c(
  solarized$blue, solarized$orange, solarized$green,
  solarized$violet, solarized$red, solarized$cyan, solarized$magenta, solarized$yellow
)

# Minimal Solarized theme that preserves your mk772 settings
theme_mk772_solarized <- function(mode = c("light","dark"), base_size = 14, base_family = "") {
  mode <- match.arg(mode)
  is_dark <- identical(mode, "dark")

  bg  <- if (is_dark) solarized$base03 else solarized$base3
  fg  <- if (is_dark) solarized$base0  else solarized$base00
  grid<- if (is_dark) solarized$base01 else "#d8d2bf"  # slightly softer than base1 for light
  strip_bg <- if (is_dark) solarized$base02 else solarized$base2

  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # keep mk772 feel
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.3, colour = grid),

      # solarized surfaces
      plot.background  = ggplot2::element_rect(fill = bg,  colour = NA),
      panel.background = ggplot2::element_rect(fill = bg,  colour = NA),

      # text & titles
      text             = ggplot2::element_text(colour = fg),
      axis.text        = ggplot2::element_text(colour = fg),
      axis.title       = ggplot2::element_text(colour = fg),
      plot.title.position = "plot",
      plot.title       = ggplot2::element_text(face = "bold", colour = if (is_dark) solarized$yellow else solarized$orange),

      # strips
      strip.background = ggplot2::element_rect(fill = strip_bg, colour = NA),
      strip.text       = ggplot2::element_text(face = "bold", colour = fg),

      # legend
      legend.background= ggplot2::element_rect(fill = bg,  colour = NA),
      legend.key       = ggplot2::element_rect(fill = bg,  colour = NA),
      legend.title     = ggplot2::element_text(colour = fg),
      legend.text      = ggplot2::element_text(colour = fg),
      legend.position  = "right"
    )
}

# Matching discrete scales
scale_color_solarized_d <- function(...){ ggplot2::scale_color_manual(values = solarized_discrete, ...) }
scale_fill_solarized_d  <- function(...){ ggplot2::scale_fill_manual(values = solarized_discrete, ...) }

# Matching continuous scales (smooth through cool→warm Solarized accents)
scale_color_solarized_c <- function(...){
  ggplot2::scale_color_gradientn(colours = c(solarized$cyan, solarized$blue, solarized$violet, solarized$orange), ...)
}
scale_fill_solarized_c <- function(...){
  ggplot2::scale_fill_gradientn(colours = c(solarized$cyan, solarized$blue, solarized$violet, solarized$orange), ...)
}
