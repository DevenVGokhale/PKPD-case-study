# app.R
suppressPackageStartupMessages({
  library(shiny)
  library(rxode2)
  library(ggplot2)
})

# --------- Defaults (edit with your "true" values) ----------
true_pars <- list(
  CL  = 5.0,   # L/h
  Vc  = 50.0,  # L
  Q2  = 10.0,  # L/h (A1 <-> A2)
  Vp1 = 100.0, # L
  Q3  = 8.0,   # L/h (A1 <-> A3)
  Vp2 = 80.0,  # L
  MTT = 2.0,   # h (N=3 transit comps => ktr = 4/MTT)
  F   = 1.0,   # fraction
  sigma = 0.0  # mg/L, optional noise overlay for plotting
)

# --------- rxode2 model: your formulation (no theta_* indirection) ----------
model <- rxode2({
  # Transit rate (N=3 → N+1 = 4)
  ktr <- 4 / MTT

  # Transit absorption
  d/dt(A0)  = -ktr * A0;
  d/dt(TR1) =  ktr * A0 - ktr * TR1;
  d/dt(TR2) =  ktr * TR1 - ktr * TR2;
  d/dt(TR3) =  ktr * TR2 - ktr * TR3;

  # 3CMT PK
  d/dt(A1) = ktr * TR3 - (CL / Vc) * A1
             - Q2 / Vc * A1 + Q2 / Vp1 * A2
             - Q3 / Vc * A1 + Q3 / Vp2 * A3;
  d/dt(A2) = Q2 / Vc * A1 - Q2 / Vp1 * A2;
  d/dt(A3) = Q3 / Vc * A1 - Q3 / Vp2 * A3;

  Cp = A1 / Vc;
})

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("3CMT with Transit Absorption (rxode2)"),
  sidebarLayout(
    sidebarPanel(
      h4("PK parameters"),
      numericInput("CL",  "CL (L/h)",  value = true_pars$CL,  min = 0, step = 0.1),
      numericInput("Vc",  "Vc (L)",    value = true_pars$Vc,  min = 1, step = 1),
      numericInput("Q2",  "Q2 (L/h)",  value = true_pars$Q2,  min = 0, step = 0.1),
      numericInput("Vp1", "Vp1 (L)",   value = true_pars$Vp1, min = 1, step = 1),
      numericInput("Q3",  "Q3 (L/h)",  value = true_pars$Q3,  min = 0, step = 0.1),
      numericInput("Vp2", "Vp2 (L)",   value = true_pars$Vp2, min = 1, step = 1),
      tags$hr(),
      h4("Transit / F"),
      numericInput("MTT", "MTT (h)",          value = true_pars$MTT, min = 0.1, step = 0.1),
      sliderInput("F",    "Bioavailability",  min = 0, max = 1, value = true_pars$F, step = 0.01),
      tags$hr(),
      h4("Dosing"),
      numericInput("dose", "Oral dose (mg)",        value = 100, min = 0, step = 10),
      numericInput("ii",   "Dosing interval τ (h)", value = 12,  min = 0.1, step = 1),
      numericInput("ndose","Number of doses",       value = 10,  min = 1, step = 1),
      tags$hr(),
      h4("Sampling grid"),
      numericInput("t_end", "Simulation end (h)", value = 240, min = 1, step = 1),
      numericInput("dt",    "Δt (h)",            value = 0.25, min = 0.01, step = 0.01),
      actionButton("apply_grid", "Apply sampling grid", class = "btn-info"),
      tags$hr(),
      h4("Plot options"),
      checkboxInput("show_points", "Overlay noisy observations", value = FALSE),
      numericInput("sigma", "Residual SD (mg/L)", value = true_pars$sigma, min = 0, step = 0.05),
      actionButton("run", "Simulate", class = "btn-primary")
    ),
    mainPanel(
      uiOutput("grid_summary"),
      plotOutput("ct_plot", height = 420),
      tags$hr(),
      radioButtons(
        "metric_window",
        label = "Exposure window",
        choices = c("Overall (0–t_end)" = "overall",
                    "Last dosing interval (≈ steady-state)" = "tau"),
        selected = "overall",
        inline = TRUE
      ),
      br(),
      tableOutput("metrics_tbl"),
      downloadButton("download_metrics", "Download metrics (CSV)")
    )
  )
)

# ---------- SERVER ----------
server <- function(input, output, session) {

  # keep sigma at 0 if overlay is off
  observe({
    if (!input$show_points && input$sigma != 0) {
      updateNumericInput(session, "sigma", value = 0)
    }
  })

  # hold the currently applied grid
  grid_times <- reactiveVal(seq(0, 240, by = 0.25))

  observeEvent(input$apply_grid, {
    validate(need(input$t_end > 0, "t_end must be > 0"),
             need(input$dt > 0,   "dt must be > 0"))
    times <- seq(0, input$t_end, by = input$dt)
    if (tail(times, 1) < input$t_end) times <- c(times, input$t_end)
    grid_times(times)
  })

  output$grid_summary <- renderUI({
    times <- grid_times()
    n <- length(times)
    div(
      tags$strong("Applied sampling grid: "),
      sprintf("n=%d, start=%.3f h, end=%.3f h, dt≈%.3f h",
              n, head(times, 1), tail(times, 1),
              ifelse(n > 1, diff(times)[1], NA_real_))
    )
  })

  # parameter vector passed directly as 'params='
  param_vec <- reactive({
    c(
      CL  = input$CL,
      Vc  = input$Vc,
      Q2  = input$Q2,
      Vp1 = input$Vp1,
      Q3  = input$Q3,
      Vp2 = input$Vp2,
      MTT = input$MTT
    )
  })

  # initial conditions (all zero)
  init_vec <- reactive({
    c(A0 = 0, TR1 = 0, TR2 = 0, TR3 = 0, A1 = 0, A2 = 0, A3 = 0)
  })

  # build event table (NOTE: use named 'time=' when adding sampling times)
  make_et <- reactive({
    validate(
      need(length(grid_times()) > 0, "Apply a sampling grid first (or keep the default)."),
      need(input$ndose >= 1, "Number of doses must be ≥ 1.")
    )
    et() |>
      et(amt = input$F * input$dose, cmt = "A0", ii = input$ii,
         addl = input$ndose - 1, time = 0) |>
      et(time = grid_times())
  })

  # simulate (noise-free curve for metrics)
  sim_df <- eventReactive(input$run, {
    rxSolve(
      model,
      params = param_vec(),   # <-- direct params (no theta)
      events = make_et(),
      inits  = init_vec(),
      atol = 1e-8, rtol = 1e-8
    ) |>
      as.data.frame() |>
      (\(d) d[, c("time", "Cp")])()
  }, ignoreInit = TRUE)

  # optional noisy overlay for plotting only
  sim_with_noise <- reactive({
    df <- sim_df()
    if (is.null(df)) return(NULL)
    if (!input$show_points || input$sigma <= 0) return(df)
    set.seed(123)
    df$Cp_obs <- pmax(0, df$Cp + rnorm(nrow(df), 0, input$sigma))
    df
  })

  output$ct_plot <- renderPlot({
    req(sim_with_noise())
    df <- sim_with_noise()
    p <- ggplot(df, aes(time, Cp)) +
      geom_line() +
      labs(x = "Time (h)", y = "Concentration (mg/L)",
           title = "3CMT + Transit (N=3): Concentration–Time") +
      theme_bw()
    if (input$show_points && "Cp_obs" %in% names(df)) {
      p <- p + geom_point(aes(y = Cp_obs))
    }
    p
  })

  # ---- exposure metrics helpers ----
  trapz <- function(x, y) {
    if (length(x) < 2) return(NA_real_)
    sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  }

  compute_metrics_overall <- function(df) {
    auc <- trapz(df$time, df$Cp)
    cmax_i <- which.max(df$Cp)
    cmax <- df$Cp[cmax_i]
    tmax <- df$time[cmax_i]
    cmin <- min(df$Cp, na.rm = TRUE)
    data.frame(
      metric = c("AUC_0_tend", "Cmax", "Tmax", "Cmin"),
      value  = c(auc, cmax, tmax, cmin)
    )
  }

  compute_metrics_tau <- function(df, ii, ndose) {
    if (ndose < 1 || ii <= 0) {
      return(data.frame(metric = "Note", value = "Need ndose ≥ 1 and τ > 0"))
    }
    t0 <- (ndose - 1) * ii
    t1 <- ndose * ii
    if (max(df$time) < t1) {
      return(data.frame(
        metric = "Note",
        value  = sprintf("Extend t_end to at least %.3f h to cover the last interval.", t1)
      ))
    }
    sel <- df[df$time >= t0 & df$time <= t1, , drop = FALSE]
    if (nrow(sel) < 2) {
      return(data.frame(metric = "Note", value = "Sampling grid too coarse in the last interval"))
    }
    auc_tau <- trapz(sel$time, sel$Cp)
    cmax_i  <- which.max(sel$Cp)
    cmax_ss <- sel$Cp[cmax_i]
    tmax_ss <- sel$time[cmax_i]
    cmin_ss <- min(sel$Cp, na.rm = TRUE)
    cav_ss  <- auc_tau / ii
    fluct   <- (cmax_ss - cmin_ss) / cav_ss
    data.frame(
      metric = c("AUC_tau", "Cmax_ss", "Tmax_ss", "Cmin_ss", "Cavg_ss", "Fluctuation"),
      value  = c(auc_tau, cmax_ss, tmax_ss, cmin_ss, cav_ss, fluct)
    )
  }

  metrics_df <- reactive({
    req(sim_df())
    df <- sim_df()
    if (identical(input$metric_window, "tau")) {
      compute_metrics_tau(df, ii = input$ii, ndose = input$ndose)
    } else {
      compute_metrics_overall(df)
    }
  })

  output$metrics_tbl <- renderTable({
    req(metrics_df())
    metrics_df() |>
      transform(value = suppressWarnings(as.numeric(value)))
  }, digits = 6)

  output$download_metrics <- downloadHandler(
    filename = function() {
      w <- if (identical(input$metric_window, "tau")) "tau" else "overall"
      paste0("exposure_metrics_", w, ".csv")
    },
    content = function(file) {
      utils::write.csv(metrics_df(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
