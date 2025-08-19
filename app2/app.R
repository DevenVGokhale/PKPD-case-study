# app.R
suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(rxode2)
})

# ----------------------------
# PK–PD MODEL (rxode2)
# ----------------------------
# PD (inhibitory on growth) with bounded resistance:
#   killFrac = KILL * (Cp/(EC50 + Cp)) * (1 - RES),  0<=KILL<=1
#   dB/dt    = KIN * B * (1 - B/KMAX) * (1 - killFrac)
#   dRES/dt  = LAMBDA * (S - RES)      [RES ∈ [0, S], S ∈ (0,1]]
model <- rxode2({
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

  # PD: inhibitory effect on growth and bounded resistance
  killFrac = KILL * (Cp / (EC50 + Cp)) * (1 - RES);
  growth   = KIN * (1 - B / KMAX);
  d/dt(B)   = B * growth * (1 - killFrac);

  d/dt(RES) = LAMBDA * (S * Cp/(EC50+Cp) - RES);

  Blog10 = log10(B + 1e-12);
})

# ----------------------------
# UI
# ----------------------------
ui <- fluidPage(
  titlePanel("3CMT + Transit PK  ➜  TB PD with Bounded Resistance (rxode2)"),
  sidebarLayout(
    sidebarPanel(
      h4("Population (typical) parameters"),
      numericInput("TVCL",  "TVCL (L/h)",  5.0,  min = 0, step = 0.1),
      numericInput("TVVc",  "TVVc (L)",   50.0, min = 1, step = 1),
      numericInput("Q2",    "Q2 (L/h)",    6.0, min = 0, step = 0.1),
      numericInput("Vp1",   "Vp1 (L)",    30.0, min = 1, step = 1),
      numericInput("Q3",    "Q3 (L/h)",    5.0, min = 0, step = 0.1),
      numericInput("Vp2",   "Vp2 (L)",    40.0, min = 1, step = 1),
      numericInput("MTT",   "MTT (h)",     2.0, min = 0.1, step = 0.1),
      sliderInput("F",      "Bioavailability (F)", min = 0, max = 1, value = 1.0, step = 0.01),

      tags$hr(),
      h4("PD (TB growth, inhibition & resistance)"),
      numericInput("TVKIN",   "TV KIN (1/h)",   0.03,  min = 0, step = 0.005),
      numericInput("KMAX",    "KMAX (CFU)",     1e8,   min = 1e3, step = 1e6),
      sliderInput("KILL",     "KILL (fractional max inhibition)", min = 0, max = 1, value = 0.6, step = 0.01),
      numericInput("EC50",    "EC50 (mg/L)",    2.0,   min = 0, step = 0.1),
      sliderInput("S",        "Resistance cap S (0–1]", min = 0.1, max = 1.0, value = 0.9, step = 0.05),
      numericInput("LAMBDA",  "LAMBDA (1/h)",   0.005, min = 0, step = 0.001),
      numericInput("B0",      "Initial bacterial load B0 (CFU)", 1e6, min = 1, step = 1e3),

      tags$hr(),
      h4("IIV (lognormal SD on ETAs)"),
      numericInput("OMEGA_CL",  "ω_CL",  0.30, min = 0, step = 0.05),
      numericInput("OMEGA_Vc",  "ω_Vc",  0.20, min = 0, step = 0.05),
      numericInput("OMEGA_KIN", "ω_KIN", 0.30, min = 0, step = 0.05),
      numericInput("Nsub", "Number of subjects", 20, min = 1, step = 1),

      tags$hr(),
      h4("Dosing"),
      numericInput("dose", "Oral dose (mg)",        600,  min = 0, step = 10),
      numericInput("ii",   "Dosing interval τ (h)", 24,   min = 0.1, step = 1),
      numericInput("ndose","Number of doses",       10,   min = 1, step = 1),

      tags$hr(),
      h4("Sampling grid"),
      numericInput("t_end", "Simulation end (h)", 240, min = 1, step = 1),
      numericInput("dt",    "Δt (h)",            0.5,  min = 0.01, step = 0.01),
      actionButton("apply_grid", "Apply sampling grid", class = "btn-info"),

      tags$hr(),
      h4("Residual error (for plotted/metric outputs)"),
      numericInput("prop_sd_cp", "Prop SD (Cp)", 0.10, min = 0, step = 0.01),
      numericInput("add_sd_cp",  "Add SD (Cp, mg/L)", 0.1, min = 0, step = 0.01),
      numericInput("add_sd_b",   "Add SD (log10 B)", 0.05, min = 0, step = 0.01),

      tags$hr(),
      actionButton("run", "Simulate", class = "btn-primary")
    ),
    mainPanel(
      uiOutput("grid_summary"),
      tabsetPanel(
        tabPanel("PK curve",
          plotOutput("pk_plot", height = 380)
        ),
        tabPanel("Exposure metrics",
          radioButtons(
            "metric_window",
            label = "Exposure window",
            choices = c("Overall (0–t_end)" = "overall",
                        "Last dosing interval (≈ steady-state)" = "tau"),
            selected = "overall",
            inline = TRUE
          ),
          tableOutput("metrics_tbl"),
          downloadButton("download_metrics", "Download metrics (CSV)")
        ),
        tabPanel("PD (TB) curve",
          plotOutput("pd_plot", height = 380),
          tableOutput("pd_tbl")
        ),
        tabPanel("Resistance diagnostics",
          plotOutput("res_plot", height = 280),
          plotOutput("kill_plot", height = 280)
        )
      )
    )
  )
)

# ----------------------------
# SERVER
# ----------------------------
server <- function(input, output, session) {

  # ---- hold applied sampling grid times ----
  grid_times <- reactiveVal(seq(0, 240, by = 0.5))

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
    HTML(sprintf(
      "<strong>Applied sampling grid:</strong> n=%d, start=%.3f h, end=%.3f h, dt≈%.3f h",
      n, head(times, 1), tail(times, 1),
      ifelse(n > 1, diff(times)[1], NA_real_)
    ))
  })

  # ---- helpers ----
  trapz <- function(x, y) {
    if (length(x) < 2) return(NA_real_)
    sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  }

  # Simulate id-level ETAs (lognormal) and return id-param frame
  make_id_params <- reactive({
    N <- input$Nsub
    set.seed(123)
    eta_cl  <- rnorm(N, 0, input$OMEGA_CL)
    eta_vc  <- rnorm(N, 0, input$OMEGA_Vc)
    eta_kin <- rnorm(N, 0, input$OMEGA_KIN)

    tibble(
      id  = 1:N,
      CL  = input$TVCL * exp(eta_cl),
      Vc  = input$TVVc * exp(eta_vc),
      Vp1 = input$Vp1,
      Vp2 = input$Vp2,
      Q2  = input$Q2,
      Q3  = input$Q3,
      MTT = input$MTT,
      F   = input$F,

      # PD with IIV on KIN; inhibitory params & resistance cap
      KIN    = input$TVKIN * exp(eta_kin),
      KMAX   = input$KMAX,
      KILL   = input$KILL,   # renamed from IMAX
      EC50   = input$EC50,
      S      = input$S,
      LAMBDA = input$LAMBDA
    )
  })

  # ---- events as a single data.frame (NO et(), NO "+") ----
  make_events <- reactive({
    N <- input$Nsub
    times <- grid_times()

    dosing <- tibble(
      id   = 1:N,
      time = 0,
      evid = 1L,
      cmt  = "A0",                      # dose into A0 (gut)
      amt  = input$F * input$dose,
      ii   = input$ii,
      addl = input$ndose - 1L
    )

    obs <- tibble(
      id   = rep(1:N, each = length(times)),
      time = rep(times, times = N),
      evid = 0L,
      cmt  = NA_character_,
      amt  = 0,
      ii   = 0,
      addl = 0
    )

    bind_rows(dosing, obs) |> arrange(id, time, desc(evid))
  })

  # Initial states (B0 now user-controlled)
  init_states <- reactive({
    c(A0 = 0, TR1 = 0, TR2 = 0, TR3 = 0,
      A1 = 0, A2 = 0, A3 = 0,
      B  = input$B0,  # was fixed at 1e6
      RES = 0)
  })

  # ---- simulate on click ----
  sim_df <- eventReactive(input$run, {
    rxSolve(
      model,
      params = make_id_params(),
      events = make_events(),
      inits  = init_states(),
      addDosing = TRUE,
      rtol = 1e-8, atol = 1e-8
    ) |>
      as.data.frame()
  }, ignoreInit = TRUE)

  # ---- add residual error for outputs (doesn't feed back into ODEs) ----
  sim_with_noise <- reactive({
    df <- sim_df()
    validate(need(!is.null(df), "Run a simulation first."))
    if (nrow(df) == 0) return(df)

    # Cp residuals
    if (input$prop_sd_cp > 0 || input$add_sd_cp > 0) {
      set.seed(202)
      df$Cp_obs <- pmax(
        0,
        df$Cp * (1 + rnorm(nrow(df), 0, input$prop_sd_cp)) +
          rnorm(nrow(df), 0, input$add_sd_cp)
      )
    } else {
      df$Cp_obs <- df$Cp
    }

    # log10(B) residuals
    if (input$add_sd_b > 0) {
      set.seed(203)
      df$Blog10_obs <- df$Blog10 + rnorm(nrow(df), 0, input$add_sd_b)
    } else {
      df$Blog10_obs <- df$Blog10
    }

    df
  })

  # ----------------------------
  # PK PLOT
  # ----------------------------
  output$pk_plot <- renderPlot({
    df <- sim_with_noise()
    ggplot(df, aes(time, Cp_obs, group = id)) +
      geom_line(alpha = 0.5) +
      labs(x = "Time (h)", y = "Concentration (mg/L)",
           title = "PK: Cp over time (population)") +
      theme_bw()
  })

  # ----------------------------
  # EXPOSURE METRICS (from noise-free Cp)
  # ----------------------------
  compute_metrics_overall <- function(df) {
    df |>
      group_by(id) |>
      summarise(
        AUC  = trapz(time, Cp),
        Cmax = max(Cp),
        Tmax = time[which.max(Cp)[1]],
        .groups = "drop"
      )
  }

  compute_metrics_tau <- function(df, ii, ndose) {
    t0 <- (ndose - 1) * ii
    t1 <- ndose * ii
    df |>
      filter(time >= t0, time <= t1) |>
      group_by(id) |>
      summarise(
        AUC_tau = trapz(time, Cp),
        Cavg_ss = AUC_tau / ii,
        Cmax_ss = max(Cp),
        Tmax_ss = time[which.max(Cp)[1]],
        Cmin_ss = min(Cp),
        Fluctuation = (Cmax_ss - Cmin_ss) / Cavg_ss,
        .groups = "drop"
      )
  }

  metrics_df <- reactive({
    df <- sim_df()
    validate(need(!is.null(df) && nrow(df) > 1, "Run a simulation first."))
    if (identical(input$metric_window, "tau")) {
      validate(need(max(df$time) >= input$ndose * input$ii,
                    "Increase t_end to cover the last dosing interval."))
      compute_metrics_tau(df, input$ii, input$ndose)
    } else {
      compute_metrics_overall(df)
    }
  })

  output$metrics_tbl <- renderTable({
    as.data.frame(metrics_df())
  }, digits = 6)

  output$download_metrics <- downloadHandler(
    filename = function() {
      w <- if (identical(input$metric_window, "tau")) "tau" else "overall"
      paste0("exposure_metrics_", w, ".csv")
    },
    content = function(file) {
      m <- metrics_df()
      req(!is.null(m), nrow(m) > 0)
      utils::write.csv(m, file, row.names = FALSE)
    }
  )

  # ----------------------------
  # PD (TB) PLOT + METRICS
  # ----------------------------
  output$pd_plot <- renderPlot({
    df <- sim_with_noise()
    ggplot(df, aes(time, Blog10_obs, group = id)) +
      geom_line(alpha = 0.5) +
      labs(x = "Time (h)", y = "log10 CFU",
           title = "PD: TB bacterial load (log10 CFU)") +
      theme_bw()
  })

  output$pd_tbl <- renderTable({
    df <- sim_df()
    validate(need(!is.null(df), "Run a simulation first."))
    out <- df |>
      group_by(id) |>
      summarise(
        B0_log10   = first(Blog10),
        Bmin_log10 = min(Blog10),
        t_Bmin     = time[which.min(Blog10)[1]],
        Bend_log10 = last(Blog10),
        .groups = "drop"
      )
    as.data.frame(out)
  }, digits = 4)

  # ----------------------------
  # Diagnostics to visualize resistance effect
  # ----------------------------
  output$res_plot <- renderPlot({
    df <- sim_df(); req(nrow(df) > 0)
    ggplot(df, aes(time, RES, group = id)) +
      geom_line(alpha = 0.5) +
      labs(x = "Time (h)", y = "RES", title = "Resistance over time") +
      theme_bw()
  })

  output$kill_plot <- renderPlot({
    df <- sim_df(); req(nrow(df) > 0)
    df$killFrac <- with(df, input$KILL * (Cp / (input$EC50 + Cp)) * (1 - RES))
    ggplot(df, aes(time, killFrac, group = id)) +
      geom_line(alpha = 0.5) +
      labs(x = "Time (h)", y = "killFrac", title = "Instantaneous kill fraction") +
      theme_bw()
  })
}

shinyApp(ui, server)
