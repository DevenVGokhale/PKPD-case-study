# Main Shiny App Entry Point
# Load required libraries
library(shiny)
library(dplyr)
library(ggplot2)
library(plotly)  # This already includes purrr as dependency
library(DT)
library(mrgsolve)
# Removed library(purrr) - it's already loaded by plotly

# Source utility functions (now local to app directory)
source("utils.R")

# Compile and load models (now local to app directory)
mod_pkpd <- mread("true_model", project = ".")

# Default parameter values
default_params <- list(
  CL = 5, VC = 50, Q2 = 10, V2 = 100, Q3 = 2, V3 = 200,
  KIN = 100, KOUT = 0.1, EC50 = 10,
  TVKG = 0.01, TVKS = 0.05, TVLAMBDA = 0.001, IC50 = 5
)

# UI
ui <- fluidPage(
  titlePanel("PK/PD Model Simulation"),
  
  fluidRow(
    column(3,
      wellPanel(
        h4("Model Parameters"),
        
        # PK Parameters
        h5("PK Parameters"),
        numericInput("CL", "Clearance (L/h)", value = default_params$CL, min = 0.1, step = 0.5),
        numericInput("VC", "Central Volume (L)", value = default_params$VC, min = 1, step = 5),
        numericInput("Q2", "Q2 (L/h)", value = default_params$Q2, min = 0.1, step = 1),
        numericInput("V2", "V2 (L)", value = default_params$V2, min = 1, step = 10),
        
        # PD Parameters
        h5("PD Parameters"),
        numericInput("KIN", "KIN", value = default_params$KIN, min = 1, step = 10),
        numericInput("KOUT", "KOUT", value = default_params$KOUT, min = 0.01, step = 0.01),
        numericInput("EC50", "EC50", value = default_params$EC50, min = 0.1, step = 1),
        
        # Dosing
        h5("Dosing"),
        numericInput("dose_amount", "Dose Amount (mg)", value = 100, min = 1, step = 10),
        radioButtons("dose_type", "Dose Type", 
                     choices = list("Flat (mg)" = "flat", "Weight-based (mg/kg)" = "weight"),
                     selected = "flat"),
        numericInput("n_subjects", "Number of Subjects", value = 10, min = 1, max = 50, step = 1),
        
        # Display options
        h5("Display Options"),
        radioButtons("display_type", "Display Type",
                     choices = list("Full Trajectory" = "full", "Sampling Window" = "sparse"),
                     selected = "full"),
        
        br(),
        actionButton("simulate", "Run Simulation", class = "btn-primary")
      )
    ),
    
    column(9,
      tabsetPanel(
        tabPanel("PK Plots", 
                 plotlyOutput("pk_plot", height = "400px"),
                 plotlyOutput("pk_obs_plot", height = "400px")),
        tabPanel("PD Plots",
                 plotlyOutput("pd_plot", height = "400px"),
                 plotlyOutput("pd_obs_plot", height = "400px")),
        tabPanel("Data Summary",
                 DT::dataTableOutput("data_summary"))
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  sim_results <- eventReactive(input$simulate, {
    
    tryCatch({
      # Create individuals
      cat("Creating individuals with n_per_arm =", input$n_subjects, "\n")
      individuals <- simulate_individuals(n_per_arm = input$n_subjects, seed = 123)
      cat("Created individuals:", nrow(individuals), "subjects\n")
      
      # Generate idata
      cat("Generating idata...\n")
      idata <- generate_idata(individuals, seed = 202)
      cat("Generated idata:", nrow(idata), "rows\n")
      
      # Update idata with user parameters
      idata$CL <- input$CL
      idata$VC <- input$VC
      idata$Q2 <- input$Q2
      idata$V2 <- input$V2
      idata$KIN <- input$KIN
      idata$KOUT <- input$KOUT
      idata$EC50 <- input$EC50
      
      # Create events
      cat("Creating event schedule...\n")
      events <- create_event_schedule(individuals, sim_days = 28)
      cat("Created events:", nrow(events), "events\n")
      
      # Adjust dose if weight-based
      if(input$dose_type == "weight") {
        events <- events |>
          left_join(individuals |> select(id, wt), by = c("ID" = "id")) |>
          mutate(amt = ifelse(evid == 1, input$dose_amount * wt, amt)) |>
          select(-wt)
      } else {
        events <- events |>
          mutate(amt = ifelse(evid == 1, input$dose_amount, amt))
      }
      
      # Run simulation
      sim_times <- if(input$display_type == "full") {
        seq(0, 28, by = 1)
      } else {
        c(0, 1, 2, 4, 8, 24, 48, 72, 168, 336, 504, 672)
      }
      
      cat("Running mrgsolve simulation...\n")
      sim_data <- mod_pkpd |>
        data_set(events) |>
        idata_set(idata) |>
        mrgsim(end = 28, delta = 1, add = sim_times) |>
        as_tibble()
      
      cat("Simulation completed:", nrow(sim_data), "rows\n")
      
      # Add observation noise
      cat("Adding observation noise...\n")
      sim_obs <- add_obs_noise(sim_data, 
                               prop_sd_cp = 0.15, add_sd_cp = 0.1,
                               prop_sd_cyt = 0.20, add_sd_cyt = 0.5)
      
      cat("Added observation noise:", nrow(sim_obs), "rows\n")
      
      result <- list(sim_data = sim_data, sim_obs = sim_obs, individuals = individuals)
      cat("Returning successful result\n")
      return(result)
      
    }, error = function(e) {
      cat("Simulation error occurred:", e$message, "\n")
      showNotification(paste("Simulation error:", e$message), type = "error", duration = 10)
      return(NULL)
    })
  })
  
  # PK plot
  output$pk_plot <- renderPlotly({
    result <- sim_results()
    if (is.null(result)) {
      cat("pk_plot: sim_results is NULL\n")
      return(NULL)
    }
    
    p <- result$sim_data |>
      filter(time > 0) |>
      ggplot(aes(x = time, y = CP, group = ID)) +
      geom_line(alpha = 0.6) +
      labs(title = "True CP Concentrations", x = "Time (hours)", y = "CP (mg/L)") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # PK observed plot
  output$pk_obs_plot <- renderPlotly({
    result <- sim_results()
    if (is.null(result)) {
      cat("pk_obs_plot: sim_results is NULL\n")
      return(NULL)
    }
    
    p <- result$sim_obs |>
      filter(time > 0) |>
      ggplot(aes(x = time, y = DV_CP, group = ID)) +
      geom_line(alpha = 0.6, color = "steelblue") +
      geom_point(alpha = 0.4, size = 0.5, color = "steelblue") +
      labs(title = "Observed CP Concentrations", x = "Time (hours)", y = "CP_obs (mg/L)") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # PD plot
  output$pd_plot <- renderPlotly({
    result <- sim_results()
    if (is.null(result)) {
      cat("pd_plot: sim_results is NULL\n")
      return(NULL)
    }
    
    p <- result$sim_data |>
      filter(time > 0) |>
      ggplot(aes(x = time, y = CYT, group = ID)) +
      geom_line(alpha = 0.6, color = "orange") +
      labs(title = "True CYT Response", x = "Time (hours)", y = "CYT") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # PD observed plot  
  output$pd_obs_plot <- renderPlotly({
    result <- sim_results()
    if (is.null(result)) {
      cat("pd_obs_plot: sim_results is NULL\n")
      return(NULL)
    }
    
    p <- result$sim_obs |>
      filter(time > 0) |>
      ggplot(aes(x = time, y = DV_CYT, group = ID)) +
      geom_line(alpha = 0.6, color = "red") +
      geom_point(alpha = 0.4, size = 0.5, color = "red") +
      labs(title = "Observed CYT Response", x = "Time (hours)", y = "CYT_obs") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # Data summary
  output$data_summary <- DT::renderDataTable({
    result <- sim_results()
    if (is.null(result)) {
      cat("data_summary: sim_results is NULL\n")
      return(data.frame(Message = "No data available. Please run simulation."))
    }
    
    result$sim_obs |>
      filter(time > 0) |>
      group_by(time) |>
      summarise(
        n_subjects = n(),
        CP_mean = round(mean(DV_CP, na.rm = TRUE), 3),
        CP_sd = round(sd(DV_CP, na.rm = TRUE), 3),
        CYT_mean = round(mean(DV_CYT, na.rm = TRUE), 3),
        CYT_sd = round(sd(DV_CYT, na.rm = TRUE), 3),
        .groups = "drop"
      )
  }, options = list(pageLength = 10))
}

# Run the application
shinyApp(ui = ui, server = server)