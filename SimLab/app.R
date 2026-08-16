library(shiny)
library(shinydashboard)
library(ggplot2)
library(data.table)
library(plotly)
library(DT)

# ============================================================
# 1. HELPER FUNCTIONS
# ============================================================

# ------------------------------------------------------------
# 1.1 Formatting helpers
# ------------------------------------------------------------

format_whole_number <- function(x) {
  format(
    round(x, 0),
    big.mark = ",",
    scientific = FALSE
  )
}

format_runtime_rate <- function(n_sims, runtime_seconds) {
  if (runtime_seconds > 0) {
    format(
      round(n_sims / runtime_seconds),
      big.mark = ",",
      scientific = FALSE
    )
  } else {
    "NA"
  }
}

# ------------------------------------------------------------
# 1.2 Reusable user interface helpers
# ------------------------------------------------------------

creator_watermark <- function() {
  fluidRow(
    column(
      width = 12,
      div(
        style = paste(
          "margin-top: 20px;",
          "padding: 18px 0 12px 0;",
          "border-top: 1px solid #d2d6de;",
          "text-align: center;",
          "color: #777;",
          "font-size: 13px;"
        ),
        tags$span("Created by Benjamin Bilyak"),
        tags$span(" \u00b7 "),
        tags$a(
          href = "https://github.com/Benjamin-Bilyak/R-Stats-Projects",
          target = "_blank",
          rel = "noopener noreferrer",
          icon("github"),
          "View project on GitHub"
        )
      )
    )
  )
}

# ------------------------------------------------------------
# 1.3 Histogram helpers
# ------------------------------------------------------------

get_histogram_breaks <- function(loss_values, bins) {
  min_loss <- min(loss_values)
  max_loss <- max(loss_values)
  
  if (min_loss == max_loss) {
    return(pretty(loss_values, n = bins))
  }
  
  seq(
    min_loss,
    max_loss,
    length.out = bins + 1
  )
}

create_loss_histogram <- function(
    sim_data,
    bins,
    fill_colour,
    border_colour,
    show_density,
    density_colour
) {
  
  histogram_border <- if (
    border_colour == "none"
  ) fill_colour else border_colour
  
  histogram_breaks <- get_histogram_breaks(
    sim_data$Loss,
    bins
  )
  
  p <- ggplot(sim_data, aes(x = Loss)) +
    geom_histogram(
      aes(y = after_stat(density)),
      breaks = histogram_breaks,
      fill = fill_colour,
      colour = histogram_border
    ) +
    labs(
      title = "Monte Carlo Aggregate Loss Distribution",
      x = "Aggregate Loss",
      y = "Density"
    )
  
  if (show_density) {
    p <- p +
      geom_density(
        colour = density_colour,
        linewidth = 1
      )
  }
  
  p
}

get_selected_bin <- function(sim_data, bins, clicked_x) {
  histogram_breaks <- get_histogram_breaks(
    sim_data$Loss,
    bins
  )
  
  bin_index <- findInterval(
    clicked_x,
    histogram_breaks,
    all.inside = TRUE
  )
  
  lower_bound <- histogram_breaks[bin_index]
  upper_bound <- histogram_breaks[bin_index + 1]
  
  is_last_bin <- bin_index == length(histogram_breaks) - 1
  
  if (is_last_bin) {
    bin_data <- sim_data[
      sim_data$Loss >= lower_bound &
        sim_data$Loss <= upper_bound,
    ]
  } else {
    bin_data <- sim_data[
      sim_data$Loss >= lower_bound &
        sim_data$Loss < upper_bound,
    ]
  }
  
  list(
    data = bin_data,
    lower_bound = lower_bound,
    upper_bound = upper_bound
  )
}

# ------------------------------------------------------------
# 1.4 Simulation engine
# ------------------------------------------------------------

run_aggregate_loss_simulation <- function(
    n_sims,
    claim_frequency,
    severity_mean,
    severity_cv
) {
  
  start_time <- Sys.time()
  
  severity_variance <- (
    severity_cv * severity_mean
  )^2
  
  # 1. Draw the number of claims in each simulation.
  claim_counts <- rpois(
    n = n_sims,
    lambda = claim_frequency
  )
  
  # 2. Count the total number of individual claims across all simulations.
  total_claims <- sum(claim_counts)
  
  # 3. Handle the edge case where every simulation has zero claims.
  if (total_claims == 0) {
    aggregate_losses <- numeric(n_sims)
  } else {
    
    # 4. Draw all individual claim severities at once.
    claim_severities <- rgamma(
      n = total_claims,
      shape = severity_mean^2 / severity_variance,
      scale = severity_variance / severity_mean
    )
    
    # 5. Match each individual claim to its simulation number.
    simulation_id <- rep(
      seq_len(n_sims),
      claim_counts
    )
    
    # 6. Add individual claims back up to simulation-level aggregate losses.
    aggregate_losses <- numeric(n_sims)
    
    loss_sums <- rowsum(
      claim_severities,
      simulation_id
    )
    
    aggregate_losses[as.integer(rownames(loss_sums))] <- loss_sums[, 1]
  }
  
  runtime <- as.numeric(
    difftime(Sys.time(), start_time, units = "secs")
  )
  
  expected_mean <- claim_frequency * severity_mean
  simulated_mean <- mean(aggregate_losses)
  mean_difference <- simulated_mean - expected_mean
  mean_difference_pct <- (
    mean_difference / expected_mean
  ) * 100
  
  expected_variance <- claim_frequency * (
    severity_variance + severity_mean^2
  )
  
  expected_sd <- sqrt(expected_variance)
  simulated_sd <- sd(aggregate_losses)
  sd_difference <- simulated_sd - expected_sd
  sd_difference_pct <- (
    sd_difference / expected_sd
  ) * 100
  
  list(
    data = data.frame(
      Simulation = seq_len(n_sims),
      Loss = aggregate_losses
    ),
    runtime = runtime,
    severity_variance = severity_variance,
    expected_mean = expected_mean,
    simulated_mean = simulated_mean,
    mean_difference = mean_difference,
    mean_difference_pct = mean_difference_pct,
    expected_sd = expected_sd,
    simulated_sd = simulated_sd,
    sd_difference = sd_difference,
    sd_difference_pct = sd_difference_pct
  )
}

# ------------------------------------------------------------
# 1.5 Galton board simulation engine
# ------------------------------------------------------------

run_galton_simulation <- function(
    balls,
    rows,
    p
) {
  
  start_time <- Sys.time()
  
  # Vectorized simulation:
  # Each ball's number of right moves is drawn directly.
  right_moves <- rbinom(
    n = balls,
    size = rows,
    prob = p
  )
  
  # Centre the board around zero.
  final_position <- right_moves - rows / 2
  
  runtime <- as.numeric(
    difftime(Sys.time(), start_time, units = "secs")
  )
  
  expected_mean <- rows * p - rows / 2
  simulated_mean <- mean(final_position)
  mean_difference <- simulated_mean - expected_mean
  
  expected_sd <- sqrt(rows * p * (1 - p))
  simulated_sd <- sd(final_position)
  sd_difference <- simulated_sd - expected_sd
  
  list(
    data = data.frame(
      Ball = seq_len(balls),
      RightMoves = right_moves,
      FinalPosition = final_position
    ),
    runtime = runtime,
    expected_mean = expected_mean,
    simulated_mean = simulated_mean,
    mean_difference = mean_difference,
    expected_sd = expected_sd,
    simulated_sd = simulated_sd,
    sd_difference = sd_difference
  )
}
# ------------------------------------------------------------
# 1.6 Galton board histogram
# ------------------------------------------------------------

create_galton_histogram <- function(
    galton_data,
    rows,
    p,
    fill_colour,
    border_colour,
    show_theory,
    theory_colour
) {
  
  slot_counts <- as.data.frame(
    table(galton_data$FinalPosition)
  )
  
  names(slot_counts) <- c(
    "FinalPosition",
    "Count"
  )
  
  slot_counts$FinalPosition <- as.numeric(
    as.character(slot_counts$FinalPosition)
  )
  
  slot_counts$SimulatedProbability <- slot_counts$Count /
    sum(slot_counts$Count)
  
  possible_right_moves <- 0:rows
  possible_positions <- possible_right_moves - rows / 2
  
  theory_data <- data.frame(
    FinalPosition = possible_positions,
    TheoreticalProbability = dbinom(
      possible_right_moves,
      size = rows,
      prob = p
    )
  )
  
  p_plot <- ggplot(
    slot_counts,
    aes(
      x = FinalPosition,
      y = SimulatedProbability
    )
  ) +
    geom_col(
      fill = fill_colour,
      colour = border_colour,
      linewidth = 0.8
    ) +
    labs(
      title = "Galton Board Final Position Distribution",
      x = "Final Position",
      y = "Probability"
    )
  
  if (show_theory) {
    p_plot <- p_plot +
      geom_point(
        data = theory_data,
        aes(
          x = FinalPosition,
          y = TheoreticalProbability
        ),
        colour = theory_colour,
        size = 2
      ) +
      geom_line(
        data = theory_data,
        aes(
          x = FinalPosition,
          y = TheoreticalProbability
        ),
        colour = theory_colour,
        linewidth = 1
      )
  }
  
  p_plot
  
}

get_selected_galton_bar <- function(galton_data, clicked_x) {
  
  possible_positions <- sort(
    unique(galton_data$FinalPosition)
  )
  
  selected_position <- possible_positions[
    which.min(abs(possible_positions - clicked_x))
  ]
  
  selected_data <- galton_data[
    galton_data$FinalPosition == selected_position,
  ]
  
  list(
    data = selected_data,
    final_position = selected_position
  )
}

# ------------------------------------------------------------
# 1.7 Normal distribution helpers
# ------------------------------------------------------------

generate_normal_data <- function(
    mean,
    sd,
    display_sd = 10,
    points = 2501
) {
  stopifnot(
    is.finite(mean),
    is.finite(sd),
    sd > 0,
    is.finite(display_sd),
    display_sd > 0,
    points >= 2
  )
  
  x_values <- seq(
    mean - display_sd * sd,
    mean + display_sd * sd,
    length.out = points
  )
  
  data.frame(
    x = x_values,
    density = dnorm(
      x_values,
      mean = mean,
      sd = sd
    )
  )
}

calculate_normal_probability <- function(
    calculation_type,
    mean,
    sd,
    x = NULL,
    lower = NULL,
    upper = NULL
) {
  if (calculation_type == "lower") {
    return(pnorm(x, mean = mean, sd = sd))
  }
  
  if (calculation_type == "upper") {
    return(pnorm(x, mean = mean, sd = sd, lower.tail = FALSE))
  }
  
  if (calculation_type == "between") {
    if (lower > upper) {
      stop("The lower boundary must not exceed the upper boundary.")
    }
    
    return(
      pnorm(upper, mean = mean, sd = sd) -
        pnorm(lower, mean = mean, sd = sd)
    )
  }
  
  stop("Unknown Normal probability calculation type.")
}

# ------------------------------------------------------------
# 1.8 Binomial distribution helpers
# ------------------------------------------------------------

generate_binomial_data <- function(
    trials,
    probability
) {
  stopifnot(
    is.finite(trials),
    trials >= 1,
    trials == round(trials),
    is.finite(probability),
    probability >= 0,
    probability <= 1
  )
  
  successes <- 0:trials
  
  data.frame(
    x = successes,
    probability = dbinom(
      successes,
      size = trials,
      prob = probability
    )
  )
}

calculate_binomial_probability <- function(
    calculation_type,
    trials,
    probability,
    x = NULL,
    lower = NULL,
    upper = NULL
) {
  if (calculation_type == "lower") {
    return(
      pbinom(
        floor(x),
        size = trials,
        prob = probability
      )
    )
  }
  
  if (calculation_type == "upper") {
    return(
      pbinom(
        ceiling(x) - 1,
        size = trials,
        prob = probability,
        lower.tail = FALSE
      )
    )
  }
  
  if (calculation_type == "between") {
    if (lower > upper) {
      stop("The lower boundary must not exceed the upper boundary.")
    }
    
    first_integer <- ceiling(lower)
    last_integer <- floor(upper)
    
    if (first_integer > last_integer) {
      return(0)
    }
    
    return(
      pbinom(
        last_integer,
        size = trials,
        prob = probability
      ) -
        pbinom(
          first_integer - 1,
          size = trials,
          prob = probability
        )
    )
  }
  
  stop("Unknown Binomial probability calculation type.")
}

# ------------------------------------------------------------
# 1.9 Poisson distribution helpers
# ------------------------------------------------------------

generate_poisson_data <- function(
    lambda,
    cumulative_probability = 0.999999
) {
  stopifnot(
    is.finite(lambda),
    lambda >= 0,
    is.finite(cumulative_probability),
    cumulative_probability > 0,
    cumulative_probability < 1
  )
  
  maximum_event_count <- qpois(
    cumulative_probability,
    lambda = lambda
  )
  
  event_counts <- 0:maximum_event_count
  
  data.frame(
    x = event_counts,
    probability = dpois(
      event_counts,
      lambda = lambda
    )
  )
}

calculate_poisson_probability <- function(
    calculation_type,
    lambda,
    x = NULL,
    lower = NULL,
    upper = NULL
) {
  if (calculation_type == "lower") {
    return(
      ppois(
        floor(x),
        lambda = lambda
      )
    )
  }
  
  if (calculation_type == "upper") {
    return(
      ppois(
        ceiling(x) - 1,
        lambda = lambda,
        lower.tail = FALSE
      )
    )
  }
  
  if (calculation_type == "between") {
    if (lower > upper) {
      stop("The lower boundary must not exceed the upper boundary.")
    }
    
    first_integer <- max(
      0,
      ceiling(lower)
    )
    
    last_integer <- floor(upper)
    
    if (first_integer > last_integer) {
      return(0)
    }
    
    return(
      ppois(
        last_integer,
        lambda = lambda
      ) -
        ppois(
          first_integer - 1,
          lambda = lambda
        )
    )
  }
  
  stop("Unknown Poisson probability calculation type.")
}

# ------------------------------------------------------------
# 1.10 Exponential distribution helpers
# ------------------------------------------------------------

generate_exponential_data <- function(
    rate,
    cumulative_probability = 0.999999,
    points = 2501
) {
  stopifnot(
    is.finite(rate),
    rate > 0,
    is.finite(cumulative_probability),
    cumulative_probability > 0,
    cumulative_probability < 1,
    points >= 2
  )
  
  maximum_x <- qexp(
    cumulative_probability,
    rate = rate
  )
  
  x_values <- seq(
    0,
    maximum_x,
    length.out = points
  )
  
  data.frame(
    x = x_values,
    density = dexp(
      x_values,
      rate = rate
    )
  )
}

calculate_exponential_probability <- function(
    calculation_type,
    rate,
    x = NULL,
    lower = NULL,
    upper = NULL
) {
  if (calculation_type == "lower") {
    return(
      pexp(
        x,
        rate = rate
      )
    )
  }
  
  if (calculation_type == "upper") {
    return(
      pexp(
        x,
        rate = rate,
        lower.tail = FALSE
      )
    )
  }
  
  if (calculation_type == "between") {
    if (lower > upper) {
      stop("The lower boundary must not exceed the upper boundary.")
    }
    
    return(
      pexp(
        upper,
        rate = rate
      ) -
        pexp(
          lower,
          rate = rate
        )
    )
  }
  
  stop("Unknown Exponential probability calculation type.")
}

# ------------------------------------------------------------
# 1.11 Gamma distribution helpers
# ------------------------------------------------------------

generate_gamma_data <- function(
    shape,
    rate,
    cumulative_probability = 0.999999,
    points = 2501
) {
  stopifnot(
    is.finite(shape),
    shape > 0,
    is.finite(rate),
    rate > 0,
    is.finite(cumulative_probability),
    cumulative_probability > 0,
    cumulative_probability < 1,
    points >= 2
  )
  
  maximum_x <- qgamma(
    cumulative_probability,
    shape = shape,
    rate = rate
  )
  
  minimum_x <- if (shape < 1) {
    max(
      qgamma(
        1 - cumulative_probability,
        shape = shape,
        rate = rate
      ),
      maximum_x / (points - 1)
    )
  } else {
    0
  }
  
  x_values <- seq(
    minimum_x,
    maximum_x,
    length.out = points
  )
  
  data.frame(
    x = x_values,
    density = dgamma(
      x_values,
      shape = shape,
      rate = rate
    )
  )
}

calculate_gamma_probability <- function(
    calculation_type,
    shape,
    rate,
    x = NULL,
    lower = NULL,
    upper = NULL
) {
  if (calculation_type == "lower") {
    return(
      pgamma(
        x,
        shape = shape,
        rate = rate
      )
    )
  }
  
  if (calculation_type == "upper") {
    return(
      pgamma(
        x,
        shape = shape,
        rate = rate,
        lower.tail = FALSE
      )
    )
  }
  
  if (calculation_type == "between") {
    if (lower > upper) {
      stop("The lower boundary must not exceed the upper boundary.")
    }
    
    return(
      pgamma(
        upper,
        shape = shape,
        rate = rate
      ) -
        pgamma(
          lower,
          shape = shape,
          rate = rate
        )
    )
  }
  
  stop("Unknown Gamma probability calculation type.")
}

# ------------------------------------------------------------
# 1.12 Lognormal distribution helpers
# ------------------------------------------------------------

generate_lognormal_data <- function(
    meanlog,
    sdlog,
    cumulative_probability = 0.999999,
    points = 2501
) {
  stopifnot(
    is.finite(meanlog),
    is.finite(sdlog),
    sdlog > 0,
    is.finite(cumulative_probability),
    cumulative_probability > 0,
    cumulative_probability < 1,
    points >= 3
  )
  
  lower_probability <- 1 - cumulative_probability
  
  log_minimum_x <- qnorm(
    lower_probability,
    mean = meanlog,
    sd = sdlog
  )
  
  log_maximum_x <- qnorm(
    cumulative_probability,
    mean = meanlog,
    sd = sdlog
  )
  
  if (log_maximum_x > log(.Machine$double.xmax)) {
    stop(
      "These parameters produce values too large to display. ",
      "Reduce the log-mean or log-standard deviation."
    )
  }
  
  minimum_positive_x <- exp(
    max(
      log_minimum_x,
      log(.Machine$double.xmin)
    )
  )
  
  maximum_x <- exp(log_maximum_x)
  
  x_values <- c(
    0,
    exp(
      seq(
        log(minimum_positive_x),
        log(maximum_x),
        length.out = points - 1
      )
    )
  )
  
  data.frame(
    x = x_values,
    density = dlnorm(
      x_values,
      meanlog = meanlog,
      sdlog = sdlog
    )
  )
}

calculate_lognormal_probability <- function(
    calculation_type,
    meanlog,
    sdlog,
    x = NULL,
    lower = NULL,
    upper = NULL
) {
  if (calculation_type == "lower") {
    return(
      plnorm(
        x,
        meanlog = meanlog,
        sdlog = sdlog
      )
    )
  }
  
  if (calculation_type == "upper") {
    return(
      plnorm(
        x,
        meanlog = meanlog,
        sdlog = sdlog,
        lower.tail = FALSE
      )
    )
  }
  
  if (calculation_type == "between") {
    if (lower > upper) {
      stop("The lower boundary must not exceed the upper boundary.")
    }
    
    return(
      plnorm(
        upper,
        meanlog = meanlog,
        sdlog = sdlog
      ) -
        plnorm(
          lower,
          meanlog = meanlog,
          sdlog = sdlog
        )
    )
  }
  
  stop("Unknown Lognormal probability calculation type.")
}

# ------------------------------------------------------------
# 1.13 Pi estimator helpers
# ------------------------------------------------------------

calculate_pi_estimate <- function(n_trapeziums) {
  x_values <- seq(
    from = 0,
    to = 1,
    length.out = n_trapeziums + 1
  )
  
  y_values <- sqrt(
    pmax(0, 1 - x_values^2)
  )
  
  trapezium_widths <- diff(x_values)
  
  trapezium_areas <- trapezium_widths * (
    head(y_values, -1) + tail(y_values, -1)
  ) / 2
  
  quarter_pi_estimate <- sum(trapezium_areas)
  pi_estimate <- 4 * quarter_pi_estimate
  
  trapezium_polygons <- data.frame(
    x = as.vector(
      rbind(
        head(x_values, -1),
        tail(x_values, -1),
        tail(x_values, -1),
        head(x_values, -1),
        head(x_values, -1),
        rep(NA_real_, n_trapeziums)
      )
    ),
    y = as.vector(
      rbind(
        rep(0, n_trapeziums),
        rep(0, n_trapeziums),
        tail(y_values, -1),
        head(y_values, -1),
        rep(0, n_trapeziums),
        rep(NA_real_, n_trapeziums)
      )
    )
  )
  
  curve_x <- seq(
    from = 0,
    to = 1,
    length.out = 1001
  )
  
  curve_data <- data.frame(
    x = curve_x,
    y = sqrt(pmax(0, 1 - curve_x^2))
  )
  
  list(
    quarter_pi_estimate = quarter_pi_estimate,
    pi_estimate = pi_estimate,
    trapezium_polygons = trapezium_polygons,
    curve_data = curve_data
  )
}

# ============================================================
# 2. USER INTERFACE
# ============================================================

ui <- dashboardPage(
  title = "SimLab",
  
  # ------------------------------------------------------------
  # 2.1 Header and sidebar
  # ------------------------------------------------------------
  
  dashboardHeader(title = "SimLab"),
  
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Aggregate Loss Simulator", tabName = "aggregateLoss", icon = icon("chart-bar")),
      menuItem("Galton Board", tabName = "galton", icon = icon("circle")),
      menuItem("Distribution Explorer", tabName = "distributions", icon = icon("chart-line")),
      menuItem("Pi Estimator", tabName = "piEstimator", icon = icon("calculator"))
    )
  ),
  
  dashboardBody(
    tabItems(
      
      # ------------------------------------------------------------
      # 2.2 Home tab
      # ------------------------------------------------------------
      
      # ------------------------------------------------------------
      # 2.2 Home tab
      # ------------------------------------------------------------
      
      tabItem(
        tabName = "home",
        
        fluidRow(
          box(
            title = "Welcome to SimLab",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            h2("A Collection of Interactive Simulations"),
            
            p("SimLab is a growing collection of interactive simulations built in R Shiny. Each module turns a mathematical, statistical, or modelling idea into a visual experiment that users can control and explore."),
            
            p("The collection currently includes an aggregate loss simulation, a Galton board, a distribution explorer, and a numerical estimator for pi."),
            
            br(),
            
            tags$ul(
              tags$li("A growing library of independent simulation modules"),
              tags$li("Adjustable parameters for hands-on experimentation"),
              tags$li("Reproducible results using controllable random seeds"),
              tags$li("Simulated results compared with theoretical models"),
              tags$li("Interactive graphs, diagnostics, and data tables"),
              tags$li("Downloadable outputs for further analysis")
            )
          )
        ),
        
        fluidRow(
          valueBox(
            "4",
            "Available Modules",
            icon = icon("flask"),
            color = "blue",
            width = 4
          ),
          
          valueBox(
            "Seed Control",
            "Reproducible Runs",
            icon = icon("redo"),
            color = "green",
            width = 4
          ),
          
          valueBox(
            "Growing Collection",
            "More Modules Planned",
            icon = icon("plus-circle"),
            color = "purple",
            width = 4
          )
        ),
        
        fluidRow(
          box(
            title = "Aggregate Loss Simulator",
            width = 6,
            status = "info",
            solidHeader = TRUE,
            
            p("The aggregate loss simulator models total insurance losses as:"),
            
            withMathJax("$$S = \\sum_{i=1}^{N} X_i$$"),
            
            p("where N is the number of claims and X is the size of each claim."),
            
            p("In the model, claim frequency follows a Poisson distribution and claim severity follows a Gamma distribution."),
            
            br(),
            
            actionButton(
              "goToSimulator",
              "Open Aggregate Loss Simulator",
              icon = icon("play")
            )
          ),
          
          box(
            title = "Galton Board Simulator",
            width = 6,
            status = "success",
            solidHeader = TRUE,
            
            p("The Galton board models repeated left/right movements. Each ball has a probability p of moving right at each row."),
            
            withMathJax("$$X \\sim \\text{Binomial}(n, p)$$"),
            
            p("The final position is centred around zero, allowing the simulated distribution to be compared with the theoretical binomial distribution."),
            
            br(),
            
            actionButton(
              "goToGalton",
              "Open Galton Board",
              icon = icon("play")
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Distribution Explorer",
            width = 6,
            status = "warning",
            solidHeader = TRUE,
            
            p("Explore six probability distributions: Normal, Binomial, Poisson, Exponential, Gamma, and Lognormal."),
            
            p("Adjust each distribution's parameters, examine its graph and properties, and calculate probabilities interactively."),
            
            br(),
            
            actionButton(
              "goToDistributions",
              "Open Distribution Explorer",
              icon = icon("chart-line")
            )
          ),
          
          box(
            title = "Pi Estimator",
            width = 6,
            status = "primary",
            solidHeader = TRUE,
            
            p("The Pi Estimator uses the trapezium rule to approximate the area under the upper-right quarter of the unit circle."),
            
            withMathJax("$$\\int_0^1 \\sqrt{1-x^2}\\,dx = \\frac{\\pi}{4}$$"),
            
            p("Change the number of trapeziums and watch the numerical estimate approach pi."),
            
            br(),
            
            actionButton(
              "goToPiEstimator",
              "Open Pi Estimator",
              icon = icon("calculator")
            )
          )
        ),
        
        fluidRow(
          box(
            title = "More Development Coming",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            
            p("SimLab is an expanding collection. Future development can add more simulations, numerical methods, probability models, and quantitative experiments."),
            
            p("Each new module will continue the same approach: adjustable inputs, clear visual outputs, and direct comparison between mathematical theory and computation.")
          )
        ),
        
        creator_watermark()
      ),
      
      # ------------------------------------------------------------
      # 2.3 Aggregate loss simulator tab
      # ------------------------------------------------------------
      
      tabItem(
        tabName = "aggregateLoss",
        
        # Main histogram and controls
        fluidRow(
          box(
            title = "Monte Carlo Aggregate Loss Distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 9,
            plotlyOutput("lossHistogram", height = "600px")
          ),
          
          box(
            title = "Model Controls",
            width = 3,
            
            sliderInput(
              "claimFrequency",
              "Frequency of Claims: lambda",
              min = 1,
              max = 10,
              value = 5
            ),
            
            sliderInput(
              "histogramBins",
              "Number of Histogram Bins",
              min = 1,
              max = 300,
              value = 50
            ),
            
            sliderInput(
              "severityMean",
              "Expected Severity: Mean",
              min = 1000000,
              max = 10000000,
              value = 5000000
            ),
            
            sliderInput(
              "severityCV",
              "Severity Coefficient of Variation",
              min = 0.05,
              max = 3,
              value = 1,
              step = 0.05
            ),
            
            selectInput(
              "histogramFillColour",
              "Histogram Colour",
              choices = c(
                "Blue" = "steelblue",
                "Green" = "green",
                "Red" = "red",
                "Orange" = "orange",
                "Purple" = "purple"
              )
            ),
            
            selectInput(
              "histogramBorderColour",
              "Border Colour",
              choices = c(
                "Black" = "black",
                "Green" = "green",
                "Red" = "red",
                "Orange" = "orange",
                "Purple" = "purple",
                "Blue" = "steelblue",
                "None" = "none"
              )
            ),
            
            checkboxInput(
              "showDensity",
              "Show Density Curve",
              value = TRUE
            ),
            
            selectInput(
              "densityColour",
              "Density Curve Colour",
              choices = c(
                "Red" = "red",
                "Green" = "green",
                "Blue" = "steelblue",
                "Orange" = "orange",
                "Purple" = "purple",
                "Black" = "black"
              )
            )
          )
        ),
        
        # Simulation controls and headline output
        fluidRow(
          box(
            title = "Run Simulation",
            width = 5,
            
            actionButton(
              "runSimulation",
              "Run Simulation",
              icon = icon("play")
            ),
            
            br(),
            br(),
            
            downloadButton(
              "downloadHistogram",
              "Download Histogram"
            ),
            br(),
            br(),
            
            downloadButton(
              "downloadLossCSV",
              "Download Full Simulation CSV"
            ),
            br(),
            br(),
            
            textOutput("simulationCountText"),
            textOutput("lossSummaryText"),
            textOutput("selectedBinSummary")
          ),
          
          box(
            title = "Simulation Settings",
            width = 7,
            
            sliderInput(
              "simulationCount",
              "Number of Simulations",
              min = 10000,
              max = 1000000,
              value = 500000
            ),
            numericInput(
              "aggregateSeed",
              "Random Seed",
              value = 123,
              min = 1,
              step = 1
            ),
            
            checkboxInput(
              "advanceAggregateSeed",
              "Advance seed after each run",
              value = TRUE
            ),
            
            helpText(
              "When enabled, the displayed seed is used for the current run, then increased by 1 for the next run. Disable this to repeat a simulation with the same seed and settings."
            ),
            
            textOutput("runtimeText"),
            textOutput("severityVarianceText")
          )
        ),
        
        # Diagnostics
        fluidRow(
          box(
            title = "Model Diagnostics",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            
            verbatimTextOutput("meanCheck"),
            verbatimTextOutput("sdCheck")
          )
        ),
        
        # Tables
        fluidRow(
          box(
            title = "Full Simulation Results",
            width = 12,
            DTOutput("lossTable")
          )
        ),
        
        fluidRow(
          box(
            title = "Simulations in Selected Histogram Bin",
            width = 12,
            
            downloadButton(
              "downloadSelectedBinCSV",
              "Download Selected Bin CSV"
            ),
            
            br(),
            br(),
            
            DTOutput("selectedBinTable")
          )
        ),
        
        creator_watermark()
      ),
      
      # ------------------------------------------------------------
      # 2.4 Galton board tab
      # ------------------------------------------------------------
      
      tabItem(
        tabName = "galton",
        
        fluidRow(
          box(
            title = "Galton Board Final Position Distribution",
            width = 9,
            status = "primary",
            solidHeader = TRUE,
            plotlyOutput("galtonHistogram", height = "600px")
          ),
          
          box(
            title = "Galton Board Controls",
            width = 3,
            
            sliderInput(
              "galtonBalls",
              "Number of Balls",
              min = 100,
              max = 100000,
              value = 10000,
              step = 100
            ),
            
            sliderInput(
              "galtonRows",
              "Number of Rows",
              min = 2,
              max = 50,
              value = 12,
              step = 2
            ),
            
            sliderInput(
              "galtonProbability",
              "Probability of Moving Right",
              min = 0,
              max = 1,
              value = 0.5,
              step = 0.01
            ),
            
            numericInput(
              "galtonSeed",
              "Random Seed",
              value = 123,
              min = 1,
              step = 1
            ),
            
            checkboxInput(
              "advanceGaltonSeed",
              "Advance seed after each run",
              value = TRUE
            ),
            
            helpText(
              "When enabled, the displayed seed is used for the current run, then increased by 1 for the next run. Disable this to repeat a simulation with the same seed and settings."
            ),
            
            selectInput(
              "galtonFillColour",
              "Bar Colour",
              choices = c(
                "Blue" = "steelblue",
                "Green" = "green",
                "Orange" = "orange",
                "Purple" = "purple",
                "Red" = "red"
              )
            ),
            
            selectInput(
              "galtonBorderColour",
              "Border Colour",
              choices = c(
                "Black" = "black",
                "Blue" = "steelblue",
                "Green" = "green",
                "Orange" = "orange",
                "Purple" = "purple",
                "Red" = "red"
              )
            ),
            
            checkboxInput(
              "showGaltonTheory",
              "Show Theoretical Binomial Distribution",
              value = TRUE
            ),
            
            selectInput(
              "galtonTheoryColour",
              "Theoretical Line Colour",
              choices = c(
                "Red" = "red",
                "Green" = "green",
                "Blue" = "steelblue",
                "Orange" = "orange",
                "Purple" = "purple",
                "Black" = "black"
              )
            ),
            
            actionButton(
              "runGalton",
              "Drop Balls",
              icon = icon("play")
            ),
            
            br(),
            br(),
            
            downloadButton(
              "downloadGaltonHistogram",
              "Download Galton Histogram"
            ),
            
            br(),
            br(),
            
            downloadButton(
              "downloadGaltonCSV",
              "Download Full Galton CSV"
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Galton Board Summary",
            width = 12,
            textOutput("galtonSummaryText"),
            textOutput("galtonRuntimeText"),
            textOutput("selectedGaltonBarSummary")
          )
        ),
        
        fluidRow(
          box(
            title = "Theoretical vs Simulated Diagnostics",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            verbatimTextOutput("galtonDiagnostics")
          )
        ),
        
        fluidRow(
          box(
            title = "Galton Board Simulation Results",
            width = 12,
            DTOutput("galtonTable")
          )
        ),
        fluidRow(
          box(
            title = "Balls in Selected Galton Bar",
            width = 12,
            
            downloadButton(
              "downloadSelectedGaltonBarCSV",
              "Download Selected Galton Bar CSV"
            ),
            
            br(),
            br(),
            
            DTOutput("selectedGaltonBarTable")
          )
        ),
        
        creator_watermark()
      ),
      
      # ------------------------------------------------------------
      # 2.5 Distribution explorer tab
      # ------------------------------------------------------------
      
      tabItem(
        tabName = "distributions",
        
        tabsetPanel(
          id = "distributionTabs",
          type = "tabs",
          
          tabPanel(
            title = "Normal",
            value = "normal",
            
            br(),
            
            fluidRow(
              box(
                title = "Normal Distribution Parameters",
                width = 3,
                status = "primary",
                solidHeader = TRUE,
                
                numericInput(
                  "normalMean",
                  "Mean (\u03bc)",
                  value = 0,
                  step = 0.5
                ),
                
                numericInput(
                  "normalSD",
                  "Standard deviation (\u03c3)",
                  value = 1,
                  min = 0.01,
                  step = 0.1
                ),
                
                sliderInput(
                  "normalSDRange",
                  "Standard deviations displayed",
                  min = 2,
                  max = 10,
                  value = 4,
                  step = 1
                ),
                
                selectInput(
                  "normalCurveColour",
                  "Curve Colour",
                  choices = c(
                    "Blue" = "steelblue",
                    "Red" = "red",
                    "Green" = "green",
                    "Orange" = "orange",
                    "Purple" = "purple",
                    "Black" = "black"
                  ),
                  selected = "steelblue"
                )
              ),
              
              box(
                title = "Normal Distribution",
                width = 9,
                status = "primary",
                solidHeader = TRUE,
                
                plotlyOutput(
                  "normalDistributionPlot",
                  height = "420px"
                )
              )
            ),
            
            fluidRow(
              box(
                title = "Distribution Information",
                width = 12,
                status = "info",
                solidHeader = TRUE,
                
                withMathJax(
                  fluidRow(
                    column(
                      width = 4,
                      h4("Definition"),
                      p(HTML("\\(X \\sim N(\\mu,\\sigma^2)\\)")),
                      p(HTML(
                        "\\(f(x)=\\frac{1}{\\sigma\\sqrt{2\\pi}}",
                        "e^{-\\frac{(x-\\mu)^2}{2\\sigma^2}}\\)"
                      ))
                    ),
                    
                    column(
                      width = 4,
                      h4("Properties"),
                      p(HTML("Support: \\(-\\infty < x < \\infty\\)")),
                      p(HTML("Mean: \\(E(X)=\\mu\\)")),
                      p(HTML("Variance: \\(\\operatorname{Var}(X)=\\sigma^2\\)"))
                    ),
                    
                    column(
                      width = 4,
                      h4("Interpretation"),
                      p(
                        "The Normal distribution models continuous values",
                        "that cluster symmetrically around a mean.",
                        "The mean controls its centre and the standard",
                        "deviation controls its spread."
                      )
                    )
                  )
                )
              )
            ),
            
            fluidRow(
              box(
                title = "Probability Calculator",
                width = 12,
                status = "warning",
                solidHeader = TRUE,
                
                fluidRow(
                  column(
                    width = 4,
                    
                    selectInput(
                      "normalProbabilityType",
                      "Calculate",
                      choices = c(
                        "P(X \u2264 x)" = "lower",
                        "P(X \u2265 x)" = "upper",
                        "P(a \u2264 X \u2264 b)" = "between"
                      ),
                      selected = "lower"
                    )
                  ),
                  
                  column(
                    width = 4,
                    
                    conditionalPanel(
                      condition = paste0(
                        "input.normalProbabilityType === 'lower' || ",
                        "input.normalProbabilityType === 'upper'"
                      ),
                      
                      numericInput(
                        "normalProbabilityX",
                        "Boundary (x)",
                        value = 0,
                        step = 0.1
                      )
                    ),
                    
                    conditionalPanel(
                      condition = "input.normalProbabilityType === 'between'",
                      
                      fluidRow(
                        column(
                          width = 6,
                          numericInput(
                            "normalProbabilityLower",
                            "Lower boundary (a)",
                            value = -1,
                            step = 0.1
                          )
                        ),
                        
                        column(
                          width = 6,
                          numericInput(
                            "normalProbabilityUpper",
                            "Upper boundary (b)",
                            value = 1,
                            step = 0.1
                          )
                        )
                      )
                    )
                  ),
                  
                  column(
                    width = 4,
                    h4("Result"),
                    textOutput("normalProbabilityResult")
                  )
                )
              )
            )
          ),
          
          tabPanel(
            title = "Binomial",
            value = "binomial",
            
            br(),
            
            fluidRow(
              box(
                title = "Binomial Distribution Parameters",
                width = 3,
                status = "primary",
                solidHeader = TRUE,
                
                numericInput(
                  "binomialTrials",
                  "Number of trials (n)",
                  value = 10,
                  min = 1,
                  max = 200,
                  step = 1
                ),
                
                sliderInput(
                  "binomialProbability",
                  "Success probability (p)",
                  min = 0,
                  max = 1,
                  value = 0.5,
                  step = 0.01
                ),
                
                selectInput(
                  "binomialBarColour",
                  "Bar Colour",
                  choices = c(
                    "Blue" = "steelblue",
                    "Red" = "red",
                    "Green" = "green",
                    "Orange" = "orange",
                    "Purple" = "purple",
                    "Black" = "black"
                  ),
                  selected = "steelblue"
                )
              ),
              
              box(
                title = "Binomial Distribution",
                width = 9,
                status = "primary",
                solidHeader = TRUE,
                
                plotlyOutput(
                  "binomialDistributionPlot",
                  height = "420px"
                )
              )
            ),
            
            fluidRow(
              box(
                title = "Distribution Information",
                width = 12,
                status = "info",
                solidHeader = TRUE,
                
                withMathJax(
                  fluidRow(
                    column(
                      width = 4,
                      h4("Definition"),
                      p(HTML("\\(X \\sim \\operatorname{Bin}(n,p)\\)")),
                      p(HTML(
                        "\\(P(X=x)=\\binom{n}{x}",
                        "p^x(1-p)^{n-x}\\)"
                      ))
                    ),
                    
                    column(
                      width = 4,
                      h4("Properties"),
                      p(HTML("Support: \\(x\\in\\{0,1,\\ldots,n\\}\\)")),
                      p(HTML("Mean: \\(E(X)=np\\)")),
                      p(HTML(
                        "Variance: ",
                        "\\(\\operatorname{Var}(X)=np(1-p)\\)"
                      ))
                    ),
                    
                    column(
                      width = 4,
                      h4("Interpretation"),
                      p(
                        "The Binomial distribution models the number of",
                        "successes in a fixed number of independent trials.",
                        "Each trial has the same probability of success."
                      )
                    )
                  )
                )
              )
            ),
            
            fluidRow(
              box(
                title = "Probability Calculator",
                width = 12,
                status = "warning",
                solidHeader = TRUE,
                
                fluidRow(
                  column(
                    width = 4,
                    
                    selectInput(
                      "binomialProbabilityType",
                      "Calculate",
                      choices = c(
                        "P(X \u2264 x)" = "lower",
                        "P(X \u2265 x)" = "upper",
                        "P(a \u2264 X \u2264 b)" = "between"
                      ),
                      selected = "lower"
                    )
                  ),
                  
                  column(
                    width = 4,
                    
                    conditionalPanel(
                      condition = paste0(
                        "input.binomialProbabilityType === 'lower' || ",
                        "input.binomialProbabilityType === 'upper'"
                      ),
                      
                      numericInput(
                        "binomialProbabilityX",
                        "Boundary (x)",
                        value = 5,
                        step = 1
                      )
                    ),
                    
                    conditionalPanel(
                      condition =
                        "input.binomialProbabilityType === 'between'",
                      
                      fluidRow(
                        column(
                          width = 6,
                          numericInput(
                            "binomialProbabilityLower",
                            "Lower boundary (a)",
                            value = 3,
                            step = 1
                          )
                        ),
                        
                        column(
                          width = 6,
                          numericInput(
                            "binomialProbabilityUpper",
                            "Upper boundary (b)",
                            value = 7,
                            step = 1
                          )
                        )
                      )
                    )
                  ),
                  
                  column(
                    width = 4,
                    h4("Result"),
                    textOutput("binomialProbabilityResult")
                  )
                )
              )
            )
          ),
          
          tabPanel(
            title = "Poisson",
            value = "poisson",
            
            br(),
            
            fluidRow(
              box(
                title = "Poisson Distribution Parameters",
                width = 3,
                status = "primary",
                solidHeader = TRUE,
                
                numericInput(
                  "poissonLambda",
                  "Expected event count (\u03bb)",
                  value = 5,
                  min = 0,
                  max = 100,
                  step = 0.1
                ),
                
                selectInput(
                  "poissonBarColour",
                  "Bar Colour",
                  choices = c(
                    "Blue" = "steelblue",
                    "Red" = "red",
                    "Green" = "green",
                    "Orange" = "orange",
                    "Purple" = "purple",
                    "Black" = "black"
                  ),
                  selected = "steelblue"
                )
              ),
              
              box(
                title = "Poisson Distribution",
                width = 9,
                status = "primary",
                solidHeader = TRUE,
                
                plotlyOutput(
                  "poissonDistributionPlot",
                  height = "420px"
                )
              )
            ),
            
            fluidRow(
              box(
                title = "Distribution Information",
                width = 12,
                status = "info",
                solidHeader = TRUE,
                
                withMathJax(
                  fluidRow(
                    column(
                      width = 4,
                      h4("Definition"),
                      p(HTML(
                        "\\(X \\sim \\operatorname{Poisson}(\\lambda)\\)"
                      )),
                      p(HTML(
                        "\\(P(X=x)=\\frac{e^{-\\lambda}",
                        "\\lambda^x}{x!}\\)"
                      ))
                    ),
                    
                    column(
                      width = 4,
                      h4("Properties"),
                      p(HTML(
                        "Support: \\(x\\in\\{0,1,2,\\ldots\\}\\)"
                      )),
                      p(HTML("Mean: \\(E(X)=\\lambda\\)")),
                      p(HTML(
                        "Variance: ",
                        "\\(\\operatorname{Var}(X)=\\lambda\\)"
                      ))
                    ),
                    
                    column(
                      width = 4,
                      h4("Interpretation"),
                      p(
                        "The Poisson distribution models the number of",
                        "events occurring within a fixed interval when",
                        "events occur independently at a constant average rate."
                      ),
                      p(
                        "The graph displays enough event counts to contain",
                        "at least 99.9999% of the total probability."
                      )
                    )
                  )
                )
              )
            ),
            
            fluidRow(
              box(
                title = "Probability Calculator",
                width = 12,
                status = "warning",
                solidHeader = TRUE,
                
                fluidRow(
                  column(
                    width = 4,
                    
                    selectInput(
                      "poissonProbabilityType",
                      "Calculate",
                      choices = c(
                        "P(X \u2264 x)" = "lower",
                        "P(X \u2265 x)" = "upper",
                        "P(a \u2264 X \u2264 b)" = "between"
                      ),
                      selected = "lower"
                    )
                  ),
                  
                  column(
                    width = 4,
                    
                    conditionalPanel(
                      condition = paste0(
                        "input.poissonProbabilityType === 'lower' || ",
                        "input.poissonProbabilityType === 'upper'"
                      ),
                      
                      numericInput(
                        "poissonProbabilityX",
                        "Boundary (x)",
                        value = 5,
                        step = 1
                      )
                    ),
                    
                    conditionalPanel(
                      condition =
                        "input.poissonProbabilityType === 'between'",
                      
                      fluidRow(
                        column(
                          width = 6,
                          numericInput(
                            "poissonProbabilityLower",
                            "Lower boundary (a)",
                            value = 2,
                            step = 1
                          )
                        ),
                        
                        column(
                          width = 6,
                          numericInput(
                            "poissonProbabilityUpper",
                            "Upper boundary (b)",
                            value = 8,
                            step = 1
                          )
                        )
                      )
                    )
                  ),
                  
                  column(
                    width = 4,
                    h4("Result"),
                    textOutput("poissonProbabilityResult")
                  )
                )
              )
            )
          ),
          
          tabPanel(
            title = "Exponential",
            value = "exponential",
            
            br(),
            
            fluidRow(
              box(
                title = "Exponential Distribution Parameters",
                width = 3,
                status = "primary",
                solidHeader = TRUE,
                
                numericInput(
                  "exponentialRate",
                  "Rate (\u03bb)",
                  value = 1,
                  min = 0.01,
                  step = 0.1
                ),
                
                selectInput(
                  "exponentialCurveColour",
                  "Curve Colour",
                  choices = c(
                    "Blue" = "steelblue",
                    "Red" = "red",
                    "Green" = "green",
                    "Orange" = "orange",
                    "Purple" = "purple",
                    "Black" = "black"
                  ),
                  selected = "steelblue"
                )
              ),
              
              box(
                title = "Exponential Distribution",
                width = 9,
                status = "primary",
                solidHeader = TRUE,
                
                plotlyOutput(
                  "exponentialDistributionPlot",
                  height = "420px"
                )
              )
            ),
            
            fluidRow(
              box(
                title = "Distribution Information",
                width = 12,
                status = "info",
                solidHeader = TRUE,
                
                withMathJax(
                  fluidRow(
                    column(
                      width = 4,
                      h4("Definition"),
                      p(HTML(
                        "\\(X \\sim \\operatorname{Exp}(\\lambda)\\)"
                      )),
                      p(HTML(
                        "\\(f(x)=\\lambda e^{-\\lambda x},\\quad x\\geq0\\)"
                      ))
                    ),
                    
                    column(
                      width = 4,
                      h4("Properties"),
                      p(HTML("Support: \\(0\\leq x<\\infty\\)")),
                      p(HTML("Mean: \\(E(X)=\\frac{1}{\\lambda}\\)")),
                      p(HTML(
                        "Variance: ",
                        "\\(\\operatorname{Var}(X)=\\frac{1}{\\lambda^2}\\)"
                      ))
                    ),
                    
                    column(
                      width = 4,
                      h4("Interpretation"),
                      p(
                        "The Exponential distribution models the waiting",
                        "time until the next event in a process where events",
                        "occur independently at a constant average rate."
                      ),
                      p(
                        "A larger rate produces shorter typical waiting times.",
                        "The graph displays at least 99.9999% of the",
                        "distribution's probability."
                      )
                    )
                  )
                )
              )
            ),
            
            fluidRow(
              box(
                title = "Probability Calculator",
                width = 12,
                status = "warning",
                solidHeader = TRUE,
                
                fluidRow(
                  column(
                    width = 4,
                    
                    selectInput(
                      "exponentialProbabilityType",
                      "Calculate",
                      choices = c(
                        "P(X \u2264 x)" = "lower",
                        "P(X \u2265 x)" = "upper",
                        "P(a \u2264 X \u2264 b)" = "between"
                      ),
                      selected = "lower"
                    )
                  ),
                  
                  column(
                    width = 4,
                    
                    conditionalPanel(
                      condition = paste0(
                        "input.exponentialProbabilityType === 'lower' || ",
                        "input.exponentialProbabilityType === 'upper'"
                      ),
                      
                      numericInput(
                        "exponentialProbabilityX",
                        "Boundary (x)",
                        value = 1,
                        step = 0.1
                      )
                    ),
                    
                    conditionalPanel(
                      condition =
                        "input.exponentialProbabilityType === 'between'",
                      
                      fluidRow(
                        column(
                          width = 6,
                          numericInput(
                            "exponentialProbabilityLower",
                            "Lower boundary (a)",
                            value = 0.5,
                            step = 0.1
                          )
                        ),
                        
                        column(
                          width = 6,
                          numericInput(
                            "exponentialProbabilityUpper",
                            "Upper boundary (b)",
                            value = 2,
                            step = 0.1
                          )
                        )
                      )
                    )
                  ),
                  
                  column(
                    width = 4,
                    h4("Result"),
                    textOutput("exponentialProbabilityResult")
                  )
                )
              )
            )
          ),
          
          tabPanel(
            title = "Gamma",
            value = "gamma",
            
            br(),
            
            fluidRow(
              box(
                title = "Gamma Distribution Parameters",
                width = 3,
                status = "primary",
                solidHeader = TRUE,
                
                numericInput(
                  "gammaShape",
                  "Shape (\u03b1)",
                  value = 2,
                  min = 0.1,
                  step = 0.1
                ),
                
                numericInput(
                  "gammaRate",
                  "Rate (\u03b2)",
                  value = 1,
                  min = 0.01,
                  step = 0.1
                ),
                
                selectInput(
                  "gammaCurveColour",
                  "Curve Colour",
                  choices = c(
                    "Blue" = "steelblue",
                    "Red" = "red",
                    "Green" = "green",
                    "Orange" = "orange",
                    "Purple" = "purple",
                    "Black" = "black"
                  ),
                  selected = "steelblue"
                )
              ),
              
              box(
                title = "Gamma Distribution",
                width = 9,
                status = "primary",
                solidHeader = TRUE,
                
                plotlyOutput(
                  "gammaDistributionPlot",
                  height = "420px"
                )
              )
            ),
            
            fluidRow(
              box(
                title = "Distribution Information",
                width = 12,
                status = "info",
                solidHeader = TRUE,
                
                withMathJax(
                  fluidRow(
                    column(
                      width = 4,
                      h4("Definition"),
                      p(HTML(
                        "\\(X \\sim \\operatorname{Gamma}(\\alpha,\\beta)\\)"
                      )),
                      p(HTML(
                        "\\(f(x)=\\frac{\\beta^\\alpha}{\\Gamma(\\alpha)}",
                        "x^{\\alpha-1}e^{-\\beta x},\\quad x\\geq0\\)"
                      ))
                    ),
                    
                    column(
                      width = 4,
                      h4("Properties"),
                      p(HTML("Support: \\(0\\leq x<\\infty\\)")),
                      p(HTML(
                        "Mean: \\(E(X)=\\frac{\\alpha}{\\beta}\\)"
                      )),
                      p(HTML(
                        "Variance: ",
                        "\\(\\operatorname{Var}(X)=",
                        "\\frac{\\alpha}{\\beta^2}\\)"
                      )),
                      p(HTML("Scale: \\(\\theta=\\frac{1}{\\beta}\\)"))
                    ),
                    
                    column(
                      width = 4,
                      h4("Interpretation"),
                      p(
                        "The Gamma distribution models positive,",
                        "right-skewed quantities such as waiting times",
                        "and claim severities."
                      ),
                      p(
                        "For integer shape, it can represent the waiting",
                        "time until the \u03b1th event in a constant-rate",
                        "Poisson process."
                      ),
                      p(
                        "The graph's upper limit contains at least 99.9999%",
                        "of the distribution's probability."
                      ),
                      p(
                        "When \u03b1 is below 1, the theoretical density is",
                        "unbounded at zero, so the displayed curve begins",
                        "just above zero."
                      )
                    )
                  )
                )
              )
            ),
            
            fluidRow(
              box(
                title = "Probability Calculator",
                width = 12,
                status = "warning",
                solidHeader = TRUE,
                
                fluidRow(
                  column(
                    width = 4,
                    
                    selectInput(
                      "gammaProbabilityType",
                      "Calculate",
                      choices = c(
                        "P(X \u2264 x)" = "lower",
                        "P(X \u2265 x)" = "upper",
                        "P(a \u2264 X \u2264 b)" = "between"
                      ),
                      selected = "lower"
                    )
                  ),
                  
                  column(
                    width = 4,
                    
                    conditionalPanel(
                      condition = paste0(
                        "input.gammaProbabilityType === 'lower' || ",
                        "input.gammaProbabilityType === 'upper'"
                      ),
                      
                      numericInput(
                        "gammaProbabilityX",
                        "Boundary (x)",
                        value = 2,
                        step = 0.1
                      )
                    ),
                    
                    conditionalPanel(
                      condition = "input.gammaProbabilityType === 'between'",
                      
                      fluidRow(
                        column(
                          width = 6,
                          numericInput(
                            "gammaProbabilityLower",
                            "Lower boundary (a)",
                            value = 1,
                            step = 0.1
                          )
                        ),
                        
                        column(
                          width = 6,
                          numericInput(
                            "gammaProbabilityUpper",
                            "Upper boundary (b)",
                            value = 4,
                            step = 0.1
                          )
                        )
                      )
                    )
                  ),
                  
                  column(
                    width = 4,
                    h4("Result"),
                    textOutput("gammaProbabilityResult")
                  )
                )
              )
            )
          ),
          
          tabPanel(
            title = "Lognormal",
            value = "lognormal",
            
            br(),
            
            fluidRow(
              box(
                title = "Lognormal Distribution Parameters",
                width = 3,
                status = "primary",
                solidHeader = TRUE,
                
                numericInput(
                  "lognormalMeanlog",
                  "Mean of ln(X) (\u03bc)",
                  value = 0,
                  step = 0.1
                ),
                
                numericInput(
                  "lognormalSdlog",
                  "Standard deviation of ln(X) (\u03c3)",
                  value = 1,
                  min = 0.01,
                  step = 0.1
                ),
                
                selectInput(
                  "lognormalCurveColour",
                  "Curve Colour",
                  choices = c(
                    "Blue" = "steelblue",
                    "Red" = "red",
                    "Green" = "green",
                    "Orange" = "orange",
                    "Purple" = "purple",
                    "Black" = "black"
                  ),
                  selected = "steelblue"
                )
              ),
              
              box(
                title = "Lognormal Distribution",
                width = 9,
                status = "primary",
                solidHeader = TRUE,
                
                plotlyOutput(
                  "lognormalDistributionPlot",
                  height = "420px"
                )
              )
            ),
            
            fluidRow(
              box(
                title = "Distribution Information",
                width = 12,
                status = "info",
                solidHeader = TRUE,
                
                withMathJax(
                  fluidRow(
                    column(
                      width = 4,
                      h4("Definition"),
                      p(HTML(
                        "\\(\\ln(X)\\sim N(\\mu,\\sigma^2)\\)"
                      )),
                      p(HTML(
                        "\\(f(x)=\\frac{1}{x\\sigma\\sqrt{2\\pi}}",
                        "\\exp\\!\\left[-\\frac{(\\ln x-\\mu)^2}",
                        "{2\\sigma^2}\\right],\\quad x>0\\)"
                      ))
                    ),
                    
                    column(
                      width = 4,
                      h4("Properties"),
                      p(HTML("Support: \\(0<x<\\infty\\)")),
                      p(HTML(
                        "Mean: \\(E(X)=e^{\\mu+\\sigma^2/2}\\)"
                      )),
                      p(HTML(
                        "Variance: ",
                        "\\(\\operatorname{Var}(X)=",
                        "(e^{\\sigma^2}-1)e^{2\\mu+\\sigma^2}\\)"
                      )),
                      p(HTML("Median: \\(e^\\mu\\)")),
                      p(HTML("Mode: \\(e^{\\mu-\\sigma^2}\\)"))
                    ),
                    
                    column(
                      width = 4,
                      h4("Interpretation"),
                      p(
                        "The Lognormal distribution models positive,",
                        "right-skewed quantities whose logarithms are",
                        "Normally distributed."
                      ),
                      p(
                        "It is useful for quantities formed by",
                        "multiplicative effects, including some asset",
                        "prices, incomes, and claim severities."
                      ),
                      p(
                        "\u03bc controls the location on the log scale, while",
                        "\u03c3 controls spread and right-skewness."
                      ),
                      p(
                        "The graph's upper limit contains at least 99.9999%",
                        "of the distribution's probability."
                      )
                    )
                  )
                )
              )
            ),
            
            fluidRow(
              box(
                title = "Probability Calculator",
                width = 12,
                status = "warning",
                solidHeader = TRUE,
                
                fluidRow(
                  column(
                    width = 4,
                    
                    selectInput(
                      "lognormalProbabilityType",
                      "Calculate",
                      choices = c(
                        "P(X \u2264 x)" = "lower",
                        "P(X \u2265 x)" = "upper",
                        "P(a \u2264 X \u2264 b)" = "between"
                      ),
                      selected = "lower"
                    )
                  ),
                  
                  column(
                    width = 4,
                    
                    conditionalPanel(
                      condition = paste0(
                        "input.lognormalProbabilityType === 'lower' || ",
                        "input.lognormalProbabilityType === 'upper'"
                      ),
                      
                      numericInput(
                        "lognormalProbabilityX",
                        "Boundary (x)",
                        value = 1,
                        step = 0.1
                      )
                    ),
                    
                    conditionalPanel(
                      condition =
                        "input.lognormalProbabilityType === 'between'",
                      
                      fluidRow(
                        column(
                          width = 6,
                          numericInput(
                            "lognormalProbabilityLower",
                            "Lower boundary (a)",
                            value = 0.5,
                            step = 0.1
                          )
                        ),
                        
                        column(
                          width = 6,
                          numericInput(
                            "lognormalProbabilityUpper",
                            "Upper boundary (b)",
                            value = 2,
                            step = 0.1
                          )
                        )
                      )
                    )
                  ),
                  
                  column(
                    width = 4,
                    h4("Result"),
                    textOutput("lognormalProbabilityResult")
                  )
                )
              )
            )
          )
        ),
        
        creator_watermark()
      ),
      
      # ------------------------------------------------------------
      # 2.6 Pi estimator tab
      # ------------------------------------------------------------
      
      tabItem(
        tabName = "piEstimator",
        
        fluidRow(
          box(
            title = "Quarter-Circle Trapezium Approximation",
            width = 9,
            status = "primary",
            solidHeader = TRUE,
            plotlyOutput("piEstimatorPlot", height = "600px")
          ),
          
          box(
            title = "Estimator Controls",
            width = 3,
            status = "primary",
            solidHeader = TRUE,
            
            sliderInput(
              "piTrapeziums",
              "Number of Trapeziums",
              min = 1,
              max = 200,
              value = 10,
              step = 1,
              sep = ""
            ),
            
            helpText(
              "The curve is y = sqrt(1 - x^2) from x = 0 to x = 1. Increasing the number of trapeziums improves the approximation."
            )
          )
        ),
        
        fluidRow(
          valueBoxOutput(
            "quarterPiEstimateBox",
            width = 12
          )
        ),
        
        fluidRow(
          valueBoxOutput(
            "piEstimateBox",
            width = 12
          )
        ),
        
        creator_watermark()
      )
    )
  )
)

# ============================================================
# 3. SERVER
# ============================================================

server <- function(input, output, session) {
  
  # ------------------------------------------------------------
  # 3.1 Page navigation
  # ------------------------------------------------------------
  
  observeEvent(input$goToSimulator, {
    updateTabItems(
      session,
      "tabs",
      selected = "aggregateLoss"
    )
  })
  
  observeEvent(input$goToGalton, {
    updateTabItems(
      session,
      "tabs",
      selected = "galton"
    )
  })
  
  observeEvent(input$goToDistributions, {
    updateTabItems(
      session,
      "tabs",
      selected = "distributions"
    )
  })
  
  observeEvent(input$goToPiEstimator, {
    updateTabItems(
      session,
      "tabs",
      selected = "piEstimator"
    )
  })
  
  # ------------------------------------------------------------
  # 3.2 Main simulation reactive
  # ------------------------------------------------------------
  
  simulationResults <- eventReactive(
    input$runSimulation,
    {
      req(input$aggregateSeed)
      
      seed_used <- as.integer(input$aggregateSeed)
      
      set.seed(seed_used)
      
      results <- run_aggregate_loss_simulation(
        n_sims = input$simulationCount,
        claim_frequency = input$claimFrequency,
        severity_mean = input$severityMean,
        severity_cv = input$severityCV
      )
      
      results$seed_used <- seed_used
      
      if (isTRUE(input$advanceAggregateSeed)) {
        next_seed <- if (
          seed_used >= .Machine$integer.max
        ) {
          1L
        } else {
          seed_used + 1L
        }
        
        updateNumericInput(
          session,
          "aggregateSeed",
          value = next_seed
        )
      }
      
      results
    },
    ignoreInit = TRUE
  )
  
  # ------------------------------------------------------------
  # 3.3 Galton board simulation reactive
  # ------------------------------------------------------------
  
  galtonResults <- eventReactive(
    input$runGalton,
    {
      req(input$galtonSeed)
      
      seed_used <- as.integer(input$galtonSeed)
      
      set.seed(seed_used)
      
      results <- run_galton_simulation(
        balls = input$galtonBalls,
        rows = input$galtonRows,
        p = input$galtonProbability
      )
      
      results$seed_used <- seed_used
      
      if (isTRUE(input$advanceGaltonSeed)) {
        next_seed <- if (
          seed_used >= .Machine$integer.max
        ) {
          1L
        } else {
          seed_used + 1L
        }
        
        updateNumericInput(
          session,
          "galtonSeed",
          value = next_seed
        )
      }
      
      results
    },
    ignoreInit = TRUE
  )
  
  # ------------------------------------------------------------
  # 3.4 Galton board histogram
  # ------------------------------------------------------------
  
  output$galtonHistogram <- renderPlotly({
    req(galtonResults())
    
    galton_data <- galtonResults()$data
    
    slot_counts <- as.data.frame(
      table(galton_data$FinalPosition)
    )
    
    names(slot_counts) <- c(
      "FinalPosition",
      "Count"
    )
    
    slot_counts$FinalPosition <- as.numeric(
      as.character(slot_counts$FinalPosition)
    )
    
    slot_counts$SimulatedProbability <- slot_counts$Count /
      sum(slot_counts$Count)
    
    possible_right_moves <- 0:input$galtonRows
    
    possible_positions <- possible_right_moves -
      input$galtonRows / 2
    
    theory_data <- data.frame(
      FinalPosition = possible_positions,
      TheoreticalProbability = dbinom(
        possible_right_moves,
        size = input$galtonRows,
        prob = input$galtonProbability
      )
    )
    
    galton_plot <- plot_ly(
      source = "galtonHistogram"
    ) %>%
      add_bars(
        data = slot_counts,
        x = ~FinalPosition,
        y = ~SimulatedProbability,
        name = "Simulated",
        marker = list(
          color = input$galtonFillColour,
          line = list(
            color = input$galtonBorderColour,
            width = 1
          )
        )
      )
    
    if (input$showGaltonTheory) {
      galton_plot <- galton_plot %>%
        add_trace(
          data = theory_data,
          x = ~FinalPosition,
          y = ~TheoreticalProbability,
          type = "scatter",
          mode = "lines+markers",
          name = "Theoretical",
          line = list(
            color = input$galtonTheoryColour,
            width = 2
          ),
          marker = list(
            color = input$galtonTheoryColour,
            size = 6
          )
        )
    }
    
    galton_plot <- galton_plot %>%
      layout(
        title = "Galton Board Final Position Distribution",
        xaxis = list(
          title = "Final Position",
          dtick = 1
        ),
        yaxis = list(
          title = "Probability"
        ),
        bargap = 0.1
      )
    
    event_register(
      galton_plot,
      "plotly_click"
    )
  })
  
  # ------------------------------------------------------------
  # 3.5 Galton board summary outputs
  # ------------------------------------------------------------
  
  selectedGaltonBarData <- reactive({
    req(galtonResults())
    
    click <- event_data(
      "plotly_click",
      source = "galtonHistogram"
    )
    
    req(click)
    
    get_selected_galton_bar(
      galton_data = galtonResults()$data,
      clicked_x = click$x
    )
  })
  
  # ------------------------------------------------------------
  # 3.6 Galton board summary outputs
  # ------------------------------------------------------------
  
  output$galtonSummaryText <- renderText({
    req(galtonResults())
    
    paste(
      "Balls dropped:",
      format_whole_number(nrow(galtonResults()$data)),
      "| Rows:",
      input$galtonRows,
      "| Probability right:",
      input$galtonProbability,
      "| Seed used:",
      galtonResults()$seed_used
    )
  })
  
  output$galtonRuntimeText <- renderText({
    req(galtonResults())
    
    paste(
      "Runtime:",
      round(galtonResults()$runtime, 4),
      "seconds",
      "| Throughput:",
      format_runtime_rate(
        n_sims = input$galtonBalls,
        runtime_seconds = galtonResults()$runtime
      ),
      "balls/sec"
    )
  })
  
  output$selectedGaltonBarSummary <- renderText({
    selected <- selectedGaltonBarData()
    
    paste(
      "Selected final position:",
      selected$final_position,
      "| Balls in this position:",
      format_whole_number(nrow(selected$data))
    )
  })
  
  # ------------------------------------------------------------
  # 3.7 Histogram plot
  # ------------------------------------------------------------
  
  output$lossHistogram <- renderPlotly({
    req(simulationResults())
    
    sim_data <- simulationResults()$data
    
    p <- create_loss_histogram(
      sim_data = sim_data,
      bins = input$histogramBins,
      fill_colour = input$histogramFillColour,
      border_colour = input$histogramBorderColour,
      show_density = input$showDensity,
      density_colour = input$densityColour
    )
    
    gg <- ggplotly(
      p,
      source = "lossHistogram"
    )
    
    event_register(
      gg,
      "plotly_click"
    )
  })
  
  # ------------------------------------------------------------
  # 3.8 Clicked histogram bin
  # ------------------------------------------------------------
  
  selectedBinData <- reactive({
    req(simulationResults())
    
    click <- event_data(
      "plotly_click",
      source = "lossHistogram"
    )
    
    req(click)
    
    get_selected_bin(
      sim_data = simulationResults()$data,
      bins = input$histogramBins,
      clicked_x = click$x
    )
  })
  
  # ------------------------------------------------------------
  # 3.9 Summary outputs
  # ------------------------------------------------------------
  
  output$simulationCountText <- renderText({
    req(simulationResults())
    
    paste(
      "Simulations completed:",
      format_whole_number(nrow(simulationResults()$data)),
      "| Seed used:",
      simulationResults()$seed_used
    )
  })
  
  output$lossSummaryText <- renderText({
    req(simulationResults())
    
    sim_data <- simulationResults()$data
    
    paste(
      "Mean Loss:",
      format_whole_number(mean(sim_data$Loss)),
      "| SD:",
      format_whole_number(sd(sim_data$Loss)),
      "| 95th Percentile:",
      format_whole_number(quantile(sim_data$Loss, 0.95)),
      "| 99th Percentile:",
      format_whole_number(quantile(sim_data$Loss, 0.99)),
      "| 99.5th Percentile:",
      format_whole_number(quantile(sim_data$Loss, 0.995))
    )
  })
  
  output$runtimeText <- renderText({
    req(simulationResults())
    
    paste(
      "Runtime:",
      round(simulationResults()$runtime, 4),
      "seconds",
      "| Throughput:",
      format_runtime_rate(
        n_sims = input$simulationCount,
        runtime_seconds = simulationResults()$runtime
      ),
      "sims/sec"
    )
  })
  
  output$severityVarianceText <- renderText({
    req(simulationResults())
    
    paste(
      "Implied Severity Variance:",
      format(
        round(simulationResults()$severity_variance),
        scientific = TRUE
      )
    )
  })
  
  output$selectedBinSummary <- renderText({
    selected <- selectedBinData()
    
    paste(
      "Selected bin:",
      format_whole_number(selected$lower_bound),
      "to",
      format_whole_number(selected$upper_bound),
      "| Simulations in bin:",
      format_whole_number(nrow(selected$data))
    )
  })
  
  
  # ------------------------------------------------------------
  # 3.10 Tables
  # ------------------------------------------------------------
  
  output$lossTable <- renderDT({
    req(simulationResults())
    
    datatable(
      simulationResults()$data,
      options = list(pageLength = 10)
    )
  })
  
  output$selectedBinTable <- renderDT({
    selected <- selectedBinData()
    
    datatable(
      selected$data,
      options = list(pageLength = 10)
    )
  })
  
  output$selectedGaltonBarTable <- renderDT({
    selected <- selectedGaltonBarData()
    
    datatable(
      selected$data,
      options = list(pageLength = 10)
    )
  })
  
  # ------------------------------------------------------------
  # 3.10 Galton board table
  # ------------------------------------------------------------
  
  output$galtonTable <- renderDT({
    req(galtonResults())
    
    datatable(
      galtonResults()$data,
      options = list(pageLength = 10)
    )
  })
  
  # ------------------------------------------------------------
  # 3.11 Model diagnostics
  # ------------------------------------------------------------
  
  output$meanCheck <- renderText({
    req(simulationResults())
    
    results <- simulationResults()
    sim_data <- results$data
    
    expected_mean <- results$expected_mean
    simulated_mean <- results$simulated_mean
    
    n_sims <- nrow(sim_data)
    sample_sd <- sd(sim_data$Loss)
    standard_error <- sample_sd / sqrt(n_sims)
    
    z_score <- (
      simulated_mean - expected_mean
    ) / standard_error
    
    lower_ci <- simulated_mean - 1.96 * standard_error
    upper_ci <- simulated_mean + 1.96 * standard_error
    
    pass_fail <- if (
      expected_mean >= lower_ci &
      expected_mean <= upper_ci
    ) {
      "PASS"
    } else {
      "WARNING"
    }
    
    paste(
      "Compound Mean Test:", pass_fail,
      "\nExpected Mean:", format_whole_number(expected_mean),
      "\nSimulated Mean:", format_whole_number(simulated_mean),
      "\nDifference:", format_whole_number(simulated_mean - expected_mean),
      "\nStandard Error:", format_whole_number(standard_error),
      "\nZ-Score:", round(z_score, 3),
      "\n",
      "\n95% Confidence Interval:",
      paste(
        format_whole_number(lower_ci),
        "to",
        format_whole_number(upper_ci)
      )
    )
  })
  
  output$sdCheck <- renderText({
    req(simulationResults())
    
    results <- simulationResults()
    
    expected_sd <- results$expected_sd
    simulated_sd <- results$simulated_sd
    sd_difference <- results$sd_difference
    sd_difference_pct <- results$sd_difference_pct
    
    pass_fail <- if (
      abs(sd_difference_pct) < 2
    ) {
      "PASS"
    } else {
      "WARNING"
    }
    
    paste(
      "Compound SD Test:", pass_fail,
      "\nExpected SD:", format_whole_number(expected_sd),
      "\nSimulated SD:", format_whole_number(simulated_sd),
      "\nDifference:", format_whole_number(sd_difference),
      "\nPercentage Difference:", round(sd_difference_pct, 3), "%"
    )
  })
  
  # ------------------------------------------------------------
  # 3.12 Galton board diagnostics
  # ------------------------------------------------------------
  
  output$galtonDiagnostics <- renderText({
    req(galtonResults())
    
    results <- galtonResults()
    
    paste(
      "Galton Board Mean Check",
      "\nExpected Mean:", round(results$expected_mean, 4),
      "\nSimulated Mean:", round(results$simulated_mean, 4),
      "\nDifference:", round(results$mean_difference, 4),
      "\n",
      "\nGalton Board SD Check",
      "\nExpected SD:", round(results$expected_sd, 4),
      "\nSimulated SD:", round(results$simulated_sd, 4),
      "\nDifference:", round(results$sd_difference, 4)
    )
  })
  
  # ------------------------------------------------------------
  # 3.13 Normal distribution explorer
  # ------------------------------------------------------------
  
  normalDistributionData <- reactive({
    req(
      input$normalMean,
      input$normalSD,
      input$normalSDRange
    )
    
    validate(
      need(
        is.finite(input$normalMean),
        "The mean must be a finite number."
      ),
      need(
        is.finite(input$normalSD) && input$normalSD > 0,
        "The standard deviation must be greater than zero."
      ),
      need(
        is.finite(input$normalSDRange) &&
          input$normalSDRange >= 2 &&
          input$normalSDRange <= 10,
        "The display range must be between 2 and 10 standard deviations."
      )
    )
    
    generate_normal_data(
      mean = input$normalMean,
      sd = input$normalSD,
      display_sd = input$normalSDRange
    )
  })
  
  normalProbability <- reactive({
    req(input$normalProbabilityType)
    normalDistributionData()
    
    if (input$normalProbabilityType == "between") {
      req(
        input$normalProbabilityLower,
        input$normalProbabilityUpper
      )
      
      validate(
        need(
          input$normalProbabilityLower <= input$normalProbabilityUpper,
          "The lower boundary must not exceed the upper boundary."
        )
      )
      
      return(
        calculate_normal_probability(
          calculation_type = "between",
          mean = input$normalMean,
          sd = input$normalSD,
          lower = input$normalProbabilityLower,
          upper = input$normalProbabilityUpper
        )
      )
    }
    
    req(input$normalProbabilityX)
    
    calculate_normal_probability(
      calculation_type = input$normalProbabilityType,
      mean = input$normalMean,
      sd = input$normalSD,
      x = input$normalProbabilityX
    )
  })
  
  output$normalDistributionPlot <- renderPlotly({
    normal_data <- normalDistributionData()
    probability_type <- input$normalProbabilityType
    
    if (probability_type == "lower") {
      shaded_data <- normal_data[
        normal_data$x <= input$normalProbabilityX,
      ]
    } else if (probability_type == "upper") {
      shaded_data <- normal_data[
        normal_data$x >= input$normalProbabilityX,
      ]
    } else {
      validate(
        need(
          input$normalProbabilityLower <= input$normalProbabilityUpper,
          "The lower boundary must not exceed the upper boundary."
        )
      )
      
      shaded_data <- normal_data[
        normal_data$x >= input$normalProbabilityLower &
          normal_data$x <= input$normalProbabilityUpper,
      ]
    }
    
    normal_plot <- plot_ly() %>%
      add_trace(
        data = normal_data,
        x = ~x,
        y = ~density,
        type = "scatter",
        mode = "lines",
        name = "Density",
        line = list(
          color = input$normalCurveColour,
          width = 3
        ),
        hovertemplate = paste(
          "x: %{x:.3f}",
          "<br>Density: %{y:.5f}",
          "<extra></extra>"
        )
      )
    
    if (nrow(shaded_data) > 0) {
      normal_plot <- normal_plot %>%
        add_trace(
          x = c(
            shaded_data$x[1],
            shaded_data$x,
            shaded_data$x[nrow(shaded_data)]
          ),
          y = c(
            0,
            shaded_data$density,
            0
          ),
          type = "scatter",
          mode = "lines",
          fill = "toself",
          name = "Selected probability",
          line = list(width = 0),
          fillcolor = "rgba(60, 141, 188, 0.30)",
          hoverinfo = "skip"
        )
    }
    
    normal_plot %>%
      layout(
        title = paste0(
          "Normal Distribution: N(",
          input$normalMean,
          ", ",
          input$normalSD,
          "\u00b2)"
        ),
        xaxis = list(
          title = "x",
          range = c(
            input$normalMean -
              input$normalSDRange * input$normalSD,
            input$normalMean +
              input$normalSDRange * input$normalSD
          )
        ),
        yaxis = list(
          title = "Probability density",
          rangemode = "tozero"
        ),
        showlegend = TRUE,
        hovermode = "x unified"
      )
  })
  
  output$normalProbabilityResult <- renderText({
    paste0(
      "Probability = ",
      formatC(
        normalProbability(),
        format = "f",
        digits = 6
      )
    )
  })
  
  # ------------------------------------------------------------
  # 3.14 Binomial distribution explorer
  # ------------------------------------------------------------
  
  binomialDistributionData <- reactive({
    req(
      input$binomialTrials,
      input$binomialProbability
    )
    
    validate(
      need(
        is.finite(input$binomialTrials) &&
          input$binomialTrials >= 1 &&
          input$binomialTrials <= 200 &&
          input$binomialTrials == round(input$binomialTrials),
        "The number of trials must be a whole number between 1 and 200."
      ),
      need(
        is.finite(input$binomialProbability) &&
          input$binomialProbability >= 0 &&
          input$binomialProbability <= 1,
        "The success probability must be between 0 and 1."
      )
    )
    
    generate_binomial_data(
      trials = input$binomialTrials,
      probability = input$binomialProbability
    )
  })
  
  binomialProbabilityResult <- reactive({
    req(input$binomialProbabilityType)
    binomialDistributionData()
    
    if (input$binomialProbabilityType == "between") {
      req(
        input$binomialProbabilityLower,
        input$binomialProbabilityUpper
      )
      
      validate(
        need(
          input$binomialProbabilityLower <=
            input$binomialProbabilityUpper,
          "The lower boundary must not exceed the upper boundary."
        )
      )
      
      return(
        calculate_binomial_probability(
          calculation_type = "between",
          trials = input$binomialTrials,
          probability = input$binomialProbability,
          lower = input$binomialProbabilityLower,
          upper = input$binomialProbabilityUpper
        )
      )
    }
    
    req(input$binomialProbabilityX)
    
    calculate_binomial_probability(
      calculation_type = input$binomialProbabilityType,
      trials = input$binomialTrials,
      probability = input$binomialProbability,
      x = input$binomialProbabilityX
    )
  })
  
  output$binomialDistributionPlot <- renderPlotly({
    binomial_data <- binomialDistributionData()
    probability_type <- input$binomialProbabilityType
    
    if (probability_type == "lower") {
      selected_bars <- binomial_data$x <=
        floor(input$binomialProbabilityX)
    } else if (probability_type == "upper") {
      selected_bars <- binomial_data$x >=
        ceiling(input$binomialProbabilityX)
    } else {
      validate(
        need(
          input$binomialProbabilityLower <=
            input$binomialProbabilityUpper,
          "The lower boundary must not exceed the upper boundary."
        )
      )
      
      selected_bars <-
        binomial_data$x >= ceiling(input$binomialProbabilityLower) &
        binomial_data$x <= floor(input$binomialProbabilityUpper)
    }
    
    unselected_data <- binomial_data[!selected_bars, , drop = FALSE]
    selected_data <- binomial_data[selected_bars, , drop = FALSE]
    
    binomial_plot <- plot_ly()
    
    if (nrow(unselected_data) > 0) {
      binomial_plot <- binomial_plot %>%
        add_bars(
          data = unselected_data,
          x = ~x,
          y = ~probability,
          name = "Outside selected probability",
          marker = list(
            color = input$binomialBarColour,
            line = list(
              color = input$binomialBarColour,
              width = 1
            )
          ),
          hovertemplate = paste(
            "Successes: %{x}",
            "<br>Probability: %{y:.6f}",
            "<extra></extra>"
          )
        )
    }
    
    if (nrow(selected_data) > 0) {
      binomial_plot <- binomial_plot %>%
        add_bars(
          data = selected_data,
          x = ~x,
          y = ~probability,
          name = "Included in probability",
          marker = list(
            color = "rgb(243, 156, 18)",
            line = list(
              color = "rgb(243, 156, 18)",
              width = 1
            )
          ),
          hovertemplate = paste(
            "Successes: %{x}",
            "<br>Probability: %{y:.6f}",
            "<extra></extra>"
          )
        )
    }
    
    binomial_plot %>%
      layout(
        title = paste0(
          "Binomial Distribution: Bin(",
          input$binomialTrials,
          ", ",
          input$binomialProbability,
          ")"
        ),
        xaxis = list(
          title = "Number of successes",
          tickmode = "linear",
          dtick = 1,
          range = c(-0.5, input$binomialTrials + 0.5)
        ),
        yaxis = list(
          title = "Probability",
          rangemode = "tozero"
        ),
        barmode = "overlay",
        bargap = 0.15,
        showlegend = FALSE,
        hovermode = "closest"
      )
  })
  
  output$binomialProbabilityResult <- renderText({
    paste0(
      "Probability = ",
      formatC(
        binomialProbabilityResult(),
        format = "f",
        digits = 6
      )
    )
  })
  
  # ------------------------------------------------------------
  # 3.15 Poisson distribution explorer
  # ------------------------------------------------------------
  
  poissonDistributionData <- reactive({
    req(input$poissonLambda)
    
    validate(
      need(
        is.finite(input$poissonLambda) &&
          input$poissonLambda >= 0 &&
          input$poissonLambda <= 100,
        "Lambda must be between 0 and 100."
      )
    )
    
    generate_poisson_data(
      lambda = input$poissonLambda
    )
  })
  
  poissonProbabilityValue <- reactive({
    req(input$poissonProbabilityType)
    poissonDistributionData()
    
    if (input$poissonProbabilityType == "between") {
      req(
        input$poissonProbabilityLower,
        input$poissonProbabilityUpper
      )
      
      validate(
        need(
          input$poissonProbabilityLower <=
            input$poissonProbabilityUpper,
          "The lower boundary must not exceed the upper boundary."
        )
      )
      
      return(
        calculate_poisson_probability(
          calculation_type = "between",
          lambda = input$poissonLambda,
          lower = input$poissonProbabilityLower,
          upper = input$poissonProbabilityUpper
        )
      )
    }
    
    req(input$poissonProbabilityX)
    
    calculate_poisson_probability(
      calculation_type = input$poissonProbabilityType,
      lambda = input$poissonLambda,
      x = input$poissonProbabilityX
    )
  })
  
  output$poissonDistributionPlot <- renderPlotly({
    poisson_data <- poissonDistributionData()
    probability_type <- input$poissonProbabilityType
    
    if (probability_type == "lower") {
      selected_bars <- poisson_data$x <=
        floor(input$poissonProbabilityX)
    } else if (probability_type == "upper") {
      selected_bars <- poisson_data$x >=
        ceiling(input$poissonProbabilityX)
    } else {
      validate(
        need(
          input$poissonProbabilityLower <=
            input$poissonProbabilityUpper,
          "The lower boundary must not exceed the upper boundary."
        )
      )
      
      selected_bars <-
        poisson_data$x >= max(
          0,
          ceiling(input$poissonProbabilityLower)
        ) &
        poisson_data$x <= floor(input$poissonProbabilityUpper)
    }
    
    unselected_data <- poisson_data[
      !selected_bars,
      ,
      drop = FALSE
    ]
    
    selected_data <- poisson_data[
      selected_bars,
      ,
      drop = FALSE
    ]
    
    poisson_plot <- plot_ly()
    
    if (nrow(unselected_data) > 0) {
      poisson_plot <- poisson_plot %>%
        add_bars(
          data = unselected_data,
          x = ~x,
          y = ~probability,
          name = "Outside selected probability",
          marker = list(
            color = input$poissonBarColour,
            line = list(
              color = input$poissonBarColour,
              width = 1
            )
          ),
          hovertemplate = paste(
            "Events: %{x}",
            "<br>Probability: %{y:.6f}",
            "<extra></extra>"
          )
        )
    }
    
    if (nrow(selected_data) > 0) {
      poisson_plot <- poisson_plot %>%
        add_bars(
          data = selected_data,
          x = ~x,
          y = ~probability,
          name = "Included in probability",
          marker = list(
            color = "rgb(243, 156, 18)",
            line = list(
              color = "rgb(243, 156, 18)",
              width = 1
            )
          ),
          hovertemplate = paste(
            "Events: %{x}",
            "<br>Probability: %{y:.6f}",
            "<extra></extra>"
          )
        )
    }
    
    poisson_plot %>%
      layout(
        title = paste0(
          "Poisson Distribution: Pois(",
          input$poissonLambda,
          ")"
        ),
        xaxis = list(
          title = "Number of events",
          tickmode = "linear",
          dtick = 1,
          range = c(
            -0.5,
            max(poisson_data$x) + 0.5
          )
        ),
        yaxis = list(
          title = "Probability",
          rangemode = "tozero"
        ),
        barmode = "overlay",
        bargap = 0.15,
        showlegend = FALSE,
        hovermode = "closest"
      )
  })
  
  output$poissonProbabilityResult <- renderText({
    paste0(
      "Probability = ",
      formatC(
        poissonProbabilityValue(),
        format = "f",
        digits = 6
      )
    )
  })
  
  # ------------------------------------------------------------
  # 3.16 Exponential distribution explorer
  # ------------------------------------------------------------
  
  exponentialDistributionData <- reactive({
    req(input$exponentialRate)
    
    validate(
      need(
        is.finite(input$exponentialRate) &&
          input$exponentialRate > 0,
        "The rate must be greater than zero."
      )
    )
    
    generate_exponential_data(
      rate = input$exponentialRate
    )
  })
  
  exponentialProbabilityValue <- reactive({
    req(input$exponentialProbabilityType)
    exponentialDistributionData()
    
    if (input$exponentialProbabilityType == "between") {
      req(
        input$exponentialProbabilityLower,
        input$exponentialProbabilityUpper
      )
      
      validate(
        need(
          is.finite(input$exponentialProbabilityLower) &&
            is.finite(input$exponentialProbabilityUpper),
          "Both boundaries must be finite numbers."
        ),
        need(
          input$exponentialProbabilityLower <=
            input$exponentialProbabilityUpper,
          "The lower boundary must not exceed the upper boundary."
        )
      )
      
      return(
        calculate_exponential_probability(
          calculation_type = "between",
          rate = input$exponentialRate,
          lower = input$exponentialProbabilityLower,
          upper = input$exponentialProbabilityUpper
        )
      )
    }
    
    req(input$exponentialProbabilityX)
    
    validate(
      need(
        is.finite(input$exponentialProbabilityX),
        "The boundary must be a finite number."
      )
    )
    
    calculate_exponential_probability(
      calculation_type = input$exponentialProbabilityType,
      rate = input$exponentialRate,
      x = input$exponentialProbabilityX
    )
  })
  
  output$exponentialDistributionPlot <- renderPlotly({
    exponential_data <- exponentialDistributionData()
    probability_type <- input$exponentialProbabilityType
    
    if (probability_type == "lower") {
      req(input$exponentialProbabilityX)
      
      shaded_data <- exponential_data[
        exponential_data$x <= input$exponentialProbabilityX,
        ,
        drop = FALSE
      ]
    } else if (probability_type == "upper") {
      req(input$exponentialProbabilityX)
      
      shaded_data <- exponential_data[
        exponential_data$x >= input$exponentialProbabilityX,
        ,
        drop = FALSE
      ]
    } else {
      req(
        input$exponentialProbabilityLower,
        input$exponentialProbabilityUpper
      )
      
      validate(
        need(
          input$exponentialProbabilityLower <=
            input$exponentialProbabilityUpper,
          "The lower boundary must not exceed the upper boundary."
        )
      )
      
      shaded_data <- exponential_data[
        exponential_data$x >= input$exponentialProbabilityLower &
          exponential_data$x <= input$exponentialProbabilityUpper,
        ,
        drop = FALSE
      ]
    }
    
    exponential_plot <- plot_ly() %>%
      add_trace(
        data = exponential_data,
        x = ~x,
        y = ~density,
        type = "scatter",
        mode = "lines",
        name = "Density",
        line = list(
          color = input$exponentialCurveColour,
          width = 3
        ),
        hovertemplate = paste(
          "x: %{x:.3f}",
          "<br>Density: %{y:.5f}",
          "<extra></extra>"
        )
      )
    
    if (nrow(shaded_data) > 0) {
      exponential_plot <- exponential_plot %>%
        add_trace(
          x = c(
            shaded_data$x[1],
            shaded_data$x,
            shaded_data$x[nrow(shaded_data)]
          ),
          y = c(
            0,
            shaded_data$density,
            0
          ),
          type = "scatter",
          mode = "lines",
          fill = "toself",
          name = "Selected probability",
          line = list(width = 0),
          fillcolor = "rgba(243, 156, 18, 0.35)",
          hoverinfo = "skip"
        )
    }
    
    exponential_plot %>%
      layout(
        title = paste0(
          "Exponential Distribution: Exp(",
          input$exponentialRate,
          ")"
        ),
        xaxis = list(
          title = "x",
          range = c(
            0,
            max(exponential_data$x)
          )
        ),
        yaxis = list(
          title = "Probability density",
          rangemode = "tozero"
        ),
        showlegend = TRUE,
        hovermode = "x unified"
      )
  })
  
  output$exponentialProbabilityResult <- renderText({
    paste0(
      "Probability = ",
      formatC(
        exponentialProbabilityValue(),
        format = "f",
        digits = 6
      )
    )
  })
  
  # ------------------------------------------------------------
  # 3.17 Gamma distribution explorer
  # ------------------------------------------------------------
  
  gammaDistributionData <- reactive({
    req(
      input$gammaShape,
      input$gammaRate
    )
    
    validate(
      need(
        is.finite(input$gammaShape) &&
          input$gammaShape > 0,
        "The shape must be greater than zero."
      ),
      need(
        is.finite(input$gammaRate) &&
          input$gammaRate > 0,
        "The rate must be greater than zero."
      )
    )
    
    generate_gamma_data(
      shape = input$gammaShape,
      rate = input$gammaRate
    )
  })
  
  gammaProbabilityValue <- reactive({
    req(input$gammaProbabilityType)
    gammaDistributionData()
    
    if (input$gammaProbabilityType == "between") {
      req(
        input$gammaProbabilityLower,
        input$gammaProbabilityUpper
      )
      
      validate(
        need(
          is.finite(input$gammaProbabilityLower) &&
            is.finite(input$gammaProbabilityUpper),
          "Both boundaries must be finite numbers."
        ),
        need(
          input$gammaProbabilityLower <=
            input$gammaProbabilityUpper,
          "The lower boundary must not exceed the upper boundary."
        )
      )
      
      return(
        calculate_gamma_probability(
          calculation_type = "between",
          shape = input$gammaShape,
          rate = input$gammaRate,
          lower = input$gammaProbabilityLower,
          upper = input$gammaProbabilityUpper
        )
      )
    }
    
    req(input$gammaProbabilityX)
    
    validate(
      need(
        is.finite(input$gammaProbabilityX),
        "The boundary must be a finite number."
      )
    )
    
    calculate_gamma_probability(
      calculation_type = input$gammaProbabilityType,
      shape = input$gammaShape,
      rate = input$gammaRate,
      x = input$gammaProbabilityX
    )
  })
  
  output$gammaDistributionPlot <- renderPlotly({
    gamma_data <- gammaDistributionData()
    probability_type <- input$gammaProbabilityType
    
    if (probability_type == "lower") {
      req(input$gammaProbabilityX)
      
      shaded_data <- gamma_data[
        gamma_data$x <= input$gammaProbabilityX,
        ,
        drop = FALSE
      ]
    } else if (probability_type == "upper") {
      req(input$gammaProbabilityX)
      
      shaded_data <- gamma_data[
        gamma_data$x >= input$gammaProbabilityX,
        ,
        drop = FALSE
      ]
    } else {
      req(
        input$gammaProbabilityLower,
        input$gammaProbabilityUpper
      )
      
      validate(
        need(
          input$gammaProbabilityLower <=
            input$gammaProbabilityUpper,
          "The lower boundary must not exceed the upper boundary."
        )
      )
      
      shaded_data <- gamma_data[
        gamma_data$x >= input$gammaProbabilityLower &
          gamma_data$x <= input$gammaProbabilityUpper,
        ,
        drop = FALSE
      ]
    }
    
    gamma_plot <- plot_ly() %>%
      add_trace(
        data = gamma_data,
        x = ~x,
        y = ~density,
        type = "scatter",
        mode = "lines",
        name = "Density",
        line = list(
          color = input$gammaCurveColour,
          width = 3
        ),
        hovertemplate = paste(
          "x: %{x:.3f}",
          "<br>Density: %{y:.5f}",
          "<extra></extra>"
        )
      )
    
    if (nrow(shaded_data) > 0) {
      gamma_plot <- gamma_plot %>%
        add_trace(
          x = c(
            shaded_data$x[1],
            shaded_data$x,
            shaded_data$x[nrow(shaded_data)]
          ),
          y = c(
            0,
            shaded_data$density,
            0
          ),
          type = "scatter",
          mode = "lines",
          fill = "toself",
          name = "Selected probability",
          line = list(width = 0),
          fillcolor = "rgba(243, 156, 18, 0.35)",
          hoverinfo = "skip"
        )
    }
    
    gamma_plot %>%
      layout(
        title = paste0(
          "Gamma Distribution: Gamma(\u03b1 = ",
          input$gammaShape,
          ", \u03b2 = ",
          input$gammaRate,
          ")"
        ),
        xaxis = list(
          title = "x",
          range = c(
            0,
            max(gamma_data$x)
          )
        ),
        yaxis = list(
          title = "Probability density",
          rangemode = "tozero"
        ),
        showlegend = TRUE,
        hovermode = "x unified"
      )
  })
  
  output$gammaProbabilityResult <- renderText({
    paste0(
      "Probability = ",
      formatC(
        gammaProbabilityValue(),
        format = "f",
        digits = 6
      )
    )
  })
  
  # ------------------------------------------------------------
  # 3.18 Lognormal distribution explorer
  # ------------------------------------------------------------
  
  lognormalDistributionData <- reactive({
    req(
      input$lognormalMeanlog,
      input$lognormalSdlog
    )
    
    validate(
      need(
        is.finite(input$lognormalMeanlog),
        "The log-mean must be a finite number."
      ),
      need(
        is.finite(input$lognormalSdlog) &&
          input$lognormalSdlog > 0,
        "The log-standard deviation must be greater than zero."
      ),
      need(
        input$lognormalMeanlog +
          qnorm(0.999999) * input$lognormalSdlog <=
          log(.Machine$double.xmax),
        paste(
          "These parameters produce values too large to display.",
          "Reduce the log-mean or log-standard deviation."
        )
      )
    )
    
    generate_lognormal_data(
      meanlog = input$lognormalMeanlog,
      sdlog = input$lognormalSdlog
    )
  })
  
  lognormalProbabilityValue <- reactive({
    req(input$lognormalProbabilityType)
    lognormalDistributionData()
    
    if (input$lognormalProbabilityType == "between") {
      req(
        input$lognormalProbabilityLower,
        input$lognormalProbabilityUpper
      )
      
      validate(
        need(
          is.finite(input$lognormalProbabilityLower) &&
            is.finite(input$lognormalProbabilityUpper),
          "Both boundaries must be finite numbers."
        ),
        need(
          input$lognormalProbabilityLower <=
            input$lognormalProbabilityUpper,
          "The lower boundary must not exceed the upper boundary."
        )
      )
      
      return(
        calculate_lognormal_probability(
          calculation_type = "between",
          meanlog = input$lognormalMeanlog,
          sdlog = input$lognormalSdlog,
          lower = input$lognormalProbabilityLower,
          upper = input$lognormalProbabilityUpper
        )
      )
    }
    
    req(input$lognormalProbabilityX)
    
    validate(
      need(
        is.finite(input$lognormalProbabilityX),
        "The boundary must be a finite number."
      )
    )
    
    calculate_lognormal_probability(
      calculation_type = input$lognormalProbabilityType,
      meanlog = input$lognormalMeanlog,
      sdlog = input$lognormalSdlog,
      x = input$lognormalProbabilityX
    )
  })
  
  output$lognormalDistributionPlot <- renderPlotly({
    lognormal_data <- lognormalDistributionData()
    probability_type <- input$lognormalProbabilityType
    
    if (probability_type == "lower") {
      req(input$lognormalProbabilityX)
      
      shaded_data <- lognormal_data[
        lognormal_data$x <= input$lognormalProbabilityX,
        ,
        drop = FALSE
      ]
    } else if (probability_type == "upper") {
      req(input$lognormalProbabilityX)
      
      shaded_data <- lognormal_data[
        lognormal_data$x >= input$lognormalProbabilityX,
        ,
        drop = FALSE
      ]
    } else {
      req(
        input$lognormalProbabilityLower,
        input$lognormalProbabilityUpper
      )
      
      validate(
        need(
          input$lognormalProbabilityLower <=
            input$lognormalProbabilityUpper,
          "The lower boundary must not exceed the upper boundary."
        )
      )
      
      shaded_data <- lognormal_data[
        lognormal_data$x >= input$lognormalProbabilityLower &
          lognormal_data$x <= input$lognormalProbabilityUpper,
        ,
        drop = FALSE
      ]
    }
    
    lognormal_plot <- plot_ly() %>%
      add_trace(
        data = lognormal_data,
        x = ~x,
        y = ~density,
        type = "scatter",
        mode = "lines",
        name = "Density",
        line = list(
          color = input$lognormalCurveColour,
          width = 3
        ),
        hovertemplate = paste(
          "x: %{x:.3f}",
          "<br>Density: %{y:.5f}",
          "<extra></extra>"
        )
      )
    
    if (nrow(shaded_data) > 0) {
      lognormal_plot <- lognormal_plot %>%
        add_trace(
          x = c(
            shaded_data$x[1],
            shaded_data$x,
            shaded_data$x[nrow(shaded_data)]
          ),
          y = c(
            0,
            shaded_data$density,
            0
          ),
          type = "scatter",
          mode = "lines",
          fill = "toself",
          name = "Selected probability",
          line = list(width = 0),
          fillcolor = "rgba(243, 156, 18, 0.35)",
          hoverinfo = "skip"
        )
    }
    
    lognormal_plot %>%
      layout(
        title = paste0(
          "Lognormal Distribution: LN(\u03bc = ",
          input$lognormalMeanlog,
          ", \u03c3 = ",
          input$lognormalSdlog,
          ")"
        ),
        xaxis = list(
          title = "x",
          range = c(
            0,
            max(lognormal_data$x)
          )
        ),
        yaxis = list(
          title = "Probability density",
          rangemode = "tozero"
        ),
        showlegend = TRUE,
        hovermode = "x unified"
      )
  })
  
  output$lognormalProbabilityResult <- renderText({
    paste0(
      "Probability = ",
      formatC(
        lognormalProbabilityValue(),
        format = "f",
        digits = 6
      )
    )
  })
  
  # ------------------------------------------------------------
  # 3.19 Pi estimator
  # ------------------------------------------------------------
  
  piEstimateResults <- reactive({
    calculate_pi_estimate(
      n_trapeziums = input$piTrapeziums
    )
  })
  
  output$piEstimatorPlot <- renderPlotly({
    pi_results <- piEstimateResults()
    
    plot_ly() %>%
      add_trace(
        data = pi_results$trapezium_polygons,
        x = ~x,
        y = ~y,
        type = "scatter",
        mode = "lines",
        fill = "toself",
        fillcolor = "rgba(52, 152, 219, 0.30)",
        line = list(
          color = "rgba(41, 128, 185, 0.85)",
          width = 1
        ),
        name = "Trapeziums",
        hoverinfo = "skip"
      ) %>%
      add_trace(
        data = pi_results$curve_data,
        x = ~x,
        y = ~y,
        type = "scatter",
        mode = "lines",
        line = list(
          color = "#c0392b",
          width = 4
        ),
        name = "y = sqrt(1 - x^2)",
        hovertemplate = paste0(
          "x = %{x:.4f}<br>",
          "y = %{y:.4f}<extra></extra>"
        )
      ) %>%
      layout(
        title = paste0(
          "Quarter Circle with ",
          input$piTrapeziums,
          if (input$piTrapeziums == 1) {
            " Trapezium"
          } else {
            " Trapeziums"
          }
        ),
        xaxis = list(
          title = "x",
          range = c(0, 1.02),
          constrain = "domain"
        ),
        yaxis = list(
          title = "y",
          range = c(0, 1.02),
          scaleanchor = "x",
          scaleratio = 1
        ),
        showlegend = TRUE,
        hovermode = "closest"
      )
  })
  
  output$quarterPiEstimateBox <- renderValueBox({
    valueBox(
      value = formatC(
        piEstimateResults()$quarter_pi_estimate,
        format = "f",
        digits = 10
      ),
      subtitle = "Estimated Area = Estimated pi / 4",
      icon = icon("chart-area"),
      color = "aqua",
      width = 12
    )
  })
  
  output$piEstimateBox <- renderValueBox({
    valueBox(
      value = formatC(
        piEstimateResults()$pi_estimate,
        format = "f",
        digits = 10
      ),
      subtitle = "Estimated pi = 4 x Estimated Area",
      icon = icon("calculator"),
      color = "green",
      width = 12
    )
  })
  
  # ------------------------------------------------------------
  # 3.20 Downloads
  # ------------------------------------------------------------
  
  output$downloadHistogram <- downloadHandler(
    
    filename = function() {
      paste0(
        "MonteCarloHistogram_",
        format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
        ".png"
      )
    },
    
    content = function(file) {
      req(simulationResults())
      
      png(
        filename = file,
        width = 1200,
        height = 800,
        res = 150
      )
      
      p <- create_loss_histogram(
        sim_data = simulationResults()$data,
        bins = input$histogramBins,
        fill_colour = input$histogramFillColour,
        border_colour = input$histogramBorderColour,
        show_density = input$showDensity,
        density_colour = input$densityColour
      )
      
      print(p)
      
      dev.off()
    }
  )
  
  output$downloadGaltonHistogram <- downloadHandler(
    
    filename = function() {
      paste0(
        "GaltonBoardHistogram_",
        format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
        ".png"
      )
    },
    
    content = function(file) {
      req(galtonResults())
      
      galton_data <- galtonResults()$data
      
      p <- create_galton_histogram(
        galton_data = galton_data,
        rows = input$galtonRows,
        p = input$galtonProbability,
        fill_colour = input$galtonFillColour,
        border_colour = input$galtonBorderColour,
        show_theory = input$showGaltonTheory,
        theory_colour = input$galtonTheoryColour
      )
      
      ggsave(
        filename = file,
        plot = p,
        width = 10,
        height = 7,
        dpi = 150
      )
    }
  )
  
  # Full aggregate loss simulation CSV
  output$downloadLossCSV <- downloadHandler(
    
    filename = function() {
      paste0(
        "AggregateLossSimulation_",
        format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
        ".csv"
      )
    },
    
    content = function(file) {
      req(simulationResults())
      
      write.csv(
        simulationResults()$data,
        file,
        row.names = FALSE
      )
    }
  )
  
  # Selected aggregate loss histogram bin CSV
  output$downloadSelectedBinCSV <- downloadHandler(
    
    filename = function() {
      selected <- selectedBinData()
      
      paste0(
        "SelectedAggregateLossBin_",
        round(selected$lower_bound),
        "_to_",
        round(selected$upper_bound),
        "_",
        format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
        ".csv"
      )
    },
    
    content = function(file) {
      selected <- selectedBinData()
      
      write.csv(
        selected$data,
        file,
        row.names = FALSE
      )
    }
  )
  
  # Full Galton board simulation CSV
  output$downloadGaltonCSV <- downloadHandler(
    
    filename = function() {
      paste0(
        "GaltonBoardSimulation_",
        format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
        ".csv"
      )
    },
    
    content = function(file) {
      req(galtonResults())
      
      write.csv(
        galtonResults()$data,
        file,
        row.names = FALSE
      )
    }
  )
  
  # Selected Galton board bar CSV
  output$downloadSelectedGaltonBarCSV <- downloadHandler(
    
    filename = function() {
      selected <- selectedGaltonBarData()
      
      paste0(
        "SelectedGaltonPosition_",
        selected$final_position,
        "_",
        format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
        ".csv"
      )
    },
    
    content = function(file) {
      selected <- selectedGaltonBarData()
      
      write.csv(
        selected$data,
        file,
        row.names = FALSE
      )
    }
  )
}

# ============================================================
# 4. RUN APP
# ============================================================

shinyApp(ui = ui, server = server)