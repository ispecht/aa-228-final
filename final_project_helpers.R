

# Adjust case counts based on policy
adjust_cases <- function(
    cases, 
    policy_number, 
    p_diagnostic,
    p_genomic,
    prev_positive_count, # By location
    prev_R_estimate, # By location
    alpha
){
  
  out <- list()
  out$cases <- cases
  out$cost <- sum(cases) * ILLNESS_COST
  
  # Test positivity count by location
  out$positive_count <- rep(0, N_LOCATIONS)
  
  # Reproductive number estimate by location
  out$R_estimate <- rep(0, N_LOCATIONS)
  
  # Number of positive cases by location
  n_cases_by_location <- rep(0, N_LOCATIONS)
  
  # Number of genomes per location
  n_genomes_by_location <- rep(0, N_LOCATIONS)
  
  # Policy 0: do nothing
  if(policy_number == 0) {
    return(out)
  }
  
  # If policy 1, we disregard the previous positive count
  if(policy_number == 1) {
    prev_positive_count <- rep(0, N_LOCATIONS)
  }
  
  # Fix alpha to produce valid probabilities
  while(p_diagnostic * max(prev_positive_count * prev_R_estimate + alpha) / mean(prev_positive_count * prev_R_estimate + alpha) > 1) {
    alpha <- alpha * 1.1
  }
  
  
  for(i in 1:length(cases)) {
    
    location = ((i - 1) %% 8) + 1
    
    n_cases_by_location[location] <- n_cases_by_location[location] + cases[i]
    
    
    effective_p_diagnostic <- p_diagnostic * ((prev_positive_count[location] * prev_R_estimate[location] + alpha) / mean(prev_positive_count * prev_R_estimate + alpha)) # Use proportion based on pseudocounts
    
    
    n_positive_tested <- rbinom(1, cases[i], effective_p_diagnostic)
    out$cost <- out$cost + (n_positive_tested * DIAGNOSTIC_COST)
    
    n_positive <- rbinom(1, n_positive_tested, DIAGNOSTIC_SENSITIVITY)
    out$positive_count[location] <- out$positive_count[location] + n_positive
    
    # Number that get tested and obey quarantine
    n_quarantined <- rbinom(1, n_positive, P_QUARANTINE)
    out$cases[i] <- out$cases[i] - n_quarantined
    out$cost <- out$cost + (n_quarantined * QUARANTINE_COST)
    
    # Genomic testing
    # Of the positive ones, get a genomic test with probability p_genomic
    if(policy_number == 3) {
      n_genomes <- rbinom(1, n_positive, p_genomic)
      
      # Get the reproductive number
      R <- R_v[ceiling(i / 8)]
      
      # Update estimated R
      out$R_estimate[location] <- out$R_estimate[location] + sum(rnorm(n_genomes, R, SEQUENCING_EPSILON))
      n_genomes_by_location[location] <- n_genomes_by_location[location] + n_genomes
      
      # Cost of this
      out$cost <- out$cost + (n_genomes * SEQUENCING_COST)
    }
  }
  
  # Negative tests by location
  for(location in 1:N_LOCATIONS) {
    effective_p_diagnostic <- p_diagnostic * ((prev_positive_count[location] * prev_R_estimate[location] + alpha) / mean(prev_positive_count * prev_R_estimate + alpha)) # Use proportion based on pseudocounts
    
    n_negative_tested <- rbinom(1, N_l[location] - n_cases_by_location[location], effective_p_diagnostic)
    
    out$cost <- out$cost + (n_negative_tested * DIAGNOSTIC_COST)
    
    # Reproductive number estimation
    # If no data, set to 1
    if(n_genomes_by_location[location] == 0) {
      out$R_estimate[location] <- 1
    }else{
      out$R_estimate[location] <- out$R_estimate[location] / n_genomes_by_location[location]
    }
  }
  
  return(out)
  
}


# Function to return the cases at the next generation
update_cases <- function(cases) {
  new_cases <- rep(0, length(cases))
  
  for(i in 1:length(cases)) {
    
    # Get the reproductive number
    R <- R_v[ceiling(i / 8)]
    
    # Number of offspring
    n_offspring <- rpois(1, R * cases[i])
    
    # Put them into their respective categories
    distribution <- rmultinom(1, n_offspring, Q[i, ])[,1]
    
    # Update number of new cases
    new_cases <- new_cases + distribution
  }
  
  new_cases
}

simulate <- function(policy_number) {
  
  cases <- INIT_CASES
  
  case_table <- matrix(0, nrow = N_GENERATIONS, ncol = N_VARIANTS * N_LOCATIONS)
  case_table[1, ] <- cases
  
  prev_positive_count <- rep(0, N_LOCATIONS)
  prev_R_estimate <- rep(1, N_LOCATIONS)
  
  total_cost <- 0
  
  
  
  for (i in 1:N_GENERATIONS) {
    
    adjusted_cases <- adjust_cases(
      cases,
      policy_number,
      P_DIAGNOSTIC,
      P_GENOMIC,
      prev_positive_count,
      prev_R_estimate,
      ALPHA
    )
    
    cases <- adjusted_cases$cases
    total_cost <- total_cost + adjusted_cases$cost
    prev_positive_count <- adjusted_cases$positive_count
    prev_R_estimate <- adjusted_cases$R_estimate
    
    if(i < N_GENERATIONS) {
      cases <- update_cases(cases)
      case_table[i + 1, ] <- cases
    }
    
  }
  
  rownames(case_table) <- NULL
  
  out <- list()
  out$case_table <- case_table
  out$total_cost <- total_cost
  out$total_cases <- sum(case_table)
  
  return(out)
}


# Function to plot cases by location over time
plot_cases_by_location <- function(case_table) {
  # Convert case_table to data frame with generation column
  n_generations <- nrow(case_table)
  
  # Extract location-level data (sum across variants)
  location_data <- matrix(0, nrow = n_generations, ncol = N_LOCATIONS)
  
  for(i in 1:n_generations) {
    for(loc in 1:N_LOCATIONS) {
      # Sum across all variants for this location
      indices <- seq(loc, N_LOCATIONS * N_VARIANTS, by = N_LOCATIONS)
      location_data[i, loc] <- sum(case_table[i, indices])
    }
  }
  
  # Convert to long format for ggplot
  df <- as.data.frame(location_data)
  colnames(df) <- paste0("Location_", 1:N_LOCATIONS)
  df$Generation <- 1:n_generations
  
  df_long <- df %>%
    pivot_longer(cols = starts_with("Location_"),
                 names_to = "Location",
                 values_to = "Cases")
  
  # Create stacked area plot
  ggplot(df_long, aes(x = Generation, y = Cases, fill = Location)) +
    geom_area(alpha = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    scale_x_continuous(breaks = 1:n_generations) +
    labs(title = "Cases by Geographic Location",
         x = "Generation",
         y = "Number of Cases",
         fill = "Location") +
    theme_minimal() +
    theme(legend.position = "right")
}

# Function to plot cases by genetic variant over time
plot_cases_by_variant <- function(case_table) {
  # Convert case_table to data frame with generation column
  n_generations <- nrow(case_table)
  
  # Extract variant-level data (sum across locations)
  variant_data <- matrix(0, nrow = n_generations, ncol = N_VARIANTS)
  
  for(i in 1:n_generations) {
    for(var in 1:N_VARIANTS) {
      # Sum across all locations for this variant
      # Variant v occupies positions (v-1)*N_LOCATIONS + 1 through v*N_LOCATIONS
      start_idx <- (var - 1) * N_LOCATIONS + 1
      end_idx <- var * N_LOCATIONS
      variant_data[i, var] <- sum(case_table[i, start_idx:end_idx])
    }
  }
  
  # Convert to long format for ggplot
  df <- as.data.frame(variant_data)
  colnames(df) <- paste0("Variant_", 1:N_VARIANTS)
  df$Generation <- 1:n_generations
  
  df_long <- df %>%
    pivot_longer(cols = starts_with("Variant_"),
                 names_to = "Variant",
                 values_to = "Cases")
  
  # Add R values to variant labels for context
  df_long$Variant <- factor(df_long$Variant,
                            levels = paste0("Variant_", 1:N_VARIANTS),
                            labels = paste0("Variant ", 1:N_VARIANTS, 
                                            " (R=", R_v, ")"))
  
  # Create stacked area plot
  ggplot(df_long, aes(x = Generation, y = Cases, fill = Variant)) +
    geom_area(alpha = 0.7) +
    scale_fill_brewer(palette = "Set1") +
    scale_x_continuous(breaks = 1:n_generations) +
    labs(title = "Cases by Genetic Variant",
         x = "Generation",
         y = "Number of Cases",
         fill = "Variant") +
    theme_minimal() +
    theme(legend.position = "right")
}



# Function to plot overlapping density plots for policy comparison
plot_policy_comparison <- function(data_by_policy, metric_name = "Total Cost", 
                                   policy_names = NULL, log_scale = TRUE) {
  
  # Default policy names if not provided
  if(is.null(policy_names)) {
    policy_names <- paste0("Policy ", 0:(length(data_by_policy) - 1))
  }
  
  # Check if data is negative (costs)
  is_negative <- all(sapply(data_by_policy, function(x) all(x <= 0)))
  
  # Convert list to data frame in long format
  df_list <- list()
  for(i in 1:length(data_by_policy)) {
    df_list[[i]] <- data.frame(
      Value = data_by_policy[[i]],
      Policy = policy_names[i]
    )
  }
  df <- bind_rows(df_list)
  
  # Convert Policy to factor to control order
  df$Policy <- factor(df$Policy, levels = policy_names)
  
  # If negative and log scale requested, take absolute value for plotting
  if(is_negative && log_scale) {
    df$PlotValue <- abs(df$Value)
    x_label <- paste0(metric_name)
  } else {
    df$PlotValue <- df$Value
    x_label <- metric_name
  }
  
  # Create overlapping density plot
  p <- ggplot(df, aes(x = PlotValue, fill = Policy, color = Policy)) +
    geom_density(alpha = 0.4, linewidth = 1) +
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1") +
    labs(title = paste0("Distribution of ", metric_name, " by Policy"),
         subtitle = paste0("Based on ", length(data_by_policy[[1]]), " replications per policy"),
         x = x_label,
         y = "Density",
         fill = "Policy",
         color = "Policy") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10))
  
  # Add log scale if requested
  if(log_scale) {
    p <- p + scale_x_log10(labels = scales::comma)
  }
  
  return(p)
}
