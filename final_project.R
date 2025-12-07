## Final project
## Ivan Specht
## AA 228 / CS 238 - Decision Making Under Uncertainty

## We are going to model an epidemic and use decision making to determine the optimal course of action given the data

# Seed
library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)
source("final_project_helpers.R")
set.seed(228)

## Genetic variants

# Number of variants we model
N_VARIANTS <- 4

# Reproductive numbers of each variant
# We make one variant have a very high reproductive numner, so that the reward of detecting it is particularly noticeable
R_v <- c(1.1, 2.5, 1.1, 1.1)

# Transition probabilities between variants per generation
P_MUTATE <- 0.01
Q_v <- matrix(P_MUTATE, nrow = N_VARIANTS, ncol = N_VARIANTS)
diag(Q_v) <- 1 - P_MUTATE * (N_VARIANTS - 1)

## Locations
N_LOCATIONS <- 8

# Total poputation by location
N_l <- rep(100000, N_LOCATIONS)
TOTAL_POP <- sum(N_l)

# Transition probabilities between locations by generation
P_MIGRATE <- 0.01
Q_l <- matrix(P_MIGRATE, nrow = N_LOCATIONS, ncol = N_LOCATIONS)
diag(Q_l) <- 1 - P_MIGRATE * (N_LOCATIONS - 1)

# Transition matrix for variants and geography
Q <- kronecker(Q_v, Q_l)

## Actions
# Can diagnostic test one or more randomly-selected person in a given geographic region
# Can genome sequence someone who tests positive
# Can quarantine someone who tests positive

## Observations
# Diagnostic tests have 95% sensitivity, 100% specificity
P_DIAGNOSTIC <- 0.2
DIAGNOSTIC_SENSITIVITY = 0.95

# Genome sequencing can determine the infectiousness (reproductive number) of the variant, plus some normal noise with SD epsilon
P_GENOMIC <- 0.2
SEQUENCING_EPSILON = 0.05

# Probability of obeying quarantine
P_QUARANTINE = 0.95

## Rewards
DIAGNOSTIC_COST = -0.5
SEQUENCING_COST = -1
QUARANTINE_COST = -1.5

# Then each infection incurs a steep cost
ILLNESS_COST = -100

## Policies
# We're going to study a few possible policies via simulation

# Policy 0: do nothing at all

# Policy 1: uninformed testing. Test random people at given testing rate, quarantine positives

# Policy 2: informed testing: Test rate depends on positivity rate from previous generation. Quarantine positives 

# Policy 3: genomically-informed testing: Allocate 800 tests to geographies based on positive rate times estimated reproductive number, inferred from 80 genomic tests. Quarantine positives

# Number of generations to simulate
N_GENERATIONS = 10


POLICY_NUMBER <- 0


# Assume the outbreak starts with 10 cases of variant 1 in geography 1
INIT_CASES <- rep(0, N_LOCATIONS * N_VARIANTS)
INIT_CASES[1:8] <- 10
INIT_CASES[3] <- 0
INIT_CASES[11] <- 10

# Basically a prior on case counts by location
ALPHA <- 1

N_REPLICATIONS <- 10000

rollout <- function(policy_number, P_DIAGNOSTIC, P_GENOMIC, N_REPLICATIONS, INIT_CASES, N_GENERATIONS, ALPHA) {
  out <- list()
  out$costs <- rep(NA, N_REPLICATIONS)
  out$counts <- rep(NA, N_REPLICATIONS)
  for (i in 1:N_REPLICATIONS) {
    sim <- simulate(policy_number, P_DIAGNOSTIC, P_GENOMIC, INIT_CASES, N_GENERATIONS, ALPHA)
    out$costs[i] <- sim$total_cost
    out$counts[i] <- sim$total_cases
  }
  return(out)
}

counts_by_policy <- list()
costs_by_policy <- list()
for(i in 1:4) {
  print(paste0("Running Policy ", i-1))
  sims <- rollout(i-1, P_DIAGNOSTIC, P_GENOMIC, N_REPLICATIONS, INIT_CASES, N_GENERATIONS, ALPHA)
  counts_by_policy[[i]] <- sims$counts
  costs_by_policy[[i]] <- sims$costs
}

# Plot total costs
p1 <- plot_policy_comparison(costs_by_policy, 
                             metric_name = "Total Cost",
                             policy_names = c("Policy 0: Do Nothing",
                                              "Policy 1: Uninformed Testing",
                                              "Policy 2: Informed Testing", 
                                              "Policy 3: Genomically-Informed"))

# Plot total case counts
p2 <- plot_policy_comparison(counts_by_policy,
                             metric_name = "Total Cases",
                             policy_names = c("Policy 0: Do Nothing",
                                              "Policy 1: Uninformed Testing",
                                              "Policy 2: Informed Testing",
                                              "Policy 3: Genomically-Informed"))



sim <- simulate(3, P_DIAGNOSTIC, P_GENOMIC, INIT_CASES, N_GENERATIONS, ALPHA)
p3 <- plot_cases_by_location(sim$case_table)
p4 <- plot_cases_by_variant(sim$case_table)

plot_grid(
  p1,
  p2,
  p3,
  p4,
  nrow = 2,
  labels = "AUTO"
)

ggsave("fig1.png", width = 12, height = 6)


# Try varying diagnostic testing rate, holding genome sequencing rate constant
costs_by_testing <- list()
p_diag_vals <- seq(0, 1, 0.1)
for(i in 1:length(p_diag_vals)) {
  print(i)
  sims <- rollout(3, p_diag_vals[i], 0.2, 10000, INIT_CASES, N_GENERATIONS, ALPHA)
  costs_by_testing[[i]] <- sims$costs
}

# Try varying genome sequencing rate, holding diagnostic rate constant
costs_by_sequencing <- c()
p_geno_vals <- seq(0, 1, 0.1)
for(i in 1:length(p_geno_vals)) {
  print(i)
  sims <- rollout(3, 0.2, p_geno_vals[i], 10000, INIT_CASES, N_GENERATIONS, ALPHA)
  costs_by_sequencing[[i]] <- sims$costs
}


# Plot effect of diagnostic testing rate on costs
p5 <- plot_parameter_violin(costs_by_testing,
                      param_values = p_diag_vals,
                      param_name = "Diagnostic Testing Rate",
                      metric_name = "Total Cost",
                      log_scale = TRUE)

# Plot effect of genome sequencing rate on costs
p6 <- plot_parameter_violin(costs_by_sequencing,
                      param_values = p_geno_vals,
                      param_name = "Genome Sequencing Rate",
                      metric_name = "Total Cost",
                      log_scale = TRUE)

plot_grid(
  p5,
  p6,
  nrow = 1,
  labels = "AUTO"
)

ggsave("fig2.png", width = 12, height = 3)

# Finally, assuming constant testing and sequencing, learn the optimal balance between random and targeted sequencing by tuning alpha
alpha_vals <- seq(0.2, 2.2, 0.2)

init_cases <- INIT_CASES
for(i in 1:N_GENERATIONS) {

  reward_by_alpha <- rep(NA, length(alpha_vals))
  for(j in 1:length(alpha_vals)) {
    print(paste0("Testing alpha ", alpha_vals[j]))
    reward_by_alpha[j] <- mean(
      rollout(3, P_DIAGNOSTIC, P_GENOMIC, 1000, init_cases, N_GENERATIONS - i + 1, alpha_vals[j])$costs
    )
    
  }
  
  best_alpha <- alpha_vals[which.max(reward_by_alpha)]
  print(paste0("Best alpha: ", best_alpha))
  
  # Simulate once to move to the next step
  if(i < N_GENERATIONS) {
    sim <- simulate(3, P_DIAGNOSTIC, P_GENOMIC, init_cases, N_GENERATIONS - i + 1, best_alpha)
    init_cases <- sim$case_table[2, ]
  }
}






