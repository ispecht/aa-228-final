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
ALPHA <- 10

N_REPLICATIONS <- 10000

rollout <- function(policy_number) {
  out <- list()
  out$costs <- rep(NA, N_REPLICATIONS)
  out$counts <- rep(NA, N_REPLICATIONS)
  for (i in 1:N_REPLICATIONS) {
    sim <- simulate(policy_number)
    out$costs[i] <- sim$total_cost
    out$counts[i] <- sim$total_cases
  }
  return(out)
}

counts_by_policy <- list()
costs_by_policy <- list()
for(i in 1:4) {
  print(paste0("Running Policy ", i-1))
  sims <- rollout(i-1)
  counts_by_policy[[i]] <- sims$counts
  costs_by_policy[[i]] <- sims$costs
}

N_REPLICATIONS <- 1000





# Usage examples:
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



sim <- simulate(3)
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

ggsave("fig.png", width = 12, height = 6)


