
# Estimates from Scenario A (see Table 2 and Supplementary Table 2), which 
# comprises the combination of parameters considered highest quality and used
# to develop the natural history model.

# Baseline Prevalence

prop_ndbe <- 0.094
prop_lgdbe <- 0.00065
prop_hgdbe <- 0.00195
prop_eac1 <- 0.0023
prop_eac2 <- 0.00038

# Annual Progression

prob_gerd_prog <- 0.01
prob_ndbe_prog <- 0.0395
prob_lgdbe_prog <- 0.1045
prob_hgdbe_prog <- 0.1905
prob_eac1_prog <- 0.33
prob_eac1_diag <- 0.024
prob_eac1_diag_s <- 0.9

# Annual Mortality

prob_eac1_death <- 0.11
prob_eac1_death_s <- 0.037
prob_eac2_y1_prog  <- 0.549
prob_eac2_y2_prog <- 0.423
prob_eac2_y3_prog <- 0.278
prob_eac2_y4_prog <- 0.25
prob_eac2_y5_prog <- 0.103
prob_death <- c(
  0.0082,0.0107,0.012,0.0136,0.0151,0.0167
  ,0.0186,0.0207,0.0230,0.0253,0.0281,0.0312
)

# Screening

sensitivity_ndbe <- 0.8
sensitivity_lgdbe <- 0.85
sensitivity_hgdbe <- 0.85
sensitivity_eac1 <- 0.9
uptake <- 0.65
efficacy <- 0.75

# Additional Parameters

start_n <- 100000
alpha1 <- 0.01
alpha2 <- 0.025
alpha3 <- 0.015
N <- 120000
