
# Simple semi-Markov natural history model of EAC, developed to understand:
# (i) the impact of novel screening strategies on EAC mortality and advanced
# stage incidence; and (ii) inform power calculations for the BEST4 trial.

# SETUP =======================================================================

# Read parameter point estimates used for model development.

# rm(list = ls()); source('parameter_point_estimates.R')

options(scipen = 999)
set.seed(35729863)

# MODEL =======================================================================

eac_model <- function(
    prop_ndbe
    ,prop_lgdbe
    ,prop_hgdbe
    ,prop_eac1
    ,prop_eac2
    ,prob_gerd_prog
    ,prob_ndbe_prog
    ,prob_lgdbe_prog
    ,prob_hgdbe_prog
    ,prob_eac1_prog
    ,prob_eac1_diag
    ,prob_eac1_diag_s
    ,prob_eac1_death
    ,prob_eac1_death_s
    ,prob_eac2_y1_prog
    ,prob_eac2_y2_prog
    ,prob_eac2_y3_prog
    ,prob_eac2_y4_prog
    ,prob_eac2_y5_prog
    ,prob_death
    ,sensitivity_ndbe
    ,sensitivity_lgdbe
    ,sensitivity_hgdbe
    ,sensitivity_eac1
    ,uptake
    ,efficacy
    ,start_n
    ,alpha1
    ,alpha2
    ,alpha3
    ,N
) {
  
  # Assign number of people in each stage of disease at time 0 in the absence
  # of screening
  
  # First, proportion of GERD; based on proportion in other stages of disease
  
  prop_gerd <- 1 - prop_ndbe - prop_lgdbe - prop_hgdbe - prop_eac1 - prop_eac2
  
  # Second, number of people in each stage of disease
  
  gerd0 <- start_n * prop_gerd
  ndbe0 <- start_n * prop_ndbe
  lgdbe0 <- start_n * prop_lgdbe
  hgdbe0 <- start_n * prop_hgdbe
  eac1_0 <- start_n * prop_eac1
  eac2_0 <- start_n * prop_eac2
  
  # Third, create holder for for-loop
  
  ns <- data.frame(
    gerd = c(gerd0,rep(NA,12))
    ,ndbe = c(ndbe0,rep(NA,12))
    ,lgdbe = c(lgdbe0,rep(NA,12))
    ,hgdbe = c(hgdbe0,rep(NA,12))
    ,und_eac1_a = c(eac1_0,rep(NA,12))
    ,und_eac1_b = c(0,rep(NA,12))
    ,und_eac1_c = c(0,rep(NA,12))
    ,und_eac1_d = c(0,rep(NA,12))
    ,d_eac1_a = c(0,rep(NA,12))
    ,d_eac1_b = c(0,rep(NA,12))
    ,d_eac1_c = c(0,rep(NA,12))
    ,d_eac1_d = c(0,rep(NA,12))
    ,d_eac1_e = c(0,rep(NA,12))
    ,d_eac1_f = c(0,rep(NA,12))
    ,d_eac1_g = c(0,rep(NA,12))
    ,dead_eac1 = c(0,rep(NA,12))
    ,eac2_a = c(eac2_0,rep(NA,12))
    ,eac2_b = c(0,rep(NA,12))
    ,eac2_c = c(0,rep(NA,12))
    ,eac2_d = c(0,rep(NA,12))
    ,eac2_e = c(0,rep(NA,12))
    ,eac2_f = c(0,rep(NA,12))
    ,eac2_g = c(0,rep(NA,12))
    ,dead_eac2 = c(0,rep(NA,12))
    ,dead_other = c(0,rep(NA,12))
  )
  
  # Finally, simulate number of people in each stage of disease each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    # GERD through HGD:
    
    ns$gerd[j] <- (ns$gerd[i]*(1-prob_death[i]))*(1-prob_gerd_prog)
    ns$ndbe[j] <- ((ns$ndbe[i]*(1-prob_death[i]))*(1-prob_ndbe_prog))+((ns$gerd[i]*(1-prob_death[i]))*prob_gerd_prog)
    ns$lgdbe[j] <- ((ns$lgdbe[i]*(1-prob_death[i]))*(1-prob_lgdbe_prog))+((ns$ndbe[i]*(1-prob_death[i]))*prob_ndbe_prog)
    ns$hgdbe[j] <- ((ns$hgdbe[i]*(1-prob_death[i]))*(1-prob_hgdbe_prog))+((ns$lgdbe[i]*(1-prob_death[i]))*prob_lgdbe_prog)
    
    # Undiagnosed EAC Stage 1:
    # Patients can be diagnosed in years 2 or 3, else they progress to Stage 2+
    # or stay undiagnosed. Patients are not allowed to die from undiagnosed EAC
    # Stage 1.
    
    ns$und_eac1_a[j] <- ((ns$hgdbe[i]*(1-prob_death[i]))*prob_hgdbe_prog)
    ns$und_eac1_b[j] <- ((ns$und_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_prog)*(1-prob_eac1_diag))
    ns$und_eac1_c[j] <- ((ns$und_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_prog)*(1-prob_eac1_diag))
    ns$und_eac1_d[j] <- ((ns$und_eac1_c[i]*(1-prob_death[i])*(1-prob_eac1_prog)))+((ns$und_eac1_d[i]*(1-prob_death[i])*(1-prob_eac1_prog)))
    
    # Diagnosed EAC Stage 1:
    
    ns$d_eac1_a[j] <- (((ns$und_eac1_a[i]+ns$und_eac1_b[i])*(1-prob_death[i])))*prob_eac1_diag                                 
    ns$d_eac1_b[j] <- (ns$d_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    ns$d_eac1_c[j] <- (ns$d_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    ns$d_eac1_d[j] <- (ns$d_eac1_c[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    ns$d_eac1_e[j] <- (ns$d_eac1_d[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    ns$d_eac1_f[j] <- (ns$d_eac1_e[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    ns$d_eac1_g[j] <- ((ns$d_eac1_f[i]*(1-prob_death[i])))+((ns$d_eac1_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 1:
    
    ns$dead_eac1[j] <- (prob_eac1_death*(ns$d_eac1_a[i]+ns$d_eac1_b[i]+ns$d_eac1_c[i]+ns$d_eac1_d[i]+ns$d_eac1_e[i])*(1-prob_death[i]))+ns$dead_eac1[i]
    
    # EAC Stage 2+:
    
    ns$eac2_a[j] <- ((((ns$und_eac1_a[i]+ns$und_eac1_b[i])*(1-prob_eac1_diag))+ns$und_eac1_c[i]+ns$und_eac1_d[i])*(1-prob_death[i])*prob_eac1_prog)
    ns$eac2_b[j] <- ((ns$eac2_a[i]*(1-prob_death[i]))*(1-prob_eac2_y1_prog))
    ns$eac2_c[j] <- ((ns$eac2_b[i]*(1-prob_death[i]))*(1-prob_eac2_y2_prog))
    ns$eac2_d[j] <- ((ns$eac2_c[i]*(1-prob_death[i]))*(1-prob_eac2_y3_prog))
    ns$eac2_e[j] <- ((ns$eac2_d[i]*(1-prob_death[i]))*(1-prob_eac2_y4_prog))
    ns$eac2_f[j] <- ((ns$eac2_e[i]*(1-prob_death[i]))*(1-prob_eac2_y5_prog))
    ns$eac2_g[j] <- ((ns$eac2_f[i]*(1-prob_death[i])))+((ns$eac2_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 2+:
    
    ns$dead_eac2[j] <- (prob_eac2_y1_prog*ns$eac2_a[i]+prob_eac2_y2_prog*ns$eac2_b[i]+prob_eac2_y3_prog*ns$eac2_c[i]+prob_eac2_y4_prog*ns$eac2_d[i]+prob_eac2_y5_prog*ns$eac2_e[i])*(1-prob_death[i])+ns$dead_eac2[i]
    
    # Cumulative death from other causes:
    
    ns$dead_other[j] <- (prob_death[i]*
                           (ns$gerd[i]+ns$ndbe[i]+ns$lgdbe[i]+ns$hgdbe[i]+
                              ns$und_eac1_a[i]+ns$und_eac1_b[i]+ns$und_eac1_c[i]+ns$und_eac1_d[i]+ns$d_eac1_a[i]+ns$d_eac1_b[i]+ns$d_eac1_c[i]+ns$d_eac1_d[i]+ns$d_eac1_e[i]+ns$d_eac1_f[i]+ns$d_eac1_g[i]+
                              ns$eac2_a[i]+ns$eac2_b[i]+ns$eac2_c[i]+ns$eac2_d[i]+ns$eac2_e[i]+ns$eac2_f[i]+ns$eac2_g[i])
    )+ns$dead_other[i]
    
  }
  
  # Simulate what happens to people with GERD at baseline, none of whom would
  # benefit from screening. Therefore, we do not need to separately model
  # what would happen in the presence and absence of screening.
  
  # First, start by resetting select parameters
  
  gerd0 <- ns[1,1]
  prop_gerd <- 1
  prop_ndbe <- prop_lgdbe <- prop_hgdbe <- prop_eac1 <- prop_eac2 <- 0
  ndbe0 <- lgdbe0 <- hgdbe0 <- eac1_0 <- eac2_0 <- 0
  
  # Second, create holder for for-loop
  
  gerd_ns <- data.frame(
    gerd = c(gerd0,rep(NA,12))
    ,ndbe = c(ndbe0,rep(NA,12))
    ,lgdbe = c(lgdbe0,rep(NA,12))
    ,hgdbe = c(hgdbe0,rep(NA,12))
    ,und_eac1_a = c(eac1_0,rep(NA,12))
    ,und_eac1_b = c(0,rep(NA,12))
    ,und_eac1_c = c(0,rep(NA,12))
    ,und_eac1_d = c(0,rep(NA,12))
    ,d_eac1_a = c(0,rep(NA,12))
    ,d_eac1_b = c(0,rep(NA,12))
    ,d_eac1_c = c(0,rep(NA,12))
    ,d_eac1_d = c(0,rep(NA,12))
    ,d_eac1_e = c(0,rep(NA,12))
    ,d_eac1_f = c(0,rep(NA,12))
    ,d_eac1_g = c(0,rep(NA,12))
    ,dead_eac1 = c(0,rep(NA,12))
    ,eac2_a = c(eac2_0,rep(NA,12))
    ,eac2_b = c(0,rep(NA,12))
    ,eac2_c = c(0,rep(NA,12))
    ,eac2_d = c(0,rep(NA,12))
    ,eac2_e = c(0,rep(NA,12))
    ,eac2_f = c(0,rep(NA,12))
    ,eac2_g = c(0,rep(NA,12))
    ,dead_eac2 = c(0,rep(NA,12))
    ,dead_other = c(0,rep(NA,12))
  )
  
  # Finally, simulate number of people in each stage of disease each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    # GERD through HGD:
    
    gerd_ns$gerd[j] <- (gerd_ns$gerd[i]*(1-prob_death[i]))*(1-prob_gerd_prog)
    gerd_ns$ndbe[j] <- ((gerd_ns$ndbe[i]*(1-prob_death[i]))*(1-prob_ndbe_prog))+((gerd_ns$gerd[i]*(1-prob_death[i]))*prob_gerd_prog)
    gerd_ns$lgdbe[j] <- ((gerd_ns$lgdbe[i]*(1-prob_death[i]))*(1-prob_lgdbe_prog))+((gerd_ns$ndbe[i]*(1-prob_death[i]))*prob_ndbe_prog)
    gerd_ns$hgdbe[j] <- ((gerd_ns$hgdbe[i]*(1-prob_death[i]))*(1-prob_hgdbe_prog))+((gerd_ns$lgdbe[i]*(1-prob_death[i]))*prob_lgdbe_prog)
    
    # Undiagnosed EAC Stage 1:
    # Patients can be diagnosed in years 2 or 3, else they progress to Stage 2
    # or stay undiagnosed. Patients are not allowed to die from undiagnosed EAC
    # Stage 1.
    
    gerd_ns$und_eac1_a[j] <- ((gerd_ns$hgdbe[i]*(1-prob_death[i]))*prob_hgdbe_prog)
    gerd_ns$und_eac1_b[j] <- ((gerd_ns$und_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_prog)*(1-prob_eac1_diag))
    gerd_ns$und_eac1_c[j] <- ((gerd_ns$und_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_prog)*(1-prob_eac1_diag))
    gerd_ns$und_eac1_d[j] <- ((gerd_ns$und_eac1_c[i]*(1-prob_death[i])*(1-prob_eac1_prog)))+((gerd_ns$und_eac1_d[i]*(1-prob_death[i])*(1-prob_eac1_prog)))
    
    # Diagnosed EAC Stage 1:
    
    gerd_ns$d_eac1_a[j] <- (((gerd_ns$und_eac1_a[i]+gerd_ns$und_eac1_b[i])*(1-prob_death[i])))*prob_eac1_diag                                 
    gerd_ns$d_eac1_b[j] <- (gerd_ns$d_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    gerd_ns$d_eac1_c[j] <- (gerd_ns$d_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    gerd_ns$d_eac1_d[j] <- (gerd_ns$d_eac1_c[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    gerd_ns$d_eac1_e[j] <- (gerd_ns$d_eac1_d[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    gerd_ns$d_eac1_f[j] <- (gerd_ns$d_eac1_e[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    gerd_ns$d_eac1_g[j] <- ((gerd_ns$d_eac1_f[i]*(1-prob_death[i])))+((gerd_ns$d_eac1_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 1:
    
    gerd_ns$dead_eac1[j] <- (prob_eac1_death*(gerd_ns$d_eac1_a[i]+gerd_ns$d_eac1_b[i]+gerd_ns$d_eac1_c[i]+gerd_ns$d_eac1_d[i]+gerd_ns$d_eac1_e[i])*(1-prob_death[i]))+gerd_ns$dead_eac1[i]
    
    # EAC Stage 2+:
    
    gerd_ns$eac2_a[j] <- ((((gerd_ns$und_eac1_a[i]+gerd_ns$und_eac1_b[i])*(1-prob_eac1_diag))+gerd_ns$und_eac1_c[i]+gerd_ns$und_eac1_d[i])*(1-prob_death[i])*prob_eac1_prog)
    gerd_ns$eac2_b[j] <- ((gerd_ns$eac2_a[i]*(1-prob_death[i]))*(1-prob_eac2_y1_prog))
    gerd_ns$eac2_c[j] <- ((gerd_ns$eac2_b[i]*(1-prob_death[i]))*(1-prob_eac2_y2_prog))
    gerd_ns$eac2_d[j] <- ((gerd_ns$eac2_c[i]*(1-prob_death[i]))*(1-prob_eac2_y3_prog))
    gerd_ns$eac2_e[j] <- ((gerd_ns$eac2_d[i]*(1-prob_death[i]))*(1-prob_eac2_y4_prog))
    gerd_ns$eac2_f[j] <- ((gerd_ns$eac2_e[i]*(1-prob_death[i]))*(1-prob_eac2_y5_prog))
    gerd_ns$eac2_g[j] <- ((gerd_ns$eac2_f[i]*(1-prob_death[i])))+((gerd_ns$eac2_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 2+:
    
    gerd_ns$dead_eac2[j] <- (prob_eac2_y1_prog*gerd_ns$eac2_a[i]+prob_eac2_y2_prog*gerd_ns$eac2_b[i]+prob_eac2_y3_prog*gerd_ns$eac2_c[i]+prob_eac2_y4_prog*gerd_ns$eac2_d[i]+prob_eac2_y5_prog*gerd_ns$eac2_e[i])*(1-prob_death[i])+gerd_ns$dead_eac2[i]
    
    # Cumulative death from other causes:
    
    gerd_ns$dead_other[j] <- (prob_death[i]*
                                (gerd_ns$gerd[i]+gerd_ns$ndbe[i]+gerd_ns$lgdbe[i]+gerd_ns$hgdbe[i]+
                                   gerd_ns$und_eac1_a[i]+gerd_ns$und_eac1_b[i]+gerd_ns$und_eac1_c[i]+gerd_ns$und_eac1_d[i]+gerd_ns$d_eac1_a[i]+gerd_ns$d_eac1_b[i]+gerd_ns$d_eac1_c[i]+gerd_ns$d_eac1_d[i]+gerd_ns$d_eac1_e[i]+gerd_ns$d_eac1_f[i]+gerd_ns$d_eac1_g[i]+
                                   gerd_ns$eac2_a[i]+gerd_ns$eac2_b[i]+gerd_ns$eac2_c[i]+gerd_ns$eac2_d[i]+gerd_ns$eac2_e[i]+gerd_ns$eac2_f[i]+gerd_ns$eac2_g[i])
    )+gerd_ns$dead_other[i]
    
  }
  
  # Simulate what happens to people with NDBE at at baseline, first among those
  # who screen positive.
  
  # First, start by resetting select parameters
  
  ndbe0 <- ns[1,2] * sensitivity_ndbe
  prop_ndbe <- 1
  prop_lgdbe <- prop_hgdbe <- prop_eac1 <- prop_eac2 <- 0
  lgdbe0 <- hgdbe0 <- eac1_0 <- eac2_0 <- 0
  
  # Second, create holder for for-loop
  
  ndbe_s <- data.frame(
    gerd = c(rep(0,13))
    ,ndbe = c(ndbe0,rep(NA,12))
    ,lgdbe = c(lgdbe0,rep(NA,12))
    ,hgdbe = c(hgdbe0,rep(NA,12))
    ,und_eac1_a = c(eac1_0,rep(NA,12))
    ,und_eac1_b = c(0,rep(NA,12))
    ,und_eac1_c = c(0,rep(NA,12))
    ,und_eac1_d = c(0,rep(NA,12))
    ,d_eac1_a = c(0,rep(NA,12))
    ,d_eac1_b = c(0,rep(NA,12))
    ,d_eac1_c = c(0,rep(NA,12))
    ,d_eac1_d = c(0,rep(NA,12))
    ,d_eac1_e = c(0,rep(NA,12))
    ,d_eac1_f = c(0,rep(NA,12))
    ,d_eac1_g = c(0,rep(NA,12))
    ,dead_eac1 = c(0,rep(NA,12))
    ,eac2_a = c(eac2_0,rep(NA,12))
    ,eac2_b = c(0,rep(NA,12))
    ,eac2_c = c(0,rep(NA,12))
    ,eac2_d = c(0,rep(NA,12))
    ,eac2_e = c(0,rep(NA,12))
    ,eac2_f = c(0,rep(NA,12))
    ,eac2_g = c(0,rep(NA,12))
    ,dead_eac2 = c(0,rep(NA,12))
    ,dead_other = c(0,rep(NA,12))
  )
  
  # Finally, simulate number of people in each stage of disease each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    # NDBE through HGD:
    
    ndbe_s$ndbe[j] <- ((ndbe_s$ndbe[i]*(1-prob_death[i]))*(1-prob_ndbe_prog))
    ndbe_s$lgdbe[j] <- ((ndbe_s$lgdbe[i]*(1-prob_death[i]))*(1-(prob_lgdbe_prog*(1-efficacy))))+((ndbe_s$ndbe[i]*(1-prob_death[i]))*prob_ndbe_prog)
    ndbe_s$hgdbe[j] <- ((ndbe_s$hgdbe[i]*(1-prob_death[i]))*(1-(prob_hgdbe_prog*(1-efficacy))))+((ndbe_s$lgdbe[i]*(1-prob_death[i]))*prob_lgdbe_prog*(1-efficacy))
    
    # Undiagnosed EAC Stage 1:
    # Patients can be diagnosed in years 2 or 3, else they progress to Stage 2
    # or stay undiagnosed. Patients are not allowed to die from undiagnosed EAC
    # Stage 1.
    
    ndbe_s$und_eac1_a[j] <- ((ndbe_s$hgdbe[i]*(1-prob_death[i]))*prob_hgdbe_prog*(1-efficacy))
    ndbe_s$und_eac1_b[j] <- ((ndbe_s$und_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_diag_s)*(1-prob_eac1_prog)) 
    ndbe_s$und_eac1_c[j] <- ((ndbe_s$und_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_diag_s)*(1-prob_eac1_prog)) 
    ndbe_s$und_eac1_d[j] <- ((ndbe_s$und_eac1_c[i]*(1-prob_death[i]))*(1-prob_eac1_prog))+((ndbe_s$und_eac1_d[i]*(1-prob_death[i]))*(1-prob_eac1_prog))
    
    # Diagnosed EAC Stage 1:
    
    ndbe_s$d_eac1_a[j] <- (((ndbe_s$und_eac1_a[i]+ndbe_s$und_eac1_b[i])*(1-prob_death[i])))*prob_eac1_diag_s                                 
    ndbe_s$d_eac1_b[j] <- (ndbe_s$d_eac1_a[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    ndbe_s$d_eac1_c[j] <- (ndbe_s$d_eac1_b[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    ndbe_s$d_eac1_d[j] <- (ndbe_s$d_eac1_c[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    ndbe_s$d_eac1_e[j] <- (ndbe_s$d_eac1_d[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    ndbe_s$d_eac1_f[j] <- (ndbe_s$d_eac1_e[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    ndbe_s$d_eac1_g[j] <- ((ndbe_s$d_eac1_f[i]*(1-prob_death[i])))+((ndbe_s$d_eac1_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 1:
    
    ndbe_s$dead_eac1[j] <- (prob_eac1_death_s*(ndbe_s$d_eac1_a[i]+ndbe_s$d_eac1_b[i]+ndbe_s$d_eac1_c[i]+ndbe_s$d_eac1_d[i]+ndbe_s$d_eac1_e[i])*(1-prob_death[i]))+ndbe_s$dead_eac1[i]
    
    # EAC Stage 2+:
    
    ndbe_s$eac2_a[j] <- ((((ndbe_s$und_eac1_a[i]+ndbe_s$und_eac1_b[i])*(1-prob_eac1_diag_s))+ndbe_s$und_eac1_c[i]+ndbe_s$und_eac1_d[i])*(1-prob_death[i])*prob_eac1_prog) 
    ndbe_s$eac2_b[j] <- ((ndbe_s$eac2_a[i]*(1-prob_death[i]))*(1-prob_eac2_y1_prog))
    ndbe_s$eac2_c[j] <- ((ndbe_s$eac2_b[i]*(1-prob_death[i]))*(1-prob_eac2_y2_prog))
    ndbe_s$eac2_d[j] <- ((ndbe_s$eac2_c[i]*(1-prob_death[i]))*(1-prob_eac2_y3_prog))
    ndbe_s$eac2_e[j] <- ((ndbe_s$eac2_d[i]*(1-prob_death[i]))*(1-prob_eac2_y4_prog))
    ndbe_s$eac2_f[j] <- ((ndbe_s$eac2_e[i]*(1-prob_death[i]))*(1-prob_eac2_y5_prog))
    ndbe_s$eac2_g[j] <- ((ndbe_s$eac2_f[i]*(1-prob_death[i])))+((ndbe_s$eac2_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 2+:
    
    ndbe_s$dead_eac2[j] <- (prob_eac2_y1_prog*ndbe_s$eac2_a[i]+prob_eac2_y2_prog*ndbe_s$eac2_b[i]+prob_eac2_y3_prog*ndbe_s$eac2_c[i]+prob_eac2_y4_prog*ndbe_s$eac2_d[i]+prob_eac2_y5_prog*ndbe_s$eac2_e[i])*(1-prob_death[i])+ndbe_s$dead_eac2[i]
    
    # Cumulative death from other causes:
    
    ndbe_s$dead_other[j] <- (prob_death[i]*
                               (ndbe_s$ndbe[i]+ndbe_s$lgdbe[i]+ndbe_s$hgdbe[i]+
                                  ndbe_s$und_eac1_a[i]+ndbe_s$und_eac1_b[i]+ndbe_s$und_eac1_c[i]+ndbe_s$und_eac1_d[i]+ndbe_s$d_eac1_a[i]+ndbe_s$d_eac1_b[i]+ndbe_s$d_eac1_c[i]+ndbe_s$d_eac1_d[i]+ndbe_s$d_eac1_e[i]+ndbe_s$d_eac1_f[i]+ndbe_s$d_eac1_g[i]+
                                  ndbe_s$eac2_a[i]+ndbe_s$eac2_b[i]+ndbe_s$eac2_c[i]+ndbe_s$eac2_d[i]+ndbe_s$eac2_e[i]+ndbe_s$eac2_f[i]+ndbe_s$eac2_g[i])
    )+ndbe_s$dead_other[i]
    
  }
  
  # Next, simulate what happens to people with NDBE at baseline who test
  # negative, which is equivalent to what would occur in the absence of
  # screening.
  
  # First, start by resetting select parameters
  
  ndbe0 <- ns[1,2] * (1-sensitivity_ndbe)
  
  # Second, create holder for for-loop
  
  ndbe_ns <- data.frame(
    gerd = c(rep(0,13))
    ,ndbe = c(ndbe0,rep(NA,12))
    ,lgdbe = c(lgdbe0,rep(NA,12))
    ,hgdbe = c(hgdbe0,rep(NA,12))
    ,und_eac1_a = c(eac1_0,rep(NA,12))
    ,und_eac1_b = c(0,rep(NA,12))
    ,und_eac1_c = c(0,rep(NA,12))
    ,und_eac1_d = c(0,rep(NA,12))
    ,d_eac1_a = c(0,rep(NA,12))
    ,d_eac1_b = c(0,rep(NA,12))
    ,d_eac1_c = c(0,rep(NA,12))
    ,d_eac1_d = c(0,rep(NA,12))
    ,d_eac1_e = c(0,rep(NA,12))
    ,d_eac1_f = c(0,rep(NA,12))
    ,d_eac1_g = c(0,rep(NA,12))
    ,dead_eac1 = c(0,rep(NA,12))
    ,eac2_a = c(eac2_0,rep(NA,12))
    ,eac2_b = c(0,rep(NA,12))
    ,eac2_c = c(0,rep(NA,12))
    ,eac2_d = c(0,rep(NA,12))
    ,eac2_e = c(0,rep(NA,12))
    ,eac2_f = c(0,rep(NA,12))
    ,eac2_g = c(0,rep(NA,12))
    ,dead_eac2 = c(0,rep(NA,12))
    ,dead_other = c(0,rep(NA,12))
  )
  
  # Finally, simulate number of people in each stage of disease each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    # NDBE through HGD:
    
    ndbe_ns$ndbe[j] <- ((ndbe_ns$ndbe[i]*(1-prob_death[i]))*(1-prob_ndbe_prog))
    ndbe_ns$lgdbe[j] <- ((ndbe_ns$lgdbe[i]*(1-prob_death[i]))*(1-prob_lgdbe_prog))+((ndbe_ns$ndbe[i]*(1-prob_death[i]))*prob_ndbe_prog)
    ndbe_ns$hgdbe[j] <- ((ndbe_ns$hgdbe[i]*(1-prob_death[i]))*(1-prob_hgdbe_prog))+((ndbe_ns$lgdbe[i]*(1-prob_death[i]))*prob_lgdbe_prog)
    
    # Undiagnosed EAC Stage 1:
    # Patients can be diagnosed in years 2 or 3, else they progress to Stage 2
    # or stay undiagnosed. Patients are not allowed to die from undiagnosed EAC
    # Stage 1.
    
    ndbe_ns$und_eac1_a[j] <- ((ndbe_ns$hgdbe[i]*(1-prob_death[i]))*prob_hgdbe_prog)
    ndbe_ns$und_eac1_b[j] <- ((ndbe_ns$und_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_prog)*(1-prob_eac1_diag))
    ndbe_ns$und_eac1_c[j] <- ((ndbe_ns$und_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_prog)*(1-prob_eac1_diag))
    ndbe_ns$und_eac1_d[j] <- ((ndbe_ns$und_eac1_c[i]*(1-prob_death[i])*(1-prob_eac1_prog)))+((ndbe_ns$und_eac1_d[i]*(1-prob_death[i])*(1-prob_eac1_prog)))
    
    # Diagnosed EAC Stage 1:
    
    ndbe_ns$d_eac1_a[j] <- (((ndbe_ns$und_eac1_a[i]+ndbe_ns$und_eac1_b[i])*(1-prob_death[i])))*prob_eac1_diag                                 
    ndbe_ns$d_eac1_b[j] <- (ndbe_ns$d_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    ndbe_ns$d_eac1_c[j] <- (ndbe_ns$d_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    ndbe_ns$d_eac1_d[j] <- (ndbe_ns$d_eac1_c[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    ndbe_ns$d_eac1_e[j] <- (ndbe_ns$d_eac1_d[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    ndbe_ns$d_eac1_f[j] <- (ndbe_ns$d_eac1_e[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    ndbe_ns$d_eac1_g[j] <- ((ndbe_ns$d_eac1_f[i]*(1-prob_death[i])))+((ndbe_ns$d_eac1_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 1:
    
    ndbe_ns$dead_eac1[j] <- (prob_eac1_death*(ndbe_ns$d_eac1_a[i]+ndbe_ns$d_eac1_b[i]+ndbe_ns$d_eac1_c[i]+ndbe_ns$d_eac1_d[i]+ndbe_ns$d_eac1_e[i])*(1-prob_death[i]))+ndbe_ns$dead_eac1[i]
    
    # EAC Stage 2+:
    
    ndbe_ns$eac2_a[j] <- ((((ndbe_ns$und_eac1_a[i]+ndbe_ns$und_eac1_b[i])*(1-prob_eac1_diag))+ndbe_ns$und_eac1_c[i]+ndbe_ns$und_eac1_d[i])*(1-prob_death[i])*prob_eac1_prog) 
    ndbe_ns$eac2_b[j] <- ((ndbe_ns$eac2_a[i]*(1-prob_death[i]))*(1-prob_eac2_y1_prog))
    ndbe_ns$eac2_c[j] <- ((ndbe_ns$eac2_b[i]*(1-prob_death[i]))*(1-prob_eac2_y2_prog))
    ndbe_ns$eac2_d[j] <- ((ndbe_ns$eac2_c[i]*(1-prob_death[i]))*(1-prob_eac2_y3_prog))
    ndbe_ns$eac2_e[j] <- ((ndbe_ns$eac2_d[i]*(1-prob_death[i]))*(1-prob_eac2_y4_prog))
    ndbe_ns$eac2_f[j] <- ((ndbe_ns$eac2_e[i]*(1-prob_death[i]))*(1-prob_eac2_y5_prog))
    ndbe_ns$eac2_g[j] <- ((ndbe_ns$eac2_f[i]*(1-prob_death[i])))+((ndbe_ns$eac2_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 2+:
    
    ndbe_ns$dead_eac2[j] <- (prob_eac2_y1_prog*ndbe_ns$eac2_a[i]+prob_eac2_y2_prog*ndbe_ns$eac2_b[i]+prob_eac2_y3_prog*ndbe_ns$eac2_c[i]+prob_eac2_y4_prog*ndbe_ns$eac2_d[i]+prob_eac2_y5_prog*ndbe_ns$eac2_e[i])*(1-prob_death[i])+ndbe_ns$dead_eac2[i]
    
    # Cumulative death from other causes:
    
    ndbe_ns$dead_other[j] <- (prob_death[i]*
                                (ndbe_ns$ndbe[i]+ndbe_ns$lgdbe[i]+ndbe_ns$hgdbe[i]+
                                   ndbe_ns$und_eac1_a[i]+ndbe_ns$und_eac1_b[i]+ndbe_ns$und_eac1_c[i]+ndbe_ns$und_eac1_d[i]+ndbe_ns$d_eac1_a[i]+ndbe_ns$d_eac1_b[i]+ndbe_ns$d_eac1_c[i]+ndbe_ns$d_eac1_d[i]+ndbe_ns$d_eac1_e[i]+ndbe_ns$d_eac1_f[i]+ndbe_ns$d_eac1_g[i]+
                                   ndbe_ns$eac2_a[i]+ndbe_ns$eac2_b[i]+ndbe_ns$eac2_c[i]+ndbe_ns$eac2_d[i]+ndbe_ns$eac2_e[i]+ndbe_ns$eac2_f[i]+ndbe_ns$eac2_g[i])
    )+ndbe_ns$dead_other[i]
    
  }
  
  # Simulate what happens to people with LGD at at baseline, first among those
  # who screen positive.
  
  # First, start by resetting select parameters
  
  lgdbe0 <- ns[1,3] * sensitivity_lgdbe
  prop_lgdbe <- 1
  prop_hgdbe <- prop_eac1 <- prop_eac2 <- 0
  hgdbe0 <- eac1_0 <- eac2_0 <- 0
  
  # Second, create holder for for-loop
  
  lgdbe_s <- data.frame(
    gerd = c(rep(0,13))
    ,ndbe = c(rep(0,13))
    ,lgdbe = c(lgdbe0,rep(NA,12))
    ,hgdbe = c(hgdbe0,rep(NA,12))
    ,und_eac1_a = c(eac1_0,rep(NA,12))
    ,und_eac1_b = c(0,rep(NA,12))
    ,und_eac1_c = c(0,rep(NA,12))
    ,und_eac1_d = c(0,rep(NA,12))
    ,d_eac1_a = c(0,rep(NA,12))
    ,d_eac1_b = c(0,rep(NA,12))
    ,d_eac1_c = c(0,rep(NA,12))
    ,d_eac1_d = c(0,rep(NA,12))
    ,d_eac1_e = c(0,rep(NA,12))
    ,d_eac1_f = c(0,rep(NA,12))
    ,d_eac1_g = c(0,rep(NA,12))
    ,dead_eac1 = c(0,rep(NA,12))
    ,eac2_a = c(eac2_0,rep(NA,12))
    ,eac2_b = c(0,rep(NA,12))
    ,eac2_c = c(0,rep(NA,12))
    ,eac2_d = c(0,rep(NA,12))
    ,eac2_e = c(0,rep(NA,12))
    ,eac2_f = c(0,rep(NA,12))
    ,eac2_g = c(0,rep(NA,12))
    ,dead_eac2 = c(0,rep(NA,12))
    ,dead_other = c(0,rep(NA,12))
  )
  
  # Finally, simulate number of people in each stage of disease each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    # LGD through HGD:
    
    lgdbe_s$lgdbe[j] <- ((lgdbe_s$lgdbe[i]*(1-prob_death[i]))*(1-(prob_lgdbe_prog*(1-efficacy))))
    lgdbe_s$hgdbe[j] <- ((lgdbe_s$hgdbe[i]*(1-prob_death[i]))*(1-(prob_hgdbe_prog*(1-efficacy))))+((lgdbe_s$lgdbe[i]*(1-prob_death[i]))*prob_lgdbe_prog*(1-efficacy))
    
    # Undiagnosed EAC Stage 1:
    # Patients can be diagnosed in years 2 or 3, else they progress to Stage 2
    # or stay undiagnosed. Patients are not allowed to die from undiagnosed EAC
    # Stage 1.
    
    lgdbe_s$und_eac1_a[j] <- ((lgdbe_s$hgdbe[i]*(1-prob_death[i]))*prob_hgdbe_prog*(1-efficacy))
    lgdbe_s$und_eac1_b[j] <- ((lgdbe_s$und_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_diag_s)*(1-prob_eac1_prog))
    lgdbe_s$und_eac1_c[j] <- ((lgdbe_s$und_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_diag_s)*(1-prob_eac1_prog)) 
    lgdbe_s$und_eac1_d[j] <- ((lgdbe_s$und_eac1_c[i]*(1-prob_death[i])*(1-prob_eac1_prog)))+((lgdbe_s$und_eac1_d[i]*(1-prob_death[i])*(1-prob_eac1_prog))) 
    
    # Diagnosed EAC Stage 1:
    
    lgdbe_s$d_eac1_a[j] <- (((lgdbe_s$und_eac1_a[i]+lgdbe_s$und_eac1_b[i])*(1-prob_death[i])))*prob_eac1_diag_s                                 
    lgdbe_s$d_eac1_b[j] <- (lgdbe_s$d_eac1_a[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    lgdbe_s$d_eac1_c[j] <- (lgdbe_s$d_eac1_b[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    lgdbe_s$d_eac1_d[j] <- (lgdbe_s$d_eac1_c[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    lgdbe_s$d_eac1_e[j] <- (lgdbe_s$d_eac1_d[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    lgdbe_s$d_eac1_f[j] <- (lgdbe_s$d_eac1_e[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    lgdbe_s$d_eac1_g[j] <- ((lgdbe_s$d_eac1_f[i]*(1-prob_death[i])))+((lgdbe_s$d_eac1_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 1:
    
    lgdbe_s$dead_eac1[j] <- (prob_eac1_death_s*(lgdbe_s$d_eac1_a[i]+lgdbe_s$d_eac1_b[i]+lgdbe_s$d_eac1_c[i]+lgdbe_s$d_eac1_d[i]+lgdbe_s$d_eac1_e[i])*(1-prob_death[i]))+lgdbe_s$dead_eac1[i]
    
    # EAC Stage 2+:
    
    lgdbe_s$eac2_a[j] <- ((((lgdbe_s$und_eac1_a[i]+lgdbe_s$und_eac1_b[i])*(1-prob_eac1_diag_s))+lgdbe_s$und_eac1_c[i]+lgdbe_s$und_eac1_d[i])*(1-prob_death[i])*prob_eac1_prog) 
    lgdbe_s$eac2_b[j] <- ((lgdbe_s$eac2_a[i]*(1-prob_death[i]))*(1-prob_eac2_y1_prog))
    lgdbe_s$eac2_c[j] <- ((lgdbe_s$eac2_b[i]*(1-prob_death[i]))*(1-prob_eac2_y2_prog))
    lgdbe_s$eac2_d[j] <- ((lgdbe_s$eac2_c[i]*(1-prob_death[i]))*(1-prob_eac2_y3_prog))
    lgdbe_s$eac2_e[j] <- ((lgdbe_s$eac2_d[i]*(1-prob_death[i]))*(1-prob_eac2_y4_prog))
    lgdbe_s$eac2_f[j] <- ((lgdbe_s$eac2_e[i]*(1-prob_death[i]))*(1-prob_eac2_y5_prog))
    lgdbe_s$eac2_g[j] <- ((lgdbe_s$eac2_f[i]*(1-prob_death[i])))+((lgdbe_s$eac2_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 2+:
    
    lgdbe_s$dead_eac2[j] <- (prob_eac2_y1_prog*lgdbe_s$eac2_a[i]+prob_eac2_y2_prog*lgdbe_s$eac2_b[i]+prob_eac2_y3_prog*lgdbe_s$eac2_c[i]+prob_eac2_y4_prog*lgdbe_s$eac2_d[i]+prob_eac2_y5_prog*lgdbe_s$eac2_e[i])*(1-prob_death[i])+lgdbe_s$dead_eac2[i]
    
    # Cumulative death from other causes:
    
    lgdbe_s$dead_other[j] <- (prob_death[i]*
                                (lgdbe_s$lgdbe[i]+lgdbe_s$hgdbe[i]+
                                   lgdbe_s$und_eac1_a[i]+lgdbe_s$und_eac1_b[i]+lgdbe_s$und_eac1_c[i]+lgdbe_s$und_eac1_d[i]+lgdbe_s$d_eac1_a[i]+lgdbe_s$d_eac1_b[i]+lgdbe_s$d_eac1_c[i]+lgdbe_s$d_eac1_d[i]+lgdbe_s$d_eac1_e[i]+lgdbe_s$d_eac1_f[i]+lgdbe_s$d_eac1_g[i]+
                                   lgdbe_s$eac2_a[i]+lgdbe_s$eac2_b[i]+lgdbe_s$eac2_c[i]+lgdbe_s$eac2_d[i]+lgdbe_s$eac2_e[i]+lgdbe_s$eac2_f[i]+lgdbe_s$eac2_g[i])
    )+lgdbe_s$dead_other[i]
    
  }

  # Next, simulate what happens to people with LGD at baseline who test
  # negative, which is equivalent to what would occur in the absence of
  # screening.
  
  # First, start by resetting select parameters
  
  lgdbe0 <- ns[1,3] * (1 - sensitivity_lgdbe)
  
  # Second, create holder for for-loop
  
  lgdbe_ns <- data.frame(
    gerd = c(rep(0,13))
    ,ndbe = c(rep(0,13))
    ,lgdbe = c(lgdbe0,rep(NA,12))
    ,hgdbe = c(hgdbe0,rep(NA,12))
    ,und_eac1_a = c(eac1_0,rep(NA,12))
    ,und_eac1_b = c(0,rep(NA,12))
    ,und_eac1_c = c(0,rep(NA,12))
    ,und_eac1_d = c(0,rep(NA,12))
    ,d_eac1_a = c(0,rep(NA,12))
    ,d_eac1_b = c(0,rep(NA,12))
    ,d_eac1_c = c(0,rep(NA,12))
    ,d_eac1_d = c(0,rep(NA,12))
    ,d_eac1_e = c(0,rep(NA,12))
    ,d_eac1_f = c(0,rep(NA,12))
    ,d_eac1_g = c(0,rep(NA,12))
    ,dead_eac1 = c(0,rep(NA,12))
    ,eac2_a = c(eac2_0,rep(NA,12))
    ,eac2_b = c(0,rep(NA,12))
    ,eac2_c = c(0,rep(NA,12))
    ,eac2_d = c(0,rep(NA,12))
    ,eac2_e = c(0,rep(NA,12))
    ,eac2_f = c(0,rep(NA,12))
    ,eac2_g = c(0,rep(NA,12))
    ,dead_eac2 = c(0,rep(NA,12))
    ,dead_other = c(0,rep(NA,12))
  )
  
  # Finally, simulate number of people in each stage of disease each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    # LGD through HGD:
    
    lgdbe_ns$lgdbe[j] <- ((lgdbe_ns$lgdbe[i]*(1-prob_death[i]))*(1-prob_lgdbe_prog))
    lgdbe_ns$hgdbe[j] <- ((lgdbe_ns$hgdbe[i]*(1-prob_death[i]))*(1-prob_hgdbe_prog))+((lgdbe_ns$lgdbe[i]*(1-prob_death[i]))*prob_lgdbe_prog)
    
    # Undiagnosed EAC Stage 1:
    # Patients can be diagnosed in years 2 or 3, else they progress to Stage 2
    # or stay undiagnosed. Patients are not allowed to die from undiagnosed EAC
    # Stage 1.
    
    lgdbe_ns$und_eac1_a[j] <- ((lgdbe_ns$hgdbe[i]*(1-prob_death[i]))*prob_hgdbe_prog)
    lgdbe_ns$und_eac1_b[j] <- ((lgdbe_ns$und_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_prog)*(1-prob_eac1_diag))
    lgdbe_ns$und_eac1_c[j] <- ((lgdbe_ns$und_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_prog)*(1-prob_eac1_diag))
    lgdbe_ns$und_eac1_d[j] <- ((lgdbe_ns$und_eac1_c[i]*(1-prob_death[i])*(1-prob_eac1_prog)))+((lgdbe_ns$und_eac1_d[i]*(1-prob_death[i])*(1-prob_eac1_prog)))
    
    # Diagnosed EAC Stage 1:
    
    lgdbe_ns$d_eac1_a[j] <- (((lgdbe_ns$und_eac1_a[i]+lgdbe_ns$und_eac1_b[i])*(1-prob_death[i])))*prob_eac1_diag                                 
    lgdbe_ns$d_eac1_b[j] <- (lgdbe_ns$d_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    lgdbe_ns$d_eac1_c[j] <- (lgdbe_ns$d_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    lgdbe_ns$d_eac1_d[j] <- (lgdbe_ns$d_eac1_c[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    lgdbe_ns$d_eac1_e[j] <- (lgdbe_ns$d_eac1_d[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    lgdbe_ns$d_eac1_f[j] <- (lgdbe_ns$d_eac1_e[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    lgdbe_ns$d_eac1_g[j] <- ((lgdbe_ns$d_eac1_f[i]*(1-prob_death[i])))+((lgdbe_ns$d_eac1_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 1:
    
    lgdbe_ns$dead_eac1[j] <- (prob_eac1_death*(lgdbe_ns$d_eac1_a[i]+lgdbe_ns$d_eac1_b[i]+lgdbe_ns$d_eac1_c[i]+lgdbe_ns$d_eac1_d[i]+lgdbe_ns$d_eac1_e[i])*(1-prob_death[i]))+lgdbe_ns$dead_eac1[i]
    
    # EAC Stage 2+:
    
    lgdbe_ns$eac2_a[j] <- ((((lgdbe_ns$und_eac1_a[i]+lgdbe_ns$und_eac1_b[i])*(1-prob_eac1_diag))+lgdbe_ns$und_eac1_c[i]+lgdbe_ns$und_eac1_d[i])*(1-prob_death[i])*prob_eac1_prog)
    lgdbe_ns$eac2_b[j] <- ((lgdbe_ns$eac2_a[i]*(1-prob_death[i]))*(1-prob_eac2_y1_prog))
    lgdbe_ns$eac2_c[j] <- ((lgdbe_ns$eac2_b[i]*(1-prob_death[i]))*(1-prob_eac2_y2_prog))
    lgdbe_ns$eac2_d[j] <- ((lgdbe_ns$eac2_c[i]*(1-prob_death[i]))*(1-prob_eac2_y3_prog))
    lgdbe_ns$eac2_e[j] <- ((lgdbe_ns$eac2_d[i]*(1-prob_death[i]))*(1-prob_eac2_y4_prog))
    lgdbe_ns$eac2_f[j] <- ((lgdbe_ns$eac2_e[i]*(1-prob_death[i]))*(1-prob_eac2_y5_prog))
    lgdbe_ns$eac2_g[j] <- ((lgdbe_ns$eac2_f[i]*(1-prob_death[i])))+((lgdbe_ns$eac2_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 2+:
    
    lgdbe_ns$dead_eac2[j] <- (prob_eac2_y1_prog*lgdbe_ns$eac2_a[i]+prob_eac2_y2_prog*lgdbe_ns$eac2_b[i]+prob_eac2_y3_prog*lgdbe_ns$eac2_c[i]+prob_eac2_y4_prog*lgdbe_ns$eac2_d[i]+prob_eac2_y5_prog*lgdbe_ns$eac2_e[i])*(1-prob_death[i])+lgdbe_ns$dead_eac2[i]
    
    # Cumulative death from other causes:
    
    lgdbe_ns$dead_other[j] <- (prob_death[i]*
                                 (lgdbe_ns$lgdbe[i]+lgdbe_ns$hgdbe[i]+
                                    lgdbe_ns$und_eac1_a[i]+lgdbe_ns$und_eac1_b[i]+lgdbe_ns$und_eac1_c[i]+lgdbe_ns$und_eac1_d[i]+lgdbe_ns$d_eac1_a[i]+lgdbe_ns$d_eac1_b[i]+lgdbe_ns$d_eac1_c[i]+lgdbe_ns$d_eac1_d[i]+lgdbe_ns$d_eac1_e[i]+lgdbe_ns$d_eac1_f[i]+lgdbe_ns$d_eac1_g[i]+
                                    lgdbe_ns$eac2_a[i]+lgdbe_ns$eac2_b[i]+lgdbe_ns$eac2_c[i]+lgdbe_ns$eac2_d[i]+lgdbe_ns$eac2_e[i]+lgdbe_ns$eac2_f[i]+lgdbe_ns$eac2_g[i])
    )+lgdbe_ns$dead_other[i]
    
  }
  
  # Simulate what happens to people with HGD at at baseline, first among those
  # who screen positive.
  
  # First, start by resetting select parameters
  
  hgdbe0 <- ns[1,4] * sensitivity_hgdbe
  prop_hgdbe <- 1
  prop_eac1 <- prop_eac2 <- 0
  eac1_0 <- eac2_0 <- 0
  
  # Second, create holder for for-loop
  
  hgdbe_s <- data.frame(
    gerd = c(rep(0,13))
    ,ndbe = c(rep(0,13))
    ,lgdbe = c(rep(0,13))
    ,hgdbe = c(hgdbe0,rep(NA,12))
    ,und_eac1_a = c(eac1_0,rep(NA,12))
    ,und_eac1_b = c(0,rep(NA,12))
    ,und_eac1_c = c(0,rep(NA,12))
    ,und_eac1_d = c(0,rep(NA,12))
    ,d_eac1_a = c(0,rep(NA,12))
    ,d_eac1_b = c(0,rep(NA,12))
    ,d_eac1_c = c(0,rep(NA,12))
    ,d_eac1_d = c(0,rep(NA,12))
    ,d_eac1_e = c(0,rep(NA,12))
    ,d_eac1_f = c(0,rep(NA,12))
    ,d_eac1_g = c(0,rep(NA,12))
    ,dead_eac1 = c(0,rep(NA,12))
    ,eac2_a = c(eac2_0,rep(NA,12))
    ,eac2_b = c(0,rep(NA,12))
    ,eac2_c = c(0,rep(NA,12))
    ,eac2_d = c(0,rep(NA,12))
    ,eac2_e = c(0,rep(NA,12))
    ,eac2_f = c(0,rep(NA,12))
    ,eac2_g = c(0,rep(NA,12))
    ,dead_eac2 = c(0,rep(NA,12))
    ,dead_other = c(0,rep(NA,12))
  )
  
  # Finally, simulate number of people in each stage of disease each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    # HGD:
    
    hgdbe_s$hgdbe[j] <- ((hgdbe_s$hgdbe[i]*(1-prob_death[i]))*(1-(prob_hgdbe_prog*(1-efficacy))))
    
    # Undiagnosed EAC Stage 1:
    # Patients can be diagnosed in years 2 or 3, else they progress to Stage 2
    # or stay undiagnosed. Patients are not allowed to die from undiagnosed EAC
    # Stage 1.
    
    hgdbe_s$und_eac1_a[j] <- ((hgdbe_s$hgdbe[i]*(1-prob_death[i]))*prob_hgdbe_prog*(1-efficacy))
    hgdbe_s$und_eac1_b[j] <- ((hgdbe_s$und_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_diag_s)*(1-prob_eac1_prog))
    hgdbe_s$und_eac1_c[j] <- ((hgdbe_s$und_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_diag_s)*(1-prob_eac1_prog)) 
    hgdbe_s$und_eac1_d[j] <- ((hgdbe_s$und_eac1_c[i]*(1-prob_death[i]))*(1-prob_eac1_prog))+((hgdbe_s$und_eac1_d[i]*(1-prob_death[i]))*(1-prob_eac1_prog)) 
    
    # Diagnosed EAC Stage 1:
    
    hgdbe_s$d_eac1_a[j] <- (((hgdbe_s$und_eac1_a[i]+hgdbe_s$und_eac1_b[i])*(1-prob_death[i])))*prob_eac1_diag_s                                 
    hgdbe_s$d_eac1_b[j] <- (hgdbe_s$d_eac1_a[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    hgdbe_s$d_eac1_c[j] <- (hgdbe_s$d_eac1_b[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    hgdbe_s$d_eac1_d[j] <- (hgdbe_s$d_eac1_c[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    hgdbe_s$d_eac1_e[j] <- (hgdbe_s$d_eac1_d[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    hgdbe_s$d_eac1_f[j] <- (hgdbe_s$d_eac1_e[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    hgdbe_s$d_eac1_g[j] <- ((hgdbe_s$d_eac1_f[i]*(1-prob_death[i])))+((hgdbe_s$d_eac1_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 1:
    
    hgdbe_s$dead_eac1[j] <- (prob_eac1_death_s*(hgdbe_s$d_eac1_a[i]+hgdbe_s$d_eac1_b[i]+hgdbe_s$d_eac1_c[i]+hgdbe_s$d_eac1_d[i]+hgdbe_s$d_eac1_e[i])*(1-prob_death[i]))+hgdbe_s$dead_eac1[i]
    
    # EAC Stage 2+:
    
    hgdbe_s$eac2_a[j] <- ((((hgdbe_s$und_eac1_a[i]+hgdbe_s$und_eac1_b[i])*(1-prob_eac1_diag_s))+hgdbe_s$und_eac1_c[i]+hgdbe_s$und_eac1_d[i])*(1-prob_death[i])*prob_eac1_prog)
    hgdbe_s$eac2_b[j] <- ((hgdbe_s$eac2_a[i]*(1-prob_death[i]))*(1-prob_eac2_y1_prog))
    hgdbe_s$eac2_c[j] <- ((hgdbe_s$eac2_b[i]*(1-prob_death[i]))*(1-prob_eac2_y2_prog))
    hgdbe_s$eac2_d[j] <- ((hgdbe_s$eac2_c[i]*(1-prob_death[i]))*(1-prob_eac2_y3_prog))
    hgdbe_s$eac2_e[j] <- ((hgdbe_s$eac2_d[i]*(1-prob_death[i]))*(1-prob_eac2_y4_prog))
    hgdbe_s$eac2_f[j] <- ((hgdbe_s$eac2_e[i]*(1-prob_death[i]))*(1-prob_eac2_y5_prog))
    hgdbe_s$eac2_g[j] <- ((hgdbe_s$eac2_f[i]*(1-prob_death[i])))+((hgdbe_s$eac2_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 2+:
    
    hgdbe_s$dead_eac2[j] <- (prob_eac2_y1_prog*hgdbe_s$eac2_a[i]+prob_eac2_y2_prog*hgdbe_s$eac2_b[i]+prob_eac2_y3_prog*hgdbe_s$eac2_c[i]+prob_eac2_y4_prog*hgdbe_s$eac2_d[i]+prob_eac2_y5_prog*hgdbe_s$eac2_e[i])*(1-prob_death[i])+hgdbe_s$dead_eac2[i]
    
    # Cumulative death from other causes:
    
    hgdbe_s$dead_other[j] <- (prob_death[i]*
                                (hgdbe_s$hgdbe[i]+
                                   hgdbe_s$und_eac1_a[i]+hgdbe_s$und_eac1_b[i]+hgdbe_s$und_eac1_c[i]+hgdbe_s$und_eac1_d[i]+hgdbe_s$d_eac1_a[i]+hgdbe_s$d_eac1_b[i]+hgdbe_s$d_eac1_c[i]+hgdbe_s$d_eac1_d[i]+hgdbe_s$d_eac1_e[i]+hgdbe_s$d_eac1_f[i]+hgdbe_s$d_eac1_g[i]+
                                   hgdbe_s$eac2_a[i]+hgdbe_s$eac2_b[i]+hgdbe_s$eac2_c[i]+hgdbe_s$eac2_d[i]+hgdbe_s$eac2_e[i]+hgdbe_s$eac2_f[i]+hgdbe_s$eac2_g[i])
    )+hgdbe_s$dead_other[i]
    
  }
  
  # Next, simulate what happens to people with HGD at baseline who test
  # negative, which is equivalent to what would occur in the absence of
  # screening.
  
  # First, start by resetting select parameters
  
  hgdbe0 <- ns[1,4] * (1 - sensitivity_hgdbe)
  
  # Second, create holder for for-loop
  
  hgdbe_ns <- data.frame(
    gerd = c(rep(0,13))
    ,ndbe = c(rep(0,13))
    ,lgdbe = c(rep(0,13))
    ,hgdbe = c(hgdbe0,rep(NA,12))
    ,und_eac1_a = c(eac1_0,rep(NA,12))
    ,und_eac1_b = c(0,rep(NA,12))
    ,und_eac1_c = c(0,rep(NA,12))
    ,und_eac1_d = c(0,rep(NA,12))
    ,d_eac1_a = c(0,rep(NA,12))
    ,d_eac1_b = c(0,rep(NA,12))
    ,d_eac1_c = c(0,rep(NA,12))
    ,d_eac1_d = c(0,rep(NA,12))
    ,d_eac1_e = c(0,rep(NA,12))
    ,d_eac1_f = c(0,rep(NA,12))
    ,d_eac1_g = c(0,rep(NA,12))
    ,dead_eac1 = c(0,rep(NA,12))
    ,eac2_a = c(eac2_0,rep(NA,12))
    ,eac2_b = c(0,rep(NA,12))
    ,eac2_c = c(0,rep(NA,12))
    ,eac2_d = c(0,rep(NA,12))
    ,eac2_e = c(0,rep(NA,12))
    ,eac2_f = c(0,rep(NA,12))
    ,eac2_g = c(0,rep(NA,12))
    ,dead_eac2 = c(0,rep(NA,12))
    ,dead_other = c(0,rep(NA,12))
  )
  
  # Finally, simulate number of people in each stage of disease each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    # HGD:
    
    hgdbe_ns$hgdbe[j] <- ((hgdbe_ns$hgdbe[i]*(1-prob_death[i]))*(1-prob_hgdbe_prog))
    
    # Undiagnosed EAC Stage 1:
    # Patients can be diagnosed in years 2 or 3, else they progress to Stage 2
    # or stay undiagnosed. Patients are not allowed to die from undiagnosed EAC
    # Stage 1.
    
    hgdbe_ns$und_eac1_a[j] <- ((hgdbe_ns$hgdbe[i]*(1-prob_death[i]))*prob_hgdbe_prog)
    hgdbe_ns$und_eac1_b[j] <- ((hgdbe_ns$und_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_prog)*(1-prob_eac1_diag))
    hgdbe_ns$und_eac1_c[j] <- ((hgdbe_ns$und_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_prog)*(1-prob_eac1_diag))
    hgdbe_ns$und_eac1_d[j] <- ((hgdbe_ns$und_eac1_c[i]*(1-prob_death[i])*(1-prob_eac1_prog)))+((hgdbe_ns$und_eac1_d[i]*(1-prob_death[i])*(1-prob_eac1_prog)))
    
    # Diagnosed EAC Stage 1:
    
    hgdbe_ns$d_eac1_a[j] <- (((hgdbe_ns$und_eac1_a[i]+hgdbe_ns$und_eac1_b[i])*(1-prob_death[i])))*prob_eac1_diag                                 
    hgdbe_ns$d_eac1_b[j] <- (hgdbe_ns$d_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    hgdbe_ns$d_eac1_c[j] <- (hgdbe_ns$d_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    hgdbe_ns$d_eac1_d[j] <- (hgdbe_ns$d_eac1_c[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    hgdbe_ns$d_eac1_e[j] <- (hgdbe_ns$d_eac1_d[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    hgdbe_ns$d_eac1_f[j] <- (hgdbe_ns$d_eac1_e[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    hgdbe_ns$d_eac1_g[j] <- ((hgdbe_ns$d_eac1_f[i]*(1-prob_death[i])))+((hgdbe_ns$d_eac1_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 1:
    
    hgdbe_ns$dead_eac1[j] <- (prob_eac1_death*(hgdbe_ns$d_eac1_a[i]+hgdbe_ns$d_eac1_b[i]+hgdbe_ns$d_eac1_c[i]+hgdbe_ns$d_eac1_d[i]+hgdbe_ns$d_eac1_e[i])*(1-prob_death[i]))+hgdbe_ns$dead_eac1[i]
    
    # EAC Stage 2+:
    
    hgdbe_ns$eac2_a[j] <- ((((hgdbe_ns$und_eac1_a[i]+hgdbe_ns$und_eac1_b[i])*(1-prob_eac1_diag))+hgdbe_ns$und_eac1_c[i]+hgdbe_ns$und_eac1_d[i])*(1-prob_death[i])*prob_eac1_prog)
    hgdbe_ns$eac2_b[j] <- ((hgdbe_ns$eac2_a[i]*(1-prob_death[i]))*(1-prob_eac2_y1_prog))
    hgdbe_ns$eac2_c[j] <- ((hgdbe_ns$eac2_b[i]*(1-prob_death[i]))*(1-prob_eac2_y2_prog))
    hgdbe_ns$eac2_d[j] <- ((hgdbe_ns$eac2_c[i]*(1-prob_death[i]))*(1-prob_eac2_y3_prog))
    hgdbe_ns$eac2_e[j] <- ((hgdbe_ns$eac2_d[i]*(1-prob_death[i]))*(1-prob_eac2_y4_prog))
    hgdbe_ns$eac2_f[j] <- ((hgdbe_ns$eac2_e[i]*(1-prob_death[i]))*(1-prob_eac2_y5_prog))
    hgdbe_ns$eac2_g[j] <- ((hgdbe_ns$eac2_f[i]*(1-prob_death[i])))+((hgdbe_ns$eac2_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 2+:
    
    hgdbe_ns$dead_eac2[j] <- (prob_eac2_y1_prog*hgdbe_ns$eac2_a[i]+prob_eac2_y2_prog*hgdbe_ns$eac2_b[i]+prob_eac2_y3_prog*hgdbe_ns$eac2_c[i]+prob_eac2_y4_prog*hgdbe_ns$eac2_d[i]+prob_eac2_y5_prog*hgdbe_ns$eac2_e[i])*(1-prob_death[i])+hgdbe_ns$dead_eac2[i]
    
    # Cumulative death from other causes:
    
    hgdbe_ns$dead_other[j] <- (prob_death[i]*
                                 (hgdbe_ns$hgdbe[i]+
                                    hgdbe_ns$und_eac1_a[i]+hgdbe_ns$und_eac1_b[i]+hgdbe_ns$und_eac1_c[i]+hgdbe_ns$und_eac1_d[i]+hgdbe_ns$d_eac1_a[i]+hgdbe_ns$d_eac1_b[i]+hgdbe_ns$d_eac1_c[i]+hgdbe_ns$d_eac1_d[i]+hgdbe_ns$d_eac1_e[i]+hgdbe_ns$d_eac1_f[i]+hgdbe_ns$d_eac1_g[i]+
                                    hgdbe_ns$eac2_a[i]+hgdbe_ns$eac2_b[i]+hgdbe_ns$eac2_c[i]+hgdbe_ns$eac2_d[i]+hgdbe_ns$eac2_e[i]+hgdbe_ns$eac2_f[i]+hgdbe_ns$eac2_g[i])
    )+hgdbe_ns$dead_other[i]
    
  }
  
  # Simulate what happens to people with EAC 1 at at baseline, first among those
  # who screen positive.
  
  # First, start by resetting select parameters
  
  eac1_0 <- ns[1,5] * sensitivity_eac1
  prop_eac1 <- 1
  prop_eac2 <- 0
  eac2_0 <- 0
  
  # Second, create holder for for-loop
  
  eac1_s <- data.frame(
    gerd = c(rep(0,13))
    ,ndbe = c(rep(0,13))
    ,lgdbe = c(rep(0,13))
    ,hgdbe = c(rep(0,13))
    ,und_eac1_a = c(eac1_0,rep(NA,12))
    ,und_eac1_b = c(0,rep(NA,12))
    ,und_eac1_c = c(0,rep(NA,12))
    ,und_eac1_d = c(0,rep(NA,12))
    ,d_eac1_a = c(0,rep(NA,12))
    ,d_eac1_b = c(0,rep(NA,12))
    ,d_eac1_c = c(0,rep(NA,12))
    ,d_eac1_d = c(0,rep(NA,12))
    ,d_eac1_e = c(0,rep(NA,12))
    ,d_eac1_f = c(0,rep(NA,12))
    ,d_eac1_g = c(0,rep(NA,12))
    ,dead_eac1 = c(0,rep(NA,12))
    ,eac2_a = c(eac2_0,rep(NA,12))
    ,eac2_b = c(0,rep(NA,12))
    ,eac2_c = c(0,rep(NA,12))
    ,eac2_d = c(0,rep(NA,12))
    ,eac2_e = c(0,rep(NA,12))
    ,eac2_f = c(0,rep(NA,12))
    ,eac2_g = c(0,rep(NA,12))
    ,dead_eac2 = c(0,rep(NA,12))
    ,dead_other = c(0,rep(NA,12))
  )
  
  # Finally, simulate number of people in each stage of disease each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    # Undiagnosed EAC Stage 1:
    # Patients can be diagnosed in years 2 or 3, else they progress to Stage 2
    # or stay undiagnosed. Patients are not allowed to die from undiagnosed EAC
    # Stage 1.
    
    eac1_s$und_eac1_a[j] <- ((eac1_s$hgdbe[i]*(1-prob_death[i]))*prob_hgdbe_prog*(1-efficacy))
    eac1_s$und_eac1_b[j] <- ((eac1_s$und_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_diag_s)*(1-prob_eac1_prog))
    eac1_s$und_eac1_c[j] <- ((eac1_s$und_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_diag_s)*(1-prob_eac1_prog)) 
    eac1_s$und_eac1_d[j] <- ((eac1_s$und_eac1_c[i]*(1-prob_death[i]))*(1-prob_eac1_prog))+((eac1_s$und_eac1_d[i]*(1-prob_death[i]))*(1-prob_eac1_prog)) 
    
    # Diagnosed EAC Stage 1:
    
    eac1_s$d_eac1_a[j] <- (((eac1_s$und_eac1_a[i]+eac1_s$und_eac1_b[i])*(1-prob_death[i])))*prob_eac1_diag_s                                 
    eac1_s$d_eac1_b[j] <- (eac1_s$d_eac1_a[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    eac1_s$d_eac1_c[j] <- (eac1_s$d_eac1_b[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    eac1_s$d_eac1_d[j] <- (eac1_s$d_eac1_c[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    eac1_s$d_eac1_e[j] <- (eac1_s$d_eac1_d[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    eac1_s$d_eac1_f[j] <- (eac1_s$d_eac1_e[i]*(1-prob_death[i]))*(1-(prob_eac1_death_s))
    eac1_s$d_eac1_g[j] <- ((eac1_s$d_eac1_f[i]*(1-prob_death[i])))+((eac1_s$d_eac1_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 1:
    
    eac1_s$dead_eac1[j] <- (prob_eac1_death_s*(eac1_s$d_eac1_a[i]+eac1_s$d_eac1_b[i]+eac1_s$d_eac1_c[i]+eac1_s$d_eac1_d[i]+eac1_s$d_eac1_e[i])*(1-prob_death[i]))+eac1_s$dead_eac1[i]
    
    # EAC Stage 2+:
    
    eac1_s$eac2_a[j] <- ((((eac1_s$und_eac1_a[i]+eac1_s$und_eac1_b[i])*(1-prob_eac1_diag_s))+eac1_s$und_eac1_c[i]+eac1_s$und_eac1_d[i])*(1-prob_death[i])*prob_eac1_prog)
    eac1_s$eac2_b[j] <- ((eac1_s$eac2_a[i]*(1-prob_death[i]))*(1-prob_eac2_y1_prog))
    eac1_s$eac2_c[j] <- ((eac1_s$eac2_b[i]*(1-prob_death[i]))*(1-prob_eac2_y2_prog))
    eac1_s$eac2_d[j] <- ((eac1_s$eac2_c[i]*(1-prob_death[i]))*(1-prob_eac2_y3_prog))
    eac1_s$eac2_e[j] <- ((eac1_s$eac2_d[i]*(1-prob_death[i]))*(1-prob_eac2_y4_prog))
    eac1_s$eac2_f[j] <- ((eac1_s$eac2_e[i]*(1-prob_death[i]))*(1-prob_eac2_y5_prog))
    eac1_s$eac2_g[j] <- ((eac1_s$eac2_f[i]*(1-prob_death[i])))+((eac1_s$eac2_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 2+:
    
    eac1_s$dead_eac2[j] <- (prob_eac2_y1_prog*eac1_s$eac2_a[i]+prob_eac2_y2_prog*eac1_s$eac2_b[i]+prob_eac2_y3_prog*eac1_s$eac2_c[i]+prob_eac2_y4_prog*eac1_s$eac2_d[i]+prob_eac2_y5_prog*eac1_s$eac2_e[i])*(1-prob_death[i])+eac1_s$dead_eac2[i]
    
    # Cumulative death from other causes:
    
    eac1_s$dead_other[j] <- (prob_death[i]*
                               (eac1_s$und_eac1_a[i]+eac1_s$und_eac1_b[i]+eac1_s$und_eac1_c[i]+eac1_s$und_eac1_d[i]+eac1_s$d_eac1_a[i]+eac1_s$d_eac1_b[i]+eac1_s$d_eac1_c[i]+eac1_s$d_eac1_d[i]+eac1_s$d_eac1_e[i]+eac1_s$d_eac1_f[i]+eac1_s$d_eac1_g[i]+
                                  eac1_s$eac2_a[i]+eac1_s$eac2_b[i]+eac1_s$eac2_c[i]+eac1_s$eac2_d[i]+eac1_s$eac2_e[i]+eac1_s$eac2_f[i]+eac1_s$eac2_g[i])
    )+eac1_s$dead_other[i]
    
  }
  
  # Next, simulate what happens to people with EAC 1 at baseline who test
  # negative, which is equivalent to what would occur in the absence of
  # screening.
  
  # First, start by resetting select parameters
  
  eac1_0 <- ns[1,5] * (1 - sensitivity_eac1)
  
  # Second, create holder for for-loop
  
  eac1_ns <- data.frame(
    gerd = c(rep(0,13))
    ,ndbe = c(rep(0,13))
    ,lgdbe = c(rep(0,13))
    ,hgdbe = c(rep(0,13))
    ,und_eac1_a = c(eac1_0,rep(NA,12))
    ,und_eac1_b = c(0,rep(NA,12))
    ,und_eac1_c = c(0,rep(NA,12))
    ,und_eac1_d = c(0,rep(NA,12))
    ,d_eac1_a = c(0,rep(NA,12))
    ,d_eac1_b = c(0,rep(NA,12))
    ,d_eac1_c = c(0,rep(NA,12))
    ,d_eac1_d = c(0,rep(NA,12))
    ,d_eac1_e = c(0,rep(NA,12))
    ,d_eac1_f = c(0,rep(NA,12))
    ,d_eac1_g = c(0,rep(NA,12))
    ,dead_eac1 = c(0,rep(NA,12))
    ,eac2_a = c(eac2_0,rep(NA,12))
    ,eac2_b = c(0,rep(NA,12))
    ,eac2_c = c(0,rep(NA,12))
    ,eac2_d = c(0,rep(NA,12))
    ,eac2_e = c(0,rep(NA,12))
    ,eac2_f = c(0,rep(NA,12))
    ,eac2_g = c(0,rep(NA,12))
    ,dead_eac2 = c(0,rep(NA,12))
    ,dead_other = c(0,rep(NA,12))
  )
  
  # Finally, simulate number of people in each stage of disease each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    # Undiagnosed EAC Stage 1:
    # Patients can be diagnosed in years 2 or 3, else they progress to Stage 2
    # or stay undiagnosed. Patients are not allowed to die from undiagnosed EAC
    # Stage 1.
    
    eac1_ns$und_eac1_a[j] <- ((eac1_ns$hgdbe[i]*(1-prob_death[i]))*prob_hgdbe_prog)
    eac1_ns$und_eac1_b[j] <- ((eac1_ns$und_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_prog)*(1-prob_eac1_diag))
    eac1_ns$und_eac1_c[j] <- ((eac1_ns$und_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_prog)*(1-prob_eac1_diag))
    eac1_ns$und_eac1_d[j] <- ((eac1_ns$und_eac1_c[i]*(1-prob_death[i])*(1-prob_eac1_prog)))+((eac1_ns$und_eac1_d[i]*(1-prob_death[i])*(1-prob_eac1_prog)))
    
    # Diagnosed EAC Stage 1:
    
    eac1_ns$d_eac1_a[j] <- (((eac1_ns$und_eac1_a[i]+eac1_ns$und_eac1_b[i])*(1-prob_death[i])))*prob_eac1_diag                                 
    eac1_ns$d_eac1_b[j] <- (eac1_ns$d_eac1_a[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    eac1_ns$d_eac1_c[j] <- (eac1_ns$d_eac1_b[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    eac1_ns$d_eac1_d[j] <- (eac1_ns$d_eac1_c[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    eac1_ns$d_eac1_e[j] <- (eac1_ns$d_eac1_d[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    eac1_ns$d_eac1_f[j] <- (eac1_ns$d_eac1_e[i]*(1-prob_death[i]))*(1-prob_eac1_death)
    eac1_ns$d_eac1_g[j] <- ((eac1_ns$d_eac1_f[i]*(1-prob_death[i])))+((eac1_ns$d_eac1_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 1:
    
    eac1_ns$dead_eac1[j] <- (prob_eac1_death*(eac1_ns$d_eac1_a[i]+eac1_ns$d_eac1_b[i]+eac1_ns$d_eac1_c[i]+eac1_ns$d_eac1_d[i]+eac1_ns$d_eac1_e[i])*(1-prob_death[i]))+eac1_ns$dead_eac1[i]
    
    # EAC Stage 2+:
    
    eac1_ns$eac2_a[j] <- ((((eac1_ns$und_eac1_a[i]+eac1_ns$und_eac1_b[i])*(1-prob_eac1_diag))+eac1_ns$und_eac1_c[i]+eac1_ns$und_eac1_d[i])*(1-prob_death[i])*prob_eac1_prog)
    eac1_ns$eac2_b[j] <- ((eac1_ns$eac2_a[i]*(1-prob_death[i]))*(1-prob_eac2_y1_prog))
    eac1_ns$eac2_c[j] <- ((eac1_ns$eac2_b[i]*(1-prob_death[i]))*(1-prob_eac2_y2_prog))
    eac1_ns$eac2_d[j] <- ((eac1_ns$eac2_c[i]*(1-prob_death[i]))*(1-prob_eac2_y3_prog))
    eac1_ns$eac2_e[j] <- ((eac1_ns$eac2_d[i]*(1-prob_death[i]))*(1-prob_eac2_y4_prog))
    eac1_ns$eac2_f[j] <- ((eac1_ns$eac2_e[i]*(1-prob_death[i]))*(1-prob_eac2_y5_prog))
    eac1_ns$eac2_g[j] <- ((eac1_ns$eac2_f[i]*(1-prob_death[i])))+((eac1_ns$eac2_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 2+:
    
    eac1_ns$dead_eac2[j] <- (prob_eac2_y1_prog*eac1_ns$eac2_a[i]+prob_eac2_y2_prog*eac1_ns$eac2_b[i]+prob_eac2_y3_prog*eac1_ns$eac2_c[i]+prob_eac2_y4_prog*eac1_ns$eac2_d[i]+prob_eac2_y5_prog*eac1_ns$eac2_e[i])*(1-prob_death[i])+eac1_ns$dead_eac2[i]
    
    # Cumulative death from other causes:
    
    eac1_ns$dead_other[j] <- (prob_death[i]*
                                (eac1_ns$und_eac1_a[i]+eac1_ns$und_eac1_b[i]+eac1_ns$und_eac1_c[i]+eac1_ns$und_eac1_d[i]+eac1_ns$d_eac1_a[i]+eac1_ns$d_eac1_b[i]+eac1_ns$d_eac1_c[i]+eac1_ns$d_eac1_d[i]+eac1_ns$d_eac1_e[i]+eac1_ns$d_eac1_f[i]+eac1_ns$d_eac1_g[i]+
                                   eac1_ns$eac2_a[i]+eac1_ns$eac2_b[i]+eac1_ns$eac2_c[i]+eac1_ns$eac2_d[i]+eac1_ns$eac2_e[i]+eac1_ns$eac2_f[i]+eac1_ns$eac2_g[i])
    )+eac1_ns$dead_other[i]
    
  }
  
  # Simulate what happens to people with EAC 2+ at at baseline, none of whom would
  # benefit from screening. Therefore, we do not need to separately model
  # what would happen in the presence and absence of screening.
  
  # First, start by resetting select parameters
  
  eac2_0 <- ns[1,17]
  prop_eac2 <- 1
  
  # Second, create holder for for-loop
  
  eac2_ns <- data.frame(
    gerd = c(rep(0,13))
    ,ndbe = c(rep(0,13))
    ,lgdbe = c(rep(0,13))
    ,hgdbe = c(rep(0,13))
    ,und_eac1_a = c(rep(0,13))
    ,und_eac1_b = c(rep(0,13))
    ,und_eac1_c = c(rep(0,13))
    ,und_eac1_d = c(rep(0,13))
    ,d_eac1_a = c(rep(0,13))
    ,d_eac1_b = c(rep(0,13))
    ,d_eac1_c = c(rep(0,13))
    ,d_eac1_d = c(rep(0,13))
    ,d_eac1_e = c(rep(0,13))
    ,d_eac1_f = c(rep(0,13))
    ,d_eac1_g = c(rep(0,13))
    ,dead_eac1 = c(rep(0,13))
    ,eac2_a = c(eac2_0,rep(NA,12))
    ,eac2_b = c(0,rep(NA,12))
    ,eac2_c = c(0,rep(NA,12))
    ,eac2_d = c(0,rep(NA,12))
    ,eac2_e = c(0,rep(NA,12))
    ,eac2_f = c(0,rep(NA,12))
    ,eac2_g = c(0,rep(NA,12))
    ,dead_eac2 = c(0,rep(NA,12))
    ,dead_other = c(0,rep(NA,12))
  )
  
  # Finally, simulate number of people in each stage of disease each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    # EAC Stage 2+:
    
    eac2_ns$eac2_a[j] <- ((((eac2_ns$und_eac1_a[i]+eac2_ns$und_eac1_b[i])*(1-prob_eac1_diag))+eac2_ns$und_eac1_c[i]+eac2_ns$und_eac1_d[i])*(1-prob_death[i])*prob_eac1_prog)
    eac2_ns$eac2_b[j] <- ((eac2_ns$eac2_a[i]*(1-prob_death[i]))*(1-prob_eac2_y1_prog))
    eac2_ns$eac2_c[j] <- ((eac2_ns$eac2_b[i]*(1-prob_death[i]))*(1-prob_eac2_y2_prog))
    eac2_ns$eac2_d[j] <- ((eac2_ns$eac2_c[i]*(1-prob_death[i]))*(1-prob_eac2_y3_prog))
    eac2_ns$eac2_e[j] <- ((eac2_ns$eac2_d[i]*(1-prob_death[i]))*(1-prob_eac2_y4_prog))
    eac2_ns$eac2_f[j] <- ((eac2_ns$eac2_e[i]*(1-prob_death[i]))*(1-prob_eac2_y5_prog))
    eac2_ns$eac2_g[j] <- ((eac2_ns$eac2_f[i]*(1-prob_death[i])))+((eac2_ns$eac2_g[i]*(1-prob_death[i])))
    
    # Cumulative death from EAC Stage 2+:
    
    eac2_ns$dead_eac2[j] <- (prob_eac2_y1_prog*eac2_ns$eac2_a[i]+prob_eac2_y2_prog*eac2_ns$eac2_b[i]+prob_eac2_y3_prog*eac2_ns$eac2_c[i]+prob_eac2_y4_prog*eac2_ns$eac2_d[i]+prob_eac2_y5_prog*eac2_ns$eac2_e[i])*(1-prob_death[i])+eac2_ns$dead_eac2[i]
    
    # Cumulative death from other causes:
    
    eac2_ns$dead_other[j] <- (prob_death[i]*
                                (eac2_ns$eac2_a[i]+eac2_ns$eac2_b[i]+eac2_ns$eac2_c[i]+eac2_ns$eac2_d[i]+eac2_ns$eac2_e[i]+eac2_ns$eac2_f[i]+eac2_ns$eac2_g[i])
    )+eac2_ns$dead_other[i]
    
  }
  
  # Calculate the number of cancer deaths (any stage) in each stage of disease, 
  # both in the presence and absence of screening (where applicable)
  
  ns$combined_cancer_deaths <- ns[ ,16] + ns[ ,24]
  
  gerd_ns$combined_cancer_deaths <- gerd_ns[ ,16] + gerd_ns[ ,24]
  ndbe_ns$combined_cancer_deaths <- ndbe_ns[ ,16] + ndbe_ns[ ,24]
  lgdbe_ns$combined_cancer_deaths <- lgdbe_ns[ ,16] + lgdbe_ns[ ,24]
  hgdbe_ns$combined_cancer_deaths <- hgdbe_ns[ ,16] + hgdbe_ns[ ,24]
  eac1_ns$combined_cancer_deaths <- eac1_ns[ ,16] + eac1_ns[ ,24]
  eac2_ns$combined_cancer_deaths <- eac2_ns[ ,16] + eac2_ns[ ,24]
  
  ndbe_s$combined_cancer_deaths <- ndbe_s[ ,16] + ndbe_s[ ,24]
  lgdbe_s$combined_cancer_deaths <- lgdbe_s[ ,16] + lgdbe_s[ ,24]
  hgdbe_s$combined_cancer_deaths <- hgdbe_s[ ,16] + hgdbe_s[ ,24]
  eac1_s$combined_cancer_deaths <- eac1_s[ ,16] + eac1_s[ ,24]
  
  # Revert back to number of people in each stage of disease before
  # accounting for sensitivity
  
  ndbe_totals_ns <- ndbe_ns / (1-sensitivity_ndbe)
  lgdbe_totals_ns <- lgdbe_ns / (1-sensitivity_lgdbe)
  hgdbe_totals_ns <- hgdbe_ns / (1-sensitivity_hgdbe)
  eac1_totals_ns <- eac1_ns / (1-sensitivity_eac1)
  
  # Calculate the number of cancer deaths (any stage) prevented by baseline
  # stage of disease. Includes people screening positive, screening negative,
  # and not undergoing screening
  
  cancer_deaths_from_ndbe_prev <- ndbe_totals_ns[,26] - (uptake*ndbe_s[ ,26]+(uptake*ndbe_ns[ ,26])+(1-uptake)*ndbe_totals_ns[ ,26])
  
  cancer_deaths_from_lgdbe_prev <- lgdbe_totals_ns[,26] - (uptake*lgdbe_s[ ,26]+(uptake*lgdbe_ns[ ,26])+(1-uptake)*lgdbe_totals_ns[ ,26])
  
  cancer_deaths_from_hgdbe_prev <- hgdbe_totals_ns[,26] - (uptake*hgdbe_s[ ,26]+(uptake*hgdbe_ns[ ,26])+(1-uptake)*hgdbe_totals_ns[ ,26])
  
  cancer_deaths_from_eac1_prev <- eac1_totals_ns[,26] - (uptake*eac1_s[ ,26]+(uptake*eac1_ns[ ,26])+(1-uptake)*eac1_totals_ns[ ,26])
  
  # Calculate total number of cancer deaths (any stage) prevented by screening
  
  cancer_deaths_prevented <- cancer_deaths_from_ndbe_prev+cancer_deaths_from_lgdbe_prev+cancer_deaths_from_hgdbe_prev+cancer_deaths_from_eac1_prev
  
  # Calculate cancer deaths (any stage) remaining among people invited for screening
  
  cancer_deaths_invited <- ns[,26]-cancer_deaths_prevented
  
  # Calculate number of cancer deaths (any stage) in absence of screening
  
  cancer_deaths <- ns[,26]
  
  # Simulate cumulative number of EAC 2+ cancers, first in the absence of
  # screening
  
  # First, start by resetting select parameters
  
  new_eac2_ns0 <- ns[1,17]
  
  # Second, create holder for for-loop
  
  n_new_eac2_ns <- data.frame(
    new_eac2_ns = c(new_eac2_ns0,rep(NA,12))
  )
  
  # Finally, simulate number of people each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    n_new_eac2_ns$new_eac2_ns[j] <- (n_new_eac2_ns$new_eac2_ns[i]+(ns[j,17]))
    
  }
  
  # Next, simulate cumulative number of EAC 2+ cancers in the presence of
  # screening

  # First, start by resetting select parameters
  
  new_eac2_s0 <- ns[1,17]
  
  # Second, create holder for for-loop
  
  n_new_eac2_s <- data.frame(
    new_eac2_s = c(new_eac2_s0,rep(NA,12))
  )
  
  # Finally, simulate number of people each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    n_new_eac2_s$new_eac2_s[j] <- (n_new_eac2_s$new_eac2_s[i]+(
      (gerd_ns[j,17]+((ndbe_s[j,17]+lgdbe_s[j,17]+hgdbe_s[j,17]+eac1_s[j,17])*uptake)+((ndbe_ns[j,17]+lgdbe_ns[j,17]+hgdbe_ns[j,17]+eac1_ns[j,17])*uptake)+eac2_ns[j,17])))
    
  }
  
  # Cumulative number of EAC 2+ cancers by stage of disease
  
  # First, create new columns for for-loop
  
  ns$eac2_cumul <- c(ns[1,17],rep(0,12))
  
  ndbe_totals_ns$eac2_cumul <- c(ndbe_totals_ns[1,17],rep(0,12))
  lgdbe_totals_ns$eac2_cumul <- c(lgdbe_totals_ns[1,17],rep(0,12))
  hgdbe_totals_ns$eac2_cumul <- c(hgdbe_totals_ns[1,17],rep(0,12))
  eac1_totals_ns$eac2_cumul <- c(eac1_totals_ns[1,17],rep(0,12))
  
  ndbe_s$eac2_cumul <-  c(ndbe_s[1,17],rep(0,12))
  lgdbe_s$eac2_cumul <-  c(lgdbe_s[1,17],rep(0,12))
  hgdbe_s$eac2_cumul <-  c(hgdbe_s[1,17],rep(0,12))
  eac1_s$eac2_cumul <-  c(eac1_s[1,17],rep(0,12))
  
  gerd_ns$eac2_cumul <- c(gerd_ns[1,17],rep(0,12))
  ndbe_ns$eac2_cumul <-c(ndbe_ns[1,17],rep(0,12))
  lgdbe_ns$eac2_cumul <- c(lgdbe_ns[1,17],rep(0,12))
  hgdbe_ns$eac2_cumul <- c(hgdbe_ns[1,17],rep(0,12))
  eac1_ns$eac2_cumul <- c(eac1_ns[1,17],rep(0,12))
  eac2_ns$eac2_cumul <- c(eac2_ns[1,17],rep(0,12))
  
  # Finally, simulate number of people in each stage of disease each year
  
  for (i in 1:12) {
    
    j <- i+1
    
    ns$eac2_cumul[j] <- ns$eac2_cumul[i] + ns$eac2_a[j]
    
    ndbe_totals_ns$eac2_cumul[j] <-  ndbe_totals_ns$eac2_cumul[i] + ndbe_totals_ns$eac2_a[j]
    lgdbe_totals_ns$eac2_cumul[j] <-  lgdbe_totals_ns$eac2_cumul[i] + lgdbe_totals_ns$eac2_a[j]
    hgdbe_totals_ns$eac2_cumul[j] <-  hgdbe_totals_ns$eac2_cumul[i] + hgdbe_totals_ns$eac2_a[j]
    eac1_totals_ns$eac2_cumul[j] <-  eac1_totals_ns$eac2_cumul[i] + eac1_totals_ns$eac2_a[j]
    
    ndbe_s$eac2_cumul[j] <-  ndbe_s$eac2_cumul[i] + ndbe_s$eac2_a[j]
    lgdbe_s$eac2_cumul[j] <-  lgdbe_s$eac2_cumul[i] + lgdbe_s$eac2_a[j]
    hgdbe_s$eac2_cumul[j] <-  hgdbe_s$eac2_cumul[i] + hgdbe_s$eac2_a[j]
    eac1_s$eac2_cumul[j] <-  eac1_s$eac2_cumul[i] + eac1_s$eac2_a[j]
    
    gerd_ns$eac2_cumul[j] <-  gerd_ns$eac2_cumul[i] + gerd_ns$eac2_a[j]
    ndbe_ns$eac2_cumul[j] <-  ndbe_ns$eac2_cumul[i] + ndbe_ns$eac2_a[j]
    lgdbe_ns$eac2_cumul[j] <-  lgdbe_ns$eac2_cumul[i] + lgdbe_ns$eac2_a[j]
    hgdbe_ns$eac2_cumul[j] <-  hgdbe_ns$eac2_cumul[i] + hgdbe_ns$eac2_a[j]
    eac1_ns$eac2_cumul[j] <-  eac1_ns$eac2_cumul[i] + eac1_ns$eac2_a[j]
    eac2_ns$eac2_cumul[j] <-  eac2_ns$eac2_cumul[i] + eac2_ns$eac2_a[j]
    
  }
  
  # Calculate the number of EAC 2+ cancers prevented by baseline
  # stage of disease. Includes people screening positive, screening negative,
  # and not undergoing screening
  
  eac2_from_ndbe_prev <- ndbe_totals_ns[,27] - (uptake*ndbe_s[ ,27]+(uptake*ndbe_ns[ ,27])+(1-uptake)*ndbe_totals_ns[ ,27])
  
  eac2_from_lgdbe_prev <- lgdbe_totals_ns[,27] - (uptake*lgdbe_s[ ,27]+(uptake*lgdbe_ns[ ,27])+(1-uptake)*lgdbe_totals_ns[ ,27])
  
  eac2_from_hgdbe_prev <- hgdbe_totals_ns[,27] - (uptake*hgdbe_s[ ,27]+(uptake*hgdbe_ns[ ,27])+(1-uptake)*hgdbe_totals_ns[ ,27])
  
  eac2_from_eac1_prev <- eac1_totals_ns[,27] - (uptake*eac1_s[ ,27]+(uptake*eac1_ns[ ,27])+(1-uptake)*eac1_totals_ns[ ,27])
  
  # Calculate total number of EAC 2+ cancers prevented by screening
  
  eac2_prevented <- eac2_from_ndbe_prev+eac2_from_lgdbe_prev+eac2_from_hgdbe_prev+eac2_from_eac1_prev
  
  # Calculate number of EAC 2+ cancers that remain among people invited for
  # screening
  
  eac2_invited <- ns[,27]-eac2_prevented
  
  # Calculate number of EAC 2+ cancers occurring in the absence of screening
  
  eac2_absent_screening <- ns[,27]
  
  # Calculate the numbers detected at baseline
  
  n_endo <- (ndbe_s[1,2])+(lgdbe_s[1,3]+hgdbe_s[1,4]+eac1_s[1,5]+eac2_ns[1,17])
  n_screen_detected_ndbe <- ndbe_s[1,2]
  n_screen_detected_lgdbe <- lgdbe_s[1,3]
  n_screen_detected_hgdbe <- hgdbe_s[1,4]
  n_screen_detected_eac <- eac1_s[1,5]+eac2_ns[1,17]
  
  # Simulate the combined incidence (stage 2+) and mortality (any) endpoint
  
  for (j in 1:13) {
    
    ns$incidence_outcome[j] <-  ns$eac2_cumul[j] + ns$dead_eac1[j] 
    
    ndbe_totals_ns$incidence_outcome[j] <-   ndbe_totals_ns$eac2_cumul[j] + ndbe_totals_ns$dead_eac1[j] 
    lgdbe_totals_ns$incidence_outcome[j] <-   lgdbe_totals_ns$eac2_cumul[j] + lgdbe_totals_ns$dead_eac1[j] 
    hgdbe_totals_ns$incidence_outcome[j] <-  hgdbe_totals_ns$eac2_cumul[j] + hgdbe_totals_ns$dead_eac1[j] 
    eac1_totals_ns$incidence_outcome[j] <-   eac1_totals_ns$eac2_cumul[j] + eac1_totals_ns$dead_eac1[j] 
    
    ndbe_s$incidence_outcome[j] <-   ndbe_s$eac2_cumul[j] + ndbe_s$dead_eac1[j] 
    lgdbe_s$incidence_outcome[j] <-   lgdbe_s$eac2_cumul[j] + lgdbe_s$dead_eac1[j] 
    hgdbe_s$incidence_outcome[j] <-   hgdbe_s$eac2_cumul[j] + hgdbe_s$dead_eac1[j] 
    eac1_s$incidence_outcome[j] <-   eac1_s$eac2_cumul[j] + eac1_s$dead_eac1[j] 
    
    gerd_ns$incidence_outcome[j] <-   gerd_ns$eac2_cumul[j] + gerd_ns$dead_eac1[j] 
    ndbe_ns$incidence_outcome[j] <-   ndbe_ns$eac2_cumul[j] + ndbe_ns$dead_eac1[j] 
    lgdbe_ns$incidence_outcome[j] <-   lgdbe_ns$eac2_cumul[j] + lgdbe_ns$dead_eac1[j] 
    hgdbe_ns$incidence_outcome[j] <-   hgdbe_ns$eac2_cumul[j] + hgdbe_ns$dead_eac1[j] 
    eac1_ns$incidence_outcome[j] <-   eac1_ns$eac2_cumul[j] + eac1_ns$dead_eac1[j] 
    eac2_ns$incidence_outcome[j] <-   eac2_ns$eac2_cumul[j] + eac2_ns$dead_eac1[j] 
    
  }
  
  # Calculate the combined incidence and mortality prevented by disease stage
  
  incidence_outcome_from_ndbe_prev <- ndbe_totals_ns[,28] - (uptake*ndbe_s[ ,28]+(uptake*ndbe_ns[ ,28])+(1-uptake)*ndbe_totals_ns[ ,28])
  
  incidence_outcome_from_lgdbe_prev <- lgdbe_totals_ns[,28] - (uptake*lgdbe_s[ ,28]+(uptake*lgdbe_ns[ ,28])+(1-uptake)*lgdbe_totals_ns[ ,28])
  
  incidence_outcome_from_hgdbe_prev <- hgdbe_totals_ns[,28] - (uptake*hgdbe_s[ ,28]+(uptake*hgdbe_ns[ ,28])+(1-uptake)*hgdbe_totals_ns[ ,28])
  
  incidence_outcome_from_eac1_prev <- eac1_totals_ns[,28] - (uptake*eac1_s[ ,28]+(uptake*eac1_ns[ ,28])+(1-uptake)*eac1_totals_ns[ ,28])
  
  # Calculate the total combined incidence and mortality prevented by screening
  
  incidence_outcome_prevented <- incidence_outcome_from_ndbe_prev+incidence_outcome_from_lgdbe_prev+incidence_outcome_from_hgdbe_prev+incidence_outcome_from_eac1_prev
  
  # Power calculation for mortality-only endpoint
  
  # Control (cancer deaths) 
  
  control_6 <- sum(ns[4:7,26])/4
  control_75 <- sum(ns[7:9,26])/3
  control_12 <- sum(ns[10:13,26])/4
  
  # Intervention (cancer deaths) 
  
  intervention_6 <- control_6-sum(cancer_deaths_prevented[4:7])/4
  intervention_75 <- control_75-sum(cancer_deaths_prevented[7:9])/3
  intervention_12 <- control_12-sum(cancer_deaths_prevented[10:13])/4
  
  # Prep for power calculations 
  
  P0_6<-control_6*0.00001
  P0_75<-control_75*0.00001
  P0_12<-control_12*0.00001
  P1_6<-intervention_6*0.00001
  P1_75<-intervention_75*0.00001
  P1_12<-intervention_12*0.00001
  delta_6<-P0_6-P1_6
  delta_75<-P0_75-P1_75
  delta_12<-P0_12-P1_12
  p_pool_6<-(2*P0_6+P1_6)/3
  p_pool_75<-(2*P0_75+P1_75)/3
  p_pool_12<-(2*P0_12+P1_12)/3
  sd_diff_6<-sqrt(P0_6*(1-P0_6)/(2*N/3)+P1_6*(1-P1_6)/(N/3))
  sd_diff_75<-sqrt(P0_75*(1-P0_75)/(2*N/3)+P1_75*(1-P1_75)/(N/3))
  sd_diff_12<-sqrt(P0_12*(1-P0_12)/(2*N/3)+P1_12*(1-P1_12)/(N/3))
  sd_pool_6<-sqrt(p_pool_6*(1-p_pool_6)*4.5/N)
  sd_pool_75<-sqrt(p_pool_75*(1-p_pool_75)*4.5/N)
  sd_pool_12<-sqrt(p_pool_12*(1-p_pool_12)*4.5/N)
  
  # Power calculations
  
  power_6<-pnorm((delta_6-qnorm(1-alpha1/2,0,1)*sd_pool_6)/sd_diff_6)-pnorm((-delta_6-qnorm(1-alpha1/2,0,1)*sd_pool_6)/sd_diff_6)
  power_75<-pnorm((delta_75-qnorm(1-alpha2/2,0,1)*sd_pool_75)/sd_diff_75)-pnorm((-delta_75-qnorm(1-alpha2/2,0,1)*sd_pool_75)/sd_diff_75)
  power_12<-pnorm((delta_12-qnorm(1-alpha3/2,0,1)*sd_pool_12)/sd_diff_12)-pnorm((-delta_12-qnorm(1-alpha3/2,0,1)*sd_pool_12)/sd_diff_12)
  power<-c(power_6,power_75,power_12)
  
  # Power calculation for combined incidence and mortality endpoint
  
  # Control (any mortality and stage 2+ incidence)
  
  i_control_4 <- sum(ns[2:5,28])/4
  i_control_55 <- sum(ns[5:7,28])/3
  i_control_85 <- sum(ns[8:10,28])/3
  
  # Intervention (any mortality and stage 2+ incidence)
  
  i_intervention_4 <- i_control_4-sum(incidence_outcome_prevented[2:5])/4
  i_intervention_55 <- i_control_55-sum(incidence_outcome_prevented[5:7])/3
  i_intervention_85 <- i_control_85-sum(incidence_outcome_prevented[8:10])/3
  
  # Prep for power calculations 
  
  P0_4<-i_control_4*0.00001
  P0_55<-i_control_55*0.00001
  P0_85<-i_control_85*0.00001
  P1_4<-i_intervention_4*0.00001
  P1_55<-i_intervention_55*0.00001
  P1_85<-i_intervention_85*0.00001
  
  delta_4<-P0_4-P1_4
  delta_55<-P0_55-P1_55
  delta_85<-P0_85-P1_85
  
  p_pool_4<-(2*P0_4+P1_4)/3
  p_pool_55<-(2*P0_55+P1_55)/3
  p_pool_85<-(2*P0_85+P1_85)/3
  
  sd_diff_4<-sqrt(P0_4*(1-P0_4)/(2*N/3)+P1_4*(1-P1_4)/(N/3))
  sd_diff_55<-sqrt(P0_55*(1-P0_55)/(2*N/3)+P1_55*(1-P1_55)/(N/3))
  sd_diff_85<-sqrt(P0_85*(1-P0_85)/(2*N/3)+P1_85*(1-P1_85)/(N/3))
  
  sd_pool_4<-sqrt(p_pool_4*(1-p_pool_4)*4.5/N)
  sd_pool_55<-sqrt(p_pool_55*(1-p_pool_55)*4.5/N)
  sd_pool_85<-sqrt(p_pool_85*(1-p_pool_85)*4.5/N)
  
  # Power calculations
  
  i_power_4 <- pnorm((delta_4-qnorm(1-alpha1/2,0,1)*sd_pool_4)/sd_diff_4)-pnorm((-delta_4-qnorm(1-alpha1/2,0,1)*sd_pool_4)/sd_diff_4)
  i_power_55<-pnorm((delta_55-qnorm(1-alpha2/2,0,1)*sd_pool_55)/sd_diff_55)-pnorm((-delta_55-qnorm(1-alpha2/2,0,1)*sd_pool_55)/sd_diff_55)
  i_power_85 <-pnorm((delta_85-qnorm(1-alpha3/2,0,1)*sd_pool_85)/sd_diff_85)-pnorm((-delta_85-qnorm(1-alpha3/2,0,1)*sd_pool_85)/sd_diff_85)
  i_power<-c(i_power_4,i_power_55,i_power_85)
  
  # Number screened
  
  n_screened <- N*uptake
  
  # Output
  
  output_list <- list(
    power
    ,i_power
    ,gerd_ns
    ,ndbe_ns
    ,ndbe_s
    ,lgdbe_ns
    ,lgdbe_s
    ,hgdbe_ns
    ,hgdbe_s
    ,eac2_absent_screening # Control
    ,eac2_invited # intervention
    ,eac2_from_ndbe_prev
    ,eac2_from_lgdbe_prev
    ,eac2_from_hgdbe_prev
    ,eac2_from_eac1_prev
    ,cancer_deaths # Control
    ,cancer_deaths_invited # intervention
    ,cancer_deaths_from_ndbe_prev
    ,cancer_deaths_from_lgdbe_prev
    ,cancer_deaths_from_hgdbe_prev
    ,cancer_deaths_from_eac1_prev
    ,incidence_outcome_from_ndbe_prev
    ,incidence_outcome_from_lgdbe_prev
    ,incidence_outcome_from_hgdbe_prev
    ,incidence_outcome_from_eac1_prev
    ,eac1_ns
    ,eac1_s
    ,eac2_ns
  )
  
}