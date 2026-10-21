source("sim/sim_func2.r")

run_simulation_study <- function(n_iters = 5, n_obs = 500) {
  all_results <- tibble()
  
  # Define the 3 Test Scenarios
  scenarios <- list(
    "1. Simple" = list(st = "t_max * exp(-(tx * log(0.5) + x1 * log(2)))", 
                       st2 = NULL),
    "2. Interaction" = list(st = "t_max * exp(-(tx * log(0.5) + x1 * log(2) + tx*x1*log(1.5)))", 
                            st2 = NULL),
    "3. Waning" = list(st = "t_max * exp(-(tx * log(0.2) + x1 * log(2)))", # Strong effect initially
                       st2 = "t_max * exp(-(tx * log(0.8) + x1 * log(2)))") # Effect wanes after switch2
  )
  
  for (s_name in names(scenarios)) {
    message("Running Scenario: ", s_name)
    params <- scenarios[[s_name]]
    print("HI")
    
    for (i in 1:n_iters) {
      df_covs <- tibble(id = 1:n_obs, x1 = rnorm(n_obs, 1, 0.5))
       sim <- generate_survival_sim(df_covs, eq_scale_t = params$st, eq_scale_t2 = params$st2)
      
      # Run Models
      res_cph  <- run_cph_suite(sim)
      res_bart <- run_bart_suite(sim)
      res_ccw  <- run_ccw_tmle_suite(sim)
      
      # Combine and compare to truth
      iter_res <- bind_rows(res_cph$estimates, res_bart$estimates, res_ccw$estimates) %>%
        left_join(sim$truth, by = "t") %>%
        mutate(bias = sd_estimate - true_sd, scenario = s_name, iter = i)
      
      all_results <- bind_rows(all_results, iter_res)
    }
  }
  
  # Final Summary
  summary <- all_results %>%
    group_by(scenario, t, method) %>%
    summarise(mean_bias = mean(bias), rmse = sqrt(mean(bias^2)), .groups = "drop")
  
  return(list(raw = all_results, summary = summary))
}

# Execute
results <- run_simulation_study(n_iters = 3, n_obs = 300)
print(results$summary)


i=1
n_obs=300
df_covs <- tibble(id = 1:n_obs, x1 = rnorm(n_obs, 1, 0.5))
sim <- generate_survival_sim(df_covs, eq_scale_t = params$st, eq_scale_t2 = params$st2)


res_cph  <- run_cph_suite(sim)
res_bart <- run_bart_suite(sim)
res_ccw  <- run_ccw_tmle_suite(sim)

# Combine and compare to truth
iter_res <- bind_rows(res_cph$estimates, res_bart$estimates, res_ccw$estimates) %>%
left_join(sim$truth, by = "t") %>%
mutate(bias = sd_estimate - true_sd, scenario = s_name, iter = i)

all_results <- bind_rows(all_results, iter_res)
