library(rstan)

outcome_reg <- rstan::stan_model("./sensitivity_change_chars.stan",
                                 model_name="time_change_sensitivity",
                                 warn_pedantic=TRUE,
                                 save_dso=TRUE)

saveRDS(outcome_reg, file="sensitivity_change_chars.RDS")

specificity_reg <- rstan::stan_model("./specificity_chars.stan",
                                 model_name="specificity",
                                 warn_pedantic=TRUE,
                                 save_dso=TRUE)

saveRDS(specificity_reg, file="specificity_chars.RDS")

