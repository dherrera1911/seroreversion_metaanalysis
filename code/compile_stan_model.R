library(rstan)

outcome_reg <- rstan::stan_model("./sensitivity_change_chars.stan",
                                 model_name="time_change_sensitivity",
                                 warn_pedantic=TRUE,
                                 save_dso=TRUE)

saveRDS(outcome_reg, file="sensitivity_change_chars.RDS")

