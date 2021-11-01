
library(survival)
library(tidyverse)
library(survminer)
library(cmdstanr)
library(rstan)

set.seed(201111)

NoCov_stan        <- cmdstan_model(paste0(getwd(),'/Stan_Surv/NoCov.stan'))
NoCov_EX_stan     <- cmdstan_model(paste0(getwd(),'/Stan_Surv/NoCov_EX.stan'))
NoCov_EXNEX_stan  <- cmdstan_model(paste0(getwd(),'/Stan_Surv/NoCov_EXNEX.stan'))
NoCov_HS_stan     <- cmdstan_model(paste0(getwd(),'/Stan_Surv/NoCov_HS.stan'))

Cov_stan        <- cmdstan_model(paste0(getwd(),'/Stan_Surv_Cov/Cov.stan'))
Cov_EX_stan     <- cmdstan_model(paste0(getwd(),'/Stan_Surv_Cov/Cov_EX.stan'))
Cov_EXNEX_stan  <- cmdstan_model(paste0(getwd(),'/Stan_Surv_Cov/Cov_EXNEX.stan'))
Cov_HS_stan     <- cmdstan_model(paste0(getwd(),'/Stan_Surv_Cov/Cov_HS.stan'))

base.haz <- 0.3
n_H1 <- 283; n_H2 <- 203; n_H3 <- 450; n_CC <- 42; n_CT <- 47
n_tot <- n_H1+n_H2+n_H3+n_CC+n_CT

ANA <- data.frame(
  TRT = c(rep(0, n_H1+n_H2+n_H3+n_CC), rep(1, n_CT)),
  SID = c(rep(1, n_H1), rep(2, n_H2), rep(3, n_H3), rep(4, n_CC+n_CT)), 
  X = c(rnorm(n_H1, 2, 1), c(rnorm(n_tot-n_H1, 0, 1)))
) %>% mutate(haz = base.haz * exp(log(0.8)*TRT + log(0.7)*X)) %>% 
  mutate(etime = 1/haz * rexp(n_tot, 1)) %>% mutate(ctime = 1/(haz/4) * rexp(n_tot, 1)) %>%
  mutate(Time = if_else(etime<ctime, etime, ctime)) %>% mutate(Event = if_else(etime<ctime, 1, 0)) %>% 
  mutate(S_TRT = if_else(TRT == 1, 1, if_else(SID == 4, 2, if_else(SID == 1, 3, if_else(SID == 2, 4, 5)))))

table(ANA$S_TRT, ANA$Event)

# Study, TRT
sf <- survfit(Surv(Time,Event)~S_TRT,data = ANA)
fortify.survfit <- function(survfit.data) {
  data.frame(time = survfit.data$time,
             n.risk = survfit.data$n.risk,
             n.event = survfit.data$n.event,
             n.censor = survfit.data$n.censor,
             surv = survfit.data$surv,
             std.err = survfit.data$std.err,
             upper = survfit.data$upper,
             lower = survfit.data$lower,
             strata = rep(names(survfit.data$strata), survfit.data$strata))
}
survdiff(Surv(Time,Event)~S_TRT,data = ANA)
p2 <- ggsurvplot(fit = sf, data = ANA, censor=T, 
                 size = 1.1,
                 censor.shape = c(124),
                 censor.size = c(4),
                 ylab = "Probability of PFS",
                 xlab = "Month", 
                 legend.title = "Group", 
                 legend.labs = c( "Current treatment","Current control", "Historical 1", "Historical 2", "Historical 3"),
                 linetype = c("strata"),
                 palette = c("black", "black", "#1976D2","#81C784","#FFB74D"),
                 font.legend = c(18 ,colour =  "black"),
                 font.x = c(18 ,colour =  "black"),
                 font.y = c(18 ,colour =  "black"),
                 font.tickslab = c(18 ,colour =  "black"),
                 legend = c(0.7, 0.7)
)
p2$plot <- p2$plot + 
  scale_linetype_manual(values = c(1,2,3,5,6)) +
  theme(
    legend.key.width = unit(2.5,"cm")
  )
p2$plot



ANA01_1 <- filter(ANA,SID==1,TRT==0)
ANA01_2 <- filter(ANA,SID==2,TRT==0)
ANA01_3 <- filter(ANA,SID==3,TRT==0)
ANA01_4 <- filter(ANA,SID==4)
ANA01 <- rbind(ANA01_1,ANA01_2,ANA01_3,ANA01_4)
ANA02 <- ANA01[, colnames(ANA01) != "Study"]

ANA_C    <- filter(ANA02, SID==4)
ANA_ALL  <- filter(ANA02)

ANA_C_E <- filter(ANA_C, Event == 1)
ANA_C_C <- filter(ANA_C, Event == 0)

ANA_ALL_E <- filter(ANA_ALL, Event == 1)
ANA_ALL_C <- filter(ANA_ALL, Event == 0)

H <- 3

# not adjusted covariate
# Current: ANA_C
I_E <- nrow(ANA_C_E); I_C <- nrow(ANA_C_C);
ANA_data <- list(I_E = I_E, I_C = I_C, TRT_E = if_else(ANA_C_E$TRT==1,0,1), TRT_C = if_else(ANA_C_C$TRT==1,0,1), Time_E = ANA_C_E$Time, Time_C = ANA_C_C$Time)
fit_Current  <- NoCov_stan$sample(data=ANA_data, seed = 954, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 4000, refresh = 0, show_messages = F)

# Pooled: ANA_ALL
I_E <- nrow(ANA_ALL_E); I_C <- nrow(ANA_ALL_C);
ANA_data <- list(I_E = I_E, I_C = I_C, TRT_E = if_else(ANA_ALL_E$TRT==1,0,1), TRT_C = if_else(ANA_ALL_C$TRT==1,0,1), Time_E = ANA_ALL_E$Time, Time_C = ANA_ALL_C$Time)
fit_Pooled  <- NoCov_stan$sample(data=ANA_data, seed = 954, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 4000, refresh = 0, show_messages = F)

# EX
I_E <- nrow(ANA_ALL_E); I_C <- nrow(ANA_ALL_C);
ANA_data <- list(I_E = I_E, I_C = I_C, H = H, SID_E = ANA_ALL_E$SID, SID_C = ANA_ALL_C$SID, Tau = 0.5,  TRT_E = if_else(ANA_ALL_E$TRT==1,0,1), TRT_C = if_else(ANA_ALL_C$TRT==1,0,1), Time_E = ANA_ALL_E$Time, Time_C = ANA_ALL_C$Time)
fit_EX   <- NoCov_EX_stan$sample(data=ANA_data, seed = 954, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 4000, adapt_delta = 0.9, refresh = 0, show_messages = F)

# EXNEX
I_E <- nrow(ANA_ALL_E); I_C <- nrow(ANA_ALL_C);
ANA_data <- list(I_E = I_E, I_C = I_C, H = H, SID_E = ANA_ALL_E$SID, SID_C = ANA_ALL_C$SID, Tau = 0.5,  TRT_E = if_else(ANA_ALL_E$TRT==1,0,1), TRT_C = if_else(ANA_ALL_C$TRT==1,0,1), Time_E = ANA_ALL_E$Time, Time_C = ANA_ALL_C$Time)
fit_EXNEX   <- NoCov_EXNEX_stan$sample(data=ANA_data, seed = 954, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 6000, adapt_delta = 0.95, refresh = 0, show_messages = F)


# HS
I_E <- nrow(ANA_ALL_E); I_C <- nrow(ANA_ALL_C);
ANA_data <- list(I_E = I_E, I_C = I_C, H = H, betascale = 1, nu = 1, SID_E = ANA_ALL_E$SID, SID_C = ANA_ALL_C$SID, Tau = 0.5,  TRT_E = if_else(ANA_ALL_E$TRT==1,0,1), TRT_C = if_else(ANA_ALL_C$TRT==1,0,1), Time_E = ANA_ALL_E$Time, Time_C = ANA_ALL_C$Time)
fit_HS   <- NoCov_HS_stan$sample(data=ANA_data, seed = 954, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 4000, adapt_delta = 0.9, refresh = 0, show_messages = F)

# Method, Mean, SD, LCI, UCI
Result_NoCov <- matrix(0, nrow = 5, ncol = 6)
Result_NoCov[1,]  <- c(1, mean(exp(-fit_Current$draws("eta_CC"))), 
                       sd(exp(-fit_Current$draws("eta_CC"))), 
                       quantile(exp(-fit_Current$draws("eta_CC")),0.025,names=F),
                       quantile(exp(-fit_Current$draws("eta_CC")),0.975,names=F), 
                       max(fit_Current$summary()$rhat)  )
Result_NoCov[2,]  <- c(2, mean(exp(-fit_Pooled$draws("eta_CC"))), 
                       sd(exp(-fit_Pooled$draws("eta_CC"))), 
                       quantile(exp(-fit_Pooled$draws("eta_CC")),0.025,names=F), 
                       quantile(exp(-fit_Pooled$draws("eta_CC")),0.975,names=F), 
                       max(fit_Pooled$summary()$rhat)  )
Result_NoCov[3,]  <- c(3, mean(exp(-fit_EX$draws("eta_CC"))), 
                       sd(exp(-fit_EX$draws("eta_CC"))), 
                       quantile(exp(-fit_EX$draws("eta_CC")),0.025,names=F),
                       quantile(exp(-fit_EX$draws("eta_CC")),0.975,names=F), 
                       max(fit_EX$summary()$rhat)  )
Result_NoCov[4,]  <- c(4, mean(exp(-fit_EXNEX$draws("eta_CC"))), 
                       sd(exp(-fit_EXNEX$draws("eta_CC"))), 
                       quantile(exp(-fit_EXNEX$draws("eta_CC")),0.025,names=F), 
                       quantile(exp(-fit_EXNEX$draws("eta_CC")),0.975,names=F), 
                       max(fit_EXNEX$summary()$rhat)  )
Result_NoCov[5,]  <- c(5, mean(exp(-fit_HS$draws("eta_CC"))), 
                       sd(exp(-fit_HS$draws("eta_CC"))), 
                       quantile(exp(-fit_HS$draws("eta_CC")),0.025,names=F), 
                       quantile(exp(-fit_HS$draws("eta_CC")),0.975,names=F), 
                       max(fit_HS$summary()$rhat)  )


Result_NoCov1 <- as.data.frame(Result_NoCov,stringsAsFactors = F) 
colnames(Result_NoCov1) <- c("MethodN", "Mean", "SD", "LCI", "UCI", "diag")
Result_NoCov1 <- Result_NoCov1 %>% dplyr::mutate(Range = (UCI-LCI)) %>%
  mutate(Method = case_when(
    MethodN == 1 ~ "Current_NoCov",
    MethodN == 2 ~ "Pooled_NoCov",
    MethodN == 3 ~ "EX_NoCov",
    MethodN == 4 ~ "EXNEX_NoCov",
    MethodN == 5 ~ "Horseshoe_NoCov",
    TRUE ~ ""
  )) 

Result_NoCov2 <- Result_NoCov1[,c(8,2,3,4,5,7,6)]
Result_NoCov2
# write.csv(Result_NoCov2, "Result_Surv.csv")


# adjusted covariate
# Current: ANA_C
I_E <- nrow(ANA_C_E); I_C <- nrow(ANA_C_C);
ANA_data <- list(I_E = I_E, I_C = I_C, Ncov = 1, X_E = matrix(ANA_C_E$X, nrow = I_E), X_C = matrix(ANA_C_C$X, nrow=I_C), TRT_E = if_else(ANA_C_E$TRT==1,0,1), TRT_C = if_else(ANA_C_C$TRT==1,0,1), Time_E = ANA_C_E$Time, Time_C = ANA_C_C$Time)
fit_Current_Cov  <- Cov_stan$sample(data=ANA_data, seed = 954, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 4000, refresh = 0, show_messages = F)

# Pooled: ANA_ALL
I_E <- nrow(ANA_ALL_E); I_C <- nrow(ANA_ALL_C);
ANA_data <- list(I_E = I_E, I_C = I_C, Ncov = 1, X_E = matrix(ANA_ALL_E$X, nrow = I_E), X_C = matrix(ANA_ALL_C$X, nrow=I_C), TRT_E = if_else(ANA_ALL_E$TRT==1,0,1), TRT_C = if_else(ANA_ALL_C$TRT==1,0,1), Time_E = ANA_ALL_E$Time, Time_C = ANA_ALL_C$Time)
fit_Pooled_Cov  <- Cov_stan$sample(data=ANA_data, seed = 954, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 4000, refresh = 0, show_messages = F)

# EX
I_E <- nrow(ANA_ALL_E); I_C <- nrow(ANA_ALL_C);
ANA_data <- list(I_E = I_E, I_C = I_C, Ncov = 1, X_E = matrix(ANA_ALL_E$X, nrow = I_E), X_C = matrix(ANA_ALL_C$X, nrow=I_C), H = H, SID_E = ANA_ALL_E$SID, SID_C = ANA_ALL_C$SID, Tau = 0.5,  TRT_E = if_else(ANA_ALL_E$TRT==1,0,1), TRT_C = if_else(ANA_ALL_C$TRT==1,0,1), Time_E = ANA_ALL_E$Time, Time_C = ANA_ALL_C$Time)
fit_EX_Cov   <- Cov_EX_stan$sample(data=ANA_data, seed = 954, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 4000, adapt_delta = 0.9, refresh = 0, show_messages = F)

# EXNEX
I_E <- nrow(ANA_ALL_E); I_C <- nrow(ANA_ALL_C);
ANA_data <- list(I_E = I_E, I_C = I_C, Ncov = 1, X_E = matrix(ANA_ALL_E$X, nrow = I_E), X_C = matrix(ANA_ALL_C$X, nrow=I_C), H = H, SID_E = ANA_ALL_E$SID, SID_C = ANA_ALL_C$SID, Tau = 0.5,  TRT_E = if_else(ANA_ALL_E$TRT==1,0,1), TRT_C = if_else(ANA_ALL_C$TRT==1,0,1), Time_E = ANA_ALL_E$Time, Time_C = ANA_ALL_C$Time)
fit_EXNEX_Cov   <- Cov_EXNEX_stan$sample(data=ANA_data, seed = 954, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 6000, adapt_delta = 0.95, refresh = 0, show_messages = F)

# HS
I_E <- nrow(ANA_ALL_E); I_C <- nrow(ANA_ALL_C);
ANA_data <- list(I_E = I_E, I_C = I_C, Ncov = 1, X_E = matrix(ANA_ALL_E$X, nrow = I_E), X_C = matrix(ANA_ALL_C$X, nrow=I_C), H = H, betascale = 1, nu = 1, SID_E = ANA_ALL_E$SID, SID_C = ANA_ALL_C$SID, Tau = 0.5,  TRT_E = if_else(ANA_ALL_E$TRT==1,0,1), TRT_C = if_else(ANA_ALL_C$TRT==1,0,1), Time_E = ANA_ALL_E$Time, Time_C = ANA_ALL_C$Time)
fit_HS_Cov   <- Cov_HS_stan$sample(data=ANA_data, seed = 954, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 4000, adapt_delta = 0.9, refresh = 0, show_messages = F)

# Method, Mean, SD, LCI, UCI
Result_Cov <- matrix(0, nrow = 5, ncol = 6)
Result_Cov[1,]  <- c(1, mean(exp(-fit_Current_Cov$draws("eta_CC"))), 
                       sd(exp(-fit_Current_Cov$draws("eta_CC"))), 
                       quantile(exp(-fit_Current_Cov$draws("eta_CC")),0.025,names=F),
                       quantile(exp(-fit_Current_Cov$draws("eta_CC")),0.975,names=F), 
                       max(fit_Current_Cov$summary()$rhat)  )
Result_Cov[2,]  <- c(2, mean(exp(-fit_Pooled_Cov$draws("eta_CC"))), 
                       sd(exp(-fit_Pooled_Cov$draws("eta_CC"))), 
                       quantile(exp(-fit_Pooled_Cov$draws("eta_CC")),0.025,names=F), 
                       quantile(exp(-fit_Pooled_Cov$draws("eta_CC")),0.975,names=F), 
                       max(fit_Pooled_Cov$summary()$rhat)  )
Result_Cov[3,]  <- c(3, mean(exp(-fit_EX_Cov$draws("eta_CC"))), 
                       sd(exp(-fit_EX_Cov$draws("eta_CC"))), 
                       quantile(exp(-fit_EX_Cov$draws("eta_CC")),0.025,names=F),
                       quantile(exp(-fit_EX_Cov$draws("eta_CC")),0.975,names=F), 
                       max(fit_EX_Cov$summary()$rhat)  )
Result_Cov[4,]  <- c(4, mean(exp(-fit_EXNEX_Cov$draws("eta_CC"))), 
                       sd(exp(-fit_EXNEX_Cov$draws("eta_CC"))), 
                       quantile(exp(-fit_EXNEX_Cov$draws("eta_CC")),0.025,names=F), 
                       quantile(exp(-fit_EXNEX_Cov$draws("eta_CC")),0.975,names=F), 
                       max(fit_EXNEX_Cov$summary()$rhat)  )
Result_Cov[5,]  <- c(5, mean(exp(-fit_HS_Cov$draws("eta_CC"))), 
                       sd(exp(-fit_HS_Cov$draws("eta_CC"))), 
                       quantile(exp(-fit_HS_Cov$draws("eta_CC")),0.025,names=F), 
                       quantile(exp(-fit_HS_Cov$draws("eta_CC")),0.975,names=F), 
                       max(fit_HS_Cov$summary()$rhat)  )


Result_Cov1 <- as.data.frame(Result_Cov,stringsAsFactors = F) 
colnames(Result_Cov1) <- c("MethodN", "Mean", "SD", "LCI", "UCI", "diag")
Result_Cov1 <- Result_Cov1 %>% dplyr::mutate(Range = (UCI-LCI)) %>%
  mutate(Method = case_when(
    MethodN == 1 ~ "Current_Cov",
    MethodN == 2 ~ "Pooled_Cov",
    MethodN == 3 ~ "EX_Cov",
    MethodN == 4 ~ "EXNEX_Cov",
    MethodN == 5 ~ "Horseshoe_Cov",
    TRUE ~ ""
  )) 

Result_Cov2 <- Result_Cov1[,c(8,2,3,4,5,7,6)]
Result_Cov2
# write.csv(Result_Cov2, "Result_Surv.csv")



result_HS <- fit_HS$output_files() %>% rstan::read_stan_csv()
plot100 <- stan_dens(result_HS, pars = c("beta"))$data
plot100_1 <- plot100 %>% filter(parameter %in% c("beta[1]","beta[2]","beta[3]")) %>%
  mutate(color = case_when(
    parameter == "beta[1]" ~ "beta1",
    parameter == "beta[2]" ~ "beta2",
    parameter == "beta[3]" ~ "beta3",
    TRUE ~ " "
  )
  )
p100 <- ggplot(data = plot100_1, aes(x = value, colour = color, linetype = color)) +
  labs(linetype="parameter", color = "parameter") +
  geom_line(stat = "density", size = 1.1) +
  scale_colour_manual(
    values = c(
      "beta1"  = "#1976D2",
      "beta2"  = "#81C784",
      "beta3"  = "#FFB74D"),
    labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]))) +
  scale_linetype_manual(
    values = c("dotted", "longdash", "twodash"),
    labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]))) +
  xlab(expression(beta)) +
  ylab("Density") +
  coord_cartesian(xlim = c(-0.5,0.5), ylim = c(0,15)) +
  scale_x_continuous(breaks=c(-0.5,-0.25,0,0.25,0.5)) +
  theme(
    plot.title = element_text(size = 18),
    axis.text = element_text(colour = "black", size = 18),
    axis.ticks=element_line(colour = "black", size=1),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.position = c(0.8,0.8),
    legend.key.width = unit(2.5,"cm")
  )
p100



result_HS_Cov <- fit_HS_Cov$output_files() %>% rstan::read_stan_csv()
plot200 <- stan_dens(result_HS_Cov, pars = c("beta"))$data
plot200_1 <- plot200 %>% filter(parameter %in% c("beta[1]","beta[2]","beta[3]")) %>%
  mutate(color = case_when(
    parameter == "beta[1]" ~ "beta1",
    parameter == "beta[2]" ~ "beta2",
    parameter == "beta[3]" ~ "beta3",
    TRUE ~ " "
  )
  )
p200 <- ggplot(data = plot200_1, aes(x = value, colour = color, linetype = color)) +
  labs(linetype="parameter", color = "parameter") +
  geom_line(stat = "density", size = 1.1) +
  scale_colour_manual(
    values = c(
      "beta1"  = "#1976D2",
      "beta2"  = "#81C784",
      "beta3"  = "#FFB74D"),
    labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]))) +
  scale_linetype_manual(
    values = c("dotted", "longdash", "twodash"),
    labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]))) +
  xlab(expression(beta)) +
  ylab("Density") +
  coord_cartesian(xlim = c(-0.5,0.5), ylim = c(0,15)) +
  scale_x_continuous(breaks=c(-0.5,-0.25,0,0.25,0.5)) +
  theme(
    plot.title = element_text(size = 18),
    axis.text = element_text(colour = "black", size = 18),
    axis.ticks=element_line(colour = "black", size=1),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.position = c(0.8,0.8),
    legend.key.width = unit(2.5,"cm")
  )
p200

