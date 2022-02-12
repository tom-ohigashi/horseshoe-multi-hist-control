# install.packages(c("VGAM", "minpack.lm"))
# install.packages("StudyPrior", repos="http://R-Forge.R-project.org")
library(stats4)
library(gtools)
library(R.utils)
library(plyr)
library(dplyr)
library(StudyPrior)
library(INLA)
library(quadprog)
library(tidyverse)
library(gridExtra)
library(rstan)
library(cmdstanr)
library(ggpubr)

n_h <- c(43,27,10,5)
x_h <- c(36,8,4,2)

n_CC <- 29
x_CC <- 9
n_CT <- 32
x_CT <- 18

Historical <- data.frame(
  N = n_h, 
  X = x_h
)

Current <- data.frame(
  TRT = c(0,1),
  N = c(n_CC, n_CT),
  X = c(x_CC, x_CT)
)

Data_all <- list(Historical=Historical, Current=Current)

n_Control <- c(Data_all$Historical[,"N"],Data_all$Current[Data_all$Current$TRT==0,"N"])
x_Control <- c(Data_all$Historical[,"X"],Data_all$Current[Data_all$Current$TRT==0,"X"])
Control <- cbind(n_Control,x_Control)
success <- c(Data_all$Historical[,"X"],Data_all$Current[,"X"])
n <- c(Data_all$Historical[,"N"],Data_all$Current[,"N"])
Data_list=list(success=success,H=length(n_h),n=n)

Current_stan <- cmdstan_model(paste0(getwd(),'/Stan_BIN/Current.stan'))
Pooled_stan  <- cmdstan_model(paste0(getwd(),'/Stan_BIN/Pooled.stan'))
EX_stan      <- cmdstan_model(paste0(getwd(),'/Stan_BIN/EX.stan'))
EXNEX_stan   <- cmdstan_model(paste0(getwd(),'/Stan_BIN/EXNEX.stan'))
DMPP_stan    <- cmdstan_model(paste0(getwd(),'/Stan_BIN/DMPP.stan'))
RDMPP_stan    <- cmdstan_model(paste0(getwd(),'/Stan_BIN/RDMPP.stan'))
HS_stan      <- cmdstan_model(paste0(getwd(),'/Stan_BIN/HS.stan'))


H <- length(Data_all$Historical[,"N"])
ANA <- list(H=length(Data_all$Historical[,"N"]),x_h=Data_all$Historical[,"X"],n_h=Data_all$Historical[,"N"],x_CC=Data_all$Current$X[1],n_CC=Data_all$Current$N[1],x_CT=Data_all$Current$X[2],n_CT=Data_all$Current$N[2]
            ,betascale=1, nu=1, C=10)
set.seed(202111)
init_list <- list(list(g =  rnorm(1,0,0.01), theta_C = rnorm(1,0.5,0.01)), list(g =  rnorm(1,0,0.01), theta_C = rnorm(1,0.5,0.01)), list(g =  rnorm(1,0,0.01), theta_C = rnorm(1,0.5,0.01)), list(g =  rnorm(1,0,0.01), theta_C = rnorm(1,0.5,0.01)))
init_list_h <- list(list(g =  rnorm(1,0,0.01)), list(g =  rnorm(1,0,0.01)), list(g =  rnorm(1,0,0.01)), list(g =  rnorm(1,0,0.01)))

# Current data, Pooled data
fit_Current  <- Current_stan$sample(data=ANA, seed = 954, init = init_list, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 4000, refresh = 0, show_messages = F)
fit_Pooled   <- Pooled_stan$sample(data=ANA, seed = 954, init = init_list, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 4000, refresh = 0, show_messages = F)

# EX, EXNEX
fit_EX     <- EX_stan$sample(data=ANA, seed = 954, init = init_list_h, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 4000, adapt_delta = 0.9, refresh = 0, show_messages = F)
fit_EXNEX  <- EXNEX_stan$sample(data=ANA, seed = 954, init = init_list_h, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 4000, adapt_delta = 0.9, refresh = 0, show_messages = F)

# MPP by Banbeta et al (2019)
fit_DMPP   <- DMPP_stan$sample(data=ANA, seed = 954, init = init_list_h, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 6000, adapt_delta = 0.99, refresh = 0, show_messages = F)
fit_RDMPP   <- RDMPP_stan$sample(data=ANA, seed = 954, init = init_list_h, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 6000, adapt_delta = 0.99, refresh = 0, show_messages = F)

# EB methods by Gravestock and Held (2019)
hx.rand <- function(hx,x) {
  y <- c(); i <- 1
  while (i <= x) {
    U <- runif(1); V <-runif(1); W <- hx(V)/10
    if (U < W){
      y <- append(y, V); i <- i+1
    }
  }
  return(y)
}
Calc.posterior <- function(prior, X, N){
  k <- adaptIntegrate(function(p) prior(p)*dbinom(X,N,p),
                      lowerLimit = 0,
                      upperLimit = 1)$integral
  function(p) prior(p)*dbinom(X,N,p)/k
}
niter <- 20000
post.theta_T <- rbeta(0.2*niter,1+Data_all$Current$X[2],1+Data_all$Current$N[2]-Data_all$Current$X[2])
unnecessary.messages <- capture.output(
  cc.eb.pp <- binom.PP.EB(Data_all$Historical[,"X"], Data_all$Historical[,"N"], N=Data_all$Current$N[1],X=Data_all$Current$X[1])
)
p.cc.eb.pp <- Calc.posterior(function(p)  cc.eb.pp(p,0), Data_all$Current$X[1], Data_all$Current$N[1])
post.theta_C.comb <- hx.rand(p.cc.eb.pp,0.2*niter)
post.delta.comb <- post.theta_T - post.theta_C.comb


# Proposed methods using Horseshoe prior
fit_HS     <- HS_stan$sample(data=ANA, seed = 954, init = init_list_h, chains = 4, parallel_chains = 4, iter_warmup = 4000, iter_sampling = 6000, adapt_delta = 0.999, refresh = 0, show_messages = F)


# Method, Mean, SD, LCI, UCI
Result_ALL <- matrix(0, nrow = 8, ncol = 6)
Result_ALL[1,]  <- c(1, mean(fit_Current$draws("g")), sd(fit_Current$draws("g")), quantile(fit_Current$draws("g"),0.025,names=F),quantile(fit_Current$draws("g"),0.975,names=F), max(fit_Current$summary()$rhat)  )
Result_ALL[2,]  <- c(2, mean(fit_Pooled$draws("g")), sd(fit_Pooled$draws("g")), quantile(fit_Pooled$draws("g"),0.025,names=F),quantile(fit_Pooled$draws("g"),0.975,names=F), max(fit_Pooled$summary()$rhat)  )
Result_ALL[3,]  <- c(3, mean(fit_EX$draws("g")), sd(fit_EX$draws("g")), quantile(fit_EX$draws("g"),0.025,names=F),quantile(fit_EX$draws("g"),0.975,names=F), max(fit_EX$summary()$rhat)  )
Result_ALL[4,]  <- c(4, mean(fit_EXNEX$draws("g")), sd(fit_EXNEX$draws("g")), quantile(fit_EXNEX$draws("g"),0.025,names=F),quantile(fit_EXNEX$draws("g"),0.975,names=F), max(fit_EXNEX$summary()$rhat)  )
Result_ALL[5,]  <- c(5, mean(fit_DMPP$draws("g")), sd(fit_DMPP$draws("g")), quantile(fit_DMPP$draws("g"),0.025,names=F),quantile(fit_DMPP$draws("g"),0.975,names=F), max(fit_DMPP$summary()$rhat)  )
Result_ALL[6,]  <- c(6, mean(fit_RDMPP$draws("g")), sd(fit_RDMPP$draws("g")), quantile(fit_RDMPP$draws("g"),0.025,names=F),quantile(fit_RDMPP$draws("g"),0.975,names=F), max(fit_RDMPP$summary()$rhat)  )
Result_ALL[7,]  <- c(7, mean(post.delta.comb),  sd(post.delta.comb), quantile(post.delta.comb,0.025,names=F),quantile(post.delta.comb,0.975,names=F),1  )
Result_ALL[8,]  <- c(8, mean(fit_HS$draws("g")), sd(fit_HS$draws("g")), quantile(fit_HS$draws("g"),0.025,names=F),quantile(fit_HS$draws("g"),0.975,names=F), max(fit_HS$summary()$rhat)  )

Result_ALL1 <- as.data.frame(Result_ALL,stringsAsFactors = F) 
colnames(Result_ALL1) <- c("MethodN", "Mean", "SD", "LCI", "UCI", "diag")
Result_ALL1 <- Result_ALL1 %>% dplyr::mutate(Range = (UCI-LCI)) %>%
  mutate(Method = case_when(
    MethodN == 1 ~ "Current",
    MethodN == 2 ~ "Pooled",
    MethodN == 3 ~ "EX",
    MethodN == 4 ~ "EXNEX",
    MethodN == 5 ~ "DMPP",
    MethodN == 6 ~ "RDMPP",
    MethodN == 7 ~ "EBPP",
    MethodN == 8 ~ "Horseshoe",
    TRUE ~ ""
  )) 
  

Result_ALL2 <- Result_ALL1[,c(8,2,3,4,5,7,6)]
Result_ALL2
# write.csv(Result_ALL2, "Result_BIN.csv")



############   posterior density for beta by proposed methods
result_HS <- fit_HS$output_files() %>% rstan::read_stan_csv()
plot100 <- stan_dens(result_HS, pars = c("beta"))$data
plot100_1 <- plot100 %>% filter(parameter %in% c("beta[1]","beta[2]","beta[3]","beta[4]")) %>% 
  mutate(color = case_when(
    parameter == "beta[1]" ~ "beta1",
    parameter == "beta[2]" ~ "beta2",
    parameter == "beta[3]" ~ "beta3",
    parameter == "beta[4]" ~ "beta4",
    TRUE ~ " "
  )
)

p100 <- ggplot(data = plot100_1, aes(x = value, linetype = color, color = color)) +
  geom_line(stat = "density", size = 1.1) + 
  scale_colour_manual(
    "beta", 
    values = c("#000000", "#616161", "#000000", "#616161"),
    labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]), expression(beta[4]))
  ) +
  scale_linetype_manual(
    name = "beta",
    values = c("solid","dashed","longdash", "dotted"),
    labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]), expression(beta[4]))
  ) +
  xlab(expression(beta)) +
  ylab("Density") +
  scale_x_continuous(breaks = seq(-2,4,by=2),limits = c(-2,4)) + 
  scale_y_continuous(breaks = seq(0,2,by=1),limits = c(0,2.5)) + 
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
    legend.key.width = unit(2.5,"cm"),
    legend.position = c(0.8,0.8)
  )
p100
# ggsave("Beta_BIN.tiff", p100, dpi = 300,width = 9, height = 6 ,bg = "transparent")

