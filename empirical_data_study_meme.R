rm(list = ls())  ## DON'T FORGET TO SET WD!!
library(ald)
library(ggplot2)
library(dplyr)
library(truncnorm)
library(mvtnorm)
library(tmvtnorm)
library(Rcpp)
library(MCMCpack)
library(quantmod)
library(RcppTN)
library(xtable)


thin <- 5 #thinning param
B <- 10000 #how many burn in draws to throw away
R <- 100000 #how many draws to keep after burn in
n_chns <- 1 #how many chains to run
#load data
BTC <- read.csv("BTC_USD_2020-07-01_2021-06-30-CoinDesk.csv")
BTC <- BTC %>% dplyr::select(!Currency)
names(BTC) <- c("Date","BTC-USD.Close","BTC-USD.Open","BTC-USD.High","BTC-USD.Low")
BTC$Date <- as.Date(BTC$Date)
getSymbols("GME",from = "2020-10-01",to = "2021-06-30")
GME <- as.data.frame(GME)
GME$Date <- as.Date(rownames(GME))
DOGE <- read.csv("DOGE_USD_2020-07-01_2021-06-30-CoinDesk.csv")
DOGE <- DOGE %>% dplyr::select(!Currency)
DOGE$Date <- as.Date(DOGE$Date)
names(DOGE) <- c("Date","DOGE-USD.Close","DOGE-USD.Open","DOGE-USD.High","DOGE-USD.Low")
getSymbols("^GSPC",from = "2020-10-01",to = "2021-06-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
getSymbols("MRNA",from = "2020-10-01",to = "2021-06-30")
MRNA <- as.data.frame(MRNA)
MRNA$Date <- as.Date(rownames(MRNA))

S <- BTC %>% merge(GME) %>% merge(DOGE) %>% merge(SP500) %>% merge(MRNA)
T <- nrow(S) - 1


######FIGURE 1
Date <- S$Date[-1]
cbind(100*(log(S[-1,c("GSPC.Close", "GME.Close", "BTC-USD.Close", "AMC.Close", "DOGE-USD.Close") ])-
      log(S[-nrow(S),c("GSPC.Close", "GME.Close", "BTC-USD.Close", "AMC.Close", "DOGE-USD.Close") ])),
  Date) %>%
  melt(id.vars = c("Date")) %>%
  mutate(variable = factor(variable, levels = c("GSPC.Close", "GME.Close", "BTC-USD.Close", "AMC.Close", "DOGE-USD.Close"),
                           labels = c("S&P", "GME", "BTC", "AMC", "DOGE"))) %>%
  ggplot() +
  geom_line(aes(x = Date, y = value)) +
  facet_grid(variable~., scales = "free_y") +
  theme_bw() +
  scale_x_date(breaks = "month")

ggsave("data_plot_meme.pdf", height = 10, width = 8)

###Data frame of model parameters
models <- data.frame(exp_jumps =  c(FALSE,   TRUE,  FALSE,     FALSE),
                     norm_jumps = c(FALSE,   FALSE, TRUE,      FALSE),
                     ind =        c(FALSE,   FALSE, FALSE,     TRUE),
                     model =      c("SVMALD", "SVLD", "SVMVN", "SVIND"))


#################################################### 
# ALL MODELS ---------- GameStop
#################################################### 
for (k in 1:nrow(models)){
  print(paste0("----- > Starting ", models$model[k], " model < -------"))
use_starting_values <- FALSE
sourceCpp("pgas_2d_threshold.cpp") #C++ updates
# #2-D MODEL MCMCb        cfv09
y <- as.matrix(100*(log(S[-1,c("GME.Close","GSPC.Close")]) - log(S[-nrow(S),c("GME.Close","GSPC.Close")])))
yprim <- array(0,dim=dim(y))
exp_jumps <- models$exp_jumps[k]
norm_jumps <- models$norm_jumps[k]
ind <- models$ind[k]
source("run_mcmc_2d_threshold.R") #R+B iterations of pgas.R and pgas.cpp updates
saveRDS(keeps,paste0("keeps_",models$model[k] ,"_GME.rds"))
}

#################################################### 
# ALL MODELS ---------- Dogecoin
#################################################### 
for (k in 1:nrow(models)){
  print(paste0("----- > Starting ", models$model[k], " model < -------"))
use_starting_values <- FALSE
sourceCpp("pgas_2d.cpp") #C++ updates
# #2-D MODEL MCMCb        cfv09
y <- as.matrix(100*(log(S[-1,c("DOGE-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("DOGE-USD.Close","GSPC.Close")])))
yprim <- array(0,dim=dim(y))
exp_jumps <- models$exp_jumps[k]
norm_jumps <- models$norm_jumps[k]
ind <- models$ind[k]
source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
saveRDS(keeps,paste0("keeps_",models$model[k] ,"_DOGE.rds"))
}

#################################################### 
# ALL MODELS ---------- BTC
#################################################### 
for (k in 1:nrow(models)){
  print(paste0("----- > Starting ", models$model[k], " model < -------"))
  use_starting_values <- FALSE
  sourceCpp("pgas_2d.cpp") #C++ updates
  # #2-D MODEL MCMCb        cfv09
  y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
  yprim <- array(0,dim=dim(y))
  exp_jumps <- models$exp_jumps[k]
  norm_jumps <- models$norm_jumps[k]
  ind <- models$ind[k]
  source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
  saveRDS(keeps,paste0("keeps_",models$model[k] ,"_BTC.rds"))
}

#################################################### 
# ALL MODELS ---------- Moderna
#################################################### 
for (k in 1:nrow(models)){
  print(paste0("----- > Starting ", models$model[k], " model < -------"))
  use_starting_values <- FALSE
  sourceCpp("pgas_2d.cpp") #C++ updates
  # #2-D MODEL MCMCb        cfv09
  y <- as.matrix(100*(log(S[-1,c("MRNA.Close","GSPC.Close")]) - log(S[-nrow(S),c("MRNA.Close","GSPC.Close")])))
  yprim <- array(0,dim=dim(y))
  exp_jumps <- models$exp_jumps[k]
  norm_jumps <- models$norm_jumps[k]
  ind <- models$ind[k]
  source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
  saveRDS(keeps,paste0("keeps_",models$model[k] ,"_MRNA.rds"))
}


#ALL OF THE ABOVE WAS FOR THE 'SHORT TIME PERIOD' ANALYSIS
