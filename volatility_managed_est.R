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

rfr <- 0.05

thin <- 5 #thinning param
B <- 10000 #how many burn in draws to throw away
R <- 50000 #how many draws to keep after burn in

#load data
getSymbols("BTC-USD",from = "2020-12-01",to = "2023-09-30")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- as.Date(rownames(BTC))
getSymbols("GME",from = "2020-12-01",to = "2023-09-30")
GME <- as.data.frame(GME)
GME$Date <- as.Date(rownames(GME))
getSymbols("DOGE-USD",from = "2020-12-01",to = "2023-09-30")
DOGE <- as.data.frame(`DOGE-USD`)
DOGE$Date <- as.Date(rownames(DOGE))
getSymbols("^GSPC",from = "2020-12-01",to = "2023-09-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
getSymbols("MRNA",from = "2020-12-01",to = "2023-09-30")
MRNA <- as.data.frame(MRNA)
MRNA$Date <- as.Date(rownames(MRNA))

S <- BTC %>% merge(GME) %>% merge(DOGE) %>% merge(SP500) %>% merge(MRNA)
T <- nrow(S) - 1

Date <- S$Date
dat <- S %>% 
  dplyr::select(Date,
                "BTC-USD.Close",
                "DOGE-USD.Close",
                "GME.Close",
                "MRNA.Close",
                "GSPC.Close") %>% 
  mutate(BTC = `BTC-USD.Close`/lag(`BTC-USD.Close`)-1,
         DOGE = `DOGE-USD.Close`/lag(`DOGE-USD.Close`)-1,
         GME = `GME.Close`/lag(`GME.Close`)-1,
         MRNA = `MRNA.Close`/lag(`MRNA.Close`)-1,
         SP500 = `GSPC.Close`/lag(`GSPC.Close`)-1)

for (i in c("BTC", "DOGE", "GME","MRNA")){
  alphas_inv <- rep(0,B)
  alphas_dir <- rep(0,B)
  # for (m in c("SVIND", "SVLD", "SVMVN", "SVMALD")){
  m <- "SVMALD"
  keeps <- readRDS(paste0("keeps_",m,"_",i,"_comb.RDS"))
  V <- keeps$v[,,1]
  assign(paste0("V_",m,"_",i),V)
  ret <- dat %>% dplyr::select(all_of(i))
  
  for (b in 1:B){
    c_star <- sd(ret[,1],na.rm=TRUE)/sd(ret[,1]/lag(V[b,]) - rfr/252/lag(V[b,]),na.rm=TRUE)
    # Inverse Variance Returns
    mod <- lm(c_star/lag(V[b,])*ret[,1] + (1 - c_star/lag(V[b,]))*rfr/252 ~ ret[,1])
    alphas_inv[b] <- rnorm(1,coef(mod)[1],sqrt(vcov(mod)[1,1]))
    
    c_star <- sd(ret[,1],na.rm=TRUE)/sd(ret[,1]*lag(V[b,]) - rfr/252*lag(V[b,]),na.rm=TRUE)
    # Inverse Variance Returns
    mod <- lm(c_star*lag(V[b,])*ret[,1] + (1 - c_star*lag(V[b,]))*rfr/252 ~ ret[,1])
    alphas_dir[b] <- rnorm(1,coef(mod)[1],sqrt(vcov(mod)[1,1]))
  }
  assign(paste0("alphas_inv_",i),alphas_inv)
  assign(paste0("alphas_dir_",i),alphas_dir)
  
  if (i == "MRNA"){
    V <- keeps$v[,,2]
    assign(paste0("V_",m,"_SP500"),V)
    ret <- dat %>% dplyr::select("SP500")
    for (b in 1:B){
      c_star <- sd(ret[,1],na.rm=TRUE)/sd(ret[,1]/lag(V[b,]) - rfr/252/lag(V[b,]),na.rm=TRUE)
      # Inverse Variance Returns
      mod <- lm(c_star/lag(V[b,])*ret[,1] + (1 - c_star/lag(V[b,]))*rfr/252 ~ ret[,1])
      alphas_inv[b] <- rnorm(1,coef(mod)[1],sqrt(vcov(mod)[1,1]))
      
      c_star <- sd(ret[,1],na.rm=TRUE)/sd(ret[,1]*lag(V[b,]) - rfr/252*lag(V[b,]),na.rm=TRUE)
      # Inverse Variance Returns
      mod <- lm(c_star*lag(V[b,])*ret[,1] + (1 - c_star*lag(V[b,]))*rfr/252 ~ ret[,1])
      alphas_dir[b] <- rnorm(1,coef(mod)[1],sqrt(vcov(mod)[1,1]))
    }
    assign(paste0("alphas_inv_SP500"),alphas_inv)
    assign(paste0("alphas_dir_SP500"),alphas_dir)
  }
}