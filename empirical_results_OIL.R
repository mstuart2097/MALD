rm(list = ls()) 
library(ald)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(truncnorm)
library(mvtnorm)
library(tmvtnorm)
library(Rcpp)
library(MCMCpack)
library(quantmod)
library(RcppTN)
library(xtable)
library(tibble)#rownames_to_column()
library(tidyr)
library(reshape2)
library(latex2exp)

domean<- function(x, total){
  R <- total
  if(length(dim(x)) == 2){ #if it's a data frame
    return(apply(x[1:R,], 2, median))
  }else if (length(dim(x)) > 2){
    return(apply(x[1:R,,], 2:3, function(x){(median(x))}))
  }else{
    return(median(x[1:R]))
  }
}

docreds <- function(x, total,q){
  R <- total
  if(length(dim(x)) == 2){ #if it's a data frame
    return(apply(x[1:R,], 2, quantile, q))
  #}#else if (length(dim(x)) > 2){
    #return(apply(x[1:R,,], 2:3, function(x){(median(x))}))
  }else{
    return(quantile(x[1:R],q))
  }
}

library(LaplacesDemon)

total <- 20000 #number of mcmc iterations saved after burn-in, thinning
doESS <- function(x, total){
  R <- total
  if(!is.null(dim(x))){ #if it's a data frame
    return(apply(x[1:R,], 2, ESS))
  }else{
    return(ESS(x[1:R]))
  }
}

get_qq <- function(Value,keeps, model){## Function to make the QQ Plots
  # delta <- rep(0,T)
  # V <- array(0, dim = c(T+1, 2))
  y <- x <- array(0, dim = c(T, 2))
  y_box <- array(0,dim = c(T, 2))
  #V[1,] <- apply(keeps$v[,1,], 2, mean)
  V <- apply(keeps$v, 2:3, mean)
  #overwrite J
  J <- apply(keeps$J, 2:3, mean)
  sim <- 0
  mu <- apply(keeps$mu, 2 ,mean)
  theta <- apply(keeps$theta,2, mean)
  kappa <- 1 - apply(keeps$phi, 2, mean)
  for (t in 1:T){
    print(t)
    set.seed(t + 4922035+ sim)
    # delta <- sample(c(0:3),prob=apply(keeps$lambda, 2, mean), 1)
    # 
    # set.seed(15866245 + t + sim)
    # if(model == "MVN"){
    #   B <- 1
    # }else{
    #   B <- rexp(1)
    # }
    # Sigma <- matrix(c(mean(keeps$sigma_c[,1])^2,
    #                   mean(keeps$rhoc)*mean(keeps$sigma_c[,1])*mean(keeps$sigma_c[,2]),
    #                   mean(keeps$rhoc)*mean(keeps$sigma_c[,1])*mean(keeps$sigma_c[,2]),
    #                   mean(keeps$sigma_c[,2]^2)),
    #                 nrow=2)
    Sigma <- matrix(c(V[t,1],
                      mean(keeps$rho[,1])*sqrt(prod(V[t,])),
                      mean(keeps$rho[,3])*mean(keeps$sigma_v[,1])*V[t,1],0,
                      mean(keeps$rho[,1])*sqrt(prod(V[t,])),V[t,2],0,mean(keeps$rho[,4])*mean(keeps$sigma_v[,2])*V[t,2],
                      mean(keeps$rho[,3])*mean(keeps$sigma_v[,1])*V[t,1],0,mean(keeps$sigma_v[,1])^2*V[t,1],mean(keeps$rho[,2])*prod(apply(keeps$sigma_v, 2, mean))*sqrt(prod(V[t,])),
                      0,mean(keeps$rho[,4])*mean(keeps$sigma_v[,2])*V[t,2],mean(keeps$rho[,2])*prod(apply(keeps$sigma_v, 2, mean))*sqrt(prod(V[t,])),mean(keeps$sigma_v[,2])^2*V[t,2]),nrow=4)
    
    Mu <- c(Value[t,] - mu - J[t,],
            V[t+1,] - V[t,] - kappa*(theta - V[t,]))
    
    set.seed(463468+t)
    temp <- rtmvnorm(n = 1,
                     mean = c(mu + J[t,],
                              V[t,] + kappa*(theta - V[t,])),
                     sigma = Sigma, lower=c(-Inf,-Inf, 0, 0))
    
    # V[t+1,] <- temp[c(3:4)]
    y[t,] <- temp[c(1:2)]
    y_box[t,] <- (solve(chol(Sigma)) %*% Mu)[1:2,]
    if( t+1 <= T){ x[t+1] <- 0 }
  }
  
  QQdat = data.frame(Value=Value,
                     y)
  
  p1 <- ggplot() +
    geom_point(aes(x=quantile(QQdat$Value.Value,seq(0.01,0.99,0.01)),y=quantile(QQdat$X1,seq(0.01,0.99,0.01)))) +
    geom_abline(slope=1,intercept=0) +
    #xlim(c(-15,15)) + ylim(c(-15,15)) +
    xlab("Actual Quantiles") + ylab("Simulated Quantiles") + theme_bw() + ggtitle(paste0("SV-", model))
  
  # p2 <- ggplot() +
  #   geom_point(aes(x=quantile(QQdat$SP,seq(0.01,0.99,0.01)),y=quantile(QQdat$X2,seq(0.01,0.99,0.01)))) +
  #   geom_abline(slope=1,intercept=0) +
  #   #xlim(c(-15,15)) + ylim(c(-15,15)) +
  #   xlab("Actual Quantiles") + ylab("Simulated Quantiles") + theme_bw() + ggtitle("S&P")
  # p3 <- grid.arrange(p1,p2, nrow = 1)
  return(list(graph = p1 + theme(plot.title = element_text(hjust = 0.5)),
              p.vals = c(ks.test(QQdat$Value.Value,QQdat$X1)$p.value,
                         Box.test(y_box[,1],lag=5,type="Ljung-Box")$p.value,
                         Box.test(y_box[,1],lag=10,type="Ljung-Box")$p.value)))
}


############################################################
#----LONG TIME SERIES ONLY---------------------------------
############################################################
OIL <- read.csv("Oil_Spot_Prices.csv")
OIL$Date <- as.Date(OIL$Date)
getSymbols("^GSPC",from = "1989-05-25",to = "2009-05-25")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
S <- merge(OIL,SP500)
S <- S[S$Date <= "2005-05-25",]
T <- nrow(S) - 1
y <- as.matrix(100*(log(S[-1,c("Value","GSPC.Close")]) - log(S[-nrow(S),c("Value","GSPC.Close")])))

keeps <- readRDS(paste0("keeps_long/keeps_SVMALD_OIL.rds"))
names <- lapply(keeps[c(4,6:17)], domean, total = 10) %>%unlist %>%names

keeps_medians <- array(dim = c(4, length(names)))
keeps_summary <- array(dim = c(4, length(names)))
keeps_creds_lower <- array(dim = c(4, length(names)))
keeps_creds_upper <- array(dim = c(4, length(names)))
keeps_v1 <- array(dim = c(4,dim(keeps$v)[2]))
keeps_v2 <- array(dim = c(4,dim(keeps$v)[2]))
keeps_j1 <- array(dim = c(4,dim(keeps$J)[2]))
keeps_j2 <- array(dim = c(4,dim(keeps$J)[2]))
rho_probs <- array(dim = c(4,4))
model <- rep(NA, 4)
colnames(keeps_summary) <- names

j = 0
for (m in c("IND", "LD", "MVN", "MALD")){# "SVMVN", "SVLD",
    j = j + 1
    keeps <- readRDS(paste0("keeps_long/keeps_SV",m , "_OIL.rds"))
    # keeps_v1[j,] <- apply(keeps$v[,,1], 2, mean,na.rm=TRUE) #alternative sv
    # keeps_v2[j,] <- apply(keeps$v[,,2], 2, mean,na.rm=TRUE) #sp sv
    # 
    # keeps_j1[j,] <- apply(keeps$J[,,1], 2, mean,na.rm=TRUE) #alternative sv
    # keeps_j2[j,] <- apply(keeps$J[,,2], 2, mean,na.rm=TRUE) #sp sv
    # 
    # resids <- (100*(S$Value[-1] - S$Value[-(T+1)]) - mean(keeps$mu[,1]) - keeps_j1[1,])/sqrt(keeps_v1[j,-(T+1)])
    # 
    # keeps_medians[j,] <- lapply(keeps[c(4,6:17)], domean, total = 10000) %>%unlist
    # keeps_creds_lower[j,] <- lapply(keeps[c(4,6:17)], docreds, q=.025, total = 10000) %>%unlist
    # keeps_creds_upper[j,] <- lapply(keeps[c(4,6:17)], docreds, q=.975, total = 10000) %>%unlist
    # 
    # keeps_medians[j,] <- ifelse(grepl("phi",names),1-keeps_medians[j,],keeps_medians[j,])
    # tmp <- ifelse(grepl("phi",names),1-keeps_creds_upper[j,],keeps_creds_lower[j,])
    # keeps_creds_upper[j,] <- ifelse(grepl("phi",names),1-keeps_creds_lower[j,],keeps_creds_upper[j,])
    # keeps_creds_lower[j,] <- tmp
    # keeps_v1[j,] <- apply(keeps$v[,,1], 2, mean) #alternative sv
    # keeps_v2[j,] <- apply(keeps$v[,,2], 2, mean) #sp sv
    # keeps_j1[j,] <- apply(keeps$J[,,1], 2, mean) #alternative sv
    # keeps_j2[j,] <- apply(keeps$J[,,2], 2, mean) #sp sv
    # model[j] <- m
    # keeps_summary[j,] <- paste0(format(round(as.numeric(keeps_medians[j,]),2),nsmall=2)," (",
    #                             format(round(as.numeric(keeps_creds_lower[j,]),2),nsmall=2),",",
    #                             format(round(as.numeric(keeps_creds_upper[j,]),2),nsmall=2),")")
    # # resids <- (100*(S$Value[-1] - S$Value[-(T+1)]) - mean(keeps$mu[,1]) - keeps_j1[1,])/sqrt(keeps_v1[j,-(T+1)])
    plot <- get_qq(y,keeps,m)
    assign(paste0("QQ_",m),plot)
    # rho_probs[j,] <- keeps$rho %>%apply(2, function(x){mean(x > 0)})
    
    # model[j] <-m
}

rownames(keeps_summary) <- model

#TABLE XXX POSTERIOR MEANS OF PARAMETERS------
t(keeps_summary)%>%
  xtable() %>%
  print(sanitize.text.function = function(x){x}, type = "latex", include.rownames=TRUE)

#FIGURE XXX STOCHASTIC VOLATILITY--------
keeps_v1_long <- sqrt(t(252*keeps_v1)) %>%as.data.frame()%>%
  rename("IND"=V1,"LD"=V2,"MVN"=V3,"MALD"=V4) %>%
  mutate(Price = S$Value,
         Date = S$Date) %>% #change when add more models
  gather("model","Volatility",-c(Price,Date)) %>%
  gather("series","value",-c(Date,model))


keeps_v2_long <- sqrt(t(252*keeps_v2)) %>%as.data.frame()%>%
  rename("IND"=V1,"LD"=V2,"MVN"=V3,"MALD"=V4) %>%
  mutate(Price = S$`GSPC.Close`,
         Date = S$Date) %>% #change when add more models
  gather("model","Volatility",-c(Price,Date)) %>%
  gather("series","value",-c(Date,model))

keeps_v1_long %>%
  mutate(Model = factor(model, levels = c("IND", "LD", "MVN", "MALD"), labels = c("SV-IND", "SV-LD", "SV-MVN", "SV-MALD"))) %>%
  ggplot() +
  geom_line(aes(x = Date, y = value, linetype = Model)) +
  facet_grid(series~., scales = "free_y") +
  theme_bw() +
  annotate("rect",
           xmin=as.Date("1990-01-01"),
           xmax=as.Date("1991-07-01"),
           ymin=-Inf,ymax=Inf,fill="blue",alpha=0.2)+
  annotate("rect",
           xmin=as.Date("2003-01-01"),
           xmax=as.Date("2003-07-01"),
           ymin=-Inf,ymax=Inf,fill="blue",alpha=0.2)+
  annotate("rect",
           xmin=as.Date("1998-01-01"),
           xmax=as.Date("1999-01-01"),
           ymin=-Inf,ymax=Inf,fill="red",alpha=0.2)+
  annotate("rect",
           xmin=as.Date("2001-07-01"),
           xmax=as.Date("2002-07-01"),
           ymin=-Inf,ymax=Inf,fill="red",alpha=0.2)+
  theme(plot.margin=unit(c(1, 8, 1, 1), 'lines'),
    legend.position = c(1.15, 0.25),
        legend.background = element_rect(fill = "white", colour = NA))+
  scale_colour_grey() +
  labs(x = "Date", y = "Value")+ theme(text = element_text(size = 20))
ggsave("Oil_Vol.pdf", height = 8, width = 11) 

keeps_v2_long %>%
  mutate(Model = factor(model, levels = c("IND", "LD", "MVN", "MALD"), labels = c("SV-IND", "SV-LD", "SV-MVN", "SV-MALD"))) %>%
  ggplot() +
  geom_line(aes(x = Date, y = value, linetype = Model)) +
  facet_grid(series~., scales = "free_y") +
  theme_bw() +
  theme(plot.margin=unit(c(1, 8, 1, 1), 'lines'),
        legend.position = c(1.15, 0.25),
        legend.background = element_rect(fill = "white", colour = NA))+
  scale_colour_grey() +
  labs(x = "Date", y = "Value")+ theme(text = element_text(size = 20))
ggsave("SP_Vol_oil.pdf", height = 8, width = 11) 
