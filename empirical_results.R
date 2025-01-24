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

# Function to calculate the median of all posterior draws for model parameters
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

# Function to calculate q% credible interval of all posterior draws for model parameters
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

# Function to calculate effective sample size for model parameters as a convergence tactic
doESS <- function(x, total){
  R <- total
  if(!is.null(dim(x))){ #if it's a data frame
    return(apply(x[1:R,], 2, ESS))
  }else{
    return(ESS(x[1:R]))
  }
}

# Function to get QQ plots for a given model
get_qq <- function(Value,keeps, model){## Function to make the QQ Plots
  #browser()
  # delta <- rep(0,T)
  # V <- array(0, dim = c(T+1, 2))
  # Fixing volatility and jump size as the posterior means
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
    geom_point(aes(x=quantile(QQdat$Value.1,seq(0.01,0.99,0.01)),y=quantile(QQdat$X1,seq(0.01,0.99,0.01)))) +
    geom_abline(slope=1,intercept=0) +
    #xlim(c(-15,15)) + ylim(c(-15,15)) +
    xlab("Actual Quantiles") + ylab("Simulated Quantiles") + theme_bw() + ggtitle(paste0("SV-", model))
  
  # Also returns p-values from a KS test for normality and
  # p-values from Ljung-Box test
  return(list(graph = p1 + theme(plot.title = element_text(hjust = 0.5)),
              p.vals = c(ks.test(QQdat$Value.1,QQdat$X1)$p.value,
                         Box.test(y_box[,1],lag=5,type="Ljung-Box")$p.value,
                         Box.test(y_box[,1],lag=10,type="Ljung-Box")$p.value)))
}


############################################################
#----SHORT TIME SERIES ONLY---------------------------------
#ALL CRYPTO, MEME STOCKS
############################################################
#load data
thin <- 5 #thinning param
B <- 10000 #how many burn in draws to throw away
R <- 50000 #how many draws to keep after burn in
n_chns <- 1 #how many chains to run
#load data
getSymbols("BTC-USD",from = "2020-10-07",to = "2021-12-31")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- as.Date(rownames(BTC))
getSymbols("GME",from = "2020-10-07",to = "2021-12-31")
GME <- as.data.frame(GME)
GME$Date <- as.Date(rownames(GME))
getSymbols("DOGE-USD",from = "2020-10-07",to = "2021-12-31")
DOGE <- as.data.frame(`DOGE-USD`)
DOGE$Date <- as.Date(rownames(DOGE))
getSymbols("^GSPC",from = "2020-10-07",to = "2021-12-31")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
getSymbols("MRNA",from = "2020-10-07",to = "2021-12-31")
MRNA <- as.data.frame(MRNA)
MRNA$Date <- as.Date(rownames(MRNA))

S <- BTC %>% merge(GME) %>% merge(DOGE) %>% merge(SP500) %>% merge(MRNA)
T <- nrow(S) - 1
Date <- S$Date

keeps <- readRDS(paste0("keeps_short/keeps_","SVMALD" ,"_","MRNA", ".rds"))
names <- lapply(keeps[c(4,6:17)], domean, total = 10) %>%unlist %>%names

# Code for obtaining plots for our meme stock study
for (i in c("MRNA","DOGE","GME","BTC")){
  # Getting correct data frame for QQ plots
  if (i == "GME"){
    S <- GME %>% merge(SP500)
    T <- nrow(S) - 1
  } else if (i == "DOGE"){
    S <- DOGE %>% merge(SP500)
    T <- nrow(S) - 1
  } else if (i == "BTC"){
    S <- BTC %>% merge(SP500)
    T <- nrow(S) - 1
  } else {
    S <- MRNA %>% merge(SP500)
    T <- nrow(S) - 1
  }
  j = 0
  # For each of the four competing models, calculate medians and credible intervals
  # for all of the model parameters
  for (m in c("SVIND", "SVLD", "SVMVN", "SVMALD")){# "SVMVN", "SVLD",
    j = j + 1
    keeps <- readRDS(paste0("keeps_short/keeps_",m,"_",i,"_comb.rds"))
    if (j == 1){
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
    }
    # keeps_v1[j,] <- apply(keeps$v[,,1], 2, mean,na.rm=TRUE) #alternative sv
    # keeps_v2[j,] <- apply(keeps$v[,,2], 2, mean,na.rm=TRUE) #sp sv
    # 
    # keeps_j1[j,] <- apply(keeps$J[,,1], 2, mean,na.rm=TRUE) #alternative sv
    # keeps_j2[j,] <- apply(keeps$J[,,2], 2, mean,na.rm=TRUE) #sp sv
    keeps_medians[j,] <- lapply(keeps[c(4,6:17)], domean, total = 10000) %>%unlist
    keeps_creds_lower[j,] <- lapply(keeps[c(4,6:17)], docreds, q=.025, total = 10000) %>%unlist
    keeps_creds_upper[j,] <- lapply(keeps[c(4,6:17)], docreds, q=.975, total = 10000) %>%unlist
    
    keeps_medians[j,] <- ifelse(grepl("phi",names),1-keeps_medians[j,],keeps_medians[j,])
    tmp <- ifelse(grepl("phi",names),1-keeps_creds_upper[j,],keeps_creds_lower[j,])
    keeps_creds_upper[j,] <- ifelse(grepl("phi",names),1-keeps_creds_lower[j,],keeps_creds_upper[j,])
    keeps_creds_lower[j,] <- tmp
    keeps_v1[j,] <- apply(keeps$v[,,1], 2, mean) #alternative sv
    keeps_v2[j,] <- apply(keeps$v[,,2], 2, mean) #sp sv
    keeps_j1[j,] <- apply(keeps$J[,,1], 2, mean) #alternative sv
    keeps_j2[j,] <- apply(keeps$J[,,2], 2, mean) #sp sv
    model[j] <- m
    keeps_summary[j,] <- paste0(format(round(as.numeric(keeps_medians[j,]),2),nsmall=2)," (",
                                format(round(as.numeric(keeps_creds_lower[j,]),2),nsmall=2),",",
                                format(round(as.numeric(keeps_creds_upper[j,]),2),nsmall=2),")")
    # rho_probs provides a table of posterior probabilities of an inverse leverage effect
    rho_probs[j,] <- keeps$rho %>%apply(2, function(x){mean(x > 0)})
    model[j] <-m
  }
  
  keeps_summary <- cbind(keeps_summary,rho_probs[,3],rho_probs[,4])
  keeps_summary[keeps_summary == 0] <- "$>0.9999$"
  rownames(keeps_summary) <- gsub("SV","SV-",model)
  colnames(keeps_summary) <- c("$\\lambda_{y1}$","$\\lambda_{y2}$","$\\lambda_{c}$","$\\lambda_0$",
                               "$\\sigma_{v1}$","$\\sigma_{v2}$","$\\nu_{y1c}$","$\\nu_{y2c}$",
                               "$\\rho_c$","$\\mu_{y1c}$","$\\mu_{y2c}$","$\\nu_{y1}$","$\\mu_{y1}$",
                               "$\\nu_{y2}$","$\\mu_{y2}$","$\\kappa_1$","$\\kappa_2$",
                               "$\\theta_1$","$\\theta_2$","$\\mu_1$","$\\mu_2$",
                               "$\\rho_y$","$\\rho_v$","$\\rho_1$","$\\rho_2$",
                               "$P(\\rho_1 > 0 | \\boldsymbol{y})$",
                               "$P(\\rho_2 < 0 | \\boldsymbol{y})$"
  )
  keeps_summary <- keeps_summary[,c("$\\lambda_{y1}$","$\\lambda_{y2}$","$\\lambda_{c}$","$\\lambda_0$",
                                    "$\\mu_1$","$\\mu_2$","$\\theta_1$","$\\theta_2$",
                                    "$\\kappa_1$","$\\kappa_2$","$\\sigma_{v1}$","$\\sigma_{v2}$",
                                    "$\\rho_y$","$\\rho_v$","$\\rho_1$","$\\rho_2$",
                                    "$\\mu_{y1c}$","$\\mu_{y2c}$","$\\nu_{y1c}$","$\\nu_{y2c}$",
                                    "$\\rho_c$","$\\mu_{y1}$","$\\nu_{y1}$","$\\mu_{y2}$",
                                    "$\\nu_{y2}$",
                                    "$P(\\rho_1 > 0 | \\boldsymbol{y})$",
                                    "$P(\\rho_2 < 0 | \\boldsymbol{y})$"
  )]
  keeps_summary[keeps_summary == " 0.00 ( 0.00, 0.00)"] <- "-"
  keeps_summary[keeps_summary == " 5.00 ( 5.00, 5.00)"] <- "-"
  
  
  #TABLE XXX POSTERIOR MEANS OF PARAMETERS------
  t(keeps_summary) %>%
    xtable() %>%
    print(sanitize.text.function = function(x){x}, 
          type = "latex", 
          include.rownames=TRUE,
          include.colnames=TRUE,
          hline.after = c(-1,0,4,16,21,25,27))
  
  # #FIGURE XXX STOCHASTIC VOLATILITY--------
  # Plot the price changes along with stochastic volatility estimates for
  # all four competing models in a wrapped plot for Bitcoin only
  # all other meme stock plots are qualitatively similar, so did not create
  if(i == "BTC"){
    sqrt(t(252*keeps_v1)) %>%as.data.frame()%>%
      rename("SV-IND"=V1,"SV-LD"=V2,"SV-MVN"=V3,"SV-MALD"=V4) %>%
      mutate(Price = S$`BTC-USD.Close`,
             Date = S$Date) %>% #change when add more models
      pivot_longer(-c(Price,Date),names_to="Model",values_to="Volatility") %>%
      pivot_longer(-c(Date,Model),names_to="series",values_to="value",) %>%
      ggplot() +
      geom_line(aes(x = Date, y = value, linetype = Model)) +
      facet_grid(series~., scales = "free_y") +
      theme_bw() +
      annotate("rect",
               xmin=as.Date("2017-11-06"),
               xmax=as.Date("2018-01-01"),
               ymin=-Inf,ymax=Inf,fill="blue",alpha=0.2)+
      annotate("rect",
               xmin=as.Date("2018-01-01"),
               xmax=as.Date("2018-07-01"),
               ymin=-Inf,ymax=Inf,fill="red",alpha=0.2)+
      annotate("rect",
               xmin=as.Date("2020-09-01"),
               xmax=as.Date("2020-12-31"),
               ymin=-Inf,ymax=Inf,fill="blue",alpha=0.2)+
      annotate("rect",
               xmin=as.Date("2021-05-01"),
               xmax=as.Date("2021-06-30"),
               ymin=-Inf,ymax=Inf,fill="red",alpha=0.2)+
      theme(plot.margin=unit(c(1, 8, 1, 1), 'lines'),
            legend.position = c(1.15, 0.25),
            legend.background = element_rect(fill = "white", colour = NA))+
      scale_colour_grey() +
      labs(x = "Date", y = "Value")+ theme(text = element_text(size = 20))
    ggsave("BTC_Vol.pdf", height = 8, width = 11)

    keeps_j1_long <- keeps_j1 %>%as.data.frame()%>%
      mutate(model = model) %>%
      melt(id.vars = c("model")) %>%
      mutate(Date = rep(S$Date[-1], each = 4), #change when add more models
             series = "BTC")

    keeps_j2_long <- keeps_j2 %>%as.data.frame()%>%
      mutate(model = model) %>%
      melt(id.vars = c("model")) %>%
      mutate(Date = rep(S$Date[-1], each = 4), #change when add more models
             series = "S&P 500")

    # Plot the posterior estimates of the jump sizes for the four competing models
    # for both the meme stock and the S&P 500
    rbind(keeps_j1_long,keeps_j2_long) %>%
      mutate(Model = factor(model, levels = c("SVIND", "SVLD", "SVMVN", "SVMALD"), labels = c("SV-IND", "SV-LD", "SV-MVN", "SV-MALD"))) %>%
      ggplot() +
      geom_line(aes(x = Date, y = value)) +
      facet_grid(series~Model, scales = "free_y") +
      theme_bw() +
      scale_colour_grey() +
      labs(x = "Date", y = "Jump size")+ theme(text = element_text(size = 20))
    ggsave("jump_sizes.pdf", height = 10, width = 14)

  }

  # keeps_v2_long <- sqrt(t(252*keeps_v2)) %>%as.data.frame()%>%
  #   rename("IND"=V1,"LD"=V2,"MVN"=V3,"MALD"=V4) %>%
  #   mutate(Price = S$`GSPC.Close`,
  #          Date = S$Date) %>% #change when add more models
  #   gather("model","Volatility",-c(Price,Date)) %>%
  #   gather("series","value",-c(Date,model))
  # 
  # keeps_v1_long %>%
  #   mutate(Model = factor(model, levels = c("IND", "LD", "MVN", "MALD"), labels = c("SV-IND", "SV-LD", "SV-MVN", "SV-MALD"))) %>%
  #   ggplot() +
  #   geom_line(aes(x = Date, y = value, linetype = Model)) +
  #   facet_grid(series~., scales = "free_y") +
  #   theme_bw() +
  #   annotate("rect",
  #            xmin=as.Date("1990-01-01"),
  #            xmax=as.Date("1991-07-01"),
  #            ymin=-Inf,ymax=Inf,fill="blue",alpha=0.2)+
  #   annotate("rect",
  #            xmin=as.Date("2003-01-01"),
  #            xmax=as.Date("2003-07-01"),
  #            ymin=-Inf,ymax=Inf,fill="blue",alpha=0.2)+
  #   annotate("rect",
  #            xmin=as.Date("1998-01-01"),
  #            xmax=as.Date("1999-01-01"),
  #            ymin=-Inf,ymax=Inf,fill="red",alpha=0.2)+
  #   annotate("rect",
  #            xmin=as.Date("2001-07-01"),
  #            xmax=as.Date("2002-07-01"),
  #            ymin=-Inf,ymax=Inf,fill="red",alpha=0.2)+
  #   theme(plot.margin=unit(c(1, 8, 1, 1), 'lines'),
  #     legend.position = c(1.15, 0.25),
  #         legend.background = element_rect(fill = "white", colour = NA))+
  #   scale_colour_grey() +
  #   labs(x = "Date", y = "Value")+ theme(text = element_text(size = 20))
  # ggsave("Oil_Vol.pdf", height = 8, width = 11)
  # 
  # keeps_v2_long %>%
  #   mutate(Model = factor(model, levels = c("IND", "LD", "MVN", "MALD"), labels = c("SV-IND", "SV-LD", "SV-MVN", "SV-MALD"))) %>%
  #   ggplot() +
  #   geom_line(aes(x = Date, y = value, linetype = Model)) +
  #   facet_grid(series~., scales = "free_y") +
  #   theme_bw() +
  #   theme(plot.margin=unit(c(1, 8, 1, 1), 'lines'),
  #         legend.position = c(1.15, 0.25),
  #         legend.background = element_rect(fill = "white", colour = NA))+
  #   scale_colour_grey() +
  #   labs(x = "Date", y = "Value")+ theme(text = element_text(size = 20))
  # ggsave("SP_Vol_oil.pdf", height = 8, width = 11)
  # 
  # keeps_v1_short <- sqrt(t(252*keeps_v1)) %>%as.data.frame()%>%
  #   rename("IND"=V1,"LD"=V2,"MVN"=V3,"MALD"=V4) %>%
  #   mutate(Date=S$Date,
  #          Series="OIL") %>%
  #   gather("Model","Value",-c(Date,Series))
  # 
  # keeps_v2_short <- sqrt(t(252*keeps_v2)) %>%as.data.frame()%>%
  #   rename("IND"=V1,"LD"=V2,"MVN"=V3,"MALD"=V4) %>%
  #   mutate(Date=S$Date,
  #          Series="S&P 500") %>%
  #   gather("Model","Value",-c(Date,Series))

  S %>%
    dplyr::select(Date,Value,`GSPC.Close`) %>%
    rename(OIL = Value,`S&P 500` = `GSPC.Close`) %>%
    gather("Series","Value",-Date) %>%
    ggplot() +
      geom_line(aes(x = Date, y = Value)) +
      facet_grid(Series~., scales = "free_y")+
      theme_bw() +
      scale_colour_grey() +
      labs(x = "Date", y = "Price")+ theme(text = element_text(size = 20))
    ggsave("Price_oil_SP.pdf", height = 8, width = 11)

  rbind(keeps_v1_short,keeps_v2_short) %>%
    ggplot() +
    geom_line(aes(x = Date, y = Value, linetype = Model)) +
    facet_grid(Series~., scales = "free_y")+
    theme_bw() +
    scale_colour_grey() +
    labs(x = "Date", y = "Annual Volatility")+ theme(text = element_text(size = 20))
  ggsave("Vol_oil_SP.pdf", height = 8, width = 11)


}

# Create four-by-four qqplots for each asset
for (i in c("MRNA","DOGE","GME","BTC")){
  if (i == "MRNA"){
    S <- MRNA %>% merge(SP500)
    T <- nrow(S) - 1
    resids <- 100*cbind((log(S$MRNA.Close[-1]) - log(S$MRNA.Close[-(T+1)])),
                        log(S$`GSPC.Close`[-1]) - log(S$`GSPC.Close`[-(T+1)]))
  } else if (i == "DOGE"){
    S <- DOGE %>% merge(SP500)
    T <- nrow(S) - 1
    resids <- 100*cbind((log(S$`DOGE-USD.Close`[-1]) - log(S$`DOGE-USD.Close`[-(T+1)])),
                        log(S$`GSPC.Close`[-1]) - log(S$`GSPC.Close`[-(T+1)]))
  } else if (i == "BTC"){
    S <- BTC %>% merge(SP500)
    T <- nrow(S) - 1
    resids <- 100*cbind((log(S$`BTC-USD.Close`[-1]) - log(S$`BTC-USD.Close`[-(T+1)])),
                        log(S$`GSPC.Close`[-1]) - log(S$`GSPC.Close`[-(T+1)]))
  } else {
    S <- GME %>% merge(SP500)
    T <- nrow(S) - 1
    resids <- 100*cbind((log(S$GME.Close[-1]) - log(S$GME.Close[-(T+1)])),
                        log(S$`GSPC.Close`[-1]) - log(S$`GSPC.Close`[-(T+1)]))
  }
  # qqplots-----------
  model <- "IND"
  keeps <- readRDS(paste0("keeps_short/keeps_SV",model,"_",i,".rds"))
  plot1 <- get_qq(resids,keeps,model)
  #ggsave(paste0("QQ_", model, ".pdf"),plot, width = 12, height = 6)
  
  model <- "LD"
  keeps <- readRDS(paste0("keeps_short/keeps_SV",model,"_",i,".rds"))
  plot2 <- get_qq(resids,keeps,model)
  #ggsave(paste0("QQ_", model, ".pdf"),plot, width = 12, height = 6)
  
  model <- "MVN"
  keeps <- readRDS(paste0("keeps_short/keeps_SV",model,"_",i,".rds"))
  plot3 <- get_qq(resids,keeps,model)
  #ggsave(paste0("QQ_", model, ".pdf"),plot, width = 12, height = 6)
  
  model <- "MALD"
  keeps <- readRDS(paste0("keeps_short/keeps_SV",model,"_",i,".rds"))
  plot4 <- get_qq(resids,keeps,model)
  
  plot5 <- grid.arrange(plot1$graph,
                        plot2$graph,
                        plot3$graph,
                        plot4$graph, nrow = 2)
  
  ggsave(paste0("QQ_",i, "_short.pdf"),plot5, width = 6, height = 6)
  assign(paste0(i,"_p_vals"),
         cbind(plot1$p.vals,plot2$p.vals,plot3$p.vals,plot4$p.vals),
  )
}

# Create trace plots for each of the model parameters for the SV-MALD model
# to detect convergence
for (i in c("GME","DOGE","MRNA","BTC")){
  m = "MALD"
  keeps <- readRDS(paste0("keeps_short/keeps_SV",m,"_",i,".rds"))
  MCMC <- data.frame(cbind(keeps$lambda[,-4],keeps$mu,keeps$theta,1-keeps$phi,
                           keeps$sigma_v,keeps$rho,keeps$xi_cw,keeps$sigma_c,keeps$rhoc,
                           keeps$xi_y1w,keeps$xi_y1eta,keeps$xi_y2w,keeps$xi_y2eta))
  colnames(MCMC) <- c("$\\lambda_{y1}$","$\\lambda_{y2}$","$\\lambda_{c}$",
                      "$\\mu_1$","$\\mu_2$","$\\theta_1$","$\\theta_2$",
                      "$\\kappa_1$","$\\kappa_2$","$\\sigma_{v1}$","$\\sigma_{v2}$",
                      "$\\rho_y$","$\\rho_v$","$\\rho_1$","$\\rho_2$",
                      "$\\mu_{y1c}$","$\\mu_{y2c}$","$\\nu_{y1c}$","$\\nu_{y2c}$",
                      "$\\rho_c$","$\\mu_{y1}$","$\\nu_{y1}$","$\\mu_{y2}$",
                      "$\\nu_{y2}$")
  
  appender <- function(string){
    TeX(string)
  }
  MCMC %>%
    mutate(Index = 1:length(keeps$rhoc)) %>%
    dplyr::select(1:8,25) %>%
    pivot_longer(cols = -Index,names_to = "Parameter",values_to = "Value") %>%
    ggplot(aes(x=Index,y=Value)) +
    geom_line() +
    facet_wrap(.~Parameter,nrow=4,scales="free",labeller=as_labeller(appender,default = label_parsed))
  ggsave(paste0("Trace_",i, "_short_1.pdf"), width = 6, height = 10)
  
  MCMC %>%
    mutate(Index = 1:length(keeps$rhoc)) %>%
    dplyr::select(9:16,25) %>%
    pivot_longer(cols = -Index,names_to = "Parameter",values_to = "Value") %>%
    ggplot(aes(x=Index,y=Value)) +
    geom_line() +
    facet_wrap(.~Parameter,nrow=4,scales="free",labeller=as_labeller(appender,default = label_parsed))
  ggsave(paste0("Trace_",i, "_short_2.pdf"), width = 6, height = 10)
  
  MCMC %>%
    mutate(Index = 1:length(keeps$rhoc)) %>%
    dplyr::select(17:24,25) %>%
    pivot_longer(cols = -Index,names_to = "Parameter",values_to = "Value") %>%
    ggplot(aes(x=Index,y=Value)) +
    geom_line() +
    facet_wrap(.~Parameter,nrow=4,scales="free",labeller=as_labeller(appender,default = label_parsed))
  ggsave(paste0("Trace_",i, "_short_3.pdf"), width = 6, height = 10)
}
