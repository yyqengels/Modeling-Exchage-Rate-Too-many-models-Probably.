library("parallel");library("foreach");library("doParallel")

#trData: train and validation data
#tsData: test data
eGARCH_Model_selection <- function(input,error_dist, order1, order2){
  library("xts");library(lubridate);library(tseries);
  library(rugarch);library(fGarch);library(forecast);
  library(stochvol)
  setwd('C:/Users/Acer S3/Dropbox/621 hw6')
  if (input == "CUR-GBP"){
    ############
    #Input first series data
    ############
    FX1 <- as.xts(read.zoo("CUR-GBP.csv", 
                           index.column=list(1),header = TRUE, 
                           sep = ",", format = "%Y-%m-%d"))
    
    FX1 <- FX1[-c(dim(FX1)[1],dim(FX1)[1]-1),]
    FX1 <- FX1[wday(FX1[,1]) %in% c(2:6),]
    logFX1 <- logret(FX1, demean = TRUE)
    
    First1 <- FX1["2005-01-01/2015-12-31"]
    Second1 <- FX1["2016-01-01/2016-03-30"]
    trData <- logFX1["2005-01-01/2015-12-31"]
    tsData <- logFX1["2016-01-01/2016-03-30"]
  }else if(input == "CUR-EUR"){
    ############
    #Input second series data
    ############
    FX2 <- as.xts(read.zoo("CUR-EUR.csv", 
                           index.column=list(1),header = TRUE, 
                           sep = ",", format = "%Y-%m-%d"))
    
    #FX2 <- FX2[-c(dim(FX2)[1],dim(FX2)[1]-1),]
    FX2 <- FX2[wday(FX2[,1]) %in% c(2:6),]
    logFX2 <- logret(FX2, demean = TRUE)
    
    trData <- logFX2["2005-01-01/2015-12-31"]
    tsData <- logFX2["2016-01-01/2016-03-30"]
  }
  
  cover_egarch90 <- matrix(NA,2,2);cover_egarch80 <- matrix(NA,2,2)
  cover_egarch50 <- matrix(NA,2,2)
  tr_val_data_len <- length(trData)
  I <- order1
  J <- order2
  
  
  
  ########
  #Confidence interval for normally distributed error
  ########
  CI_sgarch <- function(TrueValue, mu, sigma, level){
    k <- length(TrueValue)
    CI <- matrix(NA,k,2)
    cover <- 0
    for (i in 1:k){
      true_value <- TrueValue[i]
      mu_temp <- mu[i]
      sigma_temp <- sigma[i]
      upper_level <- 0.5 + level/2
      lower_level <- 0.5 - level/2
      upper <- qnorm(upper_level,mean=mu_temp,sd=sigma_temp)
      lower <- qnorm(lower_level,mean=mu_temp,sd=sigma_temp)
      CI[i,1] <- lower;CI[i,2] <- upper
      condition <- TrueValue[i]>=lower & TrueValue[i] <= upper
      if (condition) cover <- cover + 1
    }
    cover_rate <- cover/k
    list(CI = as.data.frame(CI), cover_rate = cover_rate)
  }
  
  ########
  #Confidence interval for t distributed error
  ########
  CI_t_sgarch <- function(TrueValue, mu, sigma, level, shape){
    k <- length(TrueValue)
    CI <- matrix(NA,k,2)
    cover <- 0
    for (i in 1:k){
      true_value <- TrueValue[i]
      mu_temp <- mu[i]
      sigma_temp <- sigma[i]
      upper_level <- 0.5 + level/2
      lower_level <- 0.5 - level/2
      upper <- mu[i] + sigma[i]*qt(upper_level,shape[i])
      lower <- mu[i] + sigma[i]*qt(lower_level,shape[i])
      CI[i,1] <- lower;CI[i,2] <- upper
      condition <- TrueValue[i]>=lower & TrueValue[i] <= upper
      if (condition) cover <- cover + 1
    }
    cover_rate <- cover/k
    list(CI = as.data.frame(CI), cover_rate = cover_rate)
  }
  
  
  ########################
  #Loop for model selection
  ########################
  for (i in 1:I){
    for (j in 1:J){
      if(error_dist == "norm"){
        egarch_model_temp <-ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(i,j)),distribution.model="norm",mean.model=list(armaOrder=c(0,0),include.mean=TRUE))
        mod_egarch_temp = ugarchroll(egarch_model_temp, data = trData, n.ahead = 1,refit.every = 1,
                                     n.start = round(4*tr_val_data_len/5), refit.window = "recursive",
                                     solver = "hybrid", fit.control = list(),
                                     calculate.VaR = FALSE,keep.coef = TRUE)
        sigma_egarch <- mod_egarch_temp@forecast$density$Sigma
        mu_egarch <- mod_egarch_temp@forecast$density$Mu
        package90 <- CI_sgarch(tsData, mu_sgarch,sigma_sgarch, 0.9)
        package80 <- CI_sgarch(tsData, mu_sgarch,sigma_sgarch, 0.8)
        package50 <- CI_sgarch(tsData, mu_sgarch,sigma_sgarch, 0.5)
      }else if(error_dist == "std"){
        egarch_model_temp <-ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(i,j)),distribution.model="std",mean.model=list(armaOrder=c(0,0),include.mean=TRUE))
        mod_egarch_temp = ugarchroll(egarch_model_temp, data = trData, n.ahead = 1,refit.every = 1,
                                     n.start = round(4*tr_val_data_len/5), refit.window = "recursive",
                                     solver = "hybrid", fit.control = list(),
                                     calculate.VaR = FALSE,keep.coef = TRUE)
        sigma_egarch <- mod_egarch_temp@forecast$density$Sigma
        mu_egarch <- mod_egarch_temp@forecast$density$Mu
        shape_egarch <- mod_egarch_temp@forecast$density$Shape
        package90 <- CI_t_sgarch(tsData, mu_egarch,sigma_egarch, 0.9,shape_egarch)
        package80 <- CI_t_sgarch(tsData, mu_egarch,sigma_egarch, 0.8,shape_egarch)
        package50 <- CI_t_sgarch(tsData, mu_egarch,sigma_egarch, 0.5,shape_egarch)
      }
      
      cover_egarch90[i,j] <- package90$cover_rate
      cover_egarch80[i,j] <- package80$cover_rate
      cover_egarch50[i,j] <- package50$cover_rate
    }
    final_cover_rate <- list(cover90=cover_egarch90, cover80=cover_egarch80, cover50=cover_egarch50)
  }
}




#########################

#Paralell compuating

#########################
error_disritbution <- c("norm", "std")
order <- c(2,2)

cl <- makeCluster(4)
clusterExport(cl,"eGARCH_Model_selection")

cores <- seq_along(cl)
r <- clusterApply(cl[cores], cores, function(core) {
  if (core == 1) {
    eGARCH_Model_selection(input = "CUR-GBP",error_dist=error_disritbution[1],order1 = order[1], order2 = order[2])
  } else if(core == 2){
    eGARCH_Model_selection(input = "CUR-GBP",error_dist=error_disritbution[2],order1 = order[1], order2 = order[2])
  } else if(core == 3){
    eGARCH_Model_selection(input = "CUR-EUR",error_dist=error_disritbution[1],order1 = order[1], order2 = order[2])
  } else if(core == 4){
    eGARCH_Model_selection(input = "CUR-EUR",error_dist=error_disritbution[2],order1 = order[1], order2 = order[2])
  }
})