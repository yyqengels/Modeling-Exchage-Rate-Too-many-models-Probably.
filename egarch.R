eGARCH_Model_selection2 <- function(error_dist,order1,order2){
  library("xts");library(lubridate);library(tseries);
  library(rugarch);library(fGarch);library(forecast);
  library(stochvol)
  setwd('D:\\RICE UNIVERSITY\\621\\ASSIGNMENT\\6')
  
  
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
      
      egarch_model_temp <-ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(i,j)),distribution.model=error_dist[k],mean.model=list(armaOrder=c(0,0),include.mean=TRUE))
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
      
      cover_egarch90[i,j] <- package90$cover_rate
      cover_egarch80[i,j] <- package80$cover_rate
      cover_egarch50[i,j] <- package50$cover_rate
    }
  }
  final_cover_rate <- list(cover90=cover_egarch90, cover80=cover_egarch80, cover50=cover_egarch50)
}
