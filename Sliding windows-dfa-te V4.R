##### SLIDING WINDOW USING DETRENDED FLUCTUATION ANALYSIS and transfer entropy IN R VERSION 
##### AUTHORS: Yosef Farzi


#### CHOOSE SIZE OF THE INITIAL SLIDING WINDOW
janela = 500
plan(multisession) # for using all cores of cpu instead of one
##### LOAD LIBRARIES #############################################
require(compiler)
enableJIT(level=3)
require(xts)
require(tseries)
require(lubridate)
require(nonlinearTseries)
require(PerformanceAnalytics)
library(dplyr)
library(tidyverse)
message(" \nTAREFA FINALIZADA!")

##### READ THE FILE  ##############################################
message(" \nENTRE COM A BASE DE DADOS")
dados <- read.table(file.choose(), header=T, sep=",")
dados=dados %>% filter(!is.na(Gold) & !is.na(SPY))
names(dados)[1]="Date"
n <- length(dados[,1])
nome <- names(dados)
message(" \nARQUIVO LIDO!")

seq <- seq(from = 1, to = length(dados), by=1) # CREATES A SEQUENCE
objeto = 2        # CREATES AN OBJECT OF SIZE 2
ini <- seq[1]                 
fim <- ini[1] + month(objeto)-1
base <- numeric() #CREATES A VECTOR TO RECEIVE THE CALCULATED FIELDS

for(a in 2:3){
  st1 <- ts(dados[,a])
  gpr<- ts(dados[,"GPR"])
  ##### CODE FOR SLIDING WINDOW USING DETRENDED FLUCTUATION ANALYSIS 
  message(" \nAPLICANDO SLIDING WINDOWS, AGUARDE...")
  print(Sys.time())
  expoente_dfa <- double()
  teste_t_b1 <- double()
  erro_padrao_b1 <- double()
  p_value_b1 <- double()
  start_date <- seq[1]
  end_date <- start_date[1] + (janela)-1
  skw <- double()
  kurt <- double()
  
  jarque <- double()
  N <- length(st1)
  dateStamp<-date()
  teGPRto=double()
  teGPRtoPval=double()
  teGPRFrom=double()
  teGPRFromPval=double()
  # COMPUTE THE SLIDING WINDOW 
  for (j in 1:(N-janela+1)){ 
    dateStamp[j]=dados$Date[j]
    sliding <- window(st1, start=start_date, end=end_date)
    slidingGPR <- window(gpr, start=start_date, end=end_date)
    sliding =sliding[sliding>0]
   
      # skewness and Kurtosis
    skw[j] <- skewness(sliding, method="moment") 
    kurt[j] <- kurtosis(sliding, method="moment")
    jar <- jarque.bera.test(sliding)
    jarque[j] <- jar$p.value
    
    
    slidingDFA =diff((sliding))
    # ESTIMATES VALUE OF DETRENDED FLUCTUATION IN R
    dfa1 <- dfa(time.series = slidingDFA, window.size.range=c(10,length(sliding)/4), 
                npoints=10, do.plot=F)
    expoente_dfa[j] <-  estimate(dfa1, do.plot=F)
    
    model <- lm(log10(dfa1$fluctuation.function)~log10(dfa1$window.sizes))  
    erro_padrao_b1[j] <- coef(summary(model))[2, "Std. Error"]
    teste_t_b1[j] <- coef(summary(model))[2, "t value"]
    p_value_b1[j] <- coef(summary(model))[2, "Pr(>|t|)"]
    
    # TRANSFER ENTROPY
    te=transfer_entropy(x = slidingGPR, y = sliding, seed = 12345,shuffles = 10,nboot=100)
    
    teGPRto[j]=coef(te)[1,"te"]
    teGPRtoPval[j]=coef(te)[1,"p-value"]
    teGPRFrom[j]=coef(te)[2,"te"]
    teGPRFromPval[j]=coef(te)[2,"p-value"]
    
    start_date <- start_date + 1
    end_date <- end_date + 1
  }
  
  ########SAVE DFA DATA 
  #skw = skewness,  kurt = kurtosis, expoente = DFA
  new_base2 <- cbind(dateStamp,skw, kurt, teGPRto,teGPRtoPval,teGPRFrom,teGPRFromPval, expoente_dfa, erro_padrao_b1, teste_t_b1, p_value_b1) 
  write.table(new_base2, file=paste0("Sliding_window","_",nome[a],".csv"), sep=",", row.names = F)
  print(Sys.time())
}

Sys.sleep(40)

Sys.sleep(40)