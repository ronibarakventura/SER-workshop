############################ load relevant libraries ###########################
library(readxl)
library(forecast)
library(ggplot2)
library(dplyr)
library(stats)
library(readr)

################################ define functions ##############################

# a function to seasonally adjust and detrend monthly time series
monthly_sa_dt <- function(ts) {
  ts_sa_dt <- mstl(msts(ts, seasonal.periods=c(3,4,6,12))) # Decompose the time series along several seasonal periods
  # ts_sa_dt %>% autoplot() # plot each of the decompositions
  ts_sa_dt <- ts_sa_dt[,ncol(ts_sa_dt)] # Select only the last column (seasonally adjusted and detrended time series)
  return(ts_sa_dt)
}


# a function to compute entropy
compute_entropy <- function(ts) {
  T <- nrow(ts) # length of the time series
  p_rows <- (ts %>% group_by_all %>% count)$n/T # probability of each outcome
  H <- -sum(p_rows*log2(p_rows)) # entropy
  return(H)
}


# a function to compute transfer entropy from Y to X
compute_transfer_entropy <- function(y,x,n_iter) {
  
  # compute observed transfer entropy
  TE <- compute_entropy(as.data.frame(cbind(x[2:length(x)], x[1:length(x)-1]))) - 
    compute_entropy(as.data.frame(x[1:length(x)-1])) - 
    compute_entropy(as.data.frame(cbind(x[2:length(x)], x[1:length(x)-1], y[1:length(y)-1]))) + 
    compute_entropy(as.data.frame(cbind(x[1:length(x)-1], y[1:length(y)-1])))

  # compute surrogate distribution
  SD <- c() # empty vector to be filled with TE values
  for (it in 1:n_iter) {
    y_p <- sample(y) # permutate y
    TE_p <- compute_entropy(as.data.frame(cbind(x[2:length(x)], x[1:length(x)-1]))) - 
      compute_entropy(as.data.frame(x[1:length(x)-1])) - 
      compute_entropy(as.data.frame(cbind(x[2:length(x)], x[1:length(x)-1], y_p[1:length(y_p)-1]))) + 
      compute_entropy(as.data.frame(cbind(x[1:length(x)-1], y_p[1:length(y_p)-1]))) # compute transfer entropy with permutated y
    SD <- append(SD, TE_p) # append the computed transfer entropy value to the vector
  }
  
  # compute 95th percentile of the surrogate distribution
  P <- as.numeric(quantile(SD,0.95))
  
  # compute the p-value of observed transfer entropy
  PV = as.numeric(1-which(abs(quantile(SD,seq(from=0,to=1,by=0.00001))-TE) == min(abs(quantile(SD,seq(from=0,to=1,by=0.00001))-TE)))[1]/100000)
  
  return(list("TE"=TE, "SD"=SD, "P"=P, "PV"=PV))
}


# a function to compute conditional transfer entropy from Y to X
compute_conditional_transfer_entropy <- function(y,x,z,n_iter) {
  
  # compute observed transfer entropy
  TE <- compute_entropy(as.data.frame(cbind(x[2:length(x)], x[1:length(x)-1], z[1:length(x)-1]))) - 
    compute_entropy(as.data.frame(cbind(x[1:length(x)-1], z[1:length(x)-1]))) - 
    compute_entropy(as.data.frame(cbind(x[2:length(x)], x[1:length(x)-1], y[1:length(x)-1], z[1:length(x)-1]))) + 
    compute_entropy(as.data.frame(cbind(x[1:length(x)-1], y[1:length(x)-1], z[1:length(x)-1])))
  
  # compute surrogate distribution
  SD <- c() # empty vector to be filled with TE values
  for (it in 1:n_iter) {
    
    y_p <- y # a copy pf Y to be permutated
    for (xo in unique(x)){
      for (zo in unique(z)){
        y_p[intersect(which(x==xo),which(z==zo))] <- sample(y_p[intersect(which(x==xo),which(z==zo))]) # permutate y in the subset that preserves x and z
      }
    }
    
    TE_p <- compute_entropy(as.data.frame(cbind(x[2:length(x)], x[1:length(x)-1]))) - 
      compute_entropy(as.data.frame(x[1:length(x)-1])) - 
      compute_entropy(as.data.frame(cbind(x[2:length(x)], x[1:length(x)-1], y_p[1:length(y_p)-1]))) + 
      compute_entropy(as.data.frame(cbind(x[1:length(x)-1], y_p[1:length(y_p)-1]))) # compute transfer entropy with permutated y
    SD <- append(SD, TE_p) # append the computed transfer entropy value to the vector
  }
  
  # compute 95th percentile of the surrogate distribution
  P <- as.numeric(quantile(SD,0.95))
  
  # compute the p-value of observed transfer entropy
  PV = as.numeric(1-which(abs(quantile(SD,seq(from=0,to=1,by=0.00001))-TE) == min(abs(quantile(SD,seq(from=0,to=1,by=0.00001))-TE)))[1]/100000)
  
  return(list("TE"=TE, "SD"=SD, "P"=P, "PV"=PV))
}


# a function to compute transfer entropy from Y to X with a delay d
compute_lagged_transfer_entropy <- function(y,x,d,n_iter) {
  
  # compute observed transfer entropy
  TE <- compute_entropy(as.data.frame(cbind(x[(d+2):length(x)], x[(d+1):(length(x)-1)]))) - 
    compute_entropy(as.data.frame(x[(d+1):(length(x)-1)])) - 
    compute_entropy(as.data.frame(cbind(x[(d+2):length(x)], x[(d+1):(length(x)-1)], y[1:(length(y)-1-d)]))) + 
    compute_entropy(as.data.frame(cbind(x[(d+1):(length(x)-1)], y[1:(length(y)-1-d)])))
  
  # compute surrogate distribution
  SD <- c() # empty vector to be filled with TE values
  for (it in 1:n_iter) {
    y_p <- sample(y) # permutate y
    TE_p <- compute_entropy(as.data.frame(cbind(x[(d+2):length(x)], x[(d+1):(length(x)-1)]))) - 
      compute_entropy(as.data.frame(x[(d+1):(length(x)-1)])) - 
      compute_entropy(as.data.frame(cbind(x[(d+2):length(x)], x[(d+1):(length(x)-1)], y_p[1:(length(y_p)-1-d)]))) + 
      compute_entropy(as.data.frame(cbind(x[(d+1):(length(x)-1)], y_p[1:(length(y_p)-1-d)]))) # compute transfer entropy with permutated y
    SD <- append(SD, TE_p) # append the computed transfer entropy value to the vector
  }
  
  # compute 95th percentile of the surrogate distribution
  P <- as.numeric(quantile(SD,0.95))
  
  # compute the p-value of observed transfer entropy
  PV = as.numeric(1-which(abs(quantile(SD,seq(from=0,to=1,by=0.00001))-TE) == min(abs(quantile(SD,seq(from=0,to=1,by=0.00001))-TE)))[1]/100000)
  
  return(list("TE"=TE, "SD"=SD, "P"=P, "PV"=PV))
}


################################ TOY DATA ################################ 
y = c(1,1,1,1,0,1,0,1,0,0)
x = c(1,1,0,0,1,0,1,0,0,0)
z = c(0,0,1,0,1,1,0,0,0,1)

ts <- as.data.frame(cbind(x,y))

# compute entropy
T <- nrow(ts) # 10
p_rows <- (ts %>% group_by_all %>% count)$n/T # p(0,0) = 0.2; p(1,0) = 0.2; p(0,1) = 0.4; p(1,1) = 0.2
H <- -sum(p_rows*log2(p_rows)) # 1.922 = - 0.2*(-2.322)*3 - 0.4*(-1.322) = 3*0.464 + 0.53

# compute transfer entropy from y to x
TE <- compute_entropy(as.data.frame(cbind(x[2:length(x)], x[1:length(x)-1]))) - 
  compute_entropy(as.data.frame(x[1:length(x)-1])) - 
  compute_entropy(as.data.frame(cbind(x[2:length(x)], x[1:length(x)-1], y[1:length(y)-1]))) + 
  compute_entropy(as.data.frame(cbind(x[1:length(x)-1], y[1:length(y)-1])))

# compute transfer entropy from y to x, conditioned on z
conditional_TE <- compute_entropy(as.data.frame(cbind(x[2:length(x)], x[1:length(x)-1], z[1:length(x)-1]))) - 
  compute_entropy(as.data.frame(cbind(x[1:length(x)-1], z[1:length(x)-1]))) - 
  compute_entropy(as.data.frame(cbind(x[2:length(x)], x[1:length(x)-1], y[1:length(x)-1], z[1:length(x)-1]))) + 
  compute_entropy(as.data.frame(cbind(x[1:length(x)-1], y[1:length(x)-1], z[1:length(x)-1])))


################################ REAL DATA ################################ 
setwd("/Users/ronibarakventura/Desktop")
nature <- read_csv("nature.csv")
nature$time <- time(ts(nature$Background_checks, start=c(1999, 1), end=c(2017, 12), frequency=12))

# remove seasons and trends
nature$Media_output_on_firearm_laws_and_regulations_sadt <- monthly_sa_dt(nature$Media_output_on_firearm_laws_and_regulations)
nature$Media_output_on_shootings_sadt <- monthly_sa_dt(nature$Media_output_on_shootings)
nature$Background_checks_sadt <- monthly_sa_dt(nature$Background_checks)

# symbolize time series
nature$Background_checks_symb <- c(0,as.numeric(diff(nature$Background_checks_sadt)>0))
nature$Media_output_on_shootings_symb <- c(0,as.numeric(diff(nature$Media_output_on_shootings_sadt)>0))
nature$Media_output_on_firearm_laws_and_regulations_symb <- c(0,as.numeric(diff(nature$Media_output_on_firearm_laws_and_regulations_sadt)>0))
nature$Mass_shootings_symb <- as.numeric(nature$Mass_shooting>0)

# y causes x (media output on firearm laws and regulations causes background checks)
y1 <- nature$Media_output_on_firearm_laws_and_regulations_symb
y2 <- nature$Media_output_on_shootings_symb
y3 <- nature$Mass_shooting
x <- nature$Background_checks_symb

# compute transfer entropy
compute_transfer_entropy(y1,x,10)$TE
compute_transfer_entropy(y2,x,10)$TE
compute_transfer_entropy(y3,x,10)$TE

# compute transfer entropy with permutation tests
TE_y1_x <- compute_transfer_entropy(y1,x,10000)
TE_y2_x <- compute_transfer_entropy(y2,x,10000)
TE_y3_x <- compute_transfer_entropy(y3,x,10000)

# compute conditional transfer entropy with permutation tests
cTE_y1_x_y2 <- compute_conditional_transfer_entropy(y1,x,y2,10000) 
cTE_y1_x_y3 <- compute_conditional_transfer_entropy(y1,x,y3,10000) 
cTE_y2_x_y1 <- compute_conditional_transfer_entropy(y2,x,y1,10000) 
cTE_y2_x_y3 <- compute_conditional_transfer_entropy(y2,x,y3,10000) 
cTE_y3_x_y1 <- compute_conditional_transfer_entropy(y3,x,y1,10000) 
cTE_y3_x_y2 <- compute_conditional_transfer_entropy(y3,x,y2,10000) 

# compute transfer entropy with delays
delay_results <- as.data.frame(cbind(seq(0, 12),rep(0,13),rep(0,13)))
colnames(delay_results) <- c("delay","te","pv")
for (d in 0:12){
  d_te <- compute_lagged_transfer_entropy(y1,x,d,1000)
  delay_results$te[d+1] <- d_te$TE
  delay_results$pv[d+1] <- d_te$PV
}



################################ PLOT CODES ################################

# line plot for time series
ggplot(nature,aes(x=time, y=Mass_shootings_symb)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 1)) +
  xlab("Time") +
  ylab("Mass shootings") +
  theme(text=element_text(size=12,  family="Arial", face="bold"))+
  theme_classic()

# Surrogate distribution plot
ggplot() + 
  aes(TE_y1_x$SD)+ 
  scale_y_continuous(limits = c(0, 750), breaks = seq(0, 750, 250),expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01),expand = c(0, 0)) +
  geom_histogram(binwidth=0.001, colour="black", fill="white") +
  geom_vline(xintercept = TE_y1_x$P, colour="blue", linetype = "dashed",size=1) +
  geom_vline(xintercept = TE_y1_x$TE, colour="black", linetype = "solid",size=1) +
  xlab("Transfer entropy (bits)") +
  ylab("Count") +
  theme(text=element_text(size=12,  family="Arial", face="bold"))+
  theme_classic()

# delay analysis plot
ggplot(delay_results, aes(x=delay, y=te, colour = pv >0.05)) +
  geom_point(size=6)+
  scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12,1)) +
  scale_colour_manual(name = 'pv > 0.05', values = setNames(c('grey','black'),c(T, F))) +
  labs(x="Delay", y = "Transer entropy (bits)")+
  theme(text=element_text(size=12,  family="Arial", face="bold"))+
  theme_classic()  
