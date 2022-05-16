##----------------------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(brms)
library(lubridate)
#trout <- read.csv("./results_meta_offspring.csv")
trout <- read.csv("./results_meta_hiphop_offspring.csv")
trout$tributary <- gsub("Trout ","Trout",trout$tributary)
females <- trout[which(trout$sex == "Female" & trout$tributary == "Trout"),]
males <- trout[which(trout$sex == "Male" & trout$tributary == "Trout"),]

females$day <- yday(as.Date(females$date,"%m/%d/%Y"))
females$hiphop_offspring <- sapply(females$hiphop_offspring,function(x) ifelse(is.na(x),0,x))
males$day <- yday(as.Date(males$date,"%m/%d/%Y"))
males$hiphop_offspring <- sapply(males$hiphop_offspring,function(x) ifelse(is.na(x),0,x))
inv_logit<-function(x){ exp(x)/(1+exp(x)) }
IQR <- function(y) { quantile(y, 0.75) - quantile(y, 0.25) } 

females$folded <- ifelse(females$q > 0.5,(abs(females$q - 1))*2,(females$q*2))
males$folded <- ifelse(males$q > 0.5,(abs(males$q - 1))*2,(males$q*2))
trout_temp <- read.csv("nfs_troutdailytemp2019.csv")
trout_temp$date <- as.Date(trout_temp$date,"%Y-%m-%d")
trout_temp$yday <- yday(trout_temp$date)


females$temp <- trout_temp$mean_temp[match(females$day,trout_temp$yday)]
males$temp <- trout_temp$mean_temp[match(males$day,trout_temp$yday)]
##----------------------------------------------------------------

##----------------------------------------------------------------
#Model reproductive output first
##----------------------------------------------------------------


#frequentist model to set priors
library(pscl)
library(MuMIn)
#m2 <- zeroinfl(hiphop_offspring ~ q + scale(day) + I(scale(day)^2) + scale(length) + scale(temp)*scale(day),data=females,dist="negbin")
#summary(m2)
m3 <- zeroinfl(hiphop_offspring ~ q + scale(day) + I(scale(day)^2) + scale(length),data=females,dist="negbin")
#summary(m3)
m4 <- zeroinfl(hiphop_offspring ~ q + scale(day) + I(scale(day)^2) + scale(length) + scale(temp),data=females,dist="negbin")
#summary(m4)
m5 <- zeroinfl(hiphop_offspring ~ q + scale(day) + scale(length) + scale(temp),data=females,dist="negbin")
m6 <- zeroinfl(hiphop_offspring ~ q + scale(day) + scale(length),data=females,dist="negbin")
m7 <- zeroinfl(hiphop_offspring ~ q + scale(temp) + scale(length),data=females,dist="negbin")
m8 <- zeroinfl(hiphop_offspring ~ q + I(q^2) + scale(day) + scale(temp) + scale(length),data=females,dist="negbin")
m9 <- zeroinfl(hiphop_offspring ~ q + I(q^2) + scale(day) + I(scale(day)^2) + scale(temp) + scale(length),data=females,dist="negbin")
m10 <- zeroinfl(hiphop_offspring ~ q + I(q^2) + scale(temp) + scale(length),data=females,dist="negbin")
AICc(m3,m4,m5,m6,m7,m8,m9,m10)
#m3

#try with folded q
m2 <- zeroinfl(hiphop_offspring ~ folded + scale(day) + I(scale(day)^2) + scale(length) + scale(temp)*scale(day),data=females,dist="negbin")
#summary(m2)
m3 <- zeroinfl(hiphop_offspring ~ folded + scale(day) + I(scale(day)^2) + scale(length),data=females,dist="negbin")
#summary(m3)
m4 <- zeroinfl(hiphop_offspring ~ folded + scale(day) + I(scale(day)^2) + scale(length) + scale(temp),data=females,dist="negbin")
#summary(m4)
m5 <- zeroinfl(hiphop_offspring ~ folded + scale(day) + scale(length) + scale(temp),data=females,dist="negbin")
m6 <- zeroinfl(hiphop_offspring ~ folded + scale(day) + scale(length),data=females,dist="negbin")
m7 <- zeroinfl(hiphop_offspring ~ folded + scale(temp) + scale(length),data=females,dist="negbin")
m8 <- zeroinfl(hiphop_offspring ~ folded + I(folded^2) + scale(day) + scale(temp) + scale(length),data=females,dist="negbin")
m9 <- zeroinfl(hiphop_offspring ~ folded + I(folded^2) + scale(day) + I(scale(day)^2) + scale(temp) + scale(length),data=females,dist="negbin")
m10 <- zeroinfl(hiphop_offspring ~ folded + I(folded^2) + scale(temp) + scale(length),data=females,dist="negbin")
AICc(m2,m3,m4,m5,m6,m7,m8,m9,m10) #m9 is best by big margin


##----------------------------------------------------------------
#prior predictive simulation, start with neg. binomial part of model

#q prior
alpha_samps<-rnorm(n=500, mean= -1.31, sd= 10)
b1_samps<- rnorm(n=500, mean= 2.22, sd= 10) 

plot(NULL , 
     xlim=range(females$q),
     ylim=range(females$hiphop_offspring,na.rm=TRUE)*2,
     xlab="q", ylab="Prior num offspring",main='')
abline( h=min(females$hiphop_offspring,na.rm = TRUE), lty=2,lwd=5 ) 
abline( h=max(females$hiphop_offspring,na.rm=TRUE), lty=3,lwd=5 ) 
for (i in 1:500){
  curve(alpha_samps[i] + b1_samps[i]*(x),
        from=min(females$q), 
        to=max(females$q),  
        col=alpha("blue",0.2),add=TRUE) }

#length prior
alpha_samps<-rnorm(n=500, mean= -1.31, sd= 8)
b1_samps<- rnorm(n=500, mean= 0.189, sd= 5) 

plot(NULL , 
     xlim=range(scale(females$length)), 
     ylim=range(females$hiphop_offspring,na.rm=TRUE)*2,  
     xlab="length", ylab="Prior num offspring",main='')
abline( h=min(females$hiphop_offspring,na.rm = TRUE), lty=2,lwd=5 ) 
abline( h=max(females$hiphop_offspring,na.rm=TRUE), lty=3,lwd=5 ) 
for (i in 1:500){
  curve(alpha_samps[i] + b1_samps[i]*(x), 
        from=min(scale(females$length)), 
        to=max(scale(females$length)),  
        col=alpha("blue",0.2),add=TRUE) }

#day + day^2 prior
alpha_samps<-rnorm(n=500, mean= -1.31, sd= 8)
b1_samps<- rnorm(n=500, mean= 0.47, sd= 6) 
b2_samps <- rnorm(n=500,mean=-0.79,sd=5)

plot(NULL , 
     xlim=range(scale(females$day)),  
     ylim=range(females$hiphop_offspring,na.rm=TRUE)*2,  
     xlab="day", ylab="Prior num offspring",main='')
abline( h=min(females$hiphop_offspring,na.rm = TRUE), lty=2,lwd=5 )
abline( h=max(females$hiphop_offspring,na.rm=TRUE), lty=3,lwd=5 ) 
for (i in 1:500){
  curve(alpha_samps[i] + b1_samps[i]*(x) + b2_samps[i]*(x^2), 
        from=min(scale(females$day)), 
        to=max(scale(females$day)),   
        col=alpha("blue",0.2),add=TRUE) }

#temp prior
alpha_samps<-rnorm(n=500, mean= -1.31, sd= 8)
b1_samps<- rnorm(n=500, mean= 0.634, sd= 5) 

plot(NULL , 
     xlim=range(scale(females$temp)),  
     ylim=range(females$hiphop_offspring,na.rm=TRUE)*2,   
     xlab="temp", ylab="Prior num offspring",main='')
abline( h=min(females$hiphop_offspring,na.rm = TRUE), lty=2,lwd=5 )
abline( h=max(females$hiphop_offspring,na.rm=TRUE), lty=3,lwd=5 )
for (i in 1:500){
  curve(alpha_samps[i] + b1_samps[i]*(x), 
        from=min(scale(females$temp)), 
        to=max(scale(females$temp)),   
        col=alpha("blue",0.2),add=TRUE) }


##----------------------------------------------------------------
#prior predictive simulation for bernoulli part of model

#q
alpha_samps<-rnorm(n=500, mean= 0, sd= 1.5)
b1_samps<- rnorm(n=500, mean= 0, sd= 1.5) 

plot(NULL , 
     xlim=range(females$q), 
     ylim=c(0,1),   #<- possible limits for  a probability
     xlab="q", ylab="Prior probability of 'sterile'",main='alpha ~ N(0,0.75)  beta ~ N(0,1.5)')
for (i in 1:500){
  curve( inv_logit(alpha_samps[i] + b1_samps[i]*(x)), 
         from=min(females$q), 
         to=max(females$q),  
         col=alpha("blue",0.2),add=TRUE) }

#day
alpha_samps<-rnorm(n=500, mean= 0, sd= 1.5)
b1_samps<- rnorm(n=500, mean= 0, sd= 0.75) 
b2_samps <- rnorm(n=500,mean=0,sd=0.5)

plot(NULL , 
     xlim=range(scale(females$day)), 
     ylim=c(0,1),   
     xlab="day", ylab="Prior probability of 'sterile'",main='alpha ~ N(0,1.5)  beta ~ N(0,0.75) beta2 ~ N(0,0.5)')
for (i in 1:500){
  curve( inv_logit(alpha_samps[i] + b1_samps[i]*(x) + b2_samps[i]*(x^2)), 
         from=min(scale(females$day)), 
         to=max(scale(females$day)),   #^ extent of x in the code just above
         col=alpha("blue",0.2),add=TRUE) }
##----------------------------------------------------------------
#run brms model
#try out a few and compare WAIC
library(brms)
func <- bf(hiphop_offspring ~ q + scale(day) + I(scale(day)^2) + scale(length), zi ~ q + scale(day) + I(scale(day)^2))

get_prior(func,females,family=zero_inflated_negbinomial())

priors <- c(set_prior("normal(0.132,10)",class="Intercept"),
            set_prior("normal(1.745,10)",class="b",coef="q"),
            set_prior("normal(-0.350,6)",class="b",coef="scaleday"),
            set_prior("normal(-1.447,5)",class="b",coef="IscaledayE2"),
            set_prior("normal(0.011,5)",class="b",coef="scalelength"),
            #set_prior("normal(0.634,5)",class="b",coef="scaletemp"),
            set_prior("normal(0,1.5)",class="Intercept",dpar="zi"),
            set_prior("normal(0,1.5)",class="b",coef="q",dpar="zi"),
            set_prior("normal(0,0.75)",class="b",coef="scaleday",dpar="zi"),
            set_prior("normal(0,0.5)",class="b",coef="IscaledayE2",dpar="zi"))
            #set_prior("normal(0,1.5)",class="b",coef="scaletemp",dpar="zi"))


offspring_model <-brm(formula=func, data=females, family=zero_inflated_negbinomial(), 
                      prior = priors,      
                      control=list(adapt_delta = 0.8, max_treedepth=10), 
                      warmup=1000, iter=3000, thin=1, chains=3, cores=3)   
offspring_model
save(offspring_model,file="./offspring_model.Rds")
gc()
plot(offspring_model)
plot(offspring_model,N=10)
pp_check(offspring_model,nsamples=500,type='bars')


func2 <- bf(hiphop_offspring ~ q + I(q^2) + scale(day) + I(scale(day)^2) + scale(length), zi ~ q + I(q^2) + scale(day) + I(scale(day)^2))

get_prior(func2,females,family=zero_inflated_negbinomial())

priors2 <- c(set_prior("normal(-1.75,10)",class="Intercept"),
            set_prior("normal(2.42,10)",class="b",coef="q"),
            set_prior("normal(0.7684,6)",class="b",coef="scaleday"),
            set_prior("normal(0.223,5)",class="b",coef="IscaledayE2"),
            set_prior("normal(0,5)",class="b",coef="IqE2"),
            set_prior("normal(0.28,5)",class="b",coef="scalelength"),
            set_prior("normal(0,1.5)",class="Intercept",dpar="zi"),
            set_prior("normal(0,1.5)",class="b",coef="q",dpar="zi"),
            set_prior("normal(0,0.75)",class="b",coef="IqE2",dpar="zi"),
            set_prior("normal(0,0.75)",class="b",coef="scaleday",dpar="zi"),
            set_prior("normal(0,0.5)",class="b",coef="IscaledayE2",dpar="zi"))


offspring_model2 <-brm(formula=func2, data=females, family=zero_inflated_negbinomial(), 
                      prior = priors2,      
                      control=list(adapt_delta = 0.8, max_treedepth=10), 
                      warmup=1000, iter=3000, thin=1, chains=3, cores=3)   
offspring_model2
gc()
plot(offspring_model2)
pp_check(offspring_model2,nsamples=500,type='bars')

func3 <- bf(hiphop_offspring ~ folded + I(folded^2) + scale(day) + I(scale(day)^2) + scale(length), zi ~ folded + I(folded^2) + scale(day) + I(scale(day)^2))

get_prior(func3,females,family=zero_inflated_negbinomial())

priors3 <- c(set_prior("normal(-1.75,10)",class="Intercept"),
             set_prior("normal(2.42,10)",class="b",coef="folded"),
             set_prior("normal(0.7684,6)",class="b",coef="scaleday"),
             set_prior("normal(0.223,5)",class="b",coef="IscaledayE2"),
             set_prior("normal(0,5)",class="b",coef="IfoldedE2"),
             set_prior("normal(0.28,5)",class="b",coef="scalelength"),
             set_prior("normal(0,1.5)",class="Intercept",dpar="zi"),
             set_prior("normal(0,1.5)",class="b",coef="folded",dpar="zi"),
             set_prior("normal(0,0.75)",class="b",coef="IfoldedE2",dpar="zi"),
             set_prior("normal(0,0.75)",class="b",coef="scaleday",dpar="zi"),
             set_prior("normal(0,0.5)",class="b",coef="IscaledayE2",dpar="zi"))


offspring_model3 <-brm(formula=func3, data=females, family=zero_inflated_negbinomial(), 
                       prior = priors3,      
                       control=list(adapt_delta = 0.8, max_treedepth=10), 
                       warmup=1000, iter=3000, thin=1, chains=3, cores=3)   
offspring_model3
gc()
plot(offspring_model3)
pp_check(offspring_model3,nsamples=500,type='bars')

func4 <- bf(hiphop_offspring ~ folded + scale(day) + I(scale(day)^2) + scale(length), zi ~ folded + scale(day) + I(scale(day)^2))

get_prior(func4,females,family=zero_inflated_negbinomial())

priors4 <- c(set_prior("normal(-1.75,10)",class="Intercept"),
             set_prior("normal(2.42,10)",class="b",coef="folded"),
             set_prior("normal(0.7684,6)",class="b",coef="scaleday"),
             set_prior("normal(0.223,5)",class="b",coef="IscaledayE2"),
             set_prior("normal(0.28,5)",class="b",coef="scalelength"),
             set_prior("normal(0,1.5)",class="Intercept",dpar="zi"),
             set_prior("normal(0,1.5)",class="b",coef="folded",dpar="zi"),
             set_prior("normal(0,0.75)",class="b",coef="scaleday",dpar="zi"),
             set_prior("normal(0,0.5)",class="b",coef="IscaledayE2",dpar="zi"))


offspring_model4 <-brm(formula=func4, data=females, family=zero_inflated_negbinomial(), 
                       prior = priors4,      
                       control=list(adapt_delta = 0.8, max_treedepth=10), 
                       warmup=1000, iter=3000, thin=1, chains=3, cores=3)   
offspring_model4
gc()
plot(offspring_model4)
pp_check(offspring_model4,nsamples=500,type='bars')




waic(offspring_model,offspring_model2,offspring_model3,offspring_model4) #m1 preferred



newdata <- data.frame(q=seq(0,0.9999,by=0.111),day=rep(0,5),length=rep(0,5))
test <- predict(offspring_model)
which.max(test[,1])
test_prune <- test[which(test[,2] < 10000),]
hist(test[which(test[,2] < 1000),1],labels=TRUE,breaks=seq(0,70,by=1))
plot(test_prune[,1]~test_prune[,2])

#plot conditional effect of q
#date = 19 May 2019
#size = 455.04 mm TL
alpha_samp1 <- sample(posterior_samples(offspring_model,pars="Intercept")[,1],size=6000,replace=TRUE)
q_samp1 <- sample(posterior_samples(offspring_model,pars="q")[,1],size=6000,replace=TRUE)
#q2_samp1 <- sample(posterior_samples(offspring_model3,pars="foldedE2")[,1],size=6000,replace=TRUE)
q_seq <- seq(0,1,by=0.005)

est_upper95 <- c()
est_lower95 <- c()
est_mean <- c()
for (i in 1:length(q_seq)){
  #est_dist <- exp(alpha_samp1 + q_samp1*q_seq[i] + q2_samp1*q_seq[i])  #exp() does log link to put on scale of negbin response
  est_dist <- exp(alpha_samp1 + q_samp1*q_seq[i])
  est_upper95[i] <- quantile(est_dist,0.975)
  est_lower95[i] <- quantile(est_dist,0.025)
  est_mean[i] <- mean(est_dist)
}
#plot to test
{plot(est_upper95~q_seq,type="l",col="red")
lines(est_mean~q_seq)
lines(est_lower95~q_seq,col="blue")}
#cool -- that seems to have worked

#ggplot of q conditional effect
cond_q <- data.frame("q"=q_seq,"mean"=est_mean,"upper95"=est_upper95,"lower95"=est_lower95)
p1 <-ggplot(cond_q,aes(x=q,y=upper95)) + geom_line(aes(x=q,y=mean),size=1.3,color="blue") + theme_minimal() +
  labs(x="Female prop. YCT ancestry",y="Number of offspring") +
  geom_ribbon(aes(ymin=lower95,ymax=upper95),alpha=0.4) + scale_fill_manual(values="gray30") +
  theme(panel.grid.minor.x = element_blank(),axis.text=element_text(size=12),axis.title = element_text(size=16)) +
  geom_point(data=females,aes(x=q,y=hiphop_offspring),size=3,fill="black")

#conditional effect of q on logistic regression
alpha_samp2 <- sample(posterior_samples(offspring_model,pars="Intercept")[,2],size=6000,replace=TRUE)
q_samp2 <- sample(posterior_samples(offspring_model,pars="q")[,2],size=6000,replace=TRUE)
q_seq <- seq(0,1,by=0.005)

est_upper95 <- c()
est_lower95 <- c()
est_mean <- c()
for (i in 1:length(q_seq)){
  est_dist <- inv_logit(alpha_samp2 + q_samp2*q_seq[i])
  est_upper95[i] <- quantile(est_dist,0.975) 
  est_lower95[i] <- quantile(est_dist,0.025)
  est_mean[i] <- mean(est_dist)
}
#plot to test
{plot(est_upper95~q_seq,type="l",col="red",ylim=c(0,1))
  lines(est_mean~q_seq)
  lines(est_lower95~q_seq,col="blue")}
#cool -- that seems to have worked

#ggplot of q conditional effect
#females$offspr_bin <- gsub("[1-9]+","1",females$hiphop_offspring)
cond_q <- data.frame("q"=q_seq,"mean"=est_mean,"upper95"=est_upper95,"lower95"=est_lower95)
ggplot(cond_q,aes(x=q,y=upper95)) + geom_line(aes(x=q,y=mean),size=1.3,color="blue") + theme_minimal() +
  labs(x="Female q",y="P(Certain zero offspring sampled)") + lims(y=c(0,1)) +
  geom_ribbon(aes(ymin=lower95,ymax=upper95),alpha=0.4) + scale_fill_manual(values="gray30") +
  theme(panel.grid.minor.x = element_blank(),axis.text=element_text(size=12),axis.title = element_text(size=16))# +
  #geom_point(data=females,aes(x=q,y=offspr_bin),size=3,fill="black")

#get probability of zero offspring total
alpha_samp1 <- sample(posterior_samples(offspring_model,pars="Intercept")[,1],size=6000,replace=TRUE)
alpha_samp2 <- sample(posterior_samples(offspring_model,pars="Intercept")[,2],size=6000,replace=TRUE)
q_samp1 <- sample(posterior_samples(offspring_model,pars="q")[,1],size=6000,replace=TRUE)
q_samp2 <- sample(posterior_samples(offspring_model,pars="q")[,2],size=6000,replace=TRUE)
shape_samp <- sample(posterior_samples(offspring_model,pars="shape")[,1],size=6000,replace=TRUE)
q_seq <- seq(0,1,by=0.005)

est_upper95 <- c()
est_lower95 <- c()
est_mean <- c()
for (i in 1:length(q_seq)){
  est_dist <- inv_logit(alpha_samp2 + q_samp2*q_seq[i])
  est_dist2 <- exp(alpha_samp1 + q_samp1*q_seq[i])
  prob <- est_dist + (1-est_dist)*dnbinom(0,size=shape_samp,mu=est_dist2)
  est_upper95[i] <- quantile(prob,0.975) 
  est_lower95[i] <- quantile(prob,0.025)
  est_mean[i] <- mean(prob)
}
#plot to test
{plot(est_upper95~q_seq,type="l",col="red",ylim=c(0,1))
  lines(est_mean~q_seq)
  lines(est_lower95~q_seq,col="blue")}
#cool -- that seems to have worked

#ggplot of q conditional effect
#females$offspr_bin <- gsub("[1-9]+","1",females$hiphop_offspring)
cond_q <- data.frame("q"=q_seq,"mean"=est_mean,"upper95"=est_upper95,"lower95"=est_lower95)
p2 <- ggplot(cond_q,aes(x=q,y=upper95)) + geom_line(aes(x=q,y=mean),size=1.3,color="blue") + theme_minimal() +
  labs(x="Female prop. YCT ancestry",y="P(Zero offspring sampled)") + lims(y=c(0,1)) +
  geom_ribbon(aes(ymin=lower95,ymax=upper95),alpha=0.4) + scale_fill_manual(values="gray30") +
  theme(panel.grid.minor.x = element_blank(),axis.text=element_text(size=12),axis.title = element_text(size=16))# +
#geom_point(data=females,aes(x=q,y=offspr_bin),size=3,fill="black")

library(gridExtra)
grid.arrange(p1,p2,ncol=2)





#plot num offspring per female
library(ggplot2)
ggplot(females,aes(x=hiphop_offspring)) + geom_bar() + theme_minimal()+ labs(x="Number of sampled offspring",y="Number of females") + guides()

#plot num offspring ~ q
ggplot(females,aes(x=q,y=hiphop_offspring)) + geom_point() + theme_minimal() + labs(x="Proportion YCT ancestry (q)",y="Number of sampled offspring")




#------------------------------------------------------------
#males!
#------------------------------------------------------------

#now for males
m2 <- zeroinfl(hiphop_offspring ~ q + scale(day) + I(scale(day)^2) + scale(length) + scale(temp)*scale(day),data=males,dist="negbin")
#summary(m2)
m3 <- zeroinfl(hiphop_offspring ~ q + scale(day) + I(scale(day)^2) + scale(length),data=males,dist="negbin")
#summary(m3)
m4 <- zeroinfl(hiphop_offspring ~ q + scale(day) + I(scale(day)^2) + scale(length) + scale(temp),data=males,dist="negbin")
#summary(m4)
m5 <- zeroinfl(hiphop_offspring ~ q + scale(day) + scale(length) + scale(temp),data=males,dist="negbin")
m6 <- zeroinfl(hiphop_offspring ~ q + scale(day) + scale(length),data=males,dist="negbin")
m7 <- zeroinfl(hiphop_offspring ~ q + scale(temp) + scale(length),data=males,dist="negbin")
m8 <- zeroinfl(hiphop_offspring ~ q + I(q^2) + scale(day) + scale(temp) + scale(length),data=males,dist="negbin")
m9 <- zeroinfl(hiphop_offspring ~ q + I(q^2) + scale(day) + I(scale(day)^2) + scale(temp) + scale(length),data=males,dist="negbin")
m10 <- zeroinfl(hiphop_offspring ~ q + I(q^2) + scale(temp) + scale(length),data=males,dist="negbin")
m11 <- zeroinfl(hiphop_offspring ~ q  + scale(temp) + scale(length),data=males,dist="negbin")
m12 <- zeroinfl(hiphop_offspring ~ q + scale(length),data=males,dist="negbin")
m13 <- zeroinfl(hiphop_offspring ~ q + I(q^2) + scale(length),data=males,dist="negbin")
m14 <- zeroinfl(hiphop_offspring ~ q + I(q^2) + scale(temp),data=males,dist="negbin")
AICc(m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14) #m10 358
AICc(m3,m6,m7,m9,m12,m13) #m12 371

##----------------------------------------------------------------
#prior predictive simulation, start with neg. binomial part of model

#q prior
alpha_samps<-rnorm(n=500, mean= -0.621, sd= 20)
b1_samps<- rnorm(n=500, mean= 3.771, sd= 20) 

plot(NULL , 
     xlim=range(males$q),
     ylim=range(males$hiphop_offspring,na.rm=TRUE)*2,
     xlab="q", ylab="Prior num offspring",main='alpha ~ N(-1.75,4)  beta ~ N(2.42,6)')
abline( h=min(males$hiphop_offspring,na.rm = TRUE), lty=2,lwd=5 )
abline( h=max(males$hiphop_offspring,na.rm=TRUE), lty=3,lwd=5 ) 
for (i in 1:500){
  curve(alpha_samps[i] + b1_samps[i]*(x),
        from=min(males$q), 
        to=max(males$q),  
        col=alpha("blue",0.2),add=TRUE) }

#length prior
alpha_samps<-rnorm(n=500, mean= -0.621, sd= 20)
b1_samps<- rnorm(n=500, mean= 0.189, sd= 5) 

plot(NULL , 
     xlim=range(scale(males$length)), 
     ylim=range(males$hiphop_offspring,na.rm=TRUE)*2, 
     xlab="length", ylab="Prior num offspring",main='alpha ~ N(-1.75,4)  beta ~ N(0.28,4)')
abline( h=min(males$hiphop_offspring,na.rm = TRUE), lty=2,lwd=5 )
abline( h=max(males$hiphop_offspring,na.rm=TRUE), lty=3,lwd=5 ) 
for (i in 1:500){
  curve(alpha_samps[i] + b1_samps[i]*(x), 
        from=min(scale(males$length)), 
        to=max(scale(males$length)),   #^ extent of x in the code just above
        col=alpha("blue",0.2),add=TRUE) }

#temp prior
alpha_samps<-rnorm(n=500, mean= -0.621, sd= 20)
b1_samps<- rnorm(n=500, mean= 0.188, sd= 5) 

plot(NULL , 
     xlim=range(scale(males$temp)),  
     ylim=range(males$hiphop_offspring,na.rm=TRUE)*2,   
     xlab="temp", ylab="Prior num offspring",main='alpha ~ N(-1.75,4)  beta ~ N(0.28,4)')
abline( h=min(males$hiphop_offspring,na.rm = TRUE), lty=2,lwd=5 ) 
abline( h=max(males$hiphop_offspring,na.rm=TRUE), lty=3,lwd=5 ) 
for (i in 1:500){
  curve(alpha_samps[i] + b1_samps[i]*(x), 
        from=min(scale(males$temp)), 
        to=max(scale(males$temp)),   
        col=alpha("blue",0.2),add=TRUE) }

##----------------------------------------------------------------
#run brms model
library(brms)
func <- bf(hiphop_offspring ~ q + I(q^2) + scale(length) + scale(temp), zi ~ q + I(q^2) + scale(temp) + scale(length))

get_prior(func,males,family=zero_inflated_negbinomial())

priors <- c(set_prior("normal(-0.621,10)",class="Intercept"),
            set_prior("normal(3.771,10)",class="b",coef="q"),
            set_prior("normal(0,5)",class="b",coef="IqE2"),
            set_prior("normal(0.189,5)",class="b",coef="scalelength"),
            set_prior("normal(0.188,5",class="b",coef="scaletemp"),
            set_prior("normal(0,1.5)",class="Intercept",dpar="zi"),
            set_prior("normal(0,1.5)",class="b",coef="q",dpar="zi"),
            set_prior("normal(0,1.5)",class="b",coef="IqE2",dpar="zi"),
            set_prior("normal(0,0.5)",class="b",coef="scalelength",dpar="zi"),
            set_prior("normal(0,1.5)",class="b",coef="scaletemp",dpar="zi"))


m_offspring_model <-brm(formula=func, data=males, family=zero_inflated_negbinomial(), 
                      prior = priors,      
                      control=list(adapt_delta = 0.8, max_treedepth=10), 
                      warmup=1000, iter=3000, thin=1, chains=3, cores=3)   
m_offspring_model
gc()

length(which(posterior_samples(m_offspring_model,pars="b_q")[,1] > 0))/6000 #very close to significant!

func2 <- bf(hiphop_offspring ~ q +scale(length), zi ~ q + scale(length))

priors2 <- c(set_prior("normal(-0.621,10)",class="Intercept"),
            set_prior("normal(3.771,10)",class="b",coef="q"),
            set_prior("normal(0.189,5)",class="b",coef="scalelength"),
            set_prior("normal(0,1.5)",class="Intercept",dpar="zi"),
            set_prior("normal(0,1.5)",class="b",coef="q",dpar="zi"),
            set_prior("normal(0,0.5)",class="b",coef="scalelength",dpar="zi"))


m_offspring_model2 <-brm(formula=func2, data=males, family=zero_inflated_negbinomial(), 
                        prior = priors2,      
                        control=list(adapt_delta = 0.8, max_treedepth=10), 
                        warmup=1000, iter=3000, thin=1, chains=3, cores=3)   
m_offspring_model2

waic(m_offspring_model,m_offspring_model2) #model 2 is better



save(m_offspring_model2,file="./male_offspring_model.Rds")

plot(m_offspring_model2,N=11)
pp_check(m_offspring_model2,nsamples=500,type='bars')


ggplot(males,aes(x=hiphop_offspring)) + geom_bar() + theme_minimal()+ labs(x="Number of sampled offspring",y="Number of males") + guides()


#plot conditional effect of q
#size = 419.8133 mm TL
alpha_samp1 <- sample(posterior_samples(m_offspring_model2,pars="Intercept")[,1],size=6000,replace=TRUE)
q_samp1 <- sample(posterior_samples(m_offspring_model2,pars="q")[,1],size=6000,replace=TRUE)
q_seq <- seq(0,1,by=0.005)


est_upper95 <- c()
est_lower95 <- c()
est_mean <- c()
for (i in 1:length(q_seq)){
  est_dist <- exp(alpha_samp1 + q_samp1*q_seq[i])
  est_upper95[i] <- quantile(est_dist,0.975)
  est_lower95[i] <- quantile(est_dist,0.025)
  est_mean[i] <- mean(est_dist)
}

#ggplot of q conditional effect
cond_q <- data.frame("q"=q_seq,"mean"=est_mean,"upper95"=est_upper95,"lower95"=est_lower95)
p3 <-ggplot(cond_q,aes(x=q,y=upper95)) + geom_line(aes(x=q,y=mean),size=1.3,color="blue") + theme_minimal() +
  labs(x="Male prop. YCT ancestry",y="Number of offspring") +
  geom_ribbon(aes(ymin=lower95,ymax=upper95),alpha=0.4) + scale_fill_manual(values="gray30") +
  theme(panel.grid.minor.x = element_blank(),axis.text=element_text(size=12),axis.title = element_text(size=16)) +
  geom_point(data=males,aes(x=q,y=hiphop_offspring),size=3,fill="black")


#get probability of zero offspring total
alpha_samp1 <- sample(posterior_samples(m_offspring_model2,pars="Intercept")[,1],size=6000,replace=TRUE)
alpha_samp2 <- sample(posterior_samples(m_offspring_model2,pars="Intercept")[,2],size=6000,replace=TRUE)
q_samp1 <- sample(posterior_samples(m_offspring_model2,pars="q")[,1],size=6000,replace=TRUE)
q_samp2 <- sample(posterior_samples(m_offspring_model2,pars="q")[,2],size=6000,replace=TRUE)
shape_samp <- sample(posterior_samples(m_offspring_model2,pars="shape")[,1],size=6000,replace=TRUE)
q_seq <- seq(0,1,by=0.005)

est_upper95 <- c()
est_lower95 <- c()
est_mean <- c()
for (i in 1:length(q_seq)){
  est_dist <- inv_logit(alpha_samp2 + q_samp2*q_seq[i])
  est_dist2 <- exp(alpha_samp1 + q_samp1*q_seq[i])
  prob <- est_dist + (1-est_dist)*dnbinom(0,size=shape_samp,mu=est_dist2)
  est_upper95[i] <- quantile(prob,0.975) 
  est_lower95[i] <- quantile(prob,0.025)
  est_mean[i] <- mean(prob)
}

#ggplot of q conditional effect
cond_q <- data.frame("q"=q_seq,"mean"=est_mean,"upper95"=est_upper95,"lower95"=est_lower95)
p4 <- ggplot(cond_q,aes(x=q,y=upper95)) + geom_line(aes(x=q,y=mean),size=1.3,color="blue") + theme_minimal() +
  labs(x="Male prop. YCT ancestry",y="P(Zero offspring sampled)") + lims(y=c(0,1)) +
  geom_ribbon(aes(ymin=lower95,ymax=upper95),alpha=0.4) + scale_fill_manual(values="gray30") +
  theme(panel.grid.minor.x = element_blank(),axis.text=element_text(size=12),axis.title = element_text(size=16))# +
#geom_point(data=females,aes(x=q,y=offspr_bin),size=3,fill="black")

library(gridExtra)
grid.arrange(p1,p2,p3,p4,ncol=2)



#plot of male & female reprd output by date & q
adult_df <- rbind(males,females)
adult_df$date <- as.Date(adult_df$date,"%m/%d/%Y")
pal <- rev(c("#f16d22","#ffc629","#becb70","#7dcfb6","#00b2ca","#1988af","#1d4e89"))
ggplot(data=adult_df,aes(x=date,y=hiphop_offspring)) + geom_point(shape=21,aes(fill=q, color=q,size=1.5),position = position_jitter(width=0.01,height=0.5)) +
  scale_color_gradientn(name="Prop. YCT\nancestry",colors = pal) + scale_fill_gradientn(colors=alpha(pal,0.7)) +
  labs(x="Entrance date",y="Number of sampled offspring") + guides(size="none",fill="none") +
  geom_vline(xintercept=as.Date("5/27/2019","%m/%d/%Y"),linetype=2)+ theme_minimal() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), axis.title = element_text(size=14))




