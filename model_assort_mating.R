#model assortative mating for NFS trout 2019
#WRosenthal 8/8/2020
#--------------------------------------------------------
#import data

meta <- read.csv("./results_meta_offspring.csv",stringsAsFactors = F)
meta$tributary <- gsub("Trout ","Trout",meta$tributary)
pit_data <- read.csv("./trout.pit.data.2019.csv",row.names = 1,stringsAsFactors = F)
pit_data$gen_ID <- gsub("WR","WCR",pit_data$gen_ID)
matings <- read.csv("./observed_matings_hiphop.csv",stringsAsFactors = F)

#load packages
library(lubridate)
library(pscl)
library(MASS)
library(lme4)
library(msme)
library(pROC)
#--------------------------------------------------------
#look at stay duration as a function of ancestry and in-migration date

##----------------------------------------------------------------
#mixture model to capture 'long-stayers'
mvmt_df$sclday <- scale(mvmt_df$date)
mvmt_df$sclq <- scale(mvmt_df$q)
mvmt_df_long <- mvmt_df[which(mvmt_df$stay_duration >= 100),]
mvmt_df_short <- mvmt_df[which(mvmt_df$stay_duration < 100),]

#frequentist models for short & long stay
library(MASS)
freq_m_short <- glm.nb(stay_duration ~ sclq + sex_bin + sclday, data=mvmt_df_short)
freq_m_long <- glm.nb(stay_duration ~ sclq + sex_bin + sclday, data=mvmt_df_long)

summary(freq_m_short)
summary(freq_m_long)

#prior predictive simulation -- do long & short stay separately

#q prior -- short
alpha_samps<-rnorm(n=500, mean= 3.853, sd= 0.4)
b1_samps<- rnorm(n=500, mean= 0.104, sd= 0.2) 

plot(NULL , 
     xlim=range(mvmt_df_short$sclq),  #<- range of observed age
     ylim=range(mvmt_df_short$stay_duration,na.rm=TRUE)*2,   #<- arbitrarily larger than observed range for mortality
     xlab="q", ylab="Prior stay duration",main='alpha ~ N(3.853,0.4)  beta ~ N(0.104,0.2)')
abline( h=min(mvmt_df_short$stay_duration,na.rm = TRUE), lty=2,lwd=5 ) #mark minimum obs. price
abline( h=max(mvmt_df_short$stay_duration,na.rm=TRUE), lty=3,lwd=5 ) #mark maximum obs. price
for (i in 1:500){
  curve(exp(alpha_samps[i] + b1_samps[i]*(scale(x))), #regression line: NOTE THE Z-score!
        from=min(mvmt_df_short$sclq), 
        to=max(mvmt_df_short$sclq),   #^ extent of x in the code just above
        col=alpha("blue",0.2),add=TRUE) }

#sex_bin prior -- short
alpha_samps<-rnorm(n=500, mean= 3.853, sd= 0.4)
b1_samps<- rnorm(n=500, mean= -0.418, sd= 0.5) 

plot(NULL , 
     xlim=range(as.numeric(mvmt_df_short$sex_bin)),  #<- range of observed age
     ylim=range(mvmt_df_short$stay_duration,na.rm=TRUE)*2,   #<- arbitrarily larger than observed range for mortality
     xlab="sex", ylab="Prior stay duration",main='alpha ~ N(3.853,0.4)  beta ~ N(-0.418,0.5)')
abline( h=min(mvmt_df_short$stay_duration,na.rm = TRUE), lty=2,lwd=5 ) #mark minimum obs. price
abline( h=max(mvmt_df_short$stay_duration,na.rm=TRUE), lty=3,lwd=5 ) #mark maximum obs. price
for (i in 1:500){
  curve(exp(alpha_samps[i] + b1_samps[i]*(scale(x))), #regression line: NOTE THE Z-score!
        from=min(as.numeric(mvmt_df_short$sex_bin)), 
        to=max(as.numeric(mvmt_df_short$sex_bin)),   #^ extent of x in the code just above
        col=alpha("blue",0.2),add=TRUE) }

#day prior -- short
alpha_samps<-rnorm(n=500, mean= 3.853, sd= 0.4)
b1_samps<- rnorm(n=500, mean= -0.079, sd= 0.1) 

plot(NULL , 
     xlim=range(scale(mvmt_df_short$date)),  #<- range of observed age
     ylim=range(mvmt_df_short$stay_duration,na.rm=TRUE)*2,   #<- arbitrarily larger than observed range for mortality
     xlab="date", ylab="Prior stay duration",main='alpha ~ N(4.09,30)  beta ~ N(-0.111,0.5)')
abline( h=min(mvmt_df_short$stay_duration,na.rm = TRUE), lty=2,lwd=5 ) #mark minimum obs. price
abline( h=max(mvmt_df_short$stay_duration,na.rm=TRUE), lty=3,lwd=5 ) #mark maximum obs. price
for (i in 1:500){
  curve(exp(alpha_samps[i] + b1_samps[i]*(scale(x))), #regression line: NOTE THE Z-score!
        from=min(scale(mvmt_df_short$date)), 
        to=max(scale(mvmt_df_short$date)),   #^ extent of x in the code just above
        col=alpha("blue",0.2),add=TRUE) }

#q prior -- long
alpha_samps<-rnorm(n=500, mean= 4.79, sd= 0.15)
b1_samps<- rnorm(n=500, mean= 0.014, sd= 0.08) 

plot(NULL , 
     xlim=range(mvmt_df_long$sclq),  #<- range of observed age
     ylim=c(range(mvmt_df_long$stay_duration,na.rm=TRUE)[1],range(mvmt_df_long$stay_duration,na.rm=TRUE)[2]*2),   #<- arbitrarily larger than observed range for mortality
     xlab="q", ylab="Prior stay duration",main='alpha ~ N(4.095,1)  beta ~ N(0.47,0.2)')
abline( h=min(mvmt_df_long$stay_duration,na.rm = TRUE), lty=2,lwd=5 ) #mark minimum obs. price
abline( h=max(mvmt_df_long$stay_duration,na.rm=TRUE), lty=3,lwd=5 ) #mark maximum obs. price
for (i in 1:500){
  curve(exp(alpha_samps[i] + b1_samps[i]*(scale(x))), #regression line: NOTE THE Z-score!
        from=min(mvmt_df_long$sclq), 
        to=max(mvmt_df_long$sclq),   #^ extent of x in the code just above
        col=alpha("blue",0.2),add=TRUE) }

#sex_bin prior -- long
alpha_samps<-rnorm(n=500, mean= 4.79, sd= 0.15)
b1_samps<- rnorm(n=500, mean= -0.029, sd= 0.02) 

plot(NULL , 
     xlim=range(as.numeric(mvmt_df_long$sex_bin)),  #<- range of observed age
     ylim=c(range(mvmt_df_long$stay_duration,na.rm=TRUE)[1],range(mvmt_df_long$stay_duration,na.rm=TRUE)[2]*2),   #<- arbitrarily larger than observed range for mortality
     xlab="sex", ylab="Prior stay duration",main='alpha ~ N(4.09,30)  beta ~ N(-0.639,50)')
abline( h=min(mvmt_df_long$stay_duration,na.rm = TRUE), lty=2,lwd=5 ) #mark minimum obs. price
abline( h=max(mvmt_df_long$stay_duration,na.rm=TRUE), lty=3,lwd=5 ) #mark maximum obs. price
for (i in 1:500){
  curve(exp(alpha_samps[i] + b1_samps[i]*(scale(x))), #regression line: NOTE THE Z-score!
        from=min(as.numeric(mvmt_df_long$sex_bin)), 
        to=max(as.numeric(mvmt_df_long$sex_bin)),   #^ extent of x in the code just above
        col=alpha("blue",0.2),add=TRUE) }

#day prior -- long
alpha_samps<-rnorm(n=500, mean= 4.79, sd= 0.15)
b1_samps<- rnorm(n=500, mean= -0.067, sd= 0.008) 

plot(NULL , 
     xlim=range(scale(mvmt_df_short$date)),  #<- range of observed age
     ylim=c(range(mvmt_df_long$stay_duration,na.rm=TRUE)[1],range(mvmt_df_long$stay_duration,na.rm=TRUE)[2]*2),   #<- arbitrarily larger than observed range for mortality
     xlab="date", ylab="Prior stay duration",main='alpha ~ N(4.09,30)  beta ~ N(-0.111,0.5)')
abline( h=min(mvmt_df_long$stay_duration,na.rm = TRUE), lty=2,lwd=5 ) #mark minimum obs. price
abline( h=max(mvmt_df_long$stay_duration,na.rm=TRUE), lty=3,lwd=5 ) #mark maximum obs. price
for (i in 1:500){
  curve(exp(alpha_samps[i] + b1_samps[i]*(scale(x))), #regression line: NOTE THE Z-score!
        from=min(scale(mvmt_df_short$date)), 
        to=max(scale(mvmt_df_short$date)),   #^ extent of x in the code just above
        col=alpha("blue",0.2),add=TRUE) }
##----------------------------------------------------------------


#fit mixture model

func <- bf(stay_duration ~ sex_bin + scale(q) + scale(date))
mix_fam <- mixture(negbinomial,poisson)

get_prior(func,mvmt_df,family=mix_fam)

priors <- c(set_prior("normal(3.853,0.5)",class="Intercept",dpar="mu1"),
            set_prior("normal(4.79,0.15)",class="Intercept",dpar="mu2"),
            set_prior("normal(0,1)",class="b",coef="scaleq",dpar="mu1"),
            set_prior("normal(0,1)",class="b",coef="sex_bin1",dpar="mu1"),
            set_prior("normal(0,1)",class="b",coef="scaledate",dpar="mu1"),
            set_prior("normal(0,0.5)",class="b",coef="scaleq",dpar="mu2"),
            set_prior("normal(0,0.5)",class="b",coef="sex_bin1",dpar="mu2"),
            set_prior("normal(0,0.5)",class="b",coef="scaledate",dpar="mu2"),
            set_prior("dirichlet(1)",class="theta"))



mix_model_wide <-brm(formula=func, data=mvmt_df, family=mix_fam, 
                     prior = priors,      
                     control=list(adapt_delta = 0.8, max_treedepth=12), 
                     warmup=2000, iter=6000, thin=1, chains=3, cores=3)   
mix_model_wide
plot(mix_model_wide)
pp_check(mix_model_wide,nsamples=200,type="dens_overlay") + labs(x="stay duration (days)",y="density") +
  theme(axis.title = element_text(size=16))
pp_check(mix_model_wide,nsamples=200,type="violin_grouped",group="sex_bin") + labs(y="stay duration (days)") +
  scale_x_discrete(name="Sex",labels=c("Male","Female")) + theme(axis.title = element_text(size=14),axis.text = element_text(size=12))
pp_check(mix_model_wide,nsamples=200,type="stat_grouped",group="sex_bin",stat="mean")
pp_check(mix_model_wide,nsamples=200,type="stat_grouped",group="sex_bin",stat="median")
pp_check(mix_model_wide,nsamples = 200,type="ecdf_overlay")
IQR <- function(y) { quantile(y, 0.75) - quantile(y, 0.25) } 
pp_check(mix_model_wide, type='stat_grouped',group="sex_bin", stat="IQR",nsamples=200) 

save(mix_model_wide,file="./mix_model.RData")


#-------------------------------------------------------
#begin simulation attempts
#the following section of code was run on UWyo ARCC's Teton HPC

library(brms)
library(lubridate)
#setwd("~/bayes_final/")
meta <- read.csv("/project/ysctrout/wrosenth/NFS_19/assort_mating/results_meta_offspring.csv",stringsAsFactors = F)
meta$tributary <- gsub("Trout ","Trout",meta$tributary)
pit_data <- read.csv("/project/ysctrout/wrosenth/NFS_19/assort_mating/trout.pit.data.2019.csv",row.names = 1,stringsAsFactors = F)
pit_data$gen_ID <- gsub("WR","WCR",pit_data$gen_ID)
matings <- read.csv("/project/ysctrout/wrosenth/NFS_19/assort_mating/observed_matings.csv",stringsAsFactors = F)


mvmt_df <- data.frame(meta$id,meta$tributary, meta$date, meta$sex, meta$length, meta$q,meta$Q,meta$num_offspring,rep(NA,nrow(meta)))
colnames(mvmt_df) <- c("id","trib","date","sex","length","q","Q","num_offspring","stay_duration")
mvmt_df$stay_duration <- sapply(mvmt_df$id,function(x) pit_data$difference[which(pit_data$gen_ID == x)])
mvmt_df$stay_duration <- as.integer(mvmt_df$stay_duration)
mvmt_df <- mvmt_df[which(mvmt_df$trib == "Trout" & is.na(mvmt_df$sex) == FALSE ),]
mvmt_df$date <- as.Date(mvmt_df$date,format = "%m/%d/%Y")
mvmt_df$date <- yday(mvmt_df$date)
mvmt_df$sex_bin <- gsub("Male.*",0,mvmt_df$sex)
mvmt_df$sex_bin <- gsub("Female.*",1,mvmt_df$sex_bin)

load(file="/project/ysctrout/wrosenth/NFS_19/assort_mating/mix_model.RData")

to_replace <- which(is.na(mvmt_df$stay_duration))
predicted_stays <- posterior_predict(mix_model,newdata = mvmt_df[to_replace,],nsamples = 2)


func_big1 <- bf(chosen ~ q + q_dist + scale(length) + (1|female_id))
priors_big1 <- c(set_prior("normal(-4.77,1.5)",class="Intercept"),
                 set_prior("normal(0.878,1.7)",class="b",coef="q"),
                 set_prior("normal(-0.3637,1.3)",class="b",coef="q_dist"),
                 set_prior("normal(0.156,0.8)",class="b",coef="scalelength"),
                 set_prior("normal(0,3)",class="sd",coef="Intercept",group="female_id"))

func_big2 <- bf(chosen ~ q + q_dist + scale(length) + I(q^2) + (1|female_id))
priors_big2 <- c(set_prior("normal(-4.77,1.5)",class="Intercept"),
                 set_prior("normal(0.878,1.7)",class="b",coef="q"),
                 set_prior("normal(-0.3637,1.3)",class="b",coef="q_dist"),
                 set_prior("normal(0.156,0.8)",class="b",coef="scalelength"),
                 set_prior("normal(-0.175,1.5)",class="b",coef="IqE2"),
                 set_prior("normal(0,3)",class="sd",coef="Intercept",group="female_id"))

func_big3 <- bf(chosen ~ q + q_dist + q*q_dist + scale(length) + (1|female_id))
priors_big3 <- c(set_prior("normal(-4.77,1.5)",class="Intercept"),
                 set_prior("normal(0.878,1.7)",class="b",coef="q"),
                 set_prior("normal(-0.3637,1.3)",class="b",coef="q_dist"),
                 set_prior("normal(0.156,0.8)",class="b",coef="scalelength"),
                 set_prior("normal(-0.1992,1.3)",class="b",coef="q:q_dist"),
                 set_prior("normal(0,3)",class="sd",coef="Intercept",group="female_id"))




indv_info <- mvmt_df
indv_info$stay_duration[to_replace] <- predicted_stays[1,]
females <- which(indv_info$sex == "Female")
males <- which(indv_info$sex == "Male")

indv_info$departure <- indv_info$date + indv_info$stay_duration

overlappers <- matrix(NA,nrow=1,ncol=18)
colnames(overlappers) <- c(colnames(indv_info),"chosen","female_q","female_Q","female_length","female_first","female_id","q_dist")
#female first is 1 for when the female arrives before the male, 0 if otherwise
for (f in females){
  if(indv_info$id[f] %in% matings$dam){
    days_in_trib <- seq(indv_info$date[f],indv_info$departure[f])
    male_df <- indv_info[males,]
    
    #get real mates and remove from "not chosen" options
    mates <- matings$sire[which(matings$dam == indv_info$id[f])]
    male_df <- male_df[-which(male_df$id %in% mates),]
    
    #find not chosen mates
    potential_mates <- male_df[which(male_df$date %in% days_in_trib | male_df$departure %in% days_in_trib),]
    potential_mates$chosen <- rep(0,nrow(potential_mates))
    
    #add real mates
    potential_mates <- rbind(potential_mates,data.frame(indv_info[which(indv_info$id %in% mates),],"chosen"=rep(1,length(unique(mates)))))
    potential_mates$female_q <- indv_info$q[f]
    potential_mates$female_Q <- indv_info$Q[f]
    potential_mates$female_length <- indv_info$length[f]
    potential_mates$female_id <- indv_info$id[f]
    potential_mates$female_first <- sapply(potential_mates$date,function(x) ifelse(x > indv_info$date[f],1,0))
    potential_mates$q_dist <- abs(potential_mates$q - potential_mates$female_q) 
    overlappers <- rbind(overlappers,potential_mates)
  }else{
    next()
  }
}
#make dissim and format
overlappers <- overlappers[-1,]
#overlappers$dissim <- abs(overlappers$q - overlappers$female_q)
#q_dat <- rbind(q_dat,overlappers[,c(6,11)])

#model
args <- commandArgs(trailingOnly = TRUE)
print(args)
mate_model <- brm(formula=func_big1,data=overlappers, family=bernoulli(), prior = priors_big1,
                  control=list(adapt_delta = 0.93, max_treedepth=10), 
                  warmup=1000, iter=3000, thin=1, chains=3, cores=3)
saveRDS(rbind(summary(mate_model)$fixed,summary(mate_model)$spec_pars,summary(mate_model)$random$female_id),file=paste("/project/ysctrout/wrosenth/NFS_19/assort_mating/results/model1_",args[1],".rds",sep=""))
m_waic <- waic(mate_model)$estimates
saveRDS(m_waic,file=paste("/project/ysctrout/wrosenth/NFS_19/assort_mating/results/m1_waic_",args[1],".rds",sep=""))
gc()

mate_model2 <- brm(formula=func_big2,data=overlappers, family=bernoulli(), prior = priors_big2,
                   control=list(adapt_delta = 0.93, max_treedepth=10), 
                   warmup=1000, iter=3000, thin=1, chains=3, cores=3)
saveRDS(rbind(summary(mate_model2)$fixed,summary(mate_model2)$spec_pars,summary(mate_model2)$random$female_id),file=paste("/project/ysctrout/wrosenth/NFS_19/assort_mating/results/model2_",args[1],".rds",sep=""))
m_waic2 <- waic(mate_model2)$estimates
saveRDS(m_waic2,file=paste("/project/ysctrout/wrosenth/NFS_19/assort_mating/results/m2_waic_",args[1],".rds",sep=""))
gc()

mate_model3 <- brm(formula=func_big3,data=overlappers, family=bernoulli(), prior = priors_big3,
                   control=list(adapt_delta = 0.93, max_treedepth=10), 
                   warmup=1000, iter=3000, thin=1, chains=3, cores=3)
saveRDS(rbind(summary(mate_model3)$fixed,summary(mate_model3)$spec_pars,summary(mate_model3)$random$female_id),file=paste("/project/ysctrout/wrosenth/NFS_19/assort_mating/results/model3_",args[1],".rds",sep=""))
m_waic3 <- waic(mate_model3)$estimates
saveRDS(m_waic3,file=paste("/project/ysctrout/wrosenth/NFS_19/assort_mating/results/m3_waic_",args[1],".rds",sep=""))

#----------------------------------------------------------
library(ggpubr)
q_dat$chosen <- gsub(0,"not chosen",q_dat$chosen)
q_dat$chosen <- gsub(1,"chosen",q_dat$chosen)
p <- ggdensity(q_dat, x="q",add="mean",color="chosen",fill="chosen",palette=c("#4361EE","#F72585"),xlab="Male prop. YSC Ancestry")
ggpar(p,legend.title = "")

library(ggplot2)
ggplot(mvmt_df,aes(x=stay_duration)) + geom_histogram() + theme_minimal() +
  labs(x="Stay duration (days)") + scale_y_continuous(breaks = seq(0,28,by=4))

#----------------------------------------------------------
load("./coef_data_10ksim.RData")
auc_data <- read.csv("./auc_data_10ksim.csv")

coef_data[[1]]



##----------------------------------------------------------------
#analyze output
library(ggplot2)
library(gridExtra)

models1 <- list()
models2 <- list()
models3 <- list()
m1_waic <- list()
m2_waic <- list()
m3_waic <- list()
for (i in 1:1000){
  m <- readRDS(paste0("./mating_results/",list.files("./mating_results")[grep("model1",list.files("./mating_results"))[i]]))
  models1[[i]] <- m
  m <- readRDS(paste0("./mating_results/",list.files("./mating_results")[grep("model2",list.files("./mating_results"))[i]]))
  models2[[i]] <- m
  m <- readRDS(paste0("./mating_results/",list.files("./mating_results")[grep("model3",list.files("./mating_results"))[i]]))
  models3[[i]] <- m
  
  w <- readRDS(paste0("./mating_results/",list.files("./mating_results")[grep("m1",list.files("./mating_results"))[i]]))
  m1_waic[[i]] <- w
  w <- readRDS(paste0("./mating_results/",list.files("./mating_results")[grep("m2",list.files("./mating_results"))[i]]))
  m2_waic[[i]] <- w
  w <- readRDS(paste0("./mating_results/",list.files("./mating_results")[grep("m3",list.files("./mating_results"))[i]]))
  m3_waic[[i]] <- w
}

waic1_dist <- sapply(m1_waic,function(x) x[3,1])
waic2_dist <- sapply(m2_waic,function(x) x[3,1])
waic3_dist <- sapply(m3_waic,function(x) x[3,1])

waic_df <- rbind(data.frame(value=waic1_dist,model=rep("m1",length(waic1_dist))),
                 data.frame(value=waic2_dist,model=rep("m2",length(waic2_dist))),
                 data.frame(value=waic3_dist,model=rep("m3",length(waic3_dist))))
ggplot() + geom_density(data=waic_df,aes(x=value,group=model,fill=model),alpha=0.65,adjust=2) + theme_minimal() + ylab("Density") + xlab("WAIC") + ggtitle("WAIC, 500 sims")
#they are basically equivalent
#let's just use model 1 since it's easily interpreted

estimates1 <- as.data.frame(t(sapply(models1,function(x) x[,1])))
upper_ci1 <- as.data.frame(t(sapply(models1,function(x) x[,4])))
lower_ci1 <- as.data.frame(t(sapply(models1,function(x) x[,3])))
rhat1 <- as.data.frame(t(sapply(models1,function(x) x[,5])))
bulk_ess1 <- as.data.frame(t(sapply(models1,function(x) x[,6])))

colnames(estimates1) <- rownames(models1[[1]])
colnames(upper_ci1) <- rownames(models1[[1]])
colnames(lower_ci1) <- rownames(models1[[1]])
colnames(rhat1) <- rownames(models1[[1]])
colnames(bulk_ess1) <- rownames(models1[[1]])

icc <- (0.204^2)/((0.204^2)+((3.14159^2)/3)) #0.144 -- pretty low


#convert to long format per variable for plotting
q1 <- rbind(data.frame(value=unlist(estimates1$q),type=rep("mean",nrow(estimates1))),
            data.frame(value=unlist(lower_ci1$q),type=rep("lower_ci",nrow(estimates1))),
            data.frame(value=unlist(upper_ci1$q),type=rep("upper_ci",nrow(estimates1))))

intercept1 <- rbind(data.frame(value=unlist(estimates1$Intercept),type=rep("mean",nrow(estimates1))),
                    data.frame(value=unlist(lower_ci1$Intercept),type=rep("lower_ci",nrow(estimates1))),
                    data.frame(value=unlist(upper_ci1$Intercept),type=rep("upper_ci",nrow(estimates1))))

scalelength1 <- rbind(data.frame(value=unlist(estimates1$scalelength),type=rep("mean",nrow(estimates1))),
                      data.frame(value=unlist(lower_ci1$scalelength),type=rep("lower_ci",nrow(estimates1))),
                      data.frame(value=unlist(upper_ci1$scalelength),type=rep("upper_ci",nrow(estimates1))))

q_dist1 <- rbind(data.frame(value=unlist(estimates1$q_dist),type=rep("mean",nrow(estimates1))),
                 data.frame(value=unlist(lower_ci1$q_dist),type=rep("lower_ci",nrow(estimates1))),
                 data.frame(value=unlist(upper_ci1$q_dist),type=rep("upper_ci",nrow(estimates1))))

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

{p1 <- ggplot() + geom_density(data=q1,aes(x=value,group=type,fill=type),alpha=0.5,adjust=2) + ylab("Density") +
    xlab("Estimate") + ggtitle("Prop. YCT ancestry (q)") + theme(legend.text = element_text(size=12),legend.title = element_text(size=14)) +
    scale_fill_brewer(name="Results from\n1,000 simulations",type="qual",palette = 2,
                      labels=c("Lower credible\ninterval bound","Mean estimate","Upper credible\ninterval bound")) +
    theme_minimal() + theme(plot.title=element_text(size=15)) 
  legend <- get_legend(p1)
  p1 <- p1 + guides(fill="none")
  p2 <- ggplot() + geom_density(data=intercept1,aes(x=value,group=type,fill=type),alpha=0.5,adjust=2) + ylab("Density") + 
    xlab("Estimate") + ggtitle("Intercept") + scale_fill_brewer(type="qual",palette = 2) + theme_minimal() +
    guides(fill="none")+ theme(plot.title=element_text(size=15)) 
  p3 <- ggplot() + geom_density(data=q_dist1,aes(x=value,group=type,fill=type),alpha=0.5,adjust=2) + ylab("Density") + 
    xlab("Estimate") + ggtitle("Difference in prop. YCT ancestry (q diff.)") + scale_fill_brewer(type="qual",palette = 2) + theme_minimal()+
    geom_vline(xintercept = 0,size=1.2) + guides(fill="none")+ theme(plot.title=element_text(size=15)) 
  p4 <- ggplot() + geom_density(data=scalelength1,aes(x=value,group=type,fill=type),alpha=0.5,adjust=2) + ylab("Density") + 
    xlab("Estimate") + ggtitle("Length") + scale_fill_brewer(type="qual",palette = 2) + theme_minimal()+
    geom_vline(xintercept = 0,size=1.2) + guides(fill="none")+ theme(plot.title=element_text(size=15)) 
}
grid.arrange(p2,p1,p3,p4,legend,layout_matrix=cbind(c(1,2),c(3,4),c(5)),widths=c(3.4,3.4,1.4))


length(which(q1$value < 0 & q1$type == "lower_ci")) #doesn't overlap zero!
length(which(q_dist1$value >= 0 & q_dist1$type == "upper_ci")) #only 19/1000 overlaps zero


#-------------------------------------------------------
#make heatmap figure
#use the mean estimate for each parameter, but assume mean length
library(tidyverse)
library(ggplot2)

fq_seq <- seq(0,1,by=0.01)
mq_seq <- seq(0,1,by=0.01)

inv_logit<-function(x){ exp(x)/(1+exp(x)) }
intercept_est <- mean(intercept1$value[which(intercept1$type == "mean")])  #-5.049
#technically don't need this as we'll be using relative prob.
q_est <- mean(q1$value[which(q1$type == "mean")]) #3.264
qdist_est <- mean(q_dist1$value[which(q_dist1$type == "mean")]) #-1.238



prob_df <- as.data.frame(matrix(data=NA,nrow=length(mq_seq),ncol=length(fq_seq))) 
colnames(prob_df) <- fq_seq
prob_df[,ncol(prob_df)+1] <- mq_seq
colnames(prob_df)[ncol(prob_df)] <- "male_q"
for (i in 1:length(fq_seq)){
  qdist <- abs(fq_seq[i] - mq_seq)
  prob_df[,i] <- inv_logit(intercept_est + mq_seq*q_est + qdist*qdist_est)
  #prob_df[,i] <- inv_logit(intercept_est + qdist*qdist_est) #no effect of male q because we know it's biased
}

#make into ggplot-able df
prob_df <- pivot_longer(prob_df,c(1:101),names_to = "female_q", values_to = "prob")
prob_df$rel_prob <- prob_df$prob/min(prob_df$prob)
prob_df$female_q <- as.numeric(prob_df$female_q)

#make plot
#library(wesanderson)
#pal <- wes_palette("Zissou1", 100, type = "continuous")
library(RColorBrewer)
pal <-brewer.pal(n=9,"YlOrRd")

ggplot(prob_df,x=female_q,y=male_q) + geom_raster(aes(x=female_q,y=male_q,fill=rel_prob)) +
  scale_fill_gradientn(colors=pal[2:9],name="Relative prob.\nof mating") + theme_minimal() + scale_x_continuous(breaks=seq(0,1,by=0.25)) +
  theme(panel.grid.minor = element_blank()) + labs(x="Female prop. YCT ancestry",y="Male prop. YCT ancestry")
#this doesn't work well and doesn't need the second dimension.

#plot prob. mating ~ q distance
qdist_seq <- seq(0,1,by=0.01)
prob_df2 <- as.data.frame(matrix(data=NA,nrow=length(qdist_seq),ncol=3))
prob_df2[,1] <- qdist_seq
colnames(prob_df2) <- c("q_dist","prob","rel_prob")
for (i in 1:length(qdist_seq)){
  prob_df2[i,2] <- inv_logit(intercept_est + qdist_seq[i]*qdist_est)
}
prob_df2$rel_prob <- prob_df2$prob/min(prob_df2$prob)

ggplot(prob_df2,mapping=aes(x=q_dist,y=rel_prob)) + geom_line() + theme_minimal() +
  xlab("q distance\n(difference in male and female prop. YCT ancestry)") + ylab("Relative probability of\nhaving a mating detected") +
  theme(axis.text=element_text(size=14),axis.title = element_text(size=16))









