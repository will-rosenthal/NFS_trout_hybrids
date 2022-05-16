#simulate offspring generations
#3/12/21
#update 3/30/21

library(lubridate)
library(tidyverse)
matings <- read.csv("./observed_matings.csv")
nfs19 <- read.csv("./NFS19_all_0.5_0.03_results_meta.csv")
nfs20 <- read.csv("./NFS20_results_meta.csv")

#what is the range of entry dates for adults that mated with each other?
matings <- unique(matings[,1:2])
matings$male_entry <- nfs19$date[match(matings$sire,nfs19$id)]
matings$female_entry <- nfs19$date[match(matings$dam,nfs19$id)]
matings$male_entry <- as.Date(matings$male_entry,format="%m/%d/%Y")
matings$female_entry <- as.Date(matings$female_entry,format="%m/%d/%Y")
matings$diff <- as.numeric(abs(matings$male_entry - matings$female_entry))

hist(matings$diff,labels=T)

#make a model to fit a neg.bin. dist. to these data
library(lme4)
library(MASS)
m1 <- glm.nb(diff ~ 1, data=matings)
#plot dist
x <- seq(0,36,by=1)
y <- dnbinom(x,mu=unname(coef(m1)),size=m1$theta)
plot(y~x,type="l")
#looks good?


#start simulation for "informative null"
nsim <- 5000
middle2020_sim_inform <- as.data.frame(matrix(ncol=nsim,nrow=length(which(nfs20$Location == "Middle" & is.na(nfs20$Sex)))))
middle2020_sim_assort <- as.data.frame(matrix(ncol=nsim,nrow=length(which(nfs20$Location == "Middle" & is.na(nfs20$Sex)))))
middle2020_sim_null <- as.data.frame(matrix(ncol=nsim,nrow=length(which(nfs20$Location == "Middle" & is.na(nfs20$Sex)))))
middle2020_sim_yct <- as.data.frame(matrix(ncol=nsim,nrow=length(which(nfs20$Location == "Middle" & is.na(nfs20$Sex)))))

middle2020_adults <- nfs20[which(!is.na(nfs20$Sex) & nfs20$Location == "Middle"),]
middle2020_fem <- middle2020_adults[which(middle2020_adults$Sex == "Female"),]
middle2020_male <- middle2020_adults[which(middle2020_adults$Sex == "Male"),]
middle2020_fem$Date <- as.Date(middle2020_fem$Date,format = "%m/%d/%y")
middle2020_male$Date <- as.Date(middle2020_male$Date,format = "%m/%d/%y")

num_offspr <- 30
for (s in 1:nsim){
  if (s %% 100 == 0){
    cat("Sim #",s,fill=T)
  }
  #first step: decide random(ish) matings
  #set the number of days between the female entry and male entry
  middle2020_fem$diff <- rnbinom(nrow(middle2020_fem),mu=unname(coef(m1)),size=m1$theta)  
  
  #females for informative null
  inform_offspring <- as.data.frame(matrix(ncol=1,nrow=nrow(middle2020_fem)*num_offspr))
  assort_offspring <- as.data.frame(matrix(ncol=1,nrow=nrow(middle2020_fem)*num_offspr))
  yct_offspring <- as.data.frame(matrix(ncol=1,nrow=nrow(middle2020_fem)*num_offspr))
  for (f in 1:nrow(middle2020_fem)){
    #cat("Female #",f,fill=T)
    possible_mates <- as.data.frame(matrix(nrow=0,ncol=10,data=NA))
    while (nrow(possible_mates) == 0){
      possible_mates <- middle2020_male[which(abs(middle2020_male$Date - middle2020_fem$Date[f]) == middle2020_fem$diff[f]),]
      if (nrow(possible_mates) == 0){
        #print("Resampling diff")
        middle2020_fem$diff[f] <- rnbinom(1,mu=unname(coef(m1)),size=m1$theta)
      }
    }
    
    actual_mate <- possible_mates[sample(1:nrow(possible_mates),1),]
    #once we have the mate, we need to create some offspring
    #use binomial distribution to get the number of YCT chromosomes from each parent -- update 3/30/21
    fem_q <- middle2020_fem$q[f]
    male_q <- actual_mate$q
    num_chrom <- 32
    fem_chrom <- rbinom(num_offspr,num_chrom,fem_q)
    male_chrom <- rbinom(num_offspr,num_chrom,male_q)
    juv_q <- (fem_chrom + male_chrom)/(num_chrom * 2)

    inform_offspring[c((((f - 1)*num_offspr)+1):(((f - 1)*num_offspr)+num_offspr)),1] <- juv_q
    
    possible_mates <- as.data.frame(matrix(nrow=0,ncol=10,data=NA))
    while (nrow(possible_mates) == 0){
      possible_mates <- middle2020_male[which(abs(middle2020_male$Date - middle2020_fem$Date[f]) == middle2020_fem$diff[f]),]
      if (nrow(possible_mates) == 0){
        #print("Resampling diff")
        middle2020_fem$diff[f] <- rnbinom(1,mu=unname(coef(m1)),size=m1$theta)
      }
    }
    
    #do assortative mating
    assort_mate <- possible_mates[which.min(abs(possible_mates$q - fem_q)),]
    #once we have the mate, we need to create some offspring
    #use binomial distribution to get the number of YCT chromosomes from each parent -- update 3/30/21
    male_q <- assort_mate$q
    num_chrom <- 32
    fem_chrom <- rbinom(num_offspr,num_chrom,fem_q)
    male_chrom <- rbinom(num_offspr,num_chrom,male_q)
    juv_q <- (fem_chrom + male_chrom)/(num_chrom * 2)
    
    assort_offspring[c((((f - 1)*num_offspr)+1):(((f - 1)*num_offspr)+num_offspr)),1] <- juv_q
    
    possible_mates <- as.data.frame(matrix(nrow=0,ncol=10,data=NA))
    while (nrow(possible_mates) == 0){
      possible_mates <- middle2020_male[which(abs(middle2020_male$Date - middle2020_fem$Date[f]) == middle2020_fem$diff[f]),]
      if (nrow(possible_mates) == 0){
        #print("Resampling diff")
        middle2020_fem$diff[f] <- rnbinom(1,mu=unname(coef(m1)),size=m1$theta)
      }
    }
    
    #do preference for yct 
    yct_mate <- possible_mates[which.max(possible_mates$q),]
    #once we have the mate, we need to create some offspring
    #use binomial distribution to get the number of YCT chromosomes from each parent -- update 3/30/21
    male_q <- yct_mate$q
    num_chrom <- 32
    fem_chrom <- rbinom(num_offspr,num_chrom,fem_q)
    male_chrom <- rbinom(num_offspr,num_chrom,male_q)
    juv_q <- (fem_chrom + male_chrom)/(num_chrom * 2)
    
    yct_offspring[c((((f - 1)*num_offspr)+1):(((f - 1)*num_offspr)+num_offspr)),1] <- juv_q
    
  }
  
  #now sample from big pool of juveniles to mimic real sampling
  middle2020_sim_inform[,s] <- inform_offspring[sample(c(1:nrow(inform_offspring)),length(which(nfs20$Location == "Middle" & is.na(nfs20$Sex)))),]
  middle2020_sim_assort[,s] <- assort_offspring[sample(c(1:nrow(assort_offspring)),length(which(nfs20$Location == "Middle" & is.na(nfs20$Sex)))),]
  middle2020_sim_yct[,s] <- yct_offspring[sample(c(1:nrow(yct_offspring)),length(which(nfs20$Location == "Middle" & is.na(nfs20$Sex)))),]
  #--------------------------------------------------------------------------------------------------
  #do real null model -- totally random matings
  
  null_offspring <- as.data.frame(matrix(ncol=1,nrow=nrow(middle2020_fem)*num_offspr))
  for (f in 1:nrow(middle2020_fem)){
    #cat("Female #",f,fill=T)
    actual_mate <- middle2020_male[sample(1:nrow(middle2020_male),1),]
    #once we have the mate, we need to create some offspring
    #use binomial distribution to get the number of YCT chromosomes from each parent -- update 3/30/21
    fem_q <- middle2020_fem$q[f]
    male_q <- actual_mate$q
    num_chrom <- 32
    fem_chrom <- rbinom(num_offspr,num_chrom,fem_q)
    male_chrom <- rbinom(num_offspr,num_chrom,male_q)
    juv_q <- (fem_chrom + male_chrom)/(num_chrom * 2)
    
    null_offspring[c((((f - 1)*num_offspr)+1):(((f - 1)*num_offspr)+num_offspr)),1] <- juv_q
  }
  
  #now sample from big pool of juveniles to mimic real sampling
  middle2020_sim_null[,s] <- null_offspring[sample(c(1:nrow(inform_offspring)),length(which(nfs20$Location == "Middle" & is.na(nfs20$Sex)))),]
  
}

middle2020_juv <- nfs20[which(is.na(nfs20$Sex) & nfs20$Location == "Middle"),]
#library(ggplot2)
#middle2020_sim_inform_long <- gather(middle2020_sim_inform,key="sim",value="q",c(1:ncol(middle2020_sim_inform)))
#p <- ggplot() + geom_density(data=middle2020_sim_inform,aes(x=V5), col="black",fill="seagreen", alpha=0.7)
#p + geom_density(data=middle2020_juv,aes(x=q),fill="firebrick2",col="black",alpha=0.7)

#p2 <- ggplot() + geom_density(data=middle2020_sim_inform_long,aes(x=q), col="black",fill="seagreen", alpha=0.7)
#p2 + geom_density(data=middle2020_juv,aes(x=q),fill="firebrick2",col="black",alpha=0.7)

prop_yct <- as.data.frame(matrix(ncol=4,nrow=nsim))
colnames(prop_yct) <- c("informed_random","random","informed_assort","informed_yct")
for (i in 1:ncol(middle2020_sim_inform)){
  prop_yct[i,1] <- length(which(middle2020_sim_inform[,i] >=0.9))/nrow(middle2020_sim_inform)
  prop_yct[i,2] <- length(which(middle2020_sim_null[,i] >=0.9))/nrow(middle2020_sim_null)
  prop_yct[i,3] <- length(which(middle2020_sim_assort[,i] >=0.9))/nrow(middle2020_sim_assort)
  prop_yct[i,4] <- length(which(middle2020_sim_yct[,i] >=0.9))/nrow(middle2020_sim_yct)
}
obs_prop_yct <- length(which(middle2020_juv$q >=0.9))/nrow(middle2020_juv)

prop_yct <- readRDS("./sim_prop_yct.rds")

#get quantiles
inf_rand <- ecdf(prop_yct[,1])(obs_prop_yct) #informed random 0.9048
rand <- ecdf(prop_yct[,2])(obs_prop_yct) #random 0.9998
inf_asst <- ecdf(prop_yct[,3])(obs_prop_yct) #informed assortative 0.1374
inf_yct <- ecdf(prop_yct[,4])(obs_prop_yct) #informed preference for yct 0.1358

#ggplot
library(ggplot2)
library(gridExtra)

p1 <- ggplot(data=prop_yct,aes(x=random)) + geom_histogram(binwidth = 0.01) +
  labs(x="Proportion YCT juveniles",y="Count", title = "Random mating") + theme_minimal() +
  geom_vline(xintercept = obs_prop_yct,color="firebrick2",size=1.5) + theme(panel.grid.minor = element_blank()) +
  annotate("text",x=0.13,y=900,label=paste0(sprintf("%.4f",rand),"\n quantile"))
p2 <- ggplot(data=prop_yct,aes(x=informed_random)) + geom_histogram(binwidth = 0.01) +
  labs(x="Proportion YCT juveniles",y="Count",title = "Informed random mating") + theme_minimal() +
  geom_vline(xintercept = obs_prop_yct,color="firebrick2",size=1.5) + theme(panel.grid.minor = element_blank()) +
  annotate("text",x=0.165,y=650,label=paste0(inf_rand,"\n quantile"))
p3 <- ggplot(data=prop_yct,aes(x=informed_assort)) + geom_histogram(binwidth = 0.01) +
  labs(x="Proportion YCT juveniles",y="Count", title="Informed assortative mating") + theme_minimal() +
  geom_vline(xintercept = obs_prop_yct,color="firebrick2",size=1.5) + theme(panel.grid.minor = element_blank()) +
  annotate("text",x=0.285,y=550,label=paste0(inf_asst,"\n quantile"))
p4 <- ggplot(data=prop_yct,aes(x=informed_yct)) + geom_histogram(binwidth = 0.01) +
  labs(x="Proportion YCT juveniles",y="Count",title="Informed YCT preference mating") + theme_minimal() +
  geom_vline(xintercept = obs_prop_yct,color="firebrick2",size=1.5)+ theme(panel.grid.minor = element_blank()) +
  annotate("text",x=0.25,y=525,label=paste0(inf_yct,"\n quantile"))

grid.arrange(p1,p2,p3,p4)









