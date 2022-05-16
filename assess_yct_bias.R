#calculate the probability of preference for YCT males being a product of increased YCT reproductive output
library(brms)
library(ggplot2)
library(tidyr)
library(gridExtra)
inv_logit<-function(x){ exp(x)/(1+exp(x)) }


#load in male reproductive output results
load("./../parentage/male_offspring_model.Rds")

#load in assortative mating model results
models1 <- list()

for (i in 1:1000){
  m <- readRDS(paste0("./mating_results/",list.files("./mating_results")[grep("model1_fullppdist",list.files("./mating_results"))[i]]))
  models1[[i]] <- m
}

q_ppdist <- c()
for (i in 1:1000){
  q_ppdist <- c(q_ppdist,models1[[i]]$b_q)
}
rm(models1)


#calculate probability of >0 offspring for pure YCT

alpha_samp1 <- sample(posterior_samples(m_offspring_model,pars="Intercept")[,1],size=6000,replace=TRUE)
alpha_samp2 <- sample(posterior_samples(m_offspring_model,pars="Intercept")[,2],size=6000,replace=TRUE)
q_samp1 <- sample(posterior_samples(m_offspring_model,pars="q")[,1],size=6000,replace=TRUE)
zi_q_samp1 <- sample(posterior_samples(m_offspring_model,pars="q")[,3],size=6000,replace=TRUE)
shape_samp <- sample(posterior_samples(m_offspring_model,pars="shape")[,1],size=6000,replace=TRUE)
q_seq <- c(0,1)

est_upper95 <- c()
est_lower95 <- c()
est_mean <- c()
prob_df <- as.data.frame(matrix(nrow=6000,ncol=2))
colnames(prob_df) <- c("RBT","YCT")
for (i in 1:length(q_seq)){
  est_dist <- inv_logit(alpha_samp2 + zi_q_samp1*q_seq[i])
  est_dist2 <- exp(alpha_samp1 + q_samp1*q_seq[i])
  prob <- (1-est_dist)*(1-dnbinom(0,size=shape_samp,mu=est_dist2))
  est_upper95[i] <- quantile(prob,0.975) 
  est_lower95[i] <- quantile(prob,0.025)
  est_mean[i] <- mean(prob)
  prob_df[,i] <- prob
}

prob_df <- pivot_longer(prob_df,c(1,2),names_to="spp",values_to = "Probability")
ggplot(prob_df,mapping=aes(x=Probability,color=spp)) + geom_density() + theme_minimal()

yct_prob <- prob_df$Probability[which(prob_df$spp=="YCT")]
rbt_prob <- prob_df$Probability[which(prob_df$spp=="RBT")]

#convert probability difference between RBT and YCT to a coefficient distribution
diff_logit <- logit_scaled(yct_prob) - logit_scaled(rbt_prob)

#compare this new coefficient distribution with assort mating coefficient distribution

#make them the same size first
assort_mating_coef <- sample(q_ppdist,6000,replace = T)

#amount of overlap
length(which(diff_logit >= min(assort_mating_coef)))/6000 #60%
overlap_amt <- round((length(which(diff_logit >= min(assort_mating_coef)))/6000)*100,digits = 2)
#overlap for >0.025 quantile
length(which(diff_logit >= quantile(assort_mating_coef,0.025)))/6000 #30%


#plot
coef_plot_df <- data.frame("assort_mating"=assort_mating_coef,"reprd_output"=diff_logit)
coef_plot_df <- pivot_longer(coef_plot_df,c(1,2),names_to = "Model",values_to = "coef_estimate")



ggplot(coef_plot_df,mapping=aes(x=coef_estimate,fill=Model)) + geom_density(alpha=0.7) + theme_minimal() +
  scale_color_manual("gray50") +  xlab("Coefficient estimate") + ylab("Density") + 
  scale_fill_discrete(name="Model",labels=c("Mate Choice","Male Reproductive\nOutput")) +
  annotate("text",x=-5,y=0.4,label=paste0(overlap_amt,"% of Male Reproductive\nOutput values greater than\nminimum Mate Choice value"),
            size=3.5)



