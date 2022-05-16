#simulate sampling of Middle Creek juvs
#based off of Trout Creek juvs
#goal: see if decrease in YCT is due to smaller samples + stochasticity
#6/17/21

library(lubridate)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(grid)
library(ggExtra)
library(RGraphics)
library(evmix)


meta20 <- read.csv("./NFS20_results_meta.csv")

hybrid_status <- seq(1:nrow(meta20))
for (i in 1:length(hybrid_status)){
  if (meta20$q[i] < 0.1){hybrid_status[i] <- "RBT"}
  if (meta20$q[i] > 0.9){hybrid_status[i] <- "YCT"}
  if (meta20$q[i] < 0.6 & meta20$q[i] > 0.4){
    if (meta20$big_Q[i] > 0.8){hybrid_status[i] <- "F1"}
    if (meta20$big_Q[i] > 0.4 & meta20$big_Q[i] < 0.6){hybrid_status[i] <- "F2"}
  }
  if (hybrid_status[i] == i){
    calc <- meta20$big_Q[i] - 2*meta20$q[i]
    if (calc > -0.1 & calc < 0.1){hybrid_status[i] <- "BC RBT"}
    if (hybrid_status[i] == i){
      calc <- meta20$big_Q[i] + 2*meta20$q[i]
      if (calc > 1.9 & calc < 2.1){hybrid_status[i] <- "BC YCT"}
    }
  }
  if (hybrid_status[i] == i){hybrid_status[i] <- "Other"}
}
hybrid_status <- factor(hybrid_status,c("RBT","BC RBT","F1","F2","Other","BC YCT","YCT"))
meta20$hyb_stat <- hybrid_status
adults20 <- meta20[which(!is.na(meta20$Sex) & c(1:nrow(meta20)) %in% grep("WCR2020",meta20$Gen_ID)),]
juv20 <- meta20[which(is.na(meta20$Sex) & c(1:nrow(meta20)) %in% grep("WCR2020", meta20$Gen_ID)),]

meta19 <- read.csv("./../summer_19/parentage/results_meta_offspring.csv")
hybrid_status <- seq(1:nrow(meta19))
for (i in 1:length(hybrid_status)){
  if (meta19$q[i] < 0.1){hybrid_status[i] <- "RBT"}
  if (meta19$q[i] > 0.9){hybrid_status[i] <- "YCT"}
  if (meta19$q[i] < 0.6 & meta19$q[i] > 0.4){
    if (meta19$Q[i] > 0.8){hybrid_status[i] <- "F1"}
    if (meta19$Q[i] > 0.4 & meta19$Q[i] < 0.6){hybrid_status[i] <- "F2"}
  }
  if (hybrid_status[i] == i){
    calc <- meta19$Q[i] - 2*meta19$q[i]
    if (calc > -0.1 & calc < 0.1){hybrid_status[i] <- "BC RBT"}
    if (hybrid_status[i] == i){
      calc <- meta19$Q[i] + 2*meta19$q[i]
      if (calc > 1.9 & calc < 2.1){hybrid_status[i] <- "BC YCT"}
    }
  }
  if (hybrid_status[i] == i){hybrid_status[i] <- "Other"}
}
hybrid_status <- factor(hybrid_status,c("RBT","BC RBT","F1","F2","Other","BC YCT","YCT"))
meta19$hyb_stat <- hybrid_status
adults19 <- meta19[which(!is.na(meta19$sex)),]
juv19 <- meta19[which(is.na(meta19$sex)),]


###start
kern <- "b"
bw_adj_t <- 0.5
bandwidth <- 0.025
n_spots <- 100

#reflect data first to account for boundaries
adult19_r <- c(-adults19$q, adults19$q, 2-adults19$q)
adult20_r <- c(-adults20$q, adults20$q, 2-adults20$q)
juv19_r <- c(-juv19$q, juv19$q, 2-juv19$q)
juv20_r <- c(-juv20$q, juv20$q, 2-juv20$q)

#now estimate densities
adult19_d <- density(adult19_r,kernel=kern,adj=bw_adj_t,bw=bandwidth,from=0,to=1,n=n_spots)
adult20_d <- density(adult20_r,kernel=kern,bw=bandwidth,from=0,to=1,n=n_spots)
juv19_d <- density(juv19_r,kernel=kern,adj=bw_adj_t,bw=bandwidth,from=0,to=1,n=n_spots)
juv20_d <- density(juv20_r,kernel=kern,bw=bandwidth,from=0,to=1,n=n_spots)

#use histogram for Trout Creek to avoid boundary effects
adult19_h <- hist(adults19$q,breaks=seq(0,1,by=0.01),plot=F)$density
juv19_h <- hist(juv19$q,breaks=seq(0,1,by=0.01),plot=F)$density

#multiply by 3 because we have 3x the data points
adult19_d$y <- adult19_d$y*3
adult20_d$y <- adult20_d$y*3
juv19_d$y <- juv19_d$y*3
juv20_d$y <- juv20_d$y*3

#get differences
#diff_d19 <- data.frame(x=adult19_d$x, y = juv19_d$y - adult19_d$y)
diff_d20 <- data.frame(x=adult20_d$x, y = juv20_d$y - adult20_d$y)

#plot differences in adult and juvenile densities
#p1 <- ggplot(diff_d19,aes(x=x,y=y)) + geom_line(color="blue") + theme_minimal() +
#  geom_segment(aes(x=0,xend=1,y=0,yend=0),lwd=2) + 
#  xlab(element_blank()) + ylab(element_blank()) + ggtitle("Trout Creek 2019")
#p2 <- ggplot(diff_d20,aes(x=x,y=y)) + geom_line(color="blue") + theme_minimal() +
#  geom_segment(aes(x=0,xend=1,y=0,yend=0),lwd=2) +
#  xlab(element_blank()) + ylab(element_blank()) + ggtitle("Middle Creek 2020")
#tx <- "Proportion YCT ancestry (q)"
#ty <- "Difference between\njuvenile and adult\npoint densities"

#{grid.arrange(p1,p2,textGrob(""),ncol=2,layout_matrix=cbind(c(3,3,3),c(1,2,3)),widths=c(0.6,3.5),heights=c(2.4,2.4,0.3))
#  grid.text(tx,x=unit(0.6,"npc"),y=unit(0.05,"npc"),gp=gpar(fontsize=14))
#  grid.text(ty,x=unit(0.08,"npc"),y=unit(0.55,"npc"),gp=gpar(fontsize=14))
#}
#decrease in YCT from Middle Creek adult -> juv generations is weird
#is it just due to chance?


#resample from Middle & Trout Creek adults
#see how that matches up with our observations

hyp_d <- data.frame(x=adult20_d$x, y=adult20_d$y)
hyp_d_t <- data.frame(x=adult19_d$x, y=adult19_d$y)

#sample from density "distribution" in the same way we sampled irl
nsim <- 10000
hyp_samp <- matrix(data=NA,ncol=nrow(juv20),nrow=nsim)
hyp_samp_t <- matrix(data=NA,ncol=nrow(juv19),nrow=nsim)

for (s in 1:nsim){
  hyp_samp[s,] <- sample(hyp_d$x, size=ncol(hyp_samp), replace = T, prob = hyp_d$y)
  #hyp_samp_t[s,] <- sample(adults19$q, size=ncol(hyp_samp_t), replace = T)
  hyp_samp_t[s,] <- sample(hyp_d_t$x, size=ncol(hyp_samp_t), replace = T, prob=hyp_d_t$y)
}

#now get densities for each sample
hyp_d_samp <- matrix(data=NA,nrow=length(hyp_d$x), ncol=nsim+1)
hyp_d_samp_t <- matrix(data=NA,nrow=length(hyp_d_t$x), ncol=nsim+1)
hyp_d_samp[,1] <- seq(0.01,1,by=0.01)
hyp_d_samp_t[,1] <- seq(0.01,1,by=0.01)
for (s in 1:nsim){
  #ref <- c(-hyp_samp[s,], hyp_samp[s,], 2-hyp_samp[s,])
  #ref_t <- c(-hyp_samp_t[s,], hyp_samp_t[s,], 2-hyp_samp_t[s,])
  
  #d <- (density(ref,kernel=kern,adj=bw_adj,bw=bandwidth,from=0,to=1,n=n_spots)$y)*3
  #hyp_d_samp[,s+1] <- d
  
  #d_t <- (density(ref_t,kernel=kern,adj=bw_adj,bw=bandwidth,from=0,to=1,n=n_spots)$y)*3
  #hyp_d_samp_t[,s+1] <- d_t
  
  #try using histogram bins instead
  d <- hist(hyp_samp[s,],breaks=seq(0,1,by=0.01),plot=F)
  d <- d$density
  hyp_d_samp[,s+1] <- d
  
  d_t <- hist(hyp_samp_t[s,],breaks=seq(0,1,by=0.01),plot=F)
  d_t <- d_t$density
  hyp_d_samp_t[,s+1] <- d_t
  
}
rownames(hyp_d_samp) <- hyp_d_samp[,1]
rownames(hyp_d_samp_t) <- hyp_d_samp_t[,1]
hyp_d_samp <- hyp_d_samp[,-1] #trim first column
hyp_d_samp_t <- hyp_d_samp[,-1] #trim first column

#find 0.025 and 0.975 quantiles for each x point
hyp_bounds <- data.frame(x = c(hyp_d$x,hyp_d$x), upper=NA, lower=NA)
hyp_bounds$upper <- apply(hyp_d_samp,1, function(x) quantile(x,probs=0.975))
hyp_bounds$lower <- apply(hyp_d_samp,1, function(x) quantile(x,probs=0.025))
hyp_bounds$median <- apply(hyp_d_samp,1, function(x) quantile(x,probs=0.5))

hyp_bounds_t <- data.frame(x = c(seq(0.01,1,by=0.01),seq(0.01,1,by=0.01)), upper=NA, lower=NA)
hyp_bounds_t$upper <- apply(hyp_d_samp_t,1, function(x) quantile(x,probs=0.975))
hyp_bounds_t$lower <- apply(hyp_d_samp_t,1, function(x) quantile(x,probs=0.025))
hyp_bounds_t$median <- apply(hyp_d_samp_t,1, function(x) quantile(x,probs=0.5))

#add observed
obs_dens <- data.frame(x=c(juv20_d$x,juv20_d$x))
obs_dens$obs <- c(juv20_d$y,adult20_d$y)
obs_dens$obs_type <- c(rep("juv",length(juv20_d$y)),rep("adult",length(adult20_d$y)))

obs_dens_t <- data.frame(x=c(seq(0.01,1,by=0.01),seq(0.01,1,by=0.01)))
obs_dens_t$obs <- c(juv19_d$y,adult19_d$y)
obs_dens_t$obs_type <- c(rep("juv",length(juv19_h)),rep("adult",length(adult19_h)))

#plot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


{p1 <- ggplot() + geom_line(data=obs_dens_t,aes(x=x,y=obs,color=obs_type),lwd=1.2) +
  xlab("") + ylab("") + theme_minimal() + 
  geom_ribbon(data=hyp_bounds_t,aes(x=x,ymin=lower,ymax=upper,alpha=70,fill="gray60")) +
  guides(alpha=FALSE,color=guide_legend(order=1)) + ggtitle("Trout Creek 2019") +
  scale_color_manual(name="",values=c("blue","red"),labels=c("Observed adults","Observed juveniles"))+
  scale_fill_manual(name="",values="gray60",labels="Resampled adults") +
  theme(legend.title = element_blank(),legend.spacing.y = unit(0.000001,"npc"))
legend <- get_legend(p1)
p1 <- p1 + guides(fill=FALSE,color=FALSE)


p2 <- ggplot() + geom_line(data=obs_dens,aes(x=x,y=obs,color=obs_type),lwd=1.2) +
  xlab("") + ylab("") + theme_minimal() +
  geom_ribbon(data = hyp_bounds,aes(x=x,ymin=lower,ymax=upper,alpha=65,fill="gray")) +
  guides(alpha=FALSE,color=guide_legend(order=1)) + ggtitle("Middle Creek 2020") +
  scale_color_manual(name="",values=c("blue","red"),labels=c("Observed adults","Observed juveniles"))+
  scale_fill_manual(name="",values="gray",labels="Resampled adults") + guides(fill=FALSE,color=FALSE)
tx <- "Proportion YCT ancestry (q)"
ty <- "Density"  
}

{grid.arrange(p1,p2,legend,textGrob(""),ncol=3,layout_matrix=cbind(c(4,4,4),c(1,2,4),c(3,3,3)),widths=c(0.1,3.5,0.9),heights=c(2.4,2.4,0.3))
  grid.text(tx,x=unit(0.45,"npc"),y=unit(0.05,"npc"),gp=gpar(fontsize=14))
  grid.text(ty,x=unit(0.025,"npc"),y=unit(0.55,"npc"),gp=gpar(fontsize=14),rot=90)
}

p3 <- ggplot() + geom_line(data=obs_dens,aes(x=x,y=obs,color=obs_type),lwd=1.2) +
  geom_ribbon(data = hyp_bounds,aes(x=x,ymin=lower,ymax=upper,alpha=65,fill="gray")) +
  scale_color_manual(name="",values=c("blue","red"),labels=c("Observed adults","Observed juveniles")) +
  scale_fill_manual(name="",values="gray",labels="Resampled adults") +
  guides(alpha=FALSE,color=guide_legend(order=1)) +
  xlab("Proportion YCT ancestry (q)") + ylab("Density") + theme_minimal() +
  theme(legend.title = element_blank(),legend.spacing.y = unit(0.000001,"npc"))


###stop

#calculate difference between hypothetical samples and adult generation
#hyp_d_samp_diff <- hyp_d_samp
#hyp_d_samp_diff <- apply(hyp_d_samp_diff, 2, function(x) (x - adult20_d$y))


#find 0.025 and 0.975 quantiles for each x point
#hyp_diff_bounds <- data.frame(x = hyp_d$x, upper=NA, lower=NA)
#hyp_diff_bounds$upper <- apply(hyp_d_samp_diff,1, function(x) quantile(x,probs=0.975))
#hyp_diff_bounds$lower <- apply(hyp_d_samp_diff,1, function(x) quantile(x,probs=0.025))
#hyp_diff_bounds$median <- apply(hyp_d_samp_diff,1, function(x) quantile(x,probs=0.5))

#add observed
#hyp_diff_bounds$obs <- diff_d20$y

#plot!
#ggplot(hyp_diff_bounds) + geom_line(aes(x=x,y=obs, group=1,color="blue")) +
#  geom_ribbon(aes(x=x,ymin=lower,ymax=upper, group=1, fill="gray80",alpha=65)) +
#  theme_minimal() + guides(alpha=FALSE) + xlab("Proportion YCT ancestry (q)") +
#  ylab("Density difference\nfrom observed\nadult generation") + geom_segment(aes(x=0,xend=1,y=0,yend=0),lwd=2) +
#  scale_fill_identity(name="",guide="legend",labels= c("Resampled adult generation")) +
#  scale_color_identity(name="",guide="legend",labels=c("Observed juvenile generation")) +
#  theme(axis.title.y = element_text(angle=0,vjust=0.5))

plot(adult20_d,lwd=2)
lines(juv20_d$x,juv20_d$y,col="red",lwd=2)
for (i in 1:1000){
  lines(juv20_d$x,hyp_d_samp[,runif(1,1,10000)],col=adjustcolor("green",alpha.f = 0.05))
}




#test different bw for Trout Creek
kern <- "b"
test <- density(adult19_r,bw=0.02,kernel=kern,from=0,to=1,n=n_spots)
test2 <- density(adult19_r,bw=0.025,kernel=kern,from=0,to=1,n=n_spots)
test3 <- density(adult19_r,bw=0.015,kernel=kern,from=0,to=1,n=n_spots)
test4 <- density(adult19_r,bw=0.01,kernel=kern,from=0,to=1,n=n_spots)
test5 <- density(adult19_r,bw=0.008,kernel=kern,from=0,to=1,n=n_spots)

hist(adults19$q,freq=F,breaks=seq(0,1,by=0.01))
lines(test$x,test$y*3,col="red",lwd=2)
lines(test$x,test2$y*3,col="blue",lwd=2)
lines(test$x,test3$y*3,col="green",lwd=2)
lines(test$x,test4$y*3,col="orange",lwd=2)
lines(test$x,test4$y*3,col="purple",lwd=2)

#P-splines smoothing
his = hist(adults20$q, breaks = 200, plot = F)
X = his$counts
u = his$mids

# Prepare basis (I-mat) and penalty (1st difference)
B = diag(length(X))
D1 = diff(B, diff = 1)
lambda = 20 # fixed but can be selected (e.g. AIC)
P = lambda * t(D1) %*% D1

# Smooth
tol = 1e-8
eta = log(X + 1)
for (it in 1:20) 
{
  mu = exp(eta)
  z = X - mu + mu * eta
  a = solve(t(B) %*% (c(mu) * B) + P, t(B) %*% z)
  etnew = B %*% a
  de = max(abs(etnew - eta))
  cat('Crit', it, de, '\n')
  if(de < tol) break
  eta = etnew
}

# Plot
plot(u, exp(eta), ylim = c(0, max(X)), type = 'l', col = 2)
lines(u, X, type = 'h')


test <- fpsden(adults20$q,lambdaseq = 10^seq(-5,-1),breaks=seq(0,1,length.out=100),ord = )
dens <- exp(test$bsplines %*% test$mle)
plot(test$mids,dens/test$nbinwidth,type="l")




#difference between lower bound and obs for high q
test_bound_diff <- hyp_diff_bounds[which(hyp_diff_bounds$x >=0.9),]
test_bound_diff$diff <- test_bound_diff$lower - test_bound_diff$obs
plot(test_bound_diff$x, test_bound_diff$diff,type="l")

outside_bounds <- rep(NA,nrow(hyp_diff_bounds))
for (i in 1:nrow(hyp_diff_bounds)){
  if (hyp_diff_bounds$obs[i] > hyp_diff_bounds$upper[i] | hyp_diff_bounds$obs[i] < hyp_diff_bounds$lower[i]){
    outside_bounds[i] <- 1
  } else {
    outside_bounds[i] <- 0
  }
}
















