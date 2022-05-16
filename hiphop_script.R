#run hiphop
library(hiphop)

#try new vcfs
library(vcfR)
test <- read.vcfR("./NFS19_parentage_dp5_miss0.90_maf0.01_thin10k.recode.vcf")
colnames(test@gt) <- sapply(colnames(test@gt),function(x) gsub(".*aln_(.*).sorted.bam","\\1",x)) #change names
metadata <- read.csv("./../entropy/all_samp/NFS19_all_0.5_0.03_results_meta.csv")
remove <- c(which(colnames(test@gt) %in% metadata$id == FALSE))[-1] 
keep <- c(1:8,9:(ncol(test@gt)+8))
keep <- keep[-(remove+8)]
test2 <- read.vcfR("./NFS19_parentage_dp5_miss0.90_maf0.01_thin10k.recode.vcf",cols=keep)
colnames(test2@gt) <- sapply(colnames(test2@gt),function(x) gsub(".*aln_(.*).sorted.bam","\\1",x))
write.vcf(test2,file="./NFS19_0.9_0.01_dp5_thin10k_filt.vcf.gz")
library(LEA)
vcf2geno("./NFS19_0.9_0.01_dp5_thin10k_filt.vcf",output.file = "./NFS19_parentage_0.9_0.01_dp5_thin10k_filt.geno") #geno is easier to convert to hiphop
test <- read.geno("./NFS19_parentage_0.9_0.01_dp5_thin10k_filt.geno")
missing <- apply(test,1,function(x) length(which(x == 9))/length(x))
hist(missing) #fantastic!

#create input file(s)
library(LEA)
dat2 <- read.geno("./NFS19_parentage_0.9_0.01_dp5_thin10k_filt.geno")
metadata <- read.csv("./../entropy/all_samp/NFS19_all_0.5_0.03_results_meta.csv")

geno <- read.geno("./NFS19_parentage_0.9_0.01_dp5_thin10k_filt.geno")

geno <- apply(geno,MARGIN = c(1,2),function(x) ifelse(x==9,NA,x)) #different missing
geno <- apply(geno,MARGIN = c(1,2),function(x) ifelse(x==1,2,ifelse(x==2,1,x))) #hiphop says flip hetero/homozygotes
ind <- as.data.frame(matrix(nrow=nrow(geno),ncol=5)) #indvidual metadata file
colnames(ind) <- c("brood","individual","type","social.parent","year")
ind$brood <- rep(1,nrow(ind))
ind$individual <- metadata$id
type <- metadata$sex
type <- gsub("Male","adult male",type)
type <- gsub("Female","adult female",type)
type <- sapply(type,function(x) ifelse(is.na(x),"offspring",x))
ind$type <- type
ind$social.parent <- rep(0,nrow(ind))
ind$year <- 2019
rownames(geno) <- ind$individual

#split into 10 groups but with all adults in each group
adults <- which(!is.na(metadata$sex))
juvs <- which(is.na(metadata$sex))
increment <- round(length(juvs)/10)
breaks <- seq(1,length(juvs),by=increment)
juv1 <-juvs[c(breaks[1]:breaks[2])]
juv2 <- juvs[c((breaks[2]+1):breaks[3])]
juv3 <- juvs[c((breaks[3]+1):breaks[4])]
juv4 <- juvs[c((breaks[4]+1):breaks[5])]
juv5 <- juvs[c((breaks[5]+1):breaks[6])]
juv6 <- juvs[c((breaks[6]+1):breaks[7])]
juv7 <- juvs[c((breaks[7]+1):breaks[8])]
juv8 <- juvs[c((breaks[8]+1):breaks[9])]
juv9 <- juvs[c((breaks[9]+1):breaks[10])]
juv10 <- juvs[c((breaks[10]+1):(breaks[11]+1))]


geno1 <- geno[c(adults,juv1),]
geno2 <- geno[c(adults,juv2),]
geno3 <- geno[c(adults,juv3),]
geno4 <- geno[c(adults,juv4),]
geno5 <- geno[c(adults,juv5),]
geno6 <- geno[c(adults,juv6),]
geno7 <- geno[c(adults,juv7),]
geno8 <- geno[c(adults,juv8),]
geno9 <- geno[c(adults,juv9),]
geno10 <- geno[c(adults,juv10),]

ind1 <- ind[c(adults,juv1),]
ind2 <- ind[c(adults,juv2),]
ind3 <- ind[c(adults,juv3),]
ind4 <- ind[c(adults,juv4),]
ind5 <- ind[c(adults,juv5),]
ind6 <- ind[c(adults,juv6),]
ind7 <- ind[c(adults,juv7),]
ind8 <- ind[c(adults,juv8),]
ind9 <- ind[c(adults,juv9),]
ind10 <- ind[c(adults,juv10),]

write.csv(geno1,file="./hiphop_files/ind1_geno_dp5_0.9_0.01_thin10k.csv",row.names = F,quote=F)
write.csv(geno2,file="./hiphop_files/ind2_geno_dp5_0.9_0.01_thin10k.csv",row.names = F,quote=F)
write.csv(geno3,file="./hiphop_files/ind3_geno_dp5_0.9_0.01_thin10k.csv",row.names = F,quote=F)
write.csv(geno4,file="./hiphop_files/ind4_geno_dp5_0.9_0.01_thin10k.csv",row.names = F,quote=F)
write.csv(geno5,file="./hiphop_files/ind5_geno_dp5_0.9_0.01_thin10k.csv",row.names = F,quote=F)
write.csv(geno6,file="./hiphop_files/ind6_geno_dp5_0.9_0.01_thin10k.csv",row.names = F,quote=F)
write.csv(geno7,file="./hiphop_files/ind7_geno_dp5_0.9_0.01_thin10k.csv",row.names = F,quote=F)
write.csv(geno8,file="./hiphop_files/ind8_geno_dp5_0.9_0.01_thin10k.csv",row.names = F,quote=F)
write.csv(geno9,file="./hiphop_files/ind9_geno_dp5_0.9_0.01_thin10k.csv",row.names = F,quote=F)
write.csv(geno10,file="./hiphop_files/ind10_geno_dp5_0.9_0.01_thin10k.csv",row.names = F,quote=F)

write.csv(ind1,file="./hiphop_files/ind1_ind_big.csv",row.names = F,quote=F)
write.csv(ind2,file="./hiphop_files/ind2_ind_big.csv",row.names = F,quote=F)
write.csv(ind3,file="./hiphop_files/ind3_ind_big.csv",row.names = F,quote=F)
write.csv(ind4,file="./hiphop_files/ind4_ind_big.csv",row.names = F,quote=F)
write.csv(ind5,file="./hiphop_files/ind5_ind_big.csv",row.names = F,quote=F)
write.csv(ind6,file="./hiphop_files/ind6_ind_big.csv",row.names = F,quote=F)
write.csv(ind7,file="./hiphop_files/ind7_ind_big.csv",row.names = F,quote=F)
write.csv(ind8,file="./hiphop_files/ind8_ind_big.csv",row.names = F,quote=F)
write.csv(ind9,file="./hiphop_files/ind9_ind_big.csv",row.names = F,quote=F)
write.csv(ind10,file="./hiphop_files/ind30_ind_big.csv",row.names = F,quote=F)


#read in results from running on Teton
res1 <- read.csv("./hiphop_files/hiphop_results_hothiphop_0.9_0.15_0.01_ld_filt_ind1.csv")
meta1 <- read.csv("./hiphop_files/ind1_ind_old_sequoia.csv")
res2 <- read.csv("./hiphop_files/hiphop_results_hothiphop_0.9_0.15_0.01_ld_filt_ind2.csv")
meta2 <- read.csv("./hiphop_files/ind2_ind_old_sequoia.csv")
res3 <- read.csv("./hiphop_files/hiphop_results_hothiphop_0.9_0.15_0.01_ld_filt_ind3.csv")
meta3 <- read.csv("./hiphop_files/ind3_ind_old_sequoia.csv")
res4 <- read.csv("./hiphop_files/hiphop_results_hothiphop_0.9_0.15_0.01_ld_filt_ind4.csv")
meta4 <- read.csv("./hiphop_files/ind4_ind_old_sequoia.csv")
res5 <- read.csv("./hiphop_files/hiphop_results_hothiphop_0.9_0.15_0.01_ld_filt_ind5.csv")
meta5 <- read.csv("./hiphop_files/ind5_ind_old_sequoia.csv")
res6 <- read.csv("./hiphop_files/hiphop_results_hothiphop_0.9_0.15_0.01_ld_filt_ind6.csv")
meta6 <- read.csv("./hiphop_files/ind6_ind_old_sequoia.csv")
res7 <- read.csv("./hiphop_files/hiphop_results_hothiphop_0.9_0.15_0.01_ld_filt_ind7.csv")
meta7 <- read.csv("./hiphop_files/ind7_ind_old_sequoia.csv")
res8 <- read.csv("./hiphop_files/hiphop_results_hothiphop_0.9_0.15_0.01_ld_filt_ind8.csv")
meta8 <- read.csv("./hiphop_files/ind8_ind_old_sequoia.csv")
res9 <- read.csv("./hiphop_files/hiphop_results_hothiphop_0.9_0.15_0.01_ld_filt_ind9.csv")
meta9 <- read.csv("./hiphop_files/ind9_ind_old_sequoia.csv")
res10 <- read.csv("./hiphop_files/hiphop_results_hothiphop_0.9_0.15_0.01_ld_filt_ind10.csv")
meta10 <- read.csv("./hiphop_files/ind10_ind_old_sequoia.csv")

results <- rbind(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10)
results <- results[-which(is.na(results$rank)),]
write.csv(results,file="hiphop_results_0.9_0.15_0.01_ld_filt_hothiphop.csv",row.names = F,quote=F)

 
#examine results
library(ggplot2)
metadata <- read.csv("./../entropy/all_samp/NFS19_all_0.5_0.03_results_meta.csv")
results_top2 <- read.csv("./hiphop_results_0.9_0.15_0.01_ld_filt_hothiphop.csv")
#results_top2 <- results_top2[-which(is.na(results_top2$rank)),]
results_top2$rank <- as.factor(results_top2$rank)

ggplot(results_top2,aes(x=hothiphop.parents, fill=rank,color=rank)) + geom_histogram(binwidth = 1) + theme_minimal()
length(which(results_top2$hothiphop.parents <= 3 & results_top2$rank == 1))

ggplot(results_top2,aes(x=hiphop, fill=rank,color=rank)) + geom_histogram(binwidth = 1) + theme_minimal() + geom_vline(xintercept = 20)
length(which(results_top2$hiphop <=20 & results_top2$rank==1))

ggplot(results_top2,aes(x=hothiphop.dam,y=hothiphop.sire, color=rank)) + geom_point() + theme_minimal() +
  geom_hline(yintercept = 13) + geom_vline(xintercept=12) + geom_abline(slope=1,intercept = 0)
length(which(results_top2$hothiphop.dam < 12 | results_top2$hothiphop.sire < 13))

ggplot(results_top2,aes(x=hot.dam,y=hot.sire, color=rank)) + geom_point() + theme_minimal() +
  geom_hline(yintercept = 8) + geom_vline(xintercept=8)
length(which(results_top2$rank == 1 & (results_top2$hot.dam <=8 | results_top2$hot.sire <=8)))


length(which(results_top2$rank == 1 & (results_top2$hot.dam <=8 | results_top2$hot.sire <=8 | results_top2$hothiphop.parents <=3)))
#look at first and second rank differences per offspring
offspring <- unique(results_top2$offspring)
top_dam <- c()
top_sire <- c()
second_dam <- c()
second_sire <- c()
top_hiphop <- c()
top_hothiphop <- c()
second_hiphop <- c()
second_hothiphop <- c()
top_hotdam <- c()
top_hotsire <- c()
second_hotdam <- c()
second_hotsire <- c()
top_hotparents <- c()
second_hotparents <- c()
top_loci <- c()
for (o in offspring){
  top_dam <- c(top_dam,results_top2$dam[which(results_top2$rank==1 & results_top2$offspring==o)])
  top_sire <- c(top_sire,results_top2$sire[which(results_top2$rank==1 & results_top2$offspring==o)])
  second_dam <- c(second_dam,results_top2$dam[which(results_top2$rank==2 & results_top2$offspring==o)])
  second_sire <- c(second_sire,results_top2$sire[which(results_top2$rank==2 & results_top2$offspring==o)])
  top_hiphop <- c(top_hiphop,results_top2$hiphop[which(results_top2$rank==1 & results_top2$offspring==o)])
  top_hothiphop <- c(top_hothiphop,results_top2$hothiphop.parents[which(results_top2$rank==1 & results_top2$offspring==o)])
  second_hiphop <- c(second_hiphop,results_top2$hiphop[which(results_top2$rank==2 & results_top2$offspring==o)])
  second_hothiphop <- c(second_hothiphop,results_top2$hothiphop.parents[which(results_top2$rank==2 & results_top2$offspring==o)])
  top_hotparents <- c(top_hotparents,results_top2$hot.parents[which(results_top2$rank==1 & results_top2$offspring==o)])
  second_hotparents <- c(second_hotparents,results_top2$hot.parents[which(results_top2$rank==2 & results_top2$offspring==o)])
  top_hotdam <- c(top_hotdam,results_top2$hot.dam[which(results_top2$rank==1 & results_top2$offspring==o)])
  top_hotsire <- c(top_hotsire,results_top2$hot.sire[which(results_top2$rank==1 & results_top2$offspring==o)])
  second_hotdam <- c(second_hotdam,results_top2$hot.dam[which(results_top2$rank==2 & results_top2$offspring==o)])
  second_hotsire <- c(second_hotsire,results_top2$hot.sire[which(results_top2$rank==2 & results_top2$offspring==o)])
  top_loci <- c(top_loci,results_top2$loci.triad[which(results_top2$rank==1 & results_top2$offspring == o)])
}
per_offspring_results <- data.frame(offspring,top_dam,top_sire,second_dam,second_sire,top_hiphop,top_hothiphop,
                                    second_hiphop,top_hotparents,second_hotparents,second_hothiphop,top_hotdam,top_hotsire,
                                    second_hotdam,second_hotsire,top_loci)
per_offspring_results$topdam_q <- metadata$q[match(per_offspring_results$top_dam,metadata$id)]
per_offspring_results$topsire_q <- metadata$q[match(per_offspring_results$top_sire,metadata$id)]
per_offspring_results$offspring_q <- metadata$q[match(per_offspring_results$offspring,metadata$id)]


messed_up <- rep(FALSE,nrow(per_offspring_results))
for (i in 1:nrow(per_offspring_results)){
  parent_qs <- per_offspring_results[i,c(17,18)]
  min_parent <- parent_qs[which.min(parent_qs)]
  max_parent <- parent_qs[which.max(parent_qs)]
  
  if (per_offspring_results$offspring_q[i] > max_parent | per_offspring_results$offspring_q[i] < min_parent){
    messed_up[i] <- TRUE
  }
}
per_offspring_results$messed_up <- messed_up
ggplot(per_offspring_results, mapping=aes(x=messed_up,y=top_loci)) + geom_bar(stat="identity")
messed_up <- per_offspring_results[which(messed_up==TRUE & (per_offspring_results$top_hothiphop <=50 | per_offspring_results$top_dam < 8 | per_offspring_results$top_sire < 8)),]
library(tidyverse)
messed_up_plotdf <- pivot_longer(messed_up,cols = c(offspring_q,topdam_q,topsire_q),names_to="category",values_to="q")
#ggplot(messed_up_plotdf, mapping=aes(x=index,y=q,color=category)) + geom_point() 
#none look that bad...


#compile results worth keeping
offspring <- per_offspring_results$offspring[which(per_offspring_results$top_hiphop <= 3 | per_offspring_results$top_hotdam < 9 | per_offspring_results$top_hotsire < 9)]
sire <- rep(NA,length(offspring))
dam <- rep(NA,length(offspring))
for (o in 1:length(offspring)){
  if (per_offspring_results$top_hotdam[which(per_offspring_results$offspring == offspring[o])] < 9){
    dam[o] <- per_offspring_results$top_dam[which(per_offspring_results$offspring == offspring[o])]
  }
  if (per_offspring_results$top_hotsire[which(per_offspring_results$offspring == offspring[o])] < 9){
    sire[o] <- per_offspring_results$top_sire[which(per_offspring_results$offspring == offspring[o])]
  }
  if (per_offspring_results$top_hiphop[which(per_offspring_results$offspring == offspring[o])] <= 3){
    if (is.na(sire[o])){
      sire[o] <- per_offspring_results$top_sire[which(per_offspring_results$offspring == offspring[o])]
    } 
    if (is.na(dam[o])){
      dam[o] <- per_offspring_results$top_dam[which(per_offspring_results$offspring == offspring[o])]
    }
  }
}
good_assignments <- data.frame(offspring,dam,sire)
females <- table(good_assignments$dam)
males <- table(good_assignments$sire)
female_output <- data.frame("ID"=names(females),"num_offspring"=as.numeric(females),"q"=rep(NA,length(females)))
female_output$q <- metadata$q[match(female_output$ID,metadata$id)]
male_output <- data.frame("ID"=names(males),"num_offspring"=as.numeric(males),"q"=rep(NA,length(males)))
male_output$q <- metadata$q[match(male_output$ID,metadata$id)]

#close-kin mark-recapture estimate of Trout Creek pop size
metadata$tributary <- gsub("Trout ","Trout",metadata$tributary)
adults <- length(which(!is.na(metadata$sex) & metadata$tributary == "Trout"))
juvs <- length(which(is.na(metadata$sex) & metadata$tributary == "Trout"))
assignments <- sum(length(which(!is.na(good_assignments$dam))),length(which(!is.na(good_assignments$sire))))
ckmr_est <- (2 * adults * juvs)/(assignments+1)


#logistic regression for prob(parent assigned) ~ juv q
library(lme4)
juv_q <- metadata$q[which(is.na(metadata$sex))]
juv_names <- metadata$id[which(is.na(metadata$sex))]
juv_parent_sampled <- ifelse(juv_names %in% good_assignments$offspring,1,0)
log_model <- glm(juv_parent_sampled ~ juv_q,family=binomial(link="logit"))
summary(log_model)


plot(male_output$num_offspring~male_output$q)
plot(female_output$num_offspring~female_output$q)

sequoia <- read.csv("./results_meta_offspring.csv")
length(which(sequoia$sex == "Female" & sequoia$num_offspring !=0))
length(which(sequoia$sex == "Male" & sequoia$num_offspring !=0))

indv <- unique(c(sequoia$id[which(sequoia$num_offspring !=0)],female_output$ID,male_output$ID))
sequoia_num <- sequoia$num_offspring[match(indv,sequoia$id)]
hiphop_num <- rbind(male_output,female_output)$num_offspring[match(indv,rbind(male_output,female_output)$ID)]
results_comparison <- data.frame(indv,sequoia_num,hiphop_num)
results_comparison$hiphop_num <- sapply(results_comparison$hiphop_num,function(x) ifelse(is.na(x),0,x))
results_comparison$sequoia_num <- sapply(results_comparison$sequoia_num,function(x) ifelse(is.na(x),0,x))
results_comparison <- pivot_longer(results_comparison,!indv,names_to="method",values_to = "num_offspring")
results_comparison$sex <- metadata$sex[match(results_comparison$indv,metadata$id)]
ggplot(results_comparison,mapping=aes(x=method, y=num_offspring, color=sex)) + geom_point() + geom_line(aes(group=indv)) +
  theme_minimal() + xlab("Method") + ylab("Number of offspring")

#number of fully-represented trios
length(which(!is.na(good_assignments$dam) & !is.na(good_assignments$sire)))
#about the same as sequoia, which is good.

#write new file with new offspring numbers
new_results <- data.frame(indv,hiphop_num)
metadata$hiphop_offspring <- new_results$hiphop_num[match(metadata$id,new_results$indv)]
write.csv(metadata,"./results_meta_hiphop_offspring.csv",quote = F,row.names = F)





#compare to sequoia with same data
library(tidyverse)
metadata <- read.csv("./results_meta_hiphop_offspring.csv")
load("./parentage_new_0.9.RData")
females <- table(parentage$PedigreePar$dam)
males <- table(parentage$PedigreePar$sire)
female_output <- data.frame("ID"=names(females),"num_offspring"=as.numeric(females),"q"=rep(NA,length(females)))
female_output$q <- metadata$q[match(female_output$ID,metadata$id)]
male_output <- data.frame("ID"=names(males),"num_offspring"=as.numeric(males),"q"=rep(NA,length(males)))
male_output$q <- metadata$q[match(male_output$ID,metadata$id)]

indv <- unique(c(metadata$id[which(metadata$hiphop_offspring !=0)],female_output$ID,male_output$ID))
sequoia <- rbind(male_output,female_output)$num_offspring[match(indv,rbind(male_output,female_output)$ID)]
hiphop_num <- metadata$hiphop_offspring[match(indv,metadata$id)]

results_comparison <- data.frame(indv,sequoia,hiphop_num)
results_comparison$hiphop_num <- sapply(results_comparison$hiphop_num,function(x) ifelse(is.na(x),0,x))
results_comparison$sequoia <- sapply(results_comparison$sequoia,function(x) ifelse(is.na(x),0,x))
results_comparison$diff <- results_comparison$hiphop_num - results_comparison$sequoia
results_comparison <- results_comparison[order(abs(results_comparison$diff),decreasing = T),]
results_comparison$sex <- metadata$sex[match(results_comparison$indv,metadata$id)]

results_comparison <- pivot_longer(results_comparison,!indv,names_to="method",values_to = "num_offspring")
ggplot(results_comparison,mapping=aes(x=method, y=num_offspring, color=sex)) + geom_point() + geom_line(aes(group=indv)) +
  theme_minimal() + xlab("Method") + ylab("Number of offspring")

#look at number of shared and not-shared parent-offspring relationships in sequoia & hiphop
hiphop_assignments <- read.csv("./hiphop_good_assignments.csv")
sequoia_assignments <- parentage$PedigreePar[-c(1:553),]
shared <- 0
hiphop_only <- 0
sequoia_only <- 0

for (i in 1:nrow(sequoia_assignments)){
  offspring <- sequoia_assignments$id[i]
  seq_dam <- sequoia_assignments$dam[i]
  seq_sire <- sequoia_assignments$sire[i]
  
  hiphop_row <- which(hiphop_assignments$offspring == offspring)
  if (length(hiphop_row) == 0){
    if (!is.na(seq_dam)){
      sequoia_only <- sequoia_only+1
    }
    if (!is.na(seq_sire)){
      sequoia_only <- sequoia_only+1
    }
    next
  }
  hh_dam <- hiphop_assignments$dam[hiphop_row]
  hh_sire <- hiphop_assignments$sire[hiphop_row]  
  
  if (!is.na(hh_dam)){
    if (!is.na(seq_dam)){
      if (hh_dam == seq_dam){
        shared <- shared+1
      } else {
        hiphop_only <- hiphop_only+1
        sequoia_only <- sequoia_only+1
      }
    } else {
      hiphop_only <- hiphop_only+1
    }
  } else {
    if (!is.na(seq_dam)){
      sequoia_only <- sequoia_only+1
    }
  }
  
  if (!is.na(hh_sire)){
    if (!is.na(seq_sire)){
      if (hh_sire == seq_sire){
        shared <- shared+1
      } else {
        hiphop_only <- hiphop_only+1
        sequoia_only <- sequoia_only+1
      }
    } else {
      hiphop_only <- hiphop_only+1
    }
  } else {
    if (!is.na(seq_sire)){
      sequoia_only <- sequoia_only+1
    }
  }
  
}

numbers_df <- data.frame(category=c("hiphop only","shared","sequoia only"), count = c(hiphop_only,shared,sequoia_only))
ggplot(numbers_df,mapping=aes(x=category,y=count)) + geom_bar(stat="identity") + theme_minimal()
#327 total for hiphop, 300 for sequoia

