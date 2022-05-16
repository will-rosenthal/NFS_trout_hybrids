#code for visualizing trout matings & q
#3/2/2021
library(ggplot2)
library(tidyverse)
library(gghalves)
library(ggExtra)
library(gridExtra)
library(RGraphics)

meta <- read.csv("./results_meta_offspring.csv")
meta <- meta[which(meta$sex %in% c("Male","Female")),]
meta$sex <- factor(meta$sex,levels=c("Male","Female"))
matings <- read.csv("./observed_matings_hiphop.csv")
matings <- unique(matings)

meta$jit <- jitter(as.integer(meta$sex),amount=0.03)

matings$id <- rownames(matings)
matings <- gather(matings,"dam","sire",key="sex",value="name")
matings$q <- meta$q[match(matings$name,meta$id)]
matings$jit <- meta$jit[match(matings$name,meta$id)]

pal <- rev(c("#f16d22","#ffc629","#becb70","#7dcfb6","#00b2ca","#1988af","#1d4e89"))


p <- ggplot(data=meta, aes(x=sex)) + 
  geom_point(aes(x=jit,y=q,group=sex, color=q), size=3, show.legend = FALSE) +
  scale_color_gradientn(colors=pal) +
  geom_half_violin(data=meta[which(meta$sex == "Male"),],position=position_nudge(x=-0.05), fill="gray60",aes(y=q, x=as.integer(sex), alpha=0.6), side="l", show.legend = FALSE) + 
  geom_half_violin(data=meta[which(meta$sex == "Female"),],position=position_nudge(x=0.05), fill="gray60", aes(y=q, x=as.integer(sex), alpha=0.6), side="r", show.legend = FALSE) +
  theme_minimal() + theme(axis.title = element_text(size=15),axis.text = element_text(size=14)) +
  geom_line(data=matings,aes(x=jit,y=q,group=id),color="gray40") +
  scale_x_continuous(name="Sex", breaks=c(1,2), labels=c(paste0("Male\n","(n=",as.integer(table(meta$sex)[1]),")"),paste0("Female\n","(n=",as.integer(table(meta$sex)[2]),")"))) +
  ylab("Prop. YCT ancestry (q)") + removeGridX() #+
  #geom_segment(aes(x=2.4,xend=2.6,y=0,yend=0),color="gray20",size=1.3) +
  #geom_segment(aes(x=2.4,xend=2.6,y=0.1,yend=0.1),color="gray20",size=1.3) +
  #geom_segment(aes(x=2.4,xend=2.6,y=0.4,yend=0.4),color="gray20",size=1.3) +
  #geom_segment(aes(x=2.4,xend=2.6,y=0.6,yend=0.6),color="gray20",size=1.3) +
  #geom_segment(aes(x=2.4,xend=2.6,y=0.9,yend=0.9),color="gray20",size=1.3) +
  #geom_segment(aes(x=2.4,xend=2.6,y=1,yend=1),color="gray20",size=1.3)

yct_df <- matings[which(matings$id %in% matings$id[which(matings$sex == "dam" & matings$q > 0.9)] & matings$sex == "sire"),]
hyb_df1 <- matings[which(matings$id %in% matings$id[which(matings$sex == "dam" & matings$q < 0.4 & matings$q > 0.1)] & matings$sex == "sire"),]
hyb_df2 <- matings[which(matings$id %in% matings$id[which(matings$sex == "dam" & matings$q > 0.4 & matings$q < 0.6)] & matings$sex == "sire"),]
hyb_df3 <- matings[which(matings$id %in% matings$id[which(matings$sex == "dam" & matings$q > 0.6 & matings$q < 0.9)] & matings$sex == "sire"),]
rbt_df <- matings[which(matings$id %in% matings$id[which(matings$sex == "dam" & matings$q < 0.1)] & matings$sex == "sire"),]

{p2 <- ggplot(yct_df,aes(x=q)) + geom_histogram(fill="gray50",binwidth = 0.1) + 
  scale_x_continuous(breaks=seq(0,1,by=0.25),limits = c(-0.05,1.05)) + 
  theme_minimal() + theme(panel.grid.minor.y = element_blank()) + scale_y_continuous(breaks=c(0,5,10)) +
  labs(y="Count",x="Male q",title="Mates for YCT females") + theme(axis.title = element_text(size=12),
                                                                   plot.title=element_text(size=13))
p3 <- ggplot(hyb_df1,aes(x=q)) + geom_histogram(fill="gray50",binwidth = 0.1) + 
  scale_x_continuous(breaks=seq(0,1,by=0.25),limits=c(-0.05,1.05)) + 
  theme_minimal() + theme(panel.grid.minor.y = element_blank()) + scale_y_continuous(breaks=c(0,2,4,6)) +
  labs(y="Count",x="Male q",title="Mates for females 0.1 < q < 0.4") + theme(axis.title = element_text(size=12),
                                                                   plot.title=element_text(size=13))

p4 <- ggplot(hyb_df2,aes(x=q)) + geom_histogram(fill="gray50",binwidth = 0.1) + 
  scale_x_continuous(breaks=seq(0,1,by=0.25),limits=c(-0.05,1.05)) + 
  theme_minimal() + theme(panel.grid.minor.y = element_blank()) + scale_y_continuous(breaks=c(0,5,10,15)) +
  labs(y="Count",x="Male q",title="Mates for females 0.4 < q < 0.6") + theme(axis.title = element_text(size=12),
                                                                             plot.title=element_text(size=13))
p5 <- ggplot(hyb_df3,aes(x=q)) + geom_histogram(fill="gray50",binwidth = 0.1) + 
  scale_x_continuous(breaks=seq(0,1,by=0.25),limits=c(-0.05,1.05)) + 
  theme_minimal() + theme(panel.grid.minor.y = element_blank()) + scale_y_continuous(breaks=c(0,2,4)) +
  labs(y="Count",x="Male q",title="Mates for females 0.6 < q < 0.9") + theme(axis.title = element_text(size=12),
                                                                             plot.title=element_text(size=13))

p6 <- ggplot(rbt_df,aes(x=q)) + geom_histogram(fill="gray50",binwidth = 0.1) + 
  scale_x_continuous(breaks=seq(0,1,by=0.25),limits=c(-0.05,1.05)) + 
  theme_minimal() + theme(panel.grid.minor.y = element_blank()) + scale_y_continuous(breaks=c(0,2,4,6)) +
  labs(y="Count",x="Male q",title="Mates for RBT females") + theme(axis.title = element_text(size=12),
                                                                             plot.title=element_text(size=13))
}
#grid.arrange(p,p2,p3,p4,p5,p6,textGrob(""),layout_matrix=cbind(c(1,1),c(6,3),c(4,5),c(2,7)),
#             widths=c(4.5,1.5,1.5,1.5))

grid.arrange(grobs=list(p,p2,p3,p4,p5,p6,textGrob("")),layout_matrix=cbind(c(1,1,1,1,1,1),c(1,1,1,1,1,1),c(2,5,4,3,6,7)),
             widths=c(2.7,2.7,2.5),heights=c(1,1,1,1,1,0.32))











