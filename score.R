## 评分在不同cluster之间的差异
rm(list = ls())
options(stringsAsFactors = F)

load("./data/matedata.Rdata")

my_comparisons <- list(c("cluster1","cluster2"))
library(ggplot2)
library(ggpubr)
p <- ggplot(matedata,aes(x=factor(group_name),y=senescore,fill=factor(group_name), color=factor(group_name)))+
  geom_boxplot(size=0.8,outlier.fill="black",outlier.color="white",notch = TRUE)+
  geom_jitter(aes(fill=factor(group_name)),width =0.2,shape = 21,size=1.5)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  stat_compare_means(comparisons = my_comparisons,label="p.signif") +
  stat_compare_means(label.y =1,label.x =1.4) +
  scale_fill_manual(values = c("#CD3333", "#6CA6CD"))+
  scale_color_manual(values=c("black","black"))+
  ggtitle("")+
  theme_classic()+
  theme(legend.position = "None")+
  theme(axis.text.x=element_text(colour="black",family="Times",size=14),
        axis.text.y=element_text(family="Times",size=14,face="plain"),
        axis.title.y=element_text(family="Times",size = 14,face="plain"),
        axis.title.x=element_text(family="Times",size = 14,face="plain"),
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("senescore")+xlab("cluster")

ifelse(dir.exists("outputfig/Fig4"),print("文件夹已存在"),dir.create("outputfig/Fig4"))
ggsave(filename = "./outputfig/Fig4/1_cluster_score.pdf", plot= p)
