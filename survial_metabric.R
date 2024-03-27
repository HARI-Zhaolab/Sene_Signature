### 生存
rm(list = ls())
options(stringsAsFactors = F)

load("./data/data_metabric.Rdata")

data <- matedata[order(matedata[,"senescore"]),]
table(data$vitalStat)

data$surTime <- as.numeric(data$surTime)

n <- nrow(data)
rownames(data) <- 1:n
library(survminer)
library(survival)
pvals<-c()
hrs<-c()
for(n in 1:(nrow(data)- 1)){ 
  ssdf <- cbind(data, data.frame(exp = rep(c( "Low", "High"), c(n, nrow(data)-n)))) 
  ssdf$exp <- factor(ssdf$exp, levels = c( "Low", "High")) 
  fitd <- survdiff(Surv(surTime,vitalStat) ~exp , data=ssdf,na.action = na.exclude) 
  pValue <- 1- pchisq(fitd$chisq, length(fitd$n)- 1) 
  pvals <- c(pvals, pValue) 
  hr <- fitd$obs[ 1] * fitd$exp[ 2]/(fitd$obs[ 2] * fitd$exp[ 1]) 
  hrs <- c(hrs, hr) }

sur <- data.frame(num = 1:(nrow(data)-1),
                  HR = hrs,
                  pvalue = pvals)

#取出最佳截断值
sur <- sur[(n/4):(3*n/4),]
pvals <- sur$pvalue
bestNum <-sur$num[which.min(pvals)]
ssdf <- cbind(data, data.frame(exp = rep(c( "Low", "High"), c(bestNum, nrow(data)-bestNum))))

fit<-survfit(Surv(surTime,vitalStat) ~exp , data=ssdf)
# str(data)
# fit<-survfit(Surv(SurvivalTime,status) ~ pyroscore, data=merge)
# summary(fit)

# 计算HR 和 95%cl
# ssdf$exp <- factor(ssdf$exp, levels = c("Low","High"))
# cox <- coxph(Surv(surTime, vitalStat) ~exp, data = ssdf)
# coxSummary <- summary(cox)

# 不能把coxph和surdiff混用
ssdf <- ssdf[order(ssdf$exp,decreasing = T),]
ssdf$exp <- factor(ssdf$exp,levels = c("Low","High"))
data.survdiff <- survdiff(Surv(surTime,vitalStat) ~ exp,data = ssdf)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))

label <- paste0(ifelse(round(p.val,4) == 0, "p < 0.01", 
                       paste0("p = ",round(p.val,4))),"\n",
                "Hazard Ratio = ",round(HR,1),"\n",
                "95% Cl: =  ",round(low95,2),"-",
                round(up95,2))
p <- ggsurvplot(fit,pval = label, #show p-value of log-rank test，显示log-rank分析得到的P值
                conf.int = FALSE, #添加置信区间
                conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样子
                risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)的形式展示risk table
                risk.table.y.text.col = T,###  colour risk table text annotations.
                risk.table.y.text = FALSE,###  show bars instead of names in text annotations in legend of risk table.不显示注释名字
                xlab = "Time in months", ###  customize X axis label.自定义x的标签为time in years
                #surv.median.line = "hv", #添加中位生存时间的线
                ncensor.plot = FALSE, #我这里不显示删失的图，TRUE就显示
                egend.labs =c("high risk", "low risk"),    ###  对legend的标签重新命名
                legend.labs=c("high", "low"),
                legend=c(0.8,0.8),
                palette = c("#ff0000", "#000000"), ###  自定义颜色
                ggtheme =  theme(axis.title.x = element_text(size = 15),
                                 axis.title.y = element_text(size = 15),
                                 axis.text.x = element_text(size = 12),
                                 axis.text.y = element_text(size = 12),
                                 strip.text.x = element_text(size = 12),
                                 legend.position = "none",
                                 panel.background = element_rect(fill = "white",color = "black"),
                                 panel.grid = element_blank()), #绘图主题
                #xlim= c(0,100),
                #break.x.by=20,
                break.y.by=0.2
)
pdf(paste0("./outputfig/Fig4/3_survial_score_matebric.pdf"),width = 10)
print(p)
dev.off()


# 保存一下高低分组样本信息，虽然后面分析按照的是中位数
save(ssdf,file = "./outputdata/Fig4/ssdf_metabric.Rdata")

# 多计算一个RFS
rm(list = ls())
options(stringsAsFactors = F)

load("./data/data_metabric.Rdata")

data <- matedata[order(matedata[,"senescore"]),]
table(data$RFS_STATUS)

data$surTime <- as.numeric(data$RFS_MONTHS)
library(stringr)
data$vitalStat <- as.numeric(str_split(data$RFS_STATUS,":",simplify = T)[,1])

n <- nrow(data)
rownames(data) <- 1:n
library(survminer)
library(survival)
pvals<-c()
hrs<-c()
for(n in 1:(nrow(data)- 1)){ 
  ssdf <- cbind(data, data.frame(exp = rep(c( "Low", "High"), c(n, nrow(data)-n)))) 
  ssdf$exp <- factor(ssdf$exp, levels = c( "Low", "High")) 
  fitd <- survdiff(Surv(surTime,vitalStat) ~exp , data=ssdf,na.action = na.exclude) 
  pValue <- 1- pchisq(fitd$chisq, length(fitd$n)- 1) 
  pvals <- c(pvals, pValue) 
  hr <- fitd$obs[ 1] * fitd$exp[ 2]/(fitd$obs[ 2] * fitd$exp[ 1]) 
  hrs <- c(hrs, hr) }

sur <- data.frame(num = 1:(nrow(data)-1),
                  HR = hrs,
                  pvalue = pvals)

#取出最佳截断值
sur <- sur[(n/4):(3*n/4),]
pvals <- sur$pvalue
bestNum <-sur$num[which.min(pvals)]
ssdf <- cbind(data, data.frame(exp = rep(c( "Low", "High"), c(bestNum, nrow(data)-bestNum))))

fit<-survfit(Surv(surTime,vitalStat) ~exp , data=ssdf)
# str(data)
# fit<-survfit(Surv(SurvivalTime,status) ~ pyroscore, data=merge)
# summary(fit)

# 计算HR 和 95%cl
# ssdf$exp <- factor(ssdf$exp, levels = c("Low","High"))
# cox <- coxph(Surv(surTime, vitalStat) ~exp, data = ssdf)
# coxSummary <- summary(cox)

# 不能把coxph和surdiff混用
ssdf <- ssdf[order(ssdf$exp,decreasing = T),]
ssdf$exp <- factor(ssdf$exp,levels = c("Low","High"))
data.survdiff <- survdiff(Surv(surTime,vitalStat) ~ exp,data = ssdf)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))

label <- paste0(ifelse(round(p.val,4) == 0, "p < 0.01", 
                       paste0("p = ",round(p.val,4))),"\n",
                "Hazard Ratio = ",round(HR,1),"\n",
                "95% Cl: =  ",round(low95,2),"-",
                round(up95,2))
p <- ggsurvplot(fit,pval = label, #show p-value of log-rank test，显示log-rank分析得到的P值
                conf.int = FALSE, #添加置信区间
                conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样子
                risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)的形式展示risk table
                risk.table.y.text.col = T,###  colour risk table text annotations.
                risk.table.y.text = FALSE,###  show bars instead of names in text annotations in legend of risk table.不显示注释名字
                xlab = "Time in months", ###  customize X axis label.自定义x的标签为time in years
                #surv.median.line = "hv", #添加中位生存时间的线
                ncensor.plot = FALSE, #我这里不显示删失的图，TRUE就显示
                egend.labs =c("high risk", "low risk"),    ###  对legend的标签重新命名
                legend.labs=c("high", "low"),
                legend=c(0.8,0.8),
                palette = c("#ff0000", "#000000"), ###  自定义颜色
                ggtheme =  theme(axis.title.x = element_text(size = 15),
                                 axis.title.y = element_text(size = 15),
                                 axis.text.x = element_text(size = 12),
                                 axis.text.y = element_text(size = 12),
                                 strip.text.x = element_text(size = 12),
                                 legend.position = "none",
                                 panel.background = element_rect(fill = "white",color = "black"),
                                 panel.grid = element_blank()), #绘图主题
                #xlim= c(0,100),
                #break.x.by=20,
                break.y.by=0.2
)
pdf(paste0("./outputfig/Fig4/3_survial_score_matebric.pdf"),width = 10)
print(p)
dev.off()


