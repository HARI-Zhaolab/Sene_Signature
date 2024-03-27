rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)#在调用as.data.frame的时，将stringsAsFactors设置为FALSE可以避免character类型自动转化为factor类型

load("./data_GSE2990.Rdata")
# 生存预后的影响
# 计算最佳截断值
data <- matedata[order(matedata[,"senescore"]),]
table(data$event.rfs)

data$surTime <- as.numeric(data$time.rfs)

n <- nrow(data)
rownames(data) <- 1:n
library(survminer)
library(survival)
pvals<-c()
hrs<-c()
for(n in 1:(nrow(data)- 1)){ 
  ssdf <- cbind(data, data.frame(exp = rep(c( "Low", "High"), c(n, nrow(data)-n)))) 
  ssdf$exp <- factor(ssdf$exp, levels = c( "Low", "High")) 
  fitd <- survdiff(Surv(surTime,event.rfs) ~exp , data=ssdf,na.action = na.exclude) 
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

fit<-survfit(Surv(surTime,event.rfs) ~exp , data=ssdf)
# str(data)
# fit<-survfit(Surv(SurvivalTime,status) ~ pyroscore, data=merge)
# summary(fit)

# 计算HR 和 95%cl
# ssdf$exp <- factor(ssdf$exp, levels = c("Low","High"))
# cox <- coxph(Surv(surTime, event.rfs) ~exp, data = ssdf)
# coxSummary <- summary(cox)
# 
# label <- paste0(ifelse(round(coxSummary$coefficients[,"Pr(>|z|)"],4) == 0, "p < 0.01", 
#                        paste0("p = ",round(coxSummary$coefficients[,"Pr(>|z|)"],4))),"\n",
#                 "Hazard Ratio = ",round(coxSummary$coefficients[,"exp(coef)"],1),"\n",
#                 "95% Cl: =  ",round(coxSummary$conf.int[,"lower .95"],2),"-",
#                 round(coxSummary$conf.int[,"upper .95"],2))

ssdf <- ssdf[order(ssdf$exp,decreasing = T),]
ssdf$exp <- factor(ssdf$exp,levels = c("Low","High"))
data.survdiff <- survdiff(Surv(surTime,event.rfs) ~ exp, data=ssdf)
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
                conf.int = TRUE, #添加置信区间
                conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样子
                risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)的形式展示risk table
                risk.table.y.text.col = T,###  colour risk table text annotations.
                risk.table.y.text = FALSE,###  show bars instead of names in text annotations in legend of risk table.不显示注释名字
                xlab = "Time in years", ###  customize X axis label.自定义x的标签为time in years
                #surv.median.line = "hv", #添加中位生存时间的线
                ncensor.plot = FALSE, #我这里不显示删失的图，TRUE就显示
                egend.labs =c("high risk", "low risk"),    ###  对legend的标签重新命名
                legend.labs=c("high", "low"),
                legend=c(0.9,0.9),
                size = 2,
                legend.title = "Strata",
                palette = c("#ED0000FF", "#1B1919FF"), ###  自定义颜色 #1B1919FF 黑；00468BFF，蓝；
                ggtheme =  theme(axis.title.x = element_text(size = 20),
                                 axis.title.y = element_text(size = 20),
                                 axis.text.x = element_text(size = 12),
                                 axis.text.y = element_text(size = 12),
                                 strip.text.x = element_text(size = 12),
                                 legend.text = element_text(size = 15),
                                 legend.title = element_text(size = 15),
                                 panel.background = element_rect(fill = "white",color = "black"),
                                 panel.grid = element_blank()), #绘图主题
                #xlim= c(0,100),
                break.x.by=5,
                break.y.by=0.2
)
p
pdf(paste0("./RFS_survial_score.pdf"),width = 10, height = 10, onefile = FALSE)
print(p)
dev.off()

# 计算最佳截断值
data <- matedata[order(matedata[,"senescore"]),]
table(data$event.dmfs)

data$surTime <- as.numeric(data$time.dmfs)

n <- nrow(data)
rownames(data) <- 1:n
library(survminer)
library(survival)
pvals<-c()
hrs<-c()
for(n in 1:(nrow(data)- 1)){ 
  ssdf <- cbind(data, data.frame(exp = rep(c( "Low", "High"), c(n, nrow(data)-n)))) 
  ssdf$exp <- factor(ssdf$exp, levels = c( "Low", "High")) 
  fitd <- survdiff(Surv(surTime,event.dmfs) ~exp , data=ssdf,na.action = na.exclude) 
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

fit<-survfit(Surv(surTime,event.dmfs) ~exp , data=ssdf)
# str(data)
# fit<-survfit(Surv(SurvivalTime,status) ~ pyroscore, data=merge)
# summary(fit)

# 计算HR 和 95%cl
# ssdf$exp <- factor(ssdf$exp, levels = c("Low","High"))
# cox <- coxph(Surv(surTime, event.dmfs) ~exp, data = ssdf)
# coxSummary <- summary(cox)
# 
# label <- paste0(ifelse(round(coxSummary$coefficients[,"Pr(>|z|)"],4) == 0, "p < 0.01", 
#                        paste0("p = ",round(coxSummary$coefficients[,"Pr(>|z|)"],4))),"\n",
#                 "Hazard Ratio = ",round(coxSummary$coefficients[,"exp(coef)"],1),"\n",
#                 "95% Cl: =  ",round(coxSummary$conf.int[,"lower .95"],2),"-",
#                 round(coxSummary$conf.int[,"upper .95"],2))

ssdf <- ssdf[order(ssdf$exp,decreasing = T),]
ssdf$exp <- factor(ssdf$exp,levels = c("Low","High"))
data.survdiff <- survdiff(Surv(surTime,event.dmfs) ~ exp, data=ssdf)
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
                conf.int = TRUE, #添加置信区间
                conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样子
                risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)的形式展示risk table
                risk.table.y.text.col = T,###  colour risk table text annotations.
                risk.table.y.text = FALSE,###  show bars instead of names in text annotations in legend of risk table.不显示注释名字
                xlab = "Time in years", ###  customize X axis label.自定义x的标签为time in years
                #surv.median.line = "hv", #添加中位生存时间的线
                ncensor.plot = FALSE, #我这里不显示删失的图，TRUE就显示
                egend.labs =c("high risk", "low risk"),    ###  对legend的标签重新命名
                legend.labs=c("high", "low"),
                legend=c(0.9,0.9),
                size = 2,
                legend.title = "Strata",
                palette = c("#ED0000FF", "#1B1919FF"), ###  自定义颜色 #1B1919FF 黑；00468BFF，蓝；
                ggtheme =  theme(axis.title.x = element_text(size = 20),
                                 axis.title.y = element_text(size = 20),
                                 axis.text.x = element_text(size = 12),
                                 axis.text.y = element_text(size = 12),
                                 strip.text.x = element_text(size = 12),
                                 legend.text = element_text(size = 15),
                                 legend.title = element_text(size = 15),
                                 panel.background = element_rect(fill = "white",color = "black"),
                                 panel.grid = element_blank()), #绘图主题
                #xlim= c(0,100),
                break.x.by=5,
                break.y.by=0.2
)
p
pdf(paste0("./DFS_survial_score.pdf"),width = 10, height = 10, onefile = FALSE)
print(p)
dev.off()

