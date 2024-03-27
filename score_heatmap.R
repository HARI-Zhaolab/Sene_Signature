## 衰老评分与样本信息之间的关系
rm(list = ls())
options(stringsAsFactors = F)

load("./data/matedata.Rdata")
senescencescore_matedata <- matedata[order(matedata$senescore, decreasing = F),]
senescencescore_matedata$TNBC <- ifelse(senescencescore_matedata$IHC_sub == "TNBC","TNBC","no_TNBC")
senescencescore_matedata$TNBC[is.na(senescencescore_matedata$TNBC)] <- "no_TNBC"

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
display.brewer.all()
# Lund.subtype <- c(brewer.pal(9,"BuPu"),brewer.pal(11,"PuOr")[1])
# names(Lund.subtype) <- sort(unique(pyroscore_matedata$Lund.subtype))
stage.subtype <- c(brewer.pal(7,"GnBu")[1:5])
names(stage.subtype) <- sort(unique(senescencescore_matedata$pathologic_stage)[-4])
a <- HeatmapAnnotation(senescencescore = anno_barplot(senescencescore_matedata$senescore, bar_width = 1,col = "blue",
                                                      gp = gpar(col = "blue", fill = "blue"),border = FALSE),
                       age.subtype = senescencescore_matedata$age,
                       IHC.subtype = senescencescore_matedata$IHC_sub,
                       PAM50.subtype = senescencescore_matedata$PAM50_subtype,
                       TNBC.subtype = senescencescore_matedata$TNBC,
                       stage.subtype = senescencescore_matedata$pathologic_stage,
                       cluster.subtype = senescencescore_matedata$group_name,
                       col = list(pyroscore = "bule",
                                  IHC.subtype = c("luminal A" = brewer.pal(11,"Spectral")[6], 
                                                  "luminal B" = brewer.pal(11,"Spectral")[3],
                                                  "HER2" = brewer.pal(11,"Spectral")[3],
                                                  "TNBC"= brewer.pal(11,"Spectral")[1]),
                                  PAM50.subtype = c("Basal" = brewer.pal(11,"Spectral")[9],
                                                    "Normal" = brewer.pal(11,"Spectral")[10],
                                                    "LumA" = brewer.pal(11,"Spectral")[7],
                                                    "LumB" = brewer.pal(11,"Spectral")[6],
                                                    "Her2" = brewer.pal(11,"Spectral")[8]),
                                  cluster.subtype = c("cluster1" = brewer.pal(11,"RdYlGn")[5], 
                                                      "cluster2" = brewer.pal(11,"RdYlGn")[11]),
                                  TNBC.subtype = c("no_TNBC" = brewer.pal(8,"Reds")[2], 
                                                   "TNBC" = brewer.pal(8,"Reds")[7]),
                                  stage.subtype = stage.subtype),
                       show_legend = T)
pdf("./outputfig/Fig4/8_score_matedata_heatmap.pdf")
plot(a)
dev.off()

###### 吴宜鑫绘图代码
#################################################################
names(stage.subtype) <- sort(unique(senescencescore_matedata$pathologic_stage)[-4])
senescencescore_matedata$senescore <- scale(senescencescore_matedata$senescore)

anno_col <- list(age.subtype = colorRamp2(c(60,90), c("#FFFFFF", "#F39B7FCC")),
                 IHC.subtype = c(HER2 ="#E64B35FF",`luminal A` ="#91D1C2FF",`luminal B` = "#00A087FF", TNBC = "#DC0000FF"),
                 PAM50.subtype = c(Basal="#F39B7FFF",Her2="#E64B35FF",LumB="#00A087FF",LumA="#91D1C2FF",Normal="#42B540CC"),
                 cluster.subtype =c(cluster1="#B09C85CC",cluster2="#7E6148CC"),
                 TNBC.subtype = c(no_TNBC = "#F39B7FCC", TNBC = "#DC0000FF"),
                 stage.subtype = c(`Stage I`="#3C548899",`Stage II`="#3C5488B2",`Stage III`="#3C5488CC",`Stage IV`="#3C5488E5",`Stage X`="#3C5488FF"))




a <- HeatmapAnnotation(senescencescore = anno_barplot(senescencescore_matedata$senescore, bar_width = 1,
                                                      gp = gpar(col = "#E64B35CC", fill = "#E64B35CC"),border = FALSE),       
                       age.subtype = senescencescore_matedata$age,
                       stage.subtype = senescencescore_matedata$pathologic_stage,
                       IHC.subtype = senescencescore_matedata$IHC_sub,
                       PAM50.subtype = senescencescore_matedata$PAM50_subtype,
                       TNBC.subtype = senescencescore_matedata$TNBC,
                       cluster.subtype = senescencescore_matedata$group_name,
                       col = anno_col,
                       show_legend = T)

pdf("./outputfig_new/Fig4/8_score_matedata_heatmap.pdf",width = 8, height = 8,onefile = FALSE )
plot(a)
dev.off()
