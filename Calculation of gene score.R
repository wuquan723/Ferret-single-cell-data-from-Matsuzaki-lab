####This R script is for calculation of oRG socre. Also you can change oRG with tRG to calculate tRG score
###Loading data from previous step: "scRNA for CCA analysis for release"

Human_Sci_NP <- readRDS(file = "~/Human_Sci_NP.rds")
Ferret_NP <- readRDS(file = "~/Ferret_NP.rds")
###Calculation oRG scores in new dataset--------------------
Human_NP.markers <- FindAllMarkers(Human_NP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Human_oRG_marker<-Human_NP.markers[Human_NP.markers$cluster=="oRG",]
Human_oRG_marker<-Human_oRG_marker[Human_oRG_marker$p_val_adj<0.01,]

Human_oRG_marker$genescore<-apply(Human_oRG_marker,1,function(x){as.numeric(x[2])*as.numeric(x[3])/as.numeric(x[4])})


#oRG score in Ferret_RG
Ferret_RG.raw <-GetAssayData(object = Ferret_NP, slot = "data")
Ferret_RG.raw_subtRG<-Ferret_RG.raw[is.element(rownames(Ferret_RG.raw),Human_oRG_marker$gene),]
Ferret_RG.raw_subtRG<-as.data.frame(Ferret_RG.raw_subtRG)
GeneScore_oRG <- Human_oRG_marker$genescore[match(rownames(Ferret_RG.raw_subtRG),Human_oRG_marker$gene)]

F_oRG_score<-apply(Ferret_RG.raw_subtRG,2,function(x){
  as.numeric(x)%*% as.numeric(GeneScore_oRG)
})



F_oRG_score<-as.data.frame(F_oRG_score)

Cellident_FerretRG <-as.data.frame(Idents(object = Ferret_NP))
colnames(Cellident_FerretRG)<-"Ident"
Cellident_FerretRG$Age<-Ferret_NP$Age
Cellident_FerretRG$Subtype<-Ferret_NP$Subtype
F_oRG_score$Ident<-Cellident_FerretRG$Ident[match(rownames(F_oRG_score),rownames(Cellident_FerretRG))]
F_oRG_score$Subtype<-Cellident_FerretRG$Subtype[match(rownames(F_oRG_score),rownames(Cellident_FerretRG))]
F_oRG_score$Age<-Cellident_FerretRG$Age[match(rownames(F_oRG_score),rownames(Cellident_FerretRG))]
colnames(F_oRG_score)<-c("oRG_score","ID","Subtype","Age")

F_oRG_score$anchor_oRG<-ifelse(is.element(rownames(F_oRG_score),New_ferret_anchor_cell_name),1,0)
library(beeswarm)
library(RColorBrewer)

beeswarm <- beeswarm(F_oRG_score$oRG_score ~ as.factor(F_oRG_score$anchor_oRG),cex=0.1,
                     data = F_oRG_score, method = 'swarm',ylim=c(-750, 1680),col=brewer.pal(1, "Set1"),
                     xlab="Cell type",
                     ylab="oRG Score")

boxplot(F_oRG_score$oRG_score ~ as.factor(F_oRG_score$anchor_oRG), boxwex=0.5, col=brewer.pal(1, "Set1"),
        data = F_oRG_score, add = F, at=c(1.5,2.5))

wilcox.test(F_oRG_score[F_oRG_score$anchor_oRG==1,]$oRG_score,F_oRG_score[F_oRG_score$anchor_oRG==0,]$oRG_score, data=F_oRG_score, paired = F)
##W = 859938, p-value < 2.2e-16