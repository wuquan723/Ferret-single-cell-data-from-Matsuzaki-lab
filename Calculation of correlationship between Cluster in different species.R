####This R script is for calculation of correlation of two species used in Figure 3C and 3D
###Loading data from previous step: "scRNA for CCA analysis for release"
Human_Sci_NP <- readRDS(file = "~/Human_Sci_NP.rds")
Ferret_NP <- readRDS(file = "~/Ferret_NP.rds")

#### Human markers####
Human_NP.markers <- FindAllMarkers(Human_Sci_NP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Human_NP.markers<-Human_NP.markers[Human_NP.markers$p_val_adj<0.01,]
Ferret_NP.markers <- FindAllMarkers(Ferret_NP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Ferret_NP.markers<-Ferret_NP.markers[Ferret_NP.markers$p_val_adj<0.01,]

####Here Ref_dataset is an object of seurat. Ref_markers is made by FindAllMarkers function in seurat using Ref_dataset
####Query marker is made by FindAllMarkers function in seurat using Query_dataset
####The example of these data files can be download here
Search_for_cluster<- function(Query_Cluster,Ref_dataset,Query_markers,Ref_markers)
{
  correlationship<-NULL
  correlationship1<-NULL
  for (i in 1:length(levels(Ref_dataset))) 
  {
    Query_markers_tmp<-Query_markers[Query_markers$cluster==Query_Cluster,]##Find query cluster and calculate gene score
    Ref_markers_tmp<-Ref_markers[Ref_markers$cluster==levels(Ref_dataset)[i],]##Compare query with Ref_data one by one
    Query_markers_tmp<-Query_markers_tmp[is.element(Query_markers_tmp$gene,
                                                    intersect(Query_markers_tmp$gene,Ref_markers_tmp$gene)),]
    Query_markers_tmp$Query_genescore<-apply(Query_markers_tmp,1,function(x){as.numeric(x[2])*as.numeric(x[3])/as.numeric(x[4])})
    
    Ref_markers_tmp<-Ref_markers_tmp[is.element(Ref_markers_tmp$gene,
                                                intersect(Query_markers_tmp$gene,Ref_markers_tmp$gene)),]
    Ref_markers_tmp$Ref_genescore<-apply(Ref_markers_tmp,1,function(x){as.numeric(x[2])*as.numeric(x[3])/as.numeric(x[4])})
    Query_markers_tmp$Ref_genescore<-Ref_markers_tmp$Ref_genescore[match(Query_markers_tmp$gene,Ref_markers_tmp$gene)]
    try ({cor<-cor.test(Query_markers_tmp$Query_genescore, Query_markers_tmp$Ref_genescore, method=c("pearson"))})
    correlationship1<-c(cor$p.value,cor$estimate)
    correlationship<-rbind(correlationship,correlationship1)
  }
  correlationship<-as.data.frame(correlationship)
  rownames(correlationship)<-levels(Ref_dataset)
  colnames(correlationship)<-c("p_value",paste(Query_Cluster,"correlationship"))
  return(correlationship)
}
###We can calculate how Ferret tRG cells correlate with Human clusters like this:
Ferret_tRG<-Search_for_cluster("tRG",Human_Sci_NP,Ferret_NP.markers,Human_NP.markers)

