com_kegg_go <- function(gcSample,pro){
  xx <- compareCluster(gcSample, fun="enrichKEGG",
                       organism="hsa", pvalueCutoff=0.05)
  dotplot(xx) 
  ggsave(paste0(pro,'_kegg.pdf'))
  
  # Run full GO enrichment test for BP 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont		   = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_BP_cluster_simplified.pdf') ,width = 15,height = 6)
  print(dotplot(lineage1_ego, showCategory=5))
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_BP_cluster_simplified.csv'))
  # Run full GO enrichment test for CC 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont		   = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_CC_cluster_simplified.pdf') ,width = 15,height = 6)
  print(dotplot(lineage1_ego, showCategory=5))
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_CC_cluster_simplified.csv'))
  
  # Run full GO enrichment test for MF 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont		   = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_MF_cluster_simplified.pdf') ,width = 15,height = 6)
  print(dotplot(lineage1_ego, showCategory=5))
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_MF_cluster_simplified.csv'))
  
}