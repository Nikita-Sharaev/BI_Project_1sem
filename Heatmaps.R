library(iTALK)
library(dplyr)

#counts_nonresponders0 <- readRDS("pathway")
counts_responders0 <- readRDS("pathway")

#highly_exprs_genes<-rawParse(counts_nonresponders0,top_genes=50,stats='mean')

#highly_exprs_genes_nonresponders <- readRDS("pathway/highly_exprs_genes_nonresponders.rds")
highly_exprs_genes_responders <- readRDS("pathway/highly_exprs_genes_responders.rds")

comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_col<-structure(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a', '#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a'),names=unique(counts_responders0$cell_type))
par(mfrow=c(1,2))
res<-NULL
for(comm_type in comm_list){
res_cat<-FindLR(highly_exprs_genes_responders,datatype='mean count',comm_type=comm_type)
res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
#plot by ligand category
#overall network plot
NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
#top 20 ligand-receptor pairs
LRPlot(res_cat[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:20],link.arr.width=res_cat$cell_to_mean_exprs[1:20])
title(comm_type)
res<-rbind(res,res_cat)
}

saveRDS(res, "res1.rds")
saveRDS(res_cat, "res_cat.rds")

res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:20,]
NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
LRPlot(res[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res$cell_from_mean_exprs[1:20],link.arr.width=res$cell_to_mean_exprs[1:20])
counts_responders0<-counts_responders0 %>% mutate(compare_group=sample(2,nrow(counts_responders0),replace=TRUE))
deg_t<-DEG(counts_responders0 %>% filter(cell_type=='Monocytes/macrophages'),method='Wilcox',contrast=c(2,1))
deg_nk<-DEG(counts_responders0 %>% filter(cell_type=='Monocytes/macrophages'),method='Wilcox',contrast=c(2,1))
par(mfrow=c(1,2))
res<-NULL

saveRDS(res, "res2.rds")

#non-responders top ligands in different cell types
nonresp_cells <- vector(mode = 'list', length = length(unique(res_cat_nonresp$cell_from)))
for (k in 1:length(unique(res_cat_nonresp$cell_from))){
  nonresp_onecell <- subset(res_cat_nonresp, cell_from == unique(res_cat_nonresp$cell_from)[k])
  nonresp_toplig <- data.frame(Ligand = character(), Expr = double())
  for (i in 1:length(unique(nonresp_onecell$ligand))) {
    operation <- subset(nonresp_onecell, ligand == unique(nonresp_onecell$ligand)[i])
    nonresp_toplig[i,1] <- unique(nonresp_onecell$ligand)[i]
    nonresp_toplig[i,2] <-  mean(operation$cell_from_mean_exprs)
    nonresp_cells[[k]] <- nonresp_toplig
  }
}
names(nonresp_cells) <- unique(res_cat_nonresp$cell_from)

#responders top ligands in different cell types
resp_cells <- vector(mode = 'list', length = length(unique(res_cat_resp$cell_from)))
for (k in 1:length(unique(res_cat_resp$cell_from))){
  resp_onecell <- subset(res_cat_resp, cell_from == unique(res_cat_resp$cell_from)[k])
  resp_toplig <- data.frame(Ligand = character(), Expr = double())
  for (i in 1:length(unique(resp_onecell$ligand))) {
    operation <- subset(resp_onecell, ligand == unique(resp_onecell$ligand)[i])
    resp_toplig[i,1] <- unique(resp_onecell$ligand)[i]
    resp_toplig[i,2] <-  mean(operation$cell_from_mean_exprs)
    resp_toplig <- resp_toplig %>% arrange(Ligand)
    resp_cells[[k]] <- resp_toplig
  }
}
names(resp_cells) <- unique(res_cat_resp$cell_from)

#non-responders CD274, PDCD1LG2 and CD80/86
nonresp_cells_spec <- vector(mode = 'list', length = length(unique(res_cat_nonresp$cell_from)))
for (k in 1:length(res_cat_nonresp$cell_from)){
  nonresp_onecell_spec <- subset(res_cat_nonresp, cell_from == unique(res_cat_nonresp$cell_from)[k])
  nonresp_toplig_spec <- data.frame(Ligand = character(), Expr = double())
  for (i in 1:length(c('CD86', "CD274", "CD80", "PDCD1LG2"))) {
    operation <- subset(nonresp_onecell_spec, ligand == c('CD86', "CD274", "CD80", "PDCD1LG2")[i])
    nonresp_toplig_spec[i,1] <- c('CD86', "CD274", "CD80", "PDCD1LG2")[i]
    nonresp_toplig_spec[i,2] <-  mean(operation$cell_from_mean_exprs)
    nonresp_cells_spec[[k]] <- nonresp_toplig_spec
  }
}
names(nonresp_cells_spec) <- unique(res_cat_nonresp$cell_from)

#responders CD274, PDCD1LG2 and CD80/86
resp_cells_spec <- vector(mode = 'list', length = length(unique(res_cat_resp$cell_from)))
for (k in 1:length(res_cat_resp$cell_from)){
  resp_onecell_spec <- subset(res_cat_resp, cell_from == unique(res_cat_resp$cell_from)[k])
  resp_toplig_spec <- data.frame(Ligand = character(), Expr = double())
  for (i in 1:length(c('CD86', "CD274", "CD80", "PDCD1LG2"))) {
    operation <- subset(resp_onecell_spec, ligand == c('CD86', "CD274", "CD80", "PDCD1LG2")[i])
    resp_toplig_spec[i,1] <- c('CD86', "CD274", "CD80", "PDCD1LG2")[i]
    resp_toplig_spec[i,2] <-  mean(operation$cell_from_mean_exprs)
    resp_cells_spec[[k]] <- resp_toplig_spec
  }
}
names(resp_cells_spec) <- unique(res_cat_resp$cell_from)

#Example of the heatmap
resp_heat <- setNames(data.frame(matrix(ncol = 16, nrow = 11)), resp_all_lig)
rownames(resp_heat) <- names(resp_cells)
for (r in 1:length(rownames(resp_heat))){
  print(r)
  for (coln in colnames(resp_heat)){
    orig <- subset(resp_cells[[r]], Ligand == coln)[1,2]
    resp_heat[[coln]][r] <- as.numeric(orig)
  }
}
resp_heat <- as.matrix(resp_heat)
heatmap(resp_heat, main ='Responders', dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')

#Final heatmap
res.nonresponders$test <- paste(res.nonresponders$cell_from,res.nonresponders$cell_to)
res.responders$test <- paste(res.responders$cell_from,res.responders$cell_to)
res.responders$mult <- res.responders$cell_from_mean_exprs*res.responders$cell_to_mean_exprs
res.nonresponders$mult <- res.nonresponders$cell_from_mean_exprs*res.nonresponders$cell_to_mean_exprs
test4 <- merge(res.responders, res.nonresponders, by = c("ligand","receptor","test"), all = T)
test4$multD <- test4$mult.x - test4$mult.y
test4 <- test4[order(-abs(test4$multD)),]
write.table(head(test4,500), file = "/home/v/Projects/POI/Sergei_Nikita/test091220.txt", sep = "\t", quote = F, row.names = F)

test4$absmultD <- abs(test4$multD)
test5 <- head(test4,100)
temp <- as.data.frame(unique(paste(test5$ligand,test5$receptor)))
test4check <- subset(test4 , comm_type.x == "checkpoint")
# checkpoints
SUMMARY_check <- head(test4check[,c(1,2,5,7,4,6,10,12,16)],30)
rownames(SUMMARY_check) <- paste(SUMMARY_check$ligand, SUMMARY_check$receptor, SUMMARY_check$cell_from.x, SUMMARY_check$cell_to.x)
SUMMARY_check[,1:4] <- NULL 
temp2 <- log2(SUMMARY_check+1)
library(pheatmap)
pheatmap(SUMMARY_check[,1:4],
         #scale = "row",
         fontsize_row= 6, 
         cellwidth = 10, 
         cluster_row = T, 
         cluster_cols = T,
         show_colnames = T,
         show_rownames = T,
         border_color = FALSE
)
temp2$multD <- NULL
colnames(temp2) <- c("L_Responder", "R_Responder" , "L_Non-responder", "R_Non-responder")
rownames(temp2)

pheatmap(temp2[,c(1,3,2,4)],
         #scale = "row",
         fontsize_row= 10, 
         cellwidth = 10, 
         cluster_row = F, 
         cluster_cols = F,
         show_colnames = T,
         show_rownames = T,
         border_color = FALSE
)
