load("DEGAll.rda")
sigDEG=DEGAll[abs(DEGAll$logFC)>1 & DEGAll$FDR<0.05,]
library(DOSE)
library(clusterProfiler)
eg = bitr(rownames(sigDEG), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#进行KEGG富集分析
kk <- enrichKEGG (gene = eg[[2]],    #显著差异表达基因的ENTREZID
                 organism     = 'hsa',   #富集分析的物种，
                 pAdjustMethod = "BH",   #FDR校正p值
                 pvalueCutoff = 0.01,    #p值阈值<0.05
                 qvalueCutoff = 0.01,    #q值阈值<0.2
                 minGSSize = 10 ,         #富集的GO条目至少包含10个基因
)
#将富集分析结果中的ENTREZID再转换成基因名字
kk=setReadable(kk,OrgDb=org.Hs.eg.db,keyType  ="ENTREZID")
write.csv(as.data.frame(kk),"KEGG-enrich.csv",row.names =F)
##气泡图和柱状图显示KEGG富集结果
pdf(file="KEGG_bar.pdf",width=10)
barplot(kk, showCategory=15,title="KEGG Enrichment")
dev.off()
pdf(file="KEGG_dot.pdf",width=10)
dotplot(kk,showCategory=15,title="KEGG Enrichment")
dev.off()
##在KEGG通路上显示差异表达的基因，显示高表达还是低表达
library(pathview)
#获取fold change值，用于标注基因上调还是下调
deg <- sigDEG$logFC
names(deg) <- rownames(sigDEG)
#hsa04110   cell cycle
hsa04110 <- pathview(gene.data   = deg,  #包含基因ID和fold change
                     pathway.id  = "hsa04110", #通路ID
                     species     = "hsa",   #物种缩写
                     gene.idtype = 'ENSEMBL',  #基因ID类型
                     kegg.native=TRUE,    #绘制png图
                     limit       = list(gene=max(abs(deg)), cpd=1))   #设置图注的范围

