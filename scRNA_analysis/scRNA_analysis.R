#' ---
#' title: "snRNA-seq Meninges Ananlysis"
#' author: "Nathan Boles"
#' date: "August 9th, 2023"
#' ---

###Libraries utilized
library(foreach)
library(GO.db)
library(biomaRt)
library(Seurat)
library(disgenet2r)
library(screp)
library(dplyr)
library(stringr)
library(GOfuncR)

library(hypeR)
library(polycor)
library(rcompanion)
library(ggplot2)
library(babelgene)
library(SummarizedExperiment)
library(DESeq2)
library(RCy3)

#### data was processed to a seurat object at end of mapping with STAR mapper. Using those objects for metadata
######set up gene and go annotation objects
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
x<-rownames(young.combined@assays$RNA@counts)
x<-unique(c(rownames(young.combined@assays$RNA@counts),rownames(mm126191@assays$RNA@counts)))
x<-strsplit(x,"-")
x<-unlist(x)
x<-x[seq(1,length(x),by=2)]
x<- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=x,mart= mart)
y<-x[,2]
names(y)<-x[,1]
MGI<-y
MGI.rev<-names(MGI)
names(MGI.rev)<-MGI

x <- Term(GOTERM)
goterms<-names(x)
names(goterms)<-x


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes<-str_to_title(s.genes)
g2m.genes<-str_to_title(g2m.genes)

#####add descriptors to metadata from 
mm126191@meta.data$age<-"old"
mm126191@meta.data$position<-mm126191@meta.data$sample
#######since it is single nuclei data we obtained counts that includes introns 
z<-read.csv("126191_out_incl_introns_genematrix.csv",as.is=T,row.names=1)
x<-read.csv("126191_out_genematrix.csv",as.is=T,row.names=1)
z1<-read.csv("126254_out_incl_introns_genematrix.csv",as.is=T,row.names=1)
x1<-read.csv("126254_out_genematrix.csv",as.is=T,row.names=1)
z2<-read.csv("126275_out_incl_introns_genematrix.csv",as.is=T,row.names=1)
x2<-read.csv("126275_out_genematrix.csv",as.is=T,row.names=1)

## calculate out the ratio of  exon+intron to intron counts and remove infinte/nan cases
m<-z-x
m1<-z1-x1
m2<-z2-x2

test<-apply(z,1,sum)/apply(m,1,sum)
test<-test[-(which(is.infinite(test)))]
test<-test[-(which(is.nan(test)))]

test1<-apply(z1,1,sum)/apply(m1,1,sum)
test1<-test1[-(which(is.infinite(test1)))]
test1<-test1[-(which(is.nan(test1)))]

test2<-apply(z2,1,sum)/apply(m2,1,sum)
test2<-test2[-(which(is.infinite(test2)))]
test2<-test2[-(which(is.nan(test2)))]

###remove those genes with a ratio greater than 100
### We reason that mRNA transcibed by the cell should have a fairly high level of intronic sequences vs non-nuclear mRNA. 
### By taking this approach we can weed out the greater portion of the ambient RNA

moo<-names(test[which(test>100)])
moo1<-names(test1[which(test1>100)])
moo2<-names(test2[which(test2>100)])

arna<-unique(c(moo,moo1,moo2))
genes<-setdiff(unique(c(names(test),names(test1),names(test2))),arna)

# create seurat objects for each age. Get metadata from object generated at end of mapping. 
##old(18months)
y<-mm126191@meta.data
m<-rownames(mm126191@meta.data)
m<-unlist(str_split(m,"_"))
m<-m[seq(2,length(m),2)]
rownames(y)<-m
y<-y[,10:11]
z<-z[,rownames(y)]
z<-z[genes,]
m<-rownames(z)
x<-strsplit(m,"_")
x<-unlist(x)
x<-x[seq(2,length(x),by=2)]
m<-duplicated(x)
z<-z[which(!m),]
rownames(z)<-x[which(!m)]

##young( 3months)
y0<-young.combined@meta.data
m<-rownames(young.combined@meta.data)
m<-unlist(str_split(m,"_"))
m<-m[seq(2,length(m),2)]
rownames(y0)<-m
z3<-cbind(z1,z2)
z3<-z3[,rownames(y0)]
z3<-z3[genes,]
m<-rownames(z3)
x<-strsplit(m,"_")
x<-unlist(x)
x<-x[seq(2,length(x),by=2)]
m<-duplicated(x)
z3<-z3[which(!m),]
rownames(z3)<-x[which(!m)]
y0<-y0[,9:10]

y0<-CreateSeuratObject(z3,meta.data=y0,min.cells=5,min.features=1000)
ol0<-CreateSeuratObject(z,meta.data=y,min.cells=5,min.features=1000)

## add cell cyle info to objects
y0<-CellCycleScoring(y0, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
y0<-PercentageFeatureSet(y0, "^mt-", col.name = "percent_mito")
ol0<-CellCycleScoring(ol0, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ol0<-PercentageFeatureSet(ol0, "^mt-", col.name = "percent_mito")

##young= y0, old=ol0
###analysis and cell annotation
##removing low count cells
y1<-subset(y0,subset=nCount_RNA>30000)
o1<-subset(ol0,subset=nCount_RNA>30000)
FeatureScatter(y1,feature1="nCount_RNA",feature2="nFeature_RNA",group.by = "position")
FeatureScatter(o1,feature1="nCount_RNA",feature2="nFeature_RNA",group.by = "position")
#####young analysis 
y1<-SCTransform(y1,vars.to.regress=c("S.Score", "G2M.Score","sample"))
y1<-RunPCA(y1)
ElbowPlot(y1)
y1<-RunUMAP(y1,dims=1:10)
y1<-FindNeighbors(y1,dims=1:10)
y1<-FindClusters(y1,resolution = 0.5)
DimPlot(y1,label=T)
DimPlot(y1,group.by="position",label=T)
DimPlot(y1,split.by="position",label=T)



# old

o1<-SCTransform(o1,vars.to.regress=c("S.Score", "G2M.Score"))
o1<-RunPCA(o1)
ElbowPlot(o1)
o1<-RunUMAP(o1,dims=1:10)
o1<-FindNeighbors(o1,dims=1:10)
o1<-FindClusters(o1,resolution = 0.8)
DimPlot(o1,label=T)
DimPlot(o1,group.by="position")
DimPlot(o1,split.by="position",label=T)

o1@misc$markers<-FindAllMarkers(o1,only.pos=T)
y1@misc$markers<-FindAllMarkers(y1,only.pos=T)

o1@misc$markers<-o1@misc$markers[o1@misc$markers$p_val_adj<0.1,]
y1@misc$markers<-y1@misc$markers[y1@misc$markers$p_val_adj<0.1,]
########integrate ages
z<-list(o1,y1)
z1<-SelectIntegrationFeatures(object.list=z,nfeatures=3000)
z<-PrepSCTIntegration(object.list=z,anchor.features=z1)
za<-FindIntegrationAnchors(object.list=z,normalization.method = "SCT",
                           anchor.features = z1,reduction="cca",k.anchor=30)
obj<-IntegrateData(anchorset=za,normalization.method = "SCT")

obj<-RunPCA(obj)
ElbowPlot(obj)
obj<-RunUMAP(obj,dims=1:10)
obj<-FindNeighbors(obj,dims=1:10)
obj<-FindClusters(obj,resolution = 0.5)
DimPlot(obj,label=T)
DimPlot(obj,group.by="position")
DimPlot(obj,split.by="age",label=T)
DimPlot(obj,split.by="position",label=T)

DefaultAssay(obj)<-"SCT"
obj<-PrepSCTFindMarkers(obj)

##Using the following markers for annotation
##Negr1, Dcn, Igfbp2  higher in Fibroblasts over endothelial
## Pdzd2 higher in endo over fibro
## Myh11 higher in SMC
##endo markers PMID: 33343397
##pericytes ebf1 PMID:34272603; Rpbj PMID: 31249304
###Fibroblasts
## PMID: 25013174; PMID: 32403949

x<-DotPlot(o1,features=c("Negr1","Igfbp2","Ppp1r1a","Serpine2",
                          "Runx1t1","Epha5","Cxadr","Tjp1","Nfix","Meis2",
                          "Aldh1a3",
                          "Col1a2","Dcn","Pdgfra",
                          "Pecam1","Ets1","Erg","Flt1",
                          "Angpt1","Anpep","Rbpj","Ebf1","Runx1",
                          "Abcc9",
                          "Nrp1",
                          "Efnb2","Fbln5","Vegfc"))

x+coord_flip()

x<-DotPlot(y1,features=c("Negr1","Igfbp2","Ppp1r1a","Serpine2",
                          "Runx1t1","Epha5","Cxadr","Tjp1","Nfix","Meis2",
                          "Aldh1a3",
                          "Col1a2","Dcn","Pdgfra",
                          "Pecam1","Ets1","Erg","Flt1",
                          "Angpt1","Anpep","Rbpj","Ebf1","Runx1",
                          "Abcc9",
                          "Nrp1",
                          "Efnb2","Fbln5","Vegfc"))

x+coord_flip()

x<-DotPlot(obj,features=c("Negr1","Igfbp2","Ppp1r1a","Serpine2",
                          "Runx1t1","Epha5","Cxadr","Tjp1","Nfix","Meis2",
                          "Aldh1a3",
                          "Col1a2","Dcn","Pdgfra",
                          "Pecam1","Ets1","Erg","Flt1",
                          "Angpt1","Anpep","Rbpj","Ebf1","Runx1",
                          "Abcc9",
                          "Nrp1",
                          "Efnb2","Fbln5","Vegfc"),assay="SCT")

x+coord_flip()



###young has  two classes of cells types fall out with a strong division by position
###Old is more mixed, but able to annotate cells in combined object

z<-as.vector(y1$seurat_clusters)
z<-replace(z,which(z=="0"),"Fibro1")
z<-replace(z,which(z=="1"),"Fibro2")
z<-replace(z,which(z=="2"),"Endo1")
z<-replace(z,which(z=="3"),"Endo2")
z<-replace(z,which(z=="4"),"Endo3")
z<-replace(z,which(z=="5"),"Fibro3")
z<-replace(z,which(z=="6"),"Fibro4")
y1$celltype<-z


z<-as.vector(o1$seurat_clusters)
z<-replace(z,which(z=="0"),"Fibro1")
z<-replace(z,which(z=="1"),"Endo1")
z<-replace(z,which(z=="2"),"Endo2")
z<-replace(z,which(z=="3"),"Fibro2")
z<-replace(z,which(z=="4"),"Fibro3")
o1$celltype<-z

z<-as.vector(obj$seurat_clusters)
z<-replace(z,which(z=="0"),"Fibro1")
z<-replace(z,which(z=="1"),"Fibro2")
z<-replace(z,which(z=="2"),"Endo1")
z<-replace(z,which(z=="3"),"Endo2")
z<-replace(z,which(z=="4"),"Endo3")
z<-replace(z,which(z=="5"),"Fibro3")
z<-replace(z,which(z=="6"),"Fibro4")
z<-replace(z,which(z=="7"),"Fibro5")
obj$celltype<-z

###cell type based markers
obj@misc$celltype<-FindAllMarkers(obj,only.pos=T)
obj@misc$celltype<-obj@misc$celltype[obj@misc$celltype$p_val_adj<0.1,]

###Figure panels for paper (Figure 1)
x<-DotPlot(obj,features=c("Ppp1r1a","Serpine2",
                          "Tjp1","Nfix","Meis2",
                          "Col1a2","Dcn","Ebf1","Pdzd2",
                         "Runx1","Fbln5","Vegfc","Myh11",
                          "Pecam1","Ets1","Erg"),
           group.by="celltype")

x+coord_flip()
y1$celltype<-as.factor(y1$celltype)
o1$celltype<-as.factor(o1$celltype)
obj$celltype<-as.factor(obj$celltype)
y1<-SetIdent(y1,value="celltype")
o1<-SetIdent(o1,value="celltype")
obj<-SetIdent(obj,value="celltype")
pdf("Figure_1_plots.pdf")
DimPlot(y1,label=T,pt.size=1.5)
DimPlot(o1,label=T,pt.size=1.5)
DimPlot(obj,label=T,pt.size=1.5)
DimPlot(y1,label=T,pt.size=1.5,split.by="position")
DimPlot(o1,label=T,pt.size=1.5,split.by="position")
DimPlot(obj,label=T,pt.size=1.5,split.by="position")
DimPlot(y1,label=T,pt.size=1.5,group.by="position")
DimPlot(o1,label=T,pt.size=1.5,group.by="position")
DimPlot(obj,label=T,pt.size=1.5,group.by="position")
DimPlot(obj,label=T,pt.size=1.5,split.by="age")
DimPlot(obj,label=T,pt.size=1.5,group.by="age")
dev.off()

###for Figure 1J

x<-DotPlot(obj,features=c("Ppp1r1a","Serpine2",
                          "Tjp1","Nfix","Meis2",
                          "Col1a2","Dcn",
                          "Ebf1","Runx1","Fbln5",
                          "Vegfc",
                          "Pecam1","Ets1","Erg","Myh11",
                          "Angpt1"),assay="SCT")

x+coord_flip()
pdf("Celltype_dotplot.pdf")
x
x+coord_flip()
dev.off()
###get numbers for cell type pie charts for Figure 1

x<-obj@meta.data

y<-rbind(table(x[x$group==1,]$celltype),
         table(x[x$group==2,]$celltype),
         table(x[x$group==3,]$celltype),
         table(x[x$group==4,]$celltype))
rownames(y)<-c("YA","OA","YP","OP")
write.csv(y,"CellType_age_position.csv")

###set up new metadata making groups of age and position
###1=YA,2=OA,3=YP,4=OP
### find markers for each group
x<-obj@meta.data$age
x<-replace(x,which(x=="young"),0)
x<-replace(x,which(x=="old"),1)
x<-as.numeric(x)

z<-obj@meta.data$position
z<-replace(z,which(z=="Anterior"),1)
z<-replace(z,which(z=="Posterior"),3)
z<-as.numeric(z)

z<-x+z
obj@meta.data$group<-z
obj<-SetIdent(obj,value="group")
obj@misc$group<-FindAllMarkers(obj,only.pos=T)
obj@misc$group<-obj@misc$group[obj@misc$group$p_val_adj<0.1,]

##disease enrichment using celltype and group markers
z<-intersect(obj@misc$group[obj@misc$group$p_val_adj<0.01,]$gene,
             obj@misc$celltype[obj@misc$celltype$p_val_adj<0.01,]$gene)
z<-orthologs(z,species="mouse",human=F)$human_symbol



disgenet_api_key <- get_disgenet_api_key(email="youremail",password="yourpassword")
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

res<-disease_enrichment(z)
y<-res@qresult
y<-y[y$FDR<0.01,]
x<-plot(res, class = "Enrichment", count =3,  cutoff= 0.01, nchars=50)

ze<-obj@misc$celltype[obj@misc$celltype$p_val_adj<0.01,]
ze<-intersect(z,orthologs(ze[grep("Endo",ze$cluster),]$gene,species="mouse",human=F)$human_symbol)
res_e<-disease_enrichment(ze)
x1<-plot(res_e, class = "Enrichment", count =3,  cutoff= 0.01, nchars=50)
zf<-obj@misc$celltype[obj@misc$celltype$p_val_adj<0.01,]
zf<-intersect(z,orthologs(zf[grep("Fibro",zf$cluster),]$gene,species="mouse",human=F)$human_symbol)
res_f<-disease_enrichment(zf)
x2<-plot(res_f, class = "Enrichment", count =3,  cutoff= 0.01, nchars=50)


###make supplementary tables
b<-res@qresult
b2<-res_e@qresult
b3<-res_f@qresult
b<-b[b$FDR<0.01,]
b2<-b2[b2$FDR<0.01,]
b3<-b3[b3$FDR<0.01,]
b$group<-"all"
b2$group<-"endo"
b3$group<-"fibro"
write.csv(rbind(b,b2,b3),"Disease_Enrichment.csv")

write.csv(obj@misc$celltype,"Celltype_markers.csv")
write.csv(obj@misc$group,"Group_markers.csv")

###Make panel for Figure 2
x4<-rbind(x$data,x1$data,x2$data)
x4$group<-c(rep("all",dim(x$data)[1]),rep("endo",dim(x1$data)[1]),rep("fibro",dim(x2$data)[1]))
x5<-x4[x4$gg>0.01,]

p<-ggplot(x5,aes(x=gg,y=Description,color=FDR,size=Count))
p<-p+geom_point()
p<-p+theme_bw()
p<-p+facet_grid(.~group)
p<-p+ scale_color_gradient(low="red",high="blue")

pdf("disease_panel.pdf")
p
dev.off()
####GO enrichment of groups and cell types
## by cell type 
GOBP <- msigdb_gsets(species="Mus musculus", category="C5", subcategory="BP")
x<-obj@misc$celltype
x<-x[x$p_val_adj<0.05,]
x<-x[,6:7]
y<-unique(x$cluster)
x<-foreach(i=1:length(unique(x$cluster))) %do% {
  x[x$cluster==unique(x$cluster)[i],2]
}
names(x)<-y
eres<-enrich_test(clust_list=x,species_x="Mus musculus",genome_genes=42937)
go_vis_sc<-GO_visualization(eres$Enriched_df,clust_list=x,GOcats=GOBP,goterms=goterms,numcats=10,org_db="org.Mm.eg.db")
go_vis_sc$plot

##plot panel for figure 2
y<-c("neurogenesis","tissue development","aging","cell morphogenesis","regulation of cell differentiation",
     "programmed cell death","response to abiotic stimulus","regulation of transport","cell projection organization","cytokine production")
test<-GO_viz_choose(go_vis_sc,clust_list=x,chosen_cats=y,goterms=goterms,species_x="Mus musculus")
plot(test$plot)

##by group
x<-obj@misc$group
x<-x[x$p_val_adj<0.05,]
x<-x[,6:7]
y<-unique(x$cluster)
x<-foreach(i=1:length(unique(x$cluster))) %do% {
  x[x$cluster==unique(x$cluster)[i],2]
}
names(x)<-y

eres2<-enrich_test(clust_list=x,species_x="Mus musculus",genome_genes=42937)
go_vis_sc2<-GO_visualization(eres2$Enriched_df,clust_list=x,GOcats=GOBP,goterms=goterms,numcats=10,org_db="org.Mm.eg.db")
go_vis_sc2$plot

##plot panel for figure 2
y<-c("neurogenesis","tissue development","aging","cell morphogenesis","regulation of cell differentiation",
     "programmed cell death","response to abiotic stimulus","regulation of transport","cell projection organization","cytokine production")
test1<-GO_viz_choose(go_vis_sc2,clust_list=x,chosen_cats=y,goterms=goterms,species_x="Mus musculus")
plot(test1$plot)

pdf("GO_enrich_Panel.pdf")
plot(test$plot)
plot(test1$plot)
dev.off()

###write tables for supplement
write.csv(go_vis_sc$GO_sem,"celltype_GO_table.csv")
write.csv(go_vis_sc2$GO_sem,"group_GO_table.csv")

###obtain Schizophrenia data from BrainSeq phase II
###https://eqtl.brainseq.org/phase2/
###downloaded RangedSummarizedExperiment obj : Gene along includes counts and metadata
###selcted data from dorsolateral prefrontal cortex (DLPFC)

schizo_cts<-assays(rse_gene)$counts
schizo_cts<-schizo_cts[,rownames(schizo_meta[schizo_meta$Region=="DLPFC",])]

schizo_meta<-colData(rse_gene)
schizo_meta<-schizo_meta[,30:54]
schizo_meta$Dx<-as.factor(schizo_meta$Dx)
schizo_meta<-schizo_meta[schizo_meta$Region=="DLPFC",]

###due to number of samples for control and schizo conditions chose to divide sample by sex 
###and test by wilcox rank test due to results from PMID:35292087
###first normalized counts via DESeq2 median of ratios method, got rid of low count genes, and 
###categorized samples by sex and condition into a design vector
schizo<-DESeqDataSetFromMatrix(schizo_cts,
                               colData=schizo_meta,
                               design=~Sex + Dx)
schizo<-schizo[rowSums(counts(schizo)>20)>=45,]
schizo <- estimateSizeFactors(schizo)
norm_schizo <- counts(schizo, normalized=TRUE)
design<-data.frame(sex=schizo_meta$Sex,Dx=as.character(schizo_meta$Dx))
design[design=="M"]<-0
design[design=="F"]<-1
design[design=="Control"]<-1
design[design=="Schizo"]<-3
design<-as.numeric(design$sex)+as.numeric(design$Dx)
###1=M_ctrl,2=F_ctrl,3=M_schizo,4=F_schizo
###use a wilcox-rank test to find differential genes for condition comparison within each sex
testF<-cbind(norm_schizo[,design==2],norm_schizo[,design==4])
designF<-c(rep(0,length(which(design==2))),rep(1,length(which(design==4))))

pvalsF <- sapply(1:nrow(testF),function(i){
  data<-cbind.data.frame(gene=as.numeric(testF[i,]),designF)
  p=wilcox.test(gene~designF, data)$p.value
  return(p)
})

testM<-cbind(norm_schizo[,design==1],norm_schizo[,design==3])
designM<-c(rep(0,length(which(design==1))),rep(1,length(which(design==3))))

pvalsM <- sapply(1:nrow(testM),function(i){
  data<-cbind.data.frame(gene=as.numeric(testM[i,]),designM)
  p=wilcox.test(gene~designM, data)$p.value
  return(p)
})
### clean up rownames of each counts matrix
x<-(unlist(strsplit(rownames(testF),"\\."))[seq(1,2*length(rownames(testF)),2)])
rownames(testF)<-x

x<-(unlist(strsplit(rownames(testM),"\\."))[seq(1,2*length(rownames(testF)),2)])
rownames(testM)<-x

### make data frames for each sex and select significant genes
### with a absolute logFC 0.5 threshold

resF<-data.frame(pvalues=pvalsF,
                 adj_p=p.adjust(pvalsF,"fdr"))
rownames(resF)<-rownames(testF)
resM<-data.frame(pvalues=pvalsM,
                 adj_p=p.adjust(pvalsM,"fdr"))
rownames(resM)<-rownames(testM)

resF$logFC<-log2(apply(testF[,designF==1],1,mean)/apply(testF[,designF==0],1,mean))
resF$baseMean<-apply(testF,1,mean)
resF.sig<-resF[resF$adj_p<0.05,]
resF.sig<-resF.sig[which(abs(resF.sig$logFC)>0.5),]

resM$logFC<-log2(apply(testM[,designM==1],1,mean)/apply(testM[,designM==0],1,mean))
resM$baseMean<-apply(testM,1,mean)
resM.sig<-resM[resM$adj_p<0.05,]
resM.sig<-resM.sig[which(abs(resM.sig$logFC)>0.5),]

length(intersect(rownames(resF.sig),rownames(resM.sig)))
##Females have 1532 significant genes and males have 401 with 342 shared
## add column with mouse orthologs
x<- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
          values=rownames(resM.sig),mart= mart)
x<-orthologs(rownames(resM.sig),species="mouse")
y<-x$symbol
names(y)<-x$human_ensembl
resM.sig$mus_gene<-NA
resM.sig[names(y),]$mus_gene<-y

x<- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
          values=rownames(resF.sig),mart= mart)
x<-orthologs(rownames(resF.sig),species="mouse")
y<-x$symbol
names(y)<-x$human_ensembl
resF.sig$mus_gene<-NA
resF.sig[names(y),]$mus_gene<-y

### use Schizo data to look for ligands from meninges to receptors in brain
##function for identifying ligand receptor pairs

partner<-function(y,targets,LR_db,ligand=TRUE) {
  if(isTRUE(ligand)){
    z<-foreach(i=1:length(y),.combine='rbind') %do% {
      z1<-LR_db[which(LR_db$ligand_gene_symbol==y[i]),]
      z1[,2:3]
    }
    z2<-intersect(targets,z[,2])
    z2<-data.frame(z2)
    colnames(z2)<-colnames(z)[2]
    z2<-dplyr::right_join(z,z2)
    z2<-z2[!(is.na(z2$ligand_gene_symbol)),]
    return(z2)} else{
      z<-foreach(i=1:length(y),.combine='rbind') %do% {
        z1<-LR_db[which(LR_db$receptor_gene_symbol==y[i]),]
        z1[,2:3]
      }
      z2<-intersect(targets,z[,1])
      z2<-data.frame(z2)
      colnames(z2)<-colnames(z)[1]
      z2<-dplyr::right_join(z,z2)
      z2<-z2[!(is.na(z2$receptor_gene_symbol)),]
      return(z2)
    }
  
}

###downloaded Ligand-Receptod data from http://tcm.zju.edu.cn/celltalkdb/download.php
humLR<-read.delim("human_lr_pair.txt")
musLR<-read.delim("mouse_lr_pair.txt")

x<-unique(obj@misc$group$gene)
y<-intersect(x,musLR$ligand_gene_symbol)

###against schizo DEG
z<-unique(resF.sig$mus_gene)
z<-z[-(which(is.na(z)))]
LRF_sig<-partner(y,z,LR_db=musLR)

z<-unique(resM.sig$mus_gene)
z<-z[-(which(is.na(z)))]
LRM_sig<-partner(y,z,LR_db=musLR)

###against all reasonably expressed receptors in brain data


z<-orthologs(rownames(resF[resF$baseMean>15,]),species="mouse")
z<-unique(z$symbol)
LRF<-partner(y,z,LR_db=musLR)

z<-orthologs(rownames(resM[resM$baseMean>15,]),species="mouse")
z<-unique(z$symbol)
LRM<-partner(y,z,LR_db=musLR)

#### function to help identify cell type, region, and age associated with ligand 



anno_group<-function(obj_reg,obj_ct,LR_group){
  m<-as.character(obj_reg$cluster)
  names(m)<-obj_reg$gene
  reg_age<-m[intersect(obj_reg$gene,LR_group$ligand_gene_symbol)]
  reg_age<-replace(reg_age,which(reg_age==1),"Young_anterior")
  reg_age<-replace(reg_age,which(reg_age==2),"old_anterior")
  reg_age<-replace(reg_age,which(reg_age==3),"Young_posterior")
  reg_age<-replace(reg_age,which(reg_age==4),"old_posterior")
  
  m<-as.character(obj_ct$cluster)
  names(m)<-obj_ct$gene
  celltype<-m[intersect(obj_ct$gene,LR_group$ligand_gene_symbol)]
  
  z<-setdiff(names(reg_age),names(celltype))
  z1<-setdiff(names(celltype),names(reg_age))
  if(length(z)==0) {} else {
    b<-rep(NA,length(z))
    names(b)<-z
  }
  
  if(length(z1)==0) {} else {
    b1<-rep(NA,length(z1))
    names(b1)<-z1
  }
  z3<-intersect(names(reg_age),names(celltype))
  df1<-data.frame(reg_age[z3],celltype[z3])
  colnames(df1)<-c("age_reg","celltype")
  
  if(length(z)==0) {} else{
    z<-reg_age[z]
    df2<-data.frame(z,b)
    colnames(df2)<-c("age_reg","celltype")
    df1<-rbind(df1,df2)
  }
  
  if(length(z1)==0) {} else{
    z1<-celltype[z1]
    df2<-data.frame(z1,b1)
    colnames(df2)<-c("age_reg","celltype")
    df1<-rbind(df1,df2)
  }
  df1$gene<-rownames(df1)
  return(df1)
 }

### looking at each condition


LRFsig_anno<-anno_group(obj@misc$group,obj@misc$celltype,LRF_sig)
LRF_anno<-anno_group(obj@misc$group,obj@misc$celltype,LRF)
LRMsig_anno<-anno_group(obj@misc$group,obj@misc$celltype,LRM_sig)
LRM_anno<-anno_group(obj@misc$group,obj@misc$celltype,LRM)

df1<-semi_join(LRFsig_anno,LRF_anno,by="gene")
df2<-anti_join(LRFsig_anno,LRF_anno,by="gene")
df1<-rbind(df1,df2)
df1$Schizo_Sig<-TRUE
df2<-anti_join(LRF_anno,LRFsig_anno,by="gene")
df2$Schizo_Sig<-FALSE
df1<-rbind(df1,df2)
df1$Sex<-"Female"
LRF_anno<-df1

df1<-semi_join(LRMsig_anno,LRF_anno,by="gene")
df2<-anti_join(LRMsig_anno,LRF_anno,by="gene")
df1<-rbind(df1,df2)
df1$Schizo_Sig<-TRUE
df2<-anti_join(LRM_anno,LRMsig_anno,by="gene")
df2$Schizo_Sig<-FALSE
df1<-rbind(df1,df2)
df1$Sex<-"Male"
LRM_anno<-df1

write.csv(rbind(LRF_anno,LRM_anno),"LR_schizo_all.csv")

###make networks for cytoscape and paper figure
##female
LRF_sig$mus_gene<-LRF_sig$receptor_gene_symbol 
x<-left_join(LRF_sig,resF.sig,by="mus_gene")
x$mus_gene<-x$ligand_gene_symbol
LRFsig_anno$mus_gene<-LRFsig_anno$gene
y<-left_join(x,LRFsig_anno,by="mus_gene")
y<-y[,-3]
colnames(y)<-c(colnames(LRF),
               paste("Receptor",colnames(resF.sig)[1:4],sep="_"),
               paste("Ligand",colnames(LRFsig_anno)[1:3],sep="_"))
y<-y[,-9]
z<-c(rep("ligand",dim(y)[1]),rep("receptor",dim(y)[1]))
z<-data.frame(c(y[,1],y[,2]),
              z,
              c(y$Ligand_age_reg,rep(0,dim(y)[1])),
              c(y$Ligand_celltype,rep(0,dim(y)[1])),
              c(rep(0,dim(y)[1]),y$Receptor_logFC)
              )
colnames(z)<-c("Gene","Type","Age_region","Cell","Schizo_logFC")
z<-z[!duplicated(z), ]
test<-graph_from_data_frame(y[,1:2],vertices=z)
plot(test)
createNetworkFromIgraph(test,"SFemale")

##Male
LRM_sig$mus_gene<-LRM_sig$receptor_gene_symbol 
x<-left_join(LRM_sig,resM.sig,by="mus_gene")
x$mus_gene<-x$ligand_gene_symbol
LRMsig_anno$mus_gene<-LRMsig_anno$gene
y<-left_join(x,LRMsig_anno,by="mus_gene")
y<-y[,-3]
colnames(y)<-c(colnames(LRM),
               paste("Receptor",colnames(resF.sig)[1:4],sep="_"),
               paste("Ligand",colnames(LRFsig_anno)[1:3],sep="_"))
y<-y[,-9]
z<-c(rep("ligand",dim(y)[1]),rep("receptor",dim(y)[1]))
z<-data.frame(c(y[,1],y[,2]),
              z,
              c(y$Ligand_age_reg,rep(0,dim(y)[1])),
              c(y$Ligand_celltype,rep(0,dim(y)[1])),
              c(rep(0,dim(y)[1]),y$Receptor_logFC)
)
colnames(z)<-c("Gene","Type","Age_region","Cell","Schizo_logFC")
z<-z[!duplicated(z), ]
test<-graph_from_data_frame(y[,1:2],vertices=z)
plot(test)
createNetworkFromIgraph(test,"SMale")

###check if receptors for meninges secreted ligands are enriched
###female
x<-unique(LRF_sig$receptor_gene_symbol)
z<-intersect(rownames(resF[resF$baseMean>15,]),humLR$receptor_ensembl_gene_id)
z<-orthologs(z,species="mouse")$symbol

phyper(length(x),
       length(intersect(z,resF.sig$mus_gene)),
       length(rownames(resF[resF$baseMean>15,]))-length(z),
       length(unique(resF.sig$mus_gene))-1,lower.tail=F)
###1.584996e-27

###Male

x<-unique(LRM_sig$receptor_gene_symbol)
z<-intersect(rownames(resM[resM$baseMean>15,]),humLR$receptor_ensembl_gene_id)
z<-orthologs(z,species="mouse")$symbol

phyper(length(x),
       length(intersect(z,resM.sig$mus_gene)),
       length(rownames(resM[resM$baseMean>15,]))-length(z),
       length(unique(resM.sig$mus_gene))-1,lower.tail=F)
###1.153724e-25

####test to see if general mouse fibroblasts and endothelial cells have enriched synapse categories

x<-go_vis_sc$GO_sem
y<-x[grep("Fibro",x$cluster),]
length(unique(y$term))
length(unique(y[grep("synap",y$term),]$term))
length(unique(y[grep("neur",y$term),]$term))
y<-x[grep("Endo",x$cluster),]
length(unique(y$term))
length(unique(y[grep("synap",y$term),]$term))
length(unique(y[grep("neur",y$term),]$term))

##get markers from PanglaoDB
pangdb<-data.frame(read.delim("PanglaoDB_markers_27_Mar_2020.tsv",as.is=T))
pangdb<-pangdb[grep("Mm",pangdb$species),]

test<-pangdb[grep("Fibro",pangdb$cell.type),]
x<-test$official.gene.symbol
x<-orthologs(x,species="mouse")$symbol
y<-list(x)
fibro_test<-enrich_test(clust_list=y,species_x="Mus musculus",genome_genes=42937)
x<-fibro_test$Enriched_df[[1]]$GOBP
x[grep("SYNAP",x$label),]$label
x[grep("NEUR",x$label),]$label
###no synaptic categories were enriched out, 27 enriched categories related to neuron/neural
###out of a total of 994 enriched categories 

test<-pangdb[grep("Endo",pangdb$cell.type),]
x<-test$official.gene.symbol
x<-orthologs(x,species="mouse")$symbol
y<-list(x)
endo_test<-enrich_test(clust_list=y,species_x="Mus musculus",genome_genes=42937)
x<-endo_test$Enriched_df[[1]]$GOBP
x[grep("SYNAP",x$label),]$label
x[grep("NEUR",x$label),]$label
###endo has two synapse cats and 18 neuron cats out of 1095

#####compare to embryonic meninges PMID:32634398
library(VennDiagram)
embryo_men<-read.csv("PMID32634398_Table_S3.csv")
x<-embryo_men[embryo_men$p_val_adj<0.05,]
y<-unique(x$cluster)
x<-foreach(i=1:length(unique(x$cluster))) %do% {
  x[x$cluster==unique(x$cluster)[i],]$gene
}
names(x)<-y
emen_res<-enrich_test(clust_list=x,species_x="Mus musculus",genome_genes=42937)
go_vis_emen<-GO_visualization(emen_res$Enriched_df,clust_list=x,GOcats=GOBP,goterms=goterms,numcats=10,org_db="org.Mm.eg.db")
x<-go_vis_sc$GO_sem[grep("Fibro",go_vis_sc$GO_sem$cluster),]$parentTerm
ggVennDiagram(list("Adult"=x,"Embryonic"=go_vis_emen$GO_sem$parentTerm))
draw.pairwise.venn(area1=(103+89),area2=(54+89),cross.area = 89,fill = c("blue", "red"), alpha = rep(0.4, 2),
                   category=c("Adult","Embryo"))
dev.off()
ggVennDiagram(list("Adult"=go_vis_sc2$GO_sem$parentTerm,"Embryonic"=go_vis_emen$GO_sem$parentTerm))
draw.pairwise.venn(area1=(81+81),area2=(62+81),cross.area = 81,fill = c("blue", "red"), alpha = rep(0.4, 2),
                   category=c("Adult","Embryo"))
dev.off()
y<-obj@misc$celltype
y<-y[y$p_val_adj<0.05,]
y<-unique(y[grep("Fibro",y$cluster),]$gene)
ggVennDiagram(list("Adult"=y,"Embryonic"=unique(unlist(x))))
draw.pairwise.venn(area1=(2863+205),area2=(437+205),cross.area = 205,fill = c("blue", "red"), alpha = rep(0.4, 2),
                   category=c("Adult","Embryo"))
dev.off()
x<-go_vis_emen$GO_sem$term
length(unique(x))
x[grep("synap",x)]
##19
x[grep("neur",x)]
##90 out of 1073


sessionInfo()
R version 4.2.0 (2022-04-22 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RCy3_2.16.0                 babelgene_22.3              ggplot2_3.3.6              
 [4] rcompanion_2.4.15           polycor_0.8-1               GOfuncR_1.16.0             
 [7] vioplot_0.3.7               zoo_1.8-10                  sm_2.2-5.7.1               
[10] dplyr_1.0.9                 screp_1.0.0                 disgenet2r_0.99.2          
[13] sp_1.5-0                    SeuratObject_4.1.0          Seurat_4.1.1               
[16] GO.db_3.15.0                AnnotationDbi_1.58.0        foreach_1.5.2              
[19] biomaRt_2.52.0              DESeq2_1.36.0               SummarizedExperiment_1.26.1
[22] Biobase_2.56.0              MatrixGenerics_1.8.1        matrixStats_0.62.0         
[25] GenomicRanges_1.48.0        GenomeInfoDb_1.32.2         IRanges_2.30.0             
[28] S4Vectors_0.34.0            BiocGenerics_0.42.0         STRINGdb_2.8.4             
[31] hypeR_1.12.0                stringr_1.4.0              

loaded via a namespace (and not attached):
  [1] rgl_0.109.2            ica_1.0-2              svglite_2.1.0          class_7.3-20          
  [5] lmtest_0.9-40          crayon_1.5.1           spatstat.core_2.4-4    MASS_7.3-56           
  [9] nlme_3.1-157           backports_1.4.1        GOSemSim_2.22.0        rlang_1.1.1           
 [13] XVector_0.36.0         ROCR_1.0-11            readxl_1.4.0           irlba_2.3.5           
 [17] filelock_1.0.2         proto_1.0.0            BiocParallel_1.30.3    LSAfun_0.6.2          
 [21] bit64_4.0.5            glue_1.6.2             pheatmap_1.0.12        sctransform_0.3.3     
 [25] parallel_4.2.0         spatstat.sparse_2.1-1  base64url_1.4          spatstat.geom_2.4-0   
 [29] tidyselect_1.1.2       fitdistrplus_1.1-8     XML_3.99-0.10          tidyr_1.2.0           
 [33] chron_2.3-61           xtable_1.8-4           magrittr_2.0.3         evaluate_0.15         
 [37] cli_3.3.0              zlibbioc_1.42.0        rstudioapi_0.13        miniUI_0.1.1.1        
 [41] rpart_4.1.16           wordcloud_2.6          RJSONIO_1.3-1.8        shiny_1.7.1           
 [45] xfun_0.31              tm_0.7-8               cluster_2.1.3          caTools_1.18.2        
 [49] pbdZMQ_0.3-9           KEGGREST_1.36.2        tibble_3.1.7           expm_0.999-6          
 [53] ggrepel_0.9.1          mapplots_1.5.1         listenv_0.8.0          Biostrings_2.64.0     
 [57] png_0.1-7              future_1.26.1          withr_2.5.0            lsa_0.73.3            
 [61] bitops_1.0-7           slam_0.1-50            ggforce_0.3.3          plyr_1.8.7            
 [65] cellranger_1.1.0       e1071_1.7-11           pillar_1.7.0           gplots_3.1.3          
 [69] cachem_1.0.6           multcomp_1.4-19        fs_1.5.2               NLP_0.2-1             
 [73] hash_2.2.6.3           vctrs_0.4.1            ellipsis_0.3.2         generics_0.1.2        
 [77] gsubfn_0.7             nortest_1.0-4          tools_4.2.0            munsell_0.5.0         
 [81] tweenr_1.0.2           proxy_0.4-27           DelayedArray_0.22.0    fastmap_1.1.0         
 [85] compiler_4.2.0         abind_1.4-5            httpuv_1.6.5           DescTools_0.99.45     
 [89] plotly_4.10.0          rgeos_0.5-9            GenomeInfoDbData_1.2.8 gridExtra_2.3         
 [93] lattice_0.20-45        deldir_1.0-6           visNetwork_2.1.0       utf8_1.2.2            
 [97] later_1.3.0            BiocFileCache_2.4.0    jsonlite_1.8.0         rrvgo_1.8.0           
[101] scales_1.2.0           gld_2.6.4              graph_1.74.0           pbapply_1.5-0         
[105] genefilter_1.78.0      lazyeval_0.2.2         promises_1.2.0.1       goftest_1.2-3         
[109] spatstat.utils_3.0-1   reticulate_1.25        rmarkdown_2.14         openxlsx_4.2.5        
[113] sandwich_3.0-2         cowplot_1.1.1          webshot_0.5.3          Rtsne_0.16            
[117] uchardet_1.1.1         uwot_0.1.11            treemap_2.4-3          igraph_1.3.2          
[121] survival_3.3-1         plotrix_3.8-2          systemfonts_1.0.4      htmltools_0.5.2       
[125] memoise_2.0.1          modeltools_0.2-23      locfit_1.5-9.6         viridisLite_0.4.0     
[129] digest_0.6.29          assertthat_0.2.1       mime_0.12              rappdirs_0.3.3        
[133] repr_1.1.6             RSQLite_2.2.14         sqldf_0.4-11           future.apply_1.9.0    
[137] Exact_3.1              data.table_1.14.2      blob_1.2.3             splines_4.2.0         
[141] reactable_0.3.0        RCurl_1.98-1.7         hms_1.1.1              colorspace_2.0-3      
[145] base64enc_0.1-3        libcoin_1.0-9          Rcpp_1.0.8.3           coin_1.4-2            
[149] RANN_2.6.1             mvtnorm_1.1-3          multcompView_0.1-8     fansi_1.0.3           
[153] parallelly_1.32.0      IRdisplay_1.1          SnowballC_0.7.0        R6_2.5.1              
[157] grid_4.2.0             ggridges_0.5.3         lifecycle_1.0.1        rootSolve_1.8.2.3     
[161] zip_2.2.0              curl_4.3.2             leiden_0.4.2           Matrix_1.5-1          
[165] RcppAnnoy_0.0.19       TH.data_1.1-1          RColorBrewer_1.1-3     iterators_1.0.14      
[169] htmlwidgets_1.5.4      polyclip_1.10-0        purrr_0.3.4            rvest_1.0.2           
[173] mgcv_1.8-40            globals_0.15.0         lmom_2.9               patchwork_1.1.1       
[177] spatstat.random_2.2-0  progressr_0.10.1       codetools_0.2-18       gtools_3.9.2.2        
[181] prettyunits_1.1.1      dbplyr_2.2.0           gridBase_0.4-7         gtable_0.3.0          
[185] DBI_1.1.3              tensor_1.5             httr_1.4.3             KernSmooth_2.23-20    
[189] stringi_1.7.6          progress_1.2.2         msigdbr_7.5.1          reshape2_1.4.4        
[193] farver_2.1.0           uuid_1.1-0             annotate_1.74.0        xml2_1.3.3            
[197] admisc_0.28            boot_1.3-28            IRkernel_1.3.2         kableExtra_1.3.4      
[201] geneplotter_1.74.0     scattermore_0.8        bit_4.0.4              spatstat.data_3.0-0   
[205] pkgconfig_2.0.3        knitr_1.39   




