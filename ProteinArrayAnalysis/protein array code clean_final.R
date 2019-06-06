#' ---
#' title: "Protein Array Regional Meninges Ananlysis"
#' author: "Nathan Boles"
#' date: "November 14th, 2018"
#' output: pdf_document
#' ---

library(foreach)
library(doParallel)
library(ggplot2)
library(plyr)
library(psycho)
library(dplyr)
library(cluster)
library(reshape2)

# sessionInfo()
# R version 3.5.1 (2018-07-02)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1

# Matrix products: default

# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#   [1] reshape2_1.4.3    cluster_2.0.7-1   dplyr_0.7.7       psycho_0.3.7      plyr_1.8.4        ggplot2_3.1.0    
# [7] doParallel_1.0.14 iterators_1.0.10  foreach_1.4.4    

# loaded via a namespace (and not attached):
#   [1] backports_1.1.2        Hmisc_4.1-1            BDgraph_2.52           igraph_1.2.2           lazyeval_0.2.1        
# [6] splines_3.5.1          crosstalk_1.0.0        rstantools_1.5.1       inline_0.3.15          digest_0.6.18         
# [11] htmltools_0.3.6        matrixcalc_1.0-3       rsconnect_0.8.8        lmerTest_3.0-1         magrittr_1.5          
# [16] checkmate_1.8.5        sna_2.4                matrixStats_0.54.0     xts_0.11-1             prettyunits_1.0.2     
# [21] jpeg_0.1-8             sem_3.1-9              colorspace_1.3-2       callr_3.0.0            crayon_1.3.4          
# [26] lme4_1.1-18-1          bindr_0.1.1            survival_2.42-3        zoo_1.8-4              glue_1.3.0            
# [31] gtable_0.2.0           ppcor_1.1              emmeans_1.2.4          MatrixModels_0.4-1     mi_1.0                
# [36] pkgbuild_1.0.2         rstan_2.18.1           ggm_2.3                abind_1.4-5            scales_1.0.0          
# [41] mvtnorm_1.0-8          miniUI_0.1.1.1         Rcpp_0.12.19           xtable_1.8-3           htmlTable_1.12        
# [46] foreign_0.8-70         Formula_1.2-3          stats4_3.5.1           StanHeaders_2.18.0     DT_0.4                
# [51] htmlwidgets_1.3        threejs_0.3.1          RColorBrewer_1.1-2     lavaan_0.6-3           acepack_1.4.1         
# [56] pkgconfig_2.0.2        loo_2.0.0              manipulate_1.0.1       nnet_7.3-12            tidyselect_0.2.5      
# [61] rlang_0.3.0            later_0.7.5            ggcorrplot_0.1.2       munsell_0.5.0          tools_3.5.1           
# [66] cli_1.0.1              statnet.common_4.1.4   broom_0.5.0            ggridges_0.5.1         fdrtool_1.2.15        
# [71] stringr_1.3.1          arm_1.10-1             yaml_2.2.0             processx_3.2.0         knitr_1.20            
# [76] purrr_0.2.5            bindrcpp_0.2.2         glasso_1.10            pbapply_1.3-4          nlme_3.1-137          
# [81] whisker_0.3-2          mime_0.6               rstanarm_2.18.1        compiler_3.5.1         bayesplot_1.6.0       
# [86] shinythemes_1.1.1      rstudioapi_0.8         png_0.1-7              huge_1.2.7             tibble_1.4.2          
# [91] pbivnorm_0.6.0         DescTools_0.99.25      stringi_1.2.4          ps_1.2.0               qgraph_1.5            
# [96] lattice_0.20-35        Matrix_1.2-14          psych_1.8.4            nloptr_1.2.1           markdown_0.8          
# [101] shinyjs_1.0            pillar_1.3.0           estimability_1.3       data.table_1.11.8      corpcor_1.6.9         
# [106] httpuv_1.4.5           R6_2.3.0               latticeExtra_0.6-28    nFactors_2.3.3         MuMIn_1.42.1          
# [111] promises_1.0.1         network_1.13.0.1       gridExtra_2.3          BayesFactor_0.9.12-4.2 codetools_0.2-15      
# [116] boot_1.3-20            colourpicker_1.0       MASS_7.3-50            gtools_3.8.1           assertthat_0.2.0      
# [121] rjson_0.2.20           withr_2.1.2            shinystan_2.5.0        mnormt_1.5-5           expm_0.999-3          
# [126] grid_3.5.1             rpart_4.1-13           tidyr_0.8.1            coda_0.19-2            minqa_1.2.4           
# [131] d3Network_0.5.2.1      numDeriv_2016.8-1      shiny_1.1.0            base64enc_0.1-3        ellipse_0.4.1         
# [136] dygraphs_1.1.1.6  

cl <- makeCluster(1)
registerDoParallel(cl)

######functions########

###normalization from manual Xnorm=Xy*P(ref)/Py subarray 1 will be reference (block1)###

normalizer<-function(vec.char,pro.arrays,key){
  ####vec.char is a vector of array names e.g. "M3"###
  z<-foreach(b=1:length(vec.char), .combine='rbind') %do% {
    m3<-pro.arrays[which(pro.arrays[,6]==vec.char[b]),]
    m3k<-key[which(key[,4]==vec.char[b]),]
    x<-rep(as.character(m3k[,3]),max(m3$Block))
    m3<-data.frame(m3,as.character(x))
    x<-m3[which(m3[,1]==1),]
    P3<-mean(x[grep("pos3",x[,7],ignore.case=T),5])
    x<-m3[grep("pos3",m3[,7],ignore.case=T),5]
    x<-foreach(i=1:(length(x)/2), .combine='c') %do% {
      y<-seq(2,length(x),2)
      low<-y-1
      mean(x[low[i]],x[y[i]])
    }
    scalor<-P3/x
    m3.norm<-foreach(i=1:max(m3$Block),.combine='rbind') %do% {
      x<-m3[which(m3[,1]==i),]
      y<-x[,5]*scalor[i]
      data.frame(x,y)
    }
    m3.norm<-foreach(i=1:length(unique(m3.norm[,7])),.combine='rbind') %do% {
      x<-m3.norm[which(unique(m3.norm[,7])[i]==m3.norm[,7]),]
      hi<-seq(2,dim(x)[1],2)
      low<-hi-1
      p<-foreach(m=1:(length(unique(m3.norm[,1]))/2), .combine='c') %do% {
        mean(x[low[m],8],x[hi[m],8])
      }
      data.frame(x[low,],p)
    }
    
  }
  z<-z[,1:8]
  colnames(z)<-c(colnames(pro.arrays),"Key","Norm_value")
  z
}

###normalization from manual Xnorm=Xy*P(ref)/Py subarray 1 will be reference (block1)###
#####times a scalor based on concentration of protein put on array###
normalizer.concentration<-function(vec.char,pro.arrays,key,protein.concentration){
  conc.scale<-min(protein.concentration)/protein.concentration
  ####vec.char is a vector of array names e.g. "M3"###
  z<-foreach(b=1:length(vec.char), .combine='rbind') %do% {
    m3<-pro.arrays[which(pro.arrays[,6]==vec.char[b]),]
    m3k<-key[which(key[,4]==vec.char[b]),]
    x<-rep(as.character(m3k[,3]),max(m3$Block))
    m3<-data.frame(m3,as.character(x))
    x<-m3[which(m3[,1]==1),]
    P3<-mean(x[grep("pos3",x[,7],ignore.case=T),5])
    x<-m3[grep("pos3",m3[,7],ignore.case=T),5]
    x<-foreach(i=1:(length(x)/2), .combine='c') %do% {
      y<-seq(2,length(x),2)
      low<-y-1
      mean(x[low[i]],x[y[i]])
    }
    scalor<-P3/x
    m3.norm<-foreach(i=1:max(m3$Block),.combine='rbind') %do% {
      x<-m3[which(m3[,1]==i),]
      y<-x[,5]*scalor[i]*conc.scale[i]
      data.frame(x,y)
    }
    m3.norm<-foreach(i=1:length(unique(m3.norm[,7])),.combine='rbind') %do% {
      x<-m3.norm[which(unique(m3.norm[,7])[i]==m3.norm[,7]),]
      hi<-seq(2,dim(x)[1],2)
      low<-hi-1
      p<-foreach(m=1:(length(unique(m3.norm[,1]))/2), .combine='c') %do% {
        mean(x[low[m],8],x[hi[m],8])
      }
      data.frame(x[low,],p)
    }
    
  }
  z<-z[,1:8]
  colnames(z)<-c(colnames(pro.arrays),"Key","Norm_value")
  z
}
########analyze by manufacturer's suggestion

biotech.analyzer<-function(protein.arrays,group1,group2,sample.key,background.array,ncores=1){
  require(foreach)
  require(doParallel)
  require(stringr)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  groups<-sample.key[c(group1,group2)]
  norm.val.con<-foreach(i=1:length(unique(protein.arrays$Block)),.combine="cbind") %do% 
  {protein.arrays[which(protein.arrays$Block==unique(protein.arrays$Block[i])),8]}
  colnames(norm.val.con)<-names(sample.key)
  rownames(norm.val.con)<-protein.arrays[which(protein.arrays$Block==groups[1]),7]
  norm.val.con<-norm.val.con[-(grep("pos",rownames(norm.val.con))),]
  norm.val.con<-norm.val.con[-(grep("neg",rownames(norm.val.con))),]
  
  background<-mean(background.array[which(background.array$Key=="neg"),"Norm_value"])+2*
    sd(background.array[which(background.array$Key=="neg"),"Norm_value"])
  norm.val.b<-norm.val.con-background
  norm.val.b<-data.frame(round(apply(norm.val.b,2, function(z) replace(z,which(z<0),0)),0))
  z<-which(apply(norm.val.b[,names(groups)],1, function(b) length(which(b==0)))>2)
  norm.val.b<-norm.val.b[-z,]
  colnames(norm.val.b)<-names(sample.key)
  whole.set<-norm.val.b
  norm.val.b<-norm.val.b[,names(groups)]
  x<-data.frame(rownames(norm.val.b),norm.val.b[,group2[1]]/norm.val.b[,group1[1]],norm.val.b[,group2[2]]/norm.val.b[,group1[2]])
  z<-foreach(i=1:length(group1),.combine='c') %do% {paste(group1[i],group2[i],sep=".vs.")}
  colnames(x)<-c("proteins",z)
  
  FC<- apply(x[,2:3],1,mean)
  if(length(which(is.na(FC))) >0) {
    x<-x[-(which(is.na(FC))),]
    FC<-FC[-(which(is.na(FC)))]} else {}
  
  Directionality<-foreach(i=1:dim(x)[1],.combine='rbind') %do% {
    y<-if(x[i,2]>1 & x[i,3]>1) {"Same"} else {
      if(x[i,2]<1 & x[i,3]<1) {"Same"} else {"Different"}
    }
  }
  x<-data.frame(x,FC,Directionality)
  y<-x[which(x$Directionality=="Same"),]
  ymx<-y[which(y$FC>=1.5),]
  ymn<-y[which(y$FC<=0.65),]
  biotec.diff<-rbind(ymx,ymn)
  Up_in<-sapply(biotec.diff$FC, function(b) if(b>1){str_extract(group2[1],"[aA-zZ]+")} else {str_extract(group1[1],"[aA-zZ]+")})
  biotec.diff<-data.frame(whole.set[as.character(biotec.diff$proteins),],biotec.diff,Up_in)
}

########write tables#######
set.tables<-function(data1,data2,cols) {
  interd<-data.frame(data1[intersect(data1$proteins,data2$proteins),],data2[intersect(data1$proteins,data2$proteins),cols])
  diffd1<-data1[setdiff(data1$proteins,data2$proteins),]
  diffd2<-data2[setdiff(data2$proteins,data1$proteins),]
  write.csv(interd,paste("Intersection_of_",deparse(substitute(data1)),"_in_comparison_",deparse(substitute(data2)),".csv",sep=""))
  write.csv(diffd1,paste("unique_to_",deparse(substitute(data1)),"_in_comparison_",deparse(substitute(data1)),"_to_",deparse(substitute(data2)),".csv",sep=""))
  write.csv(diffd2,paste("unique_to_",deparse(substitute(data2)),"_in_comparison_",deparse(substitute(data1)),"_to_",deparse(substitute(data2)),".csv",sep=""))
  write.csv(data1,paste(deparse(substitute(data1)),".csv",sep=""))
  write.csv(data2,paste(deparse(substitute(data2)),".csv",sep=""))
}



########data###
pro.arrays<-data.frame(read.csv("R_protein_array.csv"))
key<-data.frame(read.csv("key_protein array.csv"))
###########sample key
sample.key<-1:12
names(sample.key)<-c("AmAc1","PmPc1","AmAc2","PmPc2","Ac1","Pc1","Ac2","Media", "Am1","Am2","Pm1","Pm2")
########protein concentraions put on array##
conc<-c(536.4,1149.6,453,799.68,456.6,456.4,436.4,455.6,757.2,485.1,1200,530.7)
names(conc)<-names(sample.key)

##########comparisons and normalization

x<-c("M3","M4","M5")
norm.pro<-normalizer(x,pro.arrays,key)
norm.pro.con<-normalizer.concentration(x,pro.arrays,key,conc)


AmAc.PmPC<-biotech.analyzer(norm.pro,c("AmAc1","AmAc2"),c("PmPc1","PmPc2"),sample.key,norm.pro)
AmAc.Ac<-biotech.analyzer(norm.pro,c("AmAc1","AmAc2"),c("Ac1","Ac2"),sample.key,norm.pro)
Am.Pm<-biotech.analyzer(norm.pro,c("Am1","Am2"),c("Pm1","Pm2"),sample.key,norm.pro)
AmAc.Am<-biotech.analyzer(norm.pro,c("AmAc1","AmAc2"),c("Am1","Am2"),sample.key,norm.pro)
PmPc.Pm<-biotech.analyzer(norm.pro,c("PmPc1","PmPc2"),c("Pm1","Pm2"),sample.key,norm.pro)
Ac.Am<-biotech.analyzer(norm.pro,c("Ac1","Ac2"),c("Am1","Am2"),sample.key,norm.pro)

AmAc.PmPC.conc<-biotech.analyzer(norm.pro.con,c("AmAc1","AmAc2"),c("PmPc1","PmPc2"),sample.key,norm.pro)
AmAc.Ac.conc<-biotech.analyzer(norm.pro.con,c("AmAc1","AmAc2"),c("Ac1","Ac2"),sample.key,norm.pro)
Am.Pm.conc<-biotech.analyzer(norm.pro.con,c("Am1","Am2"),c("Pm1","Pm2"),sample.key,norm.pro)
AmAc.Am.conc<-biotech.analyzer(norm.pro.con,c("AmAc1","AmAc2"),c("Am1","Am2"),sample.key,norm.pro)
PmPc.Pm.conc<-biotech.analyzer(norm.pro.con,c("PmPc1","PmPc2"),c("Pm1","Pm2"),sample.key,norm.pro)
Ac.Am.conc<-biotech.analyzer(norm.pro.con,c("Ac1","Ac2"),c("Am1","Am2"),sample.key,norm.pro)

x<-intersect(rownames(Am.Pm),rownames(AmAc.PmPC))
y<-setdiff(rownames(Am.Pm),rownames(AmAc.PmPC))
z<-setdiff(rownames(AmAc.PmPC),rownames(Am.Pm))

x<-intersect(rownames(AmAc.Ac),rownames(AmAc.Am))
y<-intersect(x,rownames(AmAc.PmPC))
z<-data.frame(AmAc.PmPC[y,],AmAc.Ac[y,14:18],AmAc.Am[y,14:18])
write.csv(z,"intersection_AmAC.PMPC_vs_AmAc.Ac_vs_AmAc.Am.csv")
x<-intersect(rownames(AmAc.PmPC),rownames(AmAc.Am))
x<-intersect(rownames(AmAc.Ac),rownames(AmAc.PmPC))
m<-intersect(z$proteins,y)



View(Ac.Am)
View(PmPc.Pm)
View(AmAc.PmPC)
View(Am.Pm)
View(AmAc.Ac)
View(AmAc.Am)


View(AmAc.Am.conc)
View(AmAc.Ac.conc)
View(Am.Pm.conc)
View(Ac.Am.conc)
View(PmPc.Pm.conc)
View(AmAc.PmPC.conc)

#####aging arrays####
load("Protein.arrays.clean.RData")
pro.arrays.aged<-data.frame(read.csv("Protein.array.aging.csv"))

sample.aged.key<-1:8
names(sample.aged.key)<-c("AmAc1","PmPc1","AmAc2","PmPc2", "Am1","Pm1","Am2","Pm2")
x<-c("M3","M4","M5")
norm.pro.aged<-normalizer(x,pro.arrays.aged,key)

aged.AmAc.PmPC<-biotech.analyzer(norm.pro.aged,c("AmAc1","AmAc2"),c("PmPc1","PmPc2"),sample.aged.key,norm.pro.aged)
aged.Am.Pm<-biotech.analyzer(norm.pro.aged,c("Am1","Am2"),c("Pm1","Pm2"),sample.aged.key,norm.pro.aged)
aged.AmAc.Am<-biotech.analyzer(norm.pro.aged,c("AmAc1","AmAc2"),c("Am1","Am2"),sample.aged.key,norm.pro.aged)
aged.PmPc.Pm<-biotech.analyzer(norm.pro.aged,c("PmPc1","PmPc2"),c("Pm1","Pm2"),sample.aged.key,norm.pro.aged)

x<-norm.pro[which(as.numeric(norm.pro$Block)==5),]
x$Block<-rep(9,dim(x)[1])
y<-norm.pro[which(as.numeric(norm.pro$Block)==7),]
y$Block<-rep(10,dim(y)[1])

norm.pro.aged<-rbind(norm.pro.aged,x,y)
sample.aged.key<-1:10
names(sample.aged.key)<-c("AmAc1","PmPc1","AmAc2","PmPc2", "Am1","Pm1","Am2","Pm2","Ac1","Ac2")

aged.AmAc.AC<-biotech.analyzer(norm.pro.aged,c("AmAc1","AmAc2"),c("Ac1","Ac2"),sample.aged.key,norm.pro.aged)
aged.Am.AC<-biotech.analyzer(norm.pro.aged,c("Am1","Am2"),c("Ac1","Ac2"),sample.aged.key,norm.pro.aged)




write.csv(aged.Am.AC,"aged.Am.vs.Ac.csv")
write.csv(aged.AmAc.AC,"aged.AmAc.vs.Ac.csv")

save.image("Protein arrays finished.RData")

#########normalize all arrays together


unique(pro.arrays.aged$Block)


x<-pro.arrays.aged
for(i in 1:8 ) {x$Block[which(x$Block==i)]<-i+12}

all.arrays<-rbind(pro.arrays,x)
x<-13:20
names(x)<-c("AmAc1.18","PmPc1.18","AmAc2.18","PmPc2.18","Am1.18","Pm1.18","Am2.18","Pm2.18")
all.key<-c(sample.key,x)
  
  
  

x<-c("M3","M4","M5")
norm.pro.all<-normalizer(x,all.arrays,key)


yvo.Am<-biotech.analyzer(norm.pro.all,c("Am1","Am2"),c("Am1.18","Am2.18"),all.key,norm.pro.all)
yvo.AmAc<-biotech.analyzer(norm.pro.all,c("AmAc1","AmAc2"),c("AmAc1.18","AmAc2.18"),all.key,norm.pro.all)
yvo.PmPc<-biotech.analyzer(norm.pro.all,c("PmPc1","PmPc2"),c("PmPc1.18","PmPc2.18"),all.key,norm.pro.all)
write.csv(yvo.Am,"Young.vs.old.in.AM.csv")
write.csv(yvo.AmAc,"Young.vs.old.in.AmAc.csv")
write.csv(yvo.PmPc,"Young.vs.old.in.PmPc.csv")


#########make heatmaps and supplemental tables for paper
heatmap.prep<-function(protein.arrays,group1,group2,sample.key,background.array,ncores=1){
  require(foreach)
  require(doParallel)
  require(stringr)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  groups<-sample.key[c(group1,group2)]
  norm.val.con<-foreach(i=1:length(unique(protein.arrays$Block)),.combine="cbind") %do% 
  {protein.arrays[which(protein.arrays$Block==unique(protein.arrays$Block[i])),8]}
  colnames(norm.val.con)<-names(sample.key)
  rownames(norm.val.con)<-protein.arrays[which(protein.arrays$Block==groups[1]),7]
  norm.val.con<-norm.val.con[-(grep("pos",rownames(norm.val.con))),]
  norm.val.con<-norm.val.con[-(grep("neg",rownames(norm.val.con))),]
  
  background<-mean(background.array[which(background.array$Key=="neg"),"Norm_value"])+2*
    sd(background.array[which(background.array$Key=="neg"),"Norm_value"])
  norm.val.b<-norm.val.con-background
  norm.val.b<-data.frame(round(apply(norm.val.b,2, function(z) replace(z,which(z<0),0)),0))
  colnames(norm.val.b)<-names(sample.key)
  whole.set<-norm.val.b
  norm.val.b<-norm.val.b[,names(groups)]
  x<-data.frame(rownames(norm.val.b),(norm.val.b[,group2[1]]+1)/(norm.val.b[,group1[1]]+1),
                (norm.val.b[,group2[2]]+1)/(norm.val.b[,group1[2]]+1))
  z<-foreach(i=1:length(group1),.combine='c') %do% {paste(group1[i],group2[i],sep=".vs.")}
  colnames(x)<-c("proteins",z)
  
  FC<- apply(x[,2:3],1,mean)
  if(length(which(is.na(FC))) >0) {
    x<-x[-(which(is.na(FC))),]
    FC<-FC[-(which(is.na(FC)))]} else {}
  
  Directionality<-foreach(i=1:dim(x)[1],.combine='rbind') %do% {
    y<-if(x[i,2]>1 & x[i,3]>1) {"Same"} else {
      if(x[i,2]<1 & x[i,3]<1) {"Same"} else {"Different"}
    }
  }
  x<-data.frame(x,FC,Directionality)
}


z1<-heatmap.prep(norm.pro,c("AmAc1","AmAc2"),c("PmPc1","PmPc2"),sample.key,norm.pro)
z2<-heatmap.prep(norm.pro,c("Am1","Am2"),c("Pm1","Pm2"),sample.key,norm.pro)

z3<-heatmap.prep(norm.pro.aged,c("AmAc1","AmAc2"),c("PmPc1","PmPc2"),sample.aged.key,norm.pro.aged)
z4<-heatmap.prep(norm.pro.aged,c("Am1","Am2"),c("Pm1","Pm2"),sample.aged.key,norm.pro.aged)

m<-c("CC1","CC2")
m2<-c("M1","M2")

z1.1<-foreach(i=1:2,.combine='rbind') %do% {
  x<-data.frame(z1[,1],as.numeric(z1[,i+1]),rep(m[i],dim(z1)[1]),rep("3",dim(z1)[1]))
  colnames(x)<-c("Protein","FC","Samples","Age")
  return(x)
}

z2.1<-foreach(i=1:2,.combine='rbind') %do% {
  x<-data.frame(z2[,1],as.numeric(z2[,i+1]),rep(m2[i],dim(z2)[1]),rep("3",dim(z2)[1]))
  colnames(x)<-c("Protein","FC","Samples","Age")
  return(x)
}

z3.1<-foreach(i=1:2,.combine='rbind') %do% {
  x<-data.frame(z3[,1],as.numeric(z3[,i+1]),rep(m[i],dim(z3)[1]),rep("18",dim(z3)[1]))
  colnames(x)<-c("Protein","FC","Samples","Age")
  return(x)
}

z4.1<-foreach(i=1:2,.combine='rbind') %do% {
  x<-data.frame(z4[,1],as.numeric(z4[,i+1]),rep(m2[i],dim(z4)[1]),rep("18",dim(z4)[1]))
  colnames(x)<-c("Protein","FC","Samples","Age")
  return(x)
}

z.heat<-rbind(z1.1,z2.1,z3.1,z4.1)
z.heat$logFC<-log(z.heat$FC)


heat.3<-z.heat[z.heat$Age==3,]
heat.18<-z.heat[z.heat$Age==18,]


x<-as.character(read.csv("heatmaps_reg_3_1_list.csv",header=F)[,1])
heat.3<-foreach(i=1:length(x),.combine='rbind') %do% {
  blue<-which(as.character(heat.3$Protein)==x[i])
  fin<-heat.3[blue,]
}
heat.3<-rbind(heat.3[which(heat.3$Samples=="M1"),],heat.3[which(heat.3$Samples=="M2"),])

x<-as.character(read.csv("heatmaps_reg_18_2_list.csv",header=F)[,1])
heat.18<-foreach(i=1:length(x),.combine='rbind') %do% {
  blue<-which(as.character(heat.18$Protein)==x[i])
  fin<-heat.18[blue,]
}
heat.18<-rbind(heat.18[which(heat.18$Samples=="M1"),],heat.18[which(heat.18$Samples=="M2"),])

x<-as.character(read.csv("heatmaps_reg_comp_3_list.csv",header=F)[,1])
heat.comp<-foreach(i=1:length(x),.combine='rbind') %do% {
  blue<-which(as.character(z.heat$Protein)==x[i])
  fin<-z.heat[blue,]
}
heat.comp<-rbind(heat.comp[which(heat.comp$Samples=="M1"),],heat.comp[which(heat.comp$Samples=="M2"),])


ggplot(heat.3, aes(Samples, Protein)) + geom_tile(aes(fill = log(FC) )) + 
  scale_fill_gradient2(low="midnightblue",mid="mistyrose", high="tan4") 

ggplot(heat.18, aes(Samples, Protein)) + geom_tile(aes(fill = log(FC) )) + 
  scale_fill_gradient2(low="midnightblue",mid="mistyrose", high="tan4") 

ggplot(heat.comp, aes(Age, Protein)) + geom_tile(aes(fill = log(FC) )) + 
  scale_fill_gradient2(low="midnightblue",mid="mistyrose", high="tan4") 

#####function to compare across ages and make heatmap
heatmap.prep.age<-function(protein.arrays1,protein.arrays2,group1,group2,sample.key1,sample.key2,age1,age2,ncores=1){
  require(foreach)
  require(doParallel)
  require(stringr)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  group1<-sample.key[group1]
  group2<-sample.key2[group2]
  norm.val.con1<-foreach(i=1:length(unique(protein.arrays1$Block)),.combine="cbind") %do% {
    protein.arrays1[which(protein.arrays1$Block==unique(protein.arrays1$Block[i])),8]
  }
  colnames(norm.val.con1)<-names(sample.key)
  rownames(norm.val.con1)<-protein.arrays1[which(protein.arrays1$Block==group1[1]),7]
  norm.val.con1<-norm.val.con1[-(grep("pos",rownames(norm.val.con1))),]
  norm.val.con1<-norm.val.con1[-(grep("neg",rownames(norm.val.con1))),]
  
  norm.val.con2<-foreach(i=1:length(unique(protein.arrays2$Block)),.combine="cbind") %do% {
    protein.arrays2[which(protein.arrays2$Block==unique(protein.arrays2$Block[i])),8]
  }
  colnames(norm.val.con2)<-names(sample.key2)
  rownames(norm.val.con2)<-protein.arrays2[which(protein.arrays2$Block==group2[1]),7]
  norm.val.con2<-norm.val.con2[-(grep("pos",rownames(norm.val.con2))),]
  norm.val.con2<-norm.val.con2[-(grep("neg",rownames(norm.val.con2))),]
  
  
  
  background<-mean(protein.arrays1[which(protein.arrays1$Key=="neg"),"Norm_value"])+2*
    sd(protein.arrays1[which(protein.arrays1$Key=="neg"),"Norm_value"])
  background2<-mean(protein.arrays2[which(protein.arrays2$Key=="neg"),"Norm_value"])+2*
    sd(protein.arrays2[which(protein.arrays2$Key=="neg"),"Norm_value"])
  
  norm.val.b<-norm.val.con1-background
  norm.val.b<-data.frame(round(apply(norm.val.b,2, function(z) replace(z,which(z<0),0)),0))
  colnames(norm.val.b)<-names(sample.key1)
  whole.set1<-norm.val.b
  
  norm.val.b2<-norm.val.con2-background2
  norm.val.b2<-data.frame(round(apply(norm.val.b2,2, function(z) replace(z,which(z<0),0)),0))
  colnames(norm.val.b2)<-names(sample.key2)
  whole.set2<-norm.val.b2
  
  norm.val.b<-norm.val.b[,names(group1)]+1
  norm.val.b<-apply(norm.val.b,1,mean)
  norm.val.b2<-norm.val.b2[,names(group2)]+1
  norm.val.b2<-apply(norm.val.b2,1,mean)
  x<-norm.val.b/norm.val.b2
  
  return(x)
}



z5<-heatmap.prep.age(norm.pro,norm.pro.aged,c("AmAc1","AmAc2"),c("AmAc1","AmAc2"),sample.key,sample.aged.key,"3","18")
z6<-heatmap.prep.age(norm.pro,norm.pro.aged,c("PmPc1","PmPc2"),c("PmPc1","PmPc2"),sample.key,sample.aged.key,"3","18")
z7<-heatmap.prep.age(norm.pro,norm.pro.aged,c("Am1","Am2"),c("Am1","Am2"),sample.key,sample.aged.key,"3","18")
z8<-heatmap.prep.age(norm.pro,norm.pro.aged,c("Pm1","Pm2"),c("Pm1","Pm2"),sample.key,sample.aged.key,"3","18")


z5.1<-data.frame(names(z5),z5,rep("CoCulture",length(z5)),rep("Anterior",length(z5)))
colnames(z5.1)<-c("Protein","FC","Samples","Region")


z6.1<-data.frame(names(z6),z6,rep("CoCulture",length(z6)),rep("Posterior",length(z6)))
colnames(z6.1)<-c("Protein","FC","Samples","Region")

z7.1<-data.frame(names(z7),z7,rep("Mngs Only",length(z7)),rep("Anterior",length(z7)))
colnames(z7.1)<-c("Protein","FC","Samples","Region")

z8.1<-data.frame(names(z8),z8,rep("Mngs Only",length(z8)),rep("Posterior",length(z8)))
colnames(z8.1)<-c("Protein","FC","Samples","Region")

z.age<-rbind(z5.1,z6.1,z7.1,z8.1)
z.age$logFC<-log(as.numeric(z.age$FC))

x<-as.character(read.csv("heatmaps_age_3_1_list.csv",header=F)[,1])
age.1<-foreach(i=1:length(x),.combine='rbind') %do% {
  blue<-which(as.character(z.age$Protein)==x[i])
  fin<-z.age[blue,]
}
age.1<-age.1[which(age.1$Samples=="Mngs Only"),]

x<-as.character(read.csv("heatmaps_age_18_2_list.csv",header=F)[,1])
age.2<-foreach(i=1:length(x),.combine='rbind') %do% {
  blue<-which(as.character(z.age$Protein)==x[i])
  fin<-z.age[blue,]
}
age.2<-age.2[which(age.2$Samples=="Mngs Only"),]

x<-as.character(read.csv("heatmaps_age_comp_3_list.csv",header=F)[,1])
age.3<-foreach(i=1:length(x),.combine='rbind') %do% {
  blue<-which(as.character(z.age$Protein)==x[i])
  fin<-z.age[blue,]
}
age.3<-age.3[which(age.3$Samples=="Mngs Only"),]

ggplot(age.1, aes(Samples, Protein)) + geom_tile(aes(fill = logFC )) + 
  scale_fill_gradient2(low="darkorchid4",mid="slategray1", high="red") +  facet_grid(cols=vars(Region)) +
  theme(axis.text.y=element_text(size=8))

ggplot(age.2, aes(Samples, Protein)) + geom_tile(aes(fill = logFC )) + 
  scale_fill_gradient2(low="darkorchid4",mid="slategray1", high="red") +  facet_grid(cols=vars(Region)) +
  theme(axis.text.y=element_text(size=8))

ggplot(age.3, aes(Samples, Protein)) + geom_tile(aes(fill = logFC )) + 
  scale_fill_gradient2(low="darkorchid4",mid="slategray1", high="red") +  facet_grid(cols=vars(Region)) +
  theme(axis.text.y=element_text(size=8))

x<-intersect(rownames(Am.Pm),rownames(AmAc.PmPC))
y<-intersect(rownames(aged.Am.Pm),rownames(aged.AmAc.PmPC))
m<-intersect(x,y)



ggplot(z.age[m,], aes(Samples, Protein)) + geom_tile(aes(fill = logFC )) + 
  scale_fill_gradient2(low="blue",mid="lightgrey", high="red",limits=c(-8,8)) +  facet_grid(cols=vars(Region)) +
  theme(axis.text.y=element_text(size=8))

x<-foreach(i=1:length(unique(z.heat$Protein)),.combine='rbind') %do% {
  z<-z.heat[z.heat$Protein==z.heat$Protein[i],]
  dcast(z,Protein~Samples+Age,value.var="FC")
}

z5<-heatmap.prep.age(norm.pro,norm.pro.aged,c("AmAc1","AmAc2"),c("AmAc1","AmAc2"),sample.key,sample.aged.key,"3","18")
z6<-heatmap.prep.age(norm.pro,norm.pro.aged,c("PmPc1","PmPc2"),c("PmPc1","PmPc2"),sample.key,sample.aged.key,"3","18")
z7<-heatmap.prep.age(norm.pro,norm.pro.aged,c("Am1","Am2"),c("Am1","Am2"),sample.key,sample.aged.key,"3","18")
z8<-heatmap.prep.age(norm.pro,norm.pro.aged,c("Pm1","Pm2"),c("Pm1","Pm2"),sample.key,sample.aged.key,"3","18")

x<-data.frame(z5,z6,z7,z8)
colnames(x)<-c("AmAc","PMPC","AM","PM")

z1<-heatmap.prep(norm.pro,c("AmAc1","AmAc2"),c("PmPc1","PmPc2"),sample.key,norm.pro)
z2<-heatmap.prep(norm.pro,c("Am1","Am2"),c("Pm1","Pm2"),sample.key,norm.pro)

z3<-heatmap.prep(norm.pro.aged,c("AmAc1","AmAc2"),c("PmPc1","PmPc2"),sample.aged.key,norm.pro.aged)
z4<-heatmap.prep(norm.pro.aged,c("Am1","Am2"),c("Pm1","Pm2"),sample.aged.key,norm.pro.aged)

write.csv(z1,"Coculture young.csv")
write.csv(z2,"Meninges young.csv")
write.csv(z3,"Coculture old.csv")
write.csv(z4,"Meninges old.csv")


