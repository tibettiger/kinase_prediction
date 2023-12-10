

###plot RMSE,MAE,R2 benchmark results

library(ggplot2)
library(ggsci)
##barplot
cancer_type<-c("BC","LSCC","UCEC","GBM","HCC")
dat<-as.data.frame(matrix(0,5,3))
rownames(dat)<-cancer_type
colnames(dat)<-c(">0.5","0.3-0.5","<0.3")
for(i in 1:5){
    res<-read.table(paste("C:/users/90410/desktop/kinase/iscience_revision_10.9/result/",cancer_type[i],"/","XGBoost_result.csv",sep=""),sep=",")
    res<-res[-1,]
    for(j in 1:nrow(res)){
      if(res[j,2]>0.5)dat[i,1]<-dat[i,1]+1
      if(res[j,2]<0.3)dat[i,3]<-dat[i,3]+1
      if((res[j,2]>=0.3)&&(res[j,2]<=0.5))dat[i,2]<-dat[i,2]+1
    }
}
library(reshape2)
dat$cancer_type<-cancer_type
newres<-melt(dat)
newres$cancer_type<-factor(newres$cancer_type, levels=c('HCC', 'BC', 'LSCC', 'UCEC','GBM'))
ggplot(newres,aes(x=cancer_type,weight=value)) + geom_bar(aes(fill=variable),position = position_stack(reverse = TRUE),width=0.9)+
  scale_fill_lancet()+
  labs(x="Cancer type",y="Kinase number")+
  theme(title=element_text(size=8),axis.text = element_text(size = 8),axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),legend.title = element_text(size = 7),axis.text.x=element_text(angle=40,hjust=1,vjust=1))+
  theme_classic()

sample_size<-c(122,202,115,108,318)
for(i in 1:5){
  dat[i,5]<-dat[i,1]/sum(dat[i,1:3])
}

dat$sample_size<-sample_size
ggplot(dat, aes(x = sample_size, y = sum)) +
  geom_point() +
  xlab("Sample size") +
  ylab("Proportion") +
  stat_smooth(method = lm)+
  theme_classic()
cor(dat$sum,dat$sample_size,method='pearson')
pearson=0.643057


###benchmark kinase predict results(representative kinase)
### select ARAF  XGboost performance more consistent among different cancer types
model<-c("XGBoost","RandomForest","linear_regression","KNeighbors")
R2<-as.data.frame(matrix(NA,5,4))
RMSE<-as.data.frame(matrix(NA,5,4))
MAE<-as.data.frame(matrix(NA,5,4))
rownames(R2)<-cancer_type
colnames(R2)<-model
rownames(RMSE)<-cancer_type
colnames(RMSE)<-model
rownames(MAE)<-cancer_type
colnames(MAE)<-model
for(i in 1:5){
  for(j in 1:4){
  dat<-read.table(paste("C:/users/90410/desktop/kinase/iscience_revision_10.9/result/",cancer_type[i],"/",model[j],"_result.csv",sep=""),
                  sep=",")
  num<-which(dat[,1]=="CSNK1E")
  R2[i,j]<-dat[num,2]
  RMSE[i,j]<-dat[num,3]
  MAE[i,j]<-dat[num,4]
  }
}

R2$cancer_type<-rownames(R2)
newR2<-melt(R2)
ggplot(newR2,aes(x=cancer_type,y=value,fill=variable)) + geom_bar(position="dodge",stat="identity",width=0.9)+
  scale_fill_npg()+
  labs(x="Cancer type",y="R2")+
  theme(title=element_text(size=8),axis.text = element_text(size = 8),axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),legend.title = element_text(size = 7),axis.text.x=element_text(angle=40,hjust=1,vjust=1))+
  theme_classic()

RMSE$cancer_type<-rownames(RMSE)
newRMSE<-melt(RMSE)
ggplot(newRMSE,aes(x=cancer_type,y=value,fill=variable)) + geom_bar(position="dodge",stat="identity",width=0.9)+
  scale_fill_npg()+
  labs(x="Cancer type",y="RMSE")+
  theme(title=element_text(size=8),axis.text = element_text(size = 8),axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),legend.title = element_text(size = 7),axis.text.x=element_text(angle=40,hjust=1,vjust=1))+
  theme_classic()

MAE$cancer_type<-rownames(MAE)
newMAE<-melt(MAE)
ggplot(newMAE,aes(x=cancer_type,y=value,fill=variable)) + geom_bar(position="dodge",stat="identity",width=0.9)+
  scale_fill_npg()+
  labs(x="Cancer type",y="MAE")+
  theme(title=element_text(size=8),axis.text = element_text(size = 8),axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),legend.title = element_text(size = 7),axis.text.x=element_text(angle=40,hjust=1,vjust=1))+
  theme_classic()
##plot relationship
kinase<-c("TNK2","CSNK1E","SRPK1","PAK1","ABL1")
substrate_RNA<-c("AKT1","DDX3X","ASF1A","BAD","CAT")
RNA_exp<-read.table("C:/users/90410/desktop/KBPRNA/RNA/BC_RNA.csv",sep=",",header=TRUE)
kinase_act<-read.table("C:/users/90410/desktop/KBPRNA/kinase/BC_kinase.csv",sep=",",header=TRUE)
mat<-as.data.frame(matrix(NA,5,5))
rownames(mat)<-kinase
colnames(mat)<-substrate_RNA
for(i in 1:5){
  for(j in 1:5){
    num.kin<-which(kinase_act[,1]==kinase[i])
    num.RNA<-which(RNA_exp[,1]==substrate_RNA[j])
    mat[i,j]<-cor(as.numeric(kinase_act[num.kin,-1]),as.numeric(RNA_exp[num.RNA,-1]),method='pearson')
  }
}
corrplot(as.matrix(mat), method="color",outline = TRUE,tl.srt =45, tl.col = "black", tl.offset = 0.9, tl.cex = 0.9, cl.pos = 'r',number.font = 3, 
         col = col2,rect.col="white",rect.lwd=1,addrect=T,addgrid.col="white",cl.cex = 1,is.corr=FALSE,col.lim = c(-1,1),diag=TRUE)

###Plot feature importance of ARAF, ABL1, CSNK1E
##ARAF ABL1 CSNK1E             BC LSCC UCEC HCC GBM
#[(0.41473556, 'CDC25A'), (0.1310003, 'AURKB'), (0.03777455, 'MTERFD1'), (0.03501809, 'PLK1'), (0.03316619, 'EPHA3'), (0.028539868, 'MMP2'), (0.020673037, 'TIMP2'), (0.020330023, 'SUPV3L1'), (0.018990455, 'PSRC1'), (0.017968195, 'BMP4')]
#[(0.29053745, 'AURKB'), (0.15287578, 'CDC25A'), (0.08361147, 'CDH3'), (0.06814803, 'XBP1'), (0.053229537, 'CCND1'), (0.040163457, 'STAT1'), (0.02696869, 'ME2'), (0.023857126, 'NFIL3'), (0.022571746, 'E2F2'), (0.021589372, 'CDK1')]
#[(0.30987585, 'AURKB'), (0.16023886, 'CDC25A'), (0.07823381, 'CDK1'), (0.056485146, 'MVP'), (0.051649284, 'WASF3'), (0.041970734, 'CDH3'), (0.039111547, 'AURKA'), (0.03141794, 'PSIP1'), (0.024484519, 'UBE2C'), (0.02318979, 'XBP1')]

#[(0.74285966, 'UBE2C'), (0.034461495, 'E2F2'), (0.03144458, 'PYCR1'), (0.016789252, 'GAPDH'), (0.014231387, 'TLR4'), (0.012489407, 'POLR1C'), (0.010401122, 'PRSS23'), (0.010019144, 'TBX2'), (0.009295753, 'MELK'), (0.0075443764, 'WASF3')]
#[(0.6572201, 'UBE2C'), (0.13763322, 'KIF20A'), (0.03661761, 'EIF4EBP1'), (0.015945842, 'ADRB2'), (0.015158403, 'GAPDH'), (0.011750618, 'SOCS2'), (0.011443579, 'E2F2'), (0.00879897, 'HPRT1'), (0.008120913, 'ZMIZ1'), (0.0071978085, 'SLC27A3')]
#[(0.57062113, 'UBE2C'), (0.14459789, 'CCNB1'), (0.073502235, 'TIMELESS'), (0.073125, 'ORC1'), (0.03802475, 'CDC25A'), (0.019545475, 'TSPAN4'), (0.009967613, 'SOCS2'), (0.008943156, 'MYLK'), (0.008541398, 'TOP2A'), (0.00620769, 'SH3BP5')]

#[(0.25651672, 'RGS2'), (0.24838632, 'TPM1'), (0.07631669, 'PLK1'), (0.0715013, 'GHR'), (0.06856772, 'SNCA'), (0.041339662, 'TUBB6'), (0.040686417, 'PKIG'), (0.024792679, 'CSRP1'), (0.023739427, 'TCEAL4'), (0.022435332, 'CDK1')]
#[(0.46481207, 'SNCA'), (0.122242555, 'GHR'), (0.0718489, 'CGRRF1'), (0.0527083, 'ADRB2'), (0.037955273, 'LBR'), (0.035464842, 'ILK'), (0.029195663, 'DECR1'), (0.019791374, 'CRTAP'), (0.016501462, 'DCTD'), (0.014412658, 'E2F2')]
#[(0.3986306, 'RGS2'), (0.14690648, 'ILK'), (0.0857814, 'GHR'), (0.08337837, 'TCEAL4'), (0.048691712, 'TUBB6'), (0.048340164, 'RAB31'), (0.018394837, 'TSC22D3'), (0.018026903, 'CHEK1'), (0.016429896, 'MBOAT7'), (0.014341006, 'MYLK')]

#[(0.5637545, 'PHKB'), (0.14779642, 'AGL'), (0.050475936, 'PGM1'), (0.049929515, 'MSRA'), (0.020391224, 'TOP2A'), (0.017430492, 'PNP'), (0.013138151, 'ACAA1'), (0.010321407, 'IGF2BP2'), (0.008850343, 'BTK'), (0.008708257, 'GRWD1')]
#[(0.28922135, 'PGM1'), (0.28590098, 'PHKB'), (0.15005578, 'MSRA'), (0.041907303, 'PXMP2'), (0.039024502, 'GSTZ1'), (0.035707187, 'IGF2BP2'), (0.025131378, 'FAH'), (0.014648802, 'PYGL'), (0.013296606, 'CCDC86'), (0.011532077, 'PHKG2')]
#[(0.28102395, 'MSRA'), (0.23868923, 'AGL'), (0.12962556, 'MYLK'), (0.050943345, 'PGM1'), (0.03959517, 'GRWD1'), (0.031581316, 'IGF2BP2'), (0.02682722, 'CNDP2'), (0.019100044, 'NOSIP'), (0.016362512, 'TOP2A'), (0.012670117, 'ELAVL1')]

#[(0.30465952, 'SNAP25'), (0.06558413, 'CHAC1'), (0.055476196, 'SMARCA4'), (0.053609043, 'UBE2L6'), (0.040195175, 'TRAM2'), (0.038968224, 'CDK6'), (0.032269247, 'ALDOA'), (0.030546827, 'RRP8'), (0.028451888, 'ALDOC'), (0.02447965, 'NT5DC2')]
#[(0.16867863, 'ETS1'), (0.116426155, 'SNAP25'), (0.09395557, 'STXBP1'), (0.067075916, 'KIF20A'), (0.05638814, 'PYCR1'), (0.054839253, 'BIRC5'), (0.05477168, 'LAMA3'), (0.05430508, 'NOL3'), (0.035329852, 'CCND1'), (0.024148282, 'PRSS23')]
#[(0.23690271, 'SNAP25'), (0.2073464, 'SYNGR3'), (0.054448072, 'ALDOC'), (0.053817265, 'EIF4EBP1'), (0.037224494, 'GDPD5'), (0.030844353, 'IGFBP3'), (0.02402465, 'FUT1'), (0.0223478, 'ADI1'), (0.021476528, 'CHAC1'), (0.019006066, 'SNCA')]
library(readxl)
dat<-read_xlsx("C:/users/90410/desktop/kinase/iscience_revision_10.9/figure/feature_importance.xlsx")
num<-which(dat[,2]=="ARAF")
num_2<-which(dat[,1]=="GBM")
newdat<-dat[intersect(num,num_2),]
newdat<-dat[num_2,]
ggplot(newdat,aes(x=kinase,y=importance,fill=gene)) + geom_bar(position="dodge",stat="identity",width=0.9)+
  scale_fill_simpsons()+
  labs(NULL)+
  theme(title=element_text(size=8),axis.text = element_text(size = 8),
        legend.text = element_text(size = 7),legend.title = element_text(size = 7),axis.text.x=element_text(angle=40,hjust=1,vjust=1))+
  theme_classic()
###plot classfication results

dat<-read.table(paste("C:/users/90410/desktop/kinase/iscience_revision_10.9/result/",cancer_type[i],"/XGboost_result.csv",sep=""),
                sep=",")
###compare F1 score; Recall, Precision  5*6.5
library(ggplot2)
ROC<-as.data.frame(matrix(NA,0,3))
colnames(ROC)<-c("fpr","tpr","model")
model<-c("XGboost","RF","SVM","logisticRegression")
for(i in 1:4){
  total<-read.csv(paste("C:/Users/90410/Desktop/kinase/iscience_revision_10.9/figure/figure4/UCEC_",model[i],".csv",sep=""),header=TRUE)
  total<-total[,-1]
  total<-total[-1,]
  total$model<-model[i]
  new<-total[,c(1,2,5)]
  new<-na.omit(total)
  ROC<-rbind(ROC,new)}


AUC<-as.data.frame(matrix(NA,0,3))
colnames(AUC)<-c("recall","precision","model")
for(i in 1:4){
  total<-read.csv(paste("C:/Users/90410/Desktop/kinase/iscience_revision_10.9/figure/figure4/UCEC_",model[i],".csv",sep=""),header=TRUE)
  total<-total[,-1]
  total<-total[-1,]
  total$model<-model[i]
  new<-total[,c(3,4,5)]
  new<-na.omit(new)
  add<-as.data.frame(matrix(NA,1,3))
  colnames(add)<-c("recall","precision","model")
  add[1,]<-c(1,0,model[i])
  new<-rbind(add,new)
  AUC<-rbind(AUC,new)}


ggplot(ROC, mapping = aes(x = fpr, y = tpr, color = model)) + geom_line() +scale_fill_simpsons()+
  xlab('False Positive rate')+ylab('True Positive rate')+
  theme(title=element_text(size=8),axis.text = element_text(size = 8),
        legend.text = element_text(size = 7),legend.title = element_text(size = 7),axis.text.x=element_text(angle=40,hjust=1,vjust=1))+
  theme_classic()

AUC$recall<-as.numeric(AUC$recall)
AUC$precision<-as.numeric(AUC$precision)
ggplot(AUC, mapping = aes(x = recall, y = precision, color = model)) + geom_line() +scale_fill_simpsons()+
  xlab('Recall')+ylab('Precision')+
  theme(title=element_text(size=8),axis.text = element_text(size = 8),
        legend.text = element_text(size = 7),legend.title = element_text(size = 7),axis.text.x=element_text(angle=40,hjust=1,vjust=1))+
  theme_classic()

newdat<-read_xlsx("C:/Users/90410/Desktop/kinase/iscience_revision_10.9/f1_score.xlsx")

newdat$type<-factor(newdat$type, levels=c('top10', '31-40', '91-100', 'Low10'))
ggplot(newdat,aes(x=type,y=Value,fill=Metric)) + geom_bar(position="dodge",stat="identity",width=0.9)+
  scale_fill_simpsons()+
  labs(NULL)+
  theme(title=element_text(size=8),axis.text = element_text(size = 8),
        legend.text = element_text(size = 7),legend.title = element_text(size = 7),axis.text.x=element_text(angle=40,hjust=1,vjust=1))+
  theme_classic()

###plot scRNA-seq Breast cancer kinase results
##"BC","GBM","HCC","LSCC","UCEC"
library(ggsci)
library(Seurat)
###combine all five samples
dat1<-read.table("C:/Users/90410/Desktop/KBPRNA/bc_scRNA/GSE180286_RAW/GSM5457199/GSM5457199_A2019-1.expression_matrix.txt",
                sep="\t",header=TRUE,row.names=1)
dat2<-read.table("C:/Users/90410/Desktop/KBPRNA/bc_scRNA/GSE180286_RAW/GSM5457202/GSM5457202_B2019-1.expression_matrix.txt",
                 sep="\t",header=TRUE,row.names=1)
dat3<-read.table("C:/Users/90410/Desktop/KBPRNA/bc_scRNA/GSE180286_RAW/GSM5457205/GSM5457205_C2020-1.expression_matrix.txt",
                 sep="\t",header=TRUE,row.names=1)
dat4<-read.table("C:/Users/90410/Desktop/KBPRNA/bc_scRNA/GSE180286_RAW/GSM5457208/GSM5457208_D2020-1.expression_matrix.txt",
                 sep="\t",header=TRUE,row.names=1)
dat5<-read.table("C:/Users/90410/Desktop/KBPRNA/bc_scRNA/GSE180286_RAW/GSM5457211/GSM5457211_E2020-1.expression_matrix.txt",
                 sep="\t",header=TRUE,row.names=1)

###combine kinase_activity with dat
pbmc1 <- CreateSeuratObject(counts = dat1,
                           project = "kinase",
                           min.cells = 3,
                           min.features = 200)

pbmc2 <- CreateSeuratObject(counts = dat2,
                            project = "kinase",
                            min.cells = 3,
                            min.features = 200)

pbmc3 <- CreateSeuratObject(counts = dat3,
                            project = "kinase",
                            min.cells = 3,
                            min.features = 200)

pbmc4 <- CreateSeuratObject(counts = dat4,
                            project = "kinase",
                            min.cells = 3,
                            min.features = 200)

pbmc5 <- CreateSeuratObject(counts = dat5,
                            project = "kinase",
                            min.cells = 3,
                            min.features = 200)
pbmc <- merge(pbmc1, y = c(pbmc2, pbmc3,pbmc4,pbmc5), add.cell.ids = c("1","2","3","4","5"), project = "kinase")

pbmc<-pbmc1

###线粒体和红细胞基因比例
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") # 人类血液常见红细胞基因
HB_m <- match(HB.genes_total,rownames(pbmc@assays$RNA))
HB.genes <- rownames(pbmc@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
pbmc[["percent.HB"]]<-PercentageFeatureSet(pbmc,features=HB.genes)

head(pbmc@meta.data)[,c(2,3,4,5)]
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)
###查看这几个指标间关系
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")+scale_color_npg()
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+scale_color_npg()
plot3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.HB")+scale_color_npg()
CombinePlots(plots = list(plot1, plot2,plot3),legend="none")

###normalization
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)  # 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
##查找高变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2),legend="bottom")

###standardization
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

DimPlot(pbmc, reduction = "pca")

###cluster
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

table(pbmc@active.ident)

pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)
plot1<-DimPlot(pbmc, reduction = "umap",label = TRUE)+scale_color_npg()
plot2<-DimPlot(pbmc, reduction = "tsne",label = TRUE)+scale_color_npg()
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1<-VlnPlot(pbmc, features = c("PTPRC","EPCAM"))+scale_color_npg()
plot2<- VlnPlot(pbmc, features = c("MAPK1"),ncol=1, same.y.lims=T,slot = "counts", log = TRUE)+scale_color_npg()
CombinePlots(plots = list(plot1, plot2))
Idents(pbmc,cells=1:3180)<-"LUNG_N08"
Idents(pbmc,cells=3181:7145)<-"LUNG_T08"
plot1<-FeaturePlot(pbmc, features = c("ROCK1","ROCK2",
                                      "UBE2C","PTPRC"),reduction="tsne",min.cutoff = 0, max.cutoff = 4)
FeaturePlot(pbmc, features = c("ROCK1"),min.cutoff = 0,max.cutoff = 10)
Idents(pbmc)
###compute high variable kinase activities

write.csv(dat1,"C:/users/90410/desktop/KBPRNA/GSM5457199.csv")
write.csv(dat2,"C:/users/90410/desktop/KBPRNA/GSM5457202.csv")
write.csv(dat3,"C:/users/90410/desktop/KBPRNA/GSM5457205.csv")
write.csv(dat4,"C:/users/90410/desktop/KBPRNA/GSM5457208.csv")
write.csv(dat5,"C:/users/90410/desktop/KBPRNA/GSM5457211.csv")
##read in Breast cancer kinase activity

kinase_activity_BC1<-read.csv("C:/users/90410/desktop/KBPRNA/predict_result_199.csv",row.names=1)
kinase_activity_BC2<-read.csv("C:/users/90410/desktop/KBPRNA/predict_result_202.csv",row.names=1)
kinase_activity_BC3<-read.csv("C:/users/90410/desktop/KBPRNA/predict_result_205.csv",row.names=1)
kinase_activity_BC4<-read.csv("C:/users/90410/desktop/KBPRNA/predict_result_208.csv",row.names=1)
kinase_activity_BC5<-read.csv("C:/users/90410/desktop/KBPRNA/predict_result_211.csv",row.names=1)
colnames(kinase_activity_BC5)<-paste("5_",colnames(kinase_activity_BC5),sep="")

total<-cbind(kinase_activity_BC1,kinase_activity_BC2,kinase_activity_BC3,kinase_activity_BC4,kinase_activity_BC5)
newkinase_epi<-as.data.frame(matrix(NA,19,0))
rownames(newkinase_epi)<-c(0:18)
##EPCAM 2,8,9,11
##PTPRC 3,4,5,7,10,12
epithelial<-c(2,8,9,11)
immune<-c(3,4,5,7,10,12)

for(i in 1:length(epithelial)){
  epi_num<-which(Idents(pbmc)==epithelial[i])
  newkinase_epi<-cbind(newkinase_epi,total[,names(Idents(pbmc)[epi_num])])
}

newkinase_imu<-as.data.frame(matrix(NA,19,0))
rownames(newkinase_imu)<-c(0:18)

for(i in 1:length(immune)){
  imu_num<-which(Idents(pbmc)==immune[i])
  newkinase_imu<-cbind(newkinase_imu,total[,names(Idents(pbmc)[imu_num])])
}

t.test(as.numeric(newkinase_epi[1,]),as.numeric(newkinase_imu[15,]))
newepi_imu<-cbind(newkinase_epi,newkinase_imu)


dat<-as.data.frame(matrix(NA,3143,3))
colnames(dat)<-c("sample",'ARAF','ERN1','MAPK1','PLK1','CDK2',
                 'MAPK3','CDK1','TAOK3','CSNK1E','PRKCA',
                 'TNK2','ABL1','CDK5','PAK3','PAK1',
                 'MAPK13','CDK12','STK38','UHMK1')
dat[1:3143,2:20]<-t(10^newepi_imu[1:19,1:3143])

dat[1:1559,1]<-"epithelial"
dat[1560:3143,1]<-"immune"

dat$PAK1<-as.numeric(dat$PAK1)
t.test(dat$PAK1[1:1559],dat$PAK1[1560:3143])
ggplot(newdat,aes(x=sample,y=PAK1))+
  geom_violin(aes(color=sample),trim=FALSE)+
  geom_boxplot(aes(color=sample),width=0.05)+
  scale_color_brewer(palette="Set2")+theme_classic()

t.test(dat$PAK1[1:1559],dat$PAK1[1560:3143])
ten_N<-quantile(dat[1:1559,15],0.1)
ninty_N<-quantile(dat[1:1559,15],0.9)

epi<-dat[1:1559,]
num_1<-which(epi[,15]>=ten_N)
num_2<-which(epi[,15]<=ninty_N)
newepi<-epi[intersect(num_1,num_2),]

ten_T<-quantile(dat[1559:3143,15],0.1)
ninty_T<-quantile(dat[1559:3143,15],0.9)
imu<-dat[1559:3143,]

num_1<-which(imu[,15]>=ten_T)
num_2<-which(imu[,15]<=ninty_T)
newimu<-imu[intersect(num_1,num_2),]

newdat<-rbind(newepi,newimu)
###
new.cluster.ids <- c(rep("others",2), "Epithelial", 
                     rep("Immune",3), "others", "Immune",rep("Epithelial",2),
                     "Immune", "others", "Immune","others")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + scale_color_npg()

kinase_activity<-10^kinase_act
rownames(kinase_activity)<-c('ARAF','ERN1','MAPK1','PLK1','CDK2',
                             'MAPK3','CDK1','TAOK3','CSNK1E','PRKCA',
                             'TNK2','ABL1','CDK5','PAK3','PAK1',
                             'MAPK13','CDK12','STK38','UHMK1')

rownames(kinase_activity)<-paste(rownames(kinase_activity),".kinase",sep="")

dat5<-rbind(dat5,kinase_activity)



###deal with LSCC dataset
library(data.table)

LSCC<-fread("C:/Users/90410/Desktop/kinase/iscience_revision_10.9/LSCC_normal_tumor.txt",sep=",")
LSCC<-as.data.frame(LSCC)
rownames(LSCC)<-LSCC[,1]
LSCC<-LSCC[,-1]
kinase_act<-read.table("C:/Users/90410/Desktop/kinase/iscience_revision_10.9/predict_result.csv",sep=",",row.names = 1)
rownames(kinase_act)<-c("ROCK1","ROCK2")
colnames(kinase_act)<-colnames(LSCC)
kinase_act<-10^kinase_act
num<-which(gsub(colnames(kinase_act)))
newLSCC<-rbind(LSCC,kinase_act)
n<-0
for(i in 1:ncol(newLSCC)){
  if(grep("LUNG_N08",colnames(newLSCC)[i]))n<-n+1
}
##violin plot
dat<-as.data.frame(matrix(NA,7145,3))
colnames(dat)<-c("sample","kin_act_ROCK1","kin_act_ROCK2")
dat[1:7145,2]<-t(kinase_act[1,1:7145])
dat[1:7145,3]<-t(kinase_act[2,1:7145])
dat[1:3380,1]<-"LUNG_N08"
dat[3381:7145,1]<-"LUNG_T08"

ten_N<-quantile(dat[1:3380,3],0.1)
ninty_N<-quantile(dat[1:3380,3],0.9)

normal<-dat[1:3380,]
num_1<-which(normal[,3]>=ten_N)
num_2<-which(normal[,3]<=ninty_N)
newnormal<-normal[intersect(num_1,num_2),]

ten_T<-quantile(dat[3381:7145,3],0.1)
ninty_T<-quantile(dat[3381:7145,3],0.9)
tumor<-dat[3381:7145,]

num_1<-which(tumor[,3]>=ten_T)
num_2<-which(tumor[,3]<=ninty_T)
newtumor<-tumor[intersect(num_1,num_2),]

newdat<-rbind(newnormal,newtumor)
t.test(as.numeric(newtumor[,2]),as.numeric(newnormal[,2]))
ggplot(newdat,aes(x=sample,y=kin_act_ROCK2))+
  geom_boxplot(aes(color=sample),width=0.05)+
  geom_violin(aes(color=sample),trim=FALSE)+
  scale_color_brewer(palette="Set1")+theme_classic()

LSCC_seurat <- CreateSeuratObject(counts = LSCC,
                            project = "kinase",
                            min.cells = 3,
                            min.features = 200)

##assess UBE2C relationship with ROCK1
which(rownames(LSCC)=="ROCK2")
which(rownames(LSCC)=="ROCK1")
cor.test(as.numeric(LSCC[17973,]),as.numeric(LSCC[27975,]))
##(0.017081494, 'SOCS2'), (0.017979022, 'AURKB'), (0.024346163, 'CHEK2'), (0.21687329, 'CCNB1'), (0.5760872, 'UBE2C')
##(0.012541928, 'RBKS'), (0.0150713, 'SOCS2'), (0.04308557, 'PLK1'), (0.14542842, 'CHEK2'), (0.62600523, 'UBE2C')

library(readxl)
importance<-read_xlsx("C:/users/90410/desktop/kinase/iscience_revision_10.9/figure/figure5/ROCK1_ROCK2_importance.xlsx")
ggplot(importance,aes(x=Kinase,y=Value,fill=Gene)) + geom_bar(position="dodge",stat="identity",width=0.9)+
  scale_fill_npg()+
  labs(NULL)+
  theme(title=element_text(size=8),axis.text = element_text(size = 8),
        legend.text = element_text(size = 7),legend.title = element_text(size = 7),axis.text.x=element_text(angle=40,hjust=1,vjust=1))+
  theme_classic()