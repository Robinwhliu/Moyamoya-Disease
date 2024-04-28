
rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# options(BioC_mirror="https://anaconda.org/bioconda/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

getwd()
setwd('C:/Users/LWH/Downloads/bioinfo/MDD') # your path
library(dplyr)
library(plyr)
library(tidyr)
library(sva)
library(limma)
library(tidyverse)
library(magrittr)
library(set)
library(biomaRt)
library(stringr)
source("flat_violin.R")
differ = function(group_list, exp) {
  design <- model.matrix(~0+factor(group_list))#制作分组矩阵
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exp)
  
  cont.matrix=makeContrasts(contrasts=c('Case-Control'), levels = design)#制作差异比较矩阵
  
  fit <- lmFit(exp, design)#构建线性拟合模型
  
  fit2=contrasts.fit(fit, cont.matrix)#根据线性模型和比较矩阵进行差值运算
  
  fit2=eBayes(fit2)#贝叶斯检验
  
  tempOutput = topTable(fit2, coef='Case-Control', n=Inf)#生成结果报告
  
  deg = na.omit(tempOutput)
  return(deg)
}

get_log<-function(ad){
  expset<-ad
  qx <- as.numeric(quantile(expset, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) { expset[which(expset <= 0)] <- NaN
  expset<- log2(expset)
  print("log2 transform finished")}else{print("log2 transform not needed")}
  return(expset)
}
func1 = function(l, x) {
  h = c()
  for (i in c(1:l$lengths[x])) {
    p = paste0(l$values[x], i, collapse = "_")
    h = c(h,p)
  }
  return(h)
}
get_ad<-function(ad,fd){
  #ad=ad1
  index<-intersect(ad$ID_REF,fd$ID)
  rownames(ad)<-ad$ID_REF
  rownames(fd)<-fd$ID
  ad<-ad[index,-1];fd=fd[index,]
  ad%<>%as.data.frame(.) %>% mutate(max=apply(.,1,max)) %>% 
    mutate(symbol=fd[rownames(.),]$GENE_SYMBOL) %>% 
    .[order(.$symbol,.$max,decreasing = T),] %>% 
    .[!duplicated(.$symbol)&!is.na(.$symbol),]
  rownames(ad)<-ad$symbol
  return(ad[,-c(ncol(ad)-1,ncol(ad))])
}

featuredata<-read.delim("01_getdata/GPL16699-15607.txt",header = T,
                        sep = "\t",skip = 19)  # Please download GPL16699-15607 file from GEO database
colnames(featuredata)
fd<-featuredata[featuredata$GENE_SYMBOL!="",
                c("ID","GENE_SYMBOL","GENE_NAME","UNIGENE_ID","ENSEMBL_ID",
                  "DESCRIPTION","GO_ID")]

GSE<-c("GSE141022","GSE141024")
get_pd<-function(GSE){
  phenoData <- read.table(paste0("01_getdata/",GSE,"_series_matrix.txt"),
                          header = T, sep = "\t",skip = 26)
  l = rle(phenoData[duplicated(phenoData[,1]),]$X.Sample_title)
  phenoData[duplicated(phenoData[,1]),]$X.Sample_title = unlist(lapply(c(1:length(l$values)), 
                                                                       function(x) func1(l, x)))
  rownames(phenoData)=phenoData$X.Sample_title
  rownames(phenoData) <- phenoData$X.Sample_title
  rownames(phenoData)
  pd<-phenoData[,-1] %>% t() %>% as.data.frame() %>% 
    .[,c("!Sample_geo_accession","!Sample_characteristics_ch1",
         "!Sample_characteristics_ch11","!Sample_characteristics_ch12",
         "!Sample_characteristics_ch13","!Sample_characteristics_ch14")]
  colnames(pd)<-c("ID","Condition","Number","Tissue","Gender","Age")
  pd<-apply(pd, 2, function(x) gsub(".*: ","",x)) %>% as.data.frame()
  pd$Age %<>% gsub("y","",.) %>% as.numeric()
  pd$Group<-lapply(pd$Condition,
                   function(x)ifelse(grepl("Moyamoya",x),"Case","Control")) %>% 
    unlist()
  return(pd)
}
pd1<-get_pd(GSE[1])
pd2<-get_pd(GSE[2])

ad1<-read.table("01_getdata/GSE141022_array.txt",header = T,sep = "\t")
ad2<-read.table("01_getdata/GSE141024_array.txt",header = T,sep = "\t")

ad1<-get_ad(ad1,fd)
ad2<-get_ad(ad2,fd)

get_batch<-function(ad,pd){
  batch=pd$Gender %>% as.factor()
  # batch2=pd$ %>% as.factor()
  design=data.frame(row.names=rownames(pd),
                    lapply(pd$Group,function(x)ifelse(x=="Control",0,1)) %>% 
                      unlist())
  ad<-removeBatchEffect(ad,batch=batch,
                        covariates=pd$Age%>%as.numeric(),
                        design=design)%>%as.data.frame() 
  return(ad)
}
ad_GSE141022<-get_batch(ad1,pd1)
ad_GSE141024<-get_batch(ad2,pd2)


deg1<-differ(pd1$Group,ad_GSE141022)
deg2<-differ(pd2$Group,ad_GSE141024)

MALAT1 <- c("TNRC6C","HLA-A","RPL12","PODXL",
            "EPAS1","FAM118A","SOCS3","GADD45B")

get_distributions<-function(exp){
  library(ggsci)
  distributions <- 
    ggplot(data = exp, 
           aes(x=group, y=Expression,fill=Gene)) +
    geom_flat_violin(position=position_nudge(x=0.2,y=0),alpha=0.8,trim=F) +
    geom_point(aes(y=Expression, color=Gene), 
               position = position_jitter(width=0.15),size=1.8,alpha=0.5) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    stat_compare_means(aes(group=group),label="p.signif",label.x=1.5,
                       size=5,hjust=0.3,vjust=-.5,method="t.test")+
    labs(y = "Expression Level", x = NULL) +
    guides(fill='none',color='none')+scale_color_ucscgb()+scale_fill_ucscgb()+
    # scale_y_continuous(limits=c(0,5)) +
    # scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
    # scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
    # facet_grid(cols=vars(Group),rows=vars(Disease),scales="free",space="free")+
    facet_wrap(vars(Gene),ncol=3)+
    coord_flip()+theme_niwot()+
    scale_x_discrete(limits=c("Case","Control"),labels=c("MMD","Control"))+
    theme(strip.text=element_text(size=12,face="bold.italic"),
          strip.background=element_rect(fill="white"))
  return(distributions)
}

GSE141022<-ad_GSE141022[c(MALAT1,"MALAT1"),] %>% t() %>% as.data.frame()
GSE141022t<-GSE141022 %>% mutate(group=pd1$Group)
exp1 <- melt(GSE141022, id.vars = c("MALAT1"),
             variable.name="Gene",value.name="Expression") %>%
  mutate(face=paste0("MALAT1 vs ",Gene))
exp1t <- melt(GSE141022t, id.vars = c("group"),
              variable.name="Gene",value.name="Expression")

GSE141024<-ad_GSE141024[c(MALAT1,"MALAT1"),] %>% t() %>% as.data.frame()
GSE141024t<-GSE141024 %>% mutate(group=pd2$Group)
exp2 <- melt(GSE141024, id.vars = c("MALAT1"),
             variable.name="Gene",value.name="Expression") %>%
  mutate(face=paste0("MALAT1 vs ",Gene))
exp2t <- melt(GSE141024t, id.vars = c("group"),
              variable.name="Gene",value.name="Expression")

distributions1<-get_distributions(exp1t)
distributions2<-get_distributions(exp2t)
