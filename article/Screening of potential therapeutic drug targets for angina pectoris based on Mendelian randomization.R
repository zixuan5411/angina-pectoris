
library(TwoSampleMR)
library(ggplot2)
library(readxl)
library(dplyr)

#exposure<-read_excel("Zheng_Supplementary_Tables.xlsx",sheet = 2) %>%as.data.frame() #intrument 2113 protein 1699

#Tier1<-subset(exposure,Type=="Tier1") #intrument1064 protein 955

#length(unique(Tier1$Exposure))
#dup_gene<-Tier1[duplicated(Tier1$Exposure) ,]$Exposure
#exposure_dup<-Tier1[which(Tier1$Exposure %in% dup_gene) ,]
#exposure_nodup<-Tier1[-which(Tier1$Exposure %in% dup_gene) ,]
#table(exposure_nodup$Exposure %in% exposure_dup$Exposure)

#Tier2<-subset(exposure,Type=="Tier2") #intrument62 protein 58
#Tier3<-subset(exposure,Type=="Tier3") #intrument987 protein 890
#Tier12<-subset(exposure,Type %in%c("Tier1","Tier2")) #intrument 1126 protein 890 #protein 1002
exposure_data<-read.csv("filted_protein_exposure.csv",header = T,stringsAsFactors = FALSE) 
#exposure_data<-read_excel("Yi-Jun Ge-Supplementary.xlsx",sheet = 6) %>% as.data.frame()
exposure<-format_data(exposure_data,type="exposure",header=T,snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele"
                      ,other_allele_col="other_allele",pval_col="pval",eaf_col="eaf",chr_col="chr",pos_col="pos",
                      samplesize_col="sample_size",phenotype_col = "Phenotype" )


exposure<-clump_data(exposure,clump_r2=0.001,clump_kb=10000)
exposure<-subset(exposure,pval.exposure<0.00000005|pval.exposure==0.00000005)
exposure<-subset(exposure,!((chr.exposure=="chr6") & (pos.exposure>26000000 &pos.exposure<34000000)))

length(unique(exposure$SNP))==nrow(exposure) 


length(unique(exposure$exposure)) 
table(duplicated(exposure$SNP))
table(duplicated(exposure$exposure))
dup_gene<-exposure[duplicated(exposure$exposure) ,]$exposure

exposure_dup<-exposure[which(exposure$exposure %in% dup_gene),]
exposure_nodup<-exposure[-which(exposure$exposure %in% dup_gene) ,]

table(exposure_nodup$exposure %in% exposure_dup$exposure)
#exposure<-format_data(exposure_data,type="exposure",header=T,snp_col="QTL",beta_col="Beta",se_col="SE",effect_allele_col="EA"
#                         ,other_allele_col="OA",pval_col="P",eaf_col="EAF",chr_col="Chr",pos_col="Pos",
#                 samplesize_col="N",phenotype_col = "Protein" )


outcome<-extract_outcome_data(snps = exposure$SNP,outcomes="ebi-a-GCST90018793")


harmonised_nodup<- harmonise_data(exposure_nodup, outcome)
res_single_nodup<- mr_singlesnp(dat=harmonised_nodup,single_method ="mr_wald_ratio")

res_single_nodup=res_single_nodup[!is.na(res_single_nodup$p),]
res_single_nodup$Method="mr_wald_ratio"
res_single_nodup <- generate_odds_ratios(res_single_nodup)
write.csv(res_single_nodup,"ebi-a-GCST90018793_singlesnp_nodup.csv",row.names = F)

dim(res_single_nodup[res_single_nodup$p< 0.05/length(unique(exposure$exposure)),])



harmonised_dup<- harmonise_data(exposure_dup, outcome)
#res_single_dup=data.frame()
#for (i in unique(harmonised_dup$exposure)){
# a=mr_singlesnp(dat=harmonised_dup[harmonised_dup$exposure==i,],single_method ="mr_ivw")
#res_single_dup=rbind(res_single_dup,a)
#}

#res_single_dup<- mr_singlesnp(dat=harmonised_dup,single_method ="mr_ivw")
res_single_dup=mr(dat=harmonised_dup,method_list="mr_ivw")

#res_single_dup=res_single_dup[!is.na(res_single_dup$p),]
res_single_dup <- generate_odds_ratios(res_single_dup)
res_single_dup$Method<-"mr_ivw"
write.csv(res_single_dup,"ebi-a-GCST90018793_singlesnp_dup.csv",row.names = F)

dim(res_single_dup[res_single_dup$p< 0.05/length(unique(exposure$exposure)),])

#1.Heterogeneity statistics
heterogeneity=mr_heterogeneity(harmonise_data(exposure, outcome))
write.csv(heterogeneity,"heterogeneity.csv",row.names = F) 

mr_heterogeneity(harmonised_nodup)

#2.Horizontal pleiotropy 
mr_pleiotropy_test(harmonise_data(exposure, outcome)) #NA
mr_pleiotropy_test(harmonised_dup)

#3.Leave-one-out analysis
lot<- mr_leaveoneout(harmonise_data(exposure, outcome))


p1 <- mr_scatter_plot(res_single_dup, harmonised_dup)
p2<- mr_scatter_plot(res_single_nodup, harmonised_nodup)

##Forest plot
dim(res_single_nodup[res_single_nodup$p< 0.05/length(unique(exposure$exposure)),])
p3<- mr_forest_plot(mr_singlesnp(dat=harmonised_dup,single_method ="mr_ivw"))
p3[[1]]
p4<- mr_forest_plot(res_single_nodup)

sig_protein<-res_single_nodup[res_single_nodup$p<0.05/length(unique(exposure$exposure)),]
sig_protein<-sig_protein[,c(1,5,6,13,14,15,9,8)]
colnames(sig_protein)<-c("Protein","samplesize","SNP","OR","low-CI95%","high-CI95%","P","SE")
library(forestploter)
sig_protein$` ` <- paste(rep(" ", 7), collapse = " ")

tm <- forest_theme(base_size = 12,  
                   ci_pch = 15, 
                   ci_col = "black",   
                   ci_fill = "blue",     
                   ci_alpha = 1,       
                   ci_lty = 1,           
                   ci_lwd = 2,        
                   ci_Theight = 0.2,
                   vertline_lwd = 2,            
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   summary_fill = "yellow",      
                   summary_col = "#4575b4",
)
pdf("sig_protein_forest.pdf",width=10)
forest(sig_protein[,c(1,2,3,4,9,5,6,7,8)],
       est = sig_protein$OR,       
       lower = sig_protein$`low-CI95%`,     
       upper = sig_protein$`high-CI95%`,     
       sizes = sig_protein$SE,
       ci_column = 5,   
       ref_line = 1,
       arrow_lab = c("Lower", "Higher"),
       xlim = c(0.5,2 ),
       ticks_at = c(0.5, 1, 2),
       footnote = "",
       theme = tm)
dev.off()
sig_protein$Association<-ifelse(sig_protein$OR>1,"Positive","Negative")

#Leave-one-out plot
p5 <- mr_leaveoneout_plot(lot)

##
p6<- mr_funnel_plot(res_single_dup)
p7<- mr_funnel_plot(res_single_nodup)

res_single_nodup2<-res_single_nodup[,c(1,9,13)]
res_single_nodup2$Association<-ifelse(
  res_single_nodup2$or>1&res_single_nodup2$p< 0.05/length(unique(exposure$exposure)),"Positive",
  ifelse(res_single_nodup2$or<1&res_single_nodup2$p< 0.05/length(unique(exposure$exposure)),"Negative","Not sig"))


res_single_nodup2$Protein<-ifelse(res_single_nodup2$Association=="Not sig","",res_single_nodup2$exposure)
library(ggplot2)
library(ggrepel)
p1<-ggplot(data=res_single_nodup2, aes(x = or, y = -log10(p), color = Association))+
  geom_point(alpha=0.4, size=3.5)+
  geom_vline(xintercept=1,lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept= -log10(0.05/length(unique(exposure$exposure))),lty=4,col="black",lwd=0.8)+
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5),legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15),axis.title = element_text(size=15),
        panel.background = element_rect(fill = "white",colour = "black",linewidth = 0.8),
        axis.text = element_text(size=15))+
  labs(x="OR",y="-log10(Pvalue)")+
  scale_color_manual(values=c("#BA55D3", "gray","#A52A2A"))+
  geom_label_repel(data = res_single_nodup2, aes(x = or, y = -log10(p), 
                                                 label = Protein),max.overlaps =10,size = 5, box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"), segment.color = "black", 
                   show.legend = FALSE)

ggsave("sig_protein_volca.pdf",p1,width = 8,height = 10,dpi=500)

library(dplyr)  
library(coloc)

GWASdata=outcome[,c("SNP","chr","pos",'effect_allele.outcome',
                    'other_allele.outcome',"eaf.outcome","beta.outcome","se.outcome",
                    "pval.outcome","samplesize.outcome")]
GWASdata$ncase.outcome=30025
GWASdata$ncontrol.outcome=440906

GWASdata$varbeta <- GWASdata$se.outcome^2
GWASdata$s <- GWASdata$ncase.outcome/GWASdata$samplesize.outcome
GWASdata$z = GWASdata$beta.outcome/GWASdata$se.outcome

QTLdata=exposure[,-c(12,13,14)]

QTLdata$varbeta <- QTLdata$se.exposure^2
QTLdata$MAF <- ifelse(QTLdata$eaf.exposure<0.5,QTLdata$eaf.exposure,1-QTLdata$eaf.exposure)
QTLdata$z = QTLdata$beta.exposure/QTLdata$se.exposure
#
library(tidyverse)
library(patchwork)
library(locuscomparer)
coloc_res=data.frame()
for (i in sig_protein$SNP){
  GWASdata2=GWASdata[GWASdata$SNP==i,]
  QTLdata2=QTLdata[QTLdata$SNP==i,]
  coloc_result <- coloc.abf(dataset1=list(pvalues=GWASdata2$pval.outcome, snp=GWASdata2$SNP, type="cc", s=GWASdata2$s[1], N=GWASdata2$samplesize.outcome[1]), 
                            dataset2=list(pvalues=QTLdata2$pval.exposure,snp=QTLdata2$SNP, type="quant", N=QTLdata2$samplesize[1]), MAF=QTLdata2$MAF)
  
  t=data.frame(t(data.frame(coloc_result[[1]])))
  rownames(t)<-i
  coloc_res=rbind(coloc_res,t)
  chr=exposure[exposure$SNP==i,]$chr.exposure
  x=exposure[exposure$chr.exposure %in% chr,][,c("chr.exposure","SNP","pval.exposure")]
  y=outcome[outcome$chr %in% chr,][,c("chr","SNP","pval.outcome")]
  
  colnames(x)=c("chr","rsid","pval")
  colnames(y)=c("chr","rsid","pval")
  pdf(paste0(i,"_locuscompare.pdf"))
  p1=locuscompare(in_fn1 = x, in_fn2 = y, title1 = paste0(i,' pQTL'), 
                  title2 = 'Angina Outcome GWAS',snp=i,combine=T)
  print(p1)
  dev.off()
}

write.csv(coloc_res,'coloc.abf_result.csv',row.names=T) 

#x=GWASdata[GWASdata$SNP %in%sig_protein$SNP,]
#y=QTLdata[QTLdata$SNP %in% sig_protein$SNP,]

#z<- coloc.abf(dataset1=list(pvalues=x$pval.outcome, snp=x$SNP, type="cc", s=x$s[1], N=x$samplesize.outcome[1]), 
#                           dataset2=list(pvalues=y$pval.exposure,snp=y$SNP, type="quant", N=y$samplesize[1]), MAF=y$MAF)

##coloc.susie
sig_protein$SNP


c("rs1800470","rs111577916","rs415895","rs191448950","rs10950398","rs964184","rs55730499")
#LD=matrix(c(1,0.3233,0.0019,0.0831,0.0175,0.0049,0.0067,
#           0.3233,1,0.1981,1,0.3716,0.6539,0.014,
#          0.0019,0.1981,1,0.0878,0.0482,0.0435,0.3274,
#         0.0831,1,0.0878,1,0.0809,0.4376,1,
#        0.0175,0.3716,0.0482,0.0809,1,0.0528,0.2753,
##       0.0049,0.6539,0.0435,0.4376,0.0528,1,0.2017,
#     0.0067,0.014,0.3274,1,0.2753,0.2017,1),nrow=7)
rownames(LD)=sig_protein$SNP
colnames(LD)=sig_protein$SNP


GWASdata$MAF<-ifelse(GWASdata$eaf.outcome<0.5,GWASdata$eaf.outcome,1-GWASdata$eaf.outcome)
#

#QTLdata2=QTLdata[QTLdata$SNP %in%sig_protein$SNP,]
#QTLdata2=QTLdata2[c(5,6,4,3,7,1,2),]

#GWASdata2=GWASdata[GWASdata$SNP %in%sig_protein$SNP,]
#GWASdata2=GWASdata2[c(6,7,1,3,5,2,4),]
#rownames(GWASdata2)=1:nrow(GWASdata2)
#rownames(QTLdata2)=1:nrow(QTLdata2)

#QTLdata2<-runsusie(list(pvalues=QTLdata2$pval.exposure,snp=QTLdata2$SNP, type="quant", N=QTLdata2$samplesize.exposure[1], MAF=QTLdata2$MAF,
#                         LD=LD,beta=QTLdata2$beta.exposure,varbeta=QTLdata2$varbeta,position=QTLdata2$pos.exposure))

#GWASdata2<-runsusie(list(pvalues=GWASdata2$pval.outcome, snp=GWASdata2$SNP, type="cc", s=GWASdata2$s[1], N=GWASdata2$samplesize.outcome[1],
#                          LD=LD,beta=GWASdata2$beta.outcome,varbeta=GWASdata2$varbeta,MAF=GWASdata2$MAF,position=GWASdata2$pos))
#susie_res=coloc.susie(QTLdata2,GWASdata2)


#write.csv(susie_res,'coloc.susie_result.csv',row.names=T) 


direct_nodup=directionality_test(harmonised_nodup[harmonised_nodup$exposure %in%sig_protein$Protein,])
write.table(direct_nodup,"Steiger_directionality_test_sig_protein.csv",sep=",",row.names = F)
direct_dup=directionality_test(harmonised_dup)

#Steiger
steiger_nodup=steiger_filtering(harmonised_nodup[harmonised_nodup$exposure %in%sig_protein$Protein,])
steiger_dup=steiger_filtering(harmonised_dup)
write.table(steiger_nodup,"Steiger_filter_sig_protein.csv",sep=",",row.names = F)
save.image("mr1.rdata")




#sig_protein$Protein "TIRAP","SWAP70","PCSK9","TMEM106B","APOA5","LPA"
#sig_protein$SNP "rs111577916","rs415895","rs191448950","rs10950398","rs964184","rs55730499"
library(TwoSampleMR)
library(ggplot2)
library(readxl)
library(dplyr)
library(coloc)
#exposure_data<-read_excel("Maik-Pietzner-Supplementary_tables.xlsx",sheet = 4) %>%as.data.frame()
exposure_data<-read.csv("filted_protein_exposure.csv",header = T,stringsAsFactors = FALSE) 
#exposure_data<-exposure_data[,c("HGNC.symbol.protein","rsID","cis.trans","Chr","Position","EA","NEA",
#"EAF","Effect","SE","P.value")]

#colnames(exposure_data)<-c("Phenotype","SNP","cis_trans","chr","pos","effect_allele","other_allele",
# "eaf","beta","se","pval")

exposure<-format_data(exposure_data,type="exposure",header=T,snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele"
                      ,other_allele_col="other_allele",pval_col="pval",eaf_col="eaf",chr_col="chr",pos_col="pos",
                      samplesize_col="sample_size",phenotype_col = "Phenotype" )


exposure<-subset(exposure,pval.exposure<0.00000005 | pval.exposure==0.00000005 &!((chr.exposure=="chr6") & (pos.exposure>26000000 &pos.exposure<34000000)))
exposure<-clump_data(exposure, clump_kb=10000, clump_r2=0.001)

#exposure_data$sample_size=10708
#exposure_data$Phenotype<-gsub(pattern = "\'",x =exposure_data$Phenotype,replacement = "" )
#grep(pattern = "\'",x=exposure_data$Phenotype)

#exposure_data=exposure_data[exposure_data$cis_trans=="cis",]

#exposure_data$effect_allele<-gsub(pattern = "D",x=exposure_data$effect_allele,replacement = "G")
#exposure_data$other_allele<-gsub(pattern = "I",x=exposure_data$other_allele,replacement = "A")

#exposure<-format_data(exposure_data,type="exposure",header=T,snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele"
#                     ,other_allele_col="other_allele",pval_col="pval",eaf_col="eaf",chr_col="chr",pos_col="pos",
#                    samplesize_col="sample_size",phenotype_col = "Phenotype" )

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
ieugwasr::api_status()
library(ieugwasr)
#devtools::install_github("explodecomputer/plinkbinr")
library(plinkbinr)
get_plink_exe()


#clump_df<-ld_clump_local(dplyr::tibble(rsid=exposure$SNP, pval=exposure$pval.exposure, id=exposure$id.exposure), clump_kb=10000, clump_r2=0.001,clump_p=1,
#plink_bin = "D:/R-4.3.1/library/pink/plink.exe",
#bfile = "C:/Users/wangdexin/Desktop/EUR_ref/EUR") %>%as.data.frame()
#colnames(clump_df)<-c("SNP","pval.exposure","id.exposure")
#exposure<-merge(clump_df[,c(1,3)],exposure,by = c("SNP","id.exposure"))


length(unique(exposure$exposure)) 
dup_gene<-exposure[duplicated(exposure$exposure) ,]$exposure

exposure_dup<-exposure[which(exposure$exposure %in% dup_gene) ,]
exposure_nodup<-exposure[-which(exposure$exposure %in% dup_gene),]
table(exposure_nodup$exposure %in% exposure_dup$exposure)

outcome<-extract_outcome_data(snps = exposure$SNP,outcomes="finn-b-I9_ANGINA")
outcome$samplesize.outcome<-206008

harmonised_nodup<- harmonise_data(exposure_nodup, outcome)
res_single_nodup<- mr_singlesnp(dat=harmonised_nodup,single_method ="mr_wald_ratio")
res_single_nodup=res_single_nodup[!is.na(res_single_nodup$p),]
dim(res_single_nodup[res_single_nodup$p< 0.05/length(unique(exposure$exposure)),])

res_single_nodup <- generate_odds_ratios(res_single_nodup)
write.csv(res_single_nodup,"efinn-b-I9_ANGINA_singlesnp_nodup.csv",row.names = F)

harmonised_dup<- harmonise_data(exposure_dup, outcome)
res_single_dup=mr(dat=harmonised_dup,method_list="mr_ivw")
dim(res_single_dup[res_single_dup$p< 0.05/length(unique(exposure$exposure)),])
write.csv(res_single_dup,"efinn-b-I9_ANGINA_singlesnp_dup.csv",row.names = F)


heterogeneity=mr_heterogeneity(harmonise_data(exposure, outcome))
write.csv(heterogeneity,"heterogeneity_efinn-b-I9_ANGINA_validation.csv",row.names = F)
mr_heterogeneity(harmonised_nodup)

dim(res_single_nodup[res_single_nodup$p< 0.05/length(unique(exposure$exposure)),])

sig_protein<-res_single_nodup[res_single_nodup$p<0.05/length(unique(exposure$exposure)),]
sig_protein=sig_protein[-3,]
sig_protein<-sig_protein[,c(1,5,6,12,13,14,9,8)]
colnames(sig_protein)<-c("Protein","samplesize","SNP","OR","low-CI95%","high-CI95%","P","SE")
library(forestploter)
sig_protein$` ` <- paste(rep(" ", 4), collapse = " ")
tm <- forest_theme(base_size = 12,  
                   ci_pch = 15,   
                   ci_col = "black",   
                   ci_fill = "blue",    
                   ci_alpha = 1,        
                   ci_lty = 1,            
                   ci_lwd = 2,          
                   ci_Theight = 0.2,
                   vertline_lwd = 2,            
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   summary_fill = "yellow",       
                   summary_col = "#4575b4",
)

pdf("sig_protein_forest_efinn-b-I9_ANGINA_validation.pdf",width=10)
forest(sig_protein[,c(1,2,3,4,9,5,6,7,8)],
       est = sig_protein$OR,       
       lower = sig_protein$`low-CI95%`,    
       upper = sig_protein$`high-CI95%`,      
       sizes = sig_protein$SE,
       ci_column = 5,   
       ref_line = 1,
       arrow_lab = c("Lower", "Higher"),
       xlim = c(0.5,2 ),
       ticks_at = c(0.5, 1, 2),
       footnote = "",
       theme = tm)
dev.off()

sig_protein$Association<-ifelse(sig_protein$OR>1,"Positive","Negative")

res_single_nodup2<-res_single_nodup[,c(1,9,12)]
res_single_nodup2$Association<-ifelse(
  res_single_nodup2$or>1&res_single_nodup2$p< 0.05/length(unique(exposure$exposure)),"Positive",
  ifelse(res_single_nodup2$or<1&res_single_nodup2$p< 0.05/length(unique(exposure$exposure)),"Negative","Not sig"))

res_single_nodup2$Protein<-ifelse(res_single_nodup2$Association=="Not sig","",res_single_nodup2$exposure)
library(ggplot2)
library(ggrepel)
res_single_nodup2=res_single_nodup2[res_single_nodup2$Protein!="GLTPD2",]
p1<-ggplot(data=res_single_nodup2, aes(x = or, y = -log10(p), color = Association))+
  geom_point(alpha=0.4, size=3.5)+
  geom_vline(xintercept=1,lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept= -log10(0.05/length(unique(exposure$exposure))),lty=4,col="black",lwd=0.8)+
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5),legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15),axis.title = element_text(size=15),
        panel.background = element_rect(fill = "white",colour = "black",linewidth = 0.8),
        axis.text = element_text(size=15))+
  labs(x="OR",y="-log10(Pvalue)")+
  scale_color_manual(values=c("#BA55D3", "gray","#A52A2A"))+
  geom_label_repel(data = res_single_nodup2, aes(x = or, y = -log10(p), 
                                                 label = Protein),max.overlaps =10,size = 5, box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"), segment.color = "black", 
                   show.legend = FALSE)

ggsave("sig_protein_volca_efinn-b-I9_ANGINA_validation.pdf",p1,dpi=500,width=5,height = 5)

#
library(dplyr)  
library(coloc)
GWASdata=outcome[,c("SNP","chr","pos",'effect_allele.outcome',
                    'other_allele.outcome',"eaf.outcome","beta.outcome","se.outcome",
                    "pval.outcome","samplesize.outcome")]
GWASdata$ncase.outcome=18168
GWASdata$ncontrol.outcome=187840

GWASdata$varbeta <- GWASdata$se.outcome^2
GWASdata$s <- GWASdata$ncase.outcome/GWASdata$samplesize.outcome
GWASdata$z = GWASdata$beta.outcome/GWASdata$se.outcome

QTLdata=exposure[,-c(12,13,14)]

QTLdata$varbeta <- QTLdata$se.exposure^2
QTLdata$MAF <- ifelse(QTLdata$eaf.exposure<0.5,QTLdata$eaf.exposure,1-QTLdata$eaf.exposure)
QTLdata$z = QTLdata$beta.exposure/QTLdata$se.exposure
#
library(tidyverse)
library(patchwork)
library(locuscomparer)
coloc_res=data.frame()
for (i in sig_protein$SNP){
  GWASdata2=GWASdata[GWASdata$SNP==i,]
  QTLdata2=QTLdata[QTLdata$SNP==i,]
  coloc_result <- coloc.abf(dataset1=list(pvalues=GWASdata2$pval.outcome, snp=GWASdata2$SNP, type="cc", s=GWASdata2$s[1], N=GWASdata2$samplesize.outcome[1]), 
                            dataset2=list(pvalues=QTLdata2$pval.exposure,snp=QTLdata2$SNP, type="quant", N=QTLdata2$samplesize[1]), MAF=QTLdata2$MAF)
  
  t=data.frame(t(data.frame(coloc_result[[1]])))
  rownames(t)<-i
  coloc_res=rbind(coloc_res,t)
  chr=exposure[exposure$SNP==i,]$chr.exposure
  x=exposure[exposure$chr.exposure %in% chr,][,c("chr.exposure","SNP","pval.exposure")]
  y=outcome[outcome$chr %in% chr,][,c("chr","SNP","pval.outcome")]
  
  colnames(x)=c("chr","rsid","pval")
  colnames(y)=c("chr","rsid","pval")
  pdf(paste0(i,"_locuscompare_efinn-b-I9_ANGINA_validation.pdf"))
  p1=locuscompare(in_fn1 = x, in_fn2 = y, title1 = paste0(i,' pQTL'), 
                  title2 = 'Angina Outcome GWAS',snp=i,combine=T)
  print(p1)
  dev.off()
}

write.csv(coloc_res,'coloc.abf_result_efinn-b-I9_ANGINA_validation.csv',row.names=T) 


direct_nodup=directionality_test(harmonised_nodup[harmonised_nodup$exposure %in%sig_protein$Protein,])
write.table(direct_nodup,"Steiger_directionality_test_sig_protein_efinn-b-I9_ANGINA_validation.csv",sep=",",row.names = F)
direct_dup=directionality_test(harmonised_dup)

steiger_nodup=steiger_filtering(harmonised_nodup[harmonised_nodup$exposure %in%sig_protein$Protein,])
steiger_dup=steiger_filtering(harmonised_dup)
write.table(steiger_nodup,"Steiger_filter_sig_protein_efinn-b-I9_ANGINA_validation.csv",sep=",",row.names = F)
save.image("efinn-b-I9_ANGINA_validation.rdata")



#sig_protein$Protein "TGFB1","TIRAP","SWAP70","PCSK9","TMEM106B","APOA5","LPA"
#sig_protein$SNP "rs1800470","rs111577916","rs415895","rs191448950","rs10950398","rs964184","rs55730499"
library(TwoSampleMR)
library(ggplot2)
library(readxl)
library(dplyr)
library(coloc)
#exposure_data<-read_excel("Maik-Pietzner-Supplementary_tables.xlsx",sheet = 4) %>%as.data.frame()
exposure_data<-read.csv("filted_protein_exposure.csv",header = T,stringsAsFactors = FALSE) 
#exposure_data<-exposure_data[,c("HGNC.symbol.protein","rsID","cis.trans","Chr","Position","EA","NEA",
#"EAF","Effect","SE","P.value")]

#colnames(exposure_data)<-c("Phenotype","SNP","cis_trans","chr","pos","effect_allele","other_allele",
# "eaf","beta","se","pval")

exposure<-format_data(exposure_data,type="exposure",header=T,snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele"
                      ,other_allele_col="other_allele",pval_col="pval",eaf_col="eaf",chr_col="chr",pos_col="pos",
                      samplesize_col="sample_size",phenotype_col = "Phenotype" )


exposure<-subset(exposure,pval.exposure<0.00000005 |pval.exposure==0.00000005 & !((chr.exposure=="chr6") & (pos.exposure>26000000 &pos.exposure<34000000)))
exposure<-clump_data(exposure, clump_kb=10000, clump_r2=0.001)


#exposure_data$sample_size=10708
#exposure_data$Phenotype<-gsub(pattern = "\'",x =exposure_data$Phenotype,replacement = "" )
#grep(pattern = "\'",x=exposure_data$Phenotype)

#exposure_data=exposure_data[exposure_data$cis_trans=="cis",]

#exposure_data$effect_allele<-gsub(pattern = "D",x=exposure_data$effect_allele,replacement = "G")
#exposure_data$other_allele<-gsub(pattern = "I",x=exposure_data$other_allele,replacement = "A")

#exposure<-format_data(exposure_data,type="exposure",header=T,snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele"
#                     ,other_allele_col="other_allele",pval_col="pval",eaf_col="eaf",chr_col="chr",pos_col="pos",
#                    samplesize_col="sample_size",phenotype_col = "Phenotype" )

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
ieugwasr::api_status()
library(ieugwasr)
#devtools::install_github("explodecomputer/plinkbinr")
library(plinkbinr)
get_plink_exe()


#clump_df<-ld_clump_local(dplyr::tibble(rsid=exposure$SNP, pval=exposure$pval.exposure, id=exposure$id.exposure), clump_kb=10000, clump_r2=0.001,clump_p=1,
#plink_bin = "D:/R-4.3.1/library/pink/plink.exe",
#bfile = "C:/Users/wangdexin/Desktop/EUR_ref/EUR") %>%as.data.frame()
#colnames(clump_df)<-c("SNP","pval.exposure","id.exposure")
#exposure<-merge(clump_df[,c(1,3)],exposure,by = c("SNP","id.exposure"))


length(unique(exposure$exposure)) 
dup_gene<-exposure[duplicated(exposure$exposure) ,]$exposure

exposure_dup<-exposure[which(exposure$exposure %in% dup_gene) ,]
exposure_nodup<-exposure[-which(exposure$exposure %in% dup_gene),]
table(exposure_nodup$exposure %in% exposure_dup$exposure)

outcome<-extract_outcome_data(snps = exposure$SNP,outcomes="ukb-a-532")

harmonised_nodup<- harmonise_data(exposure_nodup, outcome)
res_single_nodup<- mr_singlesnp(dat=harmonised_nodup,single_method ="mr_wald_ratio")
res_single_nodup=res_single_nodup[!is.na(res_single_nodup$p),]
dim(res_single_nodup[res_single_nodup$p< 0.05/length(unique(exposure$exposure)),])

res_single_nodup <- generate_odds_ratios(res_single_nodup)
write.csv(res_single_nodup,"ukb-a-532_singlesnp_nodup.csv",row.names = F)

harmonised_dup<- harmonise_data(exposure_dup, outcome)
res_single_dup=mr(dat=harmonised_dup,method_list="mr_ivw")
dim(res_single_dup[res_single_dup$p< 0.05/length(unique(exposure$exposure)),])
write.csv(res_single_dup,"ukb-a-532_singlesnp_dup.csv",row.names = F)


heterogeneity=mr_heterogeneity(harmonise_data(exposure, outcome))
write.csv(heterogeneity,"heterogeneity_ukb-a-532_validation.csv",row.names = F) 
mr_heterogeneity(harmonised_nodup)

dim(res_single_nodup[res_single_nodup$p< 0.05/length(unique(exposure$exposure)),])

sig_protein<-res_single_nodup[res_single_nodup$p<0.05/length(unique(exposure$exposure)),]
sig_protein<-sig_protein[,c(1,5,6,12,13,14,9,8)]
colnames(sig_protein)<-c("Protein","samplesize","SNP","OR","low-CI95%","high-CI95%","P","SE")
library(forestploter)
sig_protein$` ` <- paste(rep(" ", 4), collapse = " ")
tm <- forest_theme(base_size = 12,  
                   ci_pch = 15,   
                   ci_col = "black",   
                   ci_fill = "blue",     
                   ci_alpha = 1,        
                   ci_lty = 1,            
                   ci_lwd = 2,          
                   ci_Theight = 0.2,
                   vertline_lwd = 2,            
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   summary_fill = "yellow",       
                   summary_col = "#4575b4",
)

pdf("sig_protein_forest_ukb-a-532_validation.pdf",width=10)
forest(sig_protein[,c(1,2,3,4,9,5,6,7,8)],
       est = sig_protein$OR,       
       lower = sig_protein$`low-CI95%`,     
       upper = sig_protein$`high-CI95%`,    
       sizes = sig_protein$SE,
       ci_column = 5,   
       ref_line = 1,
       arrow_lab = c("Lower", "Higher"),
       xlim = c(0.5,2 ),
       ticks_at = c(0.5, 1, 2),
       footnote = "",
       theme = tm)
dev.off()

sig_protein$Association<-ifelse(sig_protein$OR>1,"Positive","Negative")

res_single_nodup2<-res_single_nodup[,c(1,9,12)]
res_single_nodup2$Association<-ifelse(
  res_single_nodup2$or>1&res_single_nodup2$p< 0.05/length(unique(exposure$exposure)),"Positive",
  ifelse(res_single_nodup2$or<1&res_single_nodup2$p< 0.05/length(unique(exposure$exposure)),"Negative","Not sig"))

res_single_nodup2$Protein<-ifelse(res_single_nodup2$Association=="Not sig","",res_single_nodup2$exposure)
library(ggplot2)
library(ggrepel)
p1<-ggplot(data=res_single_nodup2, aes(x = or, y = -log10(p), color = Association))+
  geom_point(alpha=0.4, size=3.5)+
  geom_vline(xintercept=1,lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept= -log10(0.05/length(unique(exposure$exposure))),lty=4,col="black",lwd=0.8)+
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5),legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15),axis.title = element_text(size=15),
        panel.background = element_rect(fill = "white",colour = "black",linewidth = 0.8),
        axis.text = element_text(size=15))+
  labs(x="OR",y="-log10(Pvalue)")+
  scale_color_manual(values=c("gray", "#A52A2A","#BA55D3"))+
  geom_label_repel(data = res_single_nodup2, aes(x = or, y = -log10(p), 
                                                 label = Protein),max.overlaps =10,size = 5, box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"), segment.color = "black", 
                   show.legend = FALSE)

ggsave("sig_protein_volca_ukb-a-532_validation.pdf",p1,dpi=500,height = 5,width = 5)

#
library(dplyr)  
library(coloc)
GWASdata=outcome[,c("SNP","chr","pos",'effect_allele.outcome',
                    'other_allele.outcome',"eaf.outcome","beta.outcome","se.outcome",
                    "pval.outcome","samplesize.outcome")]
GWASdata$ncase.outcome=4837
GWASdata$ncontrol.outcome=332362

GWASdata$varbeta <- GWASdata$se.outcome^2
GWASdata$s <- GWASdata$ncase.outcome/GWASdata$samplesize.outcome
GWASdata$z = GWASdata$beta.outcome/GWASdata$se.outcome

QTLdata=exposure[,-c(12,13,14)]

QTLdata$varbeta <- QTLdata$se.exposure^2
QTLdata$MAF <- ifelse(QTLdata$eaf.exposure<0.5,QTLdata$eaf.exposure,1-QTLdata$eaf.exposure)
QTLdata$z = QTLdata$beta.exposure/QTLdata$se.exposure
#
library(tidyverse)
library(patchwork)
library(locuscomparer)
coloc_res=data.frame()
for (i in sig_protein$SNP){
  GWASdata2=GWASdata[GWASdata$SNP==i,]
  QTLdata2=QTLdata[QTLdata$SNP==i,]
  coloc_result <- coloc.abf(dataset1=list(pvalues=GWASdata2$pval.outcome, snp=GWASdata2$SNP, type="cc", s=GWASdata2$s[1], N=GWASdata2$samplesize.outcome[1]), 
                            dataset2=list(pvalues=QTLdata2$pval.exposure,snp=QTLdata2$SNP, type="quant", N=QTLdata2$samplesize[1]), MAF=QTLdata2$MAF)
  
  t=data.frame(t(data.frame(coloc_result[[1]])))
  rownames(t)<-i
  coloc_res=rbind(coloc_res,t)
  chr=exposure[exposure$SNP==i,]$chr.exposure
  x=exposure[exposure$chr.exposure %in% chr,][,c("chr.exposure","SNP","pval.exposure")]
  y=outcome[outcome$chr %in% chr,][,c("chr","SNP","pval.outcome")]
  
  colnames(x)=c("chr","rsid","pval")
  colnames(y)=c("chr","rsid","pval")
  pdf(paste0(i,"_locuscompare_ukb-a-532_validation.pdf"))
  p1=locuscompare(in_fn1 = x, in_fn2 = y, title1 = paste0(i,' pQTL'), 
                  title2 = 'Angina Outcome GWAS',snp=i,combine=T)
  print(p1)
  dev.off()
}

write.csv(coloc_res,'coloc.abf_result_ukb-a-532_validation.csv',row.names=T) 

#

direct_nodup=directionality_test(harmonised_nodup[harmonised_nodup$exposure %in%sig_protein$Protein,])
write.table(direct_nodup,"Steiger_directionality_test_sig_protein_ukb-a-532_validation.csv",sep=",",row.names = F)
direct_dup=directionality_test(harmonised_dup)

steiger_nodup=steiger_filtering(harmonised_nodup[harmonised_nodup$exposure %in%sig_protein$Protein,])
steiger_dup=steiger_filtering(harmonised_dup)
write.table(steiger_nodup,"Steiger_filter_sig_protein_ukb-a-532_validation.csv",sep=",",row.names = F)
save.image("ukb-a-532_validation.rdata")



protein<-c("TIRAP","SWAP70","PCSK9","TMEM106B","APOA5","LPA")
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
ENTREZID<-bitr(geneID = protein,OrgDb = "org.Hs.eg.db",fromType = "SYMBOL",toType = "ENTREZID")
go<-enrichGO(gene=ENTREZID$ENTREZID,OrgDb ="org.Hs.eg.db",keyType ="ENTREZID",pvalueCutoff = 0.05,ont = "ALL",
             minGSSize = 10,qvalueCutoff = 0.05)

write.table(as.data.frame(go),"protein_GO.csv",sep = ",",row.names = F)

pdf("go.pdf")
barplot(go,x="GeneRatio",color = "p.adjust",showCategory = 5,split="ONTOLOGY",font.size = 15)+
  facet_grid(ONTOLOGY~.,scales = "free")
dev.off()


kegg<-enrichKEGG(gene=ENTREZID$ENTREZID,organism = "hsa",keyType = "kegg",pvalueCutoff = 0.1,
                 minGSSize = 10,qvalueCutoff = 0.2)

write.table(as.data.frame(kegg),"protein_KEGG.csv",sep = ",",row.names = F)

pdf("kegg.pdf")
dotplot(kegg,font.size = 15)
dev.off()



expression<-read.table("GSE42148_matrix.txt",sep="\t",header = T)
probe<-read.delim2("GPL13607.txt",sep = "\t",header = T)
probe<-probe[,c(1,6)]
colnames(probe)[1]<-"ID_REF"
expression<-merge(probe,expression,by="ID_REF")
protein<-c("TIRAP","SWAP70","PCSK9","TMEM106B","APOA5","LPA")
#expression=expression[-grep(pattern = "///",x=expression$GeneName),]

table( protein %in% expression$GeneName)
library(dplyr)
#expression<-expression[-grep(x=expression$GeneName,pattern = "^[0-9]"),]
expression<-aggregate(.~GeneName,mean,data=expression[,-1])
table(is.na(expression))

expression<-expression[!expression$GeneName=="",]
#rownames(expression)<-expression$GeneName
#expression$GeneName<-NULL
group<-c(rep("control",11),rep("SA",7),rep("MI",6))
mat<-expression[expression$GeneName %in% protein,]
library(tidyverse)
rownames(mat)<-mat$GeneName
mat$GeneName<-NULL
mat=t(mat) %>%as.data.frame()
mat$group<-group

library(dplyr)
library(ggpubr) 
library(patchwork)
library(ggsci)
library(ggprism)
library(gridExtra)
#mat2<-pivot_longer(data=mat,cols=APOA5:TMEM106B,names_to = "var",values_to ="value") %>%as.data.frame()


p1=ggplot(mat[,c(1,7)], aes(x=group, y=APOA5,color=group)) +
  geom_boxplot(aes(fill=group),
               alpha=0.1)+ 
  geom_jitter()+
  scale_color_manual(values = pal_npg('nrc')(9))+
  scale_fill_manual(values = pal_npg('nrc')(9))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  stat_compare_means(comparisons=list(c("control","MI"),c("control","SA"),c("MI","SA")),
                     label.y = c(5.2,5.4,5.7),
                     method="t.test",size=8,color="black")+
  labs(x="")+
  theme_prism(palette = "flames",
              base_fontface = "plain", 
              base_family = "serif", 
              base_size = 16, 
              base_line_size = 0.8, 
              axis_text_angle = 45)+ 
  scale_fill_prism(palette = "flames")

p2=ggplot(mat[,c(2,7)], aes(x=group, y=LPA,color=group)) +
  geom_boxplot(aes(fill=group),
               alpha=0.1)+ 
  geom_jitter()+
  scale_color_manual(values = pal_npg('nrc')(9))+
  scale_fill_manual(values = pal_npg('nrc')(9))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  stat_compare_means(comparisons=list(c("control","MI"),c("control","SA"),c("MI","SA")),
                     label.y = c(7.1,7.3,7.5),
                     method="t.test",size=8,color="black")+
  labs(x="")+
  theme_prism(palette = "flames",
              base_fontface = "plain", 
              base_family = "serif", 
              base_size = 16,  
              base_line_size = 0.8, 
              axis_text_angle = 45)+ 
  scale_fill_prism(palette = "flames")

p3=ggplot(mat[,c(3,7)], aes(x=group, y=PCSK9,color=group)) +
  geom_boxplot(aes(fill=group),
               alpha=0.1)+ 
  geom_jitter()+
  scale_color_manual(values = pal_npg('nrc')(9))+
  scale_fill_manual(values = pal_npg('nrc')(9))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  stat_compare_means(comparisons=list(c("control","MI"),c("control","SA"),c("MI","SA")),
                     label.y = c(2.3,2.5,2.7),
                     method="t.test",size=8,color="black")+
  labs(x="")+
  theme_prism(palette = "flames",
              base_fontface = "plain", 
              base_family = "serif", 
              base_size = 16,  
              base_line_size = 0.8, 
              axis_text_angle = 45)+
  scale_fill_prism(palette = "flames")


p4=ggplot(mat[,c(4,7)], aes(x=group, y=SWAP70,color=group)) +
  geom_boxplot(aes(fill=group),
               alpha=0.1)+ 
  geom_jitter()+
  scale_color_manual(values = pal_npg('nrc')(9))+
  scale_fill_manual(values = pal_npg('nrc')(9))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  stat_compare_means(comparisons=list(c("control","MI"),c("control","SA"),c("MI","SA")),
                     label.y = c(9.3,9.5,9.68),
                     method="t.test",size=8,color="black")+
  labs(x="")+
  theme_prism(palette = "flames",
              base_fontface = "plain", 
              base_family = "serif",
              base_size = 16, 
              base_line_size = 0.8,
              axis_text_angle = 45)+ 
  scale_fill_prism(palette = "flames")


p6=ggplot(mat[,c(5,7)], aes(x=group, y=TIRAP,color=group)) +
  geom_boxplot(aes(fill=group),
               alpha=0.1)+ 
  geom_jitter()+
  scale_color_manual(values = pal_npg('nrc')(9))+
  scale_fill_manual(values = pal_npg('nrc')(9))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  stat_compare_means(comparisons=list(c("control","MI"),c("control","SA"),c("MI","SA")),
                     label.y = c(5.98,6.15,6.25),
                     method="t.test",size=8,color="black")+
  labs(x="")+
  theme_prism(palette = "flames",
              base_fontface = "plain", 
              base_family = "serif", 
              base_size = 16,  
              base_line_size = 0.8, 
              axis_text_angle = 45)+ 
  scale_fill_prism(palette = "flames")

p7=ggplot(mat[,c(6,7)], aes(x=group, y=TMEM106B,color=group)) +
  geom_boxplot(aes(fill=group),
               alpha=0.1)+ 
  geom_jitter()+
  scale_color_manual(values = pal_npg('nrc')(9))+
  scale_fill_manual(values = pal_npg('nrc')(9))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  stat_compare_means(comparisons=list(c("control","MI"),c("control","SA"),c("MI","SA")),
                     label.y = c(9.9,10.1,10.28),
                     method="t.test",size=8,color="black")+
  labs(x="")+
  theme_prism(palette = "flames",
              base_fontface = "plain", 
              base_family = "serif", 
              base_size = 16,  
              base_line_size = 0.8, 
              axis_text_angle = 45)+
  scale_fill_prism(palette = "flames")

p8<-grid.arrange(p1,p2,p3,p4,p6,p7,ncol=3)
ggsave("sig-protein-mRNA.pdf",p8,width=12,height = 10)
save.image("plot.rdata")




library(TwoSampleMR)
library(ggplot2)
library(readxl)
library(dplyr)
library(coloc)
library(ieugwasr)
mr_res=data.frame()
plt_res=data.frame()
het_res=data.frame()
loo_res=data.frame()

protein<-c("TGFB1","TIRAP","SWAP70","PCSK9","TMEM106B","APOA5","LPA")
snp<-c("rs1800470","rs111577916","rs415895","rs191448950","rs10950398","rs964184","rs55730499")
#apolipoprotein B,cholesterol, Systolic blood pressure,	Diastolic blood pressure,Cardiovascular disease,Coronary artery disease,Triglycerides,Myocardial infarction,height
out<-c("ebi-a-GCST90013994","ebi-a-GCST90092985","ieu-b-38","ukb-b-7992","finn-b-I9_CVD","ebi-a-GCST90013868","ieu-b-111","ebi-a-GCST011364","ebi-a-GCST90018959")


#remove MHC's snp
for (i in out){
  multi_exposure<-extract_instruments(outcomes =i)  
  multi_exposure<-subset(multi_exposure,pval.exposure<0.00000005& !((chr.exposure=="6") & (pos.exposure>26000000 &pos.exposure<34000000)))
  outcome<-extract_outcome_data(snps = multi_exposure$SNP,outcomes = "ebi-a-GCST90018793")
  harmonise <- harmonise_data(multi_exposure, outcome)
  res <- mr(harmonise)
  
  res <- generate_odds_ratios(res)
  mr_res=rbind(res,mr_res)
  #heterogeneity
  het<-mr_heterogeneity(harmonise)
  het_res=rbind(het,het_res)
  #pleiotropy test
  plt<-mr_pleiotropy_test(harmonise)
  plt_res=rbind(plt,plt_res)
  ###Leave one out analysis
  loo<- mr_leaveoneout(harmonise)

  loo_res=rbind(loo,loo_res)
}
write.table(mr_res,"mr_result_to_Angina.csv",sep = ",",row.names = F)
write.table(het_res,"het_res_to_Angina.csv",sep = ",",row.names = F)
write.table(plt_res,"plt_res_to_Angina.csv",sep = ",",row.names = F)
write.table(loo_res,"loo_res_to_Angina.csv",sep = ",",row.names = F)


library(ggplot2)
save.image("risk_factor_mr.rdata")

load("risk_factor_mr.rdata")
sub("\\s*\\|\\|.*", "", mr_res$exposure)
mr_res$exposure_trait<-c(c(rep("Height",5)),c(rep("Myocardial infarction",5)),
                         c(rep("triglycerides",5)),c(rep("Coronary artery disease",5)),c(rep("Cardiovascular diseases",5)),
                         c(rep("Diastolic blood pressure",5)),c(rep("systolic blood pressure",5)),c(rep("Total cholesterol levels",5)),
                         c(rep("Apolipoprotein B levels",5)))

my36colors <-c ('#E5D2DD','#53A85F','#F1BB72','#F3B1A0','#D6E7A3','#57C3F3','#476D87','#E95C59','#E59CC4','#AB3282','#23452F','#BD956A',
                '#8C549C','#585658','#9FA3A8','#E0D4CA','#5F3D69','#C5DEBA','#58A4C3','#E4C755','#F7F398','#AA9A59','#E63863','#E39A35',
                '#C1E6F3','#6778AE','#91D0BE','#B53E2B','#712820','#DCC1DD','#CCE0F5','#CCC9E6','#625D9E','#68A180','#3A6963','#968175')

p1=ggplot(mr_res,aes(x= -log10(pval),y=exposure_trait,color=method))+
  geom_point(mapping = aes(size= -log10(pval)))+
  labs(y="",x="",title = "Risk factors to Angina Pectoris")+
  geom_vline(xintercept= -log10(0.05),lty=4,col="black",lwd=0.8)+
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5),legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15),axis.title = element_text(size=15),
        panel.background = element_rect(fill = "white",colour = "black",linewidth = 0.8),
        axis.text = element_text(size=15,colour = "black"),plot.title = element_text(hjust=0.5,colour = "black",size = 15))+
  scale_color_manual(name="Method",values = my36colors[c(2,3,6,8,13)])

save.image("risk_factor_mr.rdata")


#Reverse causal effect 

library(TwoSampleMR)
library(ggplot2)
library(readxl)
library(dplyr)
library(coloc)
library(ieugwasr)
mr_res=data.frame()
plt_res=data.frame()
het_res=data.frame()
loo_res=data.frame()

protein<-c("TGFB1","TIRAP","SWAP70","PCSK9","TMEM106B","APOA5","LPA")
snp<-c("rs1800470","rs111577916","rs415895","rs191448950","rs10950398","rs964184","rs55730499")
#apolipoprotein B,cholesterol, Systolic blood pressure,	Diastolic blood pressure,Cardiovascular disease,Coronary artery disease,Triglycerides,Myocardial infarction,height
out<-c("ebi-a-GCST90013994","ebi-a-GCST90092985","ieu-b-38","ukb-b-7992","finn-b-I9_CVD","ebi-a-GCST90013868","ieu-b-111","ebi-a-GCST011364","ebi-a-GCST90018959")


#remove MHC's snp
for (i in out){
  multi_exposure<-extract_instruments(outcomes ="ebi-a-GCST90018793")  
  multi_exposure<-subset(multi_exposure,pval.exposure<0.00000005& !((chr.exposure=="6") & (pos.exposure>26000000 &pos.exposure<34000000)))
  outcome<-extract_outcome_data(snps = multi_exposure$SNP,outcomes = i)
  harmonise <- harmonise_data(multi_exposure, outcome)
  res <- mr(harmonise)
  loo<- mr_leaveoneout(harmonise)
  #remove the p value >0.05 snp
  loo<-subset(loo,p>0.05)
  remove_name<-loo$SNP
  #remove_name<-remove_name[remove_name!="All"]
  multi_exposure<-multi_exposure[!(multi_exposure$SNP %in% remove_name),]
  harmonise <- harmonise_data(multi_exposure, outcome)
  if (dim(multi_exposure)[1]!=0 & dim(harmonise)[1]!=0){
    res <- generate_odds_ratios(res)
    res <- mr(harmonise)
    mr_res=rbind(res,mr_res)
    #heterogeneity
    het<-mr_heterogeneity(harmonise)
    het_res=rbind(het,het_res)
    #pleiotropy test
    plt<-mr_pleiotropy_test(harmonise)
    plt_res=rbind(plt,plt_res)
    ###Leave one out analysis
    loo<- mr_leaveoneout(harmonise)
    loo_res=rbind(loo,loo_res)}
}

write.table(mr_res,"reverse_mr_result_to_Angina.csv",sep = ",",row.names = F)
write.table(het_res,"reverse_het_res_to_Angina.csv",sep = ",",row.names = F)
write.table(plt_res,"reverse_plt_res_to_Angina.csv",sep = ",",row.names = F)
write.table(loo_res,"reverse_loo_res_to_Angina.csv",sep = ",",row.names = F)


sub("\\s*\\|\\|.*", "", mr_res$outcome)
mr_res$outcome_trait<-c(c(rep("Myocardial infarction",5)),c(rep("Coronary artery disease",5)),c(rep("Cardiovascular diseases",5)),
                        c(rep("systolic blood pressure",5)))

my36colors <-c ('#E5D2DD','#53A85F','#F1BB72','#F3B1A0','#D6E7A3','#57C3F3','#476D87','#E95C59','#E59CC4','#AB3282','#23452F','#BD956A',
                '#8C549C','#585658','#9FA3A8','#E0D4CA','#5F3D69','#C5DEBA','#58A4C3','#E4C755','#F7F398','#AA9A59','#E63863','#E39A35',
                '#C1E6F3','#6778AE','#91D0BE','#B53E2B','#712820','#DCC1DD','#CCE0F5','#CCC9E6','#625D9E','#68A180','#3A6963','#968175')

p2=ggplot(mr_res,aes(x= -log10(pval),y=outcome_trait,color=method))+
  geom_point(mapping = aes(size= -log10(pval)))+
  labs(y="",x="-log10(p-value)",title = "Angina Pectoris to risk factors")+
  geom_vline(xintercept= -log10(0.05),lty=4,col="black",lwd=0.8)+
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5),legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15),axis.title = element_text(size=15),
        panel.background = element_rect(fill = "white",colour = "black",linewidth = 0.8),
        axis.text = element_text(size=15,colour = "black"),plot.title = element_text(hjust=0.5,colour = "black",size = 15),
        axis.ticks = element_blank())+
  scale_color_manual(name="Method",values = my36colors[c(2,3,6,8,13)])

save.image("reverse_risk_factor_mr.rdata")
library(gridExtra)
p3=grid.arrange(p1,p2,ncol=1)
ggsave("risk factors.pdf",p3,height = 7,width =9)


library(TwoSampleMR)
library(ggplot2)
library(readxl)
library(dplyr)
library(coloc)
library(ieugwasr)
mr_res=data.frame()
#plt_res=data.frame()
#het_res=data.frame()
#loo_res=data.frame()

out<-c("ebi-a-GCST90013994","ebi-a-GCST90092985","ieu-b-38","ukb-b-7992","finn-b-I9_CVD","ebi-a-GCST90013868","ieu-b-111","ebi-a-GCST011364","ebi-a-GCST90018959")
protein<-c("TGFB1","TIRAP","SWAP70","PCSK9","TMEM106B","APOA5","LPA")
snp<-c("rs1800470","rs111577916","rs415895","rs191448950","rs10950398","rs964184","rs55730499")

posi_exposure<-read.csv("../filted_protein_exposure.csv",header = T,stringsAsFactors = FALSE)
posi_exposure<-posi_exposure[posi_exposure$SNP %in%snp,]
posi_exposure<-format_data(posi_exposure,type="exposure",header=T,snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele"
                           ,other_allele_col="other_allele",pval_col="pval",eaf_col="eaf",chr_col="chr",pos_col="pos",
                           samplesize_col="sample_size",phenotype_col = "Phenotype" )

posi_exposure<-clump_data(posi_exposure,clump_r2=0.001,clump_kb=10000)
posi_exposure<-subset(posi_exposure,pval.exposure<0.00000005& !((chr.exposure=="chr6") & (pos.exposure>26000000 &pos.exposure<34000000)))

for (i in out){
  outcome<-extract_outcome_data(snps = posi_exposure$SNP,outcomes=i)
  harmonised<- harmonise_data(posi_exposure, outcome)
  res_single<- mr_singlesnp(dat=harmonised,single_method = "mr_wald_ratio")
  res_single=res_single[!is.na(res_single$p),]
  res_single <- generate_odds_ratios(res_single)
  
  mr_res<-rbind(res_single,mr_res)
  
  #loo<- mr_leaveoneout(harmonised)
  #remove the p value >0.05 snp
  #loo<-subset(loo,p>0.05)
  #remove_name<-loo$SNP
  #remove_name<-remove_name[remove_name!="All"]
  #multi_exposure<-multi_exposure[!(multi_exposure$SNP %in% remove_name),]
  #harmonised <- harmonise_data(multi_exposure, outcome)
  #if (dim(multi_exposure)[1]!=0 & dim(harmonise)[1]!=0){
  #
  # res <- mr(harmonised)
  #mr_res=rbind(res,mr_res)
  #heterogeneity
  #het<-mr_heterogeneity(harmonised)
  # het_res=rbind(het,het_res)
  #pleiotropy test
  # plt<-mr_pleiotropy_test(harmonised)
  # plt_res=rbind(plt,plt_res)
  ###Leave one out analysis
  # loo<- mr_leaveoneout(harmonised)
  # loo_res=rbind(loo,loo_res)}
}
write.table(mr_res,"protein_to_risk_factor_mr_result.csv",sep = ",",row.names = F)



mr_res$outcome_train<-c(rep("Height",6),c(rep("Myocardial infarction",7)),
                        c(rep("triglycerides",7),c(rep("Coronary artery disease",5)),c(rep("Cardiovascular diseases",7)),
                          c(rep("Diastolic blood pressure",7)),c(rep("systolic blood pressure",6)),c(rep("Total cholesterol levels",7)),
                          c(rep("Apolipoprotein B levels",5))))
sig_protein<-mr_res[,c(15,1,5,6,12,13,14,9,8)]
colnames(sig_protein)<-c("Trait","Protein","samplesize","SNP","OR","low-CI95%","high-CI95%","P","SE")
library(forestploter)
sig_protein$` ` <- paste(rep(" ", dim(sig_protein)[1]), collapse = " ")

tm <- forest_theme(base_size = 12,  
                   ci_pch = 15,  
                   ci_col = "black",   
                   ci_fill = "blue",     
                   ci_alpha = 1,        
                   ci_lty = 1,            
                   ci_lwd = 2,         
                   ci_Theight = 0.2,
                   vertline_lwd = 2,            
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   summary_fill = "yellow",    
                   summary_col = "#4575b4",
)

pdf("protein_to_risk_factors_forest.pdf",width=18,height = 15)
forest(sig_protein[,c(1,2,3,4,10,5,6,7,8,9)],
       est = sig_protein$OR,     
       lower = sig_protein$`low-CI95%`,     
       upper = sig_protein$`high-CI95%`,      
       sizes = sig_protein$SE,
       ci_column = 5,   
       ref_line = 1,
       arrow_lab = c("Lower", "Higher"),
       xlim = c(0.25,1.5 ),
       ticks_at = c(0.25, 1, 1.5),
       footnote = "",
       theme = tm)
dev.off()

save.image("protein_to_risk_factors.rdata")





library(TwoSampleMR)
library(ggplot2)
library(readxl)
library(dplyr)
library(coloc)
library(ieugwasr)
mr_res=data.frame()
#plt_res=data.frame()
#het_res=data.frame()
#loo_res=data.frame()

out<-c("ebi-a-GCST90013994","ebi-a-GCST90092985","ieu-b-38","ukb-b-7992","finn-b-I9_CVD","ebi-a-GCST90013868","ieu-b-111","ebi-a-GCST011364","ebi-a-GCST90018959")
protein<-c("TGFB1","TIRAP","SWAP70","PCSK9","TMEM106B","APOA5","LPA")
snp<-c("rs1800470","rs111577916","rs415895","rs191448950","rs10950398","rs964184","rs55730499")

posi_outcome<-read.csv("../filted_protein_exposure.csv",header = T,stringsAsFactors = FALSE)
posi_outcome<-posi_outcome[posi_outcome$SNP %in%snp,]
posi_outcome<-format_data(posi_outcome,type="outcome",header=T,snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele"
                          ,other_allele_col="other_allele",pval_col="pval",eaf_col="eaf",chr_col="chr",pos_col="pos",
                          samplesize_col="sample_size",phenotype_col = "Phenotype" )

direct_res<-data.frame()
steiger_res<-data.frame()
for (i in out){
  exposure<-extract_instruments(outcomes = i,clump =F)
  if(i=="finn-b-I9_CVD"){
    exposure$samplesize.exposure=c(218792)
  }
  if(i=="ieu-b-111")
  { exposure$samplesize.exposure=c(441016)}
  
  #exposure<-clump_data(exposure,clump_r2=0.001,clump_kb=10000)
  #exposure<-subset(exposure,pval.exposure<0.00000005& !((chr.exposure=="chr6") & (pos.exposure>26000000 &pos.exposure<34000000)))
  harmonised<- harmonise_data(exposure, posi_outcome)
  
  # harmonised=harmonised[harmonised$SNP %in%snp,]
  #if(dim(harmonised)[1]!=0){
  
  # res_single=mr(dat=harmonised)
  
  #}
  #mr_res=rbind(res_single,mr_res)
  direct=directionality_test(harmonised)
  direct_res=rbind(direct,direct_res)

  steiger=steiger_filtering(harmonised)
  steiger_res<-rbind(steiger,steiger_res)
  
}

write.table(direct,"Steiger_directionality_risk_factor_to_protrin.csv",sep=",",row.names = F)
write.table(harmonised,"Steiger_filter_risk_factor_to_protrin.csv",sep=",",row.names = F)
