library(maftools)
maf_BRCA= read.delim("/Users/cyu/Downloads/mutation_process/BRCAmut.txt", header =TRUE, sep ="\t",stringsAsFactors = FALSE)
BRCA = read.maf(maf = maf_BRCA)

maf_CESC= read.delim("/Users/cyu/Downloads/mutation_process/CESCmut.txt", header =TRUE, sep ="\t",stringsAsFactors = FALSE)
CESC = read.maf(maf = maf_CESC)
# MAF file summarized by Genes 
getGeneSummary(CESC)
# MAF file summarized by samples (Tumor Sample Barcode)
samplesummary<-getSampleSummary(CESC)
head(samplesummary)

maf_COAD= read.delim("/Users/cyu/Downloads/mutation_process/COADmut.txt", header =TRUE, sep ="\t",stringsAsFactors = FALSE)
COAD = read.maf(maf = maf_COAD)
# MAF file summarized by Genes 
getGeneSummary(COAD)
# MAF file summarized by samples (Tumor Sample Barcode)
samplesummary<-getSampleSummary(COAD)
head(samplesummary)

maf_COADREAD= read.delim("/Users/cyu/Downloads/mutation_process/COADREADmut.txt", header =TRUE, sep ="\t",stringsAsFactors = FALSE)
COADREAD = read.maf(maf = maf_COADREAD)
# MAF file summarized by Genes 
getGeneSummary(COADREAD)
# MAF file summarized by samples (Tumor Sample Barcode)
samplesummary<-getSampleSummary(COADREAD)
head(samplesummary)

maf_GBM= read.delim("/Users/cyu/Downloads/mutation_process/GBMmut.txt", header =TRUE, sep ="\t",stringsAsFactors = FALSE)
GBM = read.maf(maf = maf_GBM)
# MAF file summarized by Genes 
getGeneSummary(GBM)
# MAF file summarized by samples (Tumor Sample Barcode)
samplesummary<-getSampleSummary(GBM)
head(samplesummary)

maf_GBMLGG= read.delim("/Users/cyu/Downloads/mutation_process/GBMLGGmut.txt", header =TRUE, sep ="\t",stringsAsFactors = FALSE)
GBMLGG = read.maf(maf = maf_GBMLGG)
# MAF file summarized by Genes 
getGeneSummary(GBMLGG)
# MAF file summarized by samples (Tumor Sample Barcode)
samplesummary<-getSampleSummary(GBMLGG)
head(samplesummary)

maf_LUAD= read.delim("/Users/cyu/Downloads/mutation_process/LUADmut.txt", header =TRUE, sep ="\t",stringsAsFactors = FALSE)
LUAD = read.maf(maf = maf_LUAD)
# MAF file summarized by Genes 
getGeneSummary(LUAD)
# MAF file summarized by samples (Tumor Sample Barcode)
samplesummary<-getSampleSummary(LUAD)
head(samplesummary)

maf_LUSC= read.delim("/Users/cyu/Downloads/mutation_process/LUSCmut.txt", header =TRUE, sep ="\t",stringsAsFactors = FALSE)
LUSC = read.maf(maf = maf_LUSC)
# MAF file summarized by Genes 
getGeneSummary(LUSC)
# MAF file summarized by samples (Tumor Sample Barcode)
samplesummary<-getSampleSummary(LUSC)
head(samplesummary)

maf_OV= read.delim("/Users/cyu/Downloads/mutation_process/OVmut.txt", header =TRUE, sep ="\t",stringsAsFactors = FALSE)
OV = read.maf(maf = maf_OV)
# MAF file summarized by Genes 
getGeneSummary(OV)
# MAF file summarized by samples (Tumor Sample Barcode)
samplesummary<-getSampleSummary(OV)
head(samplesummary)

maf_PAAD= read.delim("/Users/cyu/Downloads/mutation_process/PAADmut.txt", header =TRUE, sep ="\t",stringsAsFactors = FALSE)
PAAD = read.maf(maf = maf_PAAD)
# MAF file summarized by Genes 
getGeneSummary(PAAD)
# MAF file summarized by samples (Tumor Sample Barcode)
samplesummary<-getSampleSummary(PAAD)
head(samplesummary)

maf_STAD= read.delim("/Users/cyu/Downloads/mutation_process/STADmut.txt", header =TRUE, sep ="\t",stringsAsFactors = FALSE)
STAD = read.maf(maf = maf_STAD)
# MAF file summarized by Genes 
getGeneSummary(STAD)
# MAF file summarized by samples (Tumor Sample Barcode)
samplesummary<-getSampleSummary(STAD)
head(samplesummary)

maf_STES= read.delim("/Users/cyu/Downloads/mutation_process/STESmut.txt", header =TRUE, sep ="\t",stringsAsFactors = FALSE)
STES = read.maf(maf = maf_STES)
# MAF file summarized by Genes 
getGeneSummary(STES)
# MAF file summarized by samples (Tumor Sample Barcode)
samplesummary<-getSampleSummary(STES)
head(samplesummary)


library(ggplot2)
library(RColorBrewer)

## summarise patient wise ##
# distribution of variant per sample #
#https://github.com/margarethannum/gnomeR/blob/19fe064fb650bdff63d3e476471bca8789d7b09c/R/maf-summary.R
#nb.cols <- length(unique(maf_GBM$Variant_Classification))
col = colorRampPalette(brewer.pal(8, "Accent"))(17)
names(col) = c("3'UTR","5'Flank","5'UTR","De_novo_Start_InFrame","De_novo_Start_OutOfFrame",
               'Frame_Shift_Del','Frame_Shift_Ins',"IGR",'In_Frame_Del','Missense_Mutation','In_Frame_Ins','Intron',
               "Indel","Missense",'Nonstop_Mutation',"Read-through",'Transaltion_Start_Site' )

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/brca_variant.pdf",width = 20, height =15)
ggplot(maf_BRCA,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
#  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
           axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
           legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/cesc_variantpsample.pdf",width = 20, height =15)
ggplot(maf_CESC,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/coad_variantpsample.pdf",width = 20, height =15)
ggplot(maf_COAD,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/COADREAD_variantpsample.pdf",width = 20, height =15)
ggplot(maf_COADREAD,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/gbm_variantpsample.pdf",width = 20, height =15)
ggplot(maf_GBM,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/gbmlgg_variantpsample.pdf",width = 20, height =15)
ggplot(maf_GBMLGG,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/luad_variantpsample.pdf",width = 20, height =15)
ggplot(maf_LUAD,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/lusc_variantpsample.pdf",width = 20, height =15)
ggplot(maf_LUSC,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/ov_variantpsample.pdf",width = 20, height =15)
ggplot(maf_OV,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/paad_variantpsample.pdf",width = 20, height =15)
ggplot(maf_PAAD,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/stad_variantpsample.pdf",width = 20, height =15)
ggplot(maf_STAD,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/stes_variantpsample.pdf",width = 20, height =15)
ggplot(maf_STES,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()





pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/brca_variantpsample.jpeg",width = 20, height =15)
ggplot(maf_BRCA,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/cesc_variantpsample.jpeg",width = 20, height =15)
ggplot(maf_CESC,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/coad_variantpsample.jpeg",width = 20, height =15)
ggplot(maf_COAD,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/COADREAD_variantpsample.jpeg",width = 20, height =15)
ggplot(maf_COADREAD,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/gbm_variantpsample.jpeg",width = 20, height =15)
ggplot(maf_GBM,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/gbmlgg_variantpsample.jpeg",width = 20, height =15)
ggplot(maf_GBMLGG,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/luad_variantpsample.jpeg",width = 20, height =15)
ggplot(maf_LUAD,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/lusc_variantpsample.jpeg",width = 20, height =15)
ggplot(maf_LUSC,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/ov_variantpsample.jpeg",width = 20, height =15)
ggplot(maf_OV,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/paad_variantpsample.jpeg",width = 20, height =15)
ggplot(maf_PAAD,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/stad_variantpsample.jpeg",width = 20, height =15)
ggplot(maf_STAD,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()

pdf("/Users/cyu/Downloads/mutation_process/variants_per_cancertpye/stes_variantpsample.jpeg",width = 20, height =15)
ggplot(maf_STES,aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
  #  ggtitle("Variants per sample") + 
  ylab("Variant Count") + xlab("Sample")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))+
  scale_fill_manual(values = col)
dev.off()








