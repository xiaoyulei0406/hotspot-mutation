tumor.type <- c(rep("BRCA",7),rep("CESC",7),rep("COAD",7),rep("COADREAD",7),rep("GBM",7),rep("GBMLGG",7),
                rep("LUAD",7),rep("LUSC",7),rep("OV",7),rep("PAAD",7),rep("STAD",7),rep("STES",7))
variant.classification <- rep(c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                      "Missense_Mutation","Nonstop_Mutation","Translation_Start_Site"),12)
value <-c(3162,	 2567, 627,  161,	55063,  133, 0,
          220,	 271,  73,   47,	26606,  84,  0,
          1230,	 803,  146,  21,	40735,  32,  2,
          1406,	 964,  176,  29,	55502,  38,  2,
          566,	 217,  214,	 28,  14213,  17,  71,
          941,	 309,  371,	 33,  20311,  24,  71,
          1701,	 1006, 259,  29,  47700,	56,  107,
          540,	 141,  98,   20,  42890,	59,  24,
          477,	 172,  165,	 44,  13579,	16,  0,
          126,	 53,	 36,   4,   19899,	6,   50,
          11082, 4107, 1021, 114, 87092,	93,	 247,
          13275, 4900, 1551, 196, 120377,	140, 247)
data <- data.frame(tumor.type,variant.classification,value)
write.table(data ,file = "/Users/yuchunchun/Documents/work/work/TCGAmuation/processed_data/2variantsamongcancer/ggplot2varianttypeamongcancer.txt",sep="\t", quote=F)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(RColorBrewer)

##colors to be selected
#col = colorRampPalette(brewer.pal(8, "Accent"))(7)
#names(col) = c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','Missense_Mutation','In_Frame_Ins',
#               'Nonstop_Mutation','Transaltion_Start_Site' )
col1 = colorRampPalette(brewer.pal(8, "Paired"))(7)
names(col1) = c('Frame_Shift_Del','Missense_Mutation','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins',
               'Nonstop_Mutation','Transaltion_Start_Site' )

pdf("/Users/yuchunchun/Documents/work/work/TCGAmuation/processed_data/2variantsamongcancer/varianttpyeamongcancertypes1.pdf",width = 20, height =15)
ggplot(data, aes(fill=variant.classification, y=value, x=tumor.type)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = col1) +
  #  ggtitle("Variants among Cancer types") + 
  ylab("Variant Count") + xlab("Cancer Type")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=23,face="bold"),legend.title = element_text( size=22,face="bold"),
        legend.text = element_text(size=30))
dev.off()
