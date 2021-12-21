library(plyr)

input_data_dir = "/Users/cyu/Documents/work/mc3data/topmutations/pancan/raw_data/"
output_data_dir <- "/Users/cyu/Documents/work/mc3data/topmutations/pancan/raw_data/"
File <-paste0(input_data_dir,"mc3.top1000mutations.txt")
maffile <-read.delim(file= File, header =TRUE, sep ="\t",stringsAsFactors = FALSE)

maffile["All_patients"]
maffile[, "new"] <- maffile[, "Freq"] / maffile[, "All_patients"]
maffile[, "mutRates"] <- maffile[, "Freq"] / maffile[, "All_patients"]
write.table(maffile,paste0(output_data_dir,"top1000muts_all.txt"),sep="\t", quote=F,row.names = FALSE)


#
input_data_dir = "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/CHOL/"
output_data_dir <- "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/CHOL/"
File <-paste0(input_data_dir,"CHOL.allaa_fre.txt")
maffile <-read.delim(file= File, header =TRUE, sep ="\t",stringsAsFactors = FALSE)
maffile[, "mutRates"] <- maffile[, "Freq"] / 36
write.table(maffile,paste0(output_data_dir,"mut_allfre_CHOL.txt"),sep="\t", quote=F,row.names = FALSE)

#
input_data_dir = "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/COAD/"
output_data_dir <- "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/COAD/"
File <-paste0(input_data_dir,"COAD.allaa_fre.txt")
maffile <-read.delim(file= File, header =TRUE, sep ="\t",stringsAsFactors = FALSE)
maffile[, "mutRates"] <- maffile[, "Freq"] / 394
write.table(maffile,paste0(output_data_dir,"mut_allfre_COAD.txt"),sep="\t", quote=F,row.names = FALSE)

#
input_data_dir = "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/LUAD/"
output_data_dir <- "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/LUAD/"
File <-paste0(input_data_dir,"LUAD.allaa_fre.txt")
maffile <-read.delim(file= File, header =TRUE, sep ="\t",stringsAsFactors = FALSE)
maffile[, "mutRates"] <- maffile[, "Freq"] / 566
write.table(maffile,paste0(output_data_dir,"mut_allfre_LUAD.txt"),sep="\t", quote=F,row.names = FALSE)

#
input_data_dir = "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/LUSC/"
output_data_dir <- "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/LUSC/"
File <-paste0(input_data_dir,"LUSC.allaa_fre.txt")
maffile <-read.delim(file= File, header =TRUE, sep ="\t",stringsAsFactors = FALSE)
maffile[, "mutRates"] <- maffile[, "Freq"] / 484
write.table(maffile,paste0(output_data_dir,"mut_allfre_LUSC.txt"),sep="\t", quote=F,row.names = FALSE)

#
input_data_dir = "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/OV/"
output_data_dir <- "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/OV/"
File <-paste0(input_data_dir,"OV.allaa_fre.txt")
maffile <-read.delim(file= File, header =TRUE, sep ="\t",stringsAsFactors = FALSE)
maffile[, "mutRates"] <- maffile[, "Freq"] / 523
write.table(maffile,paste0(output_data_dir,"mut_allfre_OV.txt"),sep="\t", quote=F,row.names = FALSE)

#
input_data_dir = "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/PAAD/"
output_data_dir <- "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/PAAD/"
File <-paste0(input_data_dir,"PAAD.allaa_fre.txt")
maffile <-read.delim(file= File, header =TRUE, sep ="\t",stringsAsFactors = FALSE)
maffile[, "mutRates"] <- maffile[, "Freq"] / 179
write.table(maffile,paste0(output_data_dir,"mut_allfre_PAAD.txt"),sep="\t", quote=F,row.names = FALSE)

#
input_data_dir = "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/READ/"
output_data_dir <- "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/READ/"
File <-paste0(input_data_dir,"READ.allaa_fre.txt")
maffile <-read.delim(file= File, header =TRUE, sep ="\t",stringsAsFactors = FALSE)
maffile[, "mutRates"] <- maffile[, "Freq"] / 140
write.table(maffile,paste0(output_data_dir,"mut_allfre_READ.txt"),sep="\t", quote=F,row.names = FALSE)

#
input_data_dir = "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/UCEC/"
output_data_dir <- "/Users/cyu/Documents/work/mc3data/topmutations/pancan/specific_cancertype/UCEC/"
File <-paste0(input_data_dir,"UCEC.allaa_fre.txt")
maffile <-read.delim(file= File, header =TRUE, sep ="\t",stringsAsFactors = FALSE)
maffile[, "mutRates"] <- maffile[, "Freq"] / 517
write.table(maffile,paste0(output_data_dir,"mut_allfre_UCEC.txt"),sep="\t", quote=F,row.names = FALSE)
