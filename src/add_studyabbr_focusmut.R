### somatic mutation data was doCurrently we mainly focus on 9 cancer type
##1 CHOL	Cholangiocarcinoma
##2 COAD Colon adenocarcinoma
##3 Colorectal adenocarcinoma
##4 Endometrial Cancer
##5 LUAD Lung adenocarcinoma (LUAD)
##6 LUSC Lung squamous cell carcinoma (LUSC)
##7 OV Ovarian Cancer
##8 PAAD Pancreatic Cancer
##9 READ Rectum adenocarcinoma
#' MC3_to_VCF
#'
#' MC3_to_VCF
#'
#'
#' @export
#' @examples
#' MC3_to_VCF()

#' @import dplyr
#' @import data.table
#' @import plyr

TSS2Study_DF <- read.table("/Users/cyu/Documents/work/hotspots/meta/TSS2studyAbbr.txt", header=TRUE,sep="\t",stringsAsFactors = FALSE)

MC3_DF <- read.table('/Users/cyu/Documents/work/hotspots/raw_data/focus_mut/select_cols_maf.txt', header=TRUE,sep="\t")
MC3_DF$TSS.Code <-gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",gsub("TCGA-","",MC3_DF$Tumor_Sample_Barcode))
head(MC3_DF$TSS.Code)
MC3_DF <- join(MC3_DF,TSS2Study_DF, by="TSS.Code")

#Save mc3-file in separate folders as VCF-files
write.table(MC3_DF,'/Users/cyu/Documents/work/hotspots/raw_data/mc3.studyAbbr.txt',sep="\t",quote=F,row.names = FALSE)

MC3_DF_CancerType <- split (MC3_DF, f=MC3_DF$Study.Abbreviation,drop=TRUE)


vcffile=vector()
subtype=vector()
sample=vector()


dir_name = "/Users/cyu/Documents/work/hotspots/raw_data/cancertype/"
#Start of loops --------------------------------------------

for (cancerType in MC3_DF_CancerType){ # Loop over all cancer types
  
  #Add a catch that prevents closing loop
  if (dim(cancerType)[1]  == 0) {
    
    next()
    
  }
  
  directory_name <- paste(dir_name,"/",unique(cancerType$Study.Abbreviation),sep="") 
  
  #crate folder if it doesnt exist
  ifelse(!dir.exists(file.path(directory_name)), dir.create(file.path(directory_name)), FALSE)
  print (directory_name)
  setwd(directory_name)
  
  write.table(cancerType, file="maf.txt", row.names=FALSE, sep="\t", quote=FALSE)
  counts<- ddply(cancerType, .(cancerType[,1], cancerType[,13]),nrow)
  names(counts) <- c("Hugo_Symbol", "amino_acid_change","Freq")
  all.AAchange <- counts[order(counts$Freq,decreasing = TRUE),]
  write.table(all.AAchange,file="aachange_Freq.txt",sep="\t", quote=F,row.names = FALSE)
  
  counts_tp53<-all.AAchange[(all.AAchange$Hugo_Symbol=="TP53"),]
  counts_kras<-all.AAchange[(all.AAchange$Hugo_Symbol=="KRAS"),]
  write.table(counts_tp53,file="aachange_tp53_Freq.txt",sep="\t", quote=F,row.names = FALSE)
  write.table(counts_kras,file="aachange_kras_Freq.txt",sep="\t", quote=F,row.names = FALSE)
  
}




counts<- ddply(MC3_DF, .(MC3_DF[,1], MC3_DF[,13], MC3_DF[,15]),nrow)
names(counts) <- c("Hugo_Symbol", "amino_acid_change", "Study.breviation","Freq")
all.AAchange <- counts[order(counts$Freq,decreasing = TRUE),]
write.table(all.AAchange,paste0(dir_name,"aachange_Freq.txt"),sep="\t", quote=F,row.names = FALSE)
counts_tp53<-all.AAchange[(all.AAchange$Hugo_Symbol=="TP53"),]
counts_kras<-all.AAchange[(all.AAchange$Hugo_Symbol=="KRAS"),]

write.table(counts_tp53,paste0(dir_name,"aachange_tp53_Freq.txt"),sep="\t", quote=F,row.names = FALSE)
write.table(counts_kras,paste0(dir_name,"aachange_kras_Freq.txt"),sep="\t", quote=F,row.names = FALSE)




##to individual vcf
#Start of loops --------------------------------------------

for (cancerType in MC3_DF_CancerType){ # Loop over all cancer types
  
  #Add a catch that prevents closing loop
  if (dim(cancerType)[1]  == 0) {
    
    next()
    
  }
  
  
  
  directory_name <- paste(dir_name,"/",unique(cancerType$Study.Abbreviation),sep="") 
  
  #crate folder if it doesnt exist
  ifelse(!dir.exists(file.path(directory_name)), dir.create(file.path(directory_name)), FALSE)
  print (directory_name)
  setwd(directory_name)
  
  
  
  
  #Save cohort level vcf file
  cohort_vcf <-  matrix(".", nrow=nrow(cancerType), ncol = 10);
  #sample_id <- unique(Sample$Tumor_Sample_Barcode);
  columns=c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT","?" )
  colnames(cohort_vcf) <- columns
  
  cohort_vcf[,1]=as.character(cancerType$Chromosome)
  cohort_vcf[,2]=as.character(cancerType$Start_Position)
  cohort_vcf[,4]=as.character(cancerType$Reference_Allele)
  cohort_vcf[,5]=as.character(cancerType$Tumor_Seq_Allele2)
  cohort_vcf[,9]="GT"
  cohort_vcf[,10]="1/0"
  
  cohort_vcf[,3]=as.character(cancerType$Tumor_Sample_Barcode)
  
  write.table(cohort_vcf, file="Cohort.vcf", row.names=FALSE, sep="\t", quote=FALSE)
  
  
  
  #Split cancer type into individual samples
  MC3_DF_Sample <- split (cancerType, f= cancerType$Tumor_Sample_Barcode,drop = TRUE) # not sure if this works
  
  
  
  
  
  
  for (Sample in MC3_DF_Sample){ #Loop over each sample in respective cancer type
    
    
    if (dim(Sample)[1]  == 0) { ########### Why does CancerType get split into empty samples??????????!!!!!!!!
      #Seem to split MC3_DF_Sample but still save old levels  # try drop = TRUE in split command
      
      next()
      
    } 
    
    
    
    sample_id <- unique(Sample$Tumor_Sample_Barcode);
    vcfdata <- matrix(".", nrow=nrow(Sample), ncol = 10);
    columns=c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",sample_id )
    colnames(vcfdata) <- columns
    
    vcfdata[,1]=as.character(Sample$Chromosome)
    vcfdata[,2]=as.character(Sample$Start_Position)
    vcfdata[,4]=as.character(Sample$Reference_Allele)
    vcfdata[,5]=as.character(Sample$Tumor_Seq_Allele2)
    vcfdata[,9]="GT"
    vcfdata[,10]="1/0"
    #vcfdata[,3]=as.character(sample_id)
    outfile <- paste(sample_id,".vcf", sep="")
    #outfile = paste(directory_name,outfile, sep="/")
    write.table(vcfdata, file=outfile, row.names=FALSE, sep="\t", quote=FALSE)
    
    subtype <- c(subtype,as.character(unique(cancerType$Study.Abbreviation)))
    sample <- c(sample,as.character(sample_id))
    vcffile <- c(vcffile, paste (unique(cancerType$Study.Abbreviation),outfile,sep = "/") )
    
  }
  
  
}



cohort <- data.frame("vcf"=vcffile, "subtype"=subtype, "sample"=sample) 
setwd(main_wd)
write.table(cohort,file="Cohort.txt")


