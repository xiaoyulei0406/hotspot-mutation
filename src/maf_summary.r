#source("https://bioconductor.org/biocLite.R")
#biocLite("maftools")
#biocLite("ComplexHeatmap")
#biocLite("VariantAnnotation")
#biocLite("Biostrings")
library(ComplexHeatmap)
library(optparse)
library(VariantAnnotation)
library(Biostrings)
library(maftools)
suppressMessages(suppressWarnings(require(optparse, quietly=TRUE)))

option_list = list(
  make_option(c("-i", "--input"),  type = "character", default = NULL,  help = "Input maf file"),
  make_option(c("-o", "--output"), type = "character", default = "vars.maflite",  help = "Output mutation analysis plot")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)
if (is.null(opt$input))
  stop(print_help(parseobj))

print('...loading maftools package')
suppressMessages(suppressWarnings(require(maftools, quietly=TRUE)))
print('...reading maf')

maf = read.delim(opt$input, header =TRUE, sep ="\t", stringsAsFactors=FALSE)
#maf = read.delim("/Users/chunleiyu/Work/data/vcf/CGU_01/maf_new/addbarcode/oscc.maf_all1.txt", header =TRUE, sep ="\t")
#clinical=read.table("/Users/chunleiyu/Work/data/vcf/CGU_01/maf_new/addbarcode/clinical.txt", header =TRUE, sep ="\t")
oscc = read.maf(maf = maf)

##Plotting MAF summary.
pdf(file=opt$output,width=8,height=5)
plotmafSummary(maf = oscc, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

###usuage:
##Rscript maf_summary.r --input /Users/chunleiyu/Work/data/vcf/CGU_01/maf_new/addbarcode/oscc.maf_all.txt -o /Users/chunleiyu/Work/data/vcf/CGU_01/maf_new/addbarcode/results_plot/mafSummary1.pdf