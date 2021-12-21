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
oscc = read.maf(maf = maf)

##Plotting MAF summary.
pdf(file=opt$output,width=8,height=5)
oncoplot(maf=oscc, top=20,annotationFontSize = 2, titleFontSize = 1, legendFontSize = 1)
dev.off()

###usuage:
##Rscript oncoplot.r --input /Users/chunleiyu/Work/data/vcf/CGU_01/maf_new/addbarcode/oscc.maf_all.txt -o /Users/chunleiyu/Work/data/vcf/CGU_01/maf_new/addbarcode/results_plot/oncoplt1.pdf