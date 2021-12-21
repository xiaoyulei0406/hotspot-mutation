library(maftools)
library(VariantAnnotation)
library(ComplexHeatmap)
library(Biostrings)
brca.mut<-read.delim(file = "/Users/cyu/Downloads/BRCA.Mutation_Packager_Calls.Level_3.2016012800.0.0/gdac.BRCA.Mutation_Packager_Calls.Level_3.2016.maf.txt",header = T,sep ="\t", stringsAsFactors = FALSE)
brca.mut_temp<-brca.mut[!(brca.mut$Variant_Classification=="Silent"),]
brca.mut_temp<-brca.mut_temp[!(brca.mut_temp$Variant_Classification=="RNA"),]
brca.mut_temp<-brca.mut_temp[!(brca.mut_temp$Variant_Classification=="Splice_Site"),]
brca.mut_temp<-brca.mut_temp[!(brca.mut_temp$Variant_Classification=="Nonsense_Mutation"),]
table(brca.mut_temp$Variant_Classification)
brca.mut<-brca.mut_temp
table(brca.mut$Variant_Classification)

write.table(brca.mut_temp,file = "/Users/cyu/Downloads/BRCA.Mutation_Packager_Calls.Level_3.2016012800.0.0/mut.txt",sep="\t", quote=F)
brcamaf<-read.maf(maf = brca.mut)
plotmafSummary(maf=brcamaf)

library(g3viz)


# load data and change names for g3viz input
colnames(brca.mut)[which(names(brca.mut) == "Hugo_Symbol")] <- "gene.symbol"
colnames(brca.mut)[which(names(brca.mut) == "Variant_Classification")] <- "variant.class"
colnames(brca.mut)[which(names(brca.mut) == "amino_acid_change_WU")] <- "protein.change"

brca.mutation.csv <- "/Users/cyu/Downloads/BRCA.Mutation_Packager_Calls.Level_3.2016012800.0.0/mut.txt"
brca.mutation.dat <- readMAF(brca.mutation.csv,
                        gene.symbol.col = "Hugo_Symbol",
                        variant.class.col = "Variant_Classification",
                        protein.change.col = "amino_acid_change_WU",
                        sep = "\t")  # column-separator of csv file

# set up chart options
plot.options <- g3Lollipop.options(
  # Chart settings
  chart.width = 600,
  chart.type = "pie",
  chart.margin = list(left = 30, right = 20, top = 20, bottom = 30),
  chart.background = "#d3d3d3",
  transition.time = 300,
  # Lollipop track settings
  lollipop.track.height = 200,
  lollipop.track.background = "#d3d3d3",
  lollipop.pop.min.size = 1,
  lollipop.pop.max.size = 8,
  lollipop.pop.info.limit = 5.5,
  lollipop.pop.info.dy = "0.24em",
  lollipop.pop.info.color = "white",
  lollipop.line.color = "#a9A9A9",
  lollipop.line.width = 3,
  lollipop.circle.color = "#ffdead",
  lollipop.circle.width = 0.4,
  lollipop.label.ratio = 2,
  lollipop.label.min.font.size = 12,
  lollipop.color.scheme = "dark2",
  highlight.text.angle = 60,
  # Domain annotation track settings
  anno.height = 16,
  anno.margin = list(top = 0, bottom = 0),
  anno.background = "#d3d3d3",
  anno.bar.fill = "#a9a9a9",
  anno.bar.margin = list(top = 4, bottom = 4),
  domain.color.scheme = "pie5",
  domain.margin = list(top = 2, bottom = 2),
  domain.text.color = "white",
  domain.text.font = "italic 8px Serif",
  # Y-axis label
  y.axis.label = "# of TP53 gene mutations",
  axis.label.color = "#303030",
  axis.label.alignment = "end",
  axis.label.font = "italic 12px Serif",
  axis.label.dy = "-1.5em",
  y.axis.line.color = "#303030",
  y.axis.line.width = 0.5,
  y.axis.line.style = "line",
  y.max.range.ratio = 1.1,
  # Chart title settings
  title.color = "#303030",
  title.text = "TP53 gene (customized chart options)",
  title.font = "bold 12px monospace",
  title.alignment = "start",
  # Chart legend settings
  legend = TRUE,
  legend.margin = list(left=20, right = 0, top = 10, bottom = 5),
  legend.interactive = TRUE,
  legend.title = "Variant classification",
  # Brush selection tool
  brush = TRUE,
  brush.selection.background = "#F8F8FF",
  brush.selection.opacity = 0.3,
  brush.border.color = "#a9a9a9",
  brush.border.width = 1,
  brush.handler.color = "#303030",
  # tooltip and zoom
  tooltip = TRUE,
  zoom = TRUE
)
g3Lollipop(brca.mutation.dat ,
           gene.symbol = "TP53",
           protein.change.col = "amino_acid_change_WU",
           btn.style = "blue", # blue-style chart download buttons
           plot.options = plot.options,
           output.filename = "customized_plot")








# load data
mutation.csv <- system.file("extdata", "ccle.csv", package = "g3viz")

# ============================================
# read in data
#   "gene.symbol.col"    : column of gene symbol
#   "variant.class.col"  : column of variant class
#   "protein.change.col" : colum of protein change column
# ============================================
mutation.dat <- readMAF(mutation.csv,
                        gene.symbol.col = "Hugo_Symbol",
                        variant.class.col = "Variant_Classification",
                        protein.change.col = "amino_acid_change",
                        sep = ",")  # column-separator of csv file



chart.options <- g3Lollipop.theme(theme.name = "default",
                                  title.text = "PIK3CA gene (default theme)")

g3Lollipop(mutation.dat,
           gene.symbol = "PIK3CA",
           protein.change.col = c("amino_acid_change"),
           plot.options = chart.options,
           output.filename = "default_theme")




