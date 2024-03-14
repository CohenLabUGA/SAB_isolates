# ----------Sourmash compare plotting---------- 
# set wd and load libraries

setwd("/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/SAB_pipeline/Sourmash")

path.fig <- "/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/Figures/bioinformatics/"

library(tidyverse)
library(Rtsne)
library(gplots)
library(circlize)
library(viridis)
#if (!require("BiocManager", quietly = TRUE)) #----NOT NEEDED EVERY TIME
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
# https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
  
# read in sourmash compare matrix

comp4 <- read_csv("04compare.csv")
labs4 <- read_tsv("04compare.labels.txt", col_names = FALSE)
labs4$short <- c("04 A","04 B","04 C")

comp6 <- read_csv("06compare.csv")
labs6 <- read_tsv("06compare.labels.txt", col_names = FALSE)
labs6$short <- c("06 A","06 B","06 C")

comp8 <- read_csv("08compare.csv")
labs8 <- read_tsv("08compare.labels.txt", col_names = FALSE)
labs8$short <- c("08 A","08 B","08 C")

comp13 <- read_csv("13compare.csv")
labs13 <- read_tsv("13compare.labels.txt", col_names = FALSE)
labs13$short <- c("13 A","13 B","13 C")

diatoms_comp <- read_csv("diatoms_compare.csv")
diatoms_labels <- read_tsv("diatoms_compare.labels.txt", col_names = FALSE)
diatoms_labels$short <- c("04 A","04 B","04 C","08 A","08 B","08 C")

coccoliths_comp <- read_csv("coccolithophores_compare.csv")
coccoliths_labels <- read_tsv("coccolithophores_compare.labels.txt", col_names = FALSE)
coccoliths_labels$short <- c("06 A","06 B","06 C","13 A","13 B","13 C")

all_taxa_comp <- read_csv("all_taxa_compare.csv")
all_taxa_labels <- read_tsv("all_taxa_compare.labels.txt", col_names = FALSE)
all_taxa_labels$short <- c("04 A","04 B","04 C",
                           "08 A","08 B","08 C",
                           "06 A","06 B","06 C",
                           "13 A","13 B","13 C")

#label the columns and turn into matrix

colnames(comp4) <- labs4$short
mat4 <- as.matrix(comp4)
rownames(mat4) <- colnames(mat4)

colnames(comp6) <- labs6$short
mat6 <- as.matrix(comp6)
rownames(mat6) <- colnames(mat6)

colnames(comp8) <- labs8$short
mat8 <- as.matrix(comp8)
rownames(mat8) <- colnames(mat8)

colnames(comp13) <- labs13$short
mat13 <- as.matrix(comp13)
rownames(mat13) <- colnames(mat13)

colnames(diatoms_comp) <- diatoms_labels$short
diatoms_mat <- as.matrix(diatoms_comp)
rownames(diatoms_mat) <- colnames(diatoms_mat)

colnames(coccoliths_comp) <- coccoliths_labels$short
coccoliths_mat <- as.matrix(coccoliths_comp)
rownames(coccoliths_mat) <- colnames(coccoliths_mat)

colnames(all_taxa_comp) <- all_taxa_labels$short
all_taxa_mat <- as.matrix(all_taxa_comp)
rownames(all_taxa_mat) <- colnames(all_taxa_mat)

####Sample-Sample Matrix####
d_s_c <- read_csv("./by_sample/diatoms_by_sample_compare.csv") #d_s_c diatoms by sample comparison matrix
c_s_c <- read_csv("./by_sample/coccoliths_by_sample_compare.csv")
a_s_c <- read_csv("./by_sample/all_samples_compare.csv")


shorten_cols <- function(df){
  # extract sample name from full path to sample character string
  # correct one sample with mislabels pFe as oFe 
  # replace column names with shortened and corrected names
 short_name <-  str_extract(colnames(df),
              "(?<=/)[[:digit:]][[:alnum:]]*.[[:alnum:]]*(?=/)")
 short_name <- str_replace(short_name, pattern="oFe", replacement="pFe")
 colnames(df) <- short_name
 df
 }

a_s_c <- shorten_cols(a_s_c)

make_mat <- function(df){
  mat <- as.matrix(df)
  rownames(mat) <- colnames(mat)
  mat
}

mat_asc <- make_mat(a_s_c)

#create heatmap annotations
#make a dataframe which maps to the data kind of like metadata with the colnames as the
#rownames and a column for each annotation I want
make.meta <- function(df){
  hh <- colnames(df)
  ha <- as.data.frame(
    str_match(hh,
            "([[:digit:]]{2})([[:alpha:]]{3}[[:digit:]]*[[:punct:]]*[[:alnum:]]*)([ABC])"))
  colnames(ha) <- c('Sample',"Organism","Treatment","Replicate")
  rownames(ha) <- hh
  ha <- ha[,-1]
  ha
}


all_meta <- make.meta(a_s_c)
all_meta$Organism <- all_meta$Organism %>%
  str_replace('04','C. closterium 4') %>% 
  str_replace('08','C. closterium 8') %>%
  str_replace( '06','G. oceanica') %>%
  str_replace( '13','G. huxleyi')

all_meta$Treatment <- all_meta$Treatment %>% 
  str_replace('pFe19', 'High Iron') %>%
  str_replace('pFe21_9', 'Low Iron') %>%
  str_replace('add_back', 'Iron Amendment') 


anno_col=list(
  #col - list of colors which contain color mapping to columns in df
  Organism=
    c("C. closterium 4"="darkred", "C. closterium 8"="lightpink", 
      "G. huxleyi"="darkgreen", "G. oceanica"="lightgreen"),
  Treatment=
    c("High Iron"="black", "Iron Amendment"="grey", "Low Iron"="white"))

make.ha <- function(df){
  #create heatmap annotation object for complex heatmaps using metadata 
  #df created in make.meta
  HeatmapAnnotation(show_legend = F,
    Organism=df$Organism,
  Treatment=df$Treatment, 
  col=anno_col)
}

all_ha <- make.ha(all_meta)
va <- rowAnnotation(Treatment=all_meta$Treatment,
              Organism=all_meta$Organism,
              col=anno_col, show_annotation_name=F)
#combine heatmap and annotation

# heatmap clustered row and columns with all annotations 

simple_heat <- function(mat, annotation_df, fig.name, text=FALSE){
  # plot heatmap with or without writing correlation values in the cells
  # and saving this plot in a defined path as a rastered png file
  
  png(filename=paste0(path.fig, fig.name, ".png", sep=""), width=800, height=800)
  ht <- Heatmap(mat, col=magma(100), right_annotation = va,
                top_annotation = annotation_df, 
                cluster_rows = TRUE, cluster_columns = TRUE,
        show_column_dend = FALSE, row_dend_side = "right",
        show_column_names = F, show_row_names = F,
        heatmap_legend_param = list(title="Correlation"),
        layer_fun = function(j,i,x,y,width,height,fill){
          if (text==TRUE){
            grid.text(sprintf("%.1f", pindex(mat, i, j)), 
                      x, y, gp = gpar(fontsize=10))}},
        use_raster = TRUE, raster_device="png")
  draw(ht)
  dev.off()
}



simple_heat(mat_asc, all_ha, "all_by_sample")
simple_heat(mat_asc, all_ha, "all_by_sample_numbered", text=TRUE)

# more traditional correlation matrix
corelation <- function(mat, fig.name){
  col_fun = colorRamp2(c(0,0.5,1),c("white","purple","red"))
  
  tiff(filename=paste0(path.fig, fig.name, ".png", sep=""), width=800, height=800)
  ht <- Heatmap(mat, col=viridis(100), #rect_gp=gpar(type='none'), 
        cluster_rows = FALSE, cluster_columns = FALSE,
        cell_fun = function(j,i,x,y,width,height,fill){
          if(i == j) { 
            # in cells where row number = column number, write the sample name
            # given in the rowname of the matrix
            grid.text(rownames(mat)[i], x=x, y=y)
          } #else if (i > j) {
            # cells below diagonal line fill with color
            #grid.circle(
            #  x=x, y=y, 
              # r=(abs(mat[i,j])/1*min(unit.c(width, height))),
            #  r=(1/2*min(unit.c(width, height))),
            #  gp=gpar(fill = viridis(100)))#col_fun(mat[i, j]),col=NA))
            
           else if (i < j){
            # cells above diagonal line, write correlation number
             grid.null()
           grid.text(sprintf("%.1f", mat[i,j]),x,y,gp=gpar(fontsize=10))
        }}, show_row_names = FALSE, show_column_names = FALSE)
  draw(ht)
  dev.off()
}



corelation(mat_asc, "all_by_sample_corr")



comp_heatmap <- function(matrix){
  ff = pheatmap(matrix, scale='none', 
         cluster_cols=TRUE,
         cluster_rows=FALSE,
         display_numbers=FALSE,
         number_color='black',
         fontsize_number = 8,
         color = hcl.colors(50, "BluYl"))
  print(ff)
         #filename = paste(x,"_sourmash.tiff", sep=""))
}

comp_heatmap(mat_dsc)
comp_heatmap(mat_csc)
comp_heatmap(mat_asc)

comp_heatmap(mat4, "4")
comp_heatmap(mat6, "6")
comp_heatmap(mat8, "8")
comp_heatmap(mat13, "13")

comp_heatmap(diatoms_mat, "diatoms")
comp_heatmap(coccoliths_mat, "coccoliths")
comp_heatmap(all_taxa_mat, "all_taxa")

comp_heatmap_simple <- function(matrix, x){
  pheatmap(matrix, scale='none', 
           cluster_cols=FALSE, 
           color = hcl.colors(50, "BluYl"),
           filename = paste(x,"_sourmash_simple.tiff", sep=""))
}

comp_heatmap_simple(mat4, "4")
comp_heatmap_simple(mat6, "6")
comp_heatmap_simple(mat8, "8")
comp_heatmap_simple(mat13, "13")

comp_heatmap_simple(diatoms_mat, "diatoms")
comp_heatmap_simple(coccoliths_mat, "coccoliths")
comp_heatmap_simple(all_taxa_mat, "all_taxa")

# ------Make a clustered heatmap---------

hc.rows <- hclust(dist(diatoms_mat))
hc.cols <- hclust(dist(t(diatoms_mat)))
heatmap(diatoms_mat[cutree(hc.rows,h=2)==1,], Colv=as.dendrogram((hc.cols), scale='none'))





transcript <- read.csv('transcriptome_compare.csv')
tc <- transcript %>% as.matrix()
rownames(tc) <- colnames(tc)

labels.new <- structure(c('C. closterium UGA4','G. oceanica',
                          'C. closterium UGA8','G. huxleyi'), names=colnames(tc))
transcriptome_compare <- Heatmap(tc, 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", tc[i, j]), x, y, 
                    gp = gpar(fontsize = 30))},
        column_names_rot=25, column_names_side='top', column_names_centered = T,
        row_labels=labels.new[rownames(tc)],
        column_labels = labels.new[colnames(tc)], column_names_gp=gpar(fontsize=20),
       row_names_gp=gpar(fontsize=20), heatmap_legend_param = gpar(title='Similarity\ncorrelation'))

png(filename=paste0(path.fig, "transcriptomes_compare.png", sep=""), width=800, height=700)
draw(transcriptome_compare)
dev.off()










