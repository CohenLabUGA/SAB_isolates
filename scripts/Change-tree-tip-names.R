#This code will subsitute more meaningful names to the tip of phylogenetic trees made 
#in mega. The tree is constructed in mega where it is saved as a .nwk file. the names
#of the organisms are saved as a .csv file, with the first column being the current tree
#tip names and second what you want to change them to. After changing the tip names, open the
#tree in mega to view the results.


setwd( "/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/Figures/Phylogenetic_Trees")
library('ape')
library('phytools')
# taken from:
#https://www.researchgate.net/post/How_to_edit_tip_labels_in_MEGA



#### author: Jinlong Zhang <jinlongzhang01@gmail.com>
#### institution: Kadoorie Farm and Botanic Garden, Hong Kong
#### package: phylotools
#### URL: http://github.com/helixcn/phylotools
#### date: 26 MAY 2015
#### Modified: 22 AUG 2018

#### Function sub.tip.label as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

sub.taxa.label <- function(tree, dat){
  
  if(!inherits(tree,"phylo")){
    stop("the input tree is not a \"phylo\" object.")
  }
  
  if(!is.data.frame(dat)){
    stop("the input dat is not a \'data.frame\'.")
  }
  
  tree2 <- tree
  nnn <- tree$tip.label
  
  if(!nrow(dat) == length(nnn)){
    warning("Number of tip labels in the phylogenetic\n tree differs from the number of rows in the reference table.\n")
  }
  
  xxx1 <- as.character(dat[,1])
  xxx2 <- as.character(dat[,2])
  
  if(!all(xxx1 %in% nnn)){
    unsub.dat <- xxx1[!xxx1 %in% nnn]
    cat("The following names in the reference data.frame \n can not be found in the phylogeny:\n", unsub.dat, "\n")
  }
  
  if(!all(nnn %in% xxx1)){
    unsub.tree <- nnn[!nnn %in% xxx1]
    cat("The following tip labels in the phylogenetic tree \ncan not be found in the reference data.frame:\n", unsub.tree, "\n")
  }
  
  label <- c()
  for(i in 1:length(dat[,1])){
    for(j in 1:length(nnn)){
      if(nnn[j] == xxx1[i]){
        label[j] <- xxx2[i]
      }
    }
  }
  
  tree2$tip.label <- label
  return(tree2)
}

tree1<-read.tree('cylindro_cox1.nwk')

#tips <- tree$tip.label
organism.names <- read.csv('cylindro_cox1.csv')

cox1_tree <- sub.taxa.label(tree1,organism.names)

write.tree(cox1_tree,file="Cylindro-cox1-alignment.nwk")



