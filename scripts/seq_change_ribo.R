setwd("/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/SAB_pipeline/Illumina")


library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rempsyc)
library(qwraps2)

#this script takes output from multiQC reports of before and after rRNA removal via
#Ribodetector to see how sequences were changed. 

#read in data
multiqc <- read.csv("seq_change_ribo.csv")[-63, ]


multiqc <- mutate(multiqc, sample = NA, treatment = NA, rep = NA) #create columns to make grouping easier later

#fill in sample column:
#loop through each row the data set, at each row, i, search for strings starting with
#04, 06, 08, or 13 and use that number to fill in the respective sample column

for(i in 1:length(multiqc$sample)){
      if (grepl("04.*", multiqc$Sample.Name[i]) == TRUE){ #when this statement is true,
  multiqc$sample[i] <- 4        #we 'write in' 04 for the sample column
}else if (grepl("06.*", multiqc$Sample.Name[i]) == TRUE){ 
  multiqc$sample[i] <- 6
}else if (grepl("08.*", multiqc$Sample.Name[i]) == TRUE){ 
  multiqc$sample[i] <- 8
}else if (grepl("13.*", multiqc$Sample.Name[i]) == TRUE){
  multiqc$sample[i] <- 13
}else{ print("no match in row")
    }
}
  
#fill in treatment column 

for(i in 1:length(multiqc$treatment)){
  if (grepl(".*add_back.*", multiqc$Sample.Name[i]) == TRUE){ 
    multiqc$treatment[i] <-  "Fe_add_back"       
  }else if (grepl(".*pFe19.*", multiqc$Sample.Name[i]) == TRUE){ 
    multiqc$treatment[i] <- "High_Fe" 
  }else if (grepl(".*pFe21_9.*", multiqc$Sample.Name[i]) == TRUE){ 
    multiqc$treatment[i] <- "Low_Fe"
  }else if (grepl(".*oFe21_9.*", multiqc$Sample.Name[i]) == TRUE){ 
    multiqc$treatment[i] <- "Low_Fe" #correcting for typo 
  }else{ 
    print("no match in row")
}}

#fill in rep column

for(i in 1:length(multiqc$rep)){
  if (grepl(".*A.*", multiqc$Sample.Name[i]) == TRUE){ 
    multiqc$rep[i] <-  "A"       
  }else if (grepl(".*B.*", multiqc$Sample.Name[i]) == TRUE){ 
    multiqc$rep[i] <- "B" 
  }else if (grepl(".*C.*", multiqc$Sample.Name[i]) == TRUE){ 
    multiqc$rep[i] <- "C"
  }else{ 
    print("no match in row")
  }}

#group data by sample and treatment 

multiqc <- group_by(multiqc, sample, treatment)


perc_seq_rem <- ggbarplot(
  multiqc, 
  x="sample", 
  y="seq_removed_perc",
  add = 'mean_sd',
  fill="treatment", 
  position=position_dodge(0.8),
  ylab=" % sequences removed",
  title= "Sequences removed by Ribodetector") +
  theme(plot.title=element_text(hjust=0.5),axis.title.x=element_blank(),
        aspect.ratio = 1, legend.position = "bottom")

ggsave("Seq_removed_ribo.pdf", plot = perc_seq_rem, width = 15, height = 10)


#make a nice table output
seq_changes <- summarize(multiqc, 
                         Mean = signif(mean(seq_removed_perc), digits = 3), 
                         SD = signif(sd(seq_removed_perc), digits = 3))

seq_changes <- unite(seq_changes, Mean, SD,col= 'Mean(SD)',sep=" ± ")

seq_changes <- pivot_wider(seq_changes, 
                           names_from = "treatment",
                           values_from = 'Mean(SD)')

colnames(seq_changes)<- c("Organism","Iron Add Back", "Low Iron", "High Iron")
seq_changes$Organism <- c("Outer Shelf Cylindrotheca",
                          "Inner Shelf G. oceanica",
                          "Inner Shelf Cylindrotheca",
                          "Outer Shelf E. huxleyi")


seq_changes <- nice_table(
  data = seq_changes, 
  title="Percent seqences removed by Ribodetector")


#save the tabel to word
flextable::save_as_docx(seq_changes, path = "seq_removed_ribo.docx")











