#### Cell Sizes ####
## Lucy Quirk 
## UGA Master's 
--------------------------------------------
## This script will run ANOVA on cell size data and create a df which will be 
## used for plotting in a script called physio_plot.R 

--------------------------------------------
## set up environment and read in data
--------------------------------------------
source('functions.R')
pkgs <- c('tidyverse','rstatix')
checkAndLoadPackages(pkgs)

path_in <- c('../rawData/')
path_out <- c('../output')


size <- read.csv(paste(path_in,'exp_cell_size.csv', sep='')) 

--------------------------------------------
#### clean up data ####
--------------------------------------------

size=pivot_longer(size, !Culture, names_to = c('Treatment', 'rep'), 
                  names_pattern = '(.*\\..*\\.)(A|B|C)', values_to = 'size')

#average by rep then by treatment
size <- size %>% group_by(Culture,rep, Treatment) %>%
  summarise(Size=mean(size), sd=sd(size))

#add empty line for missing high iron treatment in C. clost UGA4
empty.04 <- data.frame("Culture"=4, 'rep'='A', "Treatment"="High Iron", 
                       'Size'=0, 'sd'=0, 'Shelf'='Inner Shelf')
size<- bind_rows(size,empty.04)

# fix organism names and add shelf location
size$Culture <- gsub( "4", "C. closterium UGA4", size$Culture)
size$Culture <- gsub( "8", "C. closterium UGA8", size$Culture)
size$Culture <- gsub( "6", "G. oceanica", size$Culture)
size$Culture <- gsub( "13", "G. huxleyi", size$Culture)

for (i in 1:nrow(size)) {
  if (size$Culture[i] == "G. oceanica"|size$Culture[i] == "C. closterium UGA8"){
    size$Shelf[i] = "Inner shelf"} 
  else{size$Shelf[i] = "Outer shelf"}
}


size$Treatment <- str_replace_all(size$Treatment, 
                                  c('High.Iron.'='High Iron', 'Low.Iron.'='Low Iron',
                                    'Iron.Amendment.'='Iron Amendment')) 
size$Treatment <- str_remove(size$Treatment, '\\.')

size <- na.omit(size)

--------------------------------------------
#### Run ANOVA ####
--------------------------------------------
aov.size <- function(organism){
  df <- size %>% filter(Culture==organism) 
  
  #1. run anova using aov
  anova.org <- aov(Size~as.factor(Treatment), data=df)
  print(summary(anova.org))
  
  #2. create summary table of organism for plotting later
  summary.org <- df %>% 
    group_by(Treatment, Culture) %>%
    summarize(size = mean(Size),sd=sd(Size))%>%
    arrange(desc(size)) ##!!! this is important so letters from cld match right later
  
  #3. run tukey's test using anova output
  tukey.org <- TukeyHSD(anova.org)
  print(tukey.org)
  
  #4. compact letters asociated with tukey's output
  cld.org <- multcompLetters4(anova.org,tukey.org)
  print(cld.org)
  cld.org <- as.data.frame.list(cld.org$`as.factor(Treatment)`)
  
  #5. add the results to summary table in a new column
  # this will make row binding later easier
  summary.org$Tukey <- cld.org$Letters
  summary.org
}


aov.size.4 <- aov.size("C. closterium 4")
aov.size.8 <- aov.size("C. closterium 8")
aov.size.6 <- aov.size("G. oceanica")
aov.size.13 <- aov.size("G. huxleyi")

#6. once this is done for all organisms, bind rows of data frames

size_all <- rbind(aov.size.8, aov.size.4, aov.size.6, aov.size.13)
size_all <- group_by(size_all,Culture, Treatment) 


for (i in 1:nrow(size_all)) {
  if (size_all$Culture[i] == "G. oceanica"|size_all$Culture[i] == "G. huxleyi"){
    size_all$taxa[i] = "Coccolithophore"} 
  else{size_all$taxa[i] = "Diatom"}
}

size_all$Culture <- factor(size_all$Culture,
                          levels = c("G. oceanica","G. huxleyi",
                                     "C. closterium 8", "C. closterium 4"))

size_all$Treatment <- factor(size_all$Treatment,
                            levels = c("High Iron","Low Iron","Iron Amendment"))
