#### November 9 2022 Lucy Quirk ###
### This scripts combines physiological data from growth experiments into one 
### plot and runs statistical tests. Functions included in the functions.R
### script are written out for clarification, but after source('functions.R')
### these are already available. Significance and individual plots are created
### one at a time for each parameter: cell size, growth rate, chl a concentration, 
### fv/fm, tQa, sigma, and NPQ. 

### Output of the script are summary tables which are called in from the ouput/
### folder to make pretty tables in the make_tables.R script.

#------Prepare environment-----
source('functions.R')

checkAndLoadPackages('tidyverse','patchwork','rstatix','ggpubr','ggmisc','multcompView',
                     'shadowtext')
#-------Cell Sizes--------
size <- read.csv('../rawData/exp_cell_size.csv') 
size=pivot_longer(size, !Organism, names_to = c('Treatment', 'rep'), 
             names_pattern = '(.*\\..*\\.)(A|B|C)', values_to = 'size')

#average by rep then by treatment
size <- size %>% group_by(Organism,rep, Treatment) %>%
  summarise(Size=mean(size), sd=sd(size))

empty.04 <- data.frame("Organism"=4, 'rep'='A', "Treatment"="High Iron", 
                       'Size'=0, 'sd'=0, 'Shelf'='Inner Shelf')

size<- bind_rows(size,empty.04)

size$Organism <- gsub( "4", "C. closterium UGA4", size$Organism)
size$Organism <- gsub( "8", "C. closterium UGA8", size$Organism)
size$Organism <- gsub( "6", "G. oceanica", size$Organism)
size$Organism <- gsub( "13", "G. huxleyi", size$Organism)

for (i in 1:nrow(size)) {
  if (size$Organism[i] == "G. oceanica"|size$Organism[i] == "C. closterium UGA8"){
    size$Shelf[i] = "Inner shelf"} 
  else{size$Shelf[i] = "Outer shelf"}
}


size$Treatment <- str_replace_all(size$Treatment, 
                                  c('High.Iron.'='High Iron', 'Low.Iron.'='Low Iron',
                              'Iron.Amendment.'='Iron Amendment')) 
size$Treatment <- str_remove(size$Treatment, '\\.')
size <- na.omit(size)

## The anova.fire function is also located in the functions.R file. This 
## function takes two variables, an organism and a variable of interest. 
## The fire data is subset by organism and variable of interest. Then a
## one way ANOVA and Tukey's test is run on the data to find significant 
## diferences. Using the `multicompView` package multicompletters4 function
## creates letters to represent significance. 
anova.flex <- function(df,organism, var){
  aov.df <- df %>% filter(Organism==organism) 
  aov.df <- data.frame(Treatment=aov.df$Treatment, 
                       Param = aov.df[[var]], # a work around to change column of interest
                       Organism= aov.df$Organism)
  
  # 1. run anova using aov
  anova.org <- aov(Param~as.factor(Treatment), data=aov.df)
  print(summary(anova.org))
  # 2. create summary table for plotting later
  summary.org <- aov.df %>% 
    group_by(Treatment) %>%
    summarize('av.Param' = mean(Param),'sd'=sd(Param))%>%
    arrange(desc('av.Param')) ##!!! this is important so letters from cld map 
  ## to the right treatment
  
  # 3. run Tukey's test using anova output
  tukey.org <- TukeyHSD(anova.org)
  print(tukey.org)
  
  # 4. compact letters associated with Tukey's output
  cld.org <- multcompLetters4(anova.org,tukey.org)
  print(cld.org)
  # pull out letters
  cld.org <- as.data.frame.list(cld.org$`as.factor(Treatment)`)
  
  #5. add the results to summary table in a new column
  summary.org$Tukey <- cld.org$Letters
  summary.org <- mutate(summary.org, Organism=organism)
  summary.org
}

aov.size.4 <- anova.flex(df = size, organism = 'C. closterium UGA4', var = 'Size')
aov.size.8 <- anova.flex(size,"C. closterium UGA8", "Size")
aov.size.6 <- anova.flex(size,"G. oceanica", "Size")
aov.size.13 <- anova.flex(size,"G. huxleyi", "Size")

#6. once this is done for all organisms, bind rows of data frames

size.sig <- rbind(aov.size.8, aov.size.4, aov.size.6, aov.size.13)
size.sig <- group_by(size.sig,Organism, Treatment) 


for (i in 1:nrow(size.sig)) {
  if (size.sig$Organism[i] == "G. oceanica"|size.sig$Organism[i] == "G. huxleyi"){
    size.sig$taxa[i] = "Coccolithophore"} 
  else{size.sig$taxa[i] = "Diatom"}
}

size.sig$Organism <- factor(size.sig$Organism,
                           levels = c("G. oceanica","G. huxleyi",
                                      "C. closterium UGA8", "C. closterium UGA4"))

size.sig$Treatment <- factor(size.sig$Treatment,
                            levels = c("High Iron","Low Iron","Iron Amendment"))

#####-----Cell Size Plot------
barplot.aov <- function(df.sig, hist.parameter=-1, y_axis){
  df.sig <- df.sig %>% group_by(Organism, Treatment)
  plot = ggplot(df.sig, aes(Organism ,av.Param, fill=Treatment)) +
    geom_bar(color="black",stat='identity', position='dodge2') + 
    geom_errorbar(aes(ymin=av.Param-sd, ymax=av.Param+sd),width=0.2,
                  position = position_dodge(width=0.9))+
    geom_text(aes(label=Tukey, y=av.Param+sd,vjust=-1.8,group=Treatment),
              show.legend=FALSE,position=position_dodge(width=0.9)) +
    facet_grid(~taxa, scales="free_x", space="free_x", switch="x") +
    theme_pubr()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(face="italic"),legend.position='none',
          strip.placement = 'outside',
          strip.text = element_text(size=15,face='bold'),
          strip.background = element_rect(fill='white', color='white')) + 
    scale_x_discrete(labels=c(
      "C. closterium UGA8" = "C. closterium UGA8\n(Inner)",
      "G. oceanica" = "G. oceanica\n(Inner)",
      "C. closterium UGA4" = "C. closterium UGA4\n(Outer)",
      "G. huxleyi" = "G. huxleyi\n(Outer)")) +
    scale_fill_manual(labels=c('High Iron','Low Iron','Iron Amendment'),
                      values=c("black","white","grey"))+
    labs(y=y_axis, fill='Treatment') + 
    scale_y_continuous(expand=expansion(mul=c(0,0.2)))
  if( hist.parameter !=-1) {
   plot = plot + geom_shadowtext(aes(x=Organism, y=df.sig[[hist.parameter]], 
                                     label='-----',group=Treatment),bg.r=.08, 
                                 bg.colour='white',color='black', 
                                 position=position_dodge2(0.9))
  } else{plot = plot}
  plot
}

size.plot <- barplot.aov(df.sig = size.sig, y_axis='Cell size (µM^2)')
size.plot

#-----------Growth rate and Chla----------

physio<-read.csv('../rawData/physio_exp.csv')

physio$culture <- gsub( "4", "C. closterium UGA4", physio$culture)
physio$culture <- gsub( "8", "C. closterium UGA8", physio$culture)
physio$culture <- gsub( "6", "G. oceanica", physio$culture)
physio$culture <- gsub( "13", "G. huxleyi", physio$culture)

#add an empty row for missing treatment in uga.04
empty.04 <- data.frame("culture"="C. closterium UGA4", "treatment"="High Fe", 
                       "Sample"="04pfe19", "GrowthRate"=0, "chla.ug.L"=0,
                       "CellSize"=0, "Fv.Fm"=0, "pH.change"=0, "cells.L"=0)

physio<- bind_rows(physio,empty.04)
colnames(physio)[1:2] <- c('Organism','Treatment')

for (i in 1:nrow(physio)) {
  if (physio$Organism[i] == "G. oceanica"|physio$Organism[i] == "C. closterium UGA8"){
    physio$Shelf[i] = "Inner shelf"} 
  else{physio$Shelf[i] = "Outer shelf"}
}

for (i in 1:nrow(physio)) {
  if (physio$Organism[i] == "G. oceanica"|physio$Organism[i] == "G. huxleyi"){
    physio$taxa[i] = "Coccolithophore"} 
  else{physio$taxa[i] = "Diatom"}
}

## Seperating data by measurement, growth rate, chla, etc.

growth <- select(physio,c(Organism,Treatment,GrowthRate,chla.ug.L, cells.L,taxa, Shelf,Sample)) %>%
  drop_na(c(GrowthRate)) 

growth$treatment <- factor(growth$Treatment, levels=c("High Fe","Low Fe"))
avg.mu <- growth%>%group_by(Organism,Treatment,Shelf,taxa) %>%
  summarise(mu=mean(GrowthRate),mu.sd=sd(GrowthRate))

chla <- growth %>% drop_na(chla.ug.L) #%>%
 # mutate(chla.cell=chla.ug.L/cells.L)

avg.chla <- chla %>% group_by(Organism,Treatment,Shelf,taxa) %>% 
  summarise(chla=mean(chla.ug.L),chla.sd=sd(chla.ug.L))

#####------------µ/µMAX-----#####
### calculate for high and low fe treatments PER REP then 
### combine and calculate mean, SD for each organism ###

# 1. create a df for each rep and calculate mu/muMAX

a <- filter(growth, grepl("a$|a2$",growth$Sample)==TRUE) #rep a
b <- filter(growth, grepl("b$|b2$",growth$Sample)==TRUE) #rep b
c <- filter(growth, grepl("c$|c2$",growth$Sample)==TRUE) #rep c

a <- a %>% group_by(Treatment, Organism)%>%
  summarize(GrowthRate=mean(GrowthRate)) %>% 
  ungroup() %>%
  pivot_wider(names_from = "Treatment",values_from = "GrowthRate")
# avg µ for rep a
a_mu <- data.frame(mu=(a$`Low Fe`/a$`High Fe`),Organism=a$Organism)

b <- b  %>% group_by(Treatment, Organism)%>%
  summarize(GrowthRate=mean(GrowthRate)) %>% 
  ungroup() %>%
  pivot_wider(names_from = "Treatment",values_from = "GrowthRate")
b_mu <- data.frame(mu=(b$`Low Fe`/b$`High Fe`),Organism=b$Organism)

c <- c %>% group_by(Treatment, Organism)%>%
  summarize(GrowthRate=mean(GrowthRate)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Treatment",values_from = "GrowthRate")
c_mu <- data.frame(mu=(c$`Low Fe`/c$`High Fe`),Organism=c$Organism)

# 3. combine µ's into one df
mu <- bind_rows(a_mu,b_mu,c_mu) 
mu <- mu %>% group_by(Organism)%>% 
  mutate(mu_red=((1-mu)),mu_perc = mu) 

# 2. Calculate µ/µmax
mu_muMax <- summarise(mu, 
                      m = paste(signif(mean(mu_perc),2), 
                                signif(sd(mu_perc),2), sep=' ± '))
colnames(mu_muMax) <- c("Organism", "µ/µMAX")


#### run T tests

##  NOTE: in statistical tests, COLUMNS = VARIABLES. ##

# 1. test normality 
# 2. test equal variances
# 3. run appropriate test

high.vs.low <- function(df, organism, var){
  #this function tests significance for each organism given either µ or Chla as variable to test
  df <- filter(df,Organism==eval(quote(organism))) #select for organism in df
  df <- data.frame(Treatment=df$Treatment, x = df[[var]], Organism= df$Organism)
  df$treatment <- factor(df$Treatment)
        # to run a Shapiro-wilks test, separate df by organism and iron treatment
        # to test the normality WITHIN each treatment, high and low Fe
        # where H0: the sample is normal, if the p-value is > 0.05, you can run a t test,
        # otherwise if p-value < 0.05, we reject the H0 and run a Wilcox rank test
  high.fe <- filter(df, grepl("High Fe",Treatment)==TRUE)
  low.fe <- filter(df,grepl("Low Fe",Treatment)==TRUE)
  
  if (nrow(low.fe) >= 3){ #run normality check if there are enough rows
    normal.h <-shapiro_test(high.fe$x) #here x is the variable for the organism either growth rate or chla
    normal.l <- shapiro_test(low.fe$x)
    
    if (normal.h$p.value < 0.05 || normal.l$p.value < 0.05){ #if not normal,...
      print("Sample not normal, run a wilcoxon rank test")
      stat.test <- wilcox_test(df,x~Treatment) %>% 
        adjust_pvalue(method='fdr') 
      stat.test
      }else{
      print("Sample is normal, run an f test")
      f.test <- var.test(x~Treatment,df)  # to test equal variance between Treatments, we run an F test
      print(f.test)
      
      if(f.test$p.value >0.05){
        # the H0: sample variances are equal, so if p-value > 0.05, run a t test with
        # equal variances. 
        print("run a t.test with equal variance")
        stat.test <- t_test(df, x~Treatment) %>% 
          adjust_pvalue(method='fdr') 
        stat.test
        }else{
        #if p-value < 0.05, run a t test with unequal variances 
        print("run a t test with unequal variance")  
        stat.test <- t_test(df, x~Treatment, var.equal = FALSE) %>% 
          adjust_pvalue(method='fdr') 
       print(stat.test)
       }
      }
    }
  stat.test <- stat.test %>% 
    add_significance(p.col = 'p.adj', output.col = 'signif',
                     cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                     symbols = c("****", "***", "**", "*", "ns")) %>%
    mutate('Organism'=organism)
}


## Test iF growth rates significantly differ between high and low iron Treatments
mu.08 <- high.vs.low(growth,"C. closterium UGA8",'GrowthRate')  
mu.06 <- high.vs.low(growth,"G. oceanica",'GrowthRate')
mu.13 <- high.vs.low(growth,"G. huxleyi",'GrowthRate') #all samples normal
mu.signif <- bind_rows(mu.13,mu.06,mu.08)
colnames(mu.signif)[2] = c('Treatment')


mu.signif <- left_join(avg.mu, mu.signif, by=c('Organism','Treatment'))

## Test if chl. a conc significantly differs between high and low iron treatments
chla.06 <-high.vs.low(chla,"G. oceanica",'chla.ug.L') #normal, equal variance
chla.13 <-high.vs.low(chla,"G. huxleyi",'chla.ug.L') # norrmal, equal variance
chl.signif <- bind_rows(chla.13,chla.06)
colnames(chl.signif)[2] = c('Treatment')

chl.signif <- left_join(avg.chla, chl.signif, by=c('Organism','Treatment'))



#####------ Growth Plots-------

###to add bars for historical mu, reading in older data
hist.mu <- read.csv('../rawData/histData/histGrowth.csv')
hist.mu$Organism <- str_replace(hist.mu$Organism, 'E. huxleyi','G. huxleyi')
hist.av.mu <- hist.mu %>% group_by(Organism,Treatment,Taxa)%>%
  summarize('h.growth.sd'=sd(GrowthRate),'h.growth'=mean(GrowthRate))

hist.av.mu <- hist.av.mu %>% group_by(Treatment, Organism, Taxa)
colnames(hist.av.mu)[3] <- 'taxa'
mu.signif <- left_join(mu.signif, hist.av.mu, by=c('Organism','Treatment','taxa'))


mu.plot <-{
  avg.mu$Organism <- factor(avg.mu$Organism, 
                            levels=c('G. oceanica','G. huxleyi',
                                     'C. closterium UGA8','C. closterium UGA4'))
  brp <- ggplot(mu.signif, aes(Organism, mu, fill=Treatment))+
    geom_bar(color='black',stat='identity', position = position_dodge2(0.8))+
    geom_errorbar(aes(ymin=mu-mu.sd, ymax=mu+mu.sd),
                  position = position_dodge(width=0.9),
                  width=0.2) +
    geom_text(aes(label=signif, y=mu+mu.sd, 
                  vjust=-1.6),size=10,
              show.legend=F) +
    geom_shadowtext(aes(label='------', y=h.growth + h.growth.sd),
                    position = position_dodge2(0.9),size=5,
                    bg.r=.08, bg.colour='white',color='black') +
    facet_grid(~taxa, scales="free_x", space="free_x", switch="x")+
    labs(y=expression(paste('Growth Rate (',day^-1,')')), fill='Treatment')+
    theme_pubr()+
    theme(axis.title.x=element_blank(), 
          axis.text.x = element_text(face="italic"),
          legend.position = 'bottom', 
          strip.placement = 'outside',
          strip.text = element_text(size=15, face='bold'),
          strip.background = element_rect(fill='white', color='white'))  +
    scale_fill_manual(labels=c('High iron','Low iron'),values=c("black","white")) +
    scale_x_discrete(labels=c(
      "C. closterium UGA8" = "C. closterium UGA8\n(Inner)",
      "G. oceanica" = "G. oceanica\n(Inner)",
      "C. closterium UGA4" = "C. closterium UGA4\n(Outer)",
      "G. huxleyi" = "G. huxleyi\n(Outer)")) +
    scale_y_continuous(expand=expansion(mul=c(0,0.2)))
  print(brp)
}

mu.plot


#####µ/µMAX plots#####
#first make the tibble a ggtexttable

xxyy<- paste0("µ/µMAX") %>% 
  strwrap(width=25) %>%
  paste(collapse="\n")

table_plot <- mu_muMax %>%
  arrange(factor(Organism,levels=c('C. closterium UGA8','G. oceanica',
                                      'C. closterium UGA4','G. huxleyi'))) 

table_plot <- ggtexttable(table_plot,rows=NULL, theme = ttheme('blank', base_size=15))%>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(5), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2, column.side = "left", from.row = 2, linetype = 1)
table_plot <- table_cell_font(table_plot,row=2:tab_nrow(table_plot),column=1, face='italic',size = 15)
table_plot


mu_muMax <- mu_muMax %>%
  arrange(factor(Organism,levels=c(  'G. oceanica','G. huxleyi',
                                     'C. closterium UGA8','C. closterium UGA4')))

mu_muMax_df <- mutate(mu_muMax, taxa=c('Coccolithophore','Coccolithophore','Diatom','Diatom'))

tbs <- lapply(split(mu_muMax_df, mu_muMax_df$taxa), "[",-3)
df <- tibble(x = rep(-Inf, length(tbs)), 
             y = rep(Inf, length(tbs)), 
             taxa = levels(as.factor(mu_muMax_df$taxa)), 
             tbl = tbs)


mu_muMax_plot <-  mu.plot + 
  theme(legend.position='none', 
        axis.title.x=element_blank())  +
  geom_table(data = df, aes(x = x, y = y, label = tbl),
                 hjust=-.3,vjust = 1.5, 
             table.theme = ttheme_gtlight(),
             fill='white') 


mu_muMax_plot


#####Chla plot#####
chla.plot <-{
  chl.signif$Organism <- factor(avg.chla$Organism, 
                             levels=c('G. oceanica','G. huxleyi',
                                      'C. closterium UGA8','C. closterium UGA4'))
  brp <- ggplot(chl.signif, aes(Organism, chla, fill=Treatment))+
    geom_bar(color='black',stat='identity', position = position_dodge2(0.8))+
   geom_errorbar(aes(ymin=chla-chla.sd, ymax=chla+chla.sd),
                 position = position_dodge(width=0.9),
                     width=0.2) +
   geom_text(aes(label=signif, y=chla+chla.sd, 
                 vjust=-1.6),
             size=10,
             show.legend=FALSE) +
   facet_grid(~taxa, scales="free_x", space="free_x", switch="x")+
   labs(y='Chlorophyll a (µg/L)')+
   theme_pubr()+
    theme(axis.title.x = element_blank(), 
          axis.text.x=element_text(face="italic"),
          legend.position = 'none',
          strip.placement = 'outside',
          strip.text = element_text(size=15, face='bold'),
          strip.background = element_rect(fill='white', color='white'))  +
 scale_fill_manual(values=c("black","white")) +
   scale_x_discrete(labels=c(
     "C. closterium UGA8" = "C. closterium 8\n(Inner)",
     "G. oceanica" = "G. oceanica\n(Inner)",
     "C. closterium UGA4" = "C. closterium 4\n(Outer)",
     "G. huxleyi" = "G. huxleyi\n(Outer)")) +
   scale_y_continuous(expand=expansion(mul=c(0,0.2)))
  print(brp)
}
chla.plot

####---------------------Fv.Fm--------------------#####
#how much light hitting chloroplast is able to be converted into chemical energy
#indicates health of cell, is unhealthy not much light needed to fully saturate 
#cell and shut down photosystem


fire <- read.csv('../rawData/fire_exp.csv')

### The first column contains the organism, treatment, and replicate in one.
### Clean up data frame by removing replicate letter, extracting treatment  
### and organism code into new columns call Organism and Treatment.

#remove the replicate letter after the name of the file
fire$File <- gsub("_A|_B|_C","",fire$File)

# str_extract will take a pattern matching the treatment or organism codes and 
# move it to a new column via mutate
fire <- fire %>% 
  mutate(Treatment=(str_extract(fire$File, "pFe[0-9]*.[0-9]|Fe_add")),
         Organism=(str_extract(fire$File, "[0-9]{2}")))

# Replace organism and treatment codes with actual names
fire$Organism <- str_replace_all(fire$Organism, c('04'="C. closterium UGA4",
                                                    '06'="G. oceanica",
                                                    '08'="C. closterium UGA8",
                                                    '13'="G. huxleyi"))
fire$Treatment <- str_replace_all(fire$Treatment, c("pFe19"="High Iron",
                                                      "pFe21.9"="Low Iron",
                                                      "Fe_add"="Fe Ammendment"))

# set order of treatments to appear in figure
fire$Treatment <- factor(fire$Treatment,levels = c("High Iron","Low Iron","Fe Ammendment"))

##### calculate summary stats for table in output/ 
stat.fire <- {
  fire$Fv.Fm <- as.numeric(fire$Fv.Fm)
  fire %>% group_by(Organism,Treatment) %>% 
    summarise(AvgFv.Fm = mean(Fv.Fm), SdFv.Fm = sd(Fv.Fm),
              AvgtQa = mean(TauAv1),SdtQa = sd(TauAv1),
              AvgSigma = mean(Sigma), SdSigma = sd(Sigma))
}

fv.fm.anova.04 <- anova.flex(df = fire, organism = "C. closterium UGA4",var = 'Fv.Fm')
fv.fm.anova.08 <- anova.flex(fire, "C. closterium UGA8",'Fv.Fm')
fv.fm.anova.06 <- anova.flex(fire, "G. oceanica",'Fv.Fm')
fv.fm.anova.13 <- anova.flex(fire, "G. huxleyi",'Fv.Fm')

# 6. Combine anova results
fv.fm.sig <- rbind(fv.fm.anova.08, fv.fm.anova.06, fv.fm.anova.04, fv.fm.anova.13)
fv.fm.sig <- group_by(fv.fm.sig, Organism, Treatment) 

#####read in and prep historical  fire#####
histFire <- read.csv('../rawData/histData/histFire.csv')

# str_extract will take a pattern matching the treatment or organism codes and 
# move it to a new column via mutate
histFire <- histFire %>% 
  mutate(Treatment=(str_extract(histFire$File, "pFe[0-9]*.[0-9]")),
         Organism=(str_extract(histFire$File, "[0-9]{2}")))

# Replace organism and treatment codes with actual names
histFire$Organism <- str_replace_all(histFire$Organism, c('04'="C. closterium UGA4",
                                                  '06'="G. oceanica",
                                                  '08'="C. closterium UGA8",
                                                  '13'="G. huxleyi"))
histFire$Treatment <- str_replace_all(histFire$Treatment, c("pFe19"="High Iron",
                                                    "pFe21.9"="Low Iron"))

#summarize historic fire data
hist_table <- histFire %>% group_by(Organism, Treatment) %>%
    summarise(h.AvgFv.Fm=mean(Fv.Fm),h.sdFv.Fm=sd(Fv.Fm),
              h.AvgtQa=mean(TauAv1),h.sdtQa=sd(TauAv1),
              h.AvgSigma=mean(Sigma),h.sdSigma=sd(Sigma))
 
# combine anova output for fv.fm with historic data for plotting
fv.fm.sig <- full_join(fv.fm.sig, hist_table)

# set order of organisms and treatments to desired for figure
fv.fm.sig$Organism <- factor(fv.fm.sig$Organism,
                           levels = c("G. oceanica","G. huxleyi",
                                      "C. closterium UGA8", "C. closterium UGA4"))

fv.fm.sig$Treatment <- factor(fv.fm.sig$Treatment,
                            levels = c("High Iron","Low Iron","Iron Amendment"))

# add shelf zone for each isolate
for (i in 1:nrow(fv.fm.sig)) {
  if (fv.fm.sig$Organism[i] == "G. oceanica"|fv.fm.sig$Organism[i] == "C. closterium UGA8"){
    fv.fm.sig$Shelf[i] = "Inner shelf"} 
  else{fv.fm.sig$Shelf[i] = "Outer shelf"}
}

# add taxa for each isolate
for (i in 1:nrow(fv.fm.sig)) {
  if (fv.fm.sig$Organism[i] == "G. oceanica"|fv.fm.sig$Organism[i] == "G. huxleyi"){
    fv.fm.sig$taxa[i] = "Coccolithophore"} 
  else{fv.fm.sig$taxa[i] = "Diatom"}
}
######------Fv/Fm Plot----------

Fv.Fm.Plot <- barplot.aov(fv.fm.sig,hist.parameter = 'h.AvgFv.Fm',y_axis = 'Fv.Fm')
Fv.Fm.Plot

####TauAv1####
#Now I can do the same with the tQA data using the TauAv1 column

tau.04 <- anova.flex(df = fire, organism = "C. closterium UGA4",var = 'TauAv1')
tau.08 <- anova.flex(fire, "C. closterium UGA8",'TauAv1')
tau.13 <- anova.flex(fire, "G. huxleyi",'TauAv1')
tau.06 <- anova.flex(fire, "G. oceanica",'TauAv1')

#####combine exp tau with historical #####

tQa.sig <- rbind(tau.08,tau.13,tau.04,tau.06)
tQa.sig <- group_by(tQa.sig,Organism, Treatment)
tQa.sig <- full_join(tQa.sig, hist_table)
tQa.sig$Organism <- factor(tQa.sig$Organism,levels = c("C. closterium UGA8",
                                          "G. oceanica","C. closterium UGA4",
                                          "G. huxleyi"))


for (i in 1:nrow(tQa.sig)) {
  if (tQa.sig$Organism[i] == "G. oceanica"|tQa.sig$Organism[i] == "C. closterium UGA8"){
    tQa.sig$Shelf[i] = "Inner shelf"} 
  else{tQa.sig$Shelf[i] = "Outer shelf"}
}

for (i in 1:nrow(tQa.sig)) {
  if (tQa.sig$Organism[i] == "G. oceanica"|tQa.sig$Organism[i] == "G. huxleyi"){
    tQa.sig$taxa[i] = "Coccolithophore"} 
  else{tQa.sig$taxa[i] = "Diatom"}
}

tQa.sig$Treatment <- factor(tQa.sig$Treatment,
                                levels = c("High Iron","Low Iron","Iron Amendment"))

######--------tQa Plot-------

Tau.av.Plot <- barplot.aov(tQa.sig, hist.parameter = 'h.AvgtQa',y_axis=expression(paste(tau,"Q"[a],' (µs)')))

Tau.av.Plot

#### Sigma ####
## functional absorption cross-section of PSII ##

##### RUN ONE-WAY ANOVA #####

sigma.04 <- anova.flex(df = fire, organism = "C. closterium UGA4", var = 'Sigma')
sigma.08 <- anova.flex(fire, "C. closterium UGA8", "Sigma")
sigma.06 <- anova.flex(fire, "G. oceanica", 'Sigma')
sigma.13 <- anova.flex(fire, "G. huxleyi", 'Sigma')

#6. once this is done for all organisms, bind rows of data frames

Sigma.sig <- rbind(sigma.08,sigma.13,sigma.04,sigma.06)
Sigma.sig <- group_by(Sigma.sig,Organism, Treatment) 
Sigma.sig <- full_join(Sigma.sig, hist_table)
Sigma.sig$Organism <- factor(Sigma.sig$Organism,
                           levels = c("C. closterium UGA8","G. oceanica",
                                      "C. closterium UGA4","G. huxleyi"))

for (i in 1:nrow(Sigma.sig)) {
  if (Sigma.sig$Organism[i] == "G. oceanica"|Sigma.sig$Organism[i] == "C. closterium 8"){
    Sigma.sig$Shelf[i] = "Inner shelf"} 
  else{Sigma.sig$Shelf[i] = "Outer shelf"}
}

for (i in 1:nrow(Sigma.sig)) {
  if (Sigma.sig$Organism[i] == "G. oceanica"|Sigma.sig$Organism[i] == "G. huxleyi"){
    Sigma.sig$taxa[i] = "Coccolithophore"} 
  else{Sigma.sig$taxa[i] = "Diatom"}
}

##### Sigma plot #####

Sigma.Plot <- barplot.aov(df.sig = Sigma.sig,hist.parameter = 'h.AvgSigma', y_axis = expression(paste(sigma[PSII])))
Sigma.Plot

###------NPQ-------

npq <- read.csv('../rawData/als_exp.csv')

## Clean up data frame, making treatment names and shelf names consistent
npq$Treatment <- str_replace_all(npq$Treatment, c('Iron Add-Back'='Iron Amendment',
                                                  'Low Fe'='Low Iron','High Fe'='High Iron'))
npq$Shelf <- str_replace(npq$Shelf, 'Inner shelf', 'Inner Shelf ')

## Pull out NPQ from the Low Iron treatment and summarize 
npqLow <- npq %>% filter(Treatment=='Low Iron')
npqLow <- npqLow %>% group_by(Culture, Shelf, PAR, Taxa) %>%
    summarise('npq.sd'=sd(NPQ),'NPQ'=mean(NPQ))

## plot the NPQ vs Irradiance for the Low Iron treatment
npq.plot <- ggplot(npqLow, aes(PAR, NPQ))+
  geom_line(aes(linetype=Taxa))+
  geom_point(stat='identity')+
 geom_ribbon(aes(ymin=NPQ-npq.sd, ymax=NPQ+npq.sd, group=Taxa),
             alpha=0.3, show.legend = F)+
    facet_grid(~Shelf, switch='x')+
    labs(x=expression(paste('Irradiance (',µmol~quanta~m^-1~s^-1,')')),
         y='Non-Photochemical Quenching')+
    theme_pubr()+
    theme(strip.background = element_blank(), 
          strip.placement = 'outside',legend.position = 'right',
          legend.title = element_text(hjust = 0.05)
         )+
  scale_y_continuous(expand=expansion(mul=c(0,0.1)))+
  coord_cartesian(ylim=c(0, 2.2))
npq.plot


####------------plot all physiological data together-----------------------####
design <- "
1122
1122
1122
3344
3344
3344
5566
5566
5566
77##
7788
77##
"
all.physio.plot <- mu_muMax_plot + chla.plot + size.plot + Fv.Fm.Plot + 
  Tau.av.Plot + npq.plot + Sigma.Plot+
  guide_area() + plot_layout(design=design, guides='collect') &
  theme(legend.direction = 'horizontal')

all.physio.plot <- all.physio.plot + 
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(face='bold'))

all.physio.plot



#### Fluorescence after iron amendment ####

fluor <- read.csv('../rawData/fluor_add_back.csv')

fluor <- fluor %>% #group_by(Culture, rep) %>%
  mutate(ln_fluor=log(Fluorescence))

uga4=filter(fluor, Culture=='C. closterium UGA4') 
add_fluor4 <- ggplot(uga4,aes(x=Time, y=Fluorescence, group=rep)) +
  geom_point(stat='identity', aes(color=rep), alpha=0.3, size=2)+
  geom_line(stat='identity', aes(color=rep), alpha=0.3, linewidth=1)+
  geom_point(data=filter(uga4, Time>=9.9 & Fluorescence >= 1),stat='identity', 
             aes(color=rep), size=2)+
  geom_line(data=filter(uga4, Time>=9.9 & Fluorescence >= 1),stat='identity', 
            aes(color=rep), linewidth=1)+
  labs(title = 'C. closterium UGA4', tag='A') +
  theme(legend.position = 'none', plot.title = element_text(face='italic',hjust=0.5))

uga8 <- filter(fluor, Culture=='C. closterium UGA8')
add_fluor8 <- ggplot(uga8, aes(Time, Fluorescence, group=rep)) +
  geom_point(stat='identity', aes(color=rep), alpha=0.3, size=2)+
  geom_line(stat='identity', aes(color=rep), alpha=0.3, linewidth=1)+
  geom_point(data=filter(uga8, Time>=3.8),stat='identity', 
             aes(color=rep), size=2)+
  geom_line(data=filter(uga8, Time>=3.8),stat='identity', 
            aes(color=rep), linewidth=1)+
  labs(title='C. closterium UGA8', tag='B') +
  theme(legend.position = 'none',
        plot.title = element_text(face='italic',hjust=0.5))

uga6 <- filter(fluor, Culture=='G. oceanica')
uga6_highlight <- uga6 %>% group_by(rep) %>% slice_max(order_by=Time, n=2)

add_fluor6 <- ggplot(uga6, aes(Time, Fluorescence, group=rep)) +
  geom_point(stat='identity', aes(color=rep), size=2, alpha=0.3) +
  geom_line(stat='identity', aes(color=rep), linewidth = 1, alpha=0.3)+
  geom_point(data=uga6_highlight,stat='identity', 
             aes(color=rep), size=2)+
  geom_line(data=uga6_highlight,stat='identity', 
            aes(color=rep), linewidth=1)+ labs(title = 'G. oceanica', tag='C',color='Replicate')+
  theme(plot.title = element_text(face='italic',hjust=0.5))

uga13 <- filter(fluor, Culture=='G. huxleyi')
uga13_highlight <- uga13 %>% group_by(rep) %>% slice_max(order_by=Time, n=2)

add_fluor13 <- ggplot(uga13, aes(Time, Fluorescence, group=rep)) +
  geom_point(stat='identity', aes(color=rep), size=2, alpha=0.3) +
  geom_line(stat='identity', aes(color=rep), linewidth = 1, alpha=0.3)+
  geom_point(data=uga13_highlight,stat='identity', 
             aes(color=rep), size=2)+
  geom_line(data=uga13_highlight,stat='identity', 
            aes(color=rep), linewidth=1)+ labs(title = 'G. huxleyi', tag='D', color='Replicate')+
  theme(plot.title = element_text(face='italic',hjust=0.5))

add_fluor <- add_fluor4+add_fluor8+add_fluor6+add_fluor13 + 
  plot_layout(tag_level = 'new',guides='collect') &
  theme(legend.direction = 'vertical', axis.title = element_text(size=15), 
        plot.tag = element_text( hjust=11)) 

add_fluor
