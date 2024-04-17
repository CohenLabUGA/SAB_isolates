#November 9 2022 Lucy Quirk
#combining physiological data from growth experiments into one plot
#running statistical tests

library(tidyverse)
library(patchwork)
library(gridExtra)
library(rstatix)
library(ggpubr)
library(ggpmisc)
library(ggpattern)
library(multcompView)
library(RColorBrewer)
library(gt)
library(shadowtext)
library("webshot2")
library(patchwork)
library(flextable)

setwd("/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/raw.data")
path1 = "/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/Figures/Physio"
path2= "/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/output/Physio/"


#-------Cell Sizes--------
size <- read.csv('./Cell Sizes/exp_cell_size.csv') 
size=pivot_longer(size, !Culture, names_to = c('Treatment', 'rep'), 
             names_pattern = '(.*\\..*\\.)(A|B|C)', values_to = 'size')

#average by rep then by treatment
size <- size %>% group_by(Culture,rep, Treatment) %>%
  summarise(Size=mean(size), sd=sd(size))

empty.04 <- data.frame("Culture"=4, 'rep'='A', "Treatment"="High Iron", 
                       'Size'=0, 'sd'=0, 'Shelf'='Inner Shelf')

size<- bind_rows(size,empty.04)

size$Culture <- gsub( "4", "C. closterium 4", size$Culture)
size$Culture <- gsub( "8", "C. closterium 8", size$Culture)
size$Culture <- gsub( "6", "G. oceanica", size$Culture)
size$Culture <- gsub( "13", "G. huxleyi", size$Culture)

for (i in 1:nrow(size)) {
  if (size$Culture[i] == "G. oceanica"|size$Culture[i] == "C. closterium 8"){
    size$Shelf[i] = "Inner shelf"} 
  else{size$Shelf[i] = "Outer shelf"}
}


size$Treatment <- str_replace_all(size$Treatment, c('High.Iron.'='High Iron', 'Low.Iron.'='Low Iron',
                              'Iron.Amendment.'='Iron Amendment')) 
size$Treatment <- str_remove(size$Treatment, '\\.')
size <- na.omit(size)


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

cld_all <- rbind(aov.size.8, aov.size.4, aov.size.6, aov.size.13)
cld_all <- group_by(cld_all,Culture, Treatment) 


for (i in 1:nrow(cld_all)) {
  if (cld_all$Culture[i] == "G. oceanica"|cld_all$Culture[i] == "G. huxleyi"){
    cld_all$taxa[i] = "Coccolithophore"} 
  else{cld_all$taxa[i] = "Diatom"}
}

cld_all$Culture <- factor(cld_all$Culture,
                           levels = c("G. oceanica","G. huxleyi",
                                      "C. closterium 8", "C. closterium 4"))

cld_all$Treatment <- factor(cld_all$Treatment,
                            levels = c("High Iron","Low Iron","Iron Amendment"))

#####-----Cell Size Plot------

size.plot <- ggplot(
  cld_all, aes(Culture,size,fill=Treatment)) +
  geom_bar(color="black",stat='identity', size=2,position='dodge2') + 
  geom_errorbar(aes(ymin=size-sd, ymax=size+sd),
                position = position_dodge(width=0.9),
                width=0.2, size=2)+
  geom_text(aes(label=Tukey, y=size+sd, 
                vjust=-1.6,group=Treatment),
            size=c(40/.pt),
            show.legend=FALSE,
            position=position_dodge(width=0.9)) +
  facet_grid(~taxa, scales="free_x", space="free_x", switch="x") +
  theme_pubr()+
  theme(
    axis.title.x = element_blank(),
    text = element_text(size=40),
    axis.text.x = element_text(face="italic"),
    #legend.position='none',
    strip.placement = 'outside',
    strip.text = element_text(size=40, face='bold'),
    strip.background = element_rect(fill='white', color='white'), 
    legend.position = 'bottom') + 
  scale_x_discrete(labels=c(
    "C. closterium 8" = "C. closterium 8\n(Inner)",
    "G. oceanica" = "G. oceanica\n(Inner)",
    "C. closterium 4" = "C. closterium 4\n(Outer)",
    "G. huxleyi" = "G. huxleyi\n(Outer)")) +
  scale_fill_manual(values=c("black","white","grey"))+
  labs(y="Cell Size (µm)") + 
  scale_y_continuous(expand=expansion(mul=c(0,0.2)))
size.plot

ggsave('cell_size.png',size.plot, path=path1, height=30, width=30)
#-----------Growth rate and Chla----------

dat<-read.csv('physio_exp.csv')

dat$culture <- gsub( "4", "C. closterium 4", dat$culture)
dat$culture <- gsub( "8", "C. closterium 8", dat$culture)
dat$culture <- gsub( "6", "G. oceanica", dat$culture)
dat$culture <- gsub( "13", "G. huxleyi", dat$culture)

#add an empty row for missing treatment in uga.04
empty.04 <- data.frame("culture"="C. closterium 4", "treatment"="High Fe", 
                       "Sample"="04pfe19", "GrowthRate"=0, "chla.ug.L"=0,
                       "CellSize"=0, "Fv.Fm"=0, "pH.change"=0, "cells.L"=0)

dat<- bind_rows(dat,empty.04)

for (i in 1:nrow(dat)) {
  if (dat$culture[i] == "G. oceanica"|dat$culture[i] == "C. closterium 8"){
    dat$Shelf[i] = "Inner shelf"} 
  else{dat$Shelf[i] = "Outer shelf"}
}

for (i in 1:nrow(dat)) {
  if (dat$culture[i] == "G. oceanica"|dat$culture[i] == "G. huxleyi"){
    dat$taxa[i] = "Coccolithophore"} 
  else{dat$taxa[i] = "Diatom"}
}
#####PREP DATA FOR PLOTTING------
## Seperating data by measurement, growth rate, chla, etc.

growth <- select(dat,c(culture,treatment,GrowthRate,chla.ug.L, cells.L,taxa, Shelf,Sample)) %>%
  drop_na(c(GrowthRate)) 

growth$treatment <- factor(growth$treatment, levels=c("High Fe","Low Fe"))
avg.mu <- growth%>%group_by(culture,treatment,Shelf,taxa) %>%
  summarise(mu=mean(GrowthRate),mu.sd=sd(GrowthRate))

chla <- growth %>% drop_na(chla.ug.L) #%>%
 # mutate(chla.cell=chla.ug.L/cells.L)

avg.chla <- chla %>% group_by(culture,treatment,Shelf,taxa) %>% 
  summarise(chla=mean(chla.ug.L),chla.sd=sd(chla.ug.L))
avc <- avg.chla %>% group_by(culture)
filter(avc, treatment=='High Fe')

Avg_table <- {
  a <- growth %>% group_by(culture, treatment) %>% 
    summarise(GrowthRate=paste(signif(mean(GrowthRate),2), 
                             signif(sd(GrowthRate),2), sep=' ± '))
              
  b <- chla %>% group_by(culture,treatment) %>%
    summarise(Chla=paste(signif(mean(chla.ug.L),2),
                         paste(signif(sd(chla.ug.L),2)), sep=' ± '))
  ab <- left_join(a, b, by=c('culture','treatment'))
  colnames(ab) <- c("Organism","Treatment", "GrowthRate",'Chla')
  ab <- ab %>%
    arrange(factor(Organism,
                   levels=c('G. oceanica','G. huxleyi',
                            'C. closterium 8','C. closterium 4')))
}


write_csv(Avg_table, '../output/Physio/Avg_mu_chl.csv')

str_replace(Avg_table$GrowthRate, '0 ± NA','--')

#####------------µ/µMAX-----#####
### calculate for high and low fe treatments PER REP then 
### combine and calculate mean, SD for each organism ###

# 1. create a df for each rep and calculate mu/muMAX

a <- filter(growth, grepl("a$|a2$",growth$Sample)==TRUE) #rep a
b <- filter(growth, grepl("b$|b2$",growth$Sample)==TRUE) #rep b
c <- filter(growth, grepl("c$|c2$",growth$Sample)==TRUE) #rep c

a <- a %>% group_by(treatment, culture)%>%
  summarize(GrowthRate=mean(GrowthRate)) %>% 
  ungroup() %>%
  pivot_wider(names_from = "treatment",values_from = "GrowthRate")
# avg µ for rep a
a_mu <- data.frame(mu=(a$`Low Fe`/a$`High Fe`),culture=a$culture)

b <- b  %>% group_by(treatment, culture)%>%
  summarize(GrowthRate=mean(GrowthRate)) %>% 
  ungroup() %>%
  pivot_wider(names_from = "treatment",values_from = "GrowthRate")
b_mu <- data.frame(mu=(b$`Low Fe`/b$`High Fe`),culture=b$culture)

c <- c %>% group_by(treatment, culture)%>%
  summarize(GrowthRate=mean(GrowthRate)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "treatment",values_from = "GrowthRate")
c_mu <- data.frame(mu=(c$`Low Fe`/c$`High Fe`),culture=c$culture)

# 3. combine µ's into one df
mu <- bind_rows(a_mu,b_mu,c_mu) 
mu <- mu %>% group_by(culture)%>% 
  mutate(mu_red=((1-mu)),mu_perc = mu) 

# 2. Calculate µ/µmax
mu_muMax <- summarise(mu, 
                      m = paste(signif(mean(mu_perc),2), 
                                signif(sd(mu_perc),2), sep=' ± '))
colnames(mu_muMax) <- c("Organism", "µ/µMAX")
mu_muMax$Organism <- str_replace_all(mu_muMax$Organism, c('4'='UGA4','8'='UGA8'))


#### run T tests  ####

##  NOTE: in statistical tests, COLUMNS = VARIABLES. ##

# 1. test normality 
# 2. test equal variances
# 3. run appropriate test

high.vs.low <- function(df, organism, var){
  #this function tests significance for each organism given either µ or Chla as variable to test
  df <- filter(df,culture==eval(quote(organism))) #select for organism in df
  df <- data.frame(treatment=df$treatment, x = df[[var]], culture= df$culture)
  df$treatment <- factor(df$treatment)
        # to run a Shapiro-wilks test, separate df by organism and iron treatment
        # to test the normality WITHIN each treatment, high and low Fe
        # where H0: the sample is normal, if the p-value is > 0.05, you can run a t test,
        # otherwise if p-value < 0.05, we reject the H0 and run a Wilcox rank test
  high.fe <- filter(df, grepl("High Fe",treatment)==TRUE)
  low.fe <- filter(df,grepl("Low Fe",treatment)==TRUE)
  
  if (nrow(low.fe) >= 3){ #run normality check if there are enough rows
    normal.h <-shapiro_test(high.fe$x) #here x is the variable for the organism either growth rate or chla
    normal.l <- shapiro_test(low.fe$x)
    
    if (normal.h$p.value < 0.05 || normal.l$p.value < 0.05){ #if not normal,...
      print("Sample not normal, run a wilcoxon rank test")
      stat.test <- wilcox_test(df,x~treatment) %>% 
        adjust_pvalue(method='fdr') %>%
        add_significance('p.adj')
      stat.test
      }else{
      print("Sample is normal, run an f test")
      f.test <- var.test(x~treatment,df)  # to test equal variance between treatments, we run an F test
      print(f.test)
      
      if(f.test$p.value >0.05){
        # the H0: sample variances are equal, so if p-value > 0.05, run a t test with
        # equal variances. 
        print("run a t.test with equal variance")
        stat.test <- t_test(df, x~treatment) %>% 
          adjust_pvalue(method='fdr') %>%
          add_significance('p.adj')
        stat.test
        }else{
        #if p-value < 0.05, run a t test with unequal variances 
        print("run a t test with unequal variance")  
        stat.test <- t_test(df, x~treatment, var.equal = FALSE) %>% 
          adjust_pvalue(method='fdr') %>%
          add_significance('p.adj')
       print(stat.test)
       }
      }
    }else{
    print("Sample not long enough, run a single sample t test")
      
      stat.test <- high.fe %>% t_test(x~1, mu=0.01743318, alternative = "two.sided") %>% 
        adjust_pvalue(method='fdr') %>%
        add_significance('p.adj')
  }
}


## Test iF growth rates significantly differ between high and low iron treatments
mu.08 <- high.vs.low(growth,"C. closterium 8",'GrowthRate')  
mu.06 <- high.vs.low(growth,"G. oceanica",'GrowthRate')
mu.13 <- high.vs.low(growth,"G. huxleyi",'GrowthRate') #all samples normal
mu.signif <- bind_rows(mu.13,mu.06,mu.08)
mu.signif <- mutate(mu.signif,'culture'=c('G. huxleyi', 'G. oceanica', 'C. closterium 8'), 
                     '.y.'='mu','groups'=c('High Fe, Low Fe', 'High Fe, Low Fe', 'High Fe, Low Fe'),
                     'Shelf'=c('Outer Shelf', "Inner Shelf", 'Inner Shelf'))


## Test if chl. a conc significantly differs between high and low iron treatments
chla.08 <-high.vs.low(chla,"C. closterium 8",'chla.ug.L') #single sample t.test
chla.06 <-high.vs.low(chla,"G. oceanica",'chla.ug.L') #normal, equal variance
chla.13 <-high.vs.low(chla,"G. huxleyi",'chla.ug.L') # norrmal, equal variance
chl.signif <- bind_rows(chla.13,chla.06,chla.08)

chl.signif <- mutate(chl.signif,'culture'=c('G. huxleyi', 'G. oceanica', 'C. closterium 8'), 
                     '.y.'='chla.ug.L',
       'groups'=c('High Fe, Low Fe', 'High Fe, Low Fe', 'High Fe, Low Fe'),
       'x'=c(1,2,NA), 'xmin'=c(0.8,1.8,1.0), 'xmax'=c(1.2,1.2,2.0),
       'Shelf'=c('Outer Shelf', "Inner Shelf", 'Inner Shelf'))


#####------ Growth and Chla Plots-------

## use rstatix for:
#   - calculating significance 
#   - add_xy_position
##  use ggpubr to:
#   - create boxplot 
#   - add_xy_position 
#   - stat_pvalue_position
##facet by culture to compare growth rate/chla at each treatment

###to add bars for historica mu, reading in older data
hist.mu <- read.csv('./Growth-Rate/hist.growth.csv')
hist.mu$Organism <- str_replace(hist.mu$Organism, 'E. huxleyi','G. huxleyi')
hist.av.mu <- hist.mu %>% group_by(Organism,Treatment,Taxa)%>%
  summarize('sd'=sd(GrowthRate),'GrowthRate'=mean(GrowthRate))
colnames(hist.av.mu)[1:3] <- c('culture', 'treatment','taxa')
hist.av.mu <- hist.av.mu %>% group_by(treatment, culture, taxa)

hist.mu.table <- hist.mu %>% group_by(Organism, Treatment) %>%
  summarise('GrowthRate'=paste(signif(mean(GrowthRate),2), 
                           signif(sd(GrowthRate),2), sep=' ± '))

colnames(hist.mu.table) <- c("Organism",'Treatment', "GrowthRate")
write.csv(hist.mu.table, '../output/Physio/hist_mu.csv')

## mean µ for 04 high fe is 1.081544 
## add a doted line at this point
mu.plot <-{
  #factor by treatment
  growth$treatment <- as.factor(growth$treatment)
  #group by culture
  df2 <- growth %>% filter(culture != "C. closterium 4") %>% 
    group_by(culture) 
  #run t test and add p_values
  stat.test <- df2 %>% t_test(GrowthRate~treatment) %>% 
    adjust_pvalue(method='fdr') %>%
    add_significance('p.adj')
  #add the position in plot for placement of significance
  stat.test <- stat.test %>% add_xy_position(dodge=0.8)
                       
 brp <- ggbarplot(growth, x='culture',y='GrowthRate', 
                  add='mean_sd', fill='treatment',
 add.params = list(group='treatment'),
 order=c("G. oceanica","G. huxleyi",'C. closterium 8','C. closterium 4'),
 position=position_dodge(0.8), legend.title="Treatment") +
  theme(axis.title = element_text(face='bold'),
        text = element_text(size=40),
        axis.text.x=element_text(face="italic"),
        legend.position = "none")+
   scale_x_discrete(labels=c(
     "C. closterium 8" = "C. closterium 8\n(Inner)",
     "G. oceanica" = "G. oceanica\n(Inner)",
     "C. closterium 4" = "C. closterium 4\n(Outer)",
     "G. huxleyi" = "G. huxleyi\n(Outer)")) +
   labs(
     x='Coccolithophore                                           Diatom   ',
     y=expression(Growth~Rate~(day^-1))) + 
   scale_fill_manual(values=c("black","white"))+
   stat_pvalue_manual(size = 9,
     stat.test, x='culture',
     hide.ns=F,
     label = "p.adj.signif", 
     label.size=c(30/.pt),
     position=position_dodge(0.9), vjust=-0.7)+
   scale_y_continuous(expand=expansion(mul=c(0,0.3)))+
   geom_shadowtext(data = hist.av.mu, 
                   aes(culture, GrowthRate,label='------',
                       group=treatment),
             position=position_dodge(width=0.8),
             bg.r=.08, bg.colour='white',color='black',
             size=6) 
 
   brp
 
   }

mu.sig <- data.frame('culture'=c('C. closterium 8','C. closterium 4','G. oceanica','G. huxleyi'),
                       'sig'=c('*',' ','**','***'),'treatment'='High Fe')
avg.mu <- left_join(avg.mu,mu.sig)

mu.plot <-{
  avg.mu$culture <- factor(avg.mu$culture, levels=c('G. oceanica','G. huxleyi',
                                                        'C. closterium 8','C. closterium 4'))
  brp <- ggplot(avg.mu, aes(culture, mu, fill=treatment))+
    geom_bar(linewidth=2, color='black',stat='identity', position = position_dodge2(0.8))+
    geom_errorbar(aes(ymin=mu-mu.sd, ymax=mu+mu.sd),
                  position = position_dodge(width=0.9),
                  width=0.2, size=2) +
    geom_text(aes(label=sig, y=mu+mu.sd, 
                  vjust=-1.6),
              size=c(60/.pt),
              show.legend=FALSE) +
    geom_shadowtext(data = hist.av.mu, 
                    aes(culture, GrowthRate,label='------'),
                    position = position_dodge2(0.9),
                    bg.r=.08, bg.colour='white',color='black',
                    size=10) +
    facet_grid(~taxa, scales="free_x", space="free_x", switch="x")+
    labs(y=expression(paste('Growth Rate (',day^-1,')')), fill='Treatment')+
    theme_pubr()+
    theme(text = element_text(size=40),
          axis.title.x=element_blank(), 
          axis.text.x = element_text(face="italic"),
          legend.position = 'bottom', 
          strip.placement = 'outside',
          strip.text = element_text(size=40, face='bold'),
          strip.background = element_rect(fill='white', color='white'))  +
    scale_fill_manual(labels=c('High iron','Low iron'),values=c("black","white")) +
    scale_x_discrete(labels=c(
      "C. closterium 8" = "C. closterium UGA8\n(Inner)",
      "G. oceanica" = "G. oceanica\n(Inner)",
      "C. closterium 4" = "C. closterium UGA4\n(Outer)",
      "G. huxleyi" = "G. huxleyi\n(Outer)")) +
    scale_y_continuous(expand=expansion(mul=c(0,0.2)))
  print(brp)
}

mu.plot

ggsave(filename="exp_mu.png",plot=mu.plot,path=path1,width=24,height=15)                   
ggsave(filename="exp_mu.pdf",plot=mu.plot,path=path1,width=20,height=15)                   

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
ggsave('muMax_table.png',table_plot, path=path1, height=5, width=10)

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
             table.theme = ttheme_gtlight(base_size = 40, padding = unit(c(1,2),'char')),
             fill='white') 


mu_muMax_plot
ggsave("mu_muMAX.png",plot=mu_muMax_plot,path=path1, width=20, height = 15)


chla.sig <- select(chl.signif, 'culture','p.adj.signif','p.adj')%>% mutate(treatment='High Fe')
avg.chla <- left_join(avg.chla,chla.sig)
#####Chla plot#####
chla.plot <-{
  avg.chla$culture <- factor(avg.chla$culture, 
                             levels=c('G. oceanica','G. huxleyi',
                                      'C. closterium 8','C. closterium 4'))
  brp <- ggplot(avg.chla, aes(culture, chla, fill=treatment))+
    geom_bar(linewidth=2, color='black',stat='identity', position = position_dodge2(0.8))+
   geom_errorbar(aes(ymin=chla-chla.sd, ymax=chla+chla.sd),
                 position = position_dodge(width=0.9),
                     width=0.2, size=2) +
   geom_text(aes(label=p.adj.signif, y=chla+chla.sd, 
                 vjust=-1.6),
             size=c(40/.pt),
             show.legend=FALSE) +
   facet_grid(~taxa, scales="free_x", space="free_x", switch="x")+
   labs(y='Chlorophyll a (µg/L)')+
   theme_pubr()+
    theme(axis.title.x = element_blank(), 
          text = element_text(size=40),
          axis.text.x=element_text(face="italic"),
          legend.position = 'none',
          strip.placement = 'outside',
          strip.text = element_text(size=40, face='bold'),
          strip.background = element_rect(fill='white', color='white'))  +
 scale_fill_manual(values=c("black","white")) +
   scale_x_discrete(labels=c(
     "C. closterium 8" = "C. closterium 8\n(Inner)",
     "G. oceanica" = "G. oceanica\n(Inner)",
     "C. closterium 4" = "C. closterium 4\n(Outer)",
     "G. huxleyi" = "G. huxleyi\n(Outer)")) +
   scale_y_continuous(expand=expansion(mul=c(0,0.2)))
  print(brp)
}
chla.plot
ggsave(filename="exp_chla.png",plot=chla.plot,path=path1,width=10,height=5)                   

ggsave(filename="exp_chla.pdf",plot=chla.plot,path=path1,width=18,height=10)                   

mu.chla.plot = ggarrange(mu_muMax_plot +
                           theme(axis.text.x = element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.title.x = element_blank(),
                                 strip.text = element_blank()),
          chla.plot, common.legend = T, legend='bottom',
          labels="AUTO", font.label = list(size=40),
          align = "v", ncol=1)
ggsave(filename="exp_mu_chla.png",plot=mu.chla.plot,path=path1,width=12,height=12) 

####---------------------Fv.Fm--------------------#####
#how much light hitting chloroplast is able to be converted into chemical energy
#indicates health of cell, is unhealthy not much light needed to fully saturate 
#cell and shut down photosystem
path.out<-"/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/output/Fire/"
path.fig<-"/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/Figures/Physio/"


read <- function(df){
  headers <- names(read.csv(paste('./FIRe/Cultures/exp-data/',df,sep=''), nrows=1))
  data <- read.csv(paste('./FIRe/Cultures/exp-data/',df,sep=''), col.names=headers)
}


fv.fm <- read('fv.fm.csv')

###### Prep data ######

#remove the rep letter after the name of the file and summarize
fv.fm$File <- gsub("_A|_B|_C","",fv.fm$File)

fv.fm <- fv.fm %>% mutate(Treatment=(str_extract(fv.fm$File, "pFe[0-9]*.[0-9]|Fe_add")),
                          Organism=(str_extract(fv.fm$File, "[0-9]{2}")))

fv.fm$Organism <- str_replace_all(fv.fm$Organism, c('04'="C. closterium 4",
                                                    '06'="G. oceanica",
                                                    '08'="C. closterium 8",
                                                    '13'="G. huxleyi"))
fv.fm$Treatment <- str_replace_all(fv.fm$Treatment, c("pFe19"="High Fe",
                                                      "pFe21.9"="Low Fe",
                                                      "Fe_add"="Fe ammendment"))
df8 <- filter(fv.fm, Organism=='C. closterium 8')
fv.fm$Treatment <- factor(fv.fm$Treatment,levels = c("High Fe","Low Fe","Fe ammendment"))

fv.fm$Treatment <- factor(fv.fm$Treatment, levels=c("High Fe","Low Fe", "Fe ammendment"))

# factor by treatment and culture and calculate summary statistics

summary_table <- function(data){
  data$Fv.Fm <- as.numeric(data$Fv.Fm)
  a <- data %>% group_by(Organism,Treatment) %>% 
  summarise(AvgFv.Fm=paste(signif(mean(Fv.Fm),2), 
                            signif(sd(Fv.Fm),2), sep=' ± '),
             AvgTau1=paste(signif(mean(TauAv1),2),
                           signif(sd(TauAv1),2), sep=' ± '))
colnames(a) <- c("Organism",'Treatment', "Fv/Fm",'tQa')
a <- a %>%
  arrange(factor(Organism,
                 levels=c('G. oceanica','G. huxleyi',
                          'C. closterium 8','C. closterium 4')))
}

stat.fire <- summary_table(fv.fm)

write_csv(stat.fire, '../output/Physio/Avg_Fire.csv')


##### RUN ONE-WAY ANOVA #####

anova.org <- function(organism){
  df <- fv.fm %>% filter(Organism==organism) 
  
  #1. run anova using aov
  anova.org <- aov(Fv.Fm~as.factor(Treatment), data=df)
  print(summary(anova.org))
  
  #2. create summary table of organism for plotting later
  summary.org <- df %>% 
    group_by(Treatment) %>%
    summarize(fv.fm = mean(Fv.Fm),sd=sd(Fv.Fm))%>%
    arrange(desc(fv.fm)) ##!!! this is important so letters from cld match right later
  
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
  summary.org <- mutate(summary.org, Organism=organism)
  summary.org
}


anova.04 <- anova.org("C. closterium 4")
anova.08 <- anova.org("C. closterium 8")
anova.06 <- anova.org("G. oceanica")
anova.13 <- anova.org("G. huxleyi")

#6. once this is done for all organisms, bind rows of data frames

cld_all <- rbind(anova.08,anova.13,anova.04,anova.06)
cld_all <- group_by(cld_all,Organism, Treatment) 



for (i in 1:nrow(cld_all)) {
  if (cld_all$Organism[i] == "G. oceanica"|cld_all$Organism[i] == "C. closterium 8"){
    cld_all$Shelf[i] = "Inner shelf"} 
  else{cld_all$Shelf[i] = "Outer shelf"}
}

for (i in 1:nrow(cld_all)) {
  if (cld_all$Organism[i] == "G. oceanica"|cld_all$Organism[i] == "G. huxleyi"){
    cld_all$taxa[i] = "Coccolithophore"} 
  else{cld_all$taxa[i] = "Diatom"}
}


#####read in and prep historical fv.fm and tqa#####
fv.fm.04 <- read('../04.fv.fm.csv')
fv.fm.08 <- read('../08.fv.fm.csv')
fv.fm.06 <- read('../06.fv.fm.csv')
fv.fm.13 <- read('../13.fv.fm.csv')


sum.fv.fm <- function(df,id,org, shelf, taxa){
  df.sum <- df %>% 
    filter(File == paste(id, 'pFe19', sep='')|File==paste(id,'pFe21.9', sep='')) %>%
    group_by(File) %>% 
    summarize(h.fv.fm=mean(Fv.Fm), sd.fv.fm=sd(Fv.Fm),
              h.tauAv1=mean(TauAv1), sd.tauAv1=sd(TauAv1)) %>%
    mutate(Treatment=c('High Fe','Low Fe'), taxa=taxa, Organism=org, Shelf=shelf)
  df.sum[,-1]
}

hist.fv.fm.04 <- sum.fv.fm(fv.fm.04, '04', 'C. closterium 4', "Outer shelf", 'Diatom') 
hist.fv.fm.08 <- sum.fv.fm(fv.fm.08, '08', 'C. closterium 8','Inner shelf', 'Diatom')
hist.fv.fm.06 <- sum.fv.fm(fv.fm.06, '06', 'G. oceanica', 'Inner shelf','Coccolithophore')
hist.fv.fm.13 <- sum.fv.fm(fv.fm.13, '13', 'G. huxleyi', 'Outer shelf','Coccolithophore')

hist.fv.fm <- bind_rows(hist.fv.fm.04, hist.fv.fm.06, hist.fv.fm.08, hist.fv.fm.13)
hist.fv.fm <- group_by(hist.fv.fm, Organism, Treatment)

hist_table <- function(df,id,org){
  a <- df %>% filter(File == paste(id, 'pFe19', sep='')|File==paste(id,'pFe21.9', sep='')) %>%
    group_by(File) %>%
    summarise(AvgFv.Fm=paste(signif(mean(Fv.Fm),2), 
                             signif(sd(Fv.Fm),2), sep=' ± '),
              AvgTau1=paste(signif(mean(TauAv1),2),
                            signif(sd(TauAv1),2), sep=' ± ')) %>%
  mutate(Treatment=c('High Fe','Low Fe'), Organism=org)
  a <- a[,c(5,4,2,3)]
  colnames(a) <- c("Organism",'Treatment', "Fv/Fm",'tQa')
 a
}

h.fv.fm.04 <- hist_table(fv.fm.04, '04', 'C. closterium 4') 
h.fv.fm.08 <- hist_table(fv.fm.08, '08', 'C. closterium 8')
h.fv.fm.06 <- hist_table(fv.fm.06, '06', 'G. oceanica')
h.fv.fm.13 <- hist_table(fv.fm.13, '13', 'G. huxleyi')
h.fv.fm <- bind_rows(h.fv.fm.04,h.fv.fm.06,h.fv.fm.08,h.fv.fm.13)

h.fv.fm <- arrange(h.fv.fm, factor(Organism,
                 levels=c('G. oceanica','G. huxleyi',
                          'C. closterium 8','C. closterium 4')))

write_csv(h.fv.fm,'../output/Physio/hist_Fire.csv')

cld_all <- full_join(cld_all, hist.fv.fm)

cld_all$Organism <- factor(cld_all$Organism,
                           levels = c("G. oceanica","G. huxleyi",
                                      "C. closterium 8", "C. closterium 4"))

cld_all$Treatment <- factor(cld_all$Treatment,
                            levels = c("High Fe","Low Fe","Fe ammendment"))

######------Fv/Fm Plot----------

Fv.Fm.Plot <- ggplot(
  cld_all, aes(Organism ,fv.fm,fill=Treatment)) +
  geom_bar(size=2,color="black",stat='identity', position='dodge2') + 
  geom_errorbar(aes(ymin=fv.fm-sd, ymax=fv.fm+sd),
                position = position_dodge(width=0.9),
                width=0.2, size=2)+
  geom_text(aes(label=Tukey, y=fv.fm+sd, 
                vjust=-1.8,group=Treatment),
            size=c(40/.pt),
            show.legend=FALSE,
            position=position_dodge(width=0.9)) +
  geom_shadowtext(aes(x=Organism, label='-----', y=h.fv.fm,
                group=Treatment),
                size=10,bg.r=.08, bg.colour='white',color='black',
                position=position_dodge2(0.9)) +
  facet_grid(~taxa, scales="free_x", space="free_x", switch="x") +
  theme_pubr()+
  theme(
    axis.title.x = element_blank(),
    text = element_text(size=40),
    axis.text.x = element_text(face="italic"),
    legend.position='none',
    strip.placement = 'outside',
    strip.text = element_text(size=40, face='bold'),
    strip.background = element_rect(fill='white', color='white')) + 
  scale_x_discrete(labels=c(
   "C. closterium 8" = "C. closterium 8\n(Inner)",
    "G. oceanica" = "G. oceanica\n(Inner)",
   "C. closterium 4" = "C. closterium 4\n(Outer)",
    "G. huxleyi" = "G. huxleyi\n(Outer)")) +
  scale_fill_manual(labels=c('High iron','Low iron','Iron amendment'),
                    values=c("black","white","grey"))+
  labs(y="Fv/Fm", fill='Treatment') + 
  scale_y_continuous(expand=expansion(mul=c(0,0.2)))
Fv.Fm.Plot

ggsave(filename="exp_FvFm.tiff",plot=Fv.Fm.Plot,path=path.fig,width=40,height=20)                   

ggsave(filename="exp_FvFm.png",plot=Fv.Fm.Plot,path=path.fig,width=25,height=15)                   

####TauAv1####
#Now I can do the same with the tQA data using the TauAv1 column

anova.tauav1 <- function(organism){
  df <- fv.fm %>% filter(Organism==organism) 
  
  #1. run anova using aov
  anova.org <- aov(TauAv1~as.factor(Treatment), data=df)
  print(summary(anova.org))
  
  #2. create summary table of organism for plotting later
  summary.org <- df %>% 
    group_by(Treatment) %>%
    summarize(tauAv1 = mean(TauAv1), SD=sd(TauAv1)) %>%
    arrange(desc(tauAv1)) ##!!! this is important so letters from cld match right later
  
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
  summary.org <- mutate(summary.org, Organism=organism)
  summary.org
}

tau.04 <- anova.tauav1("C. closterium 4")
tau.08 <- anova.tauav1("C. closterium 8")
tau.13 <- anova.tauav1("G. huxleyi")
tau.06 <- anova.tauav1("G. oceanica")

#####combine exp tau with historical #####

tau_cld_all <- rbind(tau.08,tau.13,tau.04,tau.06)
tau_cld_all <- group_by(tau_cld_all,Organism, Treatment)
tau_cld_all <- full_join(tau_cld_all, hist.fv.fm)
tau_cld_all$Organism <- factor(tau_cld_all$Organism,
                               levels = c("C. closterium 8",
                                          "G. oceanica",
                                          "C. closterium 4",
                                          "G. huxleyi"))


for (i in 1:nrow(tau_cld_all)) {
  if (tau_cld_all$Organism[i] == "G. oceanica"|tau_cld_all$Organism[i] == "C. closterium 8"){
    tau_cld_all$Shelf[i] = "Inner shelf"} 
  else{tau_cld_all$Shelf[i] = "Outer shelf"}
}

for (i in 1:nrow(tau_cld_all)) {
  if (tau_cld_all$Organism[i] == "G. oceanica"|tau_cld_all$Organism[i] == "G. huxleyi"){
    tau_cld_all$taxa[i] = "Coccolithophore"} 
  else{tau_cld_all$taxa[i] = "Diatom"}
}

tau_cld_all$Treatment <- factor(tau_cld_all$Treatment,
                                levels = c("High Fe","Low Fe","Fe ammendment"))

tau_cld_all$Organism <- factor(tau_cld_all$Organism,
                           levels = c("G. oceanica","G. huxleyi",
                                      "C. closterium 8","C. closterium 4"))
######--------tQa Plot-------

Tau.av.Plot <- ggplot(
  tau_cld_all, aes(Organism ,tauAv1,fill=Treatment)) +
  geom_bar(size=2,color='black', stat='identity', 
           position=position_dodge(width=0.9, preserve='single')) + 
  geom_errorbar(aes(ymin=tauAv1-SD, ymax=tauAv1+SD),
                position = position_dodge(width=0.9, preserve = 'single'),
                width=0.2, size=2)+
  theme_pubr()+
  geom_shadowtext(aes(x=Organism, label='-----', y=h.tauAv1, 
                      group=Treatment),
                  size=10,bg.r=.08, bg.colour='white',color='black',
                  position=position_dodge2(0.9))+
  geom_text(aes(label=Tukey, y=tauAv1+SD, 
                vjust=-2,group=Treatment),
            show.legend=FALSE, size=c(50/.pt),
            position=position_dodge(width=0.9))+
  facet_grid(~taxa, scales="free_x", space="free_x", switch="x") +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size=40),
    axis.text.x=element_text(face="italic"),
    legend.position = 'none',
    strip.placement = 'outside',
    strip.text = element_text(size=40, face='bold'),
    strip.background = element_rect(fill='white', color='white')) + 
  scale_x_discrete(labels=c(
    "C. closterium 8" = "C. closterium 8\n(Inner)",
    "G. oceanica" = "G. oceanica\n(Inner)",
    "C. closterium 4" = "C. closterium 4\n(Outer)",
    "G. huxleyi" = "G. huxleyi\n(Outer)")) +
  labs(y=expression(paste(tau,"Q"[a],' (µs)')), fill='Treatment') + 
  scale_fill_manual(labels=c('High iron','Low iron','Iron amendment'),values=c("black","white","grey"))+
  scale_y_continuous(expand=expansion(mul=c(0,0.2)))

Tau.av.Plot

ggsave(filename="exp_TauAv1.tiff",plot=Tau.av.Plot,path=path.fig,width=40,height=20)                   
ggsave(filename="exp_TauAv1.png",plot=Tau.av.Plot,path=path.fig,width=25,height=15)                   

###------NPQ-------
#read in the files and group by culture, treatment, and par-level#
### experimental data ### 
setwd("/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/raw.data/FIRe/Cultures/exp-data/")
etr04 <- read.csv('04/etr.04.csv')
etr04 <- mutate(etr04, 'Taxa'='Diatom', 'Shelf'='Outer Shelf')
etr08 <- read.csv('08/etr08.csv')
etr08 <- mutate(etr08, 'Taxa'='Diatom', 'Shelf'='Inner Shelf')
etr06 <- read.csv('06/etr06.csv')
etr06 <- etr06 %>% mutate('Taxa'='Coccolithophore', 'Shelf'='Inner Shelf') 
etr13 <- read.csv('13/etr13.csv')
etr13 <- mutate(etr13, 'Taxa'='Coccolithophore', 'Shelf'='Outer Shelf')
#etr13 <- filter(etr13, rep!='b')


#etr13 <- etr13 %>% mutate('fv.fm0'=Fv.Fm/etr13$Fv.Fm[which(etr13$PAR==0)])


etr <- bind_rows(etr04,etr06,etr08,etr13)
etr$Shelf <- str_replace(etr$Shelf, '.Inner Shelf.','Inner Shelf')
etr <- etr %>% group_by(PAR,Culture,Treatment,Taxa)
etr$Culture <- str_replace_all(etr$Culture, c('Cylindrotheca uga.08'='C. closterium 8',
                                              'Cylindrotheca uga.04'='C. closterium 4',
                                              'Gephyrocapsa oceanica'='G. oceanica',
                                              'Emiliania huxleyi'='G. huxleyi')) 
etr$Culture <- factor(etr$Culture, levels=c('C. closterium 8',
                               'C. closterium 4',
                               'G. oceanica',
                               'G. huxleyi'))
etr$Treatment <- factor(etr$Treatment,
                        levels=c('High Fe','Low Fe','Iron Amendment'))


averaged <- summarise(etr, mETR=mean(ETR), TauAv1=mean(TauAv1), mFvFm=mean(Fv.Fm), sd=sd(ETR), max=max(ETR))
averaged$Treatment <- factor(averaged$Treatment,levels = c("High Fe","Low Fe","Iron Amendment"))

#####---NPQ Plot---------

x <- etr %>% filter(Treatment=='Low Fe')
x <- x%>% group_by(Culture, Shelf, PAR, Taxa) %>%
    summarise('npq.sd'=sd(NPQ),'NPQ'=mean(NPQ))

npq.plot <- ggplot(x, aes(PAR, NPQ))+
  geom_line(aes(linetype=Taxa), linewidth=3)+
  geom_point(stat='identity', size=5)+
 geom_ribbon(aes(ymin=NPQ-npq.sd, ymax=NPQ+npq.sd, group=Taxa),
             alpha=0.3, show.legend = F)+
    #scale_x_continuous(name='Irradiance(µmol quanta m^2 s^-1)')+
    facet_grid(~Shelf, switch='x')+
    labs(x=expression(paste('Irradiance (',µmol~quanta~m^-1~s^-1,')')),
         y='Non-Photochemical Quenching')+
    theme_pubr()+
    theme(text=element_text(size=40),
          strip.background = element_blank(), 
          strip.placement = 'outside',
          strip.text.x = element_text(size=40),
         legend.position = 'right',legend.title.align = (0.5), #legend.direction = 'vertical',
         legend.key.width = unit(60,'pt'))+
  scale_y_continuous(expand=expansion(mul=c(0,0.1)))+
  coord_cartesian(ylim=c(0, 2.2))
print(npq.plot)

npq.zoom.plot <- npq.plot + coord_cartesian(xlim=c(0,500), ylim=c(0,1.5))
print(npq.zoom.plot)

ggsave(filename='npq.plot.png',plot=npq.plot,path=path.fig,heigh=15, width=30)



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
#77#
"
all.physio.plot <- mu_muMax_plot + chla.plot + size.plot + Fv.Fm.Plot + 
  Tau.av.Plot + npq.plot +
  guide_area() + plot_layout(design=design, guides='collect') &
theme(legend.direction = 'horizontal',legend.text = element_text(size=50),
      legend.title = element_text(size=60), legend.key.size =unit(2,'cm'))

all.physio.plot <- all.physio.plot + 
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size=50, face='bold'))


all.physio.plot

ggsave('all_physiology_plot.png',plot=all.physio.plot,device = 'png',path=path1, 
       width=48, height=40, limitsize=F)

####---Summary Tables---####

set_flextable_defaults(font.size=8, font.family = 'Helvetica')

phys_table <- full_join(Avg_table,stat.fire) #%>% group_by(Treatment)
phys_table$Organism <- factor(phys_table$Organism, 
                             levels = c('C. closterium 8','C. closterium 4',
                                        'G. oceanica','G. huxleyi')) 
phys_table <- arrange(phys_table, Organism)
#phys_table[10,3:6] <- NA
phys_table <- pivot_longer(phys_table,cols=!c('Organism','Treatment'),names_to = 'Parameter',
                           values_to = 'mean')
phys_table$mean <-  str_replace_na(phys_table$mean, '--')
phys_table$mean <- str_replace(phys_table$mean, '0 ± NA','--')
phys_table$Treatment <- str_replace(phys_table$Treatment, 'Fe','Iron')
phys_table <- pivot_wider(phys_table, names_from = Organism, values_from = mean)
phys_table <- group_by(phys_table, Parameter)
phys_table <- phys_table[,c(2,1,3:6)]
phys_table <- phys_table%>%arrange(match(Parameter,c('GrowthRate','Chla','Fv/Fm','tQa')))
phys_table$Parameter <- str_replace_all(phys_table$Parameter, c('Chla'='Chl a (µg/L)',
                              'GrowthRate'='Growth Rate (day^-1)'))


ps <- as_grouped_data(phys_table, groups=c('Parameter'))
exp_tab <- as_flextable(ps, hide_grouplabel = T) %>% 
  style(j=1, i=~!is.na(Parameter), pr_t=fp_text_default(bold=TRUE),
        pr_p=officer::fp_par(text.align='left', padding=5)) %>%
  prepend_chunks(i = ~is.na(Parameter), j = 1, as_chunk("\t")) %>%
  set_header_labels(values=list(Treatment='',
                                'C. closterium 8'='C. closterium UGA8',
                                'C. closterium 4'='C. closterium UGA4',
                                'G. oceanica'='G. oceanica',
                                'G. huxleyi'='G. huxleyi')) %>%
  italic(part='header', j=c(2:5)) %>%
  hline(i=c(4,8,12), part='body') %>% autofit()
exp_tab

save_as_image(exp_tab, path='../../../../output/Physio/exp_table.png')



hist_table <- full_join(hist.mu.table, h.fv.fm)
hist_table$Treatment <- str_remove(hist_table$Treatment, ' Fe|Fe ')
hist_table <- pivot_longer(hist_table, cols=!c(Organism,Treatment), names_to = 'Parameter', 
                           values_to='xx')
hist_table <- pivot_wider(hist_table, names_from = Organism, values_from = xx)
hist_table <- group_by(hist_table, Parameter)
hist_table <- hist_table[,c(2,1,3:6)]
hist_table <- hist_table%>%arrange(match(Parameter,c('GrowthRate','Chla','Fv/Fm','tQa')))
hist_table$Parameter <- str_replace_all(hist_table$Parameter, c('Chla'='Chl a (µg/L)',
                                                                'GrowthRate'='Growth Rate (day^-1)'))

hs <- as_grouped_data(hist_table, groups=c('Parameter'))
hist_tab <- as_flextable(hs, hide_grouplabel = T) %>% 
    style(j=1, i=~!is.na(Parameter), pr_t=fp_text_default(bold=TRUE),
          pr_p=officer::fp_par(text.align='left', padding=5)) %>%
    prepend_chunks(i = ~is.na(Parameter), j = 1, as_chunk("\t")) %>%
    set_header_labels(values=list(Treatment='',
                                  'C. closterium 8'='C. closterium UGA8',
                                  'C. closterium 4'='C. closterium UGA4',
                                  'G. oceanica'='G. oceanica',
                                  'G. huxleyi'='G. huxleyi')) %>%
    italic(part='header', j=c(2:5)) %>%
    hline(i=c(3,6), part='body') %>% autofit()
  
hist_tab
save_as_image(hist_tab, path='../../../../output/Physio/hist_phys_table.png')


####-------Historical data: run stats and plot mu and Fv/Fm----####

remove.outliers <- function(data,var){
 
  quartiles <- quantile(data[[var]], probs=c(.25, .75), na.rm = FALSE)
  IQR <- IQR( data[[var]])

  Lower <- quartiles[1] - 1.5*IQR
  Upper <- quartiles[2] + 1.5*IQR 
  
  data_no_outlier <- subset(data,  data[[var]] > Lower &  data[[var]] < Upper)
  data_no_outlier
}



high.vs.low <- function(mu, organism){
  df <- filter(mu,Organism==eval(quote(organism)))
  df$Treatment <- factor(df$Treatment)
  high.fe <- filter(df,Treatment=="High Fe")
  print(str(high.fe))
  low.fe <- filter(df,Treatment=="Low Fe")
  #remove outliers
  high.fe <- remove.outliers(high.fe, 'GrowthRate')
  low.fe <- remove.outliers(low.fe, 'GrowthRate')
  if (nrow(low.fe) >= 3){ #run normality check if there are enough rows
    normal.h <-shapiro_test(high.fe$GrowthRate)
    normal.l <- shapiro_test(low.fe$GrowthRate)
    
    if (normal.h$p.value < 0.05 || normal.l$p.value < 0.05){ #if not normal,...
      print("Sample not normal, run a wilcoxon rank test")
      stat.test <- wilcox_test(df,GrowthRate~Treatment) %>% 
        adjust_pvalue(method='fdr') %>%
        add_significance('p.adj')
      stat.test
    }else{
      print("Sample is normal, run an f test")
      f.test <- var.test(GrowthRate~Treatment,df)  # to test equal variance between treatments, we run an F test
      print(f.test)
      
      if(f.test$p.value >0.05){
        # the H0: sample variances are equal, so if p-value > 0.05, run a t test with
        # equal variances. 
        print("run a t.test with equal variance")
        stat.test <- t_test(df, GrowthRate~Treatment) %>% 
          adjust_pvalue(method='fdr') %>%
          add_significance('p.adj')
        stat.test
      }else{
        #if p-value < 0.05, run a t test with unequal variances 
        print("run a t test with unequal variance")  
        stat.test <- t_test(df, GrowthRate~Treatment, var.equal = FALSE) %>% 
          adjust_pvalue(method='fdr') %>%
          add_significance('p.adj')
        print(stat.test)
      }
    }
  }else{
    print("Sample not long enough, run a single sample t test")
    stat.test <- high.fe %>% t_test(x~1, mu=0.01743318, alternative = "two.sided") %>% 
      adjust_pvalue(method='fdr') %>%
      add_significance('p.adj')
  }
  stat.test$Organism <- organism
  stat.test 
  }


hist04 <- high.vs.low(hist.mu,'C. closterium 4')
hist08 <- high.vs.low(hist.mu,'C. closterium 8')
hist06 <- high.vs.low(hist.mu,'G. oceanica')
hist13 <- high.vs.low(hist.mu,'G. huxleyi')

hist.mu.sig <- data.frame(
  'culture'=c('C. closterium 4','C. closterium 8','G. oceanica','G. huxleyi'),
  'sig'=c('**','***','***','NS'),'treatment'='High Fe')
hist.av.mu <- full_join(hist.av.mu,hist.mu.sig)

hist.av.mu$culture <- factor(hist.av.mu$culture, levels=c('G. oceanica','G. huxleyi',
                                                          'C. closterium 8','C. closterium 4'))
hist.mu.plot <- ggplot(hist.av.mu, aes(culture, GrowthRate, fill=treatment))+
  geom_bar(size=2, color='black',stat='identity', position = position_dodge2(0.8))+
  geom_errorbar(aes(ymin=GrowthRate-sd, ymax=GrowthRate+sd),
                position = position_dodge(width=0.9),
                width=0.2, size=2) +
  geom_text(aes(label=sig, y=GrowthRate+sd, 
                vjust=-1.6),
            size=c(40/.pt),
            show.legend=FALSE) +
  facet_grid(~taxa, scales="free_x", space="free_x", switch="x")+
  labs(y='Growth Rate (day^1)')+
  theme_pubr()+
  theme(text = element_text(size=40),
        axis.title.x=element_blank(), 
        axis.text.x = element_text(face="italic"),
        legend.position = 'none',
        strip.placement = 'outside',
        strip.text = element_text(size=40, face='bold'),
        strip.background = element_rect(fill='white', color='white'))  +
  scale_fill_manual(values=c("black","white")) +
  scale_x_discrete(labels=c(
    "C. closterium 8" = "C. closterium UGA8\n(Inner)",
    "G. oceanica" = "G. oceanica\n(Inner)",
    "C. closterium 4" = "C. closterium UGA4\n(Outer)",
    "G. huxleyi" = "G. huxleyi\n(Outer)")) +
  scale_y_continuous(expand=expansion(mul=c(0,0.2)))

hist.mu.plot
ggsave('../../Figures/Physio/his.mu.png',hist.mu.plot, height=10,width=15)




#calculate mu/muMax for historical growth rates
hist.mumax <- {
  ag <- hist.mu %>% group_by(Organism, Treatment) %>%
    summarise(m=mean(GrowthRate),sd=sd(GrowthRate))
  highm <- filter(ag, Treatment=="High Fe") 
  lowm <- filter(ag, Treatment=="Low Fe" ) 
  hist.muMax <- data.frame('µ/µMAX'=signif((lowm$m/highm$m),2))
  hist.muMax <- bind_cols(Organism =lowm$Organism, hist.muMax)
  hist.muMax <- hist.muMax[c(4,3,2,1),]
  print(hist.muMax)
  write.csv(hist.muMax, file="/Users/lquirk/Library/CloudStorage/OneDrive-UniversityofGeorgia/Thesis_g/output/Physio/hist_mu.muMax.csv")
  hist.muMax
}

hist_muMax_df <- mutate(hist.muMax, taxa=c('Coccolithophore','Coccolithophore','Diatom','Diatom'))

tbs <- lapply(split(hist_muMax_df, hist_muMax_df$taxa), "[",-3)
df <- tibble(x = rep(-Inf, length(tbs)), 
             y = rep(Inf, length(tbs)), 
             taxa = levels(as.factor(hist_muMax_df$taxa)), 
             tbl = tbs)


hist_muMax_plot <-  hist.mu.plot + 
  theme(legend.position='none', 
        axis.title.x=element_blank())  +
  geom_table(data = df, aes(x = x, y = y, label = tbl),
             hjust=-.3,vjust = 1.5, 
             table.theme = ttheme_gtlight(base_size = 40, padding = unit(c(1,2),'char')),
             fill='white') 
hist_muMax_plot

ggsave("hist_muMAX.png",plot=hist_muMax_plot,path=path1, width=10, height = 10)


#read in historical fire data

clean.df <- function(df, org, culture){
  high <- filter(df,File == paste(org,'pFe19', sep='')) %>% 
  bind_cols(Treatment="High Fe", Organism=culture)
  high <- remove.outliers(high,'Fv.Fm')
  low <- filter(df,File == paste(org,'pFe21.9', sep='')) %>% 
  bind_cols(Treatment="Low Fe", Organism=culture)
low <- remove.outliers(low, 'Fv.Fm')
df <- bind_rows(high,low) %>% 
  select(c('Organism', 'Treatment', 'Fv.Fm','File')) %>%
  group_by(Treatment) 
 
}

fv.fm4.h <- clean.df(fv.fm.04, '04', 'C. closterium UGA4')
fv.fm8.h <- clean.df(fv.fm.08, '08', 'C. closterium UGA8')
fv.fm6.h <- clean.df(fv.fm.06, '06', 'G. oceanica')
fv.fm13.h <- clean.df(fv.fm.13, '13', 'G. huxleyi')

fv.fm <- bind_rows(fv.fm4.h, fv.fm8.h, fv.fm6.h, fv.fm13.h)
ggplot(fv.fm2, aes(Treatment, Fv.Fm, color=Treatment))+
  geom_boxplot(stat='boxplot',outlier.fill = 'black')+
  facet_wrap('Organism',scales='free_x' )

fv.fm %>% group_by(Organism, Treatment) %>%
  summarise(fvfm=mean(Fv.Fm), sd=sd(Fv.Fm))

hist.fv.fm.t <- function(organism){
  df <- filter(fv.fm,Organism==eval(quote(organism)))
  df <- filter(fv.fm,Organism=='G. oceanica')
  df$Treatment <- factor(df$Treatment)
  high.fe <- filter(df,Treatment=="High Fe")
  print(str(high.fe))
  low.fe <- filter(df,Treatment=="Low Fe")
  if (nrow(low.fe) >= 3){ #run normality check if there are enough rows
    normal.h <-shapiro_test(high.fe$Fv.Fm)
    normal.l <- shapiro_test(low.fe$Fv.Fm)
    if (normal.h$p.value < 0.05 || normal.l$p.value < 0.05){ #if not normal,...
      print("Sample not normal, run a wilcoxon rank test")
      stat.test <- df %>% wilcox_test(Fv.Fm~Treatment) %>% 
        adjust_pvalue(method='fdr') %>%
        add_significance('p.adj')
      stat.test
    }else{
      print("Sample is normal, run an f test")
      f.test <- var.test(Fv.Fm~Treatment,df)  # to test equal variance between treatments, we run an F test
      print(f.test)
      if(f.test$p.value >0.05){
        # the H0: sample variances are equal, so if p-value > 0.05, run a t test with
        # equal variances. 
        print("run a t.test with equal variance")
        stat.test <- t_test(df, Fv.Fm~Treatment) %>% 
          adjust_pvalue(method='fdr') %>%
          add_significance('p.adj')
        stat.test
      }else{
        #if p-value < 0.05, run a t test with unequal variances 
        print("run a t test with unequal variance")  
        stat.test <- t_test(df, Fv.Fm~Treatment, var.equal = FALSE) %>% 
          adjust_pvalue(method='fdr') %>%
          add_significance('p.adj')
        print(stat.test)
      }
    }
  }
  stat.test$Organism <- organism
  stat.test 
}

hist.fv.fm.t('G. oceanica')
#6. once this is done for all organisms, bind rows of data frames

cld_all <- rbind(anova.08,anova.13,anova.04,anova.06)
cld_all <- group_by(cld_all,Organism, Treatment) 



for (i in 1:nrow(cld_all)) {
  if (cld_all$Organism[i] == "G. oceanica"|cld_all$Organism[i] == "C. closterium 8"){
    cld_all$Shelf[i] = "Inner shelf"} 
  else{cld_all$Shelf[i] = "Outer shelf"}
}

for (i in 1:nrow(cld_all)) {
  if (cld_all$Organism[i] == "G. oceanica"|cld_all$Organism[i] == "G. huxleyi"){
    cld_all$taxa[i] = "Coccolithophore"} 
  else{cld_all$taxa[i] = "Diatom"}
}

cld_all$Organism <- factor(cld_all$Organism,
                           levels = c("G. oceanica","G. huxleyi",
                                      "C. closterium 8", "C. closterium 4"))

Fv.Fm.HPlot <- ggplot(
  cld_all, aes(Organism ,fv.fm,fill=Treatment)) +
  geom_bar(color="black",stat='identity', position='dodge2') + 
  geom_errorbar(aes(ymin=fv.fm-sd, ymax=fv.fm+sd),
                position = position_dodge(width=0.9),
                width=0.2)+
  geom_text(aes(label=Tukey, y=fv.fm+sd, 
                vjust=-1.6,group=Treatment),
            size=c(40/.pt),
            show.legend=FALSE,
            position=position_dodge(width=0.9)) +
  facet_grid(~taxa, scales="free_x", space="free_x", switch="x") +
  theme_pubr()+
  theme(
    axis.title.x = element_blank(),
    text = element_text(size=40),
    axis.text.x = element_text(face="italic"),
    #legend.position='none',
    strip.placement = 'outside',
    strip.text = element_text(size=40, face='bold'),
    strip.background = element_rect(fill='white', color='white')) + 
  scale_x_discrete(labels=c(
    "C. closterium 8" = "C. closterium 8\n(Inner)",
    "G. oceanica" = "G. oceanica\n(Inner)",
    "C. closterium 4" = "C. closterium 4\n(Outer)",
    "G. huxleyi" = "G. huxleyi\n(Outer)")) +
  scale_fill_manual(values=c("black","white","grey"))+
  labs(y="Fv/Fm") + 
  scale_y_continuous(expand=expansion(mul=c(0,0.2)))
Fv.Fm.HPlot

ggsave('hist.Fv.Fm.png',Fv.Fm.HPlot, path=path1, height = 10, width = 10)


ggboxplot(fv.fm4.h, 'Treatment', 'Fv.Fm')
ggboxplot(fv.fm8.h, 'Treatment', 'TauAv1')##
ggboxplot(fv.fm6.h, 'Treatment', 'Fv.Fm')##
ggboxplot(fv.fm13.h, 'Treatment', 'Fv.Fm')


tau.04 <- anova.tauav1("C. closterium 4")
tau.08 <- anova.tauav1("C. closterium 8")
tau.06 <- anova.tauav1("G. huxleyi")
tau.13 <- anova.tauav1("G. oceanica")

tau_cld_all <- rbind(tau.08,tau.13,tau.04,tau.06)
tau_cld_all <- group_by(tau_cld_all,Organism, Treatment)
tau_cld_all$Organism <- factor(tau_cld_all$Organism,
                               levels = c("C. closterium 8",
                                          "G. oceanica",
                                          "C. closterium 4",
                                          "G. huxleyi"))


for (i in 1:nrow(tau_cld_all)) {
  if (tau_cld_all$Organism[i] == "G. oceanica"|tau_cld_all$Organism[i] == "C. closterium 8"){
    tau_cld_all$Shelf[i] = "Inner shelf"} 
  else{tau_cld_all$Shelf[i] = "Outer shelf"}
}

for (i in 1:nrow(tau_cld_all)) {
  if (tau_cld_all$Organism[i] == "G. oceanica"|tau_cld_all$Organism[i] == "G. huxleyi"){
    tau_cld_all$taxa[i] = "Coccolithophore"} 
  else{tau_cld_all$taxa[i] = "Diatom"}
}

tau_cld_all$Organism <- factor(cld_all$Organism,
                           levels = c("G. oceanica","G. huxleyi",
                                      "C. closterium 8", "C. closterium 4"))

Tau.HPlot <- ggplot(
  tau_cld_all, aes(Organism ,tauAv1,fill=Treatment)) +
  geom_bar(color='black', stat='identity', 
           position=position_dodge(width=0.9, preserve='single')) + 
  geom_errorbar(aes(ymin=tauAv1-SD, ymax=tauAv1+SD),
                position = position_dodge(width=0.9, preserve = 'single'),
                width=0.2)+
  theme_pubr()+
  geom_text(aes(label=Tukey, y=tauAv1+SD, 
                vjust=-2,group=Treatment),
            show.legend=FALSE, size=c(40/.pt),
            position=position_dodge(width=0.9))+
    facet_grid(~taxa, scales="free_x", space="free_x", switch="x") +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size=40),
    axis.text.x=element_text(face="italic"),
    legend.position = 'top',
    strip.placement = 'outside',
    strip.text = element_text(size=40, face='bold'),
    strip.background = element_rect(fill='white', color='white')) + 
  scale_x_discrete(labels=c(
    "C. closterium 8" = "C. closterium 8\n(Inner)",
    "G. oceanica" = "G. oceanica\n(Inner)",
    "C. closterium 4" = "C. closterium 4\n(Outer)",
    "G. huxleyi" = "G. huxleyi\n(Outer)")) +
  labs(y=expression(paste(tau,"Q"[a],' (µs)'))) + 
  scale_fill_manual(values=c("black","white","grey"))+
  scale_y_continuous(expand=expansion(mul=c(0,0.2)))

Tau.HPlot

ggsave('hist.tQa.png',Tau.HPlot, path=path1, height = 10, width = 10)


p <- ggarrange(common.legend = TRUE,legend = 'bottom',
               Fv.Fm.HPlot +theme(axis.text.x = element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.title.x = element_blank(), 
                                 strip.text = element_blank()),
               Tau.HPlot,
               labels="AUTO", font.label = list(size=40),
               nrow=2, align="v")



p
ggsave(filename="hist_fire.png",plot=p,path=path1,width=10,height=10)                   







####test hist vs exp significance####



df8 <- filter(fv.fm, Organism=='C. closterium 8') %>% mutate(kind='experiment') %>% 
  select(c('Fv.Fm','Treatment', 'kind'))
h.h8 <- filter(fv.fm.08, File=='08pFe19') %>% mutate(kind='hist', Treatment='High Fe')%>% 
  select(c('Fv.Fm','Treatment', 'kind'))
h.l8 <- filter(fv.fm.08, File=='08pFe21.9') %>% mutate(kind='hist', Treatment='Low Fe')%>% 
  select(c('Fv.Fm','Treatment', 'kind'))
dff8 <- bind_rows(df8,h.h8,h.l8)

aov8=aov(Fv.Fm~Treatment*kind, dff8)
t=tukey_hsd(aov8)
t #all NS

df4 <- filter(fv.fm, Organism=='C. closterium 4'& Treatment=='Low Fe') %>% mutate(kind='experiment') %>% 
  select(c('Fv.Fm','Treatment', 'kind'))

h.l4 <- filter(fv.fm.04, File=='04pFe21.9') %>% mutate(kind='hist', Treatment='Low Fe')%>% 
  select(c('Fv.Fm','Treatment', 'kind'))
dff4 <- bind_rows(df4, h.l4)

aov4=aov(Fv.Fm~kind, dff4)
t=tukey_hsd(aov4)
t # NS


df6 <- filter(fv.fm, Organism=='G. oceanica') %>% mutate(kind='experiment') %>% 
  select(c('Fv.Fm','Treatment', 'kind'))
h.h6 <- filter(fv.fm.06, File=='06pFe19') %>% mutate(kind='hist', Treatment='High Fe')%>% 
  select(c('Fv.Fm','Treatment', 'kind'))
h.l6 <- filter(fv.fm.06, File=='06pFe21.9') %>% mutate(kind='hist', Treatment='Low Fe')%>% 
  select(c('Fv.Fm','Treatment', 'kind'))
dff6 <- bind_rows(df6,h.h6,h.l6)

aov6=aov(Fv.Fm~Treatment*kind, dff6)
t=tukey_hsd(aov6)
t #all NS

df13<- filter(fv.fm, Organism=='G. huxleyi') %>% mutate(kind='experiment') %>% 
  select(c('Fv.Fm','Treatment', 'kind'))
h.h13 <- filter(fv.fm.13, File=='13pFe19') %>% mutate(kind='hist', Treatment='High Fe')%>% 
  select(c('Fv.Fm','Treatment', 'kind'))
h.l13 <- filter(fv.fm.13, File=='13pFe21.9') %>% mutate(kind='hist', Treatment='Low Fe')%>% 
  select(c('Fv.Fm','Treatment', 'kind'))
dff13 <- bind_rows(df13,h.h13,h.l13)

aov6=aov(Fv.Fm~Treatment*kind, dff13)
t=tukey_hsd(aov6)
t #all NS


####pH Changes####

pH <- dat %>% select(c('culture','treatment','Shelf','pH.starting','pH.change'))

ggbarplot(pH, 'treatment','pH.ending', color='treatment', facet.by = 'culture', add='mean_sd') +
  geom_count()

cph <- dat %>% filter(culture%in%c('G. oceanica','G. huxleyi')) 
g2 <- cph[str_detect(cph$Sample, '.*[abc]2'),] %>% mutate(run='second')
g1 <- filter(cph, (Sample%in%g2$Sample==F)) %>% mutate(run='exp')

gg=bind_rows(filter(g1, culture=='G. oceanica'),filter(g2, culture=='G. oceanica')) %>%
  group_by(treatment,run)

ggplot(gg, aes(run, pH.starting))+geom_count()+
  facet_wrap(~treatment)
ggbarplot(gg, 'run', 'pH.change', add='mean_sd', facet.by = 'treatment')

ee=bind_rows(filter(g1, culture=='G. huxleyi'),filter(g2, culture=='G. huxleyi')) %>%
  group_by(treatment,run)
ggplot(ee, aes(run, pH.starting))+geom_count()+
  facet_wrap(~treatment)
ggbarplot(ee, 'run', 'pH.change', add='mean_sd', facet.by = 'treatment')
