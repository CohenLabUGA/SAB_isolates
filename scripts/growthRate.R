#### Growth rate and Chl a ####
## Lucy Quirk
## UGA Master's
#--------------------------------------------
## This script will run t-tests on growth rate and chla data and create a df
## that can be called for plotting later
#--------------------------------------------
#### Set up environment and read in data ####
#--------------------------------------------

source('functions.R')
pkgs <- c('tidyverse','rstatix')
checkAndLoadPackages(pkgs)

path_in <- c('../rawData/')
path_out <- c('../output')

dat<-read.csv(paste(path_in,'physio_exp.csv',sep=''))

#### Clean up data ####

dat$culture <- str_replace_all(
  dat$culture, 
  c('4'='C. closterium UGA4','8'='C. closterium UGA8', 
    '6'='G. oceanica','13'='G. huxleyi'))


#add an empty row for missing treatment in uga.04
empty.04 <- data.frame("culture"="C. closterium UGA4", "treatment"="High Fe", 
                       "Sample"="04pfe19", "GrowthRate"=0, "chla.ug.L"=0,
                       "CellSize"=0, "Fv.Fm"=0, "pH.change"=0, "cells.L"=0)
dat<- bind_rows(dat,empty.04)

for (i in 1:nrow(dat)) {
  if (dat$culture[i] == "G. oceanica"|dat$culture[i] == "C. closterium UGA8"){
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
                            'C. closterium UGA8','C. closterium UGA4')))
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

